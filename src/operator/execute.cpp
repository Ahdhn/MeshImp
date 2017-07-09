/*Copyright(c) 2017, The Regents of the University of California, Davis.			*/
/*																					*/
/*																					*/
/*Redistribution and use in source and binary forms, with or without modification,	*/
/*are permitted provided that the following conditions are met :			 		*/
/*																					*/
/*1. Redistributions of source code must retain the above copyright notice, this	*/
/*list of conditions and the following disclaimer.									*/
/*2. Redistributions in binary form must reproduce the above copyright notice,		*/
/*this list of conditions and the following disclaimer in the documentation			*/
/*and / or other materials provided with the distribution.							*/
/*																					*/
/*THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND	*/
/*ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED		*/
/*WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.*/
/*IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,	*/
/*INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT	*/
/*NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR*/
/*PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,	*/
/*WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE)	*/
/*ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 		*/
/*POSSIBILITY OF SUCH DAMAGE.                                                       */
/************************************************************************************/
/************************************************************************************/

#include "execute.h"
#include "../util/Common.h"
#include "../util/KdTree.h"
#include "../util/Statistics.h"
#include "../operator/Operator.h"

//#include "DrawSpheresDebug.h" //for debugging 

#include <float.h>
#include <algorithm>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>

bool isObtuse;
bool isAcute;

void MeshImp::ErrWarnMessage(size_t lineNum, std::string message, size_t mess_id)
{
	//mess_id =0 for error (exit)
	//otherwise, it is a warning (pause)

	if (mess_id == 0){		
		fprintf(stderr, "\nError::line(%d)-->>%s", lineNum, message.c_str());				
		exit(1);
	}
	else{
		fprintf(stderr, "\nWarning::line(%d)-->>%s\n", lineNum, message.c_str());		
		system("pause");
	}
}

//** Containers **//
MeshImp::MeshImp(int numVert, double **Verts, int numTri, int**Tris) : numVert_org(numVert), numTri_org(numTri)
{
			
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "*************************** Initialize Data Structure ***************************");
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");


	//initiate containers
	MaxNumVert = size_t(3.5*numVert); //TODO only set the 1.5 factor when the injection operators are used 
	Vert_org = new vert[MaxNumVert];
	Vert_imp = new vert[MaxNumVert];
	numVert_imp = numVert_org;

	myBoundingBox.xmax = myBoundingBox.ymax = myBoundingBox.zmax = DBL_MIN;
	myBoundingBox.xmin = myBoundingBox.ymin = myBoundingBox.zmin = DBL_MAX;


	//1) store vertices coordinates and set bounding  box
	for (int i = 0; i < numVert_org; i++){
		for (int j = 0; j < 3; j++){
			Vert_imp[i].x[j] = Vert_org[i].x[j] = Verts[i][j];
		}		
		Vert_imp[i].connect[0] = Vert_org[i].connect[0] = 0;

		myBoundingBox.xmax = std::max(myBoundingBox.xmax, Verts[i][0]);
		myBoundingBox.xmin = std::min(myBoundingBox.xmin, Verts[i][0]);
		myBoundingBox.ymax = std::max(myBoundingBox.ymax, Verts[i][1]);
		myBoundingBox.ymin = std::min(myBoundingBox.ymin, Verts[i][1]);
		myBoundingBox.zmax = std::max(myBoundingBox.zmax, Verts[i][2]);
		myBoundingBox.zmin = std::min(myBoundingBox.zmin, Verts[i][2]);
	}

	myBoundingBox.lx = myBoundingBox.xmax - myBoundingBox.xmin;
	myBoundingBox.ly = myBoundingBox.ymax - myBoundingBox.ymin;
	myBoundingBox.lz = myBoundingBox.zmax - myBoundingBox.zmin;


	//2) store  tessellation  and find neighbour triangles 
	Tri_org = new tri[numTri_org];
	for (int i = 0; i < numTri_org; i++){
		for (int j = 0; j < 3; j++){
			Tri_org[i].id[j] = Tris[i][j];
		}	
	}
	int**vertex_triangles = NULL; //for each vertex i, what are the triangles i is a head on  (this array is populated in FindNeighbourTriangles)
	FindNeighbourTriangles(vertex_triangles);

	//3) scale inside a unit box
	ScaleInputMesh();

	//4) Build KdTree for the input surface 
	surfaceKdTree.BuildTree(numVert, Vert_org, 3);

	//5) get connectivity for Vert_imp/Vert_org (fan triangle)
	GetFanTriangle();

	//6) check for isolated vertices 
	for (int i = 0; i < numVert_org; i++){
		if (Vert_org[i].connect[0] == 0){
			ErrWarnMessage(__LINE__, "MeshImp::MeshImp:: isolated vertex id=" + std::to_string(i), 1);
		}
	}

	//7) sort the connectivity 
	InitialSortAroundVertex(vertex_triangles);		
	delete[]vertex_triangles;
	

	//8)Get Statistics 
	DisplayStat(numVert_org, Vert_org, 1);
	//WriteStatToFile("input_stat.txt", numVert_org, Vert_org, 1);

	//9)Draw input mesh 
	//GetMeshOBJ("input.obj", 1, "input_obtuse.obj");
}
MeshImp::~MeshImp(){};
void MeshImp::FindNeighbourTriangles(int**&myTriangles)
{
	myTriangles = new int*[numVert_org]; //for each vertex i, temporary store the triangles where i is a head on
	int common[20];
	for (int i = 0; i < numVert_org; i++){
		myTriangles[i] = new int[20];
		myTriangles[i][0] = 0;
	}

	//loop over triangles and populate myTriangles 
	for (int t = 0; t < numTri_org; t++){
		int id0(Tri_org[t].id[0]), id1(Tri_org[t].id[1]), id2(Tri_org[t].id[2]);
		myTriangles[id0][++myTriangles[id0][0]] = t;
		myTriangles[id1][++myTriangles[id1][0]] = t;
		myTriangles[id2][++myTriangles[id2][0]] = t;

		if (myTriangles[id0][0] >= 20 || myTriangles[id1][0] >= 20 || myTriangles[id2][0] >= 20){
			ErrWarnMessage(__LINE__, "MeshImp::FindNeighbourTriangles:: a vertex is shared within more than 20 triangles!!", 0);
		}
	}

	for (int t = 0; t < numTri_org; t++){
		//for each triangle t 
		//find the other triangle shared between two of its heads
		//update the t's neighbour list 

		for (int i = 0; i < 3; i++){
			int j = (i == 2) ? 0 : i + 1;


			if (FindCommonElements(myTriangles[Tri_org[t].id[i]], myTriangles[Tri_org[t].id[j]], common)){
				if (common[0] != 2){
					ErrWarnMessage(__LINE__, "MeshImp::FindNeighbourTriangles:: Non-manifold surface", 0);
				}
				else{
					Tri_org[t].neighbour[i] = (common[1] == t) ? common[2] : common[1];
				}
			}
			else{
				ErrWarnMessage(__LINE__, "MeshImp::FindNeighbourTriangles::Non-manifold surface", 0);
			}
		}
	}


	//clean up
	//for (int i = 0; i < numVert_org; i++){ delete[] myTriangles[i]; }
	
}
void MeshImp::InitialSortAroundVertex(int**myTriangles)
{
	//Sort the triangle fan around each vertex
	//only used in the initilization stage since it depends of the input triangulation 

	int common_tri[20];
	for (int ip = 0; ip < numVert_imp; ip++){
		int iq_prv;
		for (int i = 1; i < Vert_imp[ip].connect[0]; i++){
			int iq = Vert_imp[ip].connect[i];
			//get the two triangles that share the edge ip-iq
			FindCommonElements(myTriangles[ip], myTriangles[iq],common_tri);
			if (common_tri[0] != 2){
				ErrWarnMessage(__LINE__, "MeshImp::InitialSortAroundVertex:: Input mesh is not watertight", 0);
			}

			int tri1(common_tri[1]), tri2(common_tri[2]), ik1, ik2;
			//ik1 and ik2 are the third vertex in tr1 and tri2 
			for (int j = 0; j < 3; j++){
				if (Tri_org[tri1].id[j] != ip &&Tri_org[tri1].id[j] != iq){ ik1 = Tri_org[tri1].id[j]; }
			}
			for (int j = 0; j < 3; j++){
				if (Tri_org[tri2].id[j] != ip &&Tri_org[tri2].id[j] != iq){ ik2 = Tri_org[tri2].id[j]; }
			}

			int ik; //the actual replacement 
			if (i == 1){ ik = ik1; }//anyone of them would fit (this can be further utilized to make the sorting consistent i.e., CCW or CW)
			else{
				//check 
				if (iq_prv != ik1&&iq_prv != ik2){ 
					ErrWarnMessage(__LINE__, "MeshImp::InitialSortAroundVertex:: Input mesh has invalid connnectivity #1", 0);
				}
				ik = (iq_prv == ik1) ? ik2 : ik1;
			}

			int  ik_id = GetIndex(ik, Vert_imp[ip].connect); 
			
			if (ik_id < 0){
				ErrWarnMessage(__LINE__, "MeshImp::InitialSortAroundVertex:: Input mesh has invalid connnectivity #2", 0);
			}

			std::swap(Vert_org[ip].connect[ik_id], Vert_org[ip].connect[i + 1]);
			std::swap(Vert_imp[ip].connect[ik_id], Vert_imp[ip].connect[i + 1]);



			iq_prv = iq;
		}

	}
}
void MeshImp::ScaleInputMesh()
{
	//scale input mesh inside the unit box (0,0,0) (1,1,1)

	for (int i = 0; i < numVert_org; i++){
		Vert_org[i].x[0] = Vert_imp[i].x[0] = Vert_imp[i].x[0] - myBoundingBox.xmin;
		Vert_org[i].x[1] = Vert_imp[i].x[1] = Vert_imp[i].x[1] - myBoundingBox.ymin;
		Vert_org[i].x[2] = Vert_imp[i].x[2] = Vert_imp[i].x[2] - myBoundingBox.ymin;
	}

	scale_factor = 0.0;
	scale_factor = std::max(myBoundingBox.lx, myBoundingBox.ly);
	scale_factor = std::max(scale_factor, myBoundingBox.lz);

	for (int i = 0; i < numVert_org; i++){
		Vert_org[i].x[0] = Vert_imp[i].x[0] = Vert_imp[i].x[0] / scale_factor;
		Vert_org[i].x[1] = Vert_imp[i].x[1] = Vert_imp[i].x[1] / scale_factor;
		Vert_org[i].x[2] = Vert_imp[i].x[2] = Vert_imp[i].x[2] / scale_factor;
	}

	myBoundingBox.xmax -= myBoundingBox.xmin; myBoundingBox.ymax -= myBoundingBox.ymin; myBoundingBox.zmax -= myBoundingBox.zmin;
	myBoundingBox.xmax /= scale_factor; myBoundingBox.ymax /= scale_factor; myBoundingBox.zmax /= scale_factor;
	myBoundingBox.xmin = 0.0; myBoundingBox.ymin = 0.0; myBoundingBox.zmin = 0.0;
	myBoundingBox.lx = myBoundingBox.xmax - myBoundingBox.xmin; myBoundingBox.ly = myBoundingBox.ymax - myBoundingBox.ymin; myBoundingBox.lz = myBoundingBox.zmax - myBoundingBox.zmin;

}
void MeshImp::GetFanTriangle()
{
	for (int t = 0; t < numTri_org; t++){
		for (int i = 0; i < 3; i++){
			int j = (i == 2) ? 0 : i + 1;
			int vert_i(Tri_org[t].id[i]), vert_j(Tri_org[t].id[j]);
			int entry_index = GetIndex(vert_i, Vert_imp[vert_j].connect);
			if (entry_index == -1){
				Vert_imp[vert_i].connect[++Vert_imp[vert_i].connect[0]] = Vert_org[vert_i].connect[++Vert_org[vert_i].connect[0]] = vert_j;
				Vert_imp[vert_j].connect[++Vert_imp[vert_j].connect[0]] = Vert_org[vert_j].connect[++Vert_org[vert_j].connect[0]] = vert_i;
			}

			if (Vert_imp[vert_i].connect[0] >= MAX_CONNECT || Vert_imp[vert_j].connect[0] >= MAX_CONNECT){
				ErrWarnMessage(__LINE__, "MeshImp::MeshImp:: increase MAX_CONNECT in Common.h", 0);
			}
		}
	}
}


//** Sifting/Simplification **//
void MeshImp::Simp(int targetNumSamples, int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double minAngleAllow, double maxAngleAllow, bool verbose)
{
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "****************************** Simplification Started ***************************");
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	
	Statistics myStats;
	myStats.GetAllStats(numVert_imp, Vert_imp);

	double devFactor = (theta_d < 0.0) ? -1.0 : 180 - theta_d;
	Operator SimpOpt(Vert_org);
	SimpOpt.constraints.isSmooth = isSmooth;
	SimpOpt.constraints.dev = devFactor;
	SimpOpt.constraints.isDelaunay = isDelaunay;
	
	SimpOpt.constraints.isMinAngle = true;
	SimpOpt.constraints.MinAngle = (minAngleAllow<0.0) ? myStats.GetMinAngle() : minAngleAllow;

	SimpOpt.constraints.isMaxAngle = true;
	SimpOpt.constraints.MaxAngle = (maxAngleAllow<0.0) ? myStats.GetMaxAngle() : maxAngleAllow;

	SimpOpt.constraints.isNonobtuse = false;
	SimpOpt.constraints.isEdgeLen = false;
	SimpOpt.constraints.isMaximal = false;


	int verticesHandler[4];
	verticesHandler[0] = 2;


	int itter = 0;
	double vertex_void[3];

	if (verbose){
		fprintf(stdout, "\nRemoving tri-valent vertices --> ");
	}
	int numVert_before = numVert_imp;
	auto start_time = std::chrono::high_resolution_clock::now();
	SimpOpt.TriValentRemoval(numVert_imp, Vert_imp);
	auto end_time = std::chrono::high_resolution_clock::now();

	if (verbose){
		fprintf(stdout, "%d vertex removed\n", numVert_before - numVert_imp);	
	
		myStats.GetAllStats(numVert_imp, Vert_imp);
		DisplayStat(numVert_imp, Vert_imp,0);
		fprintf(stdout, "\n itter: %i",itter);
	}
	

	std::chrono::duration<double> elapsed_time_acc;

		
	while (true){
		numVert_before = numVert_imp;
		itter++;

		auto start_time1 = std::chrono::high_resolution_clock::now();

		for (int i = 0; i < numVert_imp; i++){
			indices.push_back(i);
		}
		random_shuffle(indices.begin(), indices.end());


		//******* Ejection Two
		verticesHandler[0] = 2;
		//for (int ip1 = 0; ip1 < numVert_imp; ip1++){	
		for (int id = 0; id < numVert_imp; id++){
			int ip1 = indices[id];
			if (ip1 >= numVert_imp){ continue; }


			if (Vert_imp[ip1].connect[0] == 0) { continue; }
						
			bool ejected = false;

			verticesHandler[1] = ip1;
			for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){

				verticesHandler[2] = Vert_imp[ip1].connect[i];

				vertex_void[0] = (Vert_imp[verticesHandler[1]].x[0] + 
					              Vert_imp[verticesHandler[2]].x[0])/ 2.0;

				vertex_void[1] = (Vert_imp[verticesHandler[1]].x[1] + 
					              Vert_imp[verticesHandler[2]].x[1]) / 2.0;

				vertex_void[2] = (Vert_imp[verticesHandler[1]].x[2] + 
					              Vert_imp[verticesHandler[2]].x[2]) / 2.0;


				int closestSurfID = surfaceKdTree.FindNearest(vertex_void);
				if (SimpOpt.constraints.dev > 0){ 
					SimpOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - SimpOpt.constraints.dev;
				}
								
				if (SimpOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
					ejected = true;
					break;
				}
			}

			if (numVert_imp <= targetNumSamples){ break; }
		}
		
		if (numVert_imp > targetNumSamples){

			//******* Ejection Three
			verticesHandler[0] = 3;
			for (int id = 0; id < numVert_imp; id++){
				int ip1 = indices[id];
				if (ip1 >= numVert_imp){ continue; }

				if (Vert_imp[ip1].connect[0] == 0) { continue; }

				verticesHandler[1] = ip1;

				verticesHandler[3] = Vert_imp[ip1].connect[Vert_imp[ip1].connect[0]];

				for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){

					verticesHandler[2] = Vert_imp[ip1].connect[i];

					vertex_void[0] = (Vert_imp[verticesHandler[1]].x[0] +
						Vert_imp[verticesHandler[2]].x[0] +
						Vert_imp[verticesHandler[3]].x[0]) / 3.0;

					vertex_void[1] = (Vert_imp[verticesHandler[1]].x[1] +
						Vert_imp[verticesHandler[2]].x[1] +
						Vert_imp[verticesHandler[3]].x[1]) / 3.0;

					vertex_void[2] = (Vert_imp[verticesHandler[1]].x[2] +
						Vert_imp[verticesHandler[2]].x[2] +
						Vert_imp[verticesHandler[3]].x[2]) / 3.0;


					int closestSurfID = surfaceKdTree.FindNearest(vertex_void);
					if (SimpOpt.constraints.dev > 0){
						SimpOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - SimpOpt.constraints.dev;
					}
										
					if (SimpOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
						break;
					}

					verticesHandler[3] = verticesHandler[2];
				}

				if (numVert_imp <= targetNumSamples){ break; }
			}
		}
		
		auto end_time1 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed_time1 = end_time1 - start_time1;
				
		elapsed_time_acc = elapsed_time_acc + elapsed_time1;		
				
		
		if (verbose){
			myStats.GetAllStats(numVert_imp, Vert_imp);

			DisplayStat(numVert_imp, Vert_imp,0);

			fprintf(stdout, "\n Reducation Ratio= %f %", 100.0*(double(numVert_org - numVert_imp) / double(numVert_org)));
			
			fprintf(stdout, "\n itter: %i", itter);			
		}			

		if (numVert_imp <= targetNumSamples){ break; }
		if (numVert_before == numVert_imp){ break; }
		
	}
	
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "***************************** Simplification Finished ***************************");
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
		
	std::cout << " \nTotal elapsed time: " << elapsed_time_acc.count() << " (s)\n";
		
	GetMeshOBJ("output.obj", 0, "");
}

//** Non-obtuse Remeshing **//
void MeshImp::NonobtuseRemeshing(int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double minAngleAllow, bool verbose)
{
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "************************ Non-obtuse Remeshing Started ***************************");
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

	double devFactor = (theta_d < 0.0) ? -1.0 : 180 - theta_d;
	isObtuse = false;//no partial improvement 
	isAcute = false;

	Statistics myStats;
	myStats.GetAllStats(numVert_imp, Vert_imp); 

	
	Operator NonObtuseOpt(Vert_org);	

	NonObtuseOpt.constraints.isSmooth = isSmooth;
	NonObtuseOpt.constraints.dev = devFactor; 	
	NonObtuseOpt.constraints.isDelaunay = isDelaunay;

	NonObtuseOpt.constraints.isMinAngle = true;
	NonObtuseOpt.constraints.MinAngle = (minAngleAllow<0.0) ? myStats.GetMinAngle() : minAngleAllow;
	NonObtuseOpt.constraints.isMaxAngle = true;
	NonObtuseOpt.constraints.MaxAngle = 90.0;
	
	NonObtuseOpt.constraints.isNonobtuse = true;

	NonObtuseOpt.constraints.isEdgeLen = false;	
	NonObtuseOpt.constraints.isMaximal = false;

	
	int verticesHandler[4];
	verticesHandler[0] = 2;

	int ip2, ip3;
	int itter = 0;
	double vertex_void[3];

	if (verbose){
		fprintf(stdout, "\nRemoving tri-valent vertices --> ");
	}

	int numVert_before = numVert_imp;

	auto start_time = std::chrono::high_resolution_clock::now();
	NonObtuseOpt.TriValentRemoval(numVert_imp, Vert_imp);
	auto end_time = std::chrono::high_resolution_clock::now();

	if (verbose){
		fprintf(stdout, "%d vertex removed\n", numVert_before - numVert_imp);
	}
			
	myStats.GetAllStats(numVert_imp, Vert_imp);
	
	if (verbose){
		DisplayStat(numVert_imp, Vert_imp, 1);
		fprintf(stdout, "\n itter: %i", itter);
	}
	
	
	std::chrono::duration<double> elapsed_time_acc = end_time-start_time;

	while (myStats.GetNumNonObtuse() > 0){
		itter++;
		
		auto start_time1 = std::chrono::high_resolution_clock::now();

		for (int i = 0; i < numVert_imp; i++){
			indices.push_back(i);
		}
		random_shuffle(indices.begin(), indices.end());

		//******* Relocation 		
		for (size_t id = 0; id < indices.size(); id++){
			int i = indices[id];
			if (i >= numVert_imp){ continue; }

			//printf("\nRelocating [%d]", i);
			//fflush(stdout);
					

			int closestSurfID = surfaceKdTree.FindNearest(Vert_imp[i].x);

			//only do smoothness at smooth areas
			//i.e., areas that deviates from perfect smooth (180.0 deg) by amount equal to dev factor 
			//the area is defined by closest surface vertex to i
			
			if (NonObtuseOpt.constraints.dev > 0){
				NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
			}

			//if (ObtuseHead(i, ip2, ip3)){			
				NonObtuseOpt.Relocation(i, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer);
			//}
		
		}

		
		//******* Ejection
		verticesHandler[0] = 2;		
		for (size_t id = 0; id < indices.size(); id++){
			int ip1 = indices[id];
			if (ip1 >= numVert_imp){ continue; }
			

			if (Vert_imp[ip1].connect[0] == 0) { continue; }
			if (ObtuseHead(ip1, ip2, ip3)){		
				//we could eject either ip1, ip2 or ip3 (along with another vertex)
				bool ejected = false;

				verticesHandler[1] = ip1;				
				for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){					
				
					verticesHandler[2] = Vert_imp[ip1].connect[i];										
					vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[Vert_imp[ip1].connect[i]].x[0]) / 2.0;
					vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[Vert_imp[ip1].connect[i]].x[1]) / 2.0;
					vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[Vert_imp[ip1].connect[i]].x[2]) / 2.0;
					int closestSurfID = surfaceKdTree.FindNearest(vertex_void);
					if (NonObtuseOpt.constraints.dev >0){
						NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
					}												

					if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
						ejected = true;
						break;
					}					
				}

				if (!ejected){
					verticesHandler[1] = ip2;
					for (int i = 1; i <= Vert_imp[ip2].connect[0]; i++){
					
						verticesHandler[2] = Vert_imp[ip2].connect[i];
						vertex_void[0] = (Vert_imp[ip2].x[0] + Vert_imp[Vert_imp[ip2].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip2].x[1] + Vert_imp[Vert_imp[ip2].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip2].x[2] + Vert_imp[Vert_imp[ip2].connect[i]].x[2]) / 2.0;
						int closestSurfID = surfaceKdTree.FindNearest(vertex_void);
						if (NonObtuseOpt.constraints.dev >0){
							NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
						}
												
						if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
							ejected = true;
							break;
						}
					}
				}

				if (!ejected){
					verticesHandler[1] = ip3;
					for (int i = 1; i <= Vert_imp[ip3].connect[0]; i++){
						
						verticesHandler[2] = Vert_imp[ip3].connect[i];
						vertex_void[0] = (Vert_imp[ip3].x[0] + Vert_imp[Vert_imp[ip3].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip3].x[1] + Vert_imp[Vert_imp[ip3].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip3].x[2] + Vert_imp[Vert_imp[ip3].connect[i]].x[2]) / 2.0;
						int closestSurfID = surfaceKdTree.FindNearest(vertex_void);
						if (NonObtuseOpt.constraints.dev >0){
							NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
						}
											

						if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
							ejected = true;
							break;
						}
					}
				}
			}
		}


		//******* Injection
		verticesHandler[0] = 3;		
		for (size_t id = 0; id < indices.size(); id++){
			int ip1 = indices[id];
			if (ip1 >= numVert_imp){ continue; }

			if (Vert_imp[ip1].connect[0] == 0) { continue; }
			if (ObtuseHead(ip1, ip2, ip3)){

				//call operator on triangle (ip1, ip2,ip3)
				verticesHandler[1] = ip1;
				verticesHandler[2] = ip2;
				verticesHandler[3] = ip3;

				vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
				vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
				vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

				int closestSurfID = surfaceKdTree.FindNearest(vertex_void);

				if (NonObtuseOpt.constraints.dev >0){
					NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
				}
						

				NonObtuseOpt.Injection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, 1, numSurfaceLayer, 0);
			}
		}

		//******* Attractor Ejection		
		verticesHandler[0] = 2;		
		for (size_t id = 0; id < indices.size(); id++){
			int ip1 = indices[id];
			if (ip1 >= numVert_imp){ continue; }

			if (Vert_imp[ip1].connect[0] == 0) { continue; }
			if (ObtuseHead(ip1, ip2, ip3)){
				//we could eject either ip1, ip2 or ip3 (along with another vertex)
				bool ejected = false;

				verticesHandler[1] = ip1;
				for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){

					verticesHandler[2] = Vert_imp[ip1].connect[i];
					vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
					vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
					vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

					int closestSurfID = surfaceKdTree.FindNearest(vertex_void);
					if (NonObtuseOpt.constraints.dev >0){
						NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
					}
					
					if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
						ejected = true;
						break;
					}
				}

				if (!ejected){
					verticesHandler[1] = ip2;
					for (int i = 1; i <= Vert_imp[ip2].connect[0]; i++){

						verticesHandler[2] = Vert_imp[ip2].connect[i];
						vertex_void[0] = (Vert_imp[ip2].x[0] + Vert_imp[Vert_imp[ip2].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip2].x[1] + Vert_imp[Vert_imp[ip2].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip2].x[2] + Vert_imp[Vert_imp[ip2].connect[i]].x[2]) / 2.0;

						int closestSurfID = surfaceKdTree.FindNearest(vertex_void);
						if (NonObtuseOpt.constraints.dev >0){
							NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
						}

						if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
							ejected = true;
							break;
						}
					}
				}

				if (!ejected){
					verticesHandler[1] = ip3;
					for (int i = 1; i <= Vert_imp[ip3].connect[0]; i++){

						verticesHandler[2] = Vert_imp[ip3].connect[i];
						vertex_void[0] = (Vert_imp[ip3].x[0] + Vert_imp[Vert_imp[ip3].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip3].x[1] + Vert_imp[Vert_imp[ip3].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip3].x[2] + Vert_imp[Vert_imp[ip3].connect[i]].x[2]) / 2.0;
						int closestSurfID = surfaceKdTree.FindNearest(vertex_void);
						if (NonObtuseOpt.constraints.dev >0){
							NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
						}

						if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
							ejected = true;
							break;
						}
					}
				}
			}
		}

		//******* Repeller Injection 
		verticesHandler[0] = 3;		
		for (size_t id = 0; id < indices.size(); id++){
			int ip1 = indices[id];
			if (ip1 >= numVert_imp){ continue; }

			if (Vert_imp[ip1].connect[0] == 0) { continue; }
			if (ObtuseHead(ip1, ip2, ip3)){

				//call operator on triangle (ip1, ip2,ip3)
				verticesHandler[1] = ip1;
				verticesHandler[2] = ip2;
				verticesHandler[3] = ip3;

				vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
				vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
				vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

				int closestSurfID = surfaceKdTree.FindNearest(vertex_void);

				if (NonObtuseOpt.constraints.dev >0){
					NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
				}

				NonObtuseOpt.Injection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, 1, numSurfaceLayer, 1);				

			}
		}

		auto end_time1 = std::chrono::high_resolution_clock::now();
		
		std::chrono::duration<double> elapsed_time1 = end_time1 - start_time1;
		
		elapsed_time_acc = elapsed_time_acc + elapsed_time1;
		
		myStats.GetAllStats(numVert_imp, Vert_imp);		

		if (verbose){
			DisplayStat(numVert_imp, Vert_imp,1);
			fprintf(stdout, "\n itter: %i", itter);
		}			
	}

	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "************************ Non-obtuse Remeshing Finished **************************");
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

	std::cout << "\n Total elapsed time: " << elapsed_time_acc.count() << " (s)\n";

	GetMeshOBJ("output.obj", 0, "");
}
void MeshImp::NonobtuseRemeshing_InterleaveOpt(int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double minAngleAllow, bool verbose)
{
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "******************* Non-obtuse Remeshing - Interleave Started *******************");
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

	double devFactor = (theta_d < 0.0) ? -1.0 : 180 - theta_d;

	Statistics myStats;
	myStats.GetAllStats(numVert_imp, Vert_imp);

	//Set the constriants handler 
	Operator NonObtuseOpt(Vert_org);

	NonObtuseOpt.constraints.isSmooth = isSmooth;
	NonObtuseOpt.constraints.dev = devFactor;
	NonObtuseOpt.constraints.isDelaunay = isDelaunay;

	NonObtuseOpt.constraints.isMinAngle = true;
	NonObtuseOpt.constraints.MinAngle = (minAngleAllow<0.0) ? myStats.GetMinAngle() : minAngleAllow;
	NonObtuseOpt.constraints.isMaxAngle = true;
	NonObtuseOpt.constraints.MaxAngle = 90;

	NonObtuseOpt.constraints.isNonobtuse = true;

	NonObtuseOpt.constraints.isEdgeLen = false;
	NonObtuseOpt.constraints.isMaximal = false;


	int verticesHandler[4];
	

	int ip2, ip3;
	int itter = 0;
	double vertex_void[3];

	if (verbose){
		fprintf(stdout, "\nRemoving tri-valent vertices --> ");
	}
	int numVert_before = numVert_imp;

	auto start_time = std::chrono::high_resolution_clock::now();
	NonObtuseOpt.TriValentRemoval(numVert_imp, Vert_imp);
	auto end_time = std::chrono::high_resolution_clock::now();

	if (verbose){
		fprintf(stdout, "%d vertex removed\n", numVert_before - numVert_imp);
	}


	myStats.GetAllStats(numVert_imp, Vert_imp);

	if (verbose){
		DisplayStat(numVert_imp, Vert_imp,1);
		fprintf(stdout, "\n itter: %i", itter);
	}
	
	std::chrono::duration<double> elapsed_time_acc = end_time - start_time;

	while (myStats.GetNumNonObtuse() > 0){
		itter++;

		auto start_time1 = std::chrono::high_resolution_clock::now();
		
		for (int i = 0; i < numVert_imp; i++){
			indices.push_back(i);
		}
		random_shuffle(indices.begin(), indices.end());

		// Relocation 		
		for (size_t id = 0; id < indices.size(); id++){
			
			int ip1 = indices[id];
			if (ip1 >= numVert_imp || Vert_imp[ip1].connect[0] == 0){ continue; }
			

			int closestSurfID = surfaceKdTree.FindNearest(Vert_imp[ip1].x);
						
			if (NonObtuseOpt.constraints.dev >0){
				NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
			}

			isObtuse = ObtuseHead(ip1, ip2, ip3);
			//Relocation 			
			NonObtuseOpt.Relocation(ip1, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer);
			
			
			//test if it is still an obtuse head 
			if (ObtuseHead(ip1, ip2, ip3)){
				isObtuse = true;

				bool fixed = false;

				//Ejection
				verticesHandler[0] = 2;
				
				if (!fixed){
					verticesHandler[1] = ip1;
					for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){
						verticesHandler[2] = Vert_imp[ip1].connect[i];
						vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[Vert_imp[ip1].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[Vert_imp[ip1].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[Vert_imp[ip1].connect[i]].x[2]) / 2.0;
												
						if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
							fixed = true;
							break;
						}
					}
				}

				if (!fixed){
					verticesHandler[1] = ip2;
					for (int i = 1; i <= Vert_imp[ip2].connect[0]; i++){
						verticesHandler[2] = Vert_imp[ip2].connect[i];
						vertex_void[0] = (Vert_imp[ip2].x[0] + Vert_imp[Vert_imp[ip2].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip2].x[1] + Vert_imp[Vert_imp[ip2].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip2].x[2] + Vert_imp[Vert_imp[ip2].connect[i]].x[2]) / 2.0;
												
						if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
							fixed = true;
							break;
						}
					}
				}

				if (!fixed){
					verticesHandler[1] = ip3;
					for (int i = 1; i <= Vert_imp[ip3].connect[0]; i++){
						verticesHandler[2] = Vert_imp[ip3].connect[i];
						vertex_void[0] = (Vert_imp[ip3].x[0] + Vert_imp[Vert_imp[ip3].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip3].x[1] + Vert_imp[Vert_imp[ip3].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip3].x[2] + Vert_imp[Vert_imp[ip3].connect[i]].x[2]) / 2.0;
						
						if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
							fixed = true;
							break;
						}
					}
				}

				if (!fixed){
					//Injection
					verticesHandler[0] = 3;
					verticesHandler[1] = ip1;
					verticesHandler[2] = ip2;
					verticesHandler[3] = ip3;

					vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
					vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
					vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;
										
					if (NonObtuseOpt.Injection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, 1, numSurfaceLayer, 0)){
						fixed = true;
					}
				}

				//Attracotr Ejection
				verticesHandler[0] = 2;
				if (!fixed){					
					verticesHandler[1] = ip1;
					for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){
						verticesHandler[2] = Vert_imp[ip1].connect[i];
						vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
						vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
						vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

						if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
							fixed = true;
							break;
						}
					}
				}

				if (!fixed){
					verticesHandler[1] = ip2;
					for (int i = 1; i <= Vert_imp[ip2].connect[0]; i++){

						verticesHandler[2] = Vert_imp[ip2].connect[i];
						vertex_void[0] = (Vert_imp[ip2].x[0] + Vert_imp[Vert_imp[ip2].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip2].x[1] + Vert_imp[Vert_imp[ip2].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip2].x[2] + Vert_imp[Vert_imp[ip2].connect[i]].x[2]) / 2.0;

						if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
							fixed = true;
							break;
						}
					}
				}
				if (!fixed){
					verticesHandler[1] = ip3;
					for (int i = 1; i <= Vert_imp[ip3].connect[0]; i++){
						verticesHandler[2] = Vert_imp[ip3].connect[i];
						vertex_void[0] = (Vert_imp[ip3].x[0] + Vert_imp[Vert_imp[ip3].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip3].x[1] + Vert_imp[Vert_imp[ip3].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip3].x[2] + Vert_imp[Vert_imp[ip3].connect[i]].x[2]) / 2.0;

						if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
							fixed = true;
							break;
						}
					}
				}

				//Repeller Injection 
				if (!fixed){
					verticesHandler[0] = 3;
					verticesHandler[1] = ip1;
					verticesHandler[2] = ip2;
					verticesHandler[3] = ip3;

					vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
					vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
					vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

					NonObtuseOpt.Injection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, 1, numSurfaceLayer, 1);
				}
			}
		}
		
		auto end_time1 = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> elapsed_time1 = end_time1 - start_time1;

		elapsed_time_acc = elapsed_time_acc + elapsed_time1;

		myStats.GetAllStats(numVert_imp, Vert_imp);

		if (verbose){
			DisplayStat(numVert_imp, Vert_imp,1);
			fprintf(stdout, "\n itter: %i", itter);
		}
	}

	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "****************** Non-obtuse Remeshing - Interleave Finished *******************");
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

	std::cout << "\n Total elapsed time: " << elapsed_time_acc.count() << " (s)\n";	
	
	GetMeshOBJ("output.obj", 0, "");
}
void MeshImp::AcuteRemeshing_InterleaveOpt(int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double minAngleAllow, double maxAngleAllow, bool verbose)
{
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "********************** Acute Remeshing - Interleave Started *********************");
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

	double devFactor = (theta_d < 0.0) ? -1.0 : 180 - theta_d;

	Statistics myStats;
	myStats.GetAllStats(numVert_imp, Vert_imp);

	//Set the constriants handler 
	Operator AcuteOpt(Vert_org);

	AcuteOpt.constraints.isSmooth = isSmooth;
	AcuteOpt.constraints.dev = devFactor;
	AcuteOpt.constraints.isDelaunay = isDelaunay;

	AcuteOpt.constraints.isMinAngle = true;
	AcuteOpt.constraints.MinAngle = (minAngleAllow<0.0) ? myStats.GetMinAngle() : minAngleAllow;
	AcuteOpt.constraints.isMaxAngle = true;
	AcuteOpt.constraints.MaxAngle = maxAngleAllow;

	AcuteOpt.constraints.isNonobtuse = false;

	AcuteOpt.constraints.isEdgeLen = false;
	AcuteOpt.constraints.isMaximal = false;


	int verticesHandler[4];


	int ip2, ip3;
	int itter = 0;
	double vertex_void[3];

	if (verbose){
		fprintf(stdout, "\nRemoving tri-valent vertices --> ");
	}

	int numVert_before = numVert_imp;

	auto start_time = std::chrono::high_resolution_clock::now();
	AcuteOpt.TriValentRemoval(numVert_imp, Vert_imp);
	auto end_time = std::chrono::high_resolution_clock::now();

	if (verbose){
		fprintf(stdout, "%d vertex removed\n", numVert_before - numVert_imp);
	}


	myStats.GetAllStats(numVert_imp, Vert_imp);
	int numAcute = myStats.CalcNumAcute(numVert_imp, Vert_imp, maxAngleAllow);

	if (verbose){
		DisplayStat(numVert_imp, Vert_imp,1);
		fprintf(stdout, "\n #Acute triangle = %i", numAcute);
		fprintf(stdout, "\n itter: %i", itter);
	}


	std::chrono::duration<double> elapsed_time_acc = end_time - start_time;

	while (numAcute > 0){
		itter++;

		auto start_time1 = std::chrono::high_resolution_clock::now();

		for (int i = 0; i < numVert_imp; i++){
			indices.push_back(i);
		}
		random_shuffle(indices.begin(), indices.end());

		// Relocation 		
		for (size_t id = 0; id < indices.size(); id++){
			int ip1 = indices[id];
			if (ip1 >= numVert_imp || Vert_imp[ip1].connect[0] == 0){ continue; }


			int closestSurfID = surfaceKdTree.FindNearest(Vert_imp[ip1].x);

			if (AcuteOpt.constraints.dev > 0){
				AcuteOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - AcuteOpt.constraints.dev;
			}

			isAcute = AcuteHead(ip1, ip2, ip3, maxAngleAllow);

			//Relocation 			
			AcuteOpt.Relocation(ip1, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer);

			//test if it is still an obtuse head 
			if (AcuteHead(ip1, ip2, ip3, maxAngleAllow)){
				isAcute = false;


				bool fixed = false;

				//Ejection
				verticesHandler[0] = 2;

				if (!fixed){
					verticesHandler[1] = ip1;
					for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){
						verticesHandler[2] = Vert_imp[ip1].connect[i];
						vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[Vert_imp[ip1].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[Vert_imp[ip1].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[Vert_imp[ip1].connect[i]].x[2]) / 2.0;
						if (AcuteOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
							fixed = true;
							break;
						}
					}
				}

				if (!fixed){
					verticesHandler[1] = ip2;
					for (int i = 1; i <= Vert_imp[ip2].connect[0]; i++){
						verticesHandler[2] = Vert_imp[ip2].connect[i];
						vertex_void[0] = (Vert_imp[ip2].x[0] + Vert_imp[Vert_imp[ip2].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip2].x[1] + Vert_imp[Vert_imp[ip2].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip2].x[2] + Vert_imp[Vert_imp[ip2].connect[i]].x[2]) / 2.0;
						if (AcuteOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
							fixed = true;
							break;
						}
					}
				}

				if (!fixed){
					verticesHandler[1] = ip3;
					for (int i = 1; i <= Vert_imp[ip3].connect[0]; i++){
						verticesHandler[2] = Vert_imp[ip3].connect[i];
						vertex_void[0] = (Vert_imp[ip3].x[0] + Vert_imp[Vert_imp[ip3].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip3].x[1] + Vert_imp[Vert_imp[ip3].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip3].x[2] + Vert_imp[Vert_imp[ip3].connect[i]].x[2]) / 2.0;
						if (AcuteOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
							fixed = true;
							break;
						}
					}
				}

				if (!fixed){
					//Injection
					verticesHandler[0] = 3;
					verticesHandler[1] = ip1;
					verticesHandler[2] = ip2;
					verticesHandler[3] = ip3;

					vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
					vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
					vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

					if (AcuteOpt.Injection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, 1, numSurfaceLayer, 0)){
						fixed = true;
					}
				}

				//Attracotr Ejection
				verticesHandler[0] = 2;
				if (!fixed){
					verticesHandler[1] = ip1;
					for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){
						verticesHandler[2] = Vert_imp[ip1].connect[i];
						vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
						vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
						vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;
						if (AcuteOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
							fixed = true;
							break;
						}
					}
				}

				if (!fixed){
					verticesHandler[1] = ip2;
					for (int i = 1; i <= Vert_imp[ip2].connect[0]; i++){

						verticesHandler[2] = Vert_imp[ip2].connect[i];
						vertex_void[0] = (Vert_imp[ip2].x[0] + Vert_imp[Vert_imp[ip2].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip2].x[1] + Vert_imp[Vert_imp[ip2].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip2].x[2] + Vert_imp[Vert_imp[ip2].connect[i]].x[2]) / 2.0;
						if (AcuteOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
							fixed = true;
							break;
						}
					}
				}
				if (!fixed){
					verticesHandler[1] = ip3;
					for (int i = 1; i <= Vert_imp[ip3].connect[0]; i++){
						verticesHandler[2] = Vert_imp[ip3].connect[i];
						vertex_void[0] = (Vert_imp[ip3].x[0] + Vert_imp[Vert_imp[ip3].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip3].x[1] + Vert_imp[Vert_imp[ip3].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip3].x[2] + Vert_imp[Vert_imp[ip3].connect[i]].x[2]) / 2.0;
						if (AcuteOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
							fixed = true;
							break;
						}
					}
				}

				//Repeller Injection 
				if (!fixed){
					verticesHandler[0] = 3;
					verticesHandler[1] = ip1;
					verticesHandler[2] = ip2;
					verticesHandler[3] = ip3;

					vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
					vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
					vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

					AcuteOpt.Injection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, 1, numSurfaceLayer, 1);
				}
			}
		}

		auto end_time1 = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> elapsed_time1 = end_time1 - start_time1;

		elapsed_time_acc = elapsed_time_acc + elapsed_time1;

		myStats.GetAllStats(numVert_imp, Vert_imp);
		numAcute = myStats.CalcNumAcute(numVert_imp, Vert_imp, maxAngleAllow);

		if (verbose){
			DisplayStat(numVert_imp, Vert_imp, 1);
			fprintf(stdout, "\n #Acute triangle = %i", numAcute);
			fprintf(stdout, "\n itter: %i", itter);
		}
	}

	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "********************** Acute Remeshing - Interleave Finished ********************");
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

	std::cout << "\n Total elapsed time: " << elapsed_time_acc.count() << " (s)\n";
	
	GetMeshOBJ("output.obj", 0, "");
}

//** Small Angle Elimination **//
void MeshImp::SmallAngleElimination_InterleaveOpt(double targetMinAngle, int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double maxAngleAllow, bool verbose)
{
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "****************** Small Angle Elimination - Interleave Started *****************");
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

	double devFactor = (theta_d < 0.0) ? -1.0 : 180 - theta_d;

	Statistics myStats;
	myStats.GetAllStats(numVert_imp, Vert_imp);

	//Set the constriants handler 
	Operator SmallAngOpt (Vert_org);

	SmallAngOpt.constraints.isSmooth = isSmooth;
	SmallAngOpt.constraints.dev = devFactor;
	SmallAngOpt.constraints.isDelaunay = isDelaunay;

	SmallAngOpt.constraints.isMinAngle = true;
	SmallAngOpt.constraints.MinAngle = std::max(2.0, 1.5*myStats.GetMinAngle());

	SmallAngOpt.constraints.isMaxAngle = true;
	SmallAngOpt.constraints.MaxAngle = (maxAngleAllow<0.0) ? myStats.GetMaxAngle() : maxAngleAllow;

	SmallAngOpt.constraints.isNonobtuse = false;

	SmallAngOpt.constraints.isEdgeLen = false;
	SmallAngOpt.constraints.isMaximal = false;


	int verticesHandler[4];


	int ip2, ip3;
	int itter = 0;
	double vertex_void[3];

	if (verbose){
		fprintf(stdout, "\nRemoving tri-valent vertices --> ");
	}

	int numVert_before = numVert_imp;

	auto start_time = std::chrono::high_resolution_clock::now();
	SmallAngOpt.TriValentRemoval(numVert_imp, Vert_imp);
	auto end_time = std::chrono::high_resolution_clock::now();

	if (verbose){
		fprintf(stdout, "%d vertex removed\n", numVert_before - numVert_imp);
	}
	
	myStats.GetAllStats(numVert_imp, Vert_imp);

	if (verbose){
		DisplayStat(numVert_imp, Vert_imp,0);
		fprintf(stdout, "\n itter: %i", itter);
	}

	std::chrono::duration<double> elapsed_time_acc = end_time - start_time;
	isObtuse = false;

	while (true){
		itter++;

		auto start_time1 = std::chrono::high_resolution_clock::now();

		for (int i = 0; i < numVert_imp; i++){
			indices.push_back(i);
		}
		random_shuffle(indices.begin(), indices.end());

		// Relocation 		
		for (size_t id = 0; id < indices.size(); id++){

			int ip1 = indices[id];
			if (ip1 >= numVert_imp || Vert_imp[ip1].connect[0] == 0){ continue; }


			//test if it is still an obtuse head 
			if (TooSmallHeadAngle(ip1, ip2, ip3, targetMinAngle)){
				
				int closestSurfID = surfaceKdTree.FindNearest(Vert_imp[ip1].x);

				if (SmallAngOpt.constraints.dev >0){
					SmallAngOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - SmallAngOpt.constraints.dev;
				}

				bool fixed = false;

				//Relocation 			
				fixed = SmallAngOpt.Relocation(ip1, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer);
								

				//Ejection
				verticesHandler[0] = 2;

				if (!fixed){
					verticesHandler[1] = ip1;
					for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){
						verticesHandler[2] = Vert_imp[ip1].connect[i];
						vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[Vert_imp[ip1].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[Vert_imp[ip1].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[Vert_imp[ip1].connect[i]].x[2]) / 2.0;
												

						if (SmallAngOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
							fixed = true;
							break;
						}
					}
				}

				if (!fixed){
					verticesHandler[1] = ip2;
					for (int i = 1; i <= Vert_imp[ip2].connect[0]; i++){
						verticesHandler[2] = Vert_imp[ip2].connect[i];
						vertex_void[0] = (Vert_imp[ip2].x[0] + Vert_imp[Vert_imp[ip2].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip2].x[1] + Vert_imp[Vert_imp[ip2].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip2].x[2] + Vert_imp[Vert_imp[ip2].connect[i]].x[2]) / 2.0;

						if (SmallAngOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
							fixed = true;
							break;
						}
					}
				}

				if (!fixed){
					verticesHandler[1] = ip3;
					for (int i = 1; i <= Vert_imp[ip3].connect[0]; i++){
						verticesHandler[2] = Vert_imp[ip3].connect[i];
						vertex_void[0] = (Vert_imp[ip3].x[0] + Vert_imp[Vert_imp[ip3].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip3].x[1] + Vert_imp[Vert_imp[ip3].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip3].x[2] + Vert_imp[Vert_imp[ip3].connect[i]].x[2]) / 2.0;

						if (SmallAngOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
							fixed = true;
							break;
						}
					}
				}

				if (!fixed){
					//Injection
					verticesHandler[0] = 3;
					verticesHandler[1] = ip1;
					verticesHandler[2] = ip2;
					verticesHandler[3] = ip3;

					vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
					vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
					vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

					if (SmallAngOpt.AggressiveInjection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, 1, numSurfaceLayer, 0)){
						fixed = true;
					}
				}

				//Attracotr Ejection
				verticesHandler[0] = 2;
				if (!fixed){
					verticesHandler[1] = ip1;
					for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){
						verticesHandler[2] = Vert_imp[ip1].connect[i];
						vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
						vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
						vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

						if (SmallAngOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
							fixed = true;
							break;
						}
					}
				}

				if (!fixed){
					verticesHandler[1] = ip2;
					for (int i = 1; i <= Vert_imp[ip2].connect[0]; i++){

						verticesHandler[2] = Vert_imp[ip2].connect[i];
						vertex_void[0] = (Vert_imp[ip2].x[0] + Vert_imp[Vert_imp[ip2].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip2].x[1] + Vert_imp[Vert_imp[ip2].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip2].x[2] + Vert_imp[Vert_imp[ip2].connect[i]].x[2]) / 2.0;

						if (SmallAngOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
							fixed = true;
							break;
						}
					}
				}
				if (!fixed){
					verticesHandler[1] = ip3;
					for (int i = 1; i <= Vert_imp[ip3].connect[0]; i++){
						verticesHandler[2] = Vert_imp[ip3].connect[i];
						vertex_void[0] = (Vert_imp[ip3].x[0] + Vert_imp[Vert_imp[ip3].connect[i]].x[0]) / 2.0;
						vertex_void[1] = (Vert_imp[ip3].x[1] + Vert_imp[Vert_imp[ip3].connect[i]].x[1]) / 2.0;
						vertex_void[2] = (Vert_imp[ip3].x[2] + Vert_imp[Vert_imp[ip3].connect[i]].x[2]) / 2.0;

						if (SmallAngOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
							fixed = true;
							break;
						}
					}
				}

				//Repeller Injection 
				if (!fixed){
					verticesHandler[0] = 3;
					verticesHandler[1] = ip1;
					verticesHandler[2] = ip2;
					verticesHandler[3] = ip3;

					vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
					vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
					vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

					SmallAngOpt.AggressiveInjection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, 1, numSurfaceLayer, 1);
				}
			}
		}

		auto end_time1 = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> elapsed_time1 = end_time1 - start_time1;

		elapsed_time_acc = elapsed_time_acc + elapsed_time1;

		myStats.GetAllStats(numVert_imp, Vert_imp);

		if (verbose){
			DisplayStat(numVert_imp, Vert_imp,0);
			fprintf(stdout, "\n itter: %i", itter);			
		}

		if (myStats.GetMinAngle() > targetMinAngle){ break; }
	

		SmallAngOpt.constraints.MinAngle = std::max(2.0, 1.5*myStats.GetMinAngle());		
	}

	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "****************** Small Angle Elimination - Interleave Finished ****************");
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

	std::cout << "\n Total elapsed time: " << elapsed_time_acc.count() << " (s)\n";
		
	GetMeshOBJ("output.obj", 0, "");
}

//** Evaluate head **//
bool MeshImp::ObtuseHead(int ip, int&ip1, int&ip2)
{
	
	ip2 = Vert_imp[ip].connect[Vert_imp[ip].connect[0]];

	for (int i = 1; i <= Vert_imp[ip].connect[0]; i++){ 
		ip1 = Vert_imp[ip].connect[i];

		double angle = AngleVectVect(Vert_imp[ip1].x[0] - Vert_imp[ip].x[0], Vert_imp[ip1].x[1] - Vert_imp[ip].x[1], Vert_imp[ip1].x[2] - Vert_imp[ip].x[2],
			                         Vert_imp[ip2].x[0] - Vert_imp[ip].x[0], Vert_imp[ip2].x[1] - Vert_imp[ip].x[1], Vert_imp[ip2].x[2] - Vert_imp[ip].x[2])*RadToDeg; /* 57.295779513078550 = 180.0 / PI*/
		
		if (angle > 90.0 + _tol ){
			return true;
		}
		ip2 = ip1;
	}

	return false;

}
bool MeshImp::AcuteHead(int ip, int&ip1, int&ip2, double measureAngle)
{

	ip2 = Vert_imp[ip].connect[Vert_imp[ip].connect[0]];

	for (int i = 1; i <= Vert_imp[ip].connect[0]; i++){
		ip1 = Vert_imp[ip].connect[i];

		double angle = AngleVectVect(Vert_imp[ip1].x[0] - Vert_imp[ip].x[0], Vert_imp[ip1].x[1] - Vert_imp[ip].x[1], Vert_imp[ip1].x[2] - Vert_imp[ip].x[2],
			Vert_imp[ip2].x[0] - Vert_imp[ip].x[0], Vert_imp[ip2].x[1] - Vert_imp[ip].x[1], Vert_imp[ip2].x[2] - Vert_imp[ip].x[2])*RadToDeg;

		if (angle > measureAngle + _tol){
			return true;
		}
		ip2 = ip1;
	}

	return false;

}
bool MeshImp::TooSmallHeadAngle(int ip, int&ip1, int&ip2, double measureAngle)
{

	ip2 = Vert_imp[ip].connect[Vert_imp[ip].connect[0]];

	for (int i = 1; i <= Vert_imp[ip].connect[0]; i++){
		ip1 = Vert_imp[ip].connect[i];

		double angle = AngleVectVect(Vert_imp[ip1].x[0] - Vert_imp[ip].x[0], Vert_imp[ip1].x[1] - Vert_imp[ip].x[1], Vert_imp[ip1].x[2] - Vert_imp[ip].x[2],
			Vert_imp[ip2].x[0] - Vert_imp[ip].x[0], Vert_imp[ip2].x[1] - Vert_imp[ip].x[1], Vert_imp[ip2].x[2] - Vert_imp[ip].x[2])*RadToDeg;

		if (angle + _tol< measureAngle){
			return true;
		}
		ip2 = ip1;
	}

	return false;

}

//** Postprocessing **//
void MeshImp::WriteStatToFile(std::string filename, int numV, vert*Vert, bool obt)
{
	Statistics myStats;
	myStats.GetAllStats(numV, Vert);

	FILE *fp;
	errno_t err;
	err = fopen_s(&fp, filename.c_str(), "w");
	myStats.DisplayStatistics(fp, obt);
	fclose(fp);

}
void MeshImp::DisplayStat(int numV, vert*Vert, bool obt)
{
	Statistics myStats;
	myStats.GetAllStats(numV, Vert);
	myStats.DisplayStatistics(stdout, obt);

}
void MeshImp::GetMeshOBJ(std::string filename, bool obtuse, std::string filename_obt)
{
	//write mesh in .obj file 
	fprintf(stdout, "\n Writing output file.\n");

	std::fstream file(filename.c_str(), std::ios::out);	
	std::fstream file_obt(filename_obt.c_str(), std::ios::out);

	

	file << "#V " <<numVert_imp << std::endl;
	for (int id0 = 0; id0 < numVert_imp; id0++){
		file << "v " << Vert_imp[id0].x[0] << " " << Vert_imp[id0].x[1] << " " << Vert_imp[id0].x[2] << std::endl;
		if (obtuse){
			file_obt << "v " << Vert_imp[id0].x[0] << " " << Vert_imp[id0].x[1] << " " << Vert_imp[id0].x[2] << std::endl;
		}
	}

	for (int id0 = 0; id0 < numVert_imp; id0++){
		int id1 = Vert_imp[id0].connect[Vert_imp[id0].connect[0]];

		for (int i = 1; i <= Vert_imp[id0].connect[0]; i++){
			int id2 = Vert_imp[id0].connect[i];

			if (id0 < id1 && id0 < id2){
				file << "f " << id0 + 1 << "  " << id1 + 1 << "  " << id2 + 1 << std::endl;
				if (obtuse){
					double angle1 = AngleVectVect(Vert_imp[id1].x[0] - Vert_imp[id0].x[0], Vert_imp[id1].x[1] - Vert_imp[id0].x[1], Vert_imp[id1].x[2] - Vert_imp[id0].x[2],
						                          Vert_imp[id2].x[0] - Vert_imp[id0].x[0], Vert_imp[id2].x[1] - Vert_imp[id0].x[1], Vert_imp[id2].x[2] - Vert_imp[id0].x[2])*RadToDeg;

					double angle2 = AngleVectVect(Vert_imp[id0].x[0] - Vert_imp[id1].x[0], Vert_imp[id0].x[1] - Vert_imp[id1].x[1], Vert_imp[id0].x[2] - Vert_imp[id1].x[2],
						                          Vert_imp[id2].x[0] - Vert_imp[id1].x[0], Vert_imp[id2].x[1] - Vert_imp[id1].x[1], Vert_imp[id2].x[2] - Vert_imp[id1].x[2])*RadToDeg;

					double angle3 = AngleVectVect(Vert_imp[id1].x[0] - Vert_imp[id2].x[0], Vert_imp[id1].x[1] - Vert_imp[id2].x[1], Vert_imp[id1].x[2] - Vert_imp[id2].x[2],
						                          Vert_imp[id0].x[0] - Vert_imp[id2].x[0], Vert_imp[id0].x[1] - Vert_imp[id2].x[1], Vert_imp[id0].x[2] - Vert_imp[id2].x[2])*RadToDeg;
					if (angle2>90.0 + _tol || angle3>90.0 + _tol){
						file_obt << "f " << id0 + 1 << "  " << id1 + 1 << "  " << id2 + 1 << std::endl;
					}
				}
			}
			id1 = id2;
		}
	}

	file.close();
}

void PointOnSphere(double&x1, double&y1, double&z1, double x_s, double y_s, double z_s, double r_s) // Setting point(x1,y1,z1) to be on the sphere of radius r_s centered at (x_s, y_s, z_s)
{
	double x_v(x1 - x_s), y_v(y1 - y_s), z_v(z1 - z_s);

	double n = sqrt(x_v*x_v + y_v*y_v + z_v*z_v);

	x1 = x_v*r_s / n + x_s;
	y1 = y_v*r_s / n + y_s;
	z1 = z_v*r_s / n + z_s;


}
void MeshImp::PerturbVertices()
{
	//input is vertices one a sphere 
	//move verices that are not head of obtuse angles
	//in random direction 10% of shorest edge its connected to
	
	//srand(time(NULL));
	int n1, n2;
	RndNum myRandNum;
	myRandNum.InitiateRndNum(rand());

	for (int ip = 0; ip < numVert_imp; ip++){
		if (!ObtuseHead(ip, n1, n2)){
			double shortest_edge_len = 10E10;
			for (int j = 1; j <= Vert_imp[ip].connect[0]; j++){
				int iq = Vert_imp[ip].connect[j];
				shortest_edge_len = std::min(shortest_edge_len, Dist(Vert_imp[ip].x[0], Vert_imp[ip].x[1], Vert_imp[ip].x[2],
					                                                 Vert_imp[iq].x[0], Vert_imp[iq].x[1], Vert_imp[iq].x[2]));
			}
			shortest_edge_len = sqrt(shortest_edge_len);
			shortest_edge_len /= 4.0;

			Vert_imp[ip].x[0] += (myRandNum.RandNumGenerator())* shortest_edge_len;
			Vert_imp[ip].x[1] += (myRandNum.RandNumGenerator())* shortest_edge_len;
			Vert_imp[ip].x[2] += (myRandNum.RandNumGenerator())* shortest_edge_len;

			PointOnSphere(Vert_imp[ip].x[0], Vert_imp[ip].x[1], Vert_imp[ip].x[2], 0.5, 0.5, 0.5, 0.5*sqrt(3.0));

		}
	}

	GetMeshOBJ("sphere8Perturbed.obj", 1, "input_obtuse.obj");
}