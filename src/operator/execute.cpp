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
#include <array>

bool isObtuse;
bool isAcute;
#define MAXDEL 45
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
MeshImp::MeshImp(int numSpheres, double **Spheres, int numVert, double**Verts, int numTri, int**Tris)
{	
	fprintf(stdout, "\nMeshImp:: Initialize\n");
	//1) store  tessellation  and find neighbour triangles 
	mTriangles.ids.reserve(numTri);
	mTriangles.neighbour.reserve(numTri);
	mTriangles.coord.reserve(numVert);

	for (int i = 0; i < numVert; i++){
		std::array<double, 3> x;
		x[0] = Verts[i][0];
		x[1] = Verts[i][1];
		x[2] = Verts[i][2];
		mTriangles.coord.push_back(x);
	}

	for (int i = 0; i < numTri; i++){
		std::array<int, 3> id;
		id[0] = Tris[i][0];
		id[1] = Tris[i][1];
		id[2] = Tris[i][2];
		mTriangles.ids.push_back(id);
	}		
	FindNeighbourTriangles();

	
	//2)initialize the spheres and its kd tree 
	mSphere.init(numSpheres, Spheres, mTriangles);


}
void MeshImp::FindNeighbourTriangles()
{
	int** myTriangles;
	int numTriVertices = mTriangles.coord.size();
	int numTriangles = mTriangles.ids.size();
	myTriangles = new int*[numTriVertices]; //for each vertex i, temporary store the triangles where i is a head on
	int common[20];
	for (int i = 0; i < numTriVertices; i++){
		myTriangles[i] = new int[MAXDEL];
		myTriangles[i][0] = 0;
	}

	//loop over triangles and populate myTriangles 
	for (int t = 0; t < numTriangles; t++){
		int id0(mTriangles.ids[t][0]), id1(mTriangles.ids[t][1]), id2(mTriangles.ids[t][2]);
		myTriangles[id0][++myTriangles[id0][0]] = t;
		myTriangles[id1][++myTriangles[id1][0]] = t;
		myTriangles[id2][++myTriangles[id2][0]] = t;

		if (myTriangles[id0][0] >= MAXDEL || myTriangles[id1][0] >= MAXDEL || myTriangles[id2][0] >= MAXDEL){
			ErrWarnMessage(__LINE__, "MeshImp::FindNeighbourTriangles:: a vertex is shared within more than MAXDEL triangles!!", 0);
		}
	}

	for (int t = 0; t < numTriangles; t++){
		std::array<int, 3> neigh;
		//for each triangle t 
		//find the other triangle shared between two of its heads
		//update the t's neighbour list 

		for (int i = 0; i < 3; i++){
			int j = (i == 2) ? 0 : i + 1;

			if (FindCommonElements(myTriangles[mTriangles.ids[t][i]], myTriangles[mTriangles.ids[t][j]], common)){
				if (common[0] != 2){
					ErrWarnMessage(__LINE__, "MeshImp::FindNeighbourTriangles:: Non-manifold surface", 0);
				}
				else{
					neigh[i] = (common[1] == t) ? common[2] : common[1];
				}
			}
			else{
				ErrWarnMessage(__LINE__, "MeshImp::FindNeighbourTriangles::Non-manifold surface", 0);
			}
		}

		mTriangles.neighbour.push_back(neigh);
	}


	//clean up
	for (int i = 0; i < numTriVertices; i++){ delete[] myTriangles[i]; }
	delete myTriangles;
	
}
//** Sifting/Simplification **//
void MeshImp::VC(int targetNumSamples, int samplingBudget, int numSurfaceLayer, bool isSmooth, bool verbose)
{			
	fprintf(stdout, "\nMeshImp:: Spheres Simplification\n");
	Operator SimpOpt;
	int verticesHandler[4];
	verticesHandler[0] = 2;
	
	int itter = 0;
	double vertex_void[3];
			

	auto start_time = std::chrono::high_resolution_clock::now();	
	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_time_acc;
	int numSpheres_before;
	const int numSpheres_input = mSphere.getNumSphere();
		
	while (true){
		numSpheres_before = mSphere.getNumSphere();
		itter++;

		auto start_time1 = std::chrono::high_resolution_clock::now();
				
		//******* Ejection Two
		verticesHandler[0] = 2;		

		for (int id = 0; id < mSphere.getNumSphere(); id++){
			if (mSphere.getTag(id) != 0){ continue; }
												
			bool ejected = false;
			std::vector<int> overlap;
			mSphere.getOverlap(id, overlap);
			verticesHandler[1] = id;

			double x_id, y_id, z_id, r_id2;
			int tri_id = 0;
			if (!mSphere.getCoordinates(id, x_id, y_id, z_id, r_id2, tri_id)){ continue; }
						
			for (int i = 0; i < int(overlap.size()); i++){

				verticesHandler[2] = overlap[i];

				if (mSphere.getTag(overlap[i]) != 0){ continue; }

				double x_i, y_i, z_i, r_i2;
				int tri_i = 0;


				if (!mSphere.getCoordinates(overlap[i], x_i, y_i, z_i, r_i2,tri_i)){ continue; }
				vertex_void[0] = (x_i + x_id)/ 2.0;
				vertex_void[1] = (y_i + y_id) / 2.0;
				vertex_void[2] = (z_i + z_id) / 2.0;
				

				if (false){
					mSphere.Draw_m_sphere(id);
					mSphere.Draw_m_sphere(overlap[i]);
					mSphere.Draw_m_sphere_seeds(id);
				}

				if (SimpOpt.EjectionSpheres(verticesHandler, &mTriangles, &mSphere, samplingBudget, numSurfaceLayer, isSmooth, verbose)){
					ejected = true;					
					break;
				}
				
			}			
			if (mSphere.getNumSphere() <= targetNumSamples){ break; }			
		}
		
		auto end_time1 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed_time1 = end_time1 - start_time1;				
		elapsed_time_acc = elapsed_time_acc + elapsed_time1;		
				
		
		if (verbose){
			fprintf(stdout, "\n Reducation Ratio= %f %", 100.0*(double(numSpheres_input - mSphere.getNumSphere()) / double(numSpheres_input)));			
			fprintf(stdout, "\n itter: %i", itter);		
			mSphere.vcSurface();
		}			

		if (mSphere.getNumSphere() <= targetNumSamples){ break; }
	}
	
	if (verbose){		
		std::cout << " \nTotal elapsed time: " << elapsed_time_acc.count() << " (s)\n";
	}
		

}