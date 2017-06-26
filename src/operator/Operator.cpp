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

#include <stdlib.h>
#include <algorithm>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits.h>
#include <cfloat>

#include "Operator.h"
#include "../util/Common.h"
#include "../util/RNG.h"

//#include "DrawSpheresDebug.h" //for debugging 

//#define DEBUGGING
extern bool isObtuse;
extern bool isAcute;

double OpimizerFunc_CenterAngle(vert*Verts, int*nList, double*void_verx, double xx, double yy, double zz);
double OpimizerFunc_SideAngle(vert*Verts, int*nList, double minAngle, double maxAngle, double xx, double yy, double zz);
double OpimizerFunc_Closer(vert*Verts, int*nList, double*void_verx, double xx, double yy, double zz);
double OpimizerFunc_Further(vert*Verts, int*nList, double*void_verx, double xx, double yy, double zz);

Operator::Operator(vert*Vert_org):
                    Vert_org(Vert_org)
{
	active_pool = new int *[MAXPOOL];
	tri_pool = new int *[MAXPOOL];
	tmp_active_pool = new int*[MAXPOOL];
	for (int i = 0; i < MAXPOOL; i++){
		active_pool[i] = new int[3];
		tri_pool[i] = new int[3];
		tmp_active_pool[i] = new int[3];
	}

	_tar = new double*[4];
	_tar[0] = new double[9];
	_tar[1] = new double[9];
	_tar[2] = new double[9];
	_tar[3] = new double[9];
	

	srand(time(NULL));	
	myRandNum.InitiateRndNum(rand());

}

Operator::~Operator()
{

}

///**************Special
void Operator::TriValentRemoval(int&numVert, vert*Verts)
{
	//remove all trivalent regardless to their effect on the quality
	//do it multiple times
	//everything can be fixed later
	//better to call this before applying the operators 
	//some operators (injection) quit when it finds tri-valents on a patch since it will mess up the toplogical correctness 
	//other operators can tolerate tri-valents vertices as long as topological correctness holds 
	while (true){
		bool again(false);
		for (int i = numVert; i >= 0; i--){
			if (Verts[i].connect[0] == 3){
				RemoveVertex(i, numVert, Verts);
				again = true;
			}
		}
		if (!again){ break; }
	}
}



///**************Operators
bool Operator::Relocation(int ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer)
{
	//TODO pass the optimizer function also
	//ip = the vertex to relocate 
	//numVert =  total num of vertices in mesh
	//Vert =  the mesh to update



	myConstraints.Reset(Verts, NULL);
		
	skipList[0] = 1;
	skipList[1] = ip;
	
	if (constraints.isDelaunay){
		//Delaunay constraints 
		myConstraints.Delaunay(Verts[ip].connect, skipList, Verts, true,1);
		if (!myConstraints.OverlappingInSpheres()){ 
			//empty inclusion region
			return false;
		}
	}

	myConstraints.Direction(skipList, Verts[ip].connect, Verts);

	if (constraints.isMinAngle || constraints.isMaxAngle){
		//Min max angle 
		myConstraints.MinMaxAngle(Verts[ip].connect, Verts, constraints.MinAngle, constraints.MaxAngle, Verts[ip].x, constraints.isNonobtuse,1);
	}

	if (constraints.isEdgeLen){
		//Min edge length 
		myConstraints.MinEdgeLength(Verts[ip].connect, Verts, constraints.MinEdgeLength_sq);
	}

	if (constraints.isSmooth){
		//Smoothness 
		myConstraints.Smoothness(Verts[ip].connect, skipList, Verts, Verts[ip].x, constraints.dev);
	}

	
	StartActivePool(closedtSurfaceID, numSurfaceLayer);
	if (Sampler(Verts, Verts[ip].connect, samplingBudget, -1, -1, -1, &OpimizerFunc_CenterAngle)){
		Mover(ip, Verts);
		return true;
	}
	return false;

}
bool Operator::Ejection(int* ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer, bool att)
{
	//TODO pass the optimizer function also
	//ip = pointer to vertices to eject
	//numVert =  total num of vertices in mesh
	//Vert =  the mesh to update
	
#ifdef DEBUGGING
	if (ip[0] != 2 && ip[0] != 3){
		ErrWarnMessage(__LINE__, "Operator::Ejection:: can only ejection two or three vertices", 0);
	}
#endif

	//sanity check 
	//make sure vertices in ip are not connected to the same tri valent node 
	
    skipList[0] = 0;//use this list to store the tri valent nodes connected to any of the to-be ejected vertices 
	TriValent(ip[1], Verts, skipList, ip[2], (ip[0] == 3) ? ip[3] : INT_MAX);
	TriValent(ip[2], Verts, skipList, ip[1], (ip[0] == 3) ? ip[3] : INT_MAX);
	if (ip[0] == 3){
		TriValent(ip[3], Verts, skipList, ip[1], ip[2]);
	}
	if (FindDuplication(skipList)){ return false; }
	
	//get mid vertex and copy skipList
	skipList[0] = ip[0];
	void_vertex[0] = void_vertex[1] = void_vertex[2] = 0;
	for (int i = 1; i <= ip[0]; i++){
		skipList[i] = ip[i];
		void_vertex[0] += Verts[ip[i]].x[0];
		void_vertex[1] += Verts[ip[i]].x[1];
		void_vertex[2] += Verts[ip[i]].x[2];
	}
	void_vertex[0] /= double(ip[0]); void_vertex[1] /= double(ip[0]); void_vertex[2] /= double(ip[0]);


	//get the sort neighbout list (unduplicated list of vertices connected to *ip and sorted)	
	mynList[0] = 0;
	if (ip[0] == 2){
		if (Verts[ip[1]].connect[0] == 3){
			GetListSkipVertex(ip[2], Verts, ip[1],INT_MAX, mynList);
		}
		else if (Verts[ip[2]].connect[0] == 3){
			GetListSkipVertex(ip[1], Verts, ip[2], INT_MAX, mynList);
		}
		else{
			GetEdgeSortedNeighbourList(ip[1], ip[2], Verts, mynList);
			if (mynList[0] < 0 || FindDuplication(mynList)){
				return false; 
			}	
		}
	}
	else{
		//TODO get the list correct when ejecting three vertices 
		//take care of special cases of tri-valent nodes 

		//it is not possible to have two tri-valent node conneced in a watertight manifold 
		if (Verts[ip[1]].connect[0] == 3 || Verts[ip[2]].connect[0] == 3 || Verts[ip[3]].connect[0] == 3){
			return false;
		}
		else {
			GetFaceSortedNeighbourList(ip[1], ip[2], ip[3], Verts, mynList);
			if (mynList[0] < 0 || FindDuplication(mynList)){
				return false; 
			}
		}
	}
	if (mynList[0] == 0){ 
		//could not get the neighbours right 
		return false;
	}

	///******* Attractor 
	if (att){
		SetAttractorRepeller(Verts);
		AttractorRepeller(ip, numVert, Verts, closedtSurfaceID, samplingBudget, numSurfaceLayer, true);
	}
	
	myConstraints.Reset(Verts, NULL);
	myConstraints.Direction(ip, mynList, Verts);

	if (constraints.isDelaunay){
		//Delaunay constraints 
		myConstraints.Delaunay(mynList, skipList, Verts, true, 1);
		if (!myConstraints.OverlappingInSpheres()){
			//empty inclusion region
			if (att){
				RevertAttractorRepeller(Verts);				
			}
			return false;
		}
	}

	if (constraints.isMinAngle || constraints.isMaxAngle){
		//Min max angle 
		myConstraints.MinMaxAngle(mynList, Verts, constraints.MinAngle, constraints.MaxAngle, void_vertex, constraints.isNonobtuse,1);
	}

	if (constraints.isEdgeLen){
		//Min edge length 
		myConstraints.MinEdgeLength(mynList, Verts, constraints.MinEdgeLength_sq);
	}
	if (constraints.isSmooth){
		//Smoothness 
		myConstraints.Smoothness(mynList, ip, Verts, void_vertex, constraints.dev);
	}
		
	StartActivePool(closedtSurfaceID, numSurfaceLayer);
	if (Sampler(Verts, mynList, samplingBudget, -1, -1, -1, &OpimizerFunc_CenterAngle)){
		//for update, call EdgeCollapse if it is two vertices or
		//FaceCollapse if it is three vertices 
		if (ip[0] == 2){
			EdgeCollapse(ip[1],ip[2],numVert,Verts);
		}
		else {		
			FaceCollapse(ip[1], ip[2], ip[3], numVert, Verts);
		}
		return true;
	}
		
	if (att){
		RevertAttractorRepeller(Verts);
	}

	return false;
}
bool Operator::AggressiveInjection(int*ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer, bool inj)
{
	//TODO pass the optimizer function also 
	//TODO check that three vertices in ip forms a triangle

	//ip = pointer to vertices to eject (should be three vertices)
	//numVert =  total num of vertices in mesh
	//Vert =  the mesh to update
	
#ifdef DEBUGGING
	if (ip[0] != 3){
		ErrWarnMessage(__LINE__, "Operator::Injection:: Invalud input. Correct input is three vertices (triangle)", 0);
	}
#endif

	if (Verts[ip[1]].connect[0] == 3 || Verts[ip[2]].connect[0] == 3 || Verts[ip[3]].connect[0] == 3){
		//nop, we don't do this 
		return false;
	}
		

	skipList[0] = ip[0];
	void_vertex[0] = void_vertex[1] = void_vertex[2] = 0.0;
	for (int i = 1; i <= ip[0]; i++){
		skipList[i] = ip[i];
		void_vertex[0] += Verts[ip[i]].x[0];
		void_vertex[1] += Verts[ip[i]].x[1];
		void_vertex[2] += Verts[ip[i]].x[2];
	}
	void_vertex[0] /= double(ip[0]);
	void_vertex[1] /= double(ip[0]);
	void_vertex[2] /= double(ip[0]);

	//find the correct shared vertex between each two consecutive vertices in ip
	apex[0] = ip[0];
	int i_skip = -1;
	for (int i = 1; i <= ip[0]; i++){
		int iq = (i == ip[0]) ? ip[1] : ip[i + 1];
		//apex[i] = FindCommonElement_SkipList(Verts[ip[i]].connect, Verts[iq].connect, ip);

		if (FindCommonElements_SkipList(Verts[ip[i]].connect, Verts[iq].connect, ip, temp_arr)){
			if (temp_arr[0] == 1){
				apex[i] = temp_arr[1];
			}
			else if (temp_arr[0] == 2){

				if (Verts[temp_arr[1]].connect[0] == 3 || Verts[temp_arr[2]].connect[0] == 3){
					//get the trivalen one 
					apex[i] = (Verts[temp_arr[1]].connect[0] == 3) ? temp_arr[1] : temp_arr[2];
				}
				else{

					//in this case we seek the vertex that is not connected to the other apex
					//since other apex's are not discovered yet, we skip this one and do it after getting
					//other apexs. There should not be a way that there is more than one apex that is such problomatic 
					//but we check on this anyways 
#ifdef DEBUGGING
					if (i_skip >= 0){
						ErrWarnMessage(__LINE__, "Operator::Injection:: not considered. See comment", 0);
					}
#endif
					i_skip = i;
				}
			}
			else{
#ifdef DEBUGGING
				ErrWarnMessage(__LINE__, "Operator::Injection:: 2 can not get apex", 0);
#endif
				return false;
			}

		}
		else{
			return false;
#ifdef DEBUGGING
			ErrWarnMessage(__LINE__, "Operator::Injection:: 1 can not get apex", 0);
#endif
			return false;
		}

	}

	if (i_skip >= 0){
		//we have on problomatic apex 
		//at least this apex should not be connected to 
		int iq = (i_skip == ip[0]) ? ip[1] : ip[i_skip + 1];

		FindCommonElements_SkipList(Verts[ip[i_skip]].connect, Verts[iq].connect, ip, temp_arr);

		bool one_shared(false), two_shared(false);
		for (int j = 1; j <= 3; j++){
			if (j == i_skip){ continue; }
			if (GetIndex(temp_arr[1], Verts[apex[j]].connect) >= 0){ one_shared = true; }
			if (GetIndex(temp_arr[2], Verts[apex[j]].connect) >= 0){ two_shared = true; }
		}

#ifdef DEBUGGING
		if ((one_shared&&two_shared) || (!one_shared&&!two_shared)){
			//means both of the temp_arr[1] and temp_arr[2] (candidate apex)
			//are either shared with other apex's and not shared at all
			//have not considered this case yet
			ErrWarnMessage(__LINE__, "Operator::Injection:: not considered. See comment", 0);
		}
#endif

		apex[i_skip] = (two_shared) ? temp_arr[1] : temp_arr[2];
	}
	
	if (Verts[apex[1]].connect[0] == 3 || Verts[apex[2]].connect[0] == 3 || Verts[apex[3]].connect[0] == 3){
		//nop, we don't do this either 
		return false;
	}

	


	//populate the void corner list	
	mynList[0] = 6;
	mynList[1] = ip[1];
	mynList[2] = apex[1];
	mynList[3] = ip[2];
	mynList[4] = apex[2];
	mynList[5] = ip[3];
	mynList[6] = apex[3];


	//*******1) Destory the whole triangle as if we refining it 
	//and removing the edges ip[1]-ip[2], ip[2]-ip[3] and ip[3]-ip[1]
	//thus creating a void surrounded by six vertices 
	//only if non of the three vertices in ip are 4-valent 
	//because when we insert a node in the middle in this fashion,
	//we effectively, reduce the valence of each node in ip by one 

	if (Verts[ip[1]].connect[0] > 4 && Verts[ip[2]].connect[0] > 4 && Verts[ip[3]].connect[0] > 4 /*&& InspectFeasibleRegion(mynList, Verts)*/){

		if (inj){
			SetAttractorRepeller(Verts);
			AttractorRepeller(ip, numVert, Verts, closedtSurfaceID, samplingBudget, numSurfaceLayer, false); 
		}
			
		myConstraints.Reset(Verts, NULL); 
		
		myConstraints.Direction(NULL, mynList, Verts);

		if (constraints.isDelaunay){
			//Delaunay constraints 
			myConstraints.Delaunay(mynList, skipList, Verts, false, 1);
		}

		if (myConstraints.OverlappingInSpheres()){
			if (constraints.isMinAngle || constraints.isMaxAngle){
				//Min max angle 
				myConstraints.MinMaxAngle(mynList, Verts, constraints.MinAngle, constraints.MaxAngle, void_vertex, constraints.isNonobtuse, 1);
			}

			if (constraints.isEdgeLen){
				//Min edge length 
				myConstraints.MinEdgeLength(mynList, Verts, constraints.MinEdgeLength_sq);
			}
			if (constraints.isSmooth){
				//Smoothness 
				myConstraints.Smoothness(mynList, ip, Verts, void_vertex, constraints.dev);
			}
		
			StartActivePool(closedtSurfaceID, numSurfaceLayer);
			if (Sampler(Verts, mynList, samplingBudget, -1, -1, -1, &OpimizerFunc_CenterAngle)){
				//load the six edges to break 
				//ip[1]-ip[2]
				//ip[1]-ip[3]
				//ip[2]-ip[3]
				EdgeToBreak[0] = 3;
				EdgeToBreak[1] = ip[1];
				EdgeToBreak[2] = ip[2];
				EdgeToBreak[3] = ip[1];
				EdgeToBreak[4] = ip[3];
				EdgeToBreak[5] = ip[2];
				EdgeToBreak[6] = ip[3];
				Inserter(numVert, Verts);
				return true;
			}
		}
		if (inj){
			RevertAttractorRepeller(Verts);
		}
	}


	//*******2) Destory one edge of the triangle and another one 
	//there is 9 possibilities in near, we try them all :D

	int ip3 = ip[ip[0]];
	for (int i = 1; i <= ip[0]; i++){
		int j = (i == ip[0]) ? 1 : i + 1;

		int ip1 = ip[i];
		int ip2 = ip[j];
		int iapex = apex[i];
		int japex = apex[j];

		//we take the void_vertex as the avergae of two points 
		//point one is the mid-point between ip1-ip2
		//point one is the mid-point between ip2-ip3
		void_vertex[0] = (Verts[ip1].x[0] + 2.0*Verts[ip2].x[0] + Verts[ip3].x[0]) / 4.0;
		void_vertex[1] = (Verts[ip1].x[1] + 2.0*Verts[ip2].x[1] + Verts[ip3].x[1]) / 4.0;
		void_vertex[2] = (Verts[ip1].x[2] + 2.0*Verts[ip2].x[2] + Verts[ip3].x[2]) / 4.0;

		
		//#1 first try to destroy the edges ip1-ip2 ip2-ip3
		if (Verts[ip2].connect[0] > 4){ //this prevents ip2 from being tri-valent
			//(because we remove two verices connected to it and add one, effectively reduce valence by 1)		

			mynList[0] = 5;
			mynList[1] = ip1;
			mynList[2] = iapex;
			mynList[3] = ip2;
			mynList[4] = japex;
			mynList[5] = ip3;
					
			skipList[0] = 3;

			if (inj){
				SetAttractorRepeller(Verts);
				AttractorRepeller(ip, numVert, Verts, closedtSurfaceID, samplingBudget, numSurfaceLayer, false);
			}

			myConstraints.Reset(Verts, NULL);
			myConstraints.Direction(NULL, mynList, Verts);

			if (constraints.isDelaunay){
				//Delaunay constraints 
				myConstraints.Delaunay(mynList, skipList, Verts, false, 1);
			}

			if (myConstraints.OverlappingInSpheres()){
				if (constraints.isMinAngle || constraints.isMaxAngle){
					//Min max angle 
					myConstraints.MinMaxAngle(mynList, Verts, constraints.MinAngle, constraints.MaxAngle, void_vertex, constraints.isNonobtuse, 1);
				}

				if (constraints.isEdgeLen){
					//Min edge length 
					myConstraints.MinEdgeLength(mynList, Verts, constraints.MinEdgeLength_sq);
				}

				if (constraints.isSmooth){
					//Smoothness 
					myConstraints.Smoothness(mynList, skipList, Verts, void_vertex, constraints.dev);
				}


				StartActivePool(closedtSurfaceID, numSurfaceLayer);
				if (Sampler(Verts, mynList, samplingBudget, -1, -1, -1, &OpimizerFunc_CenterAngle)){
					//load edge to break
					//ip1-ip2
					//ip2-ip3
					EdgeToBreak[0] = 2;
					EdgeToBreak[1] = ip1;
					EdgeToBreak[2] = ip2;
					EdgeToBreak[3] = ip2;
					EdgeToBreak[4] = ip3;
					Inserter(numVert, Verts);
					return true;
				}
			}

			if (inj){
				RevertAttractorRepeller(Verts);
			}
		}



		//#2 second try to destroy the edges ip1-ip2 with another edge connected to on of the apex 

		for (int id = 1; id <= 2; id++){
			if (id == 2){
				std::swap(ip1, ip2);//trust me you need this.
			}

			if (Verts[ip1].connect[0] > 4){//this prevents ip1 from being tri-valent
				//(because we remove to verices connected to it and add one, effectively reduce valence by 1)

				//the two edges to destroy are ip1-ip2  and iapex1-ip1
				int extra_apex = FindCommonElement_SkipList(Verts[ip1].connect, Verts[iapex].connect, ip);

				if (Verts[extra_apex].connect[0] > 3){//if it is 3-valent, then the void is ambigious 
					//the only solution would be creating 4-valent vertex which is equally bad for most cases
					mynList[0] = 5;
					mynList[1] = ip1;
					mynList[2] = extra_apex;
					mynList[3] = iapex;
					mynList[4] = ip2;
					mynList[5] = ip3;
					
					//here we take the void_vertex as the avergae of two points 
					//point one is the mid-point between ip1-iapex
					//point one is the mid-point between ip1-ip2
					void_vertex[0] = (2.0*Verts[ip1].x[0] + Verts[iapex].x[0] + Verts[ip2].x[0]) / 4.0;
					void_vertex[1] = (2.0*Verts[ip1].x[1] + Verts[iapex].x[1] + Verts[ip2].x[1]) / 4.0;
					void_vertex[2] = (2.0*Verts[ip1].x[2] + Verts[iapex].x[2] + Verts[ip2].x[2]) / 4.0;

					skipList[0] = 4;
					skipList[skipList[0]] = iapex;

					if (inj){
						SetAttractorRepeller(Verts);
						AttractorRepeller(ip, numVert, Verts, closedtSurfaceID, samplingBudget, numSurfaceLayer, false);
					}

					myConstraints.Reset(Verts, NULL);
					myConstraints.Direction(NULL, mynList, Verts);

					if (constraints.isDelaunay){
						//Delaunay constraints 
						myConstraints.Delaunay(mynList, skipList, Verts, false, 1);
					}

					if (myConstraints.OverlappingInSpheres()){
						if (constraints.isMinAngle || constraints.isMaxAngle){
							//Min max angle 
							myConstraints.MinMaxAngle(mynList, Verts, constraints.MinAngle, constraints.MaxAngle, void_vertex, constraints.isNonobtuse, 1);
						}

						if (constraints.isEdgeLen){
							//Min edge length 
							myConstraints.MinEdgeLength(mynList, Verts, constraints.MinEdgeLength_sq);
						}

						if (constraints.isSmooth){
							//Smoothness 
							myConstraints.Smoothness(mynList, skipList, Verts, void_vertex, constraints.dev);
						}


						StartActivePool(closedtSurfaceID, numSurfaceLayer);
						if (Sampler(Verts, mynList, samplingBudget, -1, -1, -1, &OpimizerFunc_CenterAngle)){
							EdgeToBreak[0] = 2;
							EdgeToBreak[1] = ip1;
							EdgeToBreak[2] = ip2;
							EdgeToBreak[3] = ip1;
							EdgeToBreak[4] = iapex;
							Inserter(numVert, Verts);
							return true;
						}
					}

					if (inj){
						RevertAttractorRepeller(Verts);
					}
				}
			}
		}

		ip3 = ip2;//not ip1 because we swap it 

	}



	return false;
}
bool Operator::Injection(int*ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer, bool inj)
{
	//TODO pass the optimizer function also 
	//TODO check that three vertices in ip forms a triangle

	//ip = pointer to vertices to eject (should be three vertices)
	//numVert =  total num of vertices in mesh
	//Vert =  the mesh to update

#ifdef DEBUGGING
	if (ip[0] != 3){
		ErrWarnMessage(__LINE__, "Operator::Injection:: Invalud input. Correct input is three vertices (triangle)", 0);
	}
#endif

	if (Verts[ip[1]].connect[0] == 3 || Verts[ip[2]].connect[0] == 3 || Verts[ip[3]].connect[0] == 3){
		//nop, we don't do this 
		return false;
	}


	skipList[0] = ip[0];
	void_vertex[0] = void_vertex[1] = void_vertex[2] = 0.0;
	for (int i = 1; i <= ip[0]; i++){
		skipList[i] = ip[i];
		void_vertex[0] += Verts[ip[i]].x[0];
		void_vertex[1] += Verts[ip[i]].x[1];
		void_vertex[2] += Verts[ip[i]].x[2];
	}
	void_vertex[0] /= double(ip[0]);
	void_vertex[1] /= double(ip[0]);
	void_vertex[2] /= double(ip[0]);

	//find the correct shared vertex between each two consecutive vertices in ip
	apex[0] = ip[0];
	int i_skip = -1;
	for (int i = 1; i <= ip[0]; i++){
		int iq = (i == ip[0]) ? ip[1] : ip[i + 1];
		//apex[i] = FindCommonElement_SkipList(Verts[ip[i]].connect, Verts[iq].connect, ip);

		if (FindCommonElements_SkipList(Verts[ip[i]].connect, Verts[iq].connect, ip, temp_arr)){
			if (temp_arr[0] == 1){
				apex[i] = temp_arr[1];
			}
			else if (temp_arr[0] == 2){

				if (Verts[temp_arr[1]].connect[0] == 3 || Verts[temp_arr[2]].connect[0] == 3){
					//get the trivalen one 
					apex[i] = (Verts[temp_arr[1]].connect[0] == 3) ? temp_arr[1] : temp_arr[2];
				}
				else{

					//in this case we seek the vertex that is not connected to the other apex
					//since other apex's are not discovered yet, we skip this one and do it after getting
					//other apexs. There should not be a way that there is more than one apex that is such problomatic 
					//but we check on this anyways 
#ifdef DEBUGGING
					if (i_skip >= 0){
						ErrWarnMessage(__LINE__, "Operator::Injection:: not considered. See comment", 0);
					}
#endif
					i_skip = i;
				}
			}
			else{
				return false;
				ErrWarnMessage(__LINE__, "Operator::Injection:: 2 can not get apex", 0);
			}

		}
		else{
			return false;
			ErrWarnMessage(__LINE__, "Operator::Injection:: 1 can not get apex", 0);
		}

	}

	if (i_skip >= 0){
		//we have on problomatic apex 
		//at least this apex should not be connected to 
		int iq = (i_skip == ip[0]) ? ip[1] : ip[i_skip + 1];

		FindCommonElements_SkipList(Verts[ip[i_skip]].connect, Verts[iq].connect, ip, temp_arr);

		bool one_shared(false), two_shared(false);
		for (int j = 1; j <= 3; j++){
			if (j == i_skip){ continue; }
			if (GetIndex(temp_arr[1], Verts[apex[j]].connect) >= 0){ one_shared = true; }
			if (GetIndex(temp_arr[2], Verts[apex[j]].connect) >= 0){ two_shared = true; }
		}

#ifdef DEBUGGING
		if ((one_shared&&two_shared) || (!one_shared&&!two_shared)){
			//means both of the temp_arr[1] and temp_arr[2] (candidate apex)
			//are either shared with other apex's and not shared at all
			//have not considered this case yet
			ErrWarnMessage(__LINE__, "Operator::Injection:: not considered. See comment", 0);
		}
#endif

		apex[i_skip] = (two_shared) ? temp_arr[1] : temp_arr[2];
	}

	if (Verts[apex[1]].connect[0] == 3 || Verts[apex[2]].connect[0] == 3 || Verts[apex[3]].connect[0] == 3){
		//nop, we don't do this either 
		return false;
	}



	//*******2) Destory one edge of the triangle and another one 
	//there is 9 possibilities in near, we try them all :D

	int ip3 = ip[ip[0]];
	for (int i = 1; i <= ip[0]; i++){
		int j = (i == ip[0]) ? 1 : i + 1;

		int ip1 = ip[i];
		int ip2 = ip[j];
		int iapex = apex[i];
		int japex = apex[j];

		//we take the void_vertex as the avergae of two points 
		//point one is the mid-point between ip1-ip2
		//point one is the mid-point between ip2-ip3
		void_vertex[0] = (Verts[ip1].x[0] + 2.0*Verts[ip2].x[0] + Verts[ip3].x[0]) / 4.0;
		void_vertex[1] = (Verts[ip1].x[1] + 2.0*Verts[ip2].x[1] + Verts[ip3].x[1]) / 4.0;
		void_vertex[2] = (Verts[ip1].x[2] + 2.0*Verts[ip2].x[2] + Verts[ip3].x[2]) / 4.0;


		//#1 first try to destroy the edges ip1-ip2 ip2-ip3
		if ((i == 1 || i == ip[0]) && Verts[ip2].connect[0] > 4){ //this prevents ip2 from being tri-valent
			//(because we remove two verices connected to it and add one, effectively reduce valence by 1)		

			mynList[0] = 5;
			mynList[1] = ip1;
			mynList[2] = iapex;
			mynList[3] = ip2;
			mynList[4] = japex;
			mynList[5] = ip3;

			skipList[0] = 3;

			if (inj){
				SetAttractorRepeller(Verts);
				AttractorRepeller(ip, numVert, Verts, closedtSurfaceID, samplingBudget, numSurfaceLayer, false);
			}

			myConstraints.Reset(Verts, NULL);
			myConstraints.Direction(NULL, mynList, Verts);

			if (constraints.isDelaunay){
				//Delaunay constraints 
				myConstraints.Delaunay(mynList, skipList, Verts, false, 1);
			}

			if (myConstraints.OverlappingInSpheres()){
				if (constraints.isMinAngle || constraints.isMaxAngle){
					//Min max angle 
					myConstraints.MinMaxAngle(mynList, Verts, constraints.MinAngle, constraints.MaxAngle, void_vertex, constraints.isNonobtuse, 1);
				}

				if (constraints.isEdgeLen){
					//Min edge length 
					myConstraints.MinEdgeLength(mynList, Verts, constraints.MinEdgeLength_sq);
				}

				if (constraints.isSmooth){
					//Smoothness 
					myConstraints.Smoothness(mynList, skipList, Verts, void_vertex, constraints.dev);
				}


				StartActivePool(closedtSurfaceID, numSurfaceLayer);
				if (Sampler(Verts, mynList, samplingBudget, -1, -1, -1, &OpimizerFunc_CenterAngle)){
					//load edge to break
					//ip1-ip2
					//ip2-ip3
					EdgeToBreak[0] = 2;
					EdgeToBreak[1] = ip1;
					EdgeToBreak[2] = ip2;
					EdgeToBreak[3] = ip2;
					EdgeToBreak[4] = ip3;
					Inserter(numVert, Verts);
					return true;
				}
			}

			if (inj){
				RevertAttractorRepeller(Verts);
			}
		}



		//#2 second try to destroy the edges ip1-ip2 with another edge connected to on of the apex 

		for (int id = 1; id <= 2; id++){
			if (id == 2){
				std::swap(ip1, ip2);//trust me you need this.
			}

			if (Verts[ip1].connect[0] > 4){//this prevents ip1 from being tri-valent
				//(because we remove to verices connected to it and add one, effectively reduce valence by 1)

				//the two edges to destroy are ip1-ip2  and iapex1-ip1
				int extra_apex = FindCommonElement_SkipList(Verts[ip1].connect, Verts[iapex].connect, ip);

				if (i==2 && Verts[extra_apex].connect[0] > 3){//if it is 3-valent, then the void is ambigious 
					//the only solution would be creating 4-valent vertex which is equally bad for most cases
					mynList[0] = 5;
					mynList[1] = ip1;
					mynList[2] = extra_apex;
					mynList[3] = iapex;
					mynList[4] = ip2;
					mynList[5] = ip3;

					//here we take the void_vertex as the avergae of two points 
					//point one is the mid-point between ip1-iapex
					//point one is the mid-point between ip1-ip2
					void_vertex[0] = (2.0*Verts[ip1].x[0] + Verts[iapex].x[0] + Verts[ip2].x[0]) / 4.0;
					void_vertex[1] = (2.0*Verts[ip1].x[1] + Verts[iapex].x[1] + Verts[ip2].x[1]) / 4.0;
					void_vertex[2] = (2.0*Verts[ip1].x[2] + Verts[iapex].x[2] + Verts[ip2].x[2]) / 4.0;

					skipList[0] = 4;
					skipList[skipList[0]] = iapex;

					if (inj){
						SetAttractorRepeller(Verts);
						AttractorRepeller(ip, numVert, Verts, closedtSurfaceID, samplingBudget, numSurfaceLayer, false);
					}

					myConstraints.Reset(Verts, NULL);
					myConstraints.Direction(NULL, mynList, Verts);

					if (constraints.isDelaunay){
						//Delaunay constraints 
						myConstraints.Delaunay(mynList, skipList, Verts, false, 1);
					}

					if (myConstraints.OverlappingInSpheres()){
						if (constraints.isMinAngle || constraints.isMaxAngle){
							//Min max angle 
							myConstraints.MinMaxAngle(mynList, Verts, constraints.MinAngle, constraints.MaxAngle, void_vertex, constraints.isNonobtuse, 1);
						}

						if (constraints.isEdgeLen){
							//Min edge length 
							myConstraints.MinEdgeLength(mynList, Verts, constraints.MinEdgeLength_sq);
						}

						if (constraints.isSmooth){
							//Smoothness 
							myConstraints.Smoothness(mynList, skipList, Verts, void_vertex, constraints.dev);
						}


						StartActivePool(closedtSurfaceID, numSurfaceLayer);
						if (Sampler(Verts, mynList, samplingBudget, -1, -1, -1, &OpimizerFunc_CenterAngle)){
							EdgeToBreak[0] = 2;
							EdgeToBreak[1] = ip1;
							EdgeToBreak[2] = ip2;
							EdgeToBreak[3] = ip1;
							EdgeToBreak[4] = iapex;
							Inserter(numVert, Verts);
							return true;
						}
					}

					if (inj){
						RevertAttractorRepeller(Verts);
					}
				}
			}
		}


		ip3 = ip2;//not ip1 because we swap it 

	}



	return false;
}
void Operator::AttractorRepeller(int* ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer, bool isAttractor)
{
	//for each vertex in mynList, move it in eaither away (repeller) or closer (attractor)
	//in the direction of  void_vertex
	int prv_id = mynList[0];
	skipListNN[0] = 1;
	for (int i = 1; i <= mynList[0]; i++){
		int iq = mynList[i];
		//get the neighbour list of iq
		//this neighbour list is not a conected chain since iq is connected to ip's

		int nxt_id = (i == mynList[0]) ? 1 : i + 1; 		
		int start_id = GetIndex(mynList[prv_id], Verts[iq].connect);		
		int end_id = GetIndex(mynList[nxt_id], Verts[iq].connect);

#ifdef DEBUGGING
		if (start_id < 0 || end_id < 0){
			ErrWarnMessage(__LINE__, "Operator::Attractor:: can not grab the right index", 0);
		}
#endif

		mynListNN[0] = 0;
		if (isAttractor){
			CycleList_SkipList(iq, Verts, ip, start_id, end_id, mynListNN);
		}
		else{
			int before_prv = (prv_id == 1) ? mynList[mynList[0]] : mynList[prv_id - 1];
			int after_nxt = (nxt_id == mynList[0]) ? mynList[1] : mynList[nxt_id + 1];

			int move_away_vert;
			if (GetIndex(before_prv, Verts[iq].connect) >= 0){ 
				move_away_vert = before_prv;
			}
			else if (GetIndex(after_nxt, Verts[iq].connect) >= 0){
				move_away_vert = after_nxt;
			}
			else{
				move_away_vert = mynList[nxt_id];
			}
			
			CycleList(iq, Verts, move_away_vert, start_id, end_id, mynListNN, 1); 
		}
		if (mynListNN[0] < 0){
			prv_id = i;
			continue;
		}
			

		//build your geometric primitives to maintain the quality of triangles in mynListNN

		skipListNN[1] = iq;
		myConstraints.Reset(Verts, NULL);

		//myConstraints.Direction(skipListNN, mynListNN, Verts);
		myConstraints.Direction(skipListNN, Verts[iq].connect, Verts);

		if (constraints.isDelaunay){
			if (!myConstraints.Delaunay(mynListNN, skipListNN, Verts, 1, 0)){
				prv_id = i;
				continue;
			}
			
			if (!myConstraints.DelaunayForcedNotConnected(mynList, iq, skipList, Verts)){
				prv_id = i;
				continue;
			}

			if (!myConstraints.OverlappingInSpheres()){
				prv_id = i;
				continue;
			}
		}

		if (constraints.isMinAngle || constraints.isMaxAngle){
			//Min max angle 
			myConstraints.MinMaxAngle(mynListNN, Verts, constraints.MinAngle, constraints.MaxAngle, Verts[iq].x, constraints.isNonobtuse, 0);
		}
		
		if (constraints.isEdgeLen){
			//Min edge length 
			myConstraints.MinEdgeLength(mynListNN, Verts, constraints.MinEdgeLength_sq);
		}

								
		if (isAttractor){
			myConstraints.AttractorInSphere(Verts, iq, void_vertex);
		}
		else{
			//Hybird attractor/injection
			double sum_angle = 0;
			for (int kkk = 1; kkk < mynListNN[0]; kkk++){
				int n1 = mynListNN[kkk];
				int n2 = mynListNN[kkk + 1]; 
				sum_angle += AngleVectVect(Verts[n1].x[0] - Verts[iq].x[0], Verts[n1].x[1] - Verts[iq].x[1], Verts[n1].x[2] - Verts[iq].x[2],
					                       Verts[n2].x[0] - Verts[iq].x[0], Verts[n2].x[1] - Verts[iq].x[1], Verts[n2].x[2] - Verts[iq].x[2])*RadToDeg;
			}
			sum_angle = 360.0 - sum_angle;
			if (sum_angle >= 0){
				if (sum_angle < 2.0*constraints.MinAngle){
					myConstraints.AttractorInSphere(Verts, iq, void_vertex);
				}
				else if (sum_angle>2.0*constraints.MaxAngle){
					myConstraints.RepellerExSphere(Verts, iq, void_vertex);
				}
				else{
					//check to which it is closer (repeller or attracot)
					double dif_min = sum_angle - 2.0*constraints.MinAngle;
					double dif_max = 2.0*constraints.MaxAngle - sum_angle;
					if (dif_max < dif_min){
						myConstraints.RepellerExSphere(Verts, iq, void_vertex);
					}
					else{
						myConstraints.AttractorInSphere(Verts, iq, void_vertex);
					}
				}
			}
			else{
				//the default is to repel
				myConstraints.RepellerExSphere(Verts, iq, void_vertex);
			}
		}

		//we construct a single plane such that the new point is never behind the the void vertex
		//it helps in case of large voids

		myConstraints.SinglePlane(Verts[iq].x[0] - void_vertex[0],
			                      Verts[iq].x[1] - void_vertex[1],
			                      Verts[iq].x[2] - void_vertex[2],
			                      void_vertex[0], void_vertex[1], void_vertex[2]);

		StartActivePool(closedtSurfaceID, numSurfaceLayer);		
		if (Sampler(Verts, mynListNN, samplingBudget, mynList[nxt_id], iq, mynList[prv_id], (isAttractor) ? &OpimizerFunc_Closer : &OpimizerFunc_Further)){ 
			Mover(iq, Verts);
		}

		prv_id = i;

	}
}
void Operator::SetAttractorRepeller(vert*Verts)
{
	//copy the vertices we gonna move for attracor or repeller 
	for (int i = 1; i <= mynList[0]; i++){
		original_neighbout[i][0] = Verts[mynList[i]].x[0];
		original_neighbout[i][1] = Verts[mynList[i]].x[1];
		original_neighbout[i][2] = Verts[mynList[i]].x[2];
	}
}
void Operator::RevertAttractorRepeller(vert*Verts)
{
	//reset the vertices we may have moved for attractor or repeller 
	//probably because we moved some of them and inejction or ejection failed 
	for (int i = 1; i <= mynList[0]; i++){
		Verts[mynList[i]].x[0] = original_neighbout[i][0];
		Verts[mynList[i]].x[1] = original_neighbout[i][1];
		Verts[mynList[i]].x[2] = original_neighbout[i][2];
	}
}
void Operator::TriValent(int ip, vert*Verts, int*list, int skip1, int skip2)
{
	//store in list the trivalent nodes connected to ip
	//except if it skip1 or skip2
	for (int i = 1; i <= Verts[ip].connect[0]; i++){
		int iq = Verts[ip].connect[i];
		if (iq == skip1 || iq == skip2){ continue; }
		if (Verts[iq].connect[0] == 3){
			list[++list[0]] = iq;
		}
	}
}
void Operator::GetListSkipVertex(int ip, vert*Verts, int skip1, int skip2, int*list)
{
	//store the connectivity list of ip in list with exception of skip1 and skip2
	list[0] = 0;
	for (int i = 1; i <= Verts[ip].connect[0]; i++){
		if (Verts[ip].connect[i] == skip1 || Verts[ip].connect[i] == skip2){ continue; }
		list[++list[0]] = Verts[ip].connect[i];
	}
}
void Operator::GetEdgeSortedNeighbourList(int ip1, int ip2, vert*Verts, int*list)
{
	//get the sorted list of neighbour connected to ip1 and ip2 such that 
	//ip1 and ip2 are not tri-valent and they are not connected to shared tri-valent node 

	//1) find the common node between ip1 and ip2 (should be two)
	int common[10];
	FindCommonElements(Verts[ip1].connect, Verts[ip2].connect, common);
	if (common[0] != 2){
		return;		
	}

	int common1_id = GetIndex(common[1], Verts[ip1].connect);
	int common2_id = GetIndex(common[2], Verts[ip1].connect);
		
	//2) for ip1 connectivity list, start from common1 and keep advancing till you find common2
	//   for ip2 connectivity list, start from common2 and keep advancing till you find common1 
	//   since sorting is not consistent, we could be advancing forward or backword in the connectivity list
	//   such that we are always moving away from the other ip 

	list[0] = 0;
	
	CycleList(ip1, Verts, ip2, common1_id, common2_id, list, 0);
	if (list[0] < 0){ return; }

	common1_id = GetIndex(common[1], Verts[ip2].connect);
	common2_id = GetIndex(common[2], Verts[ip2].connect);

	CycleList(ip2, Verts, ip1, common2_id, common1_id, list, 0);


#ifdef DEBUGGING
	if (list[0] != Verts[ip1].connect[0] + Verts[ip2].connect[0] - 4){
		ErrWarnMessage(__LINE__, "Operator::GetEdgeSortedNeighbourList:: incorrect number of neighbours", 0);
	}
#endif

}
void Operator::GetFaceSortedNeighbourList(int ip1, int ip2, int ip3, vert*Verts, int*list)
{
	//get the sorted list of neighbour connected to ip1, ip2 and ip3 such that 
	//ip1, ip2 and ip3 are not tri-valent and they are not connected to shared tri-valent node 

	//1) find the common node between ip1 and ip2 (should be two)
	
	int ip1_ip2_sh = FindCommonElement_SkipList(Verts[ip1].connect, Verts[ip2].connect, ip3);
	int ip2_ip3_sh = FindCommonElement_SkipList(Verts[ip2].connect, Verts[ip3].connect, ip1);
	int ip3_ip1_sh = FindCommonElement_SkipList(Verts[ip3].connect, Verts[ip1].connect, ip2);

	if (ip1_ip2_sh < 0 || ip2_ip3_sh < 0 || ip3_ip1_sh < 0){
		list[0] = -1;
		return;
	}
	
	//2) for ip1 connectivity list, start from ip3_ip1_sh and keep advancing till you find ip1_ip2_sh
	//   for ip2 connectivity list, start from ip1_ip2_sh and keep advancing till you find ip2_ip3_sh 
	//   for ip3 connectivity list, start from ip2_ip3_sh and keep advancing till you find ip3_ip1_sh 
	//   since sorting is not consistent, we could be advancing forward or backword in the connectivity list
	//   such that when loop around ip1 we move away from ip3
    //             when loop around ip2 we move away from ip1
	//             when loop around ip3 we move away from ip2 	               

	list[0] = 0;

	int sh12_id_1 = GetIndex(ip1_ip2_sh, Verts[ip1].connect);
	int sh31_id_1 = GetIndex(ip3_ip1_sh, Verts[ip1].connect);
	CycleList(ip1, Verts, ip3, sh31_id_1, sh12_id_1, list, 0);
	if (list[0] < 0){ return; }

	int sh12_id_2 = GetIndex(ip1_ip2_sh, Verts[ip2].connect);
	int sh23_id_2 = GetIndex(ip2_ip3_sh, Verts[ip2].connect);
	CycleList(ip2, Verts, ip1, sh12_id_2, sh23_id_2, list, 0);
	if (list[0] < 0){ return; }

	int sh23_id_3 = GetIndex(ip2_ip3_sh, Verts[ip3].connect);
	int sh31_id_3 = GetIndex(ip3_ip1_sh, Verts[ip3].connect);
	CycleList(ip3, Verts, ip2, sh23_id_3, sh31_id_3, list, 0);
	if (list[0] < 0){ return; }


#ifdef DEBUGGING
	//if (list[0] != Verts[ip1].connect[0] + Verts[ip2].connect[0] - 4){
	//	ErrWarnMessage(__LINE__, "Operator::GetFaceSortedNeighbourList:: incorrect number of neighbours", 0);
	//}
#endif

}
void Operator::CycleList_SkipList(int ip, vert*Verts, int*SkipList, int start_id, int end_id, int*list)
{
	//start_id and end_id are  id's in ip connectivity list
	//starting from start_id, move in ip connectivity list till we reach end_id
	//skip all vertices in SkipList which may or may not be in ip connectivity
	//store what you get in list
	//list should contain at the end a connected list of verices but only the start and end vertices are not connected

	int prv = (start_id == 1) ? Verts[ip].connect[0] : start_id - 1;
	
	prv = GetIndex(Verts[ip].connect[prv], SkipList);
	
	bool forward = (prv < 0); 

	list[++list[0]] = Verts[ip].connect[start_id];
	int init_size = list[0];
	//predicate that we gonna move forward 
	//if we meet any of ip's before reaching end_id, we than reset and go the other direction
	int my_start_id = start_id;
	while (true){		
		my_start_id = (my_start_id == Verts[ip].connect[0]) ? 1 : my_start_id + 1;
		list[++list[0]] = Verts[ip].connect[my_start_id];
		if (GetIndex(Verts[ip].connect[my_start_id], SkipList) >= 0){ break; }
		if (my_start_id == end_id){ return; }
	}
	
	//if we reach here, then that means that moving forward was wrong 
	list[0] = init_size;
	while (true){				
		start_id = (start_id == 1) ? Verts[ip].connect[0] : start_id - 1;				
		list[++list[0]] = Verts[ip].connect[start_id];
		if (start_id == end_id){ return; }
	}
}
void Operator::CycleList(int ip, vert*Verts, int iq, int start_id, int end_id, int*list, bool includeLastEle)
{
	//start_id and end_id are  id's in ip connectivity list
	//starting from start_id, move in ip connectivity list till we reach end_id
	//move such that we are moving away from iq 
	
	int prv = (start_id == 1) ? Verts[ip].connect[0] : start_id - 1;
	int nxt = (start_id == Verts[ip].connect[0]) ? 1 : start_id + 1;

	if (Verts[ip].connect[prv] != iq && Verts[ip].connect[nxt] != iq){
		list[0] = -1;
		return;
		//ErrWarnMessage(__LINE__, "Operator::CycleList:: wrong id's", 0);
	}

	bool forward = (Verts[ip].connect[prv] == iq);

	list[++list[0]] = Verts[ip].connect[start_id];
	while (true){
		if (forward){
			start_id = (start_id == Verts[ip].connect[0]) ? 1 : start_id + 1;
		}
		else{
			start_id = (start_id == 1) ? Verts[ip].connect[0] : start_id - 1;
		}

		if (start_id == end_id){ break; }
		list[++list[0]] = Verts[ip].connect[start_id];
	}
	if (includeLastEle){
		list[++list[0]] = Verts[ip].connect[start_id];
	}
}
bool Operator::InspectFeasibleRegion(int*nList, vert*Verts)
{
	//here we inspect if there is feasible region 
	//by inspect the void corners only (rather than relying on resampling to figure this out)
	
	//For the three operators, we form a void corner and then insert a vertex inside it 
	//we check if a single vertex could possibly meet all the objectives or not 
	//for example, if the void corners make an angle less than 2*min_angle 
	//then any new vertex will split this angle to two angles such that one of them is less min_angle 
	return true;
	
	if (constraints.isMinAngle || constraints.isMaxAngle){
		int n2 = nList[nList[0]];//previous
		for (int i = 1; i <= nList[0]; i++){
			int ap = nList[i];//current
			int n1 = (i == nList[0]) ? nList[1] : nList[i + 1]; //next
			double xn, yn, zn;
			//double angle = AngleVectVect(Verts[n1].x[0] - Verts[ap].x[0], Verts[n1].x[1] - Verts[ap].x[1], Verts[n1].x[2] - Verts[ap].x[2],
			//	                         Verts[n2].x[0] - Verts[ap].x[0], Verts[n2].x[1] - Verts[ap].x[1], Verts[n2].x[2] - Verts[ap].x[2])*RadToDeg;//180.0 / PI;
			double v1x = Verts[n1].x[0] - Verts[ap].x[0]; double v1y = Verts[n1].x[1] - Verts[ap].x[1]; double v1z = Verts[n1].x[2] - Verts[ap].x[2];
			double v2x = Verts[n2].x[0] - Verts[ap].x[0]; double v2y = Verts[n2].x[1] - Verts[ap].x[1]; double v2z = Verts[n2].x[2] - Verts[ap].x[2];
			 
			Cross(v1x, v1y, v1z,
				  v2x, v2y, v2z,
				  xn, yn, zn);
			NormalizeVector(xn, yn, zn);
			double angle = Angle360Vectors(v1x, v1y, v1z, v2x, v2y, v2z, xn, yn, zn);
			if (angle > 180.0){
				int sdfd = 4545;
			}
			if (constraints.isMinAngle && angle < 2.0*constraints.MinAngle){
				return false;
			}
			if (constraints.isMaxAngle && angle > 2.0*constraints.MaxAngle){
				return false;
			}
			n2 = ap;
		}
	}
	return true;
		
}

///************** Sampling 
bool Operator::Sampler(vert*Verts, int*nList, int budget, 
	                   int i_af, int i_rep, int i_bf,
					   double(*Optimizer)(vert*Verts, int*nList, double*void_ver, double xx, double yy, double zz))
{
	
	//the resampling routine, called after constructing the geometric constraints 
	//and initlize the active pool
	//it should spit out a single new vertex in case of success and return true
	//otherwise return false

	//here we should distingiush between optimizer sampling and regular sampling 
	//optimizer seek to sample number of samples and pick the best to minimize the 
	//returned value from Optimizer()	
	//regular sampling just seek one sample to solve the problem 
	bool att = i_rep >= 0;

	
	double currentBest = DBL_MAX;//we gonna minimize this 
	if (budget > 1 && !att){ 
		//if optimizer, set the best to the current sample iff it meets all the constraints 
		if (CheckNewVertex(Verts, nList, Verts[skipList[1]].x[0], Verts[skipList[1]].x[1], Verts[skipList[1]].x[2], att)){
			currentBest = Optimizer(Verts, nList, void_vertex, Verts[skipList[1]].x[0], Verts[skipList[1]].x[1], Verts[skipList[1]].x[2]);

			newVertex[0] = Verts[skipList[1]].x[0];
			newVertex[1] = Verts[skipList[1]].x[1];
			newVertex[2] = Verts[skipList[1]].x[2];
		}
	}

	int num_succ_candidate(0);

	for (int lf = 0; lf < 15; lf++){
		int attempt = int(0.8*num_active);

		//a) dart throwing 
		
		for (int i = 0; i < attempt; i++){
			int tri = int(myRandNum.RandNumGenerator()*double(num_active - 1));
			double x1, y1, z1, x2, y2, z2, x3, y3, z3;
			RetrieveCoordinates(lf, active_pool[tri], x1, y1, z1, x2, y2, z2, x3, y3, z3);
			double x_new, y_new, z_new;
			RandomPointInTri(x1, y1, z1, x2, y2, z2, x3, y3, z3, x_new, y_new, z_new);

			
			if (CheckNewVertex(Verts, nList, x_new, y_new, z_new, att)){
				num_succ_candidate++;
				if (budget == 1){
					newVertex[0] = x_new;
					newVertex[1] = y_new;
					newVertex[2] = z_new;
					return true;
				}
				else{
					//apply optimizer 
					double myVal = Optimizer(Verts, nList, void_vertex, x_new, y_new, z_new);
					if (abs(currentBest - DBL_MAX) < _tol){
						//if this is the first candidate (then take it)
						currentBest = myVal;
						newVertex[0] = x_new;
						newVertex[1] = y_new;
						newVertex[2] = z_new;
					}
					else{
						
						if (myVal < currentBest){
							currentBest = myVal;
							newVertex[0] = x_new;
							newVertex[1] = y_new;
							newVertex[2] = z_new;
						}
					}
					if (num_succ_candidate >= budget){ 
						return true;
					}
				}
			}
		}

		//b) refinement 
		int tmp_num_active = 0;
		for (int iactive = 0; iactive < num_active; iactive++){
			double x1, y1, z1, x2, y2, z2, x3, y3, z3;
			RetrieveCoordinates(lf, active_pool[iactive], x1, y1, z1, x2, y2, z2, x3, y3, z3);

			int C = 0;
			_tar[C][0] = x1;          _tar[C][1] = y1;          _tar[C][2] = z1;
			_tar[C][3] = (x1 + x2) / 2.0; _tar[C][4] = (y1 + y2) / 2.0; _tar[C][5] = (z1 + z2) / 2.0;
			_tar[C][6] = (x1 + x3) / 2.0; _tar[C][7] = (y1 + y3) / 2.0; _tar[C][8] = (z1 + z3) / 2.0;
			C = 1;
			_tar[C][0] = (x1 + x2) / 2.0; _tar[C][1] = (y1 + y2) / 2.0; _tar[C][2] = (z1 + z2) / 2.0;
			_tar[C][3] = x2;          _tar[C][4] = y2;          _tar[C][5] = z2;
			_tar[C][6] = (x2 + x3) / 2.0; _tar[C][7] = (y2 + y3) / 2.0; _tar[C][8] = (z2 + z3) / 2.0;
			C = 2;
			_tar[C][0] = (x1 + x3) / 2.0; _tar[C][1] = (y1 + y3) / 2.0; _tar[C][2] = (z1 + z3) / 2.0;
			_tar[C][3] = (x2 + x3) / 2.0; _tar[C][4] = (y2 + y3) / 2.0; _tar[C][5] = (z2 + z3) / 2.0;
			_tar[C][6] = x3;          _tar[C][7] = y3;          _tar[C][8] = z3;
			C = 3;
			_tar[C][0] = (x1 + x3) / 2.0; _tar[C][1] = (y1 + y3) / 2.0; _tar[C][2] = (z1 + z3) / 2.0;
			_tar[C][3] = (x1 + x2) / 2.0; _tar[C][4] = (y1 + y2) / 2.0; _tar[C][5] = (z1 + z2) / 2.0;
			_tar[C][6] = (x2 + x3) / 2.0; _tar[C][7] = (y2 + y3) / 2.0; _tar[C][8] = (z2 + z3) / 2.0;

			for (C = 0; C < 4; C++){
			

				if (CheckRefinement(Verts, nList, _tar[C])){

					tmp_active_pool[tmp_num_active][0] = active_pool[iactive][0];
					tmp_active_pool[tmp_num_active][1] = active_pool[iactive][1];
					tmp_active_pool[tmp_num_active][2] = active_pool[iactive][2];

					if (lf>15){
						size_t shifted = C << 2 * (lf - 16);
						tmp_active_pool[tmp_num_active][2] = tmp_active_pool[tmp_num_active][2] | shifted;
					}
					else{
						size_t shifted = C << (2 * (lf));
						tmp_active_pool[tmp_num_active][1] = tmp_active_pool[tmp_num_active][1] | shifted;
					}
					tmp_num_active += 1;
					if (tmp_num_active + 2 > MAXPOOL){
						//cout<<"Too much refining in Resampling()"<<endl;
						if (num_succ_candidate > 0){ return true; }
						return false;

					}
				}
			}
		}

		//swap
		if (tmp_num_active == 0){ 
			if (num_succ_candidate > 0){ return true; }
			return false; 
		}
		else{
			num_active = tmp_num_active;
			for (int iactive = 0; iactive < num_active; iactive++){
				active_pool[iactive][0] = tmp_active_pool[iactive][0];
				active_pool[iactive][1] = tmp_active_pool[iactive][1];
				active_pool[iactive][2] = tmp_active_pool[iactive][2];
			}
		}
	}

	return false;//why did we get here !!!

}
void Operator::StartActivePool(int id, int maxDepth)
{
	//TODO pass me the kd-tree of the original surface and let me pick the closest the surface myself

	//Start the active pool by the fan triangle of id
	//then propagate to include maxDepth of surrounding triangles of above fan triangle 

	//for relocation, it is better to let maxDepth to be 1
	//for ejection, it is better to get larger patch tho 2 would be enough 
	
	num_active = 0;
	next_layer[0] = 1;
	next_layer[1] = id;
	int start = 1;
	for (int lay = 0; lay < maxDepth; lay++){
		int end = next_layer[0];
		for (int i = start; i <=end; i++){
			int myID = next_layer[i];
			for (int j = 1; j <= Vert_org[myID].connect[0]; j++){
				int candid = Vert_org[myID].connect[j];
				if (GetIndex(candid, next_layer) < 0){
					next_layer[++next_layer[0]] = candid;
#ifdef DEBUGGING
					if (next_layer[0] + 1 > MAXPOOL){
						ErrWarnMessage(__LINE__, "Operator::StartActivePool:: increase length of next_layer", 0);
					}
#endif
				}
			}
		}
		start = end + 1; 
	}
	
	

	for (int i = 1; i <= next_layer[0]; i++){		
		GetFanTriangles(next_layer[i]);
	}
		
}
void Operator::StartActivePool(int closedtSurfaceID, double*myVert)
{
	//initialize actove pool but triangle fans of some surface vertices 
	//find the surface trianle such that the projection of myVert is min
	//get the fan triangle of the three vertices of this triangle 

	//myVert is the coordinates of the center to be projected 
	//for relocation, it is the vertex to be reolcated
	//for ejection/injection, it the mid point of the two points to be ejected, or center of triangle to be ejected
	num_active = 0;

	int closedtSurfaceID1, closedtSurfaceID2;
	int n2 = Vert_org[closedtSurfaceID].connect[Vert_org[closedtSurfaceID].connect[0]];
	double closestDist = DBL_MAX;
	for (int i = 1; i <= Vert_org[closedtSurfaceID].connect[0]; i++){
		int n1 = Vert_org[closedtSurfaceID].connect[i];
		double x_projected, y_projected, z_projected;
		double dist = PointTriangleDistance(Vert_org[closedtSurfaceID].x[0], Vert_org[closedtSurfaceID].x[1], Vert_org[closedtSurfaceID].x[2],
			                                Vert_org[n1].x[0], Vert_org[n1].x[1], Vert_org[n1].x[2],
			                                Vert_org[n2].x[0], Vert_org[n2].x[1], Vert_org[n2].x[2],
			                                myVert[0], myVert[1], myVert[2],
			                                x_projected, y_projected, z_projected);
		if (dist < closestDist){
			closestDist = dist;
			closedtSurfaceID1 = n1;
			closedtSurfaceID2 = n2;
		}
		n2 = n1;
	}


	//get the fan triangles of closedtSurfaceID, closedtSurfaceID1 and closedtSurfaceID2 
	//carful of duplication 
	GetFanTriangles(closedtSurfaceID);
	GetFanTriangles(closedtSurfaceID1);
	GetFanTriangles(closedtSurfaceID2);

}
void Operator::GetFanTriangles(int id)
{
	//get the fan triangles in id and store them in active_pool 
	//make sure there is no duplication 
	int n2 = Vert_org[id].connect[Vert_org[id].connect[0]];
	
	for (int i = 1; i <= Vert_org[id].connect[0]; i++){
		
		int n1 = Vert_org[id].connect[i];
		
		//id,n1,n2 
		//check duplication 
		bool dup = false;
		for (int j = 0; j < num_active; j++){
			if (IsDuplicated(n2, n1, id, j)){
				dup = true;
				break;
			}
		}
			
		if (!dup){
			tri_pool[num_active][0] = id;
			tri_pool[num_active][1] = n1;
			tri_pool[num_active][2] = n2;

			active_pool[num_active][0] = num_active;
			active_pool[num_active][1] = 0;
			active_pool[num_active][2] = 0;
			num_active++;

			if (num_active + 2 >= MAXPOOL){ return; }
			
		}

		n2 = n1;
	}
}
bool Operator::IsDuplicated(int ip1, int ip2, int ip3, int t2)
{
	return (ip1 == tri_pool[t2][0] || ip1 == tri_pool[t2][1] || ip1 == tri_pool[t2][2]) &&
		   (ip2 == tri_pool[t2][0] || ip2 == tri_pool[t2][1] || ip2 == tri_pool[t2][2]) &&
		   (ip3 == tri_pool[t2][0] || ip3 == tri_pool[t2][1] || ip3 == tri_pool[t2][2]);

}
void Operator::RandomPointInTri(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double&x, double&y, double&z)
{
	double r1, r2;

	r1 = myRandNum.RandNumGenerator();
	r2 = myRandNum.RandNumGenerator();

	x = (1 - sqrt(r1))*x1 + (sqrt(r1)*(1 - r2))*x2 + (r2*sqrt(r1))*x3;
	y = (1 - sqrt(r1))*y1 + (sqrt(r1)*(1 - r2))*y2 + (r2*sqrt(r1))*y3;
	z = (1 - sqrt(r1))*z1 + (sqrt(r1)*(1 - r2))*z2 + (r2*sqrt(r1))*z3;

}
bool Operator::CheckRefinement(vert*Verts, int*nList, double*tri)
{
	

	if (!myConstraints.InsideFeasibleRegion_Triangle(tri)){ return false; }

	/*if (constraints.isMinAngle || constraints.isMaxAngle){ //this is wrong 
		//TODO include the attractor here 
		int n2 = nList[nList[0]];
		for (int i = 1; i <= nList[0]; i++){
			int n1 = nList[i];
			double 	angle1 = (AngleVectVect(Verts[n1].x[0] - tri[0], Verts[n1].x[1] - tri[1], Verts[n1].x[2] - tri[2],
				                            Verts[n2].x[0] - tri[0], Verts[n2].x[1] - tri[1], Verts[n2].x[2] - tri[2]))*RadToDeg;

			double angle2 = (AngleVectVect(Verts[n1].x[0] - tri[3], Verts[n1].x[1] - tri[4], Verts[n1].x[2] - tri[5],
				                           Verts[n2].x[0] - tri[3], Verts[n2].x[1] - tri[4], Verts[n2].x[2] - tri[5]))*RadToDeg;

			double angle3 = (AngleVectVect(Verts[n1].x[0] - tri[6], Verts[n1].x[1] - tri[7], Verts[n1].x[2] - tri[8],
				                           Verts[n2].x[0] - tri[6], Verts[n2].x[1] - tri[7], Verts[n2].x[2] - tri[8]))*RadToDeg;
			if ((angle1 < constraints.MinAngle + _tol &&
				 angle2 < constraints.MinAngle + _tol && 
				 angle3 < constraints.MinAngle + _tol)||
				(angle1 > constraints.MaxAngle + _tol &&
				 angle2 > constraints.MaxAngle + _tol &&
				 angle3 > constraints.MaxAngle + _tol)){
				return false;
			}
			n2 = n1;
		}
	}*/
	return true;

}
bool Operator::CheckNewVertex(vert*Verts, int*nList, double x_new, double  y_new, double  z_new, bool att)
{
	//chech if a new vertex is an acceptable vertex (outside all exclusion regions,
	//and inside all inclusion regions)
	if (!myConstraints.InsideFeasibleRegion_Vertex(x_new, y_new, z_new)){ return false; }

	if (constraints.isMinAngle || constraints.isMaxAngle){
		//check the apex angle

		if (att){
			double accum = 0;
			for (int i = 1; i < nList[0]; i++){
				int n2 = nList[i];
				int n1 = nList[i + 1]; 
				double myAngle = AngleVectVect(Verts[n1].x[0] - x_new, Verts[n1].x[1] - y_new, Verts[n1].x[2] - z_new,
					                           Verts[n2].x[0] - x_new, Verts[n2].x[1] - y_new, Verts[n2].x[2] - z_new)*RadToDeg;

				if (myAngle < constraints.MinAngle + _tol){ return false; }
				if (myAngle + _tol > constraints.MaxAngle){ return false; }
				accum += myAngle;
			}			
		}
		else{
			int n2 = nList[nList[0]];
			for (int i = 1; i <= nList[0]; i++){
				int n1 = nList[i];

				double myAngle = AngleVectVect(Verts[n1].x[0] - x_new, Verts[n1].x[1] - y_new, Verts[n1].x[2] - z_new,
					                            Verts[n2].x[0] - x_new, Verts[n2].x[1] - y_new, Verts[n2].x[2] - z_new)*RadToDeg;

				if (myAngle < constraints.MinAngle + _tol || myAngle + _tol > constraints.MaxAngle){ return false; }				
				n2 = n1;
			}
		}
	}

	if (constraints.isSmooth){
		//TODO add smoothness to attractor 

		//check smoothnes 
		int n2 = nList[nList[0]];
		for (int i = 1; i <= nList[0]; i++){
			int n1 = nList[i];
			int n3 = (i == nList[0]) ? nList[1] : nList[i + 1];
			double angle = TriTriNormalAngle3(x_new, y_new, z_new, Verts[n1].x, Verts[n2].x, Verts[n3].x);
			//angle = 180 - angle;
			//if (angle < 180 - constraints.dev){
			if (angle > constraints.dev){ 
				return false;
			}
			n2 = n1;
		}
	}
	return true;
}
void Operator::RetrieveCoordinates(int lf, int* cell, double&x0, double&y0, double&z0, double&x1, double&y1, double&z1, double&x2, double&y2, double&z2)
{
	int po = tri_pool[cell[0]][0];
	int p1 = tri_pool[cell[0]][1];
	int p2 = tri_pool[cell[0]][2];

	x0 = Vert_org[po].x[0];
	y0 = Vert_org[po].x[1];
	z0 = Vert_org[po].x[2];
	x1 = Vert_org[p1].x[0];
	y1 = Vert_org[p1].x[1];
	z1 = Vert_org[p1].x[2];
	x2 = Vert_org[p2].x[0];
	y2 = Vert_org[p2].x[1];
	z2 = Vert_org[p2].x[2];

	size_t binary = cell[1], c_num;//c_num is the ID of refiend triangle
	//binary is the avariable holding the value of the cell[1 or 2] which will be "bit" masked to get the refined ID
	int shift = 0;
	double a0, b0, c0, a1, b1, c1, a2, b2, c2;
	for (int i = 0; i < lf; i++){ 
		//get refined tri ID
		c_num = binary& ((int(pow(2.0, (shift + 1.0))) + 1) | (int(pow(2.0, shift)) + 1));
		c_num = c_num >> shift;
		//cout<<c_num<<"\n";
		//switch and get new coordinates of C
		switch (c_num)
		{
		case 0:

			x1 = (x0 + x1) / 2.0;
			y1 = (y0 + y1) / 2.0;
			z1 = (z0 + z1) / 2.0;

			x2 = (x0 + x2) / 2.0;
			y2 = (y0 + y2) / 2.0;
			z2 = (z0 + z2) / 2.0;
			break;
		case 1:

			x0 = (x0 + x1) / 2.0; y0 = (y0 + y1) / 2.0; z0 = (z0 + z1) / 2.0;

			x2 = (x1 + x2) / 2.0; y2 = (y1 + y2) / 2.0; z2 = (z1 + z2) / 2.0;
			break;
		case 2:

			x0 = (x0 + x2) / 2.0; y0 = (y0 + y2) / 2.0; z0 = (z0 + z2) / 2.0;

			x1 = (x1 + x2) / 2.0; y1 = (y1 + y2) / 2.0; z1 = (z1 + z2) / 2.0;
			break;
		case 3://in c3 I pick Po as the (right side direction point) and then counter clockwise

			a0 = (x0 + x2) / 2.0; b0 = (y0 + y2) / 2.0; c0 = (z0 + z2) / 2.0;

			a1 = (x0 + x1) / 2.0; b1 = (y0 + y1) / 2.0; c1 = (z0 + z1) / 2.0;

			a2 = (x1 + x2) / 2.0; b2 = (y1 + y2) / 2.0; c2 = (z1 + z2) / 2.0;

			x0 = a0; y0 = b0; z0 = c0; x1 = a1; y1 = b1; z1 = c1; x2 = a2; y2 = b2; z2 = c2;
			break;
		default:
			fprintf(stderr, "\nError.\n");
		}

		if (i == 15){//switching to the second cell after 16 itaration of refinement
			binary = cell[2];
			shift = 0;
		}
		shift = shift + 2;
	}

}


///************** Optimizer Functions (for sampling) 
double OpimizerFunc_CenterAngle(vert*Verts, int*nList, double*void_verx, double xx, double yy, double zz)
{
	//return the different between myAngle and perfect_angle 
	//perfect_angle = 360.0 / nList[0];
	//myAngle is the min apex angle 

	double perfect_angle = 360.0 / nList[0];

	double myAngle = 360.0;
	int n2 = nList[nList[0]];
	for (int i = 1; i <= nList[0]; i++){
		int n1 = nList[i];
		double ang = AngleVectVect(Verts[n1].x[0] - xx, Verts[n1].x[1] - yy, Verts[n1].x[2] - zz,
			                       Verts[n2].x[0] - xx, Verts[n2].x[1] - yy, Verts[n2].x[2] - zz)*RadToDeg;
		myAngle = std::min(myAngle, ang);
		n2 = n1;
	}
	return abs(myAngle - perfect_angle);
}
double OpimizerFunc_SideAngle(vert*Verts, int*nList, double minAngle, double maxAngle, double xx, double yy, double zz)
{
	//get the angle between (xx,yy,zz), nList[1] and nList[nList[0]]
	//if it supposed to be less than 180.0
	//otherwise this method might returen wrong angle 

	//return the difference between (2minAngle - ang) + (ang - 2maxAngle)
	//we take the difference this way because we try to minimize the returned value from this function

	int n1 = nList[1];
	int n2 = nList[nList[0]];
	double ang = AngleVectVect(Verts[n1].x[0] - xx, Verts[n1].x[1] - yy, Verts[n1].x[2] - zz,
                   		       Verts[n2].x[0] - xx, Verts[n2].x[1] - yy, Verts[n2].x[2] - zz)*RadToDeg;
	
	
	double diff1 = 2.0*minAngle - ang;
	double diff2 = ang - 2.0*maxAngle;

	
	return diff2 + diff2; 
}
double OpimizerFunc_SideAngleHybird(vert*Verts, int*nList, double minAngle, double maxAngle, double xx, double yy, double zz)
{
	//TODO (smarter way for optimizing side angles)
	//use this for hybird attractor/repeller 

	//nList is discounted list of neighbours around (xx,yy,zz)
	//get sum of all apex angle at (xx,yy,zz)
	//subtract this from 360 --> sum_angle 
	//depending min(diff_min, diff_max)
	//where diff_min = sum_angle - 2.0*minAngle
	//diff_max = 2.0*maxAngle - sum_angle 
	double sum_angle = 0;
	for (int i = 1; i < nList[0]; i++){
		int n1 = nList[i];
		int n2 = nList[i + 1]; 
		sum_angle += AngleVectVect(Verts[n1].x[0] - xx, Verts[n1].x[1] - yy, Verts[n1].x[2] - zz,
			                       Verts[n2].x[0] - xx, Verts[n2].x[1] - yy, Verts[n2].x[2] - zz)*RadToDeg;
	}
	sum_angle = 360.0 - sum_angle;

	double diff_min = -(sum_angle - 2.0*minAngle); 
	double diff_max = -(2.0*maxAngle - sum_angle);

	return std::min(diff_min, diff_max);

}
double OpimizerFunc_Closer(vert*Verts, int*nList, double*void_verx, double xx, double yy, double zz)
{
	//return the distance to void_vertex 
	//since we optimizer minimize the distance, we return it as it is 

	double dist = Dist(void_verx[0], void_verx[1], void_verx[2], xx, yy, zz);
	return dist;
}
double OpimizerFunc_Further(vert*Verts, int*nList, double*void_verx, double xx, double yy, double zz)
{
	//return the distance to void_vertex 
	//since we optimizer minimize the distance, we return negative the distance 

	double dist = Dist(void_verx[0], void_verx[1], void_verx[2], xx, yy, zz);
	return -dist; 
}



///************** Update Mesh 
void Operator::Mover(int ip, vert*Verts)
{
	//move ip to newVertex (do any other neccessary toplogical changes)
	Verts[ip].x[0] = newVertex[0];
	Verts[ip].x[1] = newVertex[1];
	Verts[ip].x[2] = newVertex[2];
}
void Operator::EdgeCollapse(int ip_low, int ip_high, int&numVert, vert*Verts)
{
	//collapse the edge between ip_low and ip_high to a point (newVertex)
	//remove ip_high and then move ip_low (the order of the operation is i m p o r t a nt)
	if (ip_low > ip_high){ std::swap(ip_high, ip_low); }

	//1) remove ip_high 
	RemoveVertex(ip_high, numVert, Verts);

	//2) move ip_low
	Verts[ip_low].x[0] = newVertex[0];
	Verts[ip_low].x[1] = newVertex[1];
	Verts[ip_low].x[2] = newVertex[2];
	//remove ip_low from any connectivity list of other vertices
	for (int i = 1; i <= Verts[ip_low].connect[0]; i++){
		int iq = Verts[ip_low].connect[i];
		RemoveNodeFromList(Verts[iq].connect, ip_low);
	}

	//update the connectivity of ip_low to be mynList
	for (int i = 0; i <= mynList[0]; i++){
		Verts[ip_low].connect[i] = mynList[i];
	}

	//appropriately, add ip_low to the connectivity list of vertices on mynList
	int prv = mynList[mynList[0]];//prv
	for (int i = 1; i <= mynList[0]; i++){ 
		int nxt = (i == mynList[0]) ? mynList[1] : mynList[i + 1]; //next
		AddEntrySortedList(ip_low, Verts[mynList[i]].connect, prv, nxt);
		prv = mynList[i];
	}


}
void Operator::FaceCollapse(int ip_low, int ip_mid, int ip_high, int&numVert, vert*Verts)
{
	//collapse the triangle ip_low-ip_mid-ip_high
	//remove ip_high and then remove ip_mid and then move ip_low (the order of the operation is i m p o r t a nt)
	
	int my_ip_high, my_ip_mid, my_ip_low;
	my_ip_high = std::max(ip_low, ip_mid);
	my_ip_high = std::max(my_ip_high, ip_high);

	my_ip_low = std::min(ip_low, ip_mid);
	my_ip_low = std::min(my_ip_low, ip_high);

	my_ip_mid = (ip_mid > my_ip_low && ip_mid<my_ip_high) ? ip_mid : ((ip_low>my_ip_low && ip_low < my_ip_high) ? ip_low : ip_high);

	ip_low = my_ip_low; ip_high = my_ip_high; ip_mid = my_ip_mid;

	//1) remove ip_high, ip_mid
	RemoveVertex(ip_high, numVert, Verts);
	RemoveVertex(ip_mid, numVert, Verts);

	//2) move ip_low
	Verts[ip_low].x[0] = newVertex[0];
	Verts[ip_low].x[1] = newVertex[1];
	Verts[ip_low].x[2] = newVertex[2];
	//remove ip_low from any connectivity list of other vertices
	for (int i = 1; i <= Verts[ip_low].connect[0]; i++){
		int iq = Verts[ip_low].connect[i];
		RemoveNodeFromList(Verts[iq].connect, ip_low);
	}

	//update the connectivity of ip_low to be mynList
	for (int i = 0; i <= mynList[0]; i++){
		Verts[ip_low].connect[i] = mynList[i];
	}

	//appropriately, add ip_low to the connectivity list of vertices on mynList
	int prv = mynList[mynList[0]];//prv
	for (int i = 1; i <= mynList[0]; i++){
		int nxt = (i == mynList[0]) ? mynList[1] : mynList[i + 1]; //next
		AddEntrySortedList(ip_low, Verts[mynList[i]].connect, prv, nxt);
		prv = mynList[i];
	}

}
void Operator::RemoveVertex(int ip, int&numVert, vert*Verts)
{
	//remove ip from mesh data strcuture 

	//1) remove ip from any connectivity list of other vertices
	for (int i = 1; i <= Verts[ip].connect[0]; i++){
		int iq = Verts[ip].connect[i];
		RemoveNodeFromList(Verts[iq].connect, ip);
	}

	//2)grab the last mesh vertex and promote it to be ip
	//if ip was already the last mesh vertex, then never mind this :)
	numVert--;
	if (ip!=numVert){
		for (int i = 0; i <= Verts[numVert].connect[0]; i++){
			Verts[ip].connect[i] = Verts[numVert].connect[i];
		}
		Verts[ip].x[0] = Verts[numVert].x[0];
		Verts[ip].x[1] = Verts[numVert].x[1];
		Verts[ip].x[2] = Verts[numVert].x[2];

		//3) any vertex connected to numVert should know it new name -> ip
		for (int i = 1; i <= Verts[numVert].connect[0]; i++){
			int iq = Verts[numVert].connect[i];
			int id = GetIndex(numVert, Verts[iq].connect);
#ifdef DEBUGGING
			if (id < 0){
				ErrWarnMessage(__LINE__, "Operator::RemoveVertex error(1)", 0);
			}
#endif
			Verts[iq].connect[id] = ip;
		}
	}

	//numVert could be also there in mynList (which will be used in update another mesh connectivity in case of edge collapse)
	//so, we gotta update this one too
	int id = GetIndex(numVert, mynList);
	if (id >= 0){
		mynList[id] = ip;
	}

	
}
void Operator::Inserter(int&numVert, vert*Verts)
{
	//break edges in EdgeToBreak,
	//insert a newVertex into the data structure and connect it to mynList,
	//update mynList too, 
	
	for (int i = 1; i <= EdgeToBreak[0]; i++){
		int ip = EdgeToBreak[i * 2 - 1]; 
		int iq = EdgeToBreak[i * 2]; 
		RemoveNodeFromList(Verts[iq].connect, ip);
		RemoveNodeFromList(Verts[ip].connect, iq);
	}

	Verts[numVert].x[0] = newVertex[0];
	Verts[numVert].x[1] = newVertex[1];
	Verts[numVert].x[2] = newVertex[2];
		
	for (int i = 0; i <= mynList[0]; i++){
		Verts[numVert].connect[i] = mynList[i];
	}


	//appropriately, add numVert  to the connectivity list of vertices on mynList
	int prv = mynList[mynList[0]];//prv
	for (int i = 1; i <= mynList[0]; i++){
		int nxt = (i == mynList[0]) ? mynList[1] : mynList[i + 1]; //next
		AddEntrySortedList(numVert, Verts[mynList[i]].connect, prv, nxt);
		prv = mynList[i];
	}
		
	numVert++;

}


///************** Debug 
void Operator::ErrWarnMessage(size_t lineNum, std::string message, size_t mess_id)
{
	//mess_id =0 for error (exit)
	//otherwise, it is a warning (pause)

	if (mess_id == 0){
		fprintf(stderr, "\nError::line(%d)-->>%s", lineNum, message.c_str());
		system("pause");
	}
	else{
		fprintf(stderr, "\nWarning::line(%d)-->>%s\n", lineNum, message.c_str());
		system("pause");
	}
}
void Operator::DrawActivePool(int lf)
{
	std::cout << "\n I AM DRAWING in Operator::DrawActivePool()" << std::endl;

	std::stringstream fname;
	fname << "debug_out/active_pool.obj";
	std::fstream file(fname.str().c_str(), std::ios::out);
	double x1, y1, z1, x2, y2, z2, x3, y3, z3;

	
	for (int V = 0; V<num_active; V++){		
		RetrieveCoordinates(lf, active_pool[V], x1, y1, z1, x2, y2, z2, x3, y3, z3);
		file << "v " << x1 << " " << y1 << " " << z1 << std::endl;
		file << "v " << x2 << " " << y2 << " " << z2 << std::endl;
		file << "v " << x3 << " " << y3 << " " << z3 << std::endl;
	}
	for (int V = 1; V <= 3 * num_active; V += 3){
		file << "f " << V << " " << V + 1 << " " << V + 2 << std::endl;
	}
}
void Operator::DrawnList(vert*Verts)
{
	for (int i = 1; i <= mynList[0]; i++){
		//DrawOneSphere(mynList[i], Verts[mynList[i]].x[0], Verts[mynList[i]].x[1], Verts[mynList[i]].x[2], 0.001, 3);
	}
	
}