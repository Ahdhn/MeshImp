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

//The operators are called from the execution 
//they should be called wiht the right constraints
//the operators includes all the modification neccessary 
//i.e., sampling and update data structure  


#ifndef _OPERATOR_
#define _OPERATOR_
#define MAXPOOL 2000

#include <stdlib.h>

#include "../constraint/Constraints.h"
#include "../util/Common.h"
#include "../util/RNG.h"


class Operator
{
public:
	//the constraints are set once by the executer  
	struct AppConstraints{
		bool isDelaunay;
		bool isMinAngle; double MinAngle;
		bool isMaxAngle; double MaxAngle;
		bool isNonobtuse;
		bool isSmooth; double dev;
		bool isEdgeLen; double MinEdgeLength_sq;
		                                     //need to look into this in case of function is passed
		bool isMaximal;
	} constraints;

	Operator(vert*Vert_org);
	~Operator();

	void TriValentRemoval(int&numVert, vert*Verts);

	bool Relocation(int ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer);

	bool Ejection(int* ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer, bool att);

	bool Injection(int*ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer, bool inj);

	bool AggressiveInjection(int*ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer, bool inj);



private:
	Constraints myConstraints;  	
	RndNum myRandNum;
	
	int skipList[20];//the skip list stores the vertices that should be considered as if they are deleted 
	int skipListNN[20];//the skip list stores the vertices that should be considered as if they are deleted 
	                 //we use this when we are attracting/repelling vertices 
	double newVertex[3];//coordinates of the new vertex (updated from Sampler)
	int mynList[100];//sorted neighbour list (connected chain)
	int mynListNN[100];//sorted neighbour list for being attracted vertex (unconnected chain)
	double void_vertex[3];
	int apex[4];//for each two consecutive vertices in *ip, find their (correct) shared vertex 
	int EdgeToBreak[10];//max 3 edges (3*2 + 1)
	double original_neighbout[100][3];//keeps a copy of the original void corners to be updated in case of failure in repeller or attractor 
	int temp_arr[100];//used for various reasons when a temp array is needed 

	void AttractorRepeller(int* ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer, bool isAttractor);
	void SetAttractorRepeller(vert*Verts);
	void RevertAttractorRepeller(vert*Verts);

	void StartActivePool(int closedtSurfaceID, double*myVert);
	void StartActivePool(int id, int numLayers);

	void GetFanTriangles(int);
	bool IsDuplicated(int ip1, int ip2, int ip3, int t2);
	bool Sampler(vert*Verts, int*nList, int budget,
		         int i_af, int i_rep, int i_bf,
				 double(*Optimizer)(vert*Verts, int*nList, double*void_ver, double xx, double yy, double zz));

	void RetrieveCoordinates(int lf, int* cell, double&x0, double&y0, double&z0, double&x1, double&y1, double&z1, double&x2, double&y2, double&z2);	
	bool CheckNewVertex(vert*Verts, int*nList, double x_new, double  y_new, double  z_new, bool att);
	bool CheckRefinement(vert*Verts, int*nList, double*tri);
	void RandomPointInTri(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double&x, double&y, double&z);
	void Mover(int ip, vert*Verts);
	void EdgeCollapse(int ip_low, int ip_high, int&numVert, vert*Verts);
	void FaceCollapse(int ip_low, int ip_mid, int ip_high, int&numVert, vert*Verts);	
	void Inserter(int&numVert, vert*Verts);
	void RemoveVertex(int ip, int&numVert, vert*Verts);
	bool InspectFeasibleRegion(int*nList, vert*Verts);

	void ErrWarnMessage(size_t lineNum, std::string message, size_t mess_id);
	void TriValent(int ip, vert*Verts, int*list, int skip1, int skip2);
	void GetListSkipVertex(int ip, vert*Verts, int skip1, int skip2, int*list);
	void GetEdgeSortedNeighbourList(int ip1, int ip2, vert*Verts, int*list);
	void GetFaceSortedNeighbourList(int ip1, int ip2, int ip3, vert*Verts, int*list);
	void CycleList(int ip, vert*Verts, int iq, int start_id, int end_id, int*list, bool includeLastEle);
	void CycleList_SkipList(int ip, vert*Verts, int*SkipList, int start_id, int end_id, int*list);
	

	void DrawActivePool(int lf);
	void DrawnList(vert*Verts);
	
	vert*Vert_org;//input surface 

	//resampling grid
	int num_active;
	int**active_pool, **tmp_active_pool, **tri_pool;
	double **_tar;
	int next_layer[MAXPOOL];
};



#endif /*_OPERATOR_*/