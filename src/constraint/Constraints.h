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

//Constraints is supposed to be called from the Operators 
//it has all the allocation spaces needed and will be accessible by the 
//the Sampler. Notice that constraints are confined to single patch 

#ifndef _CONSTRAINTS_
#define _CONSTRAINTS_
#include "../util/Common.h"
#include <string>



class Constraints
{

public:
	Constraints();
	~Constraints();
	void Reset(vert*Verts, int*nList);
	bool Delaunay(int*nList, int*skipList, vert*Verts, bool loadSkipList, bool isConnected);
	bool DelaunayForcedNotConnected(int*nList, int removed, int*skipList, vert*Verts);
	void MinMaxAngle(int*nList, vert*Verts, double min_ang, double max_ang, double*void_vertex, bool nonObtuse, bool isConnected);
	void MinEdgeLength(int*nList, vert*Verts, double (*sizingfunc)(double xx, double yy, double zz));
	void MinEdgeLength(int*nList, vert*Verts, double r_min_2);
	void Smoothness(int*nList, int*skipList, vert*Verts, double*void_vertex, double dev);	
	void Direction(int*ip, int*nList, vert*Verts);
	void Maximality();
	void AttractorInSphere(vert*Verts, int ip, double*void_vertex);
	void RepellerExSphere(vert*Verts, int ip, double*void_vertex);

	bool InsideFeasibleRegion_Vertex(double xx, double yy, double zz);
	bool InsideFeasibleRegion_Triangle(double*tri);
	bool OverlappingInSpheres();
	void SinglePlane(double xn, double yn, double zn, double px, double py, double pz);
	

	//debug
	void DrawExSpheres();
	void DrawInSpheres();

private:

	int numExSphere, numInSphere, numPlane;
	int numExSphere_size, numInSphere_size, numPlane_size;
	int aux_list[20];

	sphere *ExSphere, *InSphere;
	plane *Plane;
	BoundBox myBox; //specify the bounding box using the the nList 
	                //new vertex should be inside the bounding box 
	void SetBoundingBox(vert*Verts, int*nList);
	bool DelaunayNotConnected(int*nList, int*skipList, vert*Verts, bool loadSkipList);
	void ExpandSpheres(int&currentSize, int currentNumSphere, sphere*&mySpheres);
	void ExpandPlanes(int&currentSize, int currentNumPlane, plane*&myPlanes);
	void ErrWarnMessage(size_t lineNum, std::string message, size_t mess_id);
	bool IsEmptySphere(double xc, double yc, double zc, double rc_2,
		               int*list,
		               int*skip_list,
		               vert*Verts);
	bool IsInsideBoundingBox(double xx, double yy, double zz);
	//debug
	void SpheresDrawing(std::string filename, int num, sphere *Sphere);
};



#endif 