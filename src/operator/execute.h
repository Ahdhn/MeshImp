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

//This class contains the actual implementation for mesh improvement 
//it starts with by initialize the basic data structures (containers, kd tree, etc)
//it supports the non-obtuse remeshing and simplification


#ifndef _EXECUTE_
#define _EXECUTE_

#include "../util/KdTree.h"
#include "../util/Common.h"


#include <string>
#include <vector>


class MeshImp
{
public:
	MeshImp(int numVert, double **Verts, int numTri, int**Tris);
	~MeshImp();

	
	//** Non-obtuse Remeshing **//
	void NonobtuseRemeshing(int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double minAngleAllow, bool verbose);
	void NonobtuseRemeshing_InterleaveOpt(int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double minAngleAllow, bool verbose);
	void AcuteRemeshing_InterleaveOpt(int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double minAngleAllow, double maxAngleAllow, bool verbose);
	void SmallAngleElimination_InterleaveOpt(double targetMinAngle, int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double maxAngleAllow, bool verbose);


	//** Sifting/Simplification **//
	void Simp(int targetNumSamples, int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double minAngleAllow, double maxAngleAllow, bool verbose);
	

	//** Postprocessing **//
	void WriteStatToFile(std::string, int, vert*, bool obt);
	void DisplayStat(int, vert*, bool obt);
	void GetMeshOBJ(std::string filename, bool obtuse, std::string filename_obt);
	

	void PerturbVertices();

private:
	
	//Classes 
	KdTree surfaceKdTree; // the input surface kd tree 
	BoundBox myBoundingBox;
	
	//Functions 
	void ErrWarnMessage(size_t, std::string, size_t);
	void FindNeighbourTriangles(int**&);
	bool ObtuseHead(int, int&, int&);
	bool AcuteHead(int ip, int&ip1, int&ip2, double measureAngle);
	bool TooSmallHeadAngle(int ip, int&ip1, int&ip2, double measureAngle);
	void InitialSortAroundVertex(int**);
	void ScaleInputMesh();
	void GetFanTriangle();

	//Variables 
	size_t MaxNumVert;
	int numVert_org, numTri_org, numVert_imp;
	vert *Vert_org;//the original surface 	
	tri *Tri_org;//the original surface triangles 
	vert *Vert_imp;//the improved/modified surface 
	double scale_factor;
	std::vector<int>indices; 
		
	


};



#endif /*_EXECUTE_*/