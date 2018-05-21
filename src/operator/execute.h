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
#include "../util/Sphere.h"

#include <string>
#include <vector>
#include <Eigen/Dense>

class MeshImp
{
public:
	MeshImp(){};
	~MeshImp(){};
	MeshImp(int numSpheres, double **Spheres, int numVert, double**Verts, int numTri, int**Tris);	
	

	//** Sifting **//
	void VC(int targetNumSamples, int samplingBudget, int numSurfaceLayer, bool isSmooth, bool verbose);
	

private:
	
	
	//Classes 	
	Sphere mSphere;
		
	//Functions 
	void ErrWarnMessage(size_t, std::string, size_t);
	void FindNeighbourTriangles();
	
	
	//Variables 
	Tri mTriangles;

	
};



#endif /*_EXECUTE_*/