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

#include "../util/Common.h"
#include "../util/RNG.h"
#include "../util/Sphere.h"

class Operator
{
public:
	
	Operator();	
	~Operator();
				
	bool EjectionSpheres(int*ip, Tri*Triangles, Sphere*Spheres, int samplingBudget, int numSurfaceLayer, bool isSmooth, bool verbose);
	

private:
	/*****************************************/
	int m_neighbor_tri[2000];
	double* m_two_planes1 = new double [6];	
	double* m_two_planes2 = new double[6];

	std::vector<int> m_overlap;
	std::vector<int> m_overlap_two_layers;

	inline bool InsectPointBetweenTwoSamples(Sphere* Spheres, int p1, int p2, double&x_in1, double&y_in1, double&z_in1, double&x_in2, double&y_in2, double&z_in2, Tri* Triangles, size_t&sect_num, bool get_neighbors, const int _num_layers);
	inline void GetNeighbors(Tri* Triangles, int tri, int num_layer);
	inline void GetTriNeighbors(Tri*Triangles, int tri_num);
	
	double m_insect_constraints[1000][3];
	int m_num_insect_constraints;

	int m_seeds_constraints[1000];


	int num_active;
	int**active_pool, **tmp_active_pool, **tri_pool;
	double **_tar;
	double newVertex[3];//coordinates of the new vertex (updated from Sampler)
	int newVertexTri;
	void StartActivePool(Tri* Triangles);
	bool Sampler(Sphere* Spheres, Tri*Triangles, double r_2new, int*ip);
	void RetrieveCoordinates(Tri*Triangles, int lf, int* cell, double&x0, double&y0, double&z0, double&x1, double&y1, double&z1, double&x2, double&y2, double&z2);
	void RandomPointInTri(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double&x, double&y, double&z);
	bool CheckNewVertex(Sphere* Spheres, double x_new, double  y_new, double  z_new, double r_2new, int*ip);
	bool CheckRefinement(Sphere *Spheres,double r_2new, double*tri);
	void ErrWarnMessage(size_t lineNum, std::string message, size_t mess_id);
	void DrawActivePool(Tri*Triangles, int lf);
	void DrawnList(Sphere* Spheres, std::vector<int>overlap);
	void DrawnList(Sphere* Spheres, int*ip);
	void Draw_m_neighbor_tri(Tri*Triangles);
	void DrawTriangle(Tri*Triangles, int id);
	inline bool isCoveredby(double xx, double yy, double zz, std::vector<int> overlap, Sphere* Spheres, int skip1, int skip2);
	void DrawSeedList(Sphere *Spheres);
	RndNum myRandNum;	
};



#endif /*_OPERATOR_*/