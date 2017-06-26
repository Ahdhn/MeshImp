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

#include "Statistics.h"
#include "Common.h"

#include <float.h>
#include <iostream>
#include <string.h>
#include <algorithm>

//TODO report non-delaunay triangle 
Statistics::Statistics()
{	
	min_dih_ang = DBL_MAX;
	max_dih_ang = DBL_MIN; 
	min_quality = DBL_MAX;
	max_quality = DBL_MIN;
	min_angle = DBL_MAX;
	max_angle = DBL_MIN; 
	min_edge_len = DBL_MAX;
	max_edge_len = DBL_MIN;
	num_tri = 0; 
	num_tri_obtuse = 0;
}

Statistics::~Statistics()
{
}

void Statistics::DisplayStatistics(FILE *fp, bool obt)
{
	fprintf(fp, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(fp, "********************************** Statistics ***********************************");
	fprintf(fp, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

	fprintf(fp, "\n #Triangle = %i", num_tri);
	fprintf(fp, "\n #Vertex = %i", num_vertices);
	fprintf(fp, "\n");

	if (obt){
		fprintf(fp, "\n #Obtuse triangle = %i", num_tri_obtuse);
		fprintf(fp, "\n Obtuse triangle = %f %%", 100.0*double(num_tri_obtuse) / double(num_tri));
		fprintf(fp, "\n");
	}

	fprintf(fp, "\n Minimum angle = %f", min_angle);
	fprintf(fp, "\n Maximum angle = %f", max_angle);
	fprintf(fp, "\n Average deviation from 60 deg = %f", average_div_from_60);
	fprintf(fp, "\n");

	fprintf(fp, "\n Minimum edge length = %f", sqrt(min_edge_len));
	fprintf(fp, "\n Maximum edge length = %f", sqrt(max_edge_len));
	fprintf(fp, "\n Minimum to Maximum edge length = %f  %%", 100.0*sqrt(min_edge_len / max_edge_len)); 
	fprintf(fp, "\n");
	
	fprintf(fp, "\n Minimum triangle quality = %f", min_quality);
	fprintf(fp, "\n Maximum triangle quality = %f", max_quality);
	fprintf(fp, "\n");

	fprintf(fp, "\n Minimum dihedral angle = %f", min_dih_ang);
	fprintf(fp, "\n Maximum dihedral angle = %f", max_dih_ang);
	fprintf(fp, "\n");	
}

void Statistics::GetAllStats(int numVert, vert *Vert)
{
	
	CalcDihAngle(numVert, Vert);	
	CalcAnglesEdgeLen(numVert, Vert);
	CalcTriangelQuality(numVert, Vert);
	CalcNumTriangles(numVert, Vert);
	CalcNumNonObtuse(numVert, Vert);
}

//** Num triangles **//
void Statistics::CalcNumTriangles(int numVert, vert *Vert)
{
	num_tri = 0;
	num_vertices = numVert;
	for (int id0 = 0; id0 < numVert; id0++){
		int id1 = Vert[id0].connect[Vert[id0].connect[0]];
		for (int i = 1; i <= Vert[id0].connect[0]; i++){
			int id2 = Vert[id0].connect[i];
			num_tri++;
		}
	}

	if (num_tri % 3 != 0){
		ErrWarnMessage(__LINE__, "Statistics::CalcNumTriangles num_tri%3!=0", 1);
	}

	num_tri /= 3;
}

//** Angle and edge length**//
void Statistics::CalcAnglesEdgeLen(int numVert, vert *Vert)
{
	min_angle = 180.0;
	max_angle = 0.0;
	min_edge_len = DBL_MAX;
	max_edge_len = DBL_MIN;
	num_vertices = numVert;
	average_div_from_60 = 0;
	int num_angles = 0;
	for (int id0 = 0; id0 < numVert; id0++){ 

		int id1 = Vert[id0].connect[Vert[id0].connect[0]];

		for (int i = 1; i <= Vert[id0].connect[0]; i++){
			int id2 = Vert[id0].connect[i];
						
			double angle = AngleVectVect(Vert[id2].x[0] - Vert[id0].x[0], Vert[id2].x[1] - Vert[id0].x[1], Vert[id2].x[2] - Vert[id0].x[2],
				                         Vert[id1].x[0] - Vert[id0].x[0], Vert[id1].x[1] - Vert[id0].x[1], Vert[id1].x[2] - Vert[id0].x[2])*RadToDeg;//180.0 / PI;
			
			double len = Dist(Vert[id0].x[0], Vert[id0].x[1], Vert[id0].x[2], Vert[id1].x[0], Vert[id1].x[1], Vert[id1].x[2]);

			min_angle = std::min(min_angle, angle);
			max_angle = std::max(max_angle, angle);
			min_edge_len = std::min(min_edge_len, len); 
			max_edge_len = std::max(max_edge_len, len); 

			average_div_from_60 += abs(angle - 60.0);
			num_angles++;

			id1 = id2;
		}
	}

	average_div_from_60 /= double(num_angles);

}

//** Non obtuse triangles**//
void Statistics::CalcNumNonObtuse(int numVert, vert *Vert)
{
	num_tri_obtuse = 0;
	num_vertices = numVert;
	for (int id0 = 0; id0 < numVert; id0++){

		int id1 = Vert[id0].connect[Vert[id0].connect[0]];

		for (int i = 1; i <= Vert[id0].connect[0]; i++){
			int id2 = Vert[id0].connect[i];
						 
			double angle = AngleVectVect(Vert[id2].x[0] - Vert[id0].x[0], Vert[id2].x[1] - Vert[id0].x[1], Vert[id2].x[2] - Vert[id0].x[2],
				                         Vert[id1].x[0] - Vert[id0].x[0], Vert[id1].x[1] - Vert[id0].x[1], Vert[id1].x[2] - Vert[id0].x[2])*RadToDeg; /* 57.295779513078550 = 180.0 / PI*/
			if (angle > 90.0 + _tol){ num_tri_obtuse++; }

			id1 = id2;
		}
	}
}

//** Acute triangles**//
int Statistics::CalcNumAcute(int numVert, vert *Vert, double measureAngle)
{
	int num_acute_angle = 0;	
	for (int id0 = 0; id0 < numVert; id0++){

		int id1 = Vert[id0].connect[Vert[id0].connect[0]];

		for (int i = 1; i <= Vert[id0].connect[0]; i++){
			int id2 = Vert[id0].connect[i];

			double angle = AngleVectVect(Vert[id2].x[0] - Vert[id0].x[0], Vert[id2].x[1] - Vert[id0].x[1], Vert[id2].x[2] - Vert[id0].x[2],
				                         Vert[id1].x[0] - Vert[id0].x[0], Vert[id1].x[1] - Vert[id0].x[1], Vert[id1].x[2] - Vert[id0].x[2])*RadToDeg;
			if (angle > measureAngle + _tol){ num_acute_angle++; }
			id1 = id2;
		}
	}
	return num_acute_angle;
}

//** Dih angle**//
double Statistics::getCurv(int ip, int numVert, vert *Vert)
{
	//curv is the supplementary angle of the dihderal angle between two triangles 
	//Here for vertex ip, we return the average curv between each two triangle in ip triangle fan

	double curve(0), angle; 

	int ip3 = Vert[ip].connect[Vert[ip].connect[0]]; 

	for (int i = 1; i <= Vert[ip].connect[0]; i++){
		int ip1 = Vert[ip].connect[i];
		int ip2 = (i == Vert[ip].connect[0]) ? Vert[ip].connect[1] : Vert[ip].connect[i + 1]; 
				
		angle = TriTriNormalAngle(Vert[ip].x, Vert[ip1].x, Vert[ip2].x, Vert[ip3].x);

		//curve = std::max(curve, angle);
		curve += angle;
		ip3 = ip1;
	}
	return curve / double(Vert[ip].connect[0]); 
}
void Statistics::CalcDihAngle(int numVert, vert *Vert)
{
	//Calc and store the dih angle at each vertex
	//the dih angle of a vertex is the average dih angle of the vertex triangle fan

	min_dih_ang = DBL_MAX;
	max_dih_ang = DBL_MIN;
	num_vertices = numVert;

	for (int i = 0; i < numVert; i++){
		Vert[i].dih_ang = 180.0 - getCurv(i, numVert, Vert);
		min_dih_ang = std::min(min_dih_ang, Vert[i].dih_ang);
		max_dih_ang = std::max(max_dih_ang, Vert[i].dih_ang);		
	}
}

//** Triangle quality**//
double Statistics::SingleTriangleQuality(int ip, int ip1, int ip2, vert *Vert)
{
	
	double l1, l2, l3, longest, half_perimeter, area, q;
	l1 = sqrt(Dist(Vert[ip].x[0], Vert[ip].x[1], Vert[ip].x[2], Vert[ip1].x[0], Vert[ip1].x[1], Vert[ip1].x[2]));
	l2 = sqrt(Dist(Vert[ip].x[0], Vert[ip].x[1], Vert[ip].x[2], Vert[ip2].x[0], Vert[ip2].x[1], Vert[ip2].x[2]));
	l3 = sqrt(Dist(Vert[ip2].x[0], Vert[ip2].x[1], Vert[ip2].x[2], Vert[ip1].x[0], Vert[ip1].x[1], Vert[ip1].x[2]));

	longest = std::max(l1, l2);
	longest = std::max(longest, l3);

	half_perimeter = 0.5*(l1 + l2 + l3);

	area = sqrt(half_perimeter*(half_perimeter - l1)*(half_perimeter - l2)*(half_perimeter - l3));
	q = (3.4641016151377548 * area) / (half_perimeter*longest); //3.4641016151377548 = 6.0/sqrt(3.0)
	//q = (6.0*area) / (sqrt(3.0)*half_perimeter*longest);

	return q;


}
void Statistics::CalcTriangelQuality(int numVert, vert *Vert)
{
	min_quality = DBL_MAX;
	max_quality = DBL_MIN;
	num_vertices = numVert;
	for (int ip = 0; ip < numVert; ip++){
		int ip2 = Vert[ip].connect[Vert[ip].connect[0]];
		for (int i = 1; i <= Vert[ip].connect[0]; i++){
			int ip1 = Vert[ip].connect[i];
			if (ip < ip1 && ip < ip2){//just do the triangle once 
				double qu = SingleTriangleQuality(ip, ip1, ip2, Vert);
				min_quality = std::min(min_quality, qu);
				max_quality = std::max(max_quality, qu);
			}

			ip2 = ip1;
		}
	}
}

void Statistics::ErrWarnMessage(size_t lineNum, std::string message, size_t mess_id)
{
	
	if (mess_id == 0){
		std::cerr << "\nError::line(" << lineNum << ")-->>" << message << std::endl;
		system("pause");
	}
	else{
		std::cerr << "\nWarning::line(" << lineNum << ")-->>" << message << std::endl;
		system("pause");
	}
}