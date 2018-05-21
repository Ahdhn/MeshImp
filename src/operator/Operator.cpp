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
#include "../util/Sphere.h"
//#include "../util/Draw_header.h"

//#include "DrawSpheresDebug.h" //for debugging 

#define DEBUGGING
extern bool isObtuse;
extern bool isAcute;

double OpimizerFunc_CenterAngle(vert*Verts, int*nList, double*void_verx, double xx, double yy, double zz);
double OpimizerFunc_SideAngle(vert*Verts, int*nList, double minAngle, double maxAngle, double xx, double yy, double zz);
double OpimizerFunc_Closer(vert*Verts, int*nList, double*void_verx, double xx, double yy, double zz);
double OpimizerFunc_Further(vert*Verts, int*nList, double*void_verx, double xx, double yy, double zz);

Operator::Operator()                    
{
	active_pool = new int *[MAXPOOL];
	tri_pool = new int *[MAXPOOL];
	tmp_active_pool = new int*[MAXPOOL];
	for (int i = 0; i < MAXPOOL; i++){
		active_pool[i] = new int[3];
		tri_pool[i] = new int[4];
		tmp_active_pool[i] = new int[3];
	}

	_tar = new double*[4];
	_tar[0] = new double[9];
	_tar[1] = new double[9];
	_tar[2] = new double[9];
	_tar[3] = new double[9];	

	//srand(time(NULL));	
	myRandNum.InitiateRndNum(rand());
}

Operator::~Operator()
{

}

///**************Operators
bool Operator::EjectionSpheres(int*ip, Tri* Triangles, Sphere* Spheres, int samplingBudget, int numSurfaceLayer, bool isSmooth, bool verbose)
{
	//ip = pointer to spheres to eject
	//Spheres = shperes (x,y,z,r)
#ifdef DEBUGGING
	if (ip[0] != 2 && ip[0] != 3){
		ErrWarnMessage(__LINE__, "Operator::EjectionSpheres:: can only ejection two or three vertices", 0);
	}
#endif
		
	double xp1, yp1, zp1, r_2p1, xp2, yp2, zp2, r_2p2;
	int tri1(0), tri2(0);
	Spheres->getCoordinates(ip[1], xp1, yp1, zp1, r_2p1, tri1);
	Spheres->getCoordinates(ip[2], xp2, yp2, zp2, r_2p2, tri2);


	std::vector<int> overlap1;
	Spheres->getOverlap(ip[1], overlap1);
	removeFromVector(ip[2], overlap1);
	std::vector<int> overlap2;
	Spheres->getOverlap(ip[2], overlap2);
	removeFromVector(ip[1], overlap2);
	//the overlap should only pick the unique elements
	m_overlap.clear();
	m_overlap = overlap1;
	//std::vector<int> overlap(overlap1);
	for (int l = 0; l < int(overlap2.size()); l++){
		if (std::find(overlap1.begin(), overlap1.end(), overlap2[l]) == overlap1.end()){
			m_overlap.push_back(overlap2[l]);
		}
	}
	

	if (false){
		DrawnList(Spheres, m_overlap);
		Spheres->Draw_m_sphere(ip[1]);
		Spheres->Draw_m_sphere(ip[2]);
	}

	//radius of the new sample is average of r1 and r2
	double r_2_average = 0.25*(r_2p1 + r_2p2 + 2.0*sqrt(r_2p1*r_2p1));
	

	//1st constraint -> maximality 
	//This is the inSpheres constraints where an inSpheres is centered
	//on the intersection point and has a radius of r_2_average
	m_num_insect_constraints = 0;

	for (int i = 0; i < int(m_overlap.size()) - 1; i++){
		
		if (m_overlap[i] == ip[1] || m_overlap[i] == ip[2]){ continue; }

		for (int j = i + 1; j < int(m_overlap.size()); j++){

			if (m_overlap[j] == ip[1] || m_overlap[j] == ip[2]){ continue; }

			double x_in1, y_in1, z_in1, x_in2, y_in2, z_in2;
			size_t sect_num = 0;
			if (InsectPointBetweenTwoSamples(Spheres, m_overlap.at(i), m_overlap.at(j),
				                             x_in1, y_in1, z_in1, x_in2, y_in2, z_in2,
				                             Triangles, sect_num, true, numSurfaceLayer)){
			
				
				bool in1 = (Dist(x_in1, y_in1, z_in1, xp1, yp1, zp1) < r_2p1) ||
					       (Dist(x_in1, y_in1, z_in1, xp2, yp2, zp2) < r_2p2);

				if (in1 && !isCoveredby(x_in1, y_in1, z_in1, m_overlap, Spheres, m_overlap[i], m_overlap[j])){
					//finally make sure it covered by ip[1] and only ip[1]

					m_insect_constraints[m_num_insect_constraints][0] = x_in1;
					m_insect_constraints[m_num_insect_constraints][1] = y_in1;
					m_insect_constraints[m_num_insect_constraints][2] = z_in1;
					m_num_insect_constraints++;

				}
				if (sect_num == 2){
					bool in2 = (Dist(x_in2, y_in2, z_in2, xp1, yp1, zp1) < r_2p1) ||
						       (Dist(x_in2, y_in2, z_in2, xp2, yp2, zp2) < r_2p2);
					if (in2 && !isCoveredby(x_in2, y_in2, z_in2, m_overlap, Spheres, m_overlap[i], m_overlap[j])){
						m_insect_constraints[m_num_insect_constraints][0] = x_in2;
						m_insect_constraints[m_num_insect_constraints][1] = y_in2;
						m_insect_constraints[m_num_insect_constraints][2] = z_in2;
						m_num_insect_constraints++;
					}
				}
			}
			
		}
	}
	if (m_num_insect_constraints == 0){
		//this is not good since this is the only inclusion region we have, 
		//having it empty means the new sampled is not anchored to anything and could be 
		//anywhere, so we quit
		return false;
		//std::cout << "Error  Operator::EjectionSpheres:: can detect any intersection with surface constraints!!!" << std::endl;
	}
	if (m_num_insect_constraints < 3){
		std::cout << "Error  Operator::EjectionSpheres:: detected less than three intersection with surface!!!" << std::endl;
		return false;
	}

	//check if you can cover all the seed with one sphere of radius r_2_average
	//i.e., is it feasible region
	if (true){
		for (int i = 0; i < m_num_insect_constraints; i++){
			for (int j = i + 1; j < m_num_insect_constraints; j++){
				if (Dist(m_insect_constraints[i][0], m_insect_constraints[i][1], m_insect_constraints[i][2],
					m_insect_constraints[j][0], m_insect_constraints[j][1], m_insect_constraints[j][2])>4 * r_2_average){
					return false;
				}
			}
		}
	}
	else{
		r_2_average = 0;
		for (int i = 0; i < m_num_insect_constraints; i++){
			for (int j = i + 1; j < m_num_insect_constraints; j++){
				r_2_average = std::max(r_2_average, Dist(m_insect_constraints[i][0], m_insect_constraints[i][1], m_insect_constraints[i][2],
					m_insect_constraints[j][0], m_insect_constraints[j][1], m_insect_constraints[j][2]));
			}
		}
		r_2_average *= 0.25;
		r_2_average += 2.0*_tol;

	}

	//2nd constraint -> don't cover exisiting intersection pairs 
	//these are the intersection pairs from overlapping spheres except
	//those that live over ip
	//This represents the outSphere constraints where an outSphere is represented 
	//by a seed of radius r_2_average
	m_seeds_constraints[0] = 0;
	std::vector<int>seeds;
	std::vector<int>ip1_seeds;
	std::vector<int>ip2_seeds;
	Spheres->getSeeds(ip[1], ip1_seeds);
	Spheres->getSeeds(ip[2], ip2_seeds);
	for (int i = 0; i < int(m_overlap.size()); i++){
		Spheres->getSeeds(m_overlap.at(i), seeds);
				
		for (int s = 0; s <int(seeds.size()); s++){

			if (std::find(ip1_seeds.begin(), ip1_seeds.end(), seeds[s]) != ip1_seeds.end() ||
				std::find(ip2_seeds.begin(), ip2_seeds.end(), seeds[s]) != ip2_seeds.end() ){
				//if it lives on ip1 or ip2, don't add it to the constraints 
				continue;
			}			

			bool already_taken = false;
			for (int c = 1; c <=m_seeds_constraints[0]; c++){
				if (seeds[s] == m_seeds_constraints[c]){
					already_taken = true;
					break;
				}
			}
			if (!already_taken){
				m_seeds_constraints[++m_seeds_constraints[0]] = seeds[s];
			}
		}
	}
	if (false){
		DrawSeedList(Spheres);
	}


	//we used m_overlap_two_layers for isSliverCandidate because one layer of overlapping is not 
	//enough to detect all overlap 
	m_overlap_two_layers.clear();
	m_overlap_two_layers = m_overlap;
	for (int i = 0; i < int( m_overlap.size()); i++){
		int sp = m_overlap[i];
		std::vector<int> sp_overlap;
		Spheres->getOverlap(sp, sp_overlap);
		removeFromVector(ip[1], sp_overlap);
		removeFromVector(ip[2], sp_overlap);
		for (int j = 0; j < int(sp_overlap.size()); j++){
			if (std::find(m_overlap_two_layers.begin(), m_overlap_two_layers.end(), sp_overlap[j]) == m_overlap_two_layers.end()){
				m_overlap_two_layers.push_back(sp_overlap[j]);
			}
		}
	}
	if (false){
		DrawnList(Spheres, m_overlap_two_layers);
	}
	//sampling 
	//don't need to initialize active pool since we initilize it already during 
	//computing the intersection points 
	StartActivePool(Triangles);
	if (Sampler(Spheres, Triangles, r_2_average, ip)){		
		Spheres->removeTwoAddOne(ip[1], ip[2], newVertex[0], newVertex[1], newVertex[2], r_2_average, newVertexTri);		

		std::vector<int>::iterator it = std::find(m_overlap_two_layers.begin(), m_overlap_two_layers.end(), Spheres->getNumSphere());
		if (it != m_overlap_two_layers.end()){
			//if the last sphere is in the list, change it to ip[2] since in removeTwoAddOne
			//we move the last sphere to be ip[2]
			int index = it - m_overlap_two_layers.begin();
			m_overlap_two_layers[index]= ip[2];
		}

		Spheres->re_setSeeds(m_overlap_two_layers);


		//Spheres->Draw_m_sphere(ip[1]);
		//Spheres->vcSurface();
		if (verbose){
			std::cout << "	Replaced <" << ip[1] << "," << ip[2] << "> by one" << std::endl;
		}
		if (false){
			Spheres->printSpheres_csv("out_spheres.csv");
			Spheres->printSeeds_csv("out_seeds.csv");
		}
		return true;
	}
	return false;	
}
inline bool Operator::isCoveredby(double xx, double yy, double zz, std::vector<int> overlap, Sphere* Spheres, int skip1, int skip2)
{
	for (int i = 0; i <int( overlap.size()); i++){
		if (overlap[i] == skip1 || overlap[i] == skip2){ continue; }

		double x1, y1, z1, r_2;
		int tri1;
		Spheres->getCoordinates(overlap[i], x1, y1, z1, r_2, tri1);

		double d = Dist(xx, yy, zz, x1, y1, z1);
		if (d < r_2 - _tol_sq){
			return true;
		}
	}

	return false;
}
inline void Operator::GetTriNeighbors(Tri*Triangles, int tri_num)
{
	int t1, t2, t3, V, tt;
	
	t1 = Triangles->neighbour[tri_num][0];
	t2 = Triangles->neighbour[tri_num][1];
	t3 = Triangles->neighbour[tri_num][2];
	bool in = true;

	for (V = 1; V <= m_neighbor_tri[0]; V++){
		tt = m_neighbor_tri[V];

		if (t1 == tt){ in = false; break; }
	}
	if (in){
		m_neighbor_tri[0]++;
		m_neighbor_tri[m_neighbor_tri[0]] = t1;
	}
	in = true;


	for (V = 1; V <= m_neighbor_tri[0]; V++){
		tt = m_neighbor_tri[V];

		if (t2 == tt){ in = false; break; }
	}
	if (in){
		m_neighbor_tri[0]++;
		m_neighbor_tri[m_neighbor_tri[0]] = t2;
	}
	in = true;



	for (V = 1; V <= m_neighbor_tri[0]; V++){
		tt = m_neighbor_tri[V];

		if (t3 == tt){ in = false; break; }
	}
	if (in){
		m_neighbor_tri[0]++;
		m_neighbor_tri[m_neighbor_tri[0]] = t3;
	}
	in = true;

}
inline void Operator::GetNeighbors(Tri* Triangles, int tri, int num_layer)
{
	// getting triangle neighbor for _num_layer layers of tri	
	bool in = true;
	for (int V = 1; V <= m_neighbor_tri[0]; V++){
		if (m_neighbor_tri[V] == tri){ in = false; break; }
	}
	if (in){
		m_neighbor_tri[0]++;
		m_neighbor_tri[m_neighbor_tri[0]] = tri;
	}

	GetTriNeighbors(Triangles, tri);
	
	int counter(0), start(1), num, d;

	while (counter<num_layer){
		num = m_neighbor_tri[0];
		for (d = start; d <= num; d++){
			GetTriNeighbors(Triangles,m_neighbor_tri[d]);
		}
		start = num;
		counter++;
	}
}
inline bool Operator::InsectPointBetweenTwoSamples(Sphere* Spheres,
										 int p1, int p2,
	                                     double&x_in1, double&y_in1, double&z_in1,
	                                     double&x_in2, double&y_in2, double&z_in2,
										 Tri* Triangles,	                                     
										 size_t&sect_num, bool get_neighbors, const int _num_layers)
{
	//get the intersection pair between spheres p1 and p2, and the surface 
	//sphere p1 and p2 both have access to the triangle they are attached to 
	//from which (if get_neighbors == true) we can get the neighbour triangles 
	//on which the interesection pair will lie 
	
	size_t found(0);
	/*size_t V, node1, node2, node3, tri, found(0);
	double x_c, y_c, z_c, r_c, x_v, y_v, z_v, len, dx, dy, dz, r_c_2, x_v1, y_v1, z_v1,
		line_x, line_y, line_z, line_x_v, line_y_v, line_z_v, u1, u2, x1, y1, z1, x2, y2, z2, alfa, beta, gamma, l1,
		xl, yl, zl;*/
	
	bool exist = false;
	
	double xp1, yp1, zp1, r_2p1, xp2, yp2, zp2, r_2p2;
	int tri1(0), tri2(0);
	Spheres->getCoordinates(p1, xp1, yp1, zp1, r_2p1, tri1);
	Spheres->getCoordinates(p2, xp2, yp2, zp2, r_2p2, tri2);


		
	double len_2 = Dist(xp1,yp1,zp1,xp2,yp2,zp2);
	if (len_2 + _tol_sq > r_2p1 + r_2p2 + 2.0*sqrt(r_2p1 * r_2p2)){
		return false;
	}
	double len = sqrt(len_2);

	double x_v = xp2 - xp1; // vector perpendular to the plane of the intersection
	double y_v = yp2 - yp1;
	double z_v = zp2 - zp1;
	double l1 = (r_2p1 - r_2p2 + len_2) / (2.0*len);
	double x_c = x_v*l1 / len + xp1;
	double y_c = y_v*l1 / len + yp1;
	double z_c = z_v*l1 / len + zp1;

	len_2 = Dist(xp1, yp1, zp1, x_c, y_c, z_c);
	
	double r_c_2 = r_2p1 - len_2;
	
	
	len_2 = Dist(xp2, yp2, zp2, x_c, y_c, z_c);
	

	if (abs(r_2p2 - len_2 - r_c_2)>_tol){
		std::cout << "Error (0) at Operator::InsectPointBetweenTwoSamples" << std::endl;
	}


	if (get_neighbors){
		m_neighbor_tri[0] = 0;
		GetNeighbors(Triangles, tri1, _num_layers);
		if (tri1 != tri2){ GetNeighbors(Triangles, tri2, _num_layers); }
	}

	if (false){
		Draw_m_neighbor_tri(Triangles);
		DrawTriangle(Triangles, 0);
	}



	for (int V = 1; V <= m_neighbor_tri[0]; V++){
		int tri = m_neighbor_tri[V];		
		int node1 = Triangles->ids[tri][0];
		int node2 = Triangles->ids[tri][1];
		int node3 = Triangles->ids[tri][2];

		//if(false){
		//	DrawTri2(tri,tri);
		//}
		double xl, yl, zl;		
		bool lo1 = LinePlaneIntersect(x_c, y_c, z_c, x_v, y_v, z_v,
			Triangles->coord[node1][0], Triangles->coord[node1][1], Triangles->coord[node1][2],
			Triangles->coord[node2][0], Triangles->coord[node2][1], Triangles->coord[node2][2],
			xl, yl, zl);
			

		bool lo2 = LinePlaneIntersect(x_c, y_c, z_c, x_v, y_v, z_v,
			Triangles->coord[node2][0], Triangles->coord[node2][1], Triangles->coord[node2][2],
			Triangles->coord[node3][0], Triangles->coord[node3][1], Triangles->coord[node3][2],
			xl, yl, zl);
			

		bool lo3 = LinePlaneIntersect(x_c, y_c, z_c, x_v, y_v, z_v,
			Triangles->coord[node3][0], Triangles->coord[node3][1], Triangles->coord[node3][2],
			Triangles->coord[node1][0], Triangles->coord[node1][1], Triangles->coord[node1][2],
			xl, yl, zl);
			
		if (!lo1  && !lo2 && !lo3){ continue; }


		m_two_planes1[0] = x_c;
		m_two_planes1[1] = y_c;
		m_two_planes1[2] = z_c;
		m_two_planes1[3] = x_v;
		m_two_planes1[4] = y_v;
		m_two_planes1[5] = z_v;

		double x_v1, y_v1, z_v1;
		Cross(Triangles->coord[node1][0] - Triangles->coord[node2][0], //xv1
			  Triangles->coord[node1][1] - Triangles->coord[node2][1], //yv1
			  Triangles->coord[node1][2] - Triangles->coord[node2][2], //zv1
			  Triangles->coord[node1][0] - Triangles->coord[node3][0], //xv2
			  Triangles->coord[node1][1] - Triangles->coord[node3][1], //yv2
			  Triangles->coord[node1][2] - Triangles->coord[node3][2], //zv2
			  x_v1,y_v1,z_v1 );
		
		m_two_planes2[0] = Triangles->coord[node1][0];
		m_two_planes2[1] = Triangles->coord[node1][1];
		m_two_planes2[2] = Triangles->coord[node1][2];
		m_two_planes2[3] = x_v1;
		m_two_planes2[4] = y_v1;
		m_two_planes2[5] = z_v1;
		    
			
			

		//DrawPlane(x_v1,y_v1,z_v1,_surface[node1][0],_surface[node1][1],_surface[node1][2]);
		//DrawPlane(x_v,y_v,z_v,x_c,y_c,z_c);
		double line_x(0), line_y(0), line_z(0), line_x_v(0), line_y_v(0), line_z_v(0);
		double u1(0), u2(0);
		GetInsct2Planes(line_x, line_y, line_z, line_x_v, line_y_v, line_z_v, m_two_planes1, m_two_planes2);


		double x1 = line_x ;
		double y1 = line_y ;
		double z1 = line_z ;

		double x2 = line_x + 1000*line_x_v;
		double y2 = line_y + 1000*line_y_v;
		double z2 = line_z + 1000*line_z_v;
		size_t num_sect(0);
		double scale_factor = 100.0;
		if (!SphereLineIntersection(line_x - scale_factor*line_x_v, line_y - scale_factor*line_y_v, line_z - scale_factor*line_z_v,
							        line_x + scale_factor*line_x_v, line_y + scale_factor*line_y_v, line_z + scale_factor*line_z_v,
			                        xp1, yp1, zp1, r_2p1,
			                        x1, y1, z1, x2, y2, z2,
							        num_sect)){
			continue;
		}

		if (num_sect != 2){
			continue;
			//std::cout <<"Error Operator::InsectPointBetweenTwoSamples:: Sphere line intersection does not two points!!!" << std::endl;
		}
		/*if (!(SolveQuadEqu(line_x_v*line_x_v + line_y_v*line_y_v + line_z_v*line_z_v,
			               2.0*(line_x - x_c)*line_x_v + 2.0*(line_y - y_c)*line_y_v + 2.0*(line_z - z_c)*line_z_v,
			               (line_x - x_c)*(line_x - x_c) + (line_y - y_c)*(line_y - y_c) + (line_z - z_c)*(line_z - z_c) - r_c_2,
			               u1, u2))){
			continue;
		}
		
		double x1 = line_x + u1*line_x_v;
		double y1 = line_y + u1*line_y_v;
		double z1 = line_z + u1*line_z_v;

		double x2 = line_x + u2*line_x_v;
		double y2 = line_y + u2*line_y_v;
		double z2 = line_z + u2*line_z_v;*/

		if (abs(Dist(x1, y1, z1, x_c, y_c, z_c) - r_c_2)>_tol){
			std::cout << "Error (1) at Operator::InsectPointBetweenTwoSamples" << std::endl;
		}
				
		if (abs(Dist(x2, y2, z2, x_c, y_c, z_c) - r_c_2)>_tol){
			std::cout << "Error (2) at Operator::InsectPointBetweenTwoSamples" << std::endl;			
		}


		//if(false){
		//	DrawPoint("1in1.obj",x1,y1,z1);
		//	DrawPoint("2in2.obj",x2,y2,z2);
		//}

		double alfa(0), beta(0), gamma(0);
		Barycentric(Triangles->coord[node1][0], Triangles->coord[node1][1], Triangles->coord[node1][2],
			        Triangles->coord[node2][0], Triangles->coord[node2][1], Triangles->coord[node2][2],
					Triangles->coord[node3][0], Triangles->coord[node3][1], Triangles->coord[node3][2], 
					x1, y1, z1, alfa, beta, gamma);
		if (alfa >= -_tol && beta >= -_tol && gamma >= -_tol && alfa <= 1.0 + _tol && beta <= 1.0 + _tol && gamma <= 1.0 + _tol){
			found++;
			if (found == 1){
				x_in1 = x1;
				y_in1 = y1;
				z_in1 = z1;
				tri1 = tri;
				exist = true;

			}
			else if (found == 2){
				x_in2 = x1;
				y_in2 = y1;
				z_in2 = z1;
				tri2 = tri;
				exist = true;
			}
			else{			
				double dist;
				dist = Dist(x_in2, y_in2, z_in2, x_in1, y_in1, z_in1);
				if (dist<r_2p1*10E-3){
					x_in2 = x1;
					y_in2 = y1;
					z_in2 = z1;
					tri2 = tri;
				}
			}
		}
		Barycentric(Triangles->coord[node1][0], Triangles->coord[node1][1], Triangles->coord[node1][2],
			        Triangles->coord[node2][0], Triangles->coord[node2][1], Triangles->coord[node2][2],
			        Triangles->coord[node3][0], Triangles->coord[node3][1], Triangles->coord[node3][2],
			        x2, y2, z2, alfa, beta, gamma);

		if (alfa >= -_tol && beta >= -_tol && gamma >= -_tol && alfa <= 1.0 + _tol && beta <= 1.0 + _tol && gamma <= 1.0 + _tol){
			found++;
			if (found == 1){
				x_in1 = x2;
				y_in1 = y2;
				z_in1 = z2;
				tri1 = tri;
				exist = true;
			}
			else if (found == 2){
				x_in2 = x2;
				y_in2 = y2;
				z_in2 = z2;
				tri2 = tri;
				exist = true;
			}

			else{				
				double dist;
				dist = Dist(x_in2, y_in2, z_in2, x_in1, y_in1, z_in1);
				if (dist<r_2p1*10E-3){
					x_in2 = x2;
					y_in2 = y2;
					z_in2 = z2;
					tri2 = tri;
				}
			}
		}

		if (found == 2){
			sect_num = 2;
			return exist;
		}
	}
	sect_num = found;
	return exist;
}
void Operator::StartActivePool(Tri* Triangles){
	//your active pool is just m_neighbour_tri
	num_active = 0;
	for (int i = 1; i <= m_neighbor_tri[0]; i++){
		int my_tri = m_neighbor_tri[i];
		tri_pool[num_active][0] = Triangles->ids[my_tri][0];
		tri_pool[num_active][1] = Triangles->ids[my_tri][1];
		tri_pool[num_active][2] = Triangles->ids[my_tri][2];
		tri_pool[num_active][3] = my_tri;

		active_pool[num_active][0] = num_active;
		active_pool[num_active][1] = 0;
		active_pool[num_active][2] = 0;
		num_active++;
		if (num_active + 2 >= MAXPOOL){ return; }
	}
}
bool Operator::Sampler(Sphere* Spheres, Tri*Triangles, double r_2new, int*ip)
{

	//the resampling routine, called after constructing the geometric constraints 
	//and initlize the active pool
	//it should spit out a single new vertex in case of success and return true
	//otherwise return false
	
	for (int lf = 0; lf < 5; lf++){
		int attempt = int(0.8*m_neighbor_tri[0]);

		//a) dart throwing 

		for (int i = 0; i < attempt; i++){
			int tri = int(myRandNum.RandNumGenerator()*double(num_active - 1));
			double x1, y1, z1, x2, y2, z2, x3, y3, z3;

			RetrieveCoordinates(Triangles, lf, active_pool[tri], x1, y1, z1, x2, y2, z2, x3, y3, z3);

			double x_new, y_new, z_new;

			RandomPointInTri(x1, y1, z1, x2, y2, z2, x3, y3, z3, x_new, y_new, z_new);

			if (CheckNewVertex(Spheres, x_new, y_new, z_new, r_2new, ip)){
				newVertex[0] = x_new;
				newVertex[1] = y_new;
				newVertex[2] = z_new;
				newVertexTri = tri_pool[active_pool[tri][0]][3];
				return true;
			}
		}


		//b) refinement 
		int tmp_num_active = 0;
		for (int iactive = 0; iactive < num_active; iactive++){
			double x1, y1, z1, x2, y2, z2, x3, y3, z3;
			RetrieveCoordinates(Triangles, lf, active_pool[iactive], x1, y1, z1, x2, y2, z2, x3, y3, z3);

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


				if (CheckRefinement(Spheres, r_2new,_tar[C])){

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
						return false;

					}
				}
			}
		}

		//swap
		if (tmp_num_active == 0){			
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

	return false;

}
bool Operator::CheckNewVertex(Sphere *Spheres, double x_new, double  y_new, double  z_new, double r_2new, int*ip)
{
	double xx, yy, zz, r_2;
	int tri;
	//contains all intersection pairs
	for (int i = 0; i < m_num_insect_constraints; i++){
		double dist = Dist(x_new, y_new, z_new, m_insect_constraints[i][0], m_insect_constraints[i][1], m_insect_constraints[i][2]);
		if (dist > r_2new - _tol_sq){ return false; }
	}

	//does not cover any existing seeds
	for (int i = 1; i <= m_seeds_constraints[0]; i++){
		Spheres->getSeedCoordinates(m_seeds_constraints[i], xx, yy, zz);
		double dist = Dist(xx, yy, zz, x_new, y_new, z_new);
		if (dist < r_2new + _tol_sq){ return false; }
	}
		
	//the center is not constained in another sphere
	//and it does not contain other spheres centers 
	for (int i = 0; i < int(m_overlap.size()); i++){
		Spheres->getCoordinates(m_overlap[i], xx, yy, zz,r_2,tri);
		double dist = Dist(xx, yy, zz, x_new, y_new, z_new);
		if (dist < r_2new + _tol_sq || dist < r_2 + _tol_sq){ return false; }
	}


	//won't create a sliver 
	if (Spheres->isSliverCandidate(x_new, y_new, z_new, r_2new, ip[1], ip[2], m_overlap_two_layers)){
		return false;
	}

	return true;
}
void Operator::RetrieveCoordinates(Tri*Triangles, int lf, int* cell, double&x0, double&y0, double&z0, double&x1, double&y1, double&z1, double&x2, double&y2, double&z2)
{
	int po = tri_pool[cell[0]][0];
	int p1 = tri_pool[cell[0]][1];
	int p2 = tri_pool[cell[0]][2];	

	x0 = Triangles->coord[po][0];
	y0 = Triangles->coord[po][1];
	z0 = Triangles->coord[po][2];
	x1 = Triangles->coord[p1][0];
	y1 = Triangles->coord[p1][1];
	z1 = Triangles->coord[p1][2];
	x2 = Triangles->coord[p2][0];
	y2 = Triangles->coord[p2][1];
	z2 = Triangles->coord[p2][2];

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
void Operator::RandomPointInTri(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double&x, double&y, double&z)
{
	double r1, r2;

	r1 = myRandNum.RandNumGenerator();
	r2 = myRandNum.RandNumGenerator();

	x = (1 - sqrt(r1))*x1 + (sqrt(r1)*(1 - r2))*x2 + (r2*sqrt(r1))*x3;
	y = (1 - sqrt(r1))*y1 + (sqrt(r1)*(1 - r2))*y2 + (r2*sqrt(r1))*y3;
	z = (1 - sqrt(r1))*z1 + (sqrt(r1)*(1 - r2))*z2 + (r2*sqrt(r1))*z3;

}
bool Operator::CheckRefinement(Sphere *Spheres,double r_2new, double*tri)
{	
	double xx, yy, zz, r_2;
	int t;
	for (int i = 1; i <= m_seeds_constraints[0]; i++){
		Spheres->getSeedCoordinates(m_seeds_constraints[i], xx, yy, zz);
		double dist1 = Dist(tri[0], tri[1], tri[2], xx, yy, zz);
		double dist2 = Dist(tri[3], tri[4], tri[5], xx, yy, zz);
		double dist3 = Dist(tri[6], tri[7], tri[8], xx, yy, zz);

		if (dist1<r_2new * _tol_sq_circ && dist2<r_2new * _tol_sq_circ && dist3<r_2new * _tol_sq_circ){
			//the whole triangle is in one of the exclusion sphere
			return false;
		}
	}

	
	for (int i = 0; i < int(m_overlap.size()); i++){
		Spheres->getCoordinates(m_overlap[i], xx, yy, zz, r_2, t);
		double dist1 = Dist(tri[0], tri[1], tri[2], xx, yy, zz);
		double dist2 = Dist(tri[3], tri[4], tri[5], xx, yy, zz);
		double dist3 = Dist(tri[6], tri[7], tri[8], xx, yy, zz);
		double r_2_compare = std::max(r_2, r_2new);
		if (dist1 < r_2_compare * _tol_sq_circ && dist2 < r_2_compare * _tol_sq_circ && dist3 < r_2_compare * _tol_sq_circ) {
			//the triangle is completely contained in another sphere 
			return false;
		}
	}

	//check if there is a part of the triangle inside one of the inclusion spheres 
	for (int i = 0; i < m_num_insect_constraints; i++){
		//for each inclusion sphere
		//first check if at least one vertex of the triangle is inside 
		double dist1 = Dist(tri[0], tri[1], tri[2], m_insect_constraints[i][0], m_insect_constraints[i][1], m_insect_constraints[i][2]);
		double dist2 = Dist(tri[3], tri[4], tri[5], m_insect_constraints[i][0], m_insect_constraints[i][1], m_insect_constraints[i][2]);
		double dist3 = Dist(tri[6], tri[7], tri[8], m_insect_constraints[i][0], m_insect_constraints[i][1], m_insect_constraints[i][2]);

		if (dist1 > r_2new * _tol_sq_circ && dist2 > r_2new * _tol_sq_circ && dist3 > r_2new * _tol_sq_circ){	
			//if not, then check if at least one of the triangles three line sgement cross the inclusion spheres
			double px1(0), py1(0), pz1(0), px2(0), py2(0), pz2(0);
			size_t num_sect(0);			
			bool l1 = SphereLineIntersection(tri[0], tri[1], tri[2], tri[3], tri[4], tri[5],
				m_insect_constraints[i][0], m_insect_constraints[i][1], m_insect_constraints[i][2], r_2new,
				px1, py1, pz1, px2, py2, pz2, num_sect);
			bool l2 = SphereLineIntersection(tri[6], tri[7], tri[8], tri[3], tri[4], tri[5],
				m_insect_constraints[i][0], m_insect_constraints[i][1], m_insect_constraints[i][2], r_2new,
				px1, py1, pz1, px2, py2, pz2, num_sect);
			bool l3 = SphereLineIntersection(tri[0], tri[1], tri[2], tri[6], tri[7], tri[8],
				m_insect_constraints[i][0], m_insect_constraints[i][1], m_insect_constraints[i][2], r_2new,
				px1, py1, pz1, px2, py2, pz2, num_sect);
			if (!l1 || !l2 || !l3){
				//if not, check if the inclusion sphere lie on the triangles 
				double alfa(0), beta(0), gamma(0);
				Barycentric(tri[0], tri[1], tri[2],
					        tri[3], tri[4], tri[5],
					        tri[6], tri[7], tri[8],
							m_insect_constraints[i][0], m_insect_constraints[i][1], m_insect_constraints[i][2],
							alfa, beta, gamma);
				if (alfa < -_tol || beta < -_tol || gamma < -_tol || alfa >= 1.0 + _tol || beta >= 1.0 + _tol || gamma >= 1.0 + _tol){
					return false;
				}
			}
		}
	}
		
	return true;

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
void Operator::DrawActivePool(Tri*Triangles,int lf)
{
	//std::cout << "\n I AM DRAWING in Operator::DrawActivePool()" << std::endl;

	std::stringstream fname;
	fname << "debug_out/active_pool.obj";
	std::fstream file(fname.str().c_str(), std::ios::out);
	double x1, y1, z1, x2, y2, z2, x3, y3, z3;

	
	for (int V = 0; V<num_active; V++){		
		RetrieveCoordinates(Triangles,lf, active_pool[V], x1, y1, z1, x2, y2, z2, x3, y3, z3);
		file << "v " << x1 << " " << y1 << " " << z1 << std::endl;
		file << "v " << x2 << " " << y2 << " " << z2 << std::endl;
		file << "v " << x3 << " " << y3 << " " << z3 << std::endl;
	}
	for (int V = 1; V <= 3 * num_active; V += 3){
		file << "f " << V << " " << V + 1 << " " << V + 2 << std::endl;
	}
}
void Operator::DrawnList(Sphere* Spheres, std::vector<int>overlap)
{
	for (int i = 0; i <int(overlap.size()); i++){
		double x(0), y(0), z(0), r_2(0);
		int tri = 0;
		Spheres->Draw_m_sphere(overlap[i]);
	}	
}
void Operator::DrawnList(Sphere* Spheres, int*ip)
{
	for (int i = 1; i <=ip[0]; i++){
		double x(0), y(0), z(0), r_2(0);
		int tri = 0;
		Spheres->getCoordinates(ip[i], x, y, z, r_2, tri);
	
	}
}
void Operator::Draw_m_neighbor_tri(Tri*Triangles){
	std::ofstream file("debug_out/m_neighbor_tri.obj", std::ios::out);
	for (int i = 0; i < int(Triangles->coord.size()); i++){
		file <<"v "<< Triangles->coord[i][0] << " " <<Triangles->coord[i][1] << " "<< Triangles->coord[i][2] << std::endl;
	}

	for (int i = 1; i <= m_neighbor_tri[0]; i++){
		int t = m_neighbor_tri[i];
		file << "f " << Triangles->ids[t][0] + 1 << "	" << Triangles->ids[t][1] + 1 << "	" << Triangles->ids[t][2] + 1 << std::endl;
	}

	file.close();
}
void Operator::DrawTriangle(Tri*Triangles, int id){

	std::stringstream fname;
	fname << "debug_out/Tri_" << id << ".obj";
	std::fstream file(fname.str().c_str(), std::ios::out);

	file << "v " << Triangles->coord[Triangles->ids[id][0]][0] << " " << Triangles->coord[Triangles->ids[id][0]][1] << " " << Triangles->coord[Triangles->ids[id][0]][2] << std::endl;
	file << "v " << Triangles->coord[Triangles->ids[id][1]][0] << " " << Triangles->coord[Triangles->ids[id][1]][1] << " " << Triangles->coord[Triangles->ids[id][1]][2] << std::endl;
	file << "v " << Triangles->coord[Triangles->ids[id][2]][0] << " " << Triangles->coord[Triangles->ids[id][2]][1] << " " << Triangles->coord[Triangles->ids[id][2]][2] << std::endl;
	file << "f 1 2 3" << std::endl;
	file.close();
}
void Operator::DrawSeedList(Sphere *Spheres){
	for (int i = 1; i <= m_seeds_constraints[0]; i++){
		Spheres->Draw_m_seed(m_seeds_constraints[i]);
	}
}
