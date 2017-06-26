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


#include "Constraints.h"
#include "../util/ThreeDRotation.h"
#include "../util/Common.h"

#include <stdlib.h>
#include <algorithm>

//#include "DrawSpheresDebug.h"//for debugging  

extern bool isObtuse;
extern bool isAcute;

Constraints::Constraints()
{
	numPlane = numExSphere = numInSphere = 0;
	numExSphere_size = numInSphere_size = numPlane_size = 999; 
	ExSphere = new sphere[1000];
	InSphere = new sphere[1000];
	Plane = new plane[1000];	
}
//TODO add maximality constraints as 3d points (on the original input surface) to be covered by the new vertex
//TODO add Hausdorff dist constriants 


Constraints::~Constraints(){}

void Constraints::ErrWarnMessage(size_t lineNum, std::string message, size_t mess_id)
{
	
	//mess_id =0 for error (exit)
	//otherwise, it is a warning (pause)

	if (mess_id == 0){
		fprintf(stderr, "\nError::line(%d)-->>%s", lineNum, message.c_str());
		system("pause");
	}
	else{
		fprintf(stderr, "\nWarning::line(%d)-->>%s", lineNum, message.c_str());
		system("pause");
	}
}

void Constraints::Reset(vert*Verts, int*nList)
{
	numPlane = numExSphere = numInSphere = 0;
	SetBoundingBox(Verts, nList);
}

//**** Construction
void Constraints::SetBoundingBox(vert*Verts, int*nList)
{
	//set the boundinng using nList
	//if you don't want the bounding box call this with NULL
	if (nList == NULL){
		myBox.xmin = myBox.ymin = myBox.zmin = DBL_MAX;
		myBox.xmax = myBox.ymax = myBox.zmax = DBL_MIN;
	}
	else{
		myBox.xmin = myBox.ymin = myBox.zmin = DBL_MAX;
		myBox.xmax = myBox.ymax = myBox.zmax = DBL_MIN;
		for (int i = 1; i <= nList[0]; i++){
			int ip = nList[i];
			myBox.xmin = std::min(myBox.xmin, Verts[ip].x[0]);
			myBox.ymin = std::min(myBox.ymin, Verts[ip].x[1]);
			myBox.zmin = std::min(myBox.zmin, Verts[ip].x[2]);

			myBox.xmax = std::max(myBox.xmax, Verts[ip].x[0]);
			myBox.ymax = std::max(myBox.ymax, Verts[ip].x[1]);
			myBox.zmax = std::max(myBox.zmax, Verts[ip].x[2]);
		}

		myBox.lx = myBox.xmax - myBox.xmin;
		myBox.ly = myBox.ymax - myBox.ymin;
		myBox.lz = myBox.zmax - myBox.zmin;
	}
}
bool Constraints::Delaunay(int*nList, int*skipList, vert*Verts, bool loadSkipList, bool isConnected)
{
	if (!isConnected){
		return DelaunayNotConnected(nList, skipList, Verts, loadSkipList);		
	}
	//TODO Delaunay_Extend where the list from which we calc ex sphere and in shperes are decoupled

	//nList is the sort list of neighbours around the a void 
	//1) find the exclusion spheres such that a new vertex should lie outside to maintain delaunayness of 
	//untouched/surounding triangles
	//2) find the inclusion spheres such that a new vertex should lie within to create delaunay-accurate set
	//of triangles 
	//skipList contians current vertices that should be negelected 

	//start by loading skipList to aux_list
	//for insert, we always wanna consider the vertices in skipList 
	if (loadSkipList){
		for (int i = 1; i <= skipList[0]; i++){
			aux_list[i + 3] = skipList[i];
		}
		aux_list[0] = skipList[0] + 3;
	}
	else{
		aux_list[0] = 3;
	}


	if (nList[0] == 3){
		//this is easy because this result into one InSphere (that connect these three vertices)
		//and three ExSpheres (no looping needed)

		int ap(nList[1]), n1(nList[2]), n2(nList[3]), ap2; 
		InSphere[numInSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
			                                           Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
			                                           Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
			                                           InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2]);
		numInSphere++;
		
		//have to expand skipList
		skipList[++skipList[0]] = ap;
		skipList[++skipList[0]] = n1;
		skipList[++skipList[0]] = n2;

		//first exsphere	
		ap2 = FindCommonElement_SkipList(Verts[ap].connect, Verts[n2].connect, skipList);
		ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
			                                           Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
			                                           Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
			                                           ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
		numExSphere++;

		//second exsphere
		ap2 = FindCommonElement_SkipList(Verts[ap].connect, Verts[n1].connect, skipList);
		ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
			                                           Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
			                                           Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
			                                           ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
		numExSphere++;

		//third exsphere
		ap2 = FindCommonElement_SkipList(Verts[n1].connect, Verts[n2].connect, skipList);
		ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
			                                           Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
			                                           Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
			                                           ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
		numExSphere++;
		skipList[0] -= 3;
		return  true;
	}


	

	int n2 = nList[nList[0]];
	for (int i = 1; i <= nList[0]; i++){
		int ap = nList[i];
		int n1 = (i == nList[0]) ? nList[1] : nList[i + 1];
		
		int ap2 = FindCommonElement_SkipList(Verts[ap].connect, Verts[n2].connect, skipList);
		//find the exclusion sphere (circumsphere of triangle n2-ap-ap2)
		if (ap2 >= 0){ 
			ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
								                           Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
								                           Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
								                           ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
			numExSphere++;
			if (numExSphere == numExSphere_size){ ExpandSpheres(numExSphere_size, numExSphere, ExSphere); }
		}
		else{
			ErrWarnMessage(__LINE__, "Constraints::Delaunay:: no valid apex was found", 0);
		}

		//find the inclusion sphere (circumsphere of triangle n2-ap-n1)
		InSphere[numInSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
			                                           Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
			                                           Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
													   InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2]);
		//this sphere could be too large or already contain other vertices
		//thus adding it in the containts set is useless
		if (InSphere[numInSphere].x[3] > 2){ 
			//because the domain is scaled inside the unit box 
			n2 = ap;
			continue; 
		}
		
		aux_list[1] = ap;
		aux_list[2] = n1;
		aux_list[3] = n2;

		if (!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], nList, aux_list, Verts) ||
			!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[ap].connect, aux_list, Verts) ||
			!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[n1].connect, aux_list, Verts) ||
			!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[n2].connect, aux_list, Verts)){
			//this circumsphere already contains other vertices and so it not a delaunay triangle 
			n2 = ap;
			
			continue;
		}

		numInSphere++;
		if (numInSphere == numInSphere_size){ ExpandSpheres(numInSphere_size, numInSphere, InSphere); }
		n2 = ap;
	}
	return true;	
}
bool Constraints::DelaunayNotConnected(int*nList, int*skipList, vert*Verts, bool loadSkipList)
{
	
	//nList is the sort list of neighbours around the a void but not a connected chain
	//i.e., last and first are not connected

	//1) find the exclusion spheres such that a new vertex should lie outside to maintain delaunayness of 
	//untouched/surounding triangles
	//2) find the inclusion spheres such that a new vertex should lie within to create delaunay-accurate set
	//of triangles 
	//skipList contians current vertices that should be negelected 

	//start by loading skipList to aux_list
	//for insert, we always wanna consider the vertices in skipList (then set loadSkipList to false)	
	if (loadSkipList){
		for (int i = 1; i <= skipList[0]; i++){
			aux_list[i + 3] = skipList[i];
		}
		aux_list[0] = skipList[0] + 3;
	}
	else{
		aux_list[0] = 3;
	}
	
	if (nList[0] == 3){
		//this is easy because this result into one InSphere (that connect these three vertices)
		//and two ExSpheres (no looping needed)

		int ap(nList[1]), n1(nList[2]), n2(nList[3]), ap2;
		InSphere[numInSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
			                                           Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
			                                           Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
			                                           InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2]);
		numInSphere++;

		
		//first exsphere
		ap2 = FindCommonElement_SkipList(Verts[ap].connect, Verts[n1].connect, skipList);
		ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
			                                           Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
			                                           Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
			                                           ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
		numExSphere++;

		//second exsphere
		ap2 = FindCommonElement_SkipList(Verts[n1].connect, Verts[n2].connect, skipList);
		ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
			                                           Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
			                                           Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
			                                           ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
		numExSphere++;		
		return true;
	}




	
	for (int i = 1; i < nList[0]; i++){
		int ap = nList[i];
		int n1 = nList[i + 1];

		int ap2 = FindCommonElement_SkipList(Verts[ap].connect, Verts[n1].connect, skipList);

		//find the exclusion sphere (circumsphere of triangle n2-ap-ap2)
		if (ap2 >= 0){
			ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
				                                           Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
				                                           Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
				                                           ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
			numExSphere++;
			if (numExSphere == numExSphere_size){ ExpandSpheres(numExSphere_size, numExSphere, ExSphere); }
		}
		else{
			return false;
			ErrWarnMessage(__LINE__, "Constraints::DelaunayNotConnected:: no valid apex was found", 0);
		}

		if (nList[0] - i >=2){
			int n2 = nList[i + 2];//two hubs away 

			//find the inclusion sphere (circumsphere of triangle n2-ap-n1)
			InSphere[numInSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
				                                           Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
				                                           Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
				                                           InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2]);
			//this sphere could be too large or already contain other vertices
			//thus adding it in the containts set is useless
			if (InSphere[numInSphere].x[3] > 2){
				//because the domain is scaled inside the unit box 				
				continue;
			}

			aux_list[1] = ap;
			aux_list[2] = n1;
			aux_list[3] = n2;

			if (!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], nList, aux_list, Verts) ||
				!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[ap].connect, aux_list, Verts) ||
				!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[n1].connect, aux_list, Verts) ||
				!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[n2].connect, aux_list, Verts)){
				//this circumsphere already contains other vertices and so it not a delaunay triangle 
				continue;
			}
			numInSphere++;
			if (numInSphere == numInSphere_size){ ExpandSpheres(numInSphere_size, numInSphere, InSphere); }
		}		
	}
	return true;
}
bool Constraints::DelaunayForcedNotConnected(int*nList, int removed, int*skipList, vert*Verts)
{
	//TODO add finding the circum spheres also (not well defined) 
	//force nList to be not connected by removing one vertex from it 
	//nList is the sort list of neighbours around the a void
	
	//1) find the exclusion spheres such that a new vertex should lie outside to maintain delaunayness of 
	//untouched/surounding triangles
	

	//start by loading skipList to aux_list
	//for insert, we always wanna consider the vertices in skipList 


	if (nList[0] == 3){
		return false;		
		ErrWarnMessage(__LINE__, "Constraints::DelaunayForcedNotConnected:: have not considered this case yet", 0);
		//this is easy because this result into one InSphere (that connect these three vertices)
		//and two ExSpheres (no looping needed)
		/*int ap(nList[1]), n1(nList[2]), n2(nList[3]), ap2;
		InSphere[numInSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
			Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
			Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
			InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2]);
		numInSphere++;


		//first exsphere
		ap2 = FindCommonElement_SkipList(Verts[ap].connect, Verts[n1].connect, skipList);
		ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
			Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
			Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
			ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
		numExSphere++;

		//second exsphere
		ap2 = FindCommonElement_SkipList(Verts[n1].connect, Verts[n2].connect, skipList);
		ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
			Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
			Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
			ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
		numExSphere++;
		skipList[0] -= 3;
		return;*/
	}


	int n1 = nList[nList[0]];
	for (int i = 1; i <= nList[0]; i++){ 
		int ap = nList[i];
		if (ap != removed && n1 != removed){
			int ap2 = FindCommonElement_SkipList(Verts[ap].connect, Verts[n1].connect, skipList);

			//find the exclusion sphere (circumsphere of triangle n2-ap-ap2)
			if (ap2 >= 0){
				ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
					                                           Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
					                                           Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
					                                           ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
				numExSphere++;
				if (numExSphere == numExSphere_size){ ExpandSpheres(numExSphere_size, numExSphere, ExSphere); }

			}
			else{
				return false;
				ErrWarnMessage(__LINE__, "Constraints::Delaunay:: no valid apex was found", 0);
			}
		}
		/*if (nList[0] - i > 2){
			int n2 = nList[i + 2];//two hubs away 

			//find the inclusion sphere (circumsphere of triangle n2-ap-n1)
			InSphere[numInSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
				Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
				Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
				InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2]);
			//this sphere could be too large or already contain other vertices
			//thus adding it in the containts set is useless
			if (InSphere[numInSphere].x[3] > 2){
				//because the domain is scaled inside the unit box 				
				continue;
			}

			aux_list[1] = ap;
			aux_list[2] = n1;
			aux_list[3] = n2;

			if (!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], nList, aux_list, Verts) ||
				!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[ap].connect, aux_list, Verts) ||
				!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[n1].connect, aux_list, Verts) ||
				!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[n2].connect, aux_list, Verts)){
				//this circumsphere already contains other vertices and so it not a delaunay triangle 
				continue;
			}
			numInSphere++;
			if (numInSphere == numInSphere_size){ ExpandSpheres(numInSphere_size, numInSphere, InSphere); }
		}*/
		n1 = ap;
	}
	return true;
}
bool Constraints::IsEmptySphere(double xc, double yc, double zc, double rc_2,
	                            int*list,
	                            int*skip_list,
								vert*Verts)
{
	for (int V = 1; V <= list[0]; V++){
		if (GetIndex(list[V], skip_list) >= 0){
			continue; 
		}
		double dist = Dist(Verts[list[V]].x[0], Verts[list[V]].x[1], Verts[list[V]].x[2], xc, yc, zc); 
		if (dist<rc_2 - _tol*_tol){
			return false;
		}
	}
	return true;
}
void Constraints::MinMaxAngle(int*nList, vert*Verts, double min_ang, double max_ang, double*void_vertex, bool nonObtuse, bool isConnected)
{
	//Find the inclusion and excluion region for all base angles such that a new vertex placed 
	//in the void surounded by nList and connected to all vertices in nList will not violate
	//the max_ang and min_ang constraints (just the base angles)
	//in case of nonObtuse = true, we additionally add diameter spheres such that apex angle is less than 90.0
	
	int n2 = nList[nList[0]];

	for (int i = 1; i <= nList[0]; i++){
		int n1 = nList[i];

		if (i == 1 && !isConnected){
			n2 = n1;
			continue;
		}
		if (min_ang > _tol){
			//*********** 1) Get min angle planes 
			//******* a) base angle 
			double rot_angle(90.0 - min_ang);
			double x_n, y_n, z_n;

			//normal to n1,n2, mid plane
			Cross(Verts[n2].x[0] - Verts[n1].x[0], Verts[n2].x[1] - Verts[n1].x[1], Verts[n2].x[2] - Verts[n1].x[2],
				  void_vertex[0] - Verts[n1].x[0], void_vertex[1] - Verts[n1].x[1], void_vertex[2] - Verts[n1].x[2], x_n, y_n, z_n);
			NormalizeVector(x_n, y_n, z_n);

			//rotate n2 by rot_angle about axis normal to the plane containing n1,n2 and mid
			//the axis passes through n1
			double x1, y1, z1;
			PerformRotation(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], x_n, y_n, z_n,
				            Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2], rot_angle, x1, y1, z1);
			//To be improved/revised
			//a naive way to check on the rotation direction
			//by doing it again with negative sign 
			//and check distanced min_point 
			double x_sub, y_sub, z_sub;
			rot_angle *= -1.0;
			PerformRotation(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], x_n, y_n, z_n,
				            Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2], rot_angle, x_sub, y_sub, z_sub);

			double len1(Dist(x_sub, y_sub, z_sub, void_vertex[0], void_vertex[1], void_vertex[2])), len2(Dist(x1, y1, z1, void_vertex[0], void_vertex[1], void_vertex[2]));
			if (len1 > len2){
				x1 = x_sub;
				y1 = y_sub;
				z1 = z_sub;
			}
			else{
				rot_angle *= -1.0;
			}
			

			SinglePlane(x1 - Verts[n1].x[0],
				        y1 - Verts[n1].x[1],
						z1 - Verts[n1].x[2],
						Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2]);


			//same thing is done on the other end of the edge
			Cross(Verts[n1].x[0] - Verts[n2].x[0], Verts[n1].x[1] - Verts[n2].x[1], Verts[n1].x[2] - Verts[n2].x[2],
				  void_vertex[0] - Verts[n2].x[0], void_vertex[1] - Verts[n2].x[1], void_vertex[2] - Verts[n2].x[2], x_n, y_n, z_n);
			NormalizeVector(x_n, y_n, z_n);
			double x2, y2, z2;
			PerformRotation(Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2], x_n, y_n, z_n,
				            Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], rot_angle, x2, y2, z2);
			

			SinglePlane(x2 - Verts[n2].x[0],
				        y2 - Verts[n2].x[1],
				        z2 - Verts[n2].x[2],
				        Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2]);
		}


		if (max_ang > _tol){
			//*********** 2) Get max angle planes 
			if (nonObtuse){
				//******* a) base angle 
				if (!isObtuse){
					SinglePlane(Verts[n1].x[0] - Verts[n2].x[0],
						Verts[n1].x[1] - Verts[n2].x[1],
						Verts[n1].x[2] - Verts[n2].x[2],
						Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2]);

					SinglePlane(Verts[n2].x[0] - Verts[n1].x[0],
						Verts[n2].x[1] - Verts[n1].x[1],
						Verts[n2].x[2] - Verts[n1].x[2],
						Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2]);
				}


				//******* b) apex angle 
				ExSphere[numExSphere].x[0] = (Verts[n1].x[0] + Verts[n2].x[0]) / 2.0;
				ExSphere[numExSphere].x[1] = (Verts[n1].x[1] + Verts[n2].x[1]) / 2.0;
				ExSphere[numExSphere].x[2] = (Verts[n1].x[2] + Verts[n2].x[2]) / 2.0;
				ExSphere[numExSphere].x[3] = Dist(ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2], Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2]);
								
				numExSphere++;
				if (numExSphere == numExSphere_size){ ExpandSpheres(numExSphere_size, numExSphere, ExSphere); }
			}
			else {
				if (true || !isAcute){
					//******* a) base angle 
					double rot_angle = abs(max_ang - 90.0);
					double x_n, y_n, z_n;
					//normal to n1,n2, mid plane	
					Cross(Verts[n2].x[0] - Verts[n1].x[0], Verts[n2].x[1] - Verts[n1].x[1], Verts[n2].x[2] - Verts[n1].x[2],
						  void_vertex[0] - Verts[n1].x[0], void_vertex[1] - Verts[n1].x[1], void_vertex[2] - Verts[n1].x[2], x_n, y_n, z_n);
					NormalizeVector(x_n, y_n, z_n);

					//rotate n2 by rot_angle about axis normal to the plane containing n1,n2 and mid
					//the axis passes through n1
					double x1, y1, z1;
					PerformRotation(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], x_n, y_n, z_n,
						            Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2], rot_angle, x1, y1, z1);
					//To be improved/revised
					//a naive way to check on the rotation direction
					//by doing it again with negative sign 
					//and check distanced min_point 
					double x_sub, y_sub, z_sub;
					rot_angle *= -1.0;
					PerformRotation(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], x_n, y_n, z_n,
						            Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2], rot_angle, x_sub, y_sub, z_sub);

					double len1 = Dist(x_sub, y_sub, z_sub, void_vertex[0], void_vertex[1], void_vertex[2]);
					double len2 = Dist(x1, y1, z1, void_vertex[0], void_vertex[1], void_vertex[2]);
					if (max_ang < 90.0 + _tol){ std::swap(len1, len2); }

					if (len1 < len2){
						x1 = x_sub;
						y1 = y_sub;
						z1 = z_sub;
					}
					else {
						rot_angle *= -1.0;
					}
					


					SinglePlane(Verts[n1].x[0] - x1,
								Verts[n1].x[1] - y1,
								Verts[n1].x[2] - z1,
								Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2]);


					//same thing is done on the other end of the edge
					Cross(Verts[n1].x[0] - Verts[n2].x[0], Verts[n1].x[1] - Verts[n2].x[1], Verts[n1].x[2] - Verts[n2].x[2],
						  void_vertex[0] - Verts[n2].x[0], void_vertex[1] - Verts[n2].x[1], void_vertex[2] - Verts[n2].x[2], x_n, y_n, z_n);
					NormalizeVector(x_n, y_n, z_n);
					double x2, y2, z2;
					PerformRotation(Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2], x_n, y_n, z_n,
									Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], rot_angle, x2, y2, z2);
					SinglePlane(Verts[n2].x[0] - x2,
								Verts[n2].x[1] - y2,
								Verts[n2].x[2] - z2,
								Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2]);
				}
			}
		}

		n2 = n1;
	}

}
void Constraints::MinEdgeLength(int*nList, vert*Verts, double(*sizingfunc)(double xx, double yy, double zz))
{
	//min edge lenght constraints (to respect the sizing function) 
	//it is based on the fact that new vertex in the void will be conntected to all 
	//vertices in nList 
	//sizingfunc should return the sizing function which is a protecteing sphere (ExSphere)
	//around each vertex in nList
	for (int i = 1; i <= nList[0]; i++){
		int n = nList[i];
		ExSphere[numExSphere].x[0] = Verts[n].x[0];
		ExSphere[numExSphere].x[1] = Verts[n].x[1];
		ExSphere[numExSphere].x[2] = Verts[n].x[2];
		ExSphere[numExSphere].x[3] = sizingfunc(Verts[n].x[0], Verts[n].x[1], Verts[n].x[2]);
		numExSphere++;
		if (numExSphere == numExSphere_size){ ExpandSpheres(numExSphere_size, numExSphere, ExSphere); }
	}
}
void Constraints::MinEdgeLength(int*nList, vert*Verts, double r_min_2)
{
	//min edge lenght constraints (to respect the sizing function) 
	//it is based on the fact that new vertex in the void will be conntected to all 
	//vertices in nList 
	//for unifrom sizing function, r_min_2 is square the sizing function 
	for (int i = 1; i <= nList[0]; i++){
		int n = nList[i];
		ExSphere[numExSphere].x[0] = Verts[n].x[0];
		ExSphere[numExSphere].x[1] = Verts[n].x[1];
		ExSphere[numExSphere].x[2] = Verts[n].x[2];
		ExSphere[numExSphere].x[3] = r_min_2;
				
		numExSphere++;
		if (numExSphere == numExSphere_size){ ExpandSpheres(numExSphere_size, numExSphere, ExSphere); }
	}
}
void Constraints::SinglePlane(double xn, double yn, double zn, double px, double py, double pz)
{
	//construct a single plance by specifiying the normal to the plane (xn, yn, zn)
	// and a point on it 
	NormalizeVector(xn, yn, zn);

	Plane[numPlane].x[0] = xn;
	Plane[numPlane].x[1] = yn;
	Plane[numPlane].x[2] = zn;
	Plane[numPlane].x[3] = -1.0*(px * Plane[numPlane].x[0] + py * Plane[numPlane].x[1] + pz * Plane[numPlane].x[2]); 
	
	numPlane++;
	if (numPlane == numPlane_size){ ExpandPlanes(numPlane_size, numPlane, Plane); }


}
void Constraints::Smoothness(int*nList, int*skipList, vert*Verts, double*void_vertex, double dev)
{	
	//we are looking for two planes
	//one of them is define by the mid point between n1-n2
	//and has a normal of plane containing n1,n2 and ap but tilted an angle=
	//the other has the same point but the normal tilted with -ve of previous angle
	if (nList[0] == 3){
		//expand the skipList to include nList to get the correct apex 
		skipList[++skipList[0]] = nList[1];
		skipList[++skipList[0]] = nList[2];
		skipList[++skipList[0]] = nList[3];
	}
	int n2 = nList[nList[0]];
	for (int i = 1; i <= nList[0]; i++){		
		int n1 = nList[i];
		int ap = FindCommonElement_SkipList(Verts[n1].connect, Verts[n2].connect, skipList);
		double tri_tri_angle = TriTriNormalAngle(Verts[n1].x, Verts[n2].x, Verts[ap].x, void_vertex);
		double rot_ang = dev;

		//the plane normal vector 
		double x_n, y_n, z_n;
		Cross(Verts[n1].x[0] - Verts[ap].x[0], Verts[n1].x[1] - Verts[ap].x[1], Verts[n1].x[2] - Verts[ap].x[2],
			  Verts[n2].x[0] - Verts[ap].x[0], Verts[n2].x[1] - Verts[ap].x[1], Verts[n2].x[2] - Verts[ap].x[2],
			  x_n, y_n, z_n);
		NormalizeVector(x_n, y_n, z_n);

		//point around which we will rotate the norm vector
		double x_mid = (Verts[n1].x[0] + Verts[n2].x[0]) / 2.0;
		double y_mid = (Verts[n1].x[1] + Verts[n2].x[1]) / 2.0;
		double z_mid = (Verts[n1].x[2] + Verts[n2].x[2]) / 2.0;

		//point on the vector 
		double pv_x = x_n + x_mid;
		double pv_y = y_n + y_mid;
		double pv_z = z_n + z_mid;
		double pv_x_rot, pv_y_rot, pv_z_rot;

		//first plane 
		PerformRotation(pv_x, pv_y, pv_z,
			            Verts[n1].x[0] - Verts[n2].x[0], Verts[n1].x[1] - Verts[n2].x[1], Verts[n1].x[2] - Verts[n2].x[2],
			            x_mid, y_mid, z_mid,
			            rot_ang,
			            pv_x_rot, pv_y_rot, pv_z_rot);

		Plane[numPlane].x[0] = pv_x_rot - x_mid;
		Plane[numPlane].x[1] = pv_y_rot - y_mid;
		Plane[numPlane].x[2] = pv_z_rot - z_mid;
		Plane[numPlane].x[3] = -1.0*(x_mid*Plane[numPlane].x[0] + y_mid*Plane[numPlane].x[1] + z_mid*Plane[numPlane].x[2]);
		
		if (Plane[numPlane].x[0] * void_vertex[0] + Plane[numPlane].x[1] * void_vertex[1] + Plane[numPlane].x[2] * void_vertex[2] + Plane[numPlane].x[3] > _tol){
			Plane[numPlane].x[0] *= -1.0;
			Plane[numPlane].x[1] *= -1.0;
			Plane[numPlane].x[2] *= -1.0;
			Plane[numPlane].x[3] = -1.0*(x_mid*Plane[numPlane].x[0] + y_mid*Plane[numPlane].x[1] + z_mid*Plane[numPlane].x[2]);

			if (Plane[numPlane].x[0] * void_vertex[0] + Plane[numPlane].x[1] * void_vertex[1] + Plane[numPlane].x[2] * void_vertex[2] + Plane[numPlane].x[3] > _tol){
				ErrWarnMessage(__LINE__, "Constraints::Smoothness:: error(1)", 1);				
			}
		}

		numPlane++;
		if (numPlane == numPlane_size){ ExpandPlanes(numPlane_size, numPlane, Plane); }

		//second plane 
		double rot_ang2 = -1.0*rot_ang;
		PerformRotation(pv_x, pv_y, pv_z,
			            Verts[n1].x[0] - Verts[n2].x[0],Verts[n1].x[1] -Verts[n2].x[1], Verts[n1].x[2] - Verts[n2].x[2],
			            x_mid, y_mid, z_mid,
			            rot_ang2,
			            pv_x_rot, pv_y_rot, pv_z_rot);

		Plane[numPlane].x[0] = pv_x_rot - x_mid;
		Plane[numPlane].x[1] = pv_y_rot - y_mid;
		Plane[numPlane].x[2] = pv_z_rot - z_mid;
		Plane[numPlane].x[3] = -1.0*(x_mid*Plane[numPlane].x[0] + y_mid*Plane[numPlane].x[1] + z_mid*Plane[numPlane].x[2]);

		
		if (Plane[numPlane].x[0] * void_vertex[0] + Plane[numPlane].x[1] * void_vertex[1] + Plane[numPlane].x[2] * void_vertex[2] + Plane[numPlane].x[3] > _tol){
			Plane[numPlane].x[0] *= -1.0;
			Plane[numPlane].x[1] *= -1.0;
			Plane[numPlane].x[2] *= -1.0;
			Plane[numPlane].x[3] = -1.0*(x_mid*Plane[numPlane].x[0] + y_mid*Plane[numPlane].x[1] + z_mid*Plane[numPlane].x[2]);

			if (Plane[numPlane].x[0] * void_vertex[0] + Plane[numPlane].x[1] * void_vertex[1] + Plane[numPlane].x[2] * void_vertex[2] + Plane[numPlane].x[3] > _tol){
				ErrWarnMessage(__LINE__, "Constraints::Smoothness:: error(2)", 1);
			}
		}		

		numPlane++;
		if (numPlane == numPlane_size){ ExpandPlanes(numPlane_size, numPlane, Plane); }

		n2 = n1;
	}
	if (nList[0] == 3){
		skipList[0] -= 3;
	}
}
void Constraints::Maximality()
{
	//TODO
}
void Constraints::AttractorInSphere(vert*Verts, int ip, double*void_vertex)
{
	//construct an inclusion sphere around void_vertex such that wehn ip is relocated
	//it is guarantteed to be closer to the void_vertex
	InSphere[numInSphere].x[0] = void_vertex[0];
	InSphere[numInSphere].x[1] = void_vertex[1];
	InSphere[numInSphere].x[2] = void_vertex[2];
	InSphere[numInSphere].x[3] = Dist(Verts[ip].x[0], Verts[ip].x[1], Verts[ip].x[2], void_vertex[0], void_vertex[1], void_vertex[2]);
	numInSphere++;
	if (numInSphere == numInSphere_size){ ExpandSpheres(numInSphere_size, numInSphere, InSphere); }
}
void Constraints::RepellerExSphere(vert*Verts, int ip, double*void_vertex)
{
	//construct an exclusion sphere around void_vertex such that wehn ip is relocated
	//it is guarantteed to be further away from the void_vertex
	ExSphere[numExSphere].x[0] = void_vertex[0];
	ExSphere[numExSphere].x[1] = void_vertex[1];
	ExSphere[numExSphere].x[2] = void_vertex[2];
	ExSphere[numExSphere].x[3] = Dist(Verts[ip].x[0], Verts[ip].x[1], Verts[ip].x[2], void_vertex[0], void_vertex[1], void_vertex[2]);
	numExSphere++;
	if (numExSphere == numExSphere_size){ ExpandSpheres(numExSphere_size, numExSphere, ExSphere); }
}
void Constraints::Direction(int*ip, int*nList, vert*Verts)
{
	//this builds a bounding polygon to prevent a newly relocated or created 
	//vertex from being created in a region that might create a tangled mesh
	//while preserving all other quality matric

	//loop around nList
	//for each edge in nList starting by n1-n2
	//find the common vertex of n1 and n2 in ip --> shrd
	//rotate the line between mid and n1 by 90 in direction closer to shrd 
	//construct the plane with normal as the rotate line 
	//and the point mid(n1-n2) as a point in it 

	//nList could be connected or not connected (does not make difference)
	
	int n2 = nList[nList[0]];
	for (int i = 1; i <= nList[0]; i++){
		int n1 = nList[i];

		//find the share vertex 
		if (FindCommonElements(Verts[n1].connect, Verts[n2].connect, aux_list)){
			int id = -1;
			for (int j = 1; j <= aux_list[0]; j++){						
				int sh = aux_list[j];
				id = (ip != NULL) ? GetIndex(aux_list[j], ip) : GetIndex(aux_list[j], nList); //for injection, the shared vertex is in nList

				/*if (j == aux_list[0] && id < 0 && nList[0] == 3 && ip[0] == 1){ 
					//special case of nList not connected with three vertices only
					id = 0;
					sh = ip[1];
				}*/
				if (id >= 0){
					
					double xmid(0.5*(Verts[n1].x[0] + Verts[n2].x[0])),
						   ymid(0.5*(Verts[n1].x[1] + Verts[n2].x[1])),
						   zmid(0.5*(Verts[n1].x[2] + Verts[n2].x[2]));

					//normal to n1,n2, sh
					double x_n12, y_n12, z_n12;
					Cross(Verts[n2].x[0] - Verts[n1].x[0], Verts[n2].x[1] - Verts[n1].x[1], Verts[n2].x[2] - Verts[n1].x[2],
						  Verts[sh].x[0] - Verts[n1].x[0], Verts[sh].x[1] - Verts[n1].x[1], Verts[sh].x[2] - Verts[n1].x[2], x_n12, y_n12, z_n12);
					NormalizeVector(x_n12, y_n12, z_n12);

					double x1, y1, z1;
					//PerformRotation(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], x_n12, y_n12, z_n12,
					//	            Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2], 90.0, x1, y1, z1);
					PerformRotation(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], x_n12, y_n12, z_n12,
						            xmid, ymid, zmid, 90.0, x1, y1, z1);
					double x_sub, y_sub, z_sub;
					//PerformRotation(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], x_n12, y_n12, z_n12,
					//	            Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2], -90.0, x_sub, y_sub, z_sub);
					PerformRotation(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], x_n12, y_n12, z_n12,
						            xmid, ymid, zmid, -90.0, x_sub, y_sub, z_sub);

					double len1(Dist(x_sub, y_sub, z_sub, Verts[sh].x[0], Verts[sh].x[1], Verts[sh].x[2])), len2(Dist(x1, y1, z1, Verts[sh].x[0], Verts[sh].x[1], Verts[sh].x[2]));
					if (len1 > len2){
						x1 = x_sub;
						y1 = y_sub;
						z1 = z_sub;
					}
					
					//normal
					//double xn(x1 - Verts[n1].x[0]),
					//	     yn(y1 - Verts[n1].x[1]),
					//	     zn(z1 - Verts[n1].x[2]);

					double xn(x1 - xmid),
						   yn(y1 - ymid),
						   zn(z1 - zmid);
					SinglePlane(xn, yn, zn, xmid, ymid, zmid);
					//SinglePlane(xn, yn, zn, Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2]);
					break;					
				}
			}
			if (id < 0 && nList[0] != 3){
				ErrWarnMessage(__LINE__, "Constraints::Direction edge has incorrect shared vertex!!!", 0);
			}
		}
		else{
			ErrWarnMessage(__LINE__, "Constraints::Direction edge does not have shared vertex!!!", 0);
		}
		n2 = n1;
	}
}
//*** Checking 
bool Constraints::OverlappingInSpheres()
{
	//find if all InSpheres overlaps
	//if two don't overlap, return false
	for (int i = 0; i < numInSphere - 1; i++){
		for (int j = i + 1; j < numInSphere; j++){
			if (Dist(InSphere[i].x[0], InSphere[i].x[1], InSphere[i].x[2], InSphere[j].x[0], InSphere[j].x[1], InSphere[j].x[2]) *_tol_sq_circ >
				InSphere[i].x[3] + InSphere[j].x[3] + 2.0*sqrt(InSphere[i].x[3] * InSphere[j].x[3])){
				return false;
			}
		}
	}

	return true;
}
bool Constraints::InsideFeasibleRegion_Vertex(double xx, double yy, double zz)
{
	//check if (xx,yy,zz) is insdie all the feasible region 
	//i.e., inside all inclusion regions and outside all exclusion regions 
	if (!IsInsideBoundingBox(xx, yy, zz)){
		return false;
	}

	for (int i = 0; i < numExSphere; i++){
		double dist = Dist(xx, yy, zz, ExSphere[i].x[0], ExSphere[i].x[1], ExSphere[i].x[2]);
		if (dist < ExSphere[i].x[3] + _tol_sq){ return false; }
	}

	for (int i = 0; i < numInSphere; i++){
		double dist = Dist(xx, yy, zz, InSphere[i].x[0], InSphere[i].x[1], InSphere[i].x[2]);	
		if (dist > InSphere[i].x[3] - _tol_sq){ return false; }
	}

	for (int i = 0; i < numPlane; i++){
		if (xx*Plane[i].x[0] + yy*Plane[i].x[1] + zz*Plane[i].x[2] + Plane[i].x[3]>-_tol){			
			return false;
		}
	}

	return true;
}
bool Constraints::InsideFeasibleRegion_Triangle(double*tri)
{
	//tri is the coordinates of the triangle
	//0-2 is first vertex
	//3-5 is second vertex
	//6-8 is third vertex

	for (int i = 0; i < numExSphere; i++){
		double dist1 = Dist(tri[0], tri[1], tri[2], ExSphere[i].x[0], ExSphere[i].x[1], ExSphere[i].x[2]);
		double dist2 = Dist(tri[3], tri[4], tri[5], ExSphere[i].x[0], ExSphere[i].x[1], ExSphere[i].x[2]);
		double dist3 = Dist(tri[6], tri[7], tri[8], ExSphere[i].x[0], ExSphere[i].x[1], ExSphere[i].x[2]);

		if (dist1<ExSphere[i].x[3] * _tol_sq_circ && dist2<ExSphere[i].x[3] * _tol_sq_circ && dist3<ExSphere[i].x[3] * _tol_sq_circ){
			//the whole triangle is in one of the exclusion sphere
			return false;
		}
	}
	
	for (int i = 0; i < numInSphere; i++){
		double dist1 = Dist(tri[0], tri[1], tri[2], InSphere[i].x[0], InSphere[i].x[1], InSphere[i].x[2]);
		double dist2 = Dist(tri[3], tri[4], tri[5], InSphere[i].x[0], InSphere[i].x[1], InSphere[i].x[2]);
		double dist3 = Dist(tri[6], tri[7], tri[8], InSphere[i].x[0], InSphere[i].x[1], InSphere[i].x[2]);

		if (dist1 > InSphere[i].x[3] * _tol_sq_circ && dist2 > InSphere[i].x[3] * _tol_sq_circ && dist3 > InSphere[i].x[3] * _tol_sq_circ){
			//the whole triangle is in one of the exclusion sphere
			return false;
		}
	}

	
	for (int i = 0; i < numPlane; i++){
		if (tri[0] * Plane[i].x[0] + tri[1] * Plane[i].x[1] + tri[2] * Plane[i].x[2] + Plane[i].x[3] >-_tol &&
			tri[3] * Plane[i].x[0] + tri[4] * Plane[i].x[1] + tri[5] * Plane[i].x[2] + Plane[i].x[3] >-_tol &&
			tri[6] * Plane[i].x[0] + tri[7] * Plane[i].x[1] + tri[8] * Plane[i].x[2] + Plane[i].x[3] >-_tol ){

			return false;
		}
	}
	return true;

}
bool Constraints::IsInsideBoundingBox(double xx, double yy, double zz)
{
	//check if (xx,yy,zz) is inside the bounding box myBox
	if (abs(myBox.xmin - DBL_MAX) < _tol){
		//if the box did not set, we assume it is okay
		return true;
	}
	else{
		if (xx < myBox.xmin || xx > myBox.xmax ||
			yy < myBox.ymin || yy > myBox.ymax ||
			zz < myBox.zmin || xx > myBox.zmax ){
			return false;
		}
	}
	return true;
}


void Constraints::ExpandSpheres(int&currentSize,int currentNumSphere, sphere*&mySpheres)
{
	currentSize *= 2; 
	sphere*newSpheres = new sphere[currentSize];
	for (int i = 0; i < currentNumSphere; i++){
		for (int j = 0; j < 4; j++){
			newSpheres[i].x[j] = mySpheres[i].x[j];
		}
	}

	delete[]mySpheres;
	
	mySpheres = new sphere[currentSize];
	for (int i = 0; i < currentNumSphere; i++){
		for (int j = 0; j < 4; j++){
			newSpheres[i].x[j] = mySpheres[i].x[j];
		}
	}	
	currentSize--;
}
void Constraints::ExpandPlanes(int&currentSize, int currentNumPlane, plane*&myPlanes)
{
	currentSize *= 2;
	sphere*newPlanes = new sphere[currentSize];
	for (int i = 0; i < currentNumPlane; i++){
		for (int j = 0; j < 4; j++){
			newPlanes[i].x[j] = myPlanes[i].x[j];
		}
	}

	delete[]myPlanes;

	myPlanes = new plane[currentSize];
	for (int i = 0; i < currentNumPlane; i++){
		for (int j = 0; j < 4; j++){
			newPlanes[i].x[j] = myPlanes[i].x[j];
		}
	}
	currentSize--;
}

//*** Debug
void Constraints::SpheresDrawing(std::string filename, int num, sphere *Sphere)
{
	double **drawpsheres = new double*[num];
	for (int i = 0; i < num; i++){
		drawpsheres[i] = new double[4];
		drawpsheres[i][0] = Sphere[i].x[0];
		drawpsheres[i][1] = Sphere[i].x[1];
		drawpsheres[i][2] = Sphere[i].x[2];
		drawpsheres[i][3] = Sphere[i].x[3];
	}

	//DrawManySpheres(filename, num, drawpsheres, 1);
	for (int i = 0; i < num; i++){
		delete[] drawpsheres[i];
	}
	free(drawpsheres);
}
void Constraints::DrawInSpheres()
{
	///SpheresDrawing("debug_out/InSphere.obj",numInSphere, InSphere);
}
void Constraints::DrawExSpheres()
{
	//SpheresDrawing("debug_out/ExSphere.obj", numExSphere, ExSphere);
}