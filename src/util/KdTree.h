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

//Modified version of https://rosettacode.org/wiki/K-d_tree
#ifndef _KDTREE_
#define _KDTREE_

#define MAX_DIM 3
#include "Common.h"

struct kd_node_t{
	double*x;
	int myID;//the id as in the passed point cloud 
	struct kd_node_t *left, *right;
};

class KdTree
{
public:
	KdTree();

	//Call this one time to build everything 
	void BuildTree(int numPoints, vert*Verts, int dimension);
	void BuildTree(int numPoints, double**Verts, int dimension);
	
	//Call this to get the nearest node
	int FindNearest(double*point);

	//find all points inside a sphere or radius (square) or r_2 around iPoint
	void rangeQuery(double*point, double r_2, int*inside, int&numInside);
	void rangeQuery(double point_xx, double point_yy, double point_zz, double r_2, int*inside, int&numInside);

	~KdTree();
	
private:
	double Dist(struct kd_node_t*, struct kd_node_t*, int);
	void Swap(struct kd_node_t*, struct kd_node_t*);
	struct kd_node_t* FindMedian(struct kd_node_t*, struct kd_node_t*, int);	
	struct kd_node_t* Construct(struct kd_node_t *t, int len, int i);
	void Nearest(struct kd_node_t *, struct kd_node_t *, int, struct kd_node_t**, double *);
	double Dist(struct kd_node_t *, struct kd_node_t *);
	void Range(struct kd_node_t *root, struct kd_node_t *iPoint, int i, double r_2, int*inside, int&numInside);

	int DIM;
	
	struct kd_node_t testNode;

	double*tmp; 
	struct kd_node_t *myPoints;
	struct kd_node_t *root;
	
};



#endif