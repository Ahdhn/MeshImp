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

#include "KdTree.h"

#include <iostream>
#include <algorithm>
#include "Common.h"
KdTree::KdTree()
{
	
	
}

KdTree::~KdTree()
{

}

void KdTree::BuildTree(int numPoints, double**Verts, int dimension)
{
	
	DIM = dimension;
	myPoints = new kd_node_t[numPoints];
	for (int i = 0; i < numPoints; i++){
		myPoints[i].myID = i;
		myPoints[i].x = new double[DIM];
		for (int j = 0; j < DIM; j++){
			myPoints[i].x[j] = Verts[i][j];
		}
	}
	tmp = new double[DIM];
	root = Construct(myPoints, numPoints, 0);

}

void KdTree::BuildTree(int numPoints, vert*Verts, int dimension)
{
	//Build the tree, copy data and update the root 
	/*if (dimension < 0 || numPoints <= 0){
		std::cerr << "Error in KdTree::BuildTree:: Invalid parameters" << std::endl;
		exit(1);
	}*/

	DIM = dimension;
	myPoints = new kd_node_t[numPoints];
	for (int i = 0; i < numPoints; i++){
		myPoints[i].myID = i;
		myPoints[i].x = new double[DIM];
		for (int j = 0; j < DIM; j++){
			myPoints[i].x[j] = Verts[i].x[j]; 
		}
	}
	tmp = new double[DIM];
	root = Construct(myPoints, numPoints, 0);

}
struct kd_node_t* KdTree::Construct(struct kd_node_t *t, int len, int i)
{
	//Constrct the tree and update the root
	
	struct kd_node_t *n = NULL;
	if (n = FindMedian(t, t + len, i)){ 
		i = (i + 1) % DIM;
		n->left = Construct(t, n - t, i);
		n->right = Construct(n + 1, t + len - (n + 1), i);
	}
	return n;
	
}

int KdTree::FindNearest(double*point)
{
	testNode.x = point;
	
	struct kd_node_t *found = 0;
	
	double best_dist;
	
	Nearest(root, &testNode, 0, &found, &best_dist); 

	return found->myID;

}

void KdTree::rangeQuery(double point_xx, double point_yy, double point_zz, double r_2, int*inside, int&numInside){
	
	double*point = new double[3];
	point[0] = point_xx;
	point[1] = point_yy;
	point[2] = point_zz;
	
	rangeQuery(point, r_2, inside, numInside);
	
	delete[] point;
}
void KdTree::rangeQuery(double*point, double r_2, int*inside, int&numInside){
	testNode.x = point;
	
	Range(root, &testNode, 0, r_2, inside, numInside);
}

struct kd_node_t* KdTree::FindMedian(struct kd_node_t*start, struct kd_node_t*end, int id)
{
	if (end <= start)return NULL;
	if (end == start + 1){ return start; }

	struct  kd_node_t*p, *store, *md(start + (end - start) / 2);
	double pivot;

	
	while (1){
		pivot = md->x[id];
		Swap(md, end - 1);

		

		for (store = p = start; p < end; p++){
			if (p->x[id] < pivot){
				if (p != store){
					Swap(p, store);
				}
				store++;				
			}			
		}

		Swap(store, end - 1);
				
		if (store->x[id] == md->x[id]){ return md; }

		if (store > md)end = store;
		else start = store;
	}
	
}

void KdTree::Swap(struct kd_node_t*x, struct kd_node_t*y)
{
	for (int j = 0; j < DIM; j++){
		std::swap(x->x[j], y->x[j]);
	}
	std::swap(x->myID, y->myID);
}

void KdTree::Nearest(struct kd_node_t *root, struct kd_node_t *nd, int i, struct kd_node_t**best, double *best_dist)
{	
	if (!root){ return; }
	double d = Dist(root, nd);
	double dx = root->x[i] - nd->x[i];
	double dx2 = dx*dx;

	if (!*best || d < *best_dist){
		*best_dist = d;
		*best = root;
	}

	if (++i >= DIM){ i = 0; }

	Nearest(dx > 0 ? root->left : root->right, nd, i, best, best_dist);
	if (dx2 >= *best_dist)return;
	Nearest(dx > 0 ? root->right : root->left, nd, i, best, best_dist);
	
}

double KdTree::Dist(struct kd_node_t *a, struct kd_node_t *b)
{
	double t, d = 0;
	int myDim = DIM;
	while (myDim--) {
		t = a->x[myDim] - b->x[myDim];
		d += t * t;
	}
	return d;
}

void KdTree::Range(struct kd_node_t *root, struct kd_node_t *iPoint, int i, double r_2, int*inside, int&numInside){

	if (!root){ return; }
	
	double dx = root->x[i] - iPoint->x[i];
	double dx2 = dx*dx;
		
	if (r_2 > dx2){
		double d = Dist(root, iPoint);
		if (r_2 > d){
			inside[numInside] = root->myID;
			numInside++;
		}	
	}
	
	if (++i >= DIM){ i = 0; }

	Range(dx > 0 ? root->left : root->right, iPoint, i, r_2, inside, numInside);
	if (r_2 <= dx2){ return; }
	Range(dx > 0 ? root->right : root->left, iPoint, i, r_2, inside, numInside);
	/*if (r_2 > dx2){		
		double d = Dist(root, iPoint);
		if (d < r_2){
			inside[numInside] = root->myID;
			numInside++;
		}
		Range(root->left, iPoint, i, r_2, inside, numInside);
		Range(root->right, iPoint, i, r_2, inside, numInside);
	}
	else {
		//only one side 
		Range(dx > 0 ? root->left : root->right, iPoint, i, r_2, inside, numInside);
	}*/

}