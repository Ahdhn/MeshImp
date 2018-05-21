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

//header file for constants, math and geometric functions

#ifndef _COMMON_
#define _COMMON_

#define MAX_CONNECT 20 // the max valence for a triangle 

#define PI 3.14159265358979323846264338327950288
#define TwoPI 6.2831853071795862

#define _tol 10E-8 //this tolerance is based on the unit box. may change for different configuration
#define _tol_sq _tol*_tol
#define _tol_sq_circ (_tol+1.0)*(_tol+1.0) //use this when compare distance^2 with radius^2
                                           //basically to check if the a point is inside a circle

#define RadToDeg 57.295779513078550 //180.0/PI (muliply this by the radian angle to convert to degree)

#define _ang_tol_2 5.0 //angle tolerance
#define _ang_tol 3.0 //angle tolerance
#define SIZE_T_MAX ((size_t) -1)

#include <math.h>
#include <assert.h>
#include <array>
#include <vector>
struct vert{
	double x[3];//x,y,z
	int connect[MAX_CONNECT];
	double dih_ang;
};

//Structs 
struct BoundBox{
	double xmin, xmax, ymin, ymax, zmin, zmax, lx, ly, lz;
};
struct sphere
{
	double x[4];//x,y,z and r
};
struct plane
{
	double x[4]; //a, b,c and d in the scalar equation of a plane

};
struct Tri {
	std::vector<std::array<int,3>> ids;//triangles three vertices id
	std::vector<std::array<int, 3>> neighbour; //num_tri X 3
	std::vector<std::array<double, 3>> coord;//vertices coordinates  X 3
};

struct stats
{
	double min_ang, max_ang, min_dih_ang, max_dih_ang, min_tri_quality, max_tri_quality;
};

template<typename T>
int GetIndex(T entry, T*list)
{
	//get the index of entry in array (list) of type T
	//if entry is not there, return -1
	//list store the total number of its entries in first element [0]
	for (int i = 1; i <= list[0]; i++){
		if (list[i] == entry){
			return i;
		}
	}
	return -1;
}

template <typename T>
bool FindCommonElements(T*arr1, T*arr2, T commonEle[])
{
	//store the common elements in arr1 and arr2 in commonEle
	//return true if it is more than one elements 

	commonEle[0] = 0;
	int num_arr1(arr1[0]), num_arr2(arr2[0]); 
	for (int i = 1; i <= num_arr1; i++){
		for (int j = 1; j <= num_arr2; j++){ 
			if (arr1[i] == arr2[j]){
				commonEle[++commonEle[0]] = arr1[i];				
			}
		}
	}

	if (commonEle[0] > 0){ return true; }
	else{ return false; }
}

template <typename T>
bool FindCommonElements_SkipList(T*arr1, T*arr2, T*skip_list, T commonEle[])
{
	//store the common elements in arr1 and arr2 in commonEle
	//return true if it is more than one elements 

	commonEle[0] = 0;
	int num_arr1(arr1[0]), num_arr2(arr2[0]), num_skip(skip_list[0]); 

	for (int i = 1; i <= num_arr1; i++){
		if (GetIndex(arr1[i], skip_list) >= 0){ continue; }
		for (int j = 1; j <= num_arr2; j++){
			if (arr1[i] == arr2[j]){
				commonEle[++commonEle[0]] = arr1[i];
			}
		}
	}

	if (commonEle[0] > 0){ return true; }
	else{ return false; }
}

template <typename T>
T FindCommonElement_SkipList(T*arr1, T*arr2, T*skip_list)
{
		
	//return the single common elements bteween arr1 and arr2 
	//make no check if there is more than one 
	int num_arr1(arr1[0]), num_arr2(arr2[0]), num_skip(skip_list[0]);

	for (int i = 1; i <= num_arr1; i++){
		if (GetIndex(arr1[i], skip_list) >= 0){ continue; }
		for (int j = 1; j <= num_arr2; j++){
			if (arr1[i] == arr2[j]){
				return arr1[i];
			}
		}
	}

	return -1;
}

template <typename T>
T FindCommonElement_SkipList(T*arr1, T*arr2, T skip_element)
{

	//return the single common elements bteween arr1 and arr2 
	//make no check if there is more than one 
	int num_arr1(arr1[0]), num_arr2(arr2[0]);

	for (int i = 1; i <= num_arr1; i++){
		if (arr1[i] == skip_element){ continue; }
		for (int j = 1; j <= num_arr2; j++){
			if (arr1[i] == arr2[j]){
				return arr1[i];
			}
		}
	}

	return -1;
}

template <typename T>
inline bool FindDuplication(T*arr1)
{
	//find if there is any duplication in arr1
	//arr1 lenght is stored in first element of it
	if (arr1[0] < 2){ return false; }
	for (int i = 1; i <= arr1[0] - 1; i++){
		for (int j = i + 1; j <= arr1[0]; j++){
			if (arr1[i] == arr1[j]){ return true; }
		}
	}
	return false;
}

template <typename T>
inline void RemoveNodeFromList(T*list, T entry)
{
	//remove entry from list while preseving its sorted format
	int tmp, tmp1, d;
	bool find(false);
	for (d = 1; d <= list[0]; d++){
		if (list[d] == entry){ find = true; break; }
	}
	//entry is not in the list
	if (!find){ return; }
	d = list[0];
	tmp = list[d];
	while (tmp != entry){
		tmp1 = list[d - 1];
		list[d - 1] = tmp;
		tmp = tmp1;
		d--;
	}
	list[0]--;
}

template <typename T>
void AddEntrySortedList(T entry, T* list, T prv, T nxt)
{
	//add entry to the sorted list 
	//such that the new entry fits between prv_entry and nxt_entry 
	//length of list is stored in list[0]

	if ((prv == list[list[0]] && nxt == list[1]) || (nxt == list[list[0]] && prv == list[1])){ 
		//it fits to be in the last or the start of list		
		list[++list[0]] = entry; 
		return;
	}

	if (GetIndex(prv, list) > GetIndex(nxt, list)){ 
		std::swap(prv, nxt); 
	}

	int i(++list[0]);
	while (list[i - 1] != prv){ 
		list[i] = list[i - 1];
		--i;
	}
	list[i] = entry;
}

inline void NormalizeVector(double&vector_x, double&vector_y, double&vector_z)
{
	double nn = sqrt(vector_x*vector_x + vector_y*vector_y + vector_z*vector_z);
	vector_x /= nn; vector_y /= nn; vector_z /= nn;
}

inline void Cross(double xv1, double yv1, double zv1, double xv2, double yv2, double zv2, double&xx, double&yy, double&zz)
{
	xx = yv1*zv2 - zv1*yv2;
	yy = zv1*xv2 - xv1*zv2;
	zz = xv1*yv2 - yv1*xv2;
}
inline double Dot(double xv1, double yv1, double zv1, double xv2, double yv2, double zv2)
{	
	return xv1*xv2 + yv1*yv2 + zv1*zv2;	
}
inline void PlaneNorm(double ux, double uy, double uz, double vx, double vy, double vz, double&x, double&y, double&z)
{
	//normal vector to plane defined by two vectors

	x = uy*vz - uz*vy;
	y = uz*vx - ux*vz;
	z = ux*vy - uy*vx;

	double len = sqrt(x*x + y*y + z*z);
	x /= len;
	y /= len;
	z /= len;
}

inline double Angle360Vectors(double dxba, double dyba, double dzba, double dxca, double dyca, double dzca, double norm_x, double norm_y, double norm_z)
{
	// return angle in degrees between b.a.c from (0,2PI)
	//(norm_x,norm_y,norm_z) are the perpendicular normal vector to the triangle (bac)

	double dot, pcross, theta;

	double i, j, k; // cross(ab,ac)

	dot = dxba*dxca + dyba*dyca + dzba*dzca; // dot(ab,ac)

	Cross(dxba, dyba, dzba, dxca, dyca, dzca, i, j, k);
	
	pcross = i*norm_x + j*norm_y + k*norm_z;

	theta = atan2(pcross, dot);

	if (theta < 0){ theta += TwoPI; }

	return theta *RadToDeg;//convert to deg


}

inline double AngleVectVect(double x1, double y1, double z1, double x2, double y2, double z2) 
{
	double dot = x1*x2 + y1*y2 + z1*z2;

	if (dot == 0.0) { return 1.5707963268; }//90 deg

	double angle = dot / sqrt((x1*x1 + y1*y1 + z1*z1) * (x2*x2 + y2*y2 + z2*z2));
	
	if (angle>1.0){ return 0.0; }
	if (angle<-1.0){ return PI; }

	return acos(angle);
}

inline double Dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double dx, dy, dz;
	dx = x1 - x2;
	dy = y1 - y2;
	dz = z1 - z2;
	dx *= dx;
	dy *= dy;
	dz *= dz;

	return dx + dy + dz;

}

inline double TriCircumcenter3d(double xa, double ya, double za,
                                double xb, double yb, double zb,
 						        double xc, double yc, double zc,
                                double&x_cir, double&y_cir, double&z_cir)
								/*(double*a,double*b,double*c,double*circumcenter,double*xi,double*eta)*/
{
	/*****************************************************************************/
	/*
	/*  tricircumcenter3d()   Find the circumcenter of a triangle in 3D.         */
	/*                                                                           */
	/*  The result is returned both in terms of xyz coordinates and xi-eta       */
	/*  coordinates, relative to the triangle's point `a' (that is, `a' is       */
	/*  the origin of both coordinate systems).  Hence, the xyz coordinates      */
	/*  returned are NOT absolute; one must add the coordinates of `a' to        */
	/*  find the absolute coordinates of the circumcircle.  However, this means                                              */
	/*  that the result is frequently more accurate than would be possible if    */
	/*  absolute coordinates were returned, due to limited floating-point        */
	/*  precision.  In general, the circumradius can be computed much more       */
	/*  accurately.                                                              */
	/*                                                                           */
	/*  The xi-eta coordinate system is defined in terms of the triangle.        */
	/*  Point `a' is the origin of the coordinate system.  The edge `ab' extends */
	/*  one unit along the xi axis.  The edge `ac' extends one unit along the    */
	/*  eta axis.  These coordinate values are useful for linear interpolation.  */
	/*                                                                           */
	/*  If `xi' is NULL on input, the xi-eta coordinates will not be computed.   */
	/*                                                                           */
	/*****************************************************************************/
	//http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html 
	//http://gamedev.stackexchange.com/questions/60630/how-do-i-find-the-circumcenter-of-a-triangle-in-3d
	double xba, yba, zba, xca, yca, zca;
	double balength, calength;
	double xcrossbc, ycrossbc, zcrossbc;
	double denominator;
	double xcirca, ycirca, zcirca;

	/* Use coordinates relative to point `a' of the triangle. */ 
	xba = xb - xa;
	yba = yb - ya;
	zba = zb - za;
	xca = xc - xa;
	yca = yc - ya;
	zca = zc - za;
	/* Squares of lengths of the edges incident to `a'. */
	balength = xba * xba + yba * yba + zba * zba;
	calength = xca * xca + yca * yca + zca * zca;

//	/* Cross product of these edges. */
//#ifdef EXACT
//	/* Use orient2d() from http://www.cs.cmu.edu/~quake/robust.html     */
//	/*   to ensure a correctly signed (and reasonably accurate) result, */
//	/*   avoiding any possibility of division by zero.                  */
//	xcrossbc = orient2d(b[1], b[2], c[1], c[2], a[1], a[2]);
//	ycrossbc = orient2d(b[2], b[0], c[2], c[0], a[2], a[0]);
//	zcrossbc = orient2d(b[0], b[1], c[0], c[1], a[0], a[1]);
//#else
//	/* Take your chances with floating-point roundoff. */
//	xcrossbc = yba * zca - yca * zba;
//	ycrossbc = zba * xca - zca * xba;
//	zcrossbc = xba * yca - xca * yba;
//#endif

	/* Take your chances with floating-point roundoff. */
	xcrossbc = yba * zca - yca * zba;
	ycrossbc = zba * xca - zca * xba;
	zcrossbc = xba * yca - xca * yba;

	/* Calculate the denominator of the formulae. */
	denominator = 0.5 / (xcrossbc * xcrossbc + ycrossbc * ycrossbc + zcrossbc * zcrossbc);

	/* Calculate offset (from `a') of circumcenter. */
	xcirca = ((balength * yca - calength * yba) * zcrossbc -
		(balength * zca - calength * zba) * ycrossbc) * denominator;
	ycirca = ((balength * zca - calength * zba) * xcrossbc -
		(balength * xca - calength * xba) * zcrossbc) * denominator;
	zcirca = ((balength * xca - calength * xba) * ycrossbc -
		(balength * yca - calength * yba) * xcrossbc) * denominator;
	x_cir = xcirca + xa;
	y_cir = ycirca + ya;
	z_cir = zcirca + za;

	//if (xi != (double *) NULL) {
	/* To interpolate a linear function at the circumcenter, define a     */
	/*   coordinate system with a xi-axis directed from `a' to `b' and    */
	/*   an eta-axis directed from `a' to `c'.  The values for xi and eta */
	/*   are computed by Cramer's Rule for solving systems of linear      */
	/*   equations.                                                       */

	/* There are three ways to do this calculation - using xcrossbc, */
	/*   ycrossbc, or zcrossbc.  Choose whichever has the largest    */
	/*   magnitude, to improve stability and avoid division by zero. */
	/*if (((xcrossbc >= ycrossbc) ^ (-xcrossbc > ycrossbc)) &&
	((xcrossbc >= zcrossbc) ^ (-xcrossbc > zcrossbc))) {
	*xi = (ycirca * zca - zcirca * yca) / xcrossbc;
	*eta = (zcirca * yba - ycirca * zba) / xcrossbc;
	} else if ((ycrossbc >= zcrossbc) ^ (-ycrossbc > zcrossbc)) {
	*xi = (zcirca * xca - xcirca * zca) / ycrossbc;
	*eta = (xcirca * zba - zcirca * xba) / ycrossbc;
	} else {
	*xi = (xcirca * yca - ycirca * xca) / zcrossbc;
	*eta = (ycirca * xba - xcirca * yba) / zcrossbc;
	}
	}*/

	double len1, dx, dy, dz;
	dx = xa - x_cir;
	dy = ya - y_cir;
	dz = za - z_cir;
	len1 = dx*dx + dy*dy + dz*dz;
	if (len1>2.0){ return len1; }

/*#ifdef  debug	

	dx = xb - x_cir;
	dy = yb - y_cir;
	dz = zb - z_cir;
	len2 = dx*dx + dy*dy + dz*dz;

	dx = xc - x_cir;
	dy = yc - y_cir;
	dz = zc - z_cir;
	len3 = dx*dx + dy*dy + dz*dz;

	if (abs(len1 - len2)>_tol || abs(len3 - len2)>_tol || abs(len1 - len3)>_tol){
		//return 1000;
		cout << "Error at TriCircumcenter3d()..!! " << endl;
		//system("pause");
	}
#endif*/

	return len1;


}

inline double TriTriNormalAngle3(double xip, double yip, double zip, double*ip1, double*ip2, double*ip3)
{
	//return the angles between the normal of the planes containing the two triangles sharing an edge ip-ip1
	//watch out for how you call it. the order of the points is curcial for return the correct angle

	double x_n1, y_n1, z_n1, x_n2, y_n2, z_n2, x_np, y_np, z_np, angle360;

	Cross(ip1[0] - xip, ip1[1] - yip, ip1[2] - zip,
		  ip2[0] - xip, ip2[1] - yip, ip2[2] - zip,
		  x_n1, y_n1, z_n1);
	Cross(ip3[0] - xip, ip3[1] - yip, ip3[2] - zip,
		  ip1[0] - xip, ip1[1] - yip, ip1[2] - zip,
		  x_n2, y_n2, z_n2);

	PlaneNorm(x_n1, y_n1, z_n1, x_n2, y_n2, z_n2, x_np, y_np, z_np);

	angle360 = Angle360Vectors(x_n1, y_n1, z_n1, x_n2, y_n2, z_n2, x_np, y_np, z_np);

	return angle360;
}

inline double TriTriNormalAngle(double*ip, double*ip1, double*ip2, double*ip3)
{
	return TriTriNormalAngle3(ip[0], ip[1], ip[2], ip1, ip2, ip3);
}

inline double PointTriangleDistance(double xp1, double yp1, double zp1,// triangle head1
	                                double xp2, double yp2, double zp2,// triangle head2 
	                                double xp3, double yp3, double zp3,// triangle head3
	                                double xp, double yp, double zp, // my point 
	                                double&x_new, double&y_new, double&z_new)// my projection 
{
	// http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
	// http://www.mathworks.com/matlabcentral/fileexchange/22857-distance-between-a-point-and-a-triangle-in-3d/content/pointTriangleDistance.m
	//        ^t
	//  \     |
	//   \reg2|
	//    \   |
	//     \  |
	//      \ |
	//       \|
	//        *P2
	//        |\
	//        | \
	//  reg3  |  \ reg1
	//        |   \
	//        |reg0\ 
	//        |     \ 
	//        |      \ P1
	// -------*-------*------->s
	//        |P0      \ 
	//  reg4  | reg5    \ reg6

	/*size_t p1(_surface_tri2[tri][0]),p2(_surface_tri2[tri][1]),p3(_surface_tri2[tri][2]);
	double xp1(_surface[p1][0]),yp1(_surface[p1][1]),zp1(_surface[p1][2]),
	xp2(_surface[p2][0]),yp2(_surface[p2][1]),zp2(_surface[p2][2]),
	xp3(_surface[p3][0]),yp3(_surface[p3][1]),zp3(_surface[p3][2]);*/

	double xE0(xp2 - xp1), yE0(yp2 - yp1), zE0(zp2 - zp1),
		xE1(xp3 - xp1), yE1(yp3 - yp1), zE1(zp3 - zp1),
		xD(xp1 - xp), yD(yp1 - yp), zD(zp1 - zp);

	double a = Dot(xE0, yE0, zE0, xE0, yE0, zE0);
	double b = Dot(xE0, yE0, zE0, xE1, yE1, zE1);
	double c = Dot(xE1, yE1, zE1, xE1, yE1, zE1);
	double d = Dot(xE0, yE0, zE0, xD, yD, zD);
	double e = Dot(xE1, yE1, zE1, xD, yD, zD);
	double f = Dot(xD, yD, zD, xD, yD, zD);

	double det = a*c - b*b;
	double s = b*e - c*d;
	double t = b*d - a*e;
	double invDet, dist, numer, denom, tmp0, tmp1;
	//Terible tree of conditionals to determine in which region of the diagram
	// shown above the projection of the point into the triangle-plane lies.

	if (s + t <= det){
		if (s<0){
			if (t<0){
				//region 4
				if (d < 0){
					t = 0;
					if (-1.0*d >= a){
						s = 1;
						dist = a + 2.0*d + f;
					}
					else{
						s = -1.0*d / a;
						dist = d*s + f;
					}
				}
				else{
					s = 0;
					if (e >= 0){
						t = 0;
						dist = f;
					}
					else{
						if (-e >= c){
							t = 1;
							dist = c + 2.0*e + f;
						}
						else{
							t = -1.0*e / c;
							dist = e*t + f;
						}

					}

				}

			}
			else
			{
				//region 3
				s = 0;
				if (e >= 0){
					t = 0;
					dist = f;
				}
				else{
					if (-1.0*e >= c){
						t = 1;
						dist = c + 2.0*e + f;
					}
					else{
						t = -1.0*e / c;
						dist = e*t + f;
					}
				}

			}
		}
		else if (t<0){
			//region 5
			t = 0;
			if (d >= 0){
				s = 0;
				dist = f;
			}
			else{
				if (-1.0*d >= a){
					s = 1;
					dist = a + 2.0*d + f;
				}
				else{
					s = -1.0*d / a;
					dist = d*s + f;
				}
			}

		}
		else{
			//region 0
			invDet = 1.0 / det;
			s = s*invDet;
			t = t*invDet;
			dist = s*(a*s + b*t + 2.0*d) +
				t*(b*s + c*t + 2.0*e) + f;
		}
	}
	else{
		if (s<0){
			//region 2
			tmp0 = b + d;
			tmp1 = c + e;
			if (tmp1 > tmp0){ // minimum on edge s+t=1
				numer = tmp1 - tmp0;
				denom = a - 2.0*b + c;
				if (numer >= denom){
					s = 1;
					t = 0;
					dist = a + 2.0*d + f;
				}
				else{
					s = numer / denom;
					t = 1.0 - s;
					dist = s*(a*s + b*t + 2.0*d)
						+ t*(b*s + c*t + 2 * e) + f;
				}
			}
			else{
				s = 0;
				if (tmp1 <= 0){
					t = 1;
					dist = c + 2.0*e + f;
				}
				else{
					if (e >= 0){
						t = 0;
						dist = f;
					}
					else{
						t = -1.0*e / c;
						dist = e*t + f;
					}
				}
			}
		}
		else if (t<0){
			//region 6

			tmp0 = b + e;
			tmp1 = a + d;
			if (tmp1 > tmp0){ // minimum on edge s+t=1
				numer = tmp1 - tmp0;
				denom = a - 2.0*b + c;
				if (numer >= denom){
					t = 1;
					s = 0;
					dist = c + 2.0*e + f;
				}
				else{
					t = numer / denom;
					s = 1.0 - t;
					dist = s*(a*s + b*t + 2.0*d)
						+ t*(b*s + c*t + 2.0*e) + f;
				}
			}
			else{
				t = 0;
				if (tmp1 <= 0){
					s = 1;
					dist = a + 2.0*d + f;
				}
				else{
					if (d >= 0){
						s = 0;
						dist = f;
					}
					else{
						s = -d / a;
						dist = d*s + f;
					}
				}
			}
		}
		else{
			//region 1
			numer = c + e - b - d;
			if (numer <= 0){
				s = 0;
				t = 1;
				dist = c + 2.0*e + f;
			}
			else{
				denom = a - 2.0*b + c;
				if (numer >= denom){
					s = 1;
					t = 0;
					dist = a + 2.0*d + f;
				}
				else{
					s = numer / denom;
					t = 1 - s;
					dist = s*(a*s + b*t + 2.0*d) +
						t*(b*s + c*t + 2.0*e) + f;
				}
			}
		}
	}

	x_new = xp1 + s*xE0 + t*xE1;
	y_new = yp1 + s*yE0 + t*yE1;
	z_new = zp1 + s*zE0 + t*zE1;
	if (dist<0){ dist = 0; }
	return dist;

}

template<typename T>
bool LinePlaneIntersect(T pp_x, T pp_y, T pp_z, T pv_x, T pv_y, T pv_z,
	                    T ip1_x, T ip1_y, T ip1_z, T ip2_x, T ip2_y, T ip2_z,
	                    T&point_x, T&point_y, T&point_z){
	//plane line intersection. plane define by normal vector (pv_x,pv_y,pv_z) and point on it(pp_x,pp_y,pp_z)
	// and line between point ip1 and ip2... return point (point_x,point_y,point_z)

	double ux, uy, uz; // line's vector

	ux = ip2_x - ip1_x;
	uy = ip2_y - ip1_y;
	uz = ip2_z - ip1_z;

	double dot = Dot(ux, uy, uz, pv_x, pv_y, pv_z);

	if (abs(dot) <= 0.0){
		return false;
	}

	double s;

	s = (Dot(pv_x, pv_y, pv_z, pp_x - ip1_x, pp_y - ip1_y, pp_z - ip1_z)) / (dot);


	if (s<-1.0*10E-12 || s>1.0 + 10E-12){
		return false;
	}
	point_x = ip1_x + s*ux;
	point_y = ip1_y + s*uy;
	point_z = ip1_z + s*uz;
	return true;

}

static bool GetInsct2Planes(double&point_x, double&point_y, double&point_z, double&vector_x, double&vector_y, double&vector_z, double** two_planes)
{
	//two_planes is an ainput pointer where first columne stores the vector normal to plane1 and a point on planes 1 
	// and second columne stores same stuff for plane 2 

	// double x_ip(_x_samples[_ip]), y_ip(_y_samples[_ip]), z_ip(_z_samples[_ip]), 
	//	x_ip1(_x_samples[_ip1]), y_ip1(_y_samples[_ip1]), z_ip1(_z_samples[_ip1]), 
	//	x_ip2(_x_samples[_ip2]), y_ip2(_y_samples[_ip2]), z_ip2(_z_samples[_ip2]); // all points

	double x_ip1, y_ip1, z_ip1, x_ip2, y_ip2, z_ip2, ip_ip1_x, ip_ip1_y, ip_ip1_z, ip_ip2_x, ip_ip2_y, ip_ip2_z;

	//double ip_ip1_x(x_ip1-x_ip), ip_ip1_y(y_ip1-y_ip), ip_ip1_z(z_ip1-z_ip); // vector _ip1 -> _ip 
	//double ip_ip2_x(x_ip2-x_ip), ip_ip2_y(y_ip2-y_ip), ip_ip2_z(z_ip2-z_ip); // vector _ip2 -> _ip

	ip_ip1_x = two_planes[0][0]; ip_ip1_y = two_planes[0][1]; ip_ip1_z = two_planes[0][2];
	ip_ip2_x = two_planes[1][0]; ip_ip2_y = two_planes[1][1]; ip_ip2_z = two_planes[1][2];

	x_ip1 = two_planes[0][3]; y_ip1 = two_planes[0][4]; z_ip1 = two_planes[0][5];
	x_ip2 = two_planes[1][3]; y_ip2 = two_planes[1][4]; z_ip2 = two_planes[1][5];

	double n1_cros_n2_x, n1_cros_n2_y, n1_cros_n2_z;
	Cross(ip_ip1_x, ip_ip1_y, ip_ip1_z, ip_ip2_x, ip_ip2_y, ip_ip2_z, n1_cros_n2_x, n1_cros_n2_y, n1_cros_n2_z);
	NormalizeVector(n1_cros_n2_x, n1_cros_n2_y, n1_cros_n2_z);

	if (abs(n1_cros_n2_x*n1_cros_n2_x + n1_cros_n2_y*n1_cros_n2_y + n1_cros_n2_z*n1_cros_n2_z) < _tol_sq){
		//parallel/concident planes
		return false;
		//std::cout << "Error (1) at GetInsct2Planes()..!!" << std::endl;		
	}

	double n1n1 = Dot(ip_ip1_x, ip_ip1_y, ip_ip1_z, ip_ip1_x, ip_ip1_y, ip_ip1_z); //_ip_ip1 (dot) _ip_ip1
	double n2n2 = Dot(ip_ip2_x, ip_ip2_y, ip_ip2_z, ip_ip2_x, ip_ip2_y, ip_ip2_z); //_ip_ip2 (dot) _ip_ip2
	double n1n2 = Dot(ip_ip1_x, ip_ip1_y, ip_ip1_z, ip_ip2_x, ip_ip2_y, ip_ip2_z); //_ip_ip1 (dot) _ip_ip2

	double determinant = n1n1*n2n2 - n1n2*n1n2;

	//double d1=Dot(ip_ip1_x, ip_ip1_y,ip_ip1_z, _x_sphere,_y_sphere,_z_sphere); // change 
	//double d2=Dot(ip_ip2_x, ip_ip2_y,ip_ip2_z, _x_sphere,_y_sphere,_z_sphere); //change 

	double d1 = Dot(ip_ip1_x, ip_ip1_y, ip_ip1_z, x_ip1, y_ip1, z_ip1);
	double d2 = Dot(ip_ip2_x, ip_ip2_y, ip_ip2_z, x_ip2, y_ip2, z_ip2);


	double c1 = (d1*n2n2 - d2*n1n2) / determinant;
	double c2 = (d2*n1n1 - d1*n1n2) / determinant;

	point_x = c1*ip_ip1_x + c2*ip_ip2_x;
	point_y = c1*ip_ip1_y + c2*ip_ip2_y;
	point_z = c1*ip_ip1_z + c2*ip_ip2_z;

	vector_x = n1_cros_n2_x;
	vector_y = n1_cros_n2_y;
	vector_z = n1_cros_n2_z;

	double nn = sqrt(vector_x*vector_x + vector_y*vector_y + vector_z*vector_z);
	vector_x /= nn; vector_y /= nn; vector_z /= nn;


	if (false){
		//DrawPlane(ip_ip1_x, ip_ip1_y, ip_ip1_z, x_ip1, y_ip1, z_ip1);
		//DrawPlane(ip_ip2_x, ip_ip2_y, ip_ip2_z, x_ip2, y_ip2, z_ip2);
	}

	return true;
	/*double checkx, checky, checkz;
	checkx=(_x_sphere-point_x)/vector_x;
	checky=(_y_sphere-point_y)/vector_y;
	checkz=(_z_sphere-point_z)/vector_z;

	/*if(abs(checkx-checky)>10E-6 || abs(checky-checkz)>10E-6 || abs(checkz-checkx)>10E-6){ //^^debug
	cout<<"Error(2) at GetInsct2Planes()..!! " <<endl;
	system("pause");
	}*/


}
static bool GetInsct2Planes(double&point_x, double&point_y, double&point_z, double&vector_x, double&vector_y, double&vector_z, double*C1, double*C2)
{
	//combine the two planes C1 and C2 in a double** array and call GetInsct2Planes
	double **two_planes = new double*[2];
	two_planes[0] = new double[6];
	two_planes[1] = new double[6];

	two_planes[0][3] = C1[0];//norm to plane 0
	two_planes[0][4] = C1[1];
	two_planes[0][5] = C1[2];
	two_planes[0][0] = C1[3];//point on plane 0
	two_planes[0][1] = C1[4];
	two_planes[0][2] = C1[5];


	two_planes[1][3] = C2[0];//norm to plane 1
	two_planes[1][4] = C2[1];
	two_planes[1][5] = C2[2];
	two_planes[1][0] = C2[3];//point on plane 1
	two_planes[1][1] = C2[4];
	two_planes[1][2] = C2[5];

	bool isInsect = GetInsct2Planes(point_x, point_y, point_z, vector_x, vector_y, vector_z, two_planes);

	delete two_planes[0];
	delete two_planes[1];
	delete two_planes;

	return isInsect;
}
inline bool SolveQuadEqu(double a, double b, double c, double&u1, double&u2)
{
	double delta = b*b - 4.0*a*c;
	if (delta<0.0){
		return false;
	}


	u1 = (-1.0*b + sqrt(delta)) / (2.0*a);
	u2 = (-1.0*b - sqrt(delta)) / (2.0*a);
	return true;

}

template<typename T>
inline void BarycentricLite(T x1, T y1, T z1, T x2, T y2, T z2, T x3, T y3, T z3, T x_p, T y_p, T z_p, T&alfa, T& beta, T& gamma) {

	double xx(0), yy(0), zz(0);
	
	Cross(x3 - x1, x2 - x1, x1 - x_p, 
	      y3 - y1, y2 - y1, y1 - y_p,
		  xx,yy,zz);

	Cross(y3 - y1, y2 - y1, y1 - y_p,
		  z3 - z1, z2 - z1, z1 - z_p,
		  xx, yy, zz);

	Cross(z3 - z1, z2 - z1, z1 - z_p,
		  x3 - x1, x2 - x1, x1 - x_p,
		  xx, yy, zz);
	
	if (std::abs(zz) < 1 - _tol) {
		std::cout << "Error (0) at geo::Barycentric(). Triangle is degenerate!!!" << std::endl;
		system("pause");		
		//return Vec3f(-1, 1, 1);
	}
	alfa = 1.0 - (xx + yy) / zz;
	beta = yy / zz;
	gamma = xx / zz;
	
}


template<typename T>
inline void Barycentric(T x1, T y1, T z1, T x2, T y2, T z2, T x3, T y3, T z3, T x_p, T y_p, T z_p, T&alfa, T& beta, T& gamma, int stack = 0)
{

	T u21x(x2 - x1), u21y(y2 - y1), v31x(x3 - x1), v31y(y3 - y1), wq1x(x_p - x1), wq1y(y_p - y1);
	if (u21y == 0 || u21x == 0 || v31y - v31x == 0){
		u21x = x1 - x2;
		u21y = y1 - y2;
		v31x = x3 - x2;
		v31y = y3 - y2;
		wq1x = x_p - x2;
		wq1y = y_p - y2;
	}
	if (u21y == 0 || u21x == 0 || v31y - v31x == 0){
		u21x = x2 - x3;
		u21y = y2 - y3;
		v31x = x1 - x3;
		v31y = y1 - y3;
		wq1x = x_p - x3;
		wq1y = y_p - y3;
	}
	if (u21x == 0 || (v31y - v31x*(u21y / u21x)) == 0){
		stack++;
		if (stack > 4){
			std::cout << "Error (0) at geo::Barycentric()... " << std::endl;
			system("pause");
		}
		Barycentric(z1, x1, y1, z2, x2, y2, z3, x3, y3, x_p, y_p, z_p, alfa, beta, gamma, stack);
		return;
		
	}
	if (!(u21x == 0)){
		beta = (wq1y - wq1x*(u21y / u21x)) / (v31y - v31x*(u21y / u21x));
		alfa = (wq1x - beta*v31x) / u21x;
	}
	else
	{
		if (v31x == 0){
			stack++;
			if (stack > 4){
				std::cout << "Error (0) at geo::Barycentric()... " << std::endl;
				system("pause");
			}
			Barycentric(z1, x1, y1, z2, x2, y2, z3, x3, y3, x_p, y_p, z_p, alfa, beta, gamma, stack);
			return;			
		}
		alfa = (wq1y - wq1x*(v31y / v31x)) / (u21y - v31y*(u21x / v31x));
		beta = (wq1x - alfa*u21x) / v31x;
	}
	gamma = 1.0 - alfa - beta;
}

template<typename T>
bool removeFromVector(T item, std::vector<T>&myVec)
{
	for (int i = 0; i < int(myVec.size()); i++){
		if (myVec[i] == item){
			myVec[i] = myVec.back();
			myVec.pop_back();
			return true;
		}
	}
	return false;
}

inline bool SphereLineIntersection(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z,
	                               double cx, double cy, double cz, double r_2,
	                               double&px1, double&py1, double&pz1,
	                               double&px2, double&py2, double&pz2,
	                               size_t&num_sect)
{
	//http://wiki.cgsociety.org/index.php/Ray_Sphere_Intersection
	//line(p2-p1) sphere(c) intersection
	//true if there is intersection and return (p) as the insect point 
	//return false otherwise

	double a, b, c, t1, t2;
	//a-->t^2
	//b-->t^1
	//c-->t^0
	a = (p2x - p1x)*(p2x - p1x) + (p2y - p1y)*(p2y - p1y) + (p2z - p1z)*(p2z - p1z);
	b = 2.0*(p1x - cx)*(p2x - p1x) + 2.0*(p1y - cy)*(p2y - p1y) + 2.0*(p1z - cz)*(p2z - p1z);
	c = (p1x - cx)*(p1x - cx) + (p1y - cy)*(p1y - cy) + (p1z - cz)*(p1z - cz) - r_2;

	if (!SolveQuadEqu(a, b, c, t1, t2)){ return false; }


	if (10E-8 <= t1 && t1 <= 1.0 - 10E-8 && 10E-8 <= t2 && t2 <= 1.0 - 10E-8){
		px1 = p1x + t1*(p2x - p1x);
		py1 = p1y + t1*(p2y - p1y);
		pz1 = p1z + t1*(p2z - p1z);

		px2 = p1x + t2*(p2x - p1x);
		py2 = p1y + t2*(p2y - p1y);
		pz2 = p1z + t2*(p2z - p1z);

		num_sect = 2;
		return true;
	}

	else if (10E-8 <= t1 && t1 <= 1.0 - 10E-8){
		px1 = p1x + t1*(p2x - p1x);
		py1 = p1y + t1*(p2y - p1y);
		pz1 = p1z + t1*(p2z - p1z);

		num_sect = 1;
		return true;
	}
	else if (10E-8 <= t2 && t2 <= 1.0 - 10E-8){
		px1 = p1x + t2*(p2x - p1x);
		py1 = p1y + t2*(p2y - p1y);
		pz1 = p1z + t2*(p2z - p1z);
		num_sect = 1;
		return true;
	}
	else{
		return false;
	}

}
template<typename T>
inline void Rotate3D(T x, T y, T z, T theta, size_t axis, T&xx, T&yy, T&zz)
{
	T sn, cs;
	sn = sin(theta*PI / 180.0);
	cs = cos(theta*PI / 180.0);
	if (axis == 0){
		xx = x;
		yy = y*cs - z*sn;
		zz = y*sn + z*cs;

	}
	else if (axis == 1){
		xx = x*cs + z*sn;
		yy = y;
		zz = -x*sn + z*cs;
	}
	else if (axis == 2){
		xx = x*cs - y*sn;
		yy = x*sn + y*cs;
		zz = z;
	}
	else{
		std::cout << "Error at Rotate3D(). Wrong axis entry." << std::endl;
		system("pause");
	}
}
#endif /*_COMMON_*/