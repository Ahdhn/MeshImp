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
struct tri{
	int id[3];
	int neighbour[3];
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


#endif /*_COMMON_*/