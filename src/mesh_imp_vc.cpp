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

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <limits>

#include "util/OBJ.h"
#include "util/CSV.h"
#include "operator/execute.h"
#include "util\Common.h"


void PrintHelpMessage()
{
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "************************************USE MESSAGE**********************************");
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "Command Line Syntax is: \n\n");
	fprintf(stdout, "meshimp -APP [-tar -sam -smooth -dih -ring -del -minang -maxang \n");
	fprintf(stdout, "-minedge - maxedge] INPUT.obj\n");	
	fprintf(stdout, "\n-APP could be `-sim' for Mesh Simplification or `-obt' for Non-obtuse \n");
	fprintf(stdout, "Retriangulation.\n");
	

	fprintf(stdout, "\n*************************************************************\n");
	fprintf(stdout, "Command Line Switches are:\n\n");
	fprintf(stdout, "   -sim      Invokes the Mesh Simplification algorithm on the input mesh\n");
	fprintf(stdout, "             where the complexity of the input mesh is reduced while \n");
	fprintf(stdout, "             respecting the specified mesh quality.\n");

	fprintf(stdout, "   -tar      The number followed represents the target number of samples desired\n");
	fprintf(stdout, "             in the simplified output mesh. If this number is not set, \n");
	fprintf(stdout, "             The code will attempt to reduce the samples count as much as possible.\n");

	fprintf(stdout, "   -obt      Invokes the non-obtuse retriangulation algorithm on the input\n");
	fprintf(stdout, "             mesh where we attempt to eliminate the obtuse triangles from \n");
	fprintf(stdout, "             the input mesh.\n");

	fprintf(stdout, "   -sam      The number followed represents the number of successive successful\n");
	fprintf(stdout, "             darts before termination.\n");
	fprintf(stdout, "             During the sampling process, we can sample more than one successful\n");
	fprintf(stdout, "             sample and pick the best one. The best one can have different objective.\n");
	fprintf(stdout, "             In the current implementation, the `best' here is the one with maximum.\n");
	fprintf(stdout, "             minimum apex angle across different algorihtms.\n");
	fprintf(stdout, "             The default is 10 samples.\n");

	fprintf(stdout, "   -ring     The number followed represents the number of rings (layers)\n");
	fprintf(stdout, "             taken as background grid for sampling. The default is 3 rings.\n");

	fprintf(stdout, "   -del      If set, the Delaunay property will be preserved.\n");
	fprintf(stdout, "             The default is false.\n");

	fprintf(stdout, "   -minang   The number followed represents the minimum angle preserved.\n");
	fprintf(stdout, "             The default is set to the minimum angle of the input mesh.\n");
	fprintf(stdout, "   -maxang   The number followed represents the maximum angle preserved.\n");
	fprintf(stdout, "             The default is set to the maximum angle of the input mesh.\n");

	//	fprintf(stdout, "   -minedge  The number followed represents the minimum edge preserved.\n");
	//	fprintf(stdout, "   -maxedge  The number followed represents the maximum edge preserved.\n");

	fprintf(stdout, "   -smooth   If set, the smoothness will be preserved.\n");
	fprintf(stdout, "             The default is false.\n");
	fprintf(stdout, "   -dih      The maximum dihedral angle to which the smoothness will be\n");
	fprintf(stdout, "             preserved.\n");
	fprintf(stdout, "   -v        Display various statistics throughout the execution.\n");
	fprintf(stdout, "   -h        Display the use message, and quit.\n");
	
}

void ScaleSphereAndTriangles(int numSpheres, double**Spheres, int numVert, double**Verts, double&xBase, double&yBase, double&zBase, double&scaleFactor){
	xBase = yBase = zBase = DBL_MAX;
	scaleFactor = 1.0;
	double xMax(-DBL_MAX), yMax(-DBL_MAX), zMax(-DBL_MAX);
	for (int i = 0; i < numSpheres; i++){
		xBase = std::min(xBase, Spheres[i][0]);
		yBase = std::min(yBase, Spheres[i][1]);
		zBase = std::min(zBase, Spheres[i][2]);

		xMax = std::max(xMax, Spheres[i][0]);
		yMax = std::max(yMax, Spheres[i][1]);
		zMax = std::max(zMax, Spheres[i][2]);
	}

	for (int i = 0; i < numVert; i++){
		xBase = std::min(xBase, Verts[i][0]);
		yBase = std::min(yBase, Verts[i][1]);
		zBase = std::min(zBase, Verts[i][2]);

		xMax = std::max(xMax, Verts[i][0]);
		yMax = std::max(yMax, Verts[i][1]);
		zMax = std::max(zMax, Verts[i][2]);
	}

	scaleFactor = std::max(xMax - xBase, yMax - yBase);
	scaleFactor = std::max(scaleFactor, zMax - zBase);

	for (int i = 0; i < numSpheres; i++){
		Spheres[i][0] -= xBase;
		Spheres[i][1] -= yBase;
		Spheres[i][2] -= zBase;

		Spheres[i][0] /= scaleFactor;
		Spheres[i][1] /= scaleFactor;
		Spheres[i][2] /= scaleFactor;
		Spheres[i][3] /= scaleFactor;
	}

	for (int i = 0; i < numVert; i++){
		Verts[i][0] -= xBase;
		Verts[i][1] -= yBase;
		Verts[i][2] -= zBase;

		Verts[i][0] /= scaleFactor;
		Verts[i][1] /= scaleFactor;
		Verts[i][2] /= scaleFactor;
	}

}
void RotateSphereAndTriangles(int numSpheres, double**Spheres, int numVert, double**Verts, const double xAngle, const double  yAngle, const double zAngle){

	double xx(0), yy(0), zz(0);
	for (int i = 0; i < numSpheres; i++){
		Rotate3D(Spheres[i][0], Spheres[i][1], Spheres[i][2], xAngle, 0, xx, yy, zz);
		Spheres[i][0] = xx; Spheres[i][1] = yy; Spheres[i][2] = zz;

		Rotate3D(Spheres[i][0], Spheres[i][1], Spheres[i][2], yAngle, 1, xx, yy, zz);
		Spheres[i][0] = xx; Spheres[i][1] = yy; Spheres[i][2] = zz;

		Rotate3D(Spheres[i][0], Spheres[i][1], Spheres[i][2], zAngle, 2, xx, yy, zz);
		Spheres[i][0] = xx; Spheres[i][1] = yy; Spheres[i][2] = zz;

	}

	for (int i = 0; i < numVert; i++){
		Rotate3D(Verts[i][0], Verts[i][1], Verts[i][2], xAngle, 0, xx, yy, zz);
		Verts[i][0] = xx; Verts[i][1] = yy; Verts[i][2] = zz;

		Rotate3D(Verts[i][0], Verts[i][1], Verts[i][2], yAngle, 1, xx, yy, zz);
		Verts[i][0] = xx; Verts[i][1] = yy; Verts[i][2] = zz;

		Rotate3D(Verts[i][0], Verts[i][1], Verts[i][2], zAngle, 2, xx, yy, zz);
		Verts[i][0] = xx; Verts[i][1] = yy; Verts[i][2] = zz;
	}
}
int main(int argc, char**argv)
{
	
	///char* filename = NULL;
	//char* OBJfilename = "input/fertility_vc_surface.obj";
	//char* CSVfilename = "input/fertility_surface_spheres.csv"; bool withTag = false;

	char* OBJfilename = "input/cube_vc_surface.obj";
	char* CSVfilename = "input/cube_surface_spheres_tagged.csv"; bool withTag = true;

	//char* OBJfilename = "input/scorpion_vc_surface.obj";
	//char* CSVfilename = "input/scorpion_surface_spheres_tagged.csv"; bool withTag = true;
	

	//ParseInput(argc, argv, mySwitches, filename);	

	//2) Read input mesh and spheres and build initial data structure 
	int numVert(0), numTri(0), numSpheres(0);
	double**Verts = nullptr; int**Tris = nullptr; double**Spheres = nullptr;
	
	cvsReader(CSVfilename, numSpheres, Spheres, withTag);
	objReader(OBJfilename, numVert, Verts, numTri, Tris);
	
	//3)scaling 
	double xBase(DBL_MAX), yBase(DBL_MAX), zBase(DBL_MAX), scaleFactor(1);
	ScaleSphereAndTriangles(numSpheres, Spheres, numVert, Verts, xBase, yBase, zBase, scaleFactor);
	double xAngle((double(rand()) / double(RAND_MAX))*180.0), yAngle((double(rand()) / double(RAND_MAX))*180.0), zAngle((double(rand()) / double(RAND_MAX))*180.0);
	RotateSphereAndTriangles(numSpheres, Spheres, numVert, Verts,xAngle,yAngle,zAngle);

		
	//4) Call the right application
	MeshImp myImp(numSpheres, Spheres, numVert, Verts, numTri, Tris);
	myImp.VC(-1, 10, 5, 1, 1);


		

	return 0;
}
