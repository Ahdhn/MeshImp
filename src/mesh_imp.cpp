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

#include "util/OBJ.h"
#include "operator/execute.h"


struct Switches
{
	//Input switches 
	//simp:   -sim      
	//nonobt: -nonobt  
	//samplingBudget: -samples 
	//numLayer: -ring 
	//isDelaunay: -del 

	//isSmooth: -smooth
	//devFactor: -dih
	
	//minAngle: -minang 
	//maxAngle: -maxang 

	//minEdge: -minedge 
	//maxEdge: -maxedge 

	bool simp, nonobt, isSmooth, isDelaunay, verbose;
	int  samplingBudget, numRing, targetNumSamples;
	double dih, minAngle, maxAngle, minEdge, maxEdge;	
};

void PrintHelpMessage()
{
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "************************************USE MESSAGE**********************************");
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "Command Line Syntax is: \n\n");
	fprintf(stdout, "mesh_imp.exe -APP [-tar -sam -smooth -dih -ring -del -minang -maxang \n");
	fprintf(stdout, "-minedge - maxedge] INPUT.obj\n");	
	fprintf(stdout, "\n-APP could be `-sim' for Mesh Simplification or `-obt' for Non-obtuse \n");
	fprintf(stdout, "Retriangulation.\n");
	

	fprintf(stdout, "\n*************************************************************\n");
	fprintf(stdout, "Command Line Switches are:\n\n");
	fprintf(stdout, "   -sim      Invokes the Mesh Simplification algorithm on the input mesh\n");
	fprintf(stdout, "             where the complexity of the input mesh is reduced while \n");
	fprintf(stdout, "             respecting the specified mesh quality.\n");

	fprintf(stdout, "   -tar      The number followed represents the target number of samples desired\n");
	fprintf(stdout, "             of the simplified in the output mesh. If this number is not set, \n");	
	fprintf(stdout, "             The code will attempt to reduce the samples count as much as possible.\n");
	
	fprintf(stdout, "   -obt      Invokes the non-obtuse retriangulation algorithm on the input\n");
	fprintf(stdout, "             mesh where we attempt to eliminate the obtuse triangles from \n");
	fprintf(stdout, "             the input mesh.\n");

	fprintf(stdout, "   -sam      The number followed represents the number of successful samples.\n");
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
	fprintf(stdout, "   -h        Display the use message.\n");
	
}

bool ParseInteger(int startingIndex, char*searchString, int&IntegerVal)
{
	//parse the integer from searchString starting from startingIndex
	//return false if it is not integer 
	//return true if it is integer and store the value in IntegerVal

	char numString[2048];

	IntegerVal = 0;
	size_t myStrLen = strlen(searchString);
	int k = 0;

	for (size_t i = startingIndex; i < myStrLen; i++){
		if (((searchString[i] >= '0') && (searchString[i] <= '9'))){
			numString[k] = searchString[i];
			k++;
			startingIndex++;
		}
		else{
			return false;
		}
	}

	IntegerVal = strtoul(numString, (char**)NULL, 0);

	return true;

}
bool ParseDouble(int startingIndex, char*searchString, double&DoubleVal)
{
	//parse the integer from searchString starting from startingIndex
	//return false if it is not integer 
	//return true if it is integer and store the value in IntegerVal
	char numString[2048];

	DoubleVal = 0;
	size_t myStrLen = strlen(searchString);
	int k = 0;

	for (size_t i = startingIndex; i < myStrLen; i++){
		if (((searchString[i] >= '0') && (searchString[i] <= '9')) || (i == startingIndex && searchString[i] == '.')){
			numString[k] = searchString[i];
			k++;
			startingIndex++;
		}
		else {
			return false;
		}
	}

	DoubleVal = strtod(numString, (char**)NULL);
		
	return true;
}
void ParseInput(int argc, char**argv, struct Switches *mySwitches, char*&filename )
{
	mySwitches->dih = -170.0;
	mySwitches->isDelaunay = false;
	mySwitches->isSmooth = false;
	mySwitches->minAngle = -1.0;
	mySwitches->maxAngle = -1.0;
	mySwitches->nonobt = false;
	mySwitches->simp = false;
	mySwitches->numRing = 3;
	mySwitches->samplingBudget = 10;
	mySwitches->verbose = false;
	mySwitches->targetNumSamples = -1;

	
	//filename = NULL;

	if (argc < 2){
		fprintf(stderr, "\nERROR:Invalid input. No switches or input file are provided!!!\n");
		PrintHelpMessage();
		exit(1);
	}

	for (int i = 1; i < argc; i++){
		//std::cout << argv[i] << std::endl;
		if (argv[i][0] == '-'){

			//Help
			if (argv[i][1] == 'h' && argv[i][2] == '\0'){
				PrintHelpMessage();		
				exit(0);
			}
			
			//Simplification
			else if (argv[i][1] == 's' && argv[i][2] == 'i'&& argv[i][3] == 'm'&& argv[i][4] == '\0'){
				mySwitches->simp = true;
			}	

			//Target number of samples
			else if (argv[i][1] == 't' && argv[i][2] == 'a'&& argv[i][3] == 'r'){
				if (!ParseInteger(4, argv[i], mySwitches->targetNumSamples)){
					fprintf(stderr, "\nERROR:Invalid input. Target number of samples should be an integer!!!\n");
					PrintHelpMessage();
					exit(1);
				}
			}

			//Non-obtuse
			else if (argv[i][1] == 'o' && argv[i][2] == 'b'&& argv[i][3] == 't'&& argv[i][4] == '\0'){
				mySwitches->nonobt = true;
			}

			//isDelaunay
			else if (argv[i][1] == 'd' && argv[i][2] == 'e' && argv[i][3] == 'l'&& argv[i][4] == '\0'){
				mySwitches->isDelaunay = true;
			}

			//smoothness 
			else if (argv[i][1] == 's' && argv[i][2] == 'm' && argv[i][3] == 'o'&& argv[i][4] == 'o' &&
				argv[i][5] == 't'&& argv[i][6] == 'h' && argv[i][7] == '\0'){
				mySwitches->isSmooth = true;
			}

			//verbose
			else if (argv[i][1] == 'v' && argv[i][2] == '\0'){
				mySwitches->verbose = true;
			}

			//numRing 
			else if (argv[i][1] == 'r' && argv[i][2] == 'i' && argv[i][3] == 'n' && argv[i][4] == 'g'){
				if (!ParseInteger(5, argv[i], mySwitches->numRing)){
					fprintf(stderr, "\nERROR:Invalid input. Number of rings should be an integer!!!\n");
					PrintHelpMessage();
					exit(1);
				}
			}

			//number of samples
			else if (argv[i][1] == 's' && argv[i][2] == 'a' && argv[i][3] == 'm'){
				if (!ParseInteger(4, argv[i], mySwitches->samplingBudget)){
					fprintf(stderr, "\nERROR:Invalid input. Number of samples should be an integer!!!\n");
					PrintHelpMessage();
					exit(1);
				}
			}

			//dih angle 
			else if (argv[i][1] == 'd'&&argv[i][2] == 'i'&&argv[i][3] == 'h'){
				if (!ParseDouble(4, argv[i], mySwitches->dih)){
					fprintf(stderr, "\nERROR:Invalid input. Dihedral angle should be double!!!\n");
					PrintHelpMessage();
					exit(1);
				}
			}
			
			//min angle 
			else if (argv[i][1] == 'm'&&argv[i][2] == 'i'&&argv[i][3] == 'n' &&
				argv[i][4] == 'a'&&argv[i][5] == 'n'&&argv[i][6] == 'g'){

				if (strlen(argv[i]) == 7){
					//the default
					mySwitches->minAngle = -1.0;
				}
				else{
					if (!ParseDouble(7, argv[i], mySwitches->minAngle)){
						fprintf(stderr, "\nERROR:Invalid input. Minimum angle should be double!!!\n");
						PrintHelpMessage();
						exit(1);
					}
				}
			}

			//max angle 
			else if (argv[i][1] == 'm'&&argv[i][2] == 'a'&&argv[i][3] == 'x' &&
				argv[i][4] == 'a'&&argv[i][5] == 'n'&&argv[i][6] == 'g'){

				if (strlen(argv[i]) == 7){
					//the default
					mySwitches->maxAngle = -1.0;
				}
				else{
					if (!ParseDouble(7, argv[i], mySwitches->maxAngle)){
						fprintf(stderr, "\nERROR:Invalid input. Maximum angle should be double!!!\n");
						PrintHelpMessage();
						exit(1);
					}
				}
			}

		} else {
			//this should be the input file			
			for (size_t j = 0; j < strlen(argv[i]); j++){
				if (argv[i][j] == '.'){
					if (!((argv[i][j + 1] == 'o' || argv[i][j + 1] == 'O') &&
						  (argv[i][j + 2] == 'b' || argv[i][j + 2] == 'B') &&
						  (argv[i][j + 3] == 'j' || argv[i][j + 3] == 'J'))){
						fprintf(stderr, "\nERROR:Invalid input file!!!\n");
						PrintHelpMessage();
						exit(1);
					}
					break;
				}
			}			
			filename = new char[5000];
			filename = argv[i];
		}
	}

	if (filename == NULL){
		fprintf(stderr, "\nERROR:No input file provided!!!\n");
		PrintHelpMessage();
		exit(1);
	}

	if ((!mySwitches->nonobt && !mySwitches->simp) || (mySwitches->nonobt && mySwitches->simp)){
		fprintf(stderr, "\nERROR:Invalid input!!!\n");
		PrintHelpMessage();
		exit(1);
	}
}

int main(int argc, char**argv)
{
	//1) Parse Input command
	struct Switches *mySwitches =  new Switches;	
	char* filename = NULL;
	ParseInput(argc, argv, mySwitches, filename);
	

	//2) Read input mesh and build initial data structure 
	int numVert(0), numTri(0);
	double**Verts = NULL; 
	int**Tris = NULL; 	
	
	objReader(filename, numVert, Verts, numTri, Tris);
		
	//3) Call the right application
	MeshImp myImp(numVert, Verts, numTri, Tris);
	
	if (mySwitches->nonobt){
		
		//For acute >85
		/*myImp.AcuteRemeshing_InterleaveOpt(mySwitches->samplingBudget,
										   mySwitches->numRing,
										   mySwitches->isSmooth,
										   mySwitches->dih,
										   mySwitches->isDelaunay,
										   mySwitches->minAngle,
										   mySwitches->maxAngle,
										   mySwitches->verbose);*/

		//For removal of tiny small angles <20
		/*myImp.SmallAngleElimination_InterleaveOpt(mySwitches->minAngle,
												  mySwitches->samplingBudget, 
												  mySwitches->numRing,
												  mySwitches->isSmooth,
												  mySwitches->dih,
												  mySwitches->isDelaunay,
												  mySwitches->maxAngle,
												  mySwitches->verbose);*/

		//For typical non-obtuse remeshing
		myImp.NonobtuseRemeshing(mySwitches->samplingBudget,
			                     mySwitches->numRing, 
								 mySwitches->isSmooth, 
								 mySwitches->dih, 
								 mySwitches->isDelaunay,
								 mySwitches->minAngle,
								 mySwitches->verbose);

		//For improved non-obtuse remeshing
		/*myImp.NonobtuseRemeshing_InterleaveOpt(mySwitches->samplingBudget, 
			                                   mySwitches->numRing, 
											   mySwitches->isSmooth, 
											   mySwitches->dih, 
											   mySwitches->isDelaunay,
											   mySwitches->minAngle,
											   mySwitches->verbose);*/
	}
	else if (mySwitches->simp){
		myImp.Simp(mySwitches->targetNumSamples,
				   mySwitches->samplingBudget,
				   mySwitches->numRing,
				   mySwitches->isSmooth,
				   mySwitches->dih,
				   mySwitches->isDelaunay,
				   mySwitches->minAngle,
				   mySwitches->maxAngle,
				   mySwitches->verbose);
	}

	return 0;
}
