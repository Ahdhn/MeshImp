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

#ifndef _OBJREADER_
#define _OBJREADER_

#include <stdlib.h>  
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

/*
	Read the vertices and triangles from a input obj file
	WARNING: read would return garbbage if the texture is included in the obj file
	Two passes; first one to count number of verices and faces, second to read them
*/
void objReader(char*filename, int&numVert, double**&Verts, int&numTri, int**&Tris)
{
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "*************************** Reading OBJ File ************************************");
	fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	fprintf(stdout, "\nReading %s\n", filename);

	std::ifstream myObjFile;
	myObjFile.open(filename, std::ifstream::in);
	std::string lineStr;

	numVert = 0;
	numTri = 0;

	if (myObjFile.is_open()){
		//pass one 
		while (std::getline(myObjFile, lineStr)){
			if (lineStr[0] == 'V' || lineStr[0] == 'v'){
				numVert++;
			}
			else if (lineStr[0] == 'F' || lineStr[0] == 'f'){
				numTri++;
			}
		}
		Verts = new double*[numVert];
		Tris = new int *[numTri];
		myObjFile.close();
		myObjFile.open(filename, std::ifstream::in);
		//pass two 		
		int vert_read(0), tri_read(0);

		char myChar[10];
		while (!myObjFile.eof()){
			myObjFile >> myChar;
			if (strcmp(myChar, "V") == 0 || strcmp(myChar, "v") == 0){
				Verts[vert_read] = new double[3];
				myObjFile >> Verts[vert_read][0] >> Verts[vert_read][1] >> Verts[vert_read][2];
				vert_read++;
			}
			else if (strcmp(myChar, "F") == 0 || strcmp(myChar, "f") == 0){
				Tris[tri_read] = new int[3];
				myObjFile >> Tris[tri_read][0] >> Tris[tri_read][1] >> Tris[tri_read][2];
				Tris[tri_read][0]--;
				Tris[tri_read][1]--;
				Tris[tri_read][2]--;
				tri_read++;
			}
			else {
				std::getline(myObjFile, lineStr);//read the rest of the line 
			}


		}

		if (vert_read != numVert || tri_read != numTri){
			std::cerr << "Error at objReader::Error in read " << filename << std::endl;
			exit(1);
		}
	}
	else{
		std::cerr << "Error at objReader::Can not open file " << filename << std::endl;
		exit(1);
	}
}



#endif /*_OBJREADER_*/