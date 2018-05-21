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

#ifndef _CSVREADER_
#define _CSVREADER_

#include <stdlib.h>  
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

/*
	Read shperes (X,Y,Z,R)
	WARNING: read would return garbbage if the texture is included in the obj file
	Two passes; first one to count number of verices and faces, second to read them
*/
unsigned fileNumLines(char*filename)
{
	//https://stackoverflow.com/questions/3482064/counting-the-number-of-lines-in-a-text-file

	std::ifstream myfile(filename);
	// new lines will be skipped unless we stop it from happening:    
	myfile.unsetf(std::ios_base::skipws);
	// count the newlines with an algorithm specialized for counting:
	unsigned line_count = std::count(
		std::istream_iterator<char>(myfile),
		std::istream_iterator<char>(),
		'\n');
	myfile.close();

	return line_count;
}
void cvsReader(char*filename, int&numSpheres, double**&Spheres, bool isTagged = false)
{
	
	fprintf(stdout, "\nReading %s\n", filename);

	
	std::ifstream myObjFile;
	myObjFile.open(filename, std::ifstream::in);
	std::string lineStr;

	numSpheres = fileNumLines(filename);
	fprintf(stdout, " NumSpheres = %d\n", numSpheres);

	int num_current = 0;
	if (myObjFile.is_open()){			
		if (num_current == 0){ std::getline(myObjFile, lineStr); }
		Spheres = new double*[numSpheres];
		char myChar;

		while (true){
			
			Spheres[num_current] = new double[5];
			
			myObjFile >> Spheres[num_current][0] >> myChar;
			myObjFile >> Spheres[num_current][1] >> myChar;
			myObjFile >> Spheres[num_current][2] >> myChar;
			if (isTagged){
				myObjFile >> Spheres[num_current][3] >> myChar;
				myObjFile >> Spheres[num_current][4];
			}
			else{
				myObjFile >> Spheres[num_current][3];
				Spheres[num_current][4] = 0;
			}			
			num_current++;
			if (num_current >= numSpheres){ break; }
		}

	}
	else{
		std::cerr << "Error at cvsReader::Can not open file " << filename << std::endl;
		exit(1);
	}

	myObjFile.close();
}

#endif /*_CSVREADER_*/