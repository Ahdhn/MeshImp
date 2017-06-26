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

#ifndef _STAT_
#define _STAT_


#include <string.h>
#include <stdio.h>
#include <string>

#include "Common.h"

//TODO get diiferent stats at single vertex 
//TODO report valence 

class Statistics
{
public:
	Statistics();
	~Statistics();

	//First call this to calc all statistics
	void GetAllStats(int numVert, vert *Vert);

	//Quary functions for the calculated statistics 
	double GetMinAngle(){ return min_angle; }
	double GetMaxAngle(){ return max_angle; }
	void GetMinMaxEdgeLen(double&myMinEdgeLen, double&myMaxEdgeLen){ myMinEdgeLen = min_edge_len; myMaxEdgeLen = max_edge_len; }
	void GetMinMaxTriQuality(double&myMinQuality, double&myMaxQuality){ myMinQuality = min_quality; myMaxQuality = max_quality; }
	int GetNumTri(){ return num_tri; }
	int GetNumNonObtuse(){ return num_tri_obtuse; }
	int CalcNumAcute(int numVert, vert *Vert, double measureAngle);
	void DisplayStatistics(FILE *fp, bool obt);

private:
	double min_dih_ang, max_dih_ang,
		   min_quality, max_quality,
		   min_angle, max_angle,
		   min_edge_len, max_edge_len,
		   average_div_from_60; //average deviation of interior angles from 60 degree 
	
	int num_tri, num_tri_obtuse, num_vertices;
	
	
	double getCurv(int, int, vert *);
	void CalcDihAngle(int, vert *);
	void CalcAnglesEdgeLen(int, vert *);
	void CalcTriangelQuality(int, vert *);
	void CalcNumTriangles(int, vert *);
	void CalcNumNonObtuse(int, vert *);
	double Statistics::SingleTriangleQuality(int, int, int, vert *);	
	void ErrWarnMessage(size_t lineNum, std::string message, size_t mess_id);
};





#endif /*_STAT_*/