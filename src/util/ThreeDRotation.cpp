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

#include "ThreeDRotation.h"

#include <iostream>
#include <cmath> 
#define _USE_MATH_DEFINES
#include <math.h>
 
using namespace std; 
 
 
double rotationMatrix[4][4];
double inputMatrix[4][1] = {0.0, 0.0, 0.0, 0.0};
double outputMatrix[4][1] = {0.0, 0.0, 0.0, 0.0}; 
 
void showPoint(){
    cout<<"("<<outputMatrix[0][0]<<","<<outputMatrix[1][0]<<","<<outputMatrix[2][0]<<")"<<endl;
}  
void multiplyMatrix()
{
    for(int i = 0; i < 4; i++ ){
        for(int j = 0; j < 1; j++){
            outputMatrix[i][j] = 0;
            for(int k = 0; k < 4; k++){
                outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
            }
        }
    }
}
void setUpRotationMatrix(double angle, double u, double v, double w, double a, double b, double c)
{
    double L = (u*u + v * v + w * w);
	double L_root=sqrt(L);
    angle = angle * M_PI / 180.0; //converting to radian value
	double cs(cos(angle)),sn(sin(angle));
    double u2 = u * u;
    double v2 = v * v;
    double w2 = w * w; 
 
    rotationMatrix[0][0] = (u2 + (v2 + w2) * cs) / L;
    rotationMatrix[0][1] = (u * v * (1 - cs) - w * L_root * sn) / L;
    rotationMatrix[0][2] = (u * w * (1 - cs) + v * L_root * sn) / L;
    rotationMatrix[0][3] = (((a*(v2+w2)-u*(b*v+c*w))*(1-cs)) + ((b*w-c*v)*L_root*sn))/L;
 
    rotationMatrix[1][0] = (u * v * (1 - cs) + w * L_root * sn) / L;
    rotationMatrix[1][1] = (v2 + (u2 + w2) * cs) / L;
    rotationMatrix[1][2] = (v * w * (1 - cs) - u * L_root * sn) / L;
    rotationMatrix[1][3] = (((b*(u2+w2)-v*(a*u+c*w))*(1-cs)) + ((c*u-a*w)*L_root*sn))/L;
 
    rotationMatrix[2][0] = (u * w * (1 - cs) - v * L_root * sn) / L;
    rotationMatrix[2][1] = (v * w * (1 - cs) + u * L_root * sn) / L;
    rotationMatrix[2][2] = (w2 + (u2 + v2) * cs) / L;
    rotationMatrix[2][3] = (((c*(u2+v2)-w*(a*u+b*v))*(1-cs)) + ((a*v-b*u)*L_root*sn))/L;; 
 
    rotationMatrix[3][0] = 0.0;
    rotationMatrix[3][1] = 0.0;
    rotationMatrix[3][2] = 0.0;
    rotationMatrix[3][3] = 1.0;
}  
void PerformRotation(double point_x,double point_y,double point_z,//point to rotate
   	                        double u_axis, double v_axis, double w_axis,//axis to rotate about
					        double a, double b, double c, //point the axis of rotation passes through
	                        double angle,//angle of rotation
					        double&rot_point_x,double&rot_point_y,double&rot_point_z)//point after roation
{
                
    inputMatrix[0][0] = point_x;
    inputMatrix[1][0] = point_y;
    inputMatrix[2][0] = point_z;
    inputMatrix[3][0] = 1.0;  

    setUpRotationMatrix(angle, u_axis, v_axis, w_axis,a,b,c);
    multiplyMatrix();
	rot_point_x=outputMatrix[0][0];
	rot_point_y=outputMatrix[0][1];
	rot_point_z=outputMatrix[0][2];

	// showPoint(); 

}