#include "DrawSpheres.h"
#include <iostream>
#include <string.h>  
#include <fstream>
#include <sstream>
#include <math.h>

inline void PointOnSphere(double&x1, double&y1, double&z1,double rad);






int DrawSpheres::Draw(std::string inputfilename, std::string outputfilename,size_t refine_level)
{
	size_t _num_samples;// number of samples to draw
	double **_samples; // samples coordinates
	size_t _sphere_points; // number of points of the sphere  
	double **_sphere; // sphere points coordinates 
	size_t _num_faces;// number of face (triangles) making up the sphere mesh
	size_t **_sphere_mesh; // sphere mesh connectivity
	

	size_t V,d;

	//sample_file.open("samples.txt");
		
	//char inputfilename[50];
	std::ifstream sample_file;
	sample_file.open(inputfilename);

	/*while(!sample_file.is_open()){
		cout<<"File doesn't exist!!!"<<endl;
		cout<<"Enter the file name again"<<endl;
		cin.getline(inputfilename,50);
		sample_file.open(inputfilename);
	}*/

	sample_file>>_num_samples;

	_samples=new double*[_num_samples];
	for(V=0;V<_num_samples;V++){
		_samples[V]=new double[4];
		sample_file>>_samples[V][0]>>_samples[V][1]>>_samples[V][2]>>_samples[V][3];
	}

	sample_file.close();


	// read base sphere	
	std::ifstream sphere_file;	


	if(refine_level==1){
		sphere_file.open("sphere1.obj");
	}
	else if(refine_level==2){
		sphere_file.open("sphere2.obj");
	}
	else if(refine_level==3){
		sphere_file.open("sphere3.obj");
	}
	else if(refine_level==4){
		sphere_file.open("sphere4.obj");
	}
	else if(refine_level==5){
		sphere_file.open("sphere5.obj");
	}
	else if(refine_level==6){
		sphere_file.open("sphere6.obj");
	}
	else {
		std::cout<<"Error at  DrawSpheres()"<<std::endl;
		std::cout << "Level of refined sphere was wrong " << std::endl;
		std::cout << "Please re-enter a value between 1 to 6" << std::endl;
		std::cout << "I'll terminate now!! " << std::endl;
		return 0;

	}

	

	char tuna[3];
	
	sphere_file>>tuna>>_sphere_points;
	sphere_file>>tuna>>_num_faces;

	_sphere=new double*[_sphere_points];

	for(V=0;V<_sphere_points;V++){ // read coordinates 
		_sphere[V]=new double[3];
		sphere_file>>tuna>>_sphere[V][0]>>_sphere[V][1]>>_sphere[V][2];
	}

	_sphere_mesh=new size_t*[_num_faces];
	for(V=0;V<_num_faces;V++){ // read connectivity 
		_sphere_mesh[V]=new size_t[3];
		sphere_file>>tuna>>_sphere_mesh[V][0]>>_sphere_mesh[V][1]>>_sphere_mesh[V][2];
	}


	_x_sphere=0.5; // sphere center and radius
	_y_sphere=0.5;
	_z_sphere=0.5;
	_r_sphere=0.5*sqrt(3.0);

	for(V=0;V<_sphere_points;V++){
		_sphere[V][0]-=0.5;
		_sphere[V][1]-=0.5;
		_sphere[V][2]-=0.5;
	}

	_x_sphere=0.0; // new center
	_y_sphere=0.0;
	_z_sphere=0.0;


	for(V=0;V<_sphere_points;V++){
		PointOnSphere(_sphere[V][0],_sphere[V][1],_sphere[V][2],1.0);
	}

	_r_sphere=1.0; // new center

	std::ofstream Out(outputfilename,std::ios::out);

	//fstream Out("Output.obj", ios::out);
	Out<<"#v "<<_num_samples*_sphere_points<<std::endl;
	Out << "#f " << _num_samples*_num_faces << std::endl;;

	for(V=0;V<_num_samples;V++){

		for(d=0;d<_sphere_points;d++){
			PointOnSphere(_sphere[d][0],_sphere[d][1],_sphere[d][2],_samples[V][3]);
			Out << "v " << _sphere[d][0] + _samples[V][0] << " " << _sphere[d][1] + _samples[V][1] << " " << _sphere[d][2] + _samples[V][2] << std::endl;;

		}
	}

	size_t huh=0;
	for(V=0;V<_num_samples;V++){

		for(d=0;d<_num_faces;d++){
			Out << "f " << _sphere_mesh[d][0] + huh << " " << _sphere_mesh[d][1] + huh << " " << _sphere_mesh[d][2] + huh << std::endl;;

		}

		huh+=_sphere_points;
	}


	delete[]_samples;
	delete[]_sphere_mesh;
	delete[]_sphere;

	return 0;

}

void DrawSpheres::PointOnSphere(double&x1, double&y1, double&z1, double rad)
{	
	// changing the point location so that it lies on the new radius 

	double x_v(x1-_x_sphere), y_v(y1-_y_sphere), z_v(z1-_z_sphere);
	
	double n=sqrt(x_v*x_v + y_v*y_v + z_v*z_v);
		
	x1=x_v*rad/n + _x_sphere;
	y1=y_v*rad/n + _y_sphere;
	z1=z_v*rad/n + _z_sphere;


}