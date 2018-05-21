#ifndef _DRAWDEBUG_
#define _DRAWDEBUG_

#include "DrawSpheres.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

void DrawOneSphere(int ip, double xx, double yy, double zz, double rr, size_t res)
{
	//std::cout << "\n I AM DRAWING in DrawOneSphere()" << std::endl;
	std::fstream file("ip.dat", std::ios::out);
	file << 1 << std::endl;
	file << xx << " " << yy << " " << zz << " " << rr << std::endl;

	DrawSpheres dr;
	dr.Draw("ip.dat", "debug_out/one_" + std::to_string(ip) + ".obj", res);	
	remove("ip.dat");
}
void DrawManySpheres(std::string filename, size_t num, double**spheres, bool have_radius)
{
	//std::cout << "\n I AM DRAWING in DrawManySpheres()" << std::endl;
	std::fstream file("ip.dat", std::ios::out);
	file << num << std::endl;
	size_t V;
	for (V = 0; V<num; V++){
		if (have_radius){
			file << spheres[V][0] << " " << spheres[V][1] << " " << spheres[V][2] << " " << spheres[V][3] << std::endl;
		}
		else{
			file << spheres[V][0] << " " << spheres[V][1] << " " << spheres[V][2] << " " << 0.005 << std::endl;
		}
	}
	DrawSpheres dr;
	dr.Draw("ip.dat", filename, 3);
	remove("ip.dat");
}

void DrawManySpheres(std::string filename, size_t num, double**spheres, double fixed_radius)
{
	//std::cout << "\n I AM DRAWING in DrawManySpheres()" << std::endl;
	std::fstream file("ip.dat", std::ios::out);
	file << num << std::endl;
	size_t V;
	for (V = 0; V<num; V++){		
		file << spheres[V][0] << " " << spheres[V][1] << " " << spheres[V][2] << " " << fixed_radius << std::endl;
	}
	DrawSpheres dr;
	dr.Draw("ip.dat", filename, 3);
	remove("ip.dat");
}

void DrawManySpheres(std::string filename, size_t num, std::vector<std::vector<double>>spheres, double fixed_radius)
{
	//std::cout << "\n I AM DRAWING in DrawManySpheres()" << std::endl;
	std::fstream file("ip.dat", std::ios::out);
	file << num << std::endl;
	size_t V;
	for (V = 0; V<num; V++){
		file << spheres[V][0] << " " << spheres[V][1] << " " << spheres[V][2] << " " << fixed_radius << std::endl;
	}
	DrawSpheres dr;
	dr.Draw("ip.dat", filename, 4);
	remove("ip.dat");
}

#ifdef VERT_DEFINED
void DrawManySpheres(std::string filename, std::vector<vert>spheres, double fixed_radius, size_t res)
{
	std::cout << "\n I AM DRAWING in DrawManySpheres()" << std::endl;
	std::fstream file("ip.dat", std::ios::out);
	file << spheres.size() << std::endl;
	size_t V;
	for (V = 0; V<spheres.size(); V++){
		if (fixed_radius < 0){
			file << spheres[V].x[0] << " " << spheres[V].x[1] << " " << spheres[V].x[2] << " " << sqrt(spheres[V].r) << std::endl;
		}
		else{
			file << spheres[V].x[0] << " " << spheres[V].x[1] << " " << spheres[V].x[2] << " " << fixed_radius << std::endl;
		}
	}
	DrawSpheres dr;
	dr.Draw("ip.dat", filename, res);
	remove("ip.dat");
}
#endif

void DrawManySpheresSpecial(int ip, size_t num, double**spheres, double flag, double frac, size_t res)
{
	//std::cout << "\n I AM DRAWING in DrawManySpheresSpecial()" << std::endl;
	std::fstream file("ip.dat", std::ios::out);
	file.precision(30);

	size_t real_num(0);
	for (size_t V = 0; V<num; V++){
		if (abs(spheres[V][4] - flag)>0.0000001){ continue; }
		real_num++;
	}
	file << real_num << std::endl;

	for (size_t V = 0; V<num; V++){
		if (abs(spheres[V][4] - flag)>0.0000001){ continue; }
		file << spheres[V][0] << " " << spheres[V][1] << " " << spheres[V][2] << " " << spheres[V][3] / frac << std::endl;
	}

	DrawSpheres dr;	
	dr.Draw("ip.dat", "debug_out/many_" + std::to_string(ip) + ".obj", res);
	remove("ip.dat");
}

void DrawManySpheresSpecial(int ip, size_t num, std::vector<std::vector<double>>spheres, double flag, double frac, size_t res)
{
	//std::cout << "\n I AM DRAWING in DrawManySpheresSpecial()" << std::endl;
	std::fstream file("ip.dat", std::ios::out);
	file.precision(30);

	size_t real_num(0);
	for (size_t V = 0; V<num; V++){
		if (abs(spheres[V][4] - flag)>0.0000001){ continue; }
		real_num++;
	}
	file << real_num << std::endl;

	for (size_t V = 0; V<num; V++){
		if (abs(spheres[V][4] - flag)>0.0000001){ continue; }
		file << spheres[V][0] << " " << spheres[V][1] << " " << spheres[V][2] << " " << spheres[V][3] / frac << std::endl;
	}

	DrawSpheres dr;
	dr.Draw("ip.dat", "debug_out/many_" + std::to_string(ip) + ".obj", res);
	remove("ip.dat");
}

void DrawPlane(double n_x, double n_y, double n_z, double pointx, double pointy, double pointz)
{
	//std::cout << "\n I AM DRAWING in DrawPlane()" << std::endl;

	double dd = -1.0*n_x*pointx - n_y*pointy - n_z*pointz;
	double scale(1.9);
	double a, b, c;
	a = scale; b = scale;
	c = (-1.0*dd - a*n_x - b*n_y) / n_z;
	std::string ff = "debug_out/plane.obj";

	std::fstream plot(ff, std::ios::out);

	plot << "v " << a << " " << b << " " << c << std::endl;

	a = scale; b = -scale;
	c = (-1.0*dd - a*n_x - b*n_y) / n_z;
	plot << "v " << a << " " << b << " " << c << std::endl;

	a = -scale; b = scale;
	c = (-1.0*dd - a*n_x - b*n_y) / n_z;
	plot << "v " << a << " " << b << " " << c << std::endl;

	a = -scale; b = -scale;
	c = (-1.0*dd - a*n_x - b*n_y) / n_z;
	plot << "v " << a << " " << b << " " << c << std::endl;

	plot << "f 1 3 4 2" << std::endl;

}
void DrawOneTri(int ip, 
	                   double xip1, double yip1, double zip1, 
					   double xip2, double yip2, double zip2, 
					   double xip3, double yip3, double zip3)
{
	//std::cout << "\n I AM DRAWING in DrawOneTri()" << std::endl;

	std::stringstream fname;
	fname << "debug_out/Tri_" << ip << ".obj";
	std::fstream file(fname.str().c_str(), std::ios::out);
	file << "v " << xip1 << " " << yip1 << " " << zip1 << std::endl;
	file << "v " << xip2 << " " << yip2 << " " << zip2 << std::endl;
	file << "v " << xip3 << " " << yip3 << " " << zip3 << std::endl;
	file << "f 1 2 3" << std::endl;
}

#endif 