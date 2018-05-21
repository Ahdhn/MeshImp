#ifndef __DRAWSPHERES_
#define __DRAWSPHERES_
#include <string> 
#include <string.h>  
using std::string;

class DrawSpheres
{
public:

	int Draw(string, string, size_t);
private :
	void PointOnSphere(double&x1, double&y1, double&z1, double rad);
	double _x_sphere, _y_sphere, _z_sphere, _r_sphere; // sphere center and radius 
};
#endif