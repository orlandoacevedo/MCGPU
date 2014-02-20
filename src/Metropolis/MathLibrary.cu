#include "MathLibrary.cuh"

using namespace std;

//Point
Point createPoint(double x, double y, double z)
{
	Point point;
	point.x = x;
	point.y = y;
	point.z = z;
	return point;
}
void printPoint(Point point);


//Vector
Vector createVector(Point origin, double deltaX, double deltaY, double deltaZ)
{
	Point end;
	end = origin;
	end.x += deltaX;
	end.y += deltaY;
	end.z += deltaZ;
	
	Vector vector;
	vector.origin = origin;
	vector.end = end;

	return vector;
}


//Plane
Plane createPlane(Point point, Vector normal)
{
	Plane plane;
	plane.point = point;
	plane.normal = normal;
	return plane;
}