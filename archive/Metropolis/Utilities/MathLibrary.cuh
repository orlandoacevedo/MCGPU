/*
	Contains all mathematical resources for the simulation

	Author: Nathan Coleman
	Last Changed: February 19, 2014
*/

#ifndef MATHLIBRARY_H
#define MATHLIBRARY_H

#include <math.h>

//Forward declaration so that each can be used in methods below
struct Point
{
	double x, y, z;
};
struct Vector
{
	Point origin, end;
};
struct Plane
{
	Point point;
	Vector normal;
};

//Point
Point createPoint(double x, double y, double z);
void printPoint(Point point);


//Vector
Vector createVector(Point origin, double deltaX, double deltaY, double deltaZ);


//Plane
Plane createPlane(Point point, Vector normal);
Vector getNormal(Plane plane){return plane.normal;};

#endif