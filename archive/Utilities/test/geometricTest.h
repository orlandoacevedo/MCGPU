#ifndef GEOMETRICTEST_H
#define GEOMETRICTEST_H

#include "../geometricUtil.h"
#include "../metroUtil.h"
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <time.h>

#define PRECISION .0001

/*
    tests the getNormal function
*/
void testGetNormal();

/*
    tests the getAngle(Plane, Plane) function
*/
void testGetAngleBetweenPlanes();

/*
    tests the getBond function
*/
void testGetBond();

/*
    tests the getAllBonds function
*/
void testGetAllBonds();

/*
    tests the getIntersection function
*/
void testGetIntersection();

/*
    tests the isMember function
*/
void testIsMember();


/**
 test degrees to radians and radians to degrees
*/
void testD2RandR2D();

/**
 test getOpositeAtom
*/
void testGetOppositeAtom();

/**
 test getOpositeAtom
*/
void testGetCommonAtom();

/**
 test getDistance
*/
void testGetDistance();

/**
 test getAngle
*/
void testGetAngle();

/**
    tests the translateAtom function
*/
void testTranslateAtom();

/**
    tests the rotateAboutX function
*/
void testRotateAboutX();

/**
    tests the rotateAboutY function
*/
void testRotateAboutY();

/**
    tests the rotateAboutZ function
*/
void testRotateAboutZ();

/**
  tests rotating an atom about a vector.
*/
void testRotateAboutVector();

/**
  tests rotating within a plane
*/
void testRotateInPlane();

/**
  runs all the test in this file
*/

void testGeometric();

#endif //GEOMETRICTEST_H
