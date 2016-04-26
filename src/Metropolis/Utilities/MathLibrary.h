/*!\file
  \brief Structures and functions used to do geometric transformations.
  \author Alexander Luchs, Riley Spahn, Seth Wooten
 
 */


#ifndef GEOMETRICUTIL_H
#define GEOMETRICUTIL_H

//using namespace std;

//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//  INCLUDE statements
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________

#include <math.h> //AlbertExcludes
//#include "math.h" //for Albert's testing //AlbertIncludes

#include <assert.h>
#include <errno.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>


#include "StructLibrary.h" //AlbertExcludes
//#include "StructLibrary.cpp" //AlbertIncludes
//include Metropolis/DataTypes.h

//_________________________________________________________________________________________________________________
//  Specific namespace/using requirements
//_________________________________________________________________________________________________________________
		//to avoid "using namespace std;"
using std::vector;
using std::string;
using std::map;  


#define DOUBPREC .000001
#define PI 3.14159265

void seed(int seed);
Real randomReal(const Real start, const Real end);
void adjustAtomIDs(Molecule* molec, int m);
void revertAtomIDs(Molecule* molec, int m);

/**
  Structure representing a geometic point.
*/
struct Point{
    /**
      x coordinate of a point.
    */
    double x;
    /**
      y coordinate of a point.
    */
    double y;
    /**
      z coordinate of a point.
    */
    double z;
};

/**
  @param x - the x coord of the point.
  @param y - the y coord of the point.
  @param z - the z coord of the point.
  @return - a point at x,y,z
*/
Point createPoint(double X, double Y, double Z);

/**
  Structure representing a vector with an origin and end point.
*/
struct Vector{
    /**
      Point representing the start point fo the vector.
    */
    Point origin;
    /**
      Point representing the ending point for the vector.
    */
    Point end;
};

/**
  @param startPoint - the starting point of the vector.
  @param endPoint - the ending point of the vector.
  @return - a vector that starts at startPoint and ends at endPoint
*/
Vector createVector(Point startPoint, Point endPoint);

/**
  @param startPoint - the starting point for the vector.
  @param deltaX - the change in the x direction from the startPoint to the
  end point of the vector.
  @param deltaY - the change in the y direction from the startPoint to the
  end point of the vector.
  @param deltaZ - the change in the z direction from the startPoint to the
  end point of the vectoro.
  @return - a vector that starts at startPoint and ends at
  (startPoint.x + deltaX, startPoint.y + deltaY, startPoint.z + deltaZ)
*/
Vector createVector(Point startPoint, double deltaX, double deltaY, double deltaZ);

/**
  Structure representing a geometric plane.
  A plane is described by three points.
*/
struct Plane{
    /**
      The first atom (point) to define the plane.
    */
    Atom atom1;
    /**
      The second atom (point) to define the plane.
    */
    Atom atom2;
    /**
      The third atom (point) to define the plane.
    */
    Atom atom3;
};

/**
  Prints position information about point p.
  @param p - the point to be printed.
*/
void printPoint(Point p);

/**
  @param p1 - the first point in the plane.
  @param p2 - the second point in the plane.
  @param p3 - the third point in the plane.
  @return - a plane represented by points p1, p2 and p3
*/
Plane createPlane(Atom p1, Atom p2, Atom p3);

/**
  Returns a specific atom from a vector of atoms.
  @param atoms - vector of atoms from which to get the Atom
  @param atomID - the ID of the atom to find.
  @return - the atom with atomID or Atom with parameters = -1 if not found.
*/
Atom getAtom(vector<Atom> atoms, unsigned long atomID);

/**
  Returns a bond (if it exists) containing a1 and a2.
  @param bonds - the vector of bonds to be tested.
  @param a1 - the first atom needed in the bond
  @param a2 - the second atom needed in the bond
  @return - the bond containing a1 and a2
          - Bond(-1,-1,false) if a suitable bond has not been found
 */
Bond getBond(const vector<Bond>& bonds, unsigned long a1, unsigned long a2);

/**
  Returns a vector of all of the atoms bonded to atom id
  @param bonds - the vector of bonds to be checked.
  @param atomID - the atom that must be involved in the bonds
  @return - a vector of atomIDs that are bonded to  atomID
*/
vector<unsigned long> getAllBonds(vector<Bond> bonds, unsigned long atomID);

/**
  Returns the intersection of the two vectors.
  @param v1 - the first vector  set.
  @param v2 - the seconds vector set.
  @return - vector containing the elements of the intersection of the two sets.
*/
vector<unsigned long> getIntersection(vector<unsigned long> v1, vector<unsigned long> v2);

/**
  Determines wheter the atom is a member of the vector.
  @param atoms - the list of atom ids to check.
  @param toCheck - the atomID to be checked for in the vector.
  @return - true if toCheck is a member of the vector.
            false if toCheck is not a member of the vector.
*/
bool isMember(vector<unsigned long> atoms, unsigned long toCheck);

/**
  @param degrees - the measure of an angle in degrees
  @return - the value of the same angle in radians.
*/
double degreesToRadians(double degrees);

/**
  @param radians - the measure of an angle in radians
  @return - the value of the same angle in degrees
*/
double radiansToDegrees(double radians);

/**
  Returns the id in the bond that is not atomID
  @param bond - the bond to be checked.
  @param atomID - the atomID that you want the opposite of.
  @return - the atom id of the atom opposite atomID
            -1 if atomID is not found in the bond
*/
unsigned long getOppositeAtom(Bond bond, unsigned long atomID);

/**
  Returns the id in the angle that is not atomID
  @param angle - the angle to be checked.
  @param atomID - the atomID that you want the opposite of.
  @return - the atom id of the atom opposite atomID
            -1 if atomID is not found in the angle
*/
unsigned long getOppositeAtom(Angle angle, unsigned long atomID);

/**
  Returns the id of the atom in the dihedral that is not atomID
  @param dihedral - the dihedral to be checked.
  @param atomID - the atomID that do not want to return.
  @return - the atom id of the atom opposite atomID
            -1 if atomID is not found in the dihedral 
*/
unsigned long getOppositeAtom(Dihedral dihedral, unsigned long atomID);

/**
  Returns the id of the third atom in the angle between atom1 and atom2.
  Returns -1 if no common angle is found.
  @param bonds - vector of bonds to be checked.
  @param atom1 - the first atom in the angle.
  @param atom2 - the second atom in the angle.
  @return - the id of the third atom in the angle or -1 if none exists.
*/
unsigned long getCommonAtom(vector<Bond> bonds, unsigned long atom1,
        unsigned long atom2);

/**
  @param atom1 - the first atom.
  @param atom2 - the second atom.
  @return - the distance between atom1 and atom2.
*/
double getDistance(Atom atom1, Atom atom2);

/**
  Determines if an atom is on the xz plane
  @atom - the atom to examine.
  @return - true if the atom is in the xz plane
*/
bool inXZPlane(Atom atom);

/**
  Returns the angle made by the three atoms.
  @param atom1 - edge of the angle. Bonded to atom2
  @param atom2 - the corner of the angle, Bonded to atom1 and atom3.
  @param atom3 - the other edge of the angle.  Bonded to atom2.
  @return - the angle in degrees created by the three atoms.
*/
double getAngle(Atom atom1, Atom atom2, Atom atom3);

/**
  Returns the acute angle between the two planes.
  @param p1 - the first plane that creates the angle.
  @param p2 - the second plane that creates the angle.
  @return - the acute angle in degrees between the planes.
          - returns 90 degrees if the planes are perdenicular.
*/
double getAngle(Plane p1, Plane p2);

/**
  Returns the normal vector of a plane.
  The vector has a start point of 0,0,0 and an end point of @return.
  @param p - the plane used to find the normal.
  @return - the end point of the normal vector
*/
Point getNormal(Plane p);

/**
  Returns a point that creates a normal vector.
  @param atom1 - the first atom in line (start).
  @param atom2 - the middle atom in the line.
  @param atom3 - the end of the line.
  @return - point that can be used to create an normal vector with the line)
*/
Point getNormal(Atom atom1, Atom atom2, Atom atom3);

/**
  Returns the distance between an atom and a line
  @param atom0 - base atom
  @param atom1 - first atom in the defined line
  @param atom2 - second atom in the defined line
  @param atom3 - third atom in the defined line
  @return - the distance between atom0 and the line degined by atoms 1, 2, and 3
*/
double getDistance(Atom atom0, Atom atom1, Atom atom2, Atom atom3);

/**
  Translates the atom by x, y, z
  @param atom - the atom to be translated.
  @param x - the distance in the x direction to be translated.
  @param y - the distance in the y direction to be translated.
  @param z - the distance in the z direction to be translated.
  @return - a copy of atom that has been translated accordingly.
*/
Atom translateAtom(Atom atom, double x, double y, double z);

/**
  @param atom - the atom to be rotated
  @param theta - the distance in degrees to be rotated.
  @return - atom rotated theta degrees about the x axis.
*/
Atom rotateAboutX(Atom atom, double theta);

/**
  @param atom - the atom to be rotated
  @param theta - the distance in degrees to be rotated.
  @return - atom rotated theta degrees about the y axis.
*/
Atom rotateAboutY(Atom atom, double theta);

/**
  @param atom - the atom to be rotated
  @param theta - the distance in degrees to be rotated.
  @return - atom rotated theta degrees about the z axis.
*/
Atom rotateAboutZ(Atom atom, double theta);


/**
  Rotates atom1 about atom2 in the plane defined by atom1, atom2 and atom3
  theta degrees.
  @param atom1 - the atom to be rotated.
  @param atom2 - the atom about which atom1 will be rotated.
  @param atom3 - the atom used to complete the plane.
  @param theta - the number of degrees to be rotated.
  @return - atom that has been rotated theta degrees in the plane defined
  by the locations of atom1, atom2 and atom3.
  If these points are collinear then the atom will be rotated about the 
  vector created by atom2 and (atom2.x + 1, atom2.y, atom2.z)
*/
Atom rotateAtomInPlane(Atom atom1, Atom atom2, Atom atom3, double theta);

/**
  Rotates atom1 about the vector defined by atom3 and atom2.
  @param atom1 - the atom to be rotated.
  @param atom2 - the atom defining the start point of the vector.
  @param atom3 - the atom defining the end point of the vector.
  @param theta - number of degrees to rotate
  @return - atom1 post-rotation
*/
Atom rotateAtomAboutVector(Atom atom1, Atom atom2, Atom atom3, double theta);

/**
  @param a - the first double to compare.
  @param b - the second double to compare.
  @param precision - the threshold for equality
  @return if the difference in a & b is less than precision then returns true
         - else returns false.
*/
bool compareDoubleDifference(double a, double b, double precision);

/**
  Translates a molecule based on the xTrans, yTrans and zTrans parameters
  and then rotates the atoms in the molecule about an atom.
  @param molec - the molecule to be translated and rotated
  @param pivot - the atom that will be the point about which all of the otehr
  atoms are rotated.
  @param xTrans - the distance in angstroms that the molecule will be translated
  in the x direction.
  @param yTrans - the distance in angstroms that the molecule will be translated
  in the y direction.
  @param zTrans - the distance in angstroms that the molecule will be translated
  in the z direction.
  @param xRot - the amount to be rotated about the vector through the pivot atom 
  parallel to the x axis. (degrees)
  @param yRot - the amount to be rotated about the vector through the pivot atom 
  parallel to the y axis. (degrees)
  @param zRot - the amount to be rotated about the vector through the pivot atom 
  parallel to the z axis. (degrees)
  @return - the molecule after the transformations
*/
Molecule moveMolecule(Molecule molec, Atom pivot, double xTrans, double yTrans,
        double zTrans, double xRot, double yRot, double zRot);

/**
	Returns a random number
	@param start - lowest value possible
	@param end - largest value possible
*/
double randomNUM(const double start, const double end);

/**
  Translates a molecule a random amount and returns it
  @param molecule - the molecule to be translated
  @param maxTranslation - the maximum distance (Ang) for each axis of translation
*/
Molecule translateMolecule(Molecule molec, double maxTranslation);

/**
  Rotates a molecule about a given atom a random amount and returns it
  @param molecule - the molecule to be rotated
  @param pivot - the atom that the molecule is rotated about
  @param maxRotation - the maximum number of degrees for each axis of rotation
*/
Molecule rotateMolec(Molecule molec, Atom pivot, double maxRotation);

/**
  Checks the vector of bonds to see if the atom has a bond with another atom in the Molecule.
  If the atomId is 1, and X is other atoms it looks for 1--X and not X--1.
  @param Bonds - the vector of Bond structs that it seaches.
  @param atomId - the atoms id to check for.
  @return -  the index in the Bond vector that matches the Bond structure. 
  Returns -1 if not found. 
*/
int hasBond(vector<Bond> Bonds, unsigned long atomId);

/**
  Checks the vector of angles to see if the atom has a angle with another atom in the Molecule.
  If the atomId is 1, and X is other atoms it looks for 1--X and not X--1.
  @param Angles - the vector of Angle structs that it seaches.
  @param atomId - the atoms id to check for.
  @return -  the index in the Angle vector that matches the Angle structure. 
  Returns -1 if not found. 
*/
int hasAngle(vector<Angle> Angles, unsigned long atomId);

/**
  Checks the vector of dihedrals to see if the atom has a dihedral with another atom in the Molecule.
  If the atomId is 1, and X is other atoms it looks for 1--X and not X--1.
  @param Dihedrals - the vector of Dihedral structs that it seaches.
  @param atomId - the atoms id to check for.
  @return -  the index in the Dihedrals vector that matches the Dihedral structure. 
  Returns -1 if not found. 
*/
int hasDihedral(vector<Dihedral> Dihedrals, unsigned long atomId);

/**
  Sets the Bond, Angle, and Dihedral vectors with their corresponding values from the Molecule
  where the the structures have an atom1 id <= lineAtomId. Used in buildMoleculeXYZ() to
  fill the vectors with defined Structs from the Molecule that contain relationships between
  atoms that are >= to the lineAtom.
  @param *molec - the molecule structure to copy the data from.
  @param &lineAtomId -  the atom id used to find defined structs
  @param &bondVector - the vector of bonds to fill
  @param &angleVector - the vector of angles to fill
  @param &dihedralVector -  the vector of dihedrals to fill 
*/
void setMoleculeVectors(Molecule *molec, int numBonded, unsigned long lineAtomId, vector<Bond> &bondVector, 
    vector<Angle> &angleVector, vector<Dihedral> &dihedralVector);

/**
  Calculates and sets the Atom positions in a Molecule based on the atom to atom relationships
  defined in the Bond, Angle, and Dihedral structures. Builds the Molecule structure one atom at
  a time, starting with the first atom in the Atom array. Builds the Molecule around the first
  atoms position.  
  @param *molec - an array of bonded molecules to be set
  @param numBounded - the number of bonded molecules being passed in
*/
void buildMoleculeInSpace(Molecule *molec,int numBonded=1);

/**
  Converts z-matrix to Cartesian coordinates.
  Calculates and sets the Atom positions in a Molecule based on the atom to atom relationships
  defined in the Bond, Angle, and Dihedral structures. 
  @param *molec - an array of bonded molecules to be set
  @param numBounded - the number of bonded molecules being passed in
*/
void buildMoleculeXYZ(Molecule *molec,int numBonded=1);

/**
  Calculates the cross product for the transformation (direction cosine) matrix
  in buildMoleculeXYZ(). Returns normalized values.
  @param *A - accepts an array of vectors for bond, angle, dihedrals coordinates
  @param *B - accepts an array of vectors for bond, angle, dihedrals coordinates
  @param *C - returns array of normalized coordinates
*/
void cross(double *A, double *B, double *C);

#endif //GEOMETRICUTIL_H
