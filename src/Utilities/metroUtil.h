/*!\file
  \brief Useful abstractions of for use in linear and parallel simulations.
  \author Alexander Luchs, Riley Spahn, Seth Wooten
 
 */

#ifndef METROUTIL_H
#define METROUTIL_H


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <time.h>

using namespace std;

//potentially to be used for an atom, to avoid passing x y and z separately. Supports doubles, for current implementation.
struct Vector3
{
    /*!x coordinate of the atom*/
    double x;
    /*!y coordinate of the atom*/
    double y;
    /*!z coordinate of the atom*/
    double z;
};

//potentially to be used for an atom, to avoid passing x y and z separately. Supports floats, for current implementation.
struct Vector3f
{
    /*!x coordinate of the atom*/
    float x;
    /*!y coordinate of the atom*/
    float y;
    /*!z coordinate of the atom*/
    float z;
};

/*!
Structure that represents an atom in the simulation.
*/
struct Atom
{
	//the following coordinates can be modified to use the newly created structs above
    /*!x coordinate of the atom*/
    double x;
    /*!y coordinate of the atom*/
    double y;
    /*!z coordinate of the atom*/
    double z;

    /*!unique id of the atom*/
    unsigned long id; 

    /*!sigma value for the atom for the LJ calculations.*/
    double sigma;
    /*!  epsilon value for the atom for the LJ calculation*/
    double epsilon;
    /*!epsilon value for the atom for the LJ calculation*/
    double charge;
    /*!name for atom in z-matrix.*/
    char name;

};



/*!
  @param id - the id of the atom to be created.
  @param x - the x coordinate of the atom.
  @param y - the y coordinate of the atom.
  @param z - the z coordinate of the atom.
  @param sigma - the sigma value of the atom to be used in the LJ calculation.
  @param epsilon - the epsilon value of the atom to be used in the LJ calculation.
  @param charge - the charge of the atom used in the energy calculation.
  @param name - the name of the atomtype in the OPLS file.
  @return - an instance of the an atom structure.
*/
Atom createAtom(unsigned long id, double x, double y, double z, double sigma, double epsilon, double charge, char name);

/*!
  @param id - the id of the atom to be created.
  @param x - the x coordinate of the atom.
  @param y - the y coordinate of the atom.
  @param z - the z coordinate of the atom.
  @param sigma - the sigma value of the atom to be used in the LJ calculation.
  @param epsilon - the epsilon value of the atom to be used in the LJ calculation.
  @return - an instance of the an atom structure.
*/
Atom createAtom(unsigned long id, double x, double y, double z, double sigma, double epsilon);

/*!
  @param id - the id of the atom to be created.
  @param x - the x coordinate of the atom.
  @param y - the y coordinate of the atom.
  @param z - the z coordinate of the atom.
  @return - an instance of the an atom structure.
*/
Atom createAtom(unsigned long id, double x, double y, double z);

/*!
  Structure representing the simulation's environment.
*/
struct Environment
{
    /*!
    length of the box in the x direction
    */
    double x;
    
    /*!
    length of the box in the y direction
    */
    double y; 
    
    /*!
    length of the box in the z direction
    */
    double z; 

    /*!
    the maximum distance in angstroms that an atom can move
    */
    double maxTranslation; 
    
    /*!
    the temperature of the box in kelvin
    */
    double temperature;
    
    /*!
    the number of atoms in the environment
    */
    int numOfAtoms; 
    
    /*!
    the number of molecues in the environment
    */
    int numOfMolecules;
    
    /*!
    the nonbonded cutoff distance in the environment
    */
    double cutoff; 
    
    /*!
    the maximum rotation in degrees that an atom can rotate
    */
    double maxRotation;

    /*!
    seed used in rand() function, make it repeatable 
    */
    unsigned int randomseed;
	
    /*!
    Which atom(s) from molecule to use as primaries when cutting off atoms. 
    */
	int primaryAtomIndex;
    
};

/*!
   @param x - the x dimension of the environment (box).
   @param y - the y dimension of the environment (box).
   @param z - the z dimension of the environment (box).
   @param maxTrans - the maximum translation of an atom/molecule in the simulation.
   @param temp - the temperature in kelvin of the environment.
   @param numOfAtoms - the number of atoms in the environment.
   @param cutoff - the nonbonded cutoff distance in the environment.
   @param maxRot - the maximum rotation of an atom/molecule in the simulation.
   @return - an instance of an environment structure.
*/
Environment createEnvironment(double x, double y, double z, double maxTrans, double temp, int numOfAtoms, double cutoff, double maxRot);

/**
  Writes the list of atoms to the file named filename
  @param atoms - list of atoms to be written
  @param enviro - environmental information
  @param filename - name of the file to be written
  @pararm accepts - the number of moves accepted by the system
  @param rejects - the number of moves rejected by the system
  @param totalEnergy - the totalEnergy in the system
*/
void writeOutAtoms(Atom *atoms, Environment *enviro, string filename, int accepts, int rejects, double totalEnergy);

/*!
  Structure to respresent bonds between atoms in a molecule.
*/
struct Bond
{
    /*!
       The first atom in the bond.
    */
    int atom1;
    /*!
       The second atom in the bond.
    */
    int atom2;

    /*!
       The length of the bond.
    */
    double distance;
    /*!
       Bool indicating if the bond can be varied.
    */
    bool variable;
};

/**!
  @param atom1 - the id of the first atom in the bond
  @param atom2 - the id of the second atom in the bond
  @param distance - the distance between the the two atoms
  @param variable - boolean if distance is variable
  @return - instance of the new bond
*/
Bond createBond(int atom1, int atom2, double distance, bool variable);

/*!
  Structure representing an angle between two atoms in a molecule.
  A geometric angle is defined using three points(atoms).  Only the end
  points are listed here.  It is up to the programmer to find the third
  atom in the angle.
*/
struct Angle
{
    /*!
    the first atom in the angle
    */
    int atom1; 
    
    /*!
    the second atom in the angle
    */
    int atom2;
    
    /*!
    the angle between the atoms
    */
    double value; 
    
    /*!
    if the angle is variable
    */
    bool variable; 
};

/**
  @param atom1 - the first atom in the angle
  @param atom2 - the second atom in the angle
  @param value - the value of the angle in degrees
  @param variable - if the angle between the atoms can change
  @return - instance of the new Angle
*/
Angle createAngle(int atom1, int atom2, double value, bool variable);

/**
  Structure representing a dihedral in the atom.  A dihedral is the angle created
  by two planes.  The structure is defined using two atoms.  It is up to the programmer
  to find the other two atoms needed to define two planes.
*/
struct Dihedral
{
    /**
    the first atom in the dihedral
    */
    int atom1; 
    
    /**
      the second atom in the dihedral
    */
    int atom2; 
    
    /**
    the distance between the atoms
    */
    double value; 
    
    /**
    if the distance between atoms is variable
    */
    bool variable; 
};

/**
  @param atom1 - the first atom in the dihedral
  @param atom2 - the second atom in the dihedral
  @param value - the distance between the atoms
  @param variable - if the dihedral is variable
  @return - instance of the new dihedral
*/
Dihedral createDihedral(int atom1, int atom2, double value, bool variable);

/**
  Atom pairs and their node distance(hops) away from each other
  used in the fudge factor, for total energy calculations
*/
struct Hop
{
    /**
    the starting atom
    */
    int atom1;
    /**
     the ending atom
    */
	int atom2;
    /**
    the number of nodes between start and finish
    */
	int hop;
};

/**
  @param atom1 - the starting atom
  @param atom2 - the ending atom
  @param hops - the number of nodes between the start and finish
  @return - an instance of the Hop structure.
*/
Hop createHop(int atom1, int atom2, int hops);

/**
	Used in box consturction. The tables are used to look up hop values between
	two atom pairs. A table is created for each different type of molecule in the simulation
*/
struct Table
{

    int **hopTable;

};

/**
	@param table - takes a pointer to a multidimentional array to create the table
	@return - an instance of the Table structure.
*/

Table * createTable(int **table);

/**
  Structure to represent a molecule in the simulation.
*/
struct Molecule
{
    /**
    the name of the molecule
    */
    int id;
    /**
    array of atoms in the molecule
    */
    Atom *atoms;
    /**
    array of bonds of the atoms in the molecule.
    */
    Bond *bonds;
    /**
    angles in the molecule between atoms
    */
    Angle *angles;
    /**
    array of dihedrals in the molecule
    */
    Dihedral *dihedrals;
    /**
    array containing a list of atoms that are less than 4 nodes away
    */
    Hop *hops;
    /**
     the number of atoms in the molecule
    */
    int numOfAtoms;
    /**
    the number of bonds in the molecule
    */
    int numOfBonds;
    /**
    the number of angles in the molecule
    */
    int numOfAngles;
    /**
     the number of dihedrals in the atom
    */
    int numOfDihedrals;
    /**
    the number of Hops or pairs of atoms that are less than 4 nodes away
    */
	 int numOfHops;
};

/**
    @param id - the integer id of the molecule
    @param atoms - an array of the atoms in the molecule
    @param angles - an array of the angles in the molecule
    @param bonds - an array of the bonds in the molecule
    @pararm dihedrals - array of dihedrals in the molecule
    @param atomCount - the number of atoms in the molecule
    @param angleCount - the number of angles in the molecule
    @param bondCount - the number of bonds in the molecule
    @param dihedralCount - the number of dihedrals in the molecule
    @return - an instance of the new molecule
*/
Molecule createMolecule(int id,
                        Atom *atoms, Angle *angles, Bond *bonds, Dihedral *dihedrals,
                        int atomCount, int angleCount, int bondCount, int dihedralCount );

/**
    @param id - the integer id of the molecule
    @param atoms - an array of the atoms in the molecule
    @param angles - an array of the angles in the molecule
    @param bonds - an array of the bonds in the molecule
    @param dihedrals - array of dihedrals in the molecule
    @param hops - array of hops in the molecule
    @param atomCount - the number of atoms in the molecule
    @param angleCount - the number of angles in the molecule
    @param bondCount - the number of bonds in the molecule
    @param dihedralCount - the number of dihedrals in the molecule
    @param hopsCount - the number of hops in the molecule
    @return - an instance of the new molecule
*/
Molecule createMolecule(int id,
                        Atom *atoms, Angle *angles, Bond *bonds, Dihedral *dihedrals, Hop *hops,
                        int atomCount, int angleCount, int bondCount, int dihedralCount, int hopsCount );

								
/**
  @param id - the integer id of the molecule
  @param atoms - an array of the atoms in the molecule
  @param atomCount - the number of atoms in the molecule
  @return - an instance of the new molecule
*/
Molecule createMolecule(int id, 
                        Atom *atoms,
                        int atomCount);

/**
  Copies by value molec2 into molec1
  @param molec1 - the destination molecule
  @param molec2 - the source molecule
 */
void copyMolecule(Molecule *molec1, Molecule *molec2);

/**
  Prints the position of all the atoms in an array.
  @param atoms - array of atoms
  @param count - number of atoms
*/
void printAtoms(Atom *atoms, int count);

/**
  @param atoms - array of atoms to be written to the file
  @param enviro - structure holding environmental data.
  @param filename - the name of the file that is to be written
*/
void writePDB(Atom *atoms, Environment enviro, string filename);

/**
  @param molecules - array of molecules to be written to the file
  @param enviro - structure holding environmental data.
  @param filename - the name of the file that is to be written
*/
void writePDB(Molecule *molecules, Environment enviro, string filename);

#define DEFAULT 0
#define START 1
#define END 2
#define OPLS 3
#define Z_MATRIX 4
#define GEOM 5
/**
  Logs output to the OutputLog file
  @param text - the text to be written to the output log
  @param stamp - optional flag used to add the specified tag see above for MACROSs
*/
void writeToLog(string text, int stamp=0 );

/**
  Wrapper that takes in a stringstream instead of a string.
  flushes the stringstream after writing to the OutputLog file
  @param ss - the stringstream to be written to the output log and flushed
  @param stamp - optional flag used to add the specified tag see above for MACROS
*/
void writeToLog(stringstream& ss, int stamp=0 );

/**
  Prints a molecule and its fields
  @param *molec - the molecule to be printed
*/
void printMolecule(Molecule *molec);

/**
  Calulates the percent difference between two doubles by finding the percentage
  of the difference between the two over their adverage.
  @param d1 - the first double
  @param d2 - the second double
  @return - true if the percent difference is less than 3%, else returns false
*/
bool percentDifference(double d1, double d2);

/**
  Checks to see if two boolean values are equal.
  @param b1 - the first boolean value
  @param b2 - the second boolean value
  @return - returns true if they are the same, else returns false
*/
//bool asserTwoBool(bool b1, bool b2);


#endif //METROUTIL_H
