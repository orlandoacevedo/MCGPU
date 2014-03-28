/*
	Contains all structs and methods for simulation objects (atom,molecule,etc.)

	Author: Nathan Coleman
	Last Changed: March 10, 2014 by Albert Wallace
	Previously Changed: Feb 19, Feb 26, Feb 28, Mar 03, Mar 05, Mar 06, Mar 07. [By Albert]
*/

#ifndef STRUCTLIBRARY_H
#define STRUCTLIBRARY_H

#include <string>
#include "../DataTypes.h" //use this statement on OS X, maybe for individual testing
//#include "Metropolis/DataTypes.h" //use this statement on Linux (?), or whenever the makefile is being applied
			///this should be included to allow all DOUBLE *and* FLOAT value types to be dynamically determined; ...
			/// ...redefine each DOUBLE *and* FLOAT using the new type defined in the DataTypes file.


//Forward declaration so that each can be used in methods below

/**
	Used in box consturction. The tables are used to look up hop values between
	two atom pairs. A table is created for each different type of molecule in the simulation
*/
struct Table
{

    int **hopTable;
    
    Table()
    {
    }
    
    Table(int **table) //constructor with parameters
	{
	hopTable = table;
	}
};


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
    double lengthOfBond;
    /*!
       Bool indicating if the bond can be varied.
    */
    bool canBeVaried;
    
    Bond()
    	{
    	atom1 = 0;
    	atom2 = 0;
    	lengthOfBond = 0;
    	canBeVaried = true;
    	}
    	
    Bond(int atom1, int atom2, double distance, bool variable)
	{
    atom1 = atom1;
    atom2 = atom2;
    lengthOfBond = distance;
    canBeVaried = variable;
	}

};


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
    the angle between the atoms; used to be "value"
    */
    double magnitudeOfAngle; 
    
    /*!
    if the angle is variable
    */
    bool canBeVaried;
    
    Angle()
    	{
    	atom1 = 0;
    	atom2 = 0;
    	magnitudeOfAngle = 0;
    	canBeVaried = true;
    	}
};

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
    the distance between the atoms; used to be "value"
    */
    double distanceBetweenAtoms; 
    
    /**
    if the distance between atoms is variable
    */
    bool canBeVaried; 
    
    Dihedral()
    	{
    	atom1 = 0;
    	atom2 = 0;
    	distanceBetweenAtoms = 0;
    	canBeVaried = true;
    	}
};


/**
  Atom pairs and their node distance(hops) away from each other
  used in the fudge factor, for total energy calculations
*/
//This HOP structure is only used as a placeholder; we don't do anything with it.
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
	
	Hop() {
		atom1 = 0;
		atom2 = 0;
		hop = 0;
	}
	Hop(int atom1InAsInt, int atom2InAsInt, int numNodesBetweenStartAndFinish)
	{
		atom1 = atom1InAsInt;
		atom2 = atom2InAsInt;
		hop = numNodesBetweenStartAndFinish;
	}
};

/**
  Structure used to represent the 4 Fourier Coeficients
*/
struct Fourier
{
    double vValues[4];
};

struct Atom
{
	char name;
	double x, y, z, sigma, epsilon, charge;
	unsigned long id;
	Atom(){
		name = '0';
		x = 0;
		y = 0;
		z = 0;
		sigma = 0;
		epsilon = 0;
		charge = 0;
		id = 0;
	}
	
	Atom(unsigned long inID, double inX, double inY, double inZ, double inSigma, double inEpsilon, double inCharge, char inName)
	{
		x = inX;
		y = inY;
		z = inZ;
		sigma = inSigma;
		epsilon = inEpsilon;
		charge = inCharge;
		name = inName;
        id = inID;
	}
	
};

struct Environment
{
	Real x, y, z, cutoff, temp, maxTranslation, maxRotation;
	int numOfAtoms;
	int numOfMolecules; //this line was added in by Albert to make IOUtilities compile
	int primaryAtomIndex;

	int randomseed; //--Albert
	
	Environment() //constructor/initialize all values to 0 or some other default, where applicable
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
		cutoff = 0.0;
		temp = 0.0;
		maxTranslation = 0.0;
		maxRotation = 0.0;
		numOfAtoms = 0;
		numOfMolecules = 0;
		primaryAtomIndex = 0;
		randomseed = 0;
	}

};

struct Molecule
{
	/*
	The number of atoms in the molecule.
	*/
	int numOfAtoms;
	/*
	The identification number for a given molecule. Previously referred to as simply "ID".
	*/
	int moleculeIdentificationNumber;
	/*
	The array representing the collection of atoms in the molecule.
	*/
	Atom *atoms;
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
    
    /*
    Constructor(s) for the molecule
    */
    Molecule(int id, Atom *atomsIn, Angle *anglesIn, Bond *bondsIn, Dihedral *dihedralsIn, Hop *hopsIn, int atomCount, int angleCount, int bondCount, int dihedralCount, int hopCount) 
    	{
		moleculeIdentificationNumber = id;

		atoms = atomsIn;
		angles = anglesIn;
		bonds = bondsIn;
		dihedrals = dihedralsIn;
		hops = hopsIn;

		numOfAtoms = atomCount;
		numOfAngles = angleCount;
		numOfBonds = bondCount;
		numOfDihedrals = dihedralCount;
		numOfHops = hopCount; 	
		}
		
	Molecule() {
		moleculeIdentificationNumber = 0;

		atoms = new Atom();
		angles = new Angle();
		bonds = new Bond();
		dihedrals = new Dihedral();
		hops = new Hop();

		numOfAtoms = 0;
		numOfAngles = 0;
		numOfBonds = 0;
		numOfDihedrals = 0;
		numOfHops = 0; 	
		}
    
};


//Atom
Atom createAtom(unsigned long id, Real x, Real y, Real z);
Atom createAtom(unsigned long id, Real x, Real y, Real z, Real sigma, Real epsilon);
Atom createAtom(unsigned long id, Real x, Real y, Real z, Real sigma, Real epsilon, Real charge, char name);
void printAtoms(Atom *atoms, int count);
void writeOutAtoms(Atom *atoms, Environment *environment, std::string filename, int accepts, int rejects, Real totalEnergy);

//Environment
Environment createEnvironment(Real x, Real y, Real z, Real maxTranslation, Real temp, int numAtoms, Real cutoff, Real maxRotation);

//Molecule
Molecule createMolecule(int id, Atom *atoms, int atomCount);
void copyMolecule(Molecule *destination, Molecule *source);
void printMolecule(Molecule *molecule);

#endif