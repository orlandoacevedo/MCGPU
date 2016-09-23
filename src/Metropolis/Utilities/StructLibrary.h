/** Contains all structs and methods for simulation objects */

#ifndef STRUCTLIBRARY_H
#define STRUCTLIBRARY_H

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include "../DataTypes.h"
//Forward declaration so that each can be used in methods below



struct NLC_Node {
  int index; // The molecule index that this node holds.
  struct NLC_Node* next; // The next node in this NLC chain.

  NLC_Node() {};
};

/** This struct contains data about the bond between two atoms. */
struct BondData {
  /** The force constant of this bond. */
  Real kBond;

  /** The equilibrium bond distance of this bond. */
  Real eqBondDist;

  /** Default constructor. Used when creating a map to bondData. */
  BondData() {}

  /**
   * Preferred constructor. Initializes variables.
   * @param init_kBond The value to initialize the bond's force constant to.
   * @param init_eqBondDist The value to initialize the bond's equilibrium
   * distance to.
   */
  BondData(Real init_kBond, Real init_eqBondDist) {
    kBond = init_kBond;
    eqBondDist = init_eqBondDist;
  }
};

/**
 * This struct holds data about the angle formed by three atoms.
 */
struct AngleData {

  /** The force constant of the angle. */
  Real kAngle;

  /**
   * The equilibrium size of the angle formed by the three atoms.
   * This is measured in degrees.
   */
  Real eqAngle;

  /**
   * Default constructor for angleData. Used when creating a map to angleData.
   */
  AngleData() {}

  /**
   * Preferred constructor for angleData. Initializes variables.
   * @param init_kAngle The value to initialize the angle's force constant to.
   * @param init_eqAngle The value to initialize the angle's size to.
   */
  AngleData(Real init_kAngle, Real init_eqAngle) {
    kAngle = init_kAngle;
    eqAngle = init_eqAngle;
  }
};

/**
 * Used in box consturction. The tables are used to look up hop values between
 * two atom pairs. A table is created for each different type of molecule in
 * the simulation
*/
struct Table {

  int **hopTable;

  Table() { }

  Table(int **table) {
    hopTable = table;
  }
};


/** Structure to respresent bonds between atoms in a molecule.  */
struct Bond {
  /** The first atom in the bond. */
  int atom1;

  /** The second atom in the bond. */
  int atom2;

  /** The length of the bond. */
  Real distance;

  /** Bool indicating if the bond can be varied.  */
  bool variable;

  /** The force constant K of the bond. */
  Real forceConstant;

  /** The equilibrium bond distance of this bond. */
  Real eqBondDist;

  /** Parameterless constructor for the Bond struct. */
  Bond() {
    atom1 = 0;
    atom2 = 0;
    distance = 0;
    variable = true;
  }

  /**
   * Standard constructor for the Bond struct.
   * @param a1 The first atom in the Bond.
   * @param a2 The second atom in thee Bond.
   * @param distance The bond distacne.
   * @param varied True if the bond can be varied, false otherwise.
   */
  Bond(int a1, int a2, Real dist, bool varied) {
    atom1 = a1;
    atom2 = a2;
    distance = dist;
    variable = varied;
  }
};


/**
 * Structure representing an angle between two atoms in a molecule.
 * A geometric angle is defined using three points(atoms).  Only the end
 * points are listed here.  It is up to the programmer to find the third
 * atom in the angle.
 */
struct Angle {
    /** The first atom in the angle */
    int atom1;

    /** The second atom in the angle */
    int atom2;

    /** * The index of the common atom in the angle.  */
    int commonAtom;

    /** The angle between the atoms; used to be "value" */
    Real value;

    /** True if the angle is variable */
    bool variable;

  /** The equilibrium angle of this set of atoms. */
  Real eqAngle;

  /**
   * The force constant k of this angle.
   */
  Real forceConstant;

  Angle() {
    atom1 = 0;
    atom2 = 0;
    commonAtom = 0;
    value = 0;
    variable = true;
  }


  Angle(int a1, int a2, Real val, bool varied) {
    atom1 = a1;
    atom2 = a2;
    commonAtom = 0;
    value = val;
    variable = varied;
  }
};

/**
 * Structure representing a dihedral in the atom.  A dihedral is the angle
 * created by two planes.  The structure is defined using two atoms. 
 * It is up to the programmer to find the other two atoms needed to define two
 * planes.
 */
struct Dihedral {
    /** The first atom in the dihedral */
    int atom1;

    /** The second atom in the dihedral */
    int atom2;

    /** The distance between the atoms; used to be "value" */
    Real value;

    /** True if the distance between atoms is variable */
    bool variable;

    Dihedral() {
      atom1 = 0;
      atom2 = 0;
      value = 0;
      variable = true;
  }


    Dihedral(int a1, int a2, Real val, bool varied) {
      atom1 = a1;
      atom2 = a2;
      value = val;
      variable = varied;
  }
};


/**
  Atom pairs and their node distance(hops) away from each other
  used in the fudge factor, for total energy calculations
*/
//This HOP structure is only used as a placeholder; we don't do anything with it.
struct Hop {
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

    Hop(int atom1In, int atom2In, int hopsIn) {
        atom1 = atom1In;
        atom2 = atom2In;
        hop = hopsIn;
    }
};

/**
  Structure used to represent the 4 Fourier Coeficients
*/
struct Fourier {
    Real vValues[4];
};

struct Atom {
  std::string *name;
  Real x, y, z, sigma, epsilon, charge;
  unsigned long id;
  //std::map<std::string, BondData> bonds;
  //std::map<std::string, std::map<std::string, AngleData> > angles;

  Atom()
    {
    name = new std::string("000");
    x = 0;
    y = 0;
    z = 0;
    sigma = 0;
    epsilon = 0;
    charge = 0;
    id = 0;
  }

  Atom(unsigned long inID, Real inX, Real inY, Real inZ, Real inSigma = 0, Real inEpsilon = 0, Real inCharge = 0, std::string inName = "")
  {
    x = inX;
    y = inY;
    z = inZ;
    sigma = inSigma;
    epsilon = inEpsilon;
    charge = inCharge;
    name = new std::string(inName);
        id = inID;
  }


};

struct Environment
{
  Real x, y, z, cutoff, temp, maxTranslation, maxRotation;
  Real maxBondDelta, maxAngleDelta;
  int numOfAtoms;
  int numOfMolecules; //this line was added in by Albert to make IOUtilities compile
  std::string* primaryAtomIndexConfigLine;
  int primaryAtomIndexDefinitions;
  int primaryAtomIndex;
  std::vector< std::vector<int>* >* primaryAtomIndexArray;
  int randomseed;

  /**
   * Constructor/initialize all values to 0 or some other default, where
   * applicable
   */
  Environment() {
    x = 0.0;
    y = 0.0;
    z = 0.0;
    cutoff = 0.0;
    temp = 0.0;
    maxTranslation = 0.0;
    maxRotation = 0.0;
    numOfAtoms = 0;
    numOfMolecules = 0;
    primaryAtomIndexConfigLine = new std::string("0");
    primaryAtomIndexDefinitions = 0;
    primaryAtomIndex = 0;
    primaryAtomIndexArray = new std::vector< std::vector<int>* >;
    randomseed = 0;

    // Taken from Chrys Woods' slide deck
    maxBondDelta = 0.06;
    maxAngleDelta = 3.1;
  }

  Environment(Environment* environment) {
    if (!environment) {
      std::cerr << "Error: Environment(): Copy constructor given NULL argument"
                << std::endl;
    }

    x = environment->x;
    y = environment->y;
    z = environment->z;
    cutoff = environment->cutoff;
    temp = environment->temp;
    maxTranslation = environment->maxTranslation;
    maxRotation = environment->maxRotation;
    maxBondDelta = environment->maxBondDelta;
    maxAngleDelta = environment->maxAngleDelta;
    numOfAtoms = environment->numOfAtoms;
    numOfMolecules = environment->numOfMolecules;
    primaryAtomIndexConfigLine = environment->primaryAtomIndexConfigLine;
    primaryAtomIndexDefinitions = environment->primaryAtomIndexDefinitions;
    primaryAtomIndex = environment->primaryAtomIndex;
    primaryAtomIndexArray = (environment->primaryAtomIndexArray);
    randomseed = environment->randomseed;
  }

};

struct Molecule
{
  int type;
  /*
  The number of atoms in the molecule.
  */
  int numOfAtoms;
  /*
  The identification number for a given molecule. Previously referred to as simply "ID".
  */
  int id;
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
    Molecule(int idIn, int typeIn, Atom *atomsIn, Angle *anglesIn, Bond *bondsIn, Dihedral *dihedralsIn, Hop *hopsIn, int atomCount, int angleCount, int bondCount, int dihedralCount, int hopCount)
    {
    id = idIn;
    type = typeIn;

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

  Molecule()
    {
    id = 0;
    type = 0;

    atoms = new Atom[0];
    angles = new Angle[0];
    bonds = new Bond[0];
    dihedrals = new Dihedral[0];
    hops = new Hop[0];

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
Atom createAtom(unsigned long id, Real x, Real y, Real z, Real sigma, Real epsilon, Real charge, std::string name);
void printAtoms(Atom *atoms, int count);
void writeOutAtoms(Atom *atoms, Environment *environment, std::string filename, int accepts, int rejects, Real totalEnergy);

//Environment
Environment createEnvironment(Real x, Real y, Real z, Real maxTranslation, Real temp, int numAtoms, Real cutoff, Real maxRotation);

//Molecule
Molecule createMolecule(int id, Atom *atoms, int atomCount);
Molecule createMolecule(int id, int type,  Atom *atoms, Angle *angles, Bond *bonds, Dihedral *dihedrals,
                        int atomCount, int angleCount, int bondCount, int dihedralCount);
void copyMolecule(Molecule *destination, Molecule *source);
void printMolecule(Molecule *molecule);

#endif
