#ifndef STRUCTLIBRARY_CPP
#define STRUCTLIBRARY_CPP

/*
  Contains all structs and methods for simulation objects (atom,molecule,etc.)

  Author: Nathan Coleman
  Last Changed: February 27, 2014 by Albert Wallace
  Previously Changed: February 19; February 26, 2014
*/

#include "StructLibrary.h"
#include "Metropolis/DataTypes.h" //AlbertExcludes

using namespace std;

//Atom
Atom createAtom(unsigned long id, Real x, Real y, Real z)
{
  Atom atom;
  atom.id = id;
  atom.x = x;
  atom.y = y;
  atom.z = z;
  return atom;
}

Atom createAtom(unsigned long id, Real x, Real y, Real z, Real sigma, Real epsilon)
{
  Atom atom;
  atom.id = id;
  atom.x = x;
  atom.y = y;
  atom.z = z;
  atom.sigma = sigma;
  atom.epsilon = epsilon;
  return atom;
}

Atom createAtom(unsigned long id, Real x, Real y, Real z, Real sigma, Real epsilon, Real charge, std::string name)
{
  Atom atom;
  atom.id = id;
  atom.x = x;
  atom.y = y;
  atom.z = z;
  atom.sigma = sigma;
  atom.epsilon = epsilon;
  atom.charge = charge;
  *atom.name = name;
  return atom;  
}

void printAtoms(Atom *atoms, int count){}

void writeOutAtoms(Atom *atoms, Environment *environment, std::string filename, int accepts, int rejects, Real totalEnergy){}

//Environment
Environment createEnvironment(Real x, Real y, Real z, Real maxTranslation,
                              Real temp, int numOfAtoms, Real cutoff,
                              Real maxRotation, Real maxBondDelta,
                              Real maxAngleDelta) {
  Environment environment;
  environment.x = x;
  environment.y = y;
  environment.z = z;
  environment.maxTranslation = maxTranslation;
  environment.maxRotation = maxRotation;
  environment.maxBondDelta = maxBondDelta;
  environment.maxAngleDelta = maxAngleDelta;
  environment.cutoff = cutoff;
  environment.temp = temp;
  environment.numOfAtoms = numOfAtoms;
  return environment;
}

//Molecule
Molecule createMolecule(Atom *atoms, int numOfAtoms)
{
  Molecule molecule;
  molecule.atoms = atoms;
  molecule.numOfAtoms = numOfAtoms;
  return molecule;
}

Molecule createMolecule(int id, int type, Atom *atoms, Angle *angles, Bond *bonds, Dihedral *dihedrals, 
                        int atomCount, int angleCount, int bondCount, int dihedralCount)
{
    Molecule molecule;
    molecule.id = id;
    molecule.type = type;

    molecule.atoms = atoms;
    molecule.angles = angles;
    molecule.bonds = bonds;
    molecule.dihedrals = dihedrals;

    molecule.numOfAtoms = atomCount;
    molecule.numOfAngles = angleCount;
    molecule.numOfBonds = bondCount;
    molecule.numOfDihedrals = dihedralCount;

    return molecule;
}

void copyMolecule(Molecule *destination, Molecule *source){}

void printMolecule(Molecule *molecule){}


#endif
