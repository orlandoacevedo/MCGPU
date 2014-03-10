#ifndef STRUCTLIBRARY_CPP
#define STRUCTLIBRARY_CPP

/*
	Contains all structs and methods for simulation objects (atom,molecule,etc.)

	Author: Nathan Coleman
	Last Changed: February 27, 2014 by Albert Wallace
	Previously Changed: February 19; February 26, 2014
*/

#include "StructLibrary.h"
#include "Metropolis/DataTypes.h"

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

Atom createAtom(unsigned long id, Real x, Real y, Real z, Real sigma, Real epsilon, Real charge, char name)
{
	Atom atom;
	atom.id = id;
	atom.x = x;
	atom.y = y;
	atom.z = z;
	atom.sigma = sigma;
	atom.epsilon = epsilon;
	atom.charge = charge;
	atom.name = name;
	return atom;	
}

void printAtoms(Atom *atoms, int count){}

void writeOutAtoms(Atom *atoms, Environment *environment, std::string filename, int accepts, int rejects, Real totalEnergy){}



//Environment
Environment createEnvironment(Real x, Real y, Real z, Real maxTranslation, Real temp, int numAtoms, Real cutoff, Real maxRotation)
{
	Environment environment;
	environment.x = x;
	environment.y = y;
	environment.z = z;
	environment.maxTranslation = maxTranslation;
	environment.maxRotation = maxRotation;
	environment.cutoff = cutoff;
	environment.temp = temp;
	environment.numAtoms = numAtoms;
	return environment;
}



//Molecule
Molecule createMolecule(Atom *atoms, int numAtoms)
{
	Molecule molecule;
	molecule.atoms = atoms;
	molecule.numAtoms = numAtoms;
	return molecule;
}

void copyMolecule(Molecule *destination, Molecule *source){}

void printMolecule(Molecule *molecule){}

#endif