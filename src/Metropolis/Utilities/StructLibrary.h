/*
	Contains all structs and methods for simulation objects (atom,molecule,etc.)

	Author: Nathan Coleman
	Last Changed: February 19, 2014
*/

#ifndef STRUCTLIBRARY_H
#define STRUCTLIBRARY_H

#include <string>

//Forward declaration so that each can be used in methods below
struct Atom
{
	double x, y, z, sigma, epsilon, charge;
	unsigned long id;
};

struct Environment
{
	double x, y, z, cutoff;
	int numOfMolecules;
	int primaryAtomIndex;
};

struct Molecule
{
	int numAtoms;
	Atom *atoms;
};

//Atom
Atom createAtom(unsigned long id, double x, double y, double z);
Atom createAtom(unsigned long id, double x, double y, double z, double sigma, double epsilon);
Atom createAtom(unsigned long id, double x, double y, double z, double sigma, double epsilon, double charge, char name);
void printAtoms(Atom *atoms, int count);
void writeOutAtoms(Atom *atoms, Environment *environment, std::string filename, int accepts, int rejects, double totalEnergy);

//Environment
Environment createEnvironment(double x, double y, double z, double maxTrans, double temp, int numAtoms, double cutoff, double maxRot);

//Molecule
Molecule createMolecule(int id, Atom *atoms, int atomCount);
void copyMolecule(Molecule *destination, Molecule *source);
void printMolecule(Molecule *molecule);

#endif