/*
	Contains all structs and methods for simulation objects (atom,molecule,etc.)

	Author: Nathan Coleman
	Last Changed: February 26, 2014 by Albert Wallace
	Previously Changed: February 19, 2014
*/

#ifndef STRUCTLIBRARY_H
#define STRUCTLIBRARY_H

#include <string>

//Forward declaration so that each can be used in methods below
struct Atom
{
	std::string name; //this line was added in by albert to make IOUtilities compile
	double x, y, z, sigma, epsilon, charge;
	unsigned long id;
	
};

struct Environment
{
	double x, y, z, cutoff, temp, maxTranslation, maxRotation;
	int numAtoms;
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
		numAtoms = 0;
		numOfMolecules = 0;
		primaryAtomIndex = 0;
		randomseed = 0;
	}

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
Environment createEnvironment(double x, double y, double z, double maxTranslation, double temp, int numAtoms, double cutoff, double maxRotation);

//Molecule
Molecule createMolecule(int id, Atom *atoms, int atomCount);
void copyMolecule(Molecule *destination, Molecule *source);
void printMolecule(Molecule *molecule);

#endif