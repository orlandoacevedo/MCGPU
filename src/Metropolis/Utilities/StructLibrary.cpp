#include "StructLibrary.h"

using namespace std;

//Atom
Atom createAtom(unsigned long id, double x, double y, double z)
{
	Atom atom;
	atom.id = id;
	atom.x = x;
	atom.y = y;
	atom.z = z;
	return atom;
}

Atom createAtom(unsigned long id, double x, double y, double z, double sigma, double epsilon)
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

Atom createAtom(unsigned long id, double x, double y, double z, double sigma, double epsilon, double charge, char name)
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

void writeOutAtoms(Atom *atoms, Environment *environment, std::string filename, int accepts, int rejects, double totalEnergy){}



//Environment
Environment createEnvironment(double x, double y, double z, double maxTranslation, double temp, int numAtoms, double cutoff, double maxRotation)
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
Molecule createMolecule(int id, Atom *atoms, int numAtoms)
{
	Molecule molecule;
	molecule.id = id;
	molecule.atoms = atoms;
	molecule.numAtoms = numAtoms;
	return molecule;
}

void copyMolecule(Molecule *destination, Molecule *source){}

void printMolecule(Molecule *molecule){}