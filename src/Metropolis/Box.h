#ifndef BOX_H
#define BOX_H

#include "Metropolis/Utilities/StructLibrary.h"

class Box
{
	Atom *atoms;
	Environment *environment;
	Molecule *molecules;
	double *energies;
	int atomCount, moleculeCount, energyCount, maxMolSize;

	//Utility
	int copyMolecule(Molecule *destination, Molecule *source);
	int saveMolecule(int moleculeIndex);

	Atom *getAtoms(){return atoms;};
	Environment *getEnvironment(){return environment;};
	Molecule *getMolecules(){return molecules;};
};

#endif