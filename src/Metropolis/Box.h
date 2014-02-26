#ifndef BOX_H
#define BOX_H

#include "Metropolis/Utilities/StructLibrary.h"

class Box
{
	Atom *atomPool;
	Environment *environment;
	Molecule changedMolecule;
	Molecule *molecules;
	int moleculeType;

	//Utility
	int copyMolecule(Molecule *destination, Molecule *source);
	int saveMolecule(int moleculeIndex);

	Atom *getAtoms(){return atomPool;};
	Environment *getEnvironment(){return environment;};
	Molecule *getMolecules(){return molecules;};
};

#endif