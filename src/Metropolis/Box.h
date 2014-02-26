#ifndef BOX_H
#define BOX_H

#include "Metropolis/Utilities/StructLibrary.h"

#define FREE(ptr) if(ptr!=NULL) { free(ptr);ptr=NULL;}

class Box
{
	private:
		Atom *atoms;
		Environment *environment;
		Molecule *molecules;
		double *energies;
		int atomCount, moleculeCount, energyCount, maxMolSize;

	public:
		Box(){atoms=NULL;environment=NULL;molecules=NULL;energies=NULL;};
		~Box(){FREE(atoms);FREE(environment);FREE(molecules);FREE(energies);};
		Atom *getAtoms(){return atoms;};
		Environment *getEnvironment(){return environment;};
		Molecule *getMolecules(){return molecules;};
		double *getEnergies(){return energies;};
		int *getAtomCount(){return atomCount;};
		int *getMoleculeCount(){return moleculeCount;};
		int *getEnergyCount(){return energyCount;};
		int *getMaxMolSize(){return maxMolSize;};
};

#endif