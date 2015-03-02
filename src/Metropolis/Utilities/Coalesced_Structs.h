#ifndef COALESCED_STRUCTS_H
#define COALESCED_STRUCTS_H

struct AtomData
{
	Real *x, *y, *z;
	Real *sigma, *epsilon, *charge;
	int atomCount;
	
	AtomData(Atom *atoms, int numA)
	{
		x = (Real*) malloc(numA * sizeof(Real));
		y = (Real*) malloc(numA * sizeof(Real));
		z = (Real*) malloc(numA * sizeof(Real));
		sigma = (Real*) malloc(numA * sizeof(Real));
		epsilon = (Real*) malloc(numA * sizeof(Real));
		charge = (Real*) malloc(numA * sizeof(Real));
		
		for (int i = 0; i < numA; i++)
		{
			x[i] = atoms[i].x;
			y[i] = atoms[i].y;
			z[i] = atoms[i].z;
			sigma[i] = atoms[i].sigma;
			epsilon[i] = atoms[i].epsilon;
			charge[i] = atoms[i].charge;
		}
		
		atomCount = numA;
	}
};

struct MoleculeData
{
	char *type;
	int *atomsIdx, *numOfAtoms;
	int moleculeCount;
	
	MoleculeData(Molecule *molecules, int numM)
	{
		int idx = 0;
		
		type = (char*) malloc(numM * sizeof(char));
		atomsIdx = (int*) malloc(numM * sizeof(int));
		numOfAtoms = (int*) malloc(numM * sizeof(int));
		
		for (int i = 0; i < numM; i++)
		{
			numOfAtoms[i] = molecules[i].numOfAtoms;
			type[i] = molecules[i].type;
			atomsIdx[i] = idx;
			idx += numOfAtoms[i];
		}
		
		moleculeCount = numM;
	}
};

#endif
