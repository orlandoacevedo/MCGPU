#ifndef COALESCED_STRUCTS_H
#define COALESCED_STRUCTS_H
#include "../SerialSim/NeighborList.h"
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
	int *atomsIdx, *numOfAtoms, *type, *primaryIndexes;
	int moleculeCount, totalPrimaryIndexSize;

	MoleculeData(Molecule *molecules, int numM, Environment *enviro)
	{
		int idx = 0;
		int sizeOfPrimaryIndexVector = enviro->primaryAtomIndexDefinitions;
 		int totalPrimaryIndexCount = 0;
		
		for (int i = 0; i < sizeOfPrimaryIndexVector; i++)
		{
		    totalPrimaryIndexCount += (*(*enviro->primaryAtomIndexArray)[i]).size();
		}		
		totalPrimaryIndexCount += 2 * sizeOfPrimaryIndexVector;

		type = (int*) malloc(numM * sizeof(int));
		atomsIdx = (int*) malloc(numM * sizeof(int));
		numOfAtoms = (int*) malloc(numM * sizeof(int));
		primaryIndexes = (int*) malloc(totalPrimaryIndexCount * sizeof(int));

		for (int i = 0; i < numM; i++)
		{
			numOfAtoms[i] = molecules[i].numOfAtoms;
			type[i] = molecules[i].type;
			atomsIdx[i] = idx;
			idx += numOfAtoms[i];
		}
		
		int index = 0;
		for (int i = 0; i < sizeOfPrimaryIndexVector; i++)
		{
		    primaryIndexes[index] = (*(*enviro->primaryAtomIndexArray)[i]).size() + 1;
		    primaryIndexes[index+1] = i;
		    
		    for(int j = 0; j < (*(*enviro->primaryAtomIndexArray)[i]).size(); j++)
		    {
			primaryIndexes[j + index + 2] = (*(*enviro->primaryAtomIndexArray)[i])[j];
		    }
			
		    index += 2 + (*(*enviro->primaryAtomIndexArray)[i]).size();
		}

		printf("TotalPrimaryIndexSize: %d\n", totalPrimaryIndexCount);
		for (int i = 0; i < totalPrimaryIndexCount; i++)
		{
		    printf("%d ", primaryIndexes[i]);
		}
	
		moleculeCount = numM;
		totalPrimaryIndexSize = totalPrimaryIndexCount;
	}
};

struct NeighborListData
{
	int *numCells, *head, *linkedCellList;
	int nmax, nclmax, empty, numCellsYZ, numCellsXYZ;
	Real *lengthCell, *region;
	Real rrCut;

	NeighborListData(NeighborList *nl)
	{
		int nmaxIn= NMAX;
		int nclmaxIn = NCLMAX;
		int emptyIn = EMPTY;
	
		int numCellsYZIn = nl->numCellsYZ;
		int numCellsXYZIn = nl->numCellsXYZ;
		Real rrCutIn = nl->rrCut;
	
		numCells = (int*) malloc(3 * sizeof(int));
		head = (int*) malloc(NCLMAX * sizeof(int));
		linkedCellList = (int*) malloc(NMAX * sizeof(int));
		lengthCell = (Real*) malloc(3 * sizeof(Real));
		region = (Real*) malloc(3 * sizeof(Real));

		for (int i = 0; i < 3; i++)
		{
			numCells[i] = nl->numCells[i];
			lengthCell[i] = nl->lengthCell[i];
			region[i] = nl->region[i];
		}

		for (int i = 0; i < nclmaxIn; i++)
		{
			head[i] = nl->head[i];
		}
		
		for (int i = 0; i < nmaxIn; i++)
		{
			linkedCellList[i] = nl->linkedCellList[i];
		}

		nmax = nmaxIn;
		nclmax = nclmaxIn;
		empty = emptyIn;
		numCellsYZ = numCellsYZIn;
		numCellsXYZ = numCellsXYZIn;
		rrCut = rrCutIn;
	}

};
#endif
