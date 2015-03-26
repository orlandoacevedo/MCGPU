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
	int *type;
	int *atomsIdx, *numOfAtoms;
	int moleculeCount;
	int *primaryIndexes;
	int totalPrimaryIndexSize;
	
	MoleculeData(Molecule *molecules, int numM, Environment *enviro)
	{
		int idx = 0;
                
		//int sizeOfPrimaryIndexVector = (*enviro->primaryAtomIndexArray).size();
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

#endif
