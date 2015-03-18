/*
	Implements methods related to managing data between the host and device.
	Subclass of Box.

	Created: February 21, 2014
	
	-> February 26, by Albert Wallace
	-> March 28, by Joshua Mosby
	-> April 21, by Nathan Coleman
*/

#include "ParallelBox.cuh"

using namespace std;

ParallelBox::ParallelBox(): Box()
{
	//Is anything really needed here?
}

ParallelBox::~ParallelBox()
{
	// TODO: free device memory
}

int ParallelBox::changeMolecule(int molIdx)
{
	Box::changeMolecule(molIdx);
	writeChangeToDevice(molIdx);
	
	return molIdx;
}

int ParallelBox::rollback(int molIdx)
{
	Box::rollback(molIdx);
	writeChangeToDevice(molIdx);
	
	return molIdx;
}

void ParallelBox::copyDataToDevice()
{
	//create AtomData on host, and fill atomic data arrays on device
	atomsH = new AtomData(atoms, atomCount);
	cudaMalloc(&xD, atomCount * sizeof(Real));
	cudaMalloc(&yD, atomCount * sizeof(Real));
	cudaMalloc(&zD, atomCount * sizeof(Real));
	cudaMalloc(&sigmaD, atomCount * sizeof(Real));
	cudaMalloc(&epsilonD, atomCount * sizeof(Real));
	cudaMalloc(&chargeD, atomCount * sizeof(Real));
	cudaMemcpy(xD, atomsH->x, atomCount * sizeof(Real), cudaMemcpyHostToDevice);
	cudaMemcpy(yD, atomsH->y, atomCount * sizeof(Real), cudaMemcpyHostToDevice);
	cudaMemcpy(zD, atomsH->z, atomCount * sizeof(Real), cudaMemcpyHostToDevice);
	cudaMemcpy(sigmaD, atomsH->sigma, atomCount * sizeof(Real), cudaMemcpyHostToDevice);
	cudaMemcpy(epsilonD, atomsH->epsilon, atomCount * sizeof(Real), cudaMemcpyHostToDevice);
	cudaMemcpy(chargeD, atomsH->charge, atomCount * sizeof(Real), cudaMemcpyHostToDevice);
	
	//create device AtomData struct with pointers to filled-in atomic data arrays
	AtomData *tempAD = (AtomData*) malloc(sizeof(AtomData));
	tempAD->x = xD;
	tempAD->y = yD;
	tempAD->z = zD;
	tempAD->sigma = sigmaD;
	tempAD->epsilon = epsilonD;
	tempAD->charge = chargeD;
	tempAD->atomCount = atomsH->atomCount;
	cudaMalloc(&atomsD, sizeof(AtomData));
	cudaMemcpy(atomsD, tempAD, sizeof(AtomData), cudaMemcpyHostToDevice);
	
	//create MoleculeData on host, and fill molecular data arrays on device
//	printf("TotalPrimaryIndexSize: %d\n", moleculesH->totalPrimaryIndexSize);
	moleculesH = new MoleculeData(molecules, moleculeCount, environment);
	cudaMalloc(&atomsIdxD, moleculeCount * sizeof(int));
	cudaMalloc(&numOfAtomsD, moleculeCount * sizeof(int));
	cudaMalloc(&typeD, moleculeCount * sizeof(int));
	cudaMalloc(&primaryIndexesD, moleculesH->totalPrimaryIndexSize * sizeof(int));
	cudaMemcpy(atomsIdxD, moleculesH->atomsIdx, moleculeCount * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(numOfAtomsD, moleculesH->numOfAtoms, moleculeCount * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(typeD, moleculesH->type, moleculeCount * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(primaryIndexesD, moleculesH->primaryIndexes, moleculesH->totalPrimaryIndexSize * sizeof(int), cudaMemcpyHostToDevice);
	
	//create device MoleculeData struct with pointers to filled-in molecular data arrays
	MoleculeData *tempMD = (MoleculeData*) malloc(sizeof(MoleculeData));
	tempMD->atomsIdx = atomsIdxD;
	tempMD->numOfAtoms = numOfAtomsD;
	tempMD->type = typeD;
	tempMD->primaryIndexes = primaryIndexesD;
	tempMD->moleculeCount = moleculesH->moleculeCount;
	tempMD->totalPrimaryIndexSize = moleculesH->totalPrimaryIndexSize;
	cudaMalloc(&moleculesD, sizeof(MoleculeData));
	cudaMemcpy(moleculesD, tempMD, sizeof(MoleculeData), cudaMemcpyHostToDevice);
	
	//data structures for neighbor batch in energy calculation
	nbrMolsH = (int*) malloc(moleculeCount * sizeof(int));
	molBatchH = (int*) malloc(moleculeCount * sizeof(int));
	cudaMalloc(&(nbrMolsD), moleculeCount * sizeof(int));
	cudaMalloc(&(molBatchD), moleculeCount * sizeof(int));
	
	//upper bound on number of atoms in any molecule
	maxMolSize = 0;
	for (int i = 0; i < moleculesH->moleculeCount; i++)
	{
		if (moleculesH->numOfAtoms[i] > maxMolSize)
		{
			maxMolSize = moleculesH->numOfAtoms[i];
		}
	}
	
	//energies array on device has one segment for each molecule
	//where each segment has the maximum number of
	//possible interatomic energies for one pair of molecules
	energyCount = moleculesH->moleculeCount * maxMolSize * maxMolSize;
	cudaMalloc(&(energiesD), energyCount * sizeof(Real));
	
	//initialize energies to 0
	cudaMemset(energiesD, 0, energyCount * sizeof(Real));
	
	//copy Environment to device
	cudaMalloc(&(environmentD), sizeof(Environment));
	cudaMemcpy(environmentD, environment, sizeof(Environment), cudaMemcpyHostToDevice);
}

void ParallelBox::writeChangeToDevice(int changeIdx)
{
	//update AtomData atomsH (MoleculeData will not change)
	int startIdx = moleculesH->atomsIdx[changeIdx];
	for (int i = 0; i < molecules[changeIdx].numOfAtoms; i++)
	{
		atomsH->x[startIdx + i] = molecules[changeIdx].atoms[i].x;
		atomsH->y[startIdx + i] = molecules[changeIdx].atoms[i].y;
		atomsH->z[startIdx + i] = molecules[changeIdx].atoms[i].z;
		//sigma, epsilon, and charge will not change, so there is no need to update those arrays
	}

	//copy changed atom data to device
	cudaMemcpy(xD + startIdx, atomsH->x + startIdx, moleculesH->numOfAtoms[changeIdx] * sizeof(Real), cudaMemcpyHostToDevice);
	cudaMemcpy(yD + startIdx, atomsH->y + startIdx, moleculesH->numOfAtoms[changeIdx] * sizeof(Real), cudaMemcpyHostToDevice);
	cudaMemcpy(zD + startIdx, atomsH->z + startIdx, moleculesH->numOfAtoms[changeIdx] * sizeof(Real), cudaMemcpyHostToDevice);
	//sigma, epsilon, and charge will not change, so there is no need to update those arrays
}
