/*
	New version of GPUSimBox
	Serves as a wrapper for SimBox

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

#include "ParallelBox.cuh"

using namespace std;

//Constructor & Destructor
ParallelBox::ParallelBox(): Box()
{
	
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

int ParallelBox::rollback(int moleno)
{
	Box::rollback(moleno);
	writeChangeToDevice(moleno);
	
	return moleno;
}

void ParallelBox::copyDataToDevice()
{
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
	
	moleculesH = new MoleculeData(molecules, moleculeCount);
	cudaMalloc(&atomsIdxD, moleculeCount * sizeof(int));
	cudaMalloc(&numOfAtomsD, moleculeCount * sizeof(int));
	cudaMemcpy(atomsIdxD, moleculesH->atomsIdx, moleculeCount * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(numOfAtomsD, moleculesH->numOfAtoms, moleculeCount * sizeof(int), cudaMemcpyHostToDevice);
	
	MoleculeData *tempMD = (MoleculeData*) malloc(sizeof(MoleculeData));
	tempMD->atomsIdx = atomsIdxD;
	tempMD->numOfAtoms = numOfAtomsD;
	tempMD->moleculeCount = moleculesH->moleculeCount;
	cudaMalloc(&moleculesD, sizeof(MoleculeData));
	cudaMemcpy(moleculesD, tempMD, sizeof(MoleculeData), cudaMemcpyHostToDevice);
	
	nbrMolsH = (int*) malloc(moleculeCount * sizeof(int));
	molBatchH = (int*) malloc(moleculeCount * sizeof(int));
	
	cudaMalloc(&(environmentD), sizeof(Environment));
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
	
	energyCount = moleculesH->moleculeCount * maxMolSize * maxMolSize;
	cudaMalloc(&(energiesD), energyCount * sizeof(Real));
	
	//initialize energies
	cudaMemset(energiesD, 0, sizeof(Real));
	
	//copy data to device
	cudaMemcpy(environmentD, environment, sizeof(Environment), cudaMemcpyHostToDevice);
}

void ParallelBox::writeChangeToDevice(int changeIdx)
{
	//This is temporary until we convert the host data structures as well.
	//TEMP START
	//update AtomData atomsH (MoleculeData will not change)
	int startIdx = moleculesH->atomsIdx[changeIdx];
	for (int i = 0; i < molecules[changeIdx].numOfAtoms; i++)
	{
		atomsH->x[startIdx + i] = molecules[changeIdx].atoms[i].x;
		atomsH->y[startIdx + i] = molecules[changeIdx].atoms[i].y;
		atomsH->z[startIdx + i] = molecules[changeIdx].atoms[i].z;
		//atomsH->sigma[startIdx + i] = molecules[changeIdx].atoms[i].sigma;
		//atomsH->epsilon[startIdx + i] = molecules[changeIdx].atoms[i].epsilon;
		//atomsH->charge[startIdx + i] = molecules[changeIdx].atoms[i].charge;
	}
	//TEMP FINISH

	//copy changed atom data to device
	cudaMemcpy(xD + startIdx, atomsH->x + startIdx, moleculesH->numOfAtoms[changeIdx] * sizeof(Real), cudaMemcpyHostToDevice);
	cudaMemcpy(yD + startIdx, atomsH->y + startIdx, moleculesH->numOfAtoms[changeIdx] * sizeof(Real), cudaMemcpyHostToDevice);
	cudaMemcpy(zD + startIdx, atomsH->z + startIdx, moleculesH->numOfAtoms[changeIdx] * sizeof(Real), cudaMemcpyHostToDevice);
	//cudaMemcpy(sigmaD + startIdx, atomsH->sigma + startIdx, moleculesH->numOfAtoms[changeIdx] * sizeof(Real), cudaMemcpyHostToDevice);
	//cudaMemcpy(epsilonD + startIdx, atomsH->epsilon + startIdx, moleculesH->numOfAtoms[changeIdx] * sizeof(Real), cudaMemcpyHostToDevice);
	//cudaMemcpy(chargeD + startIdx, atomsH->charge + startIdx, moleculesH->numOfAtoms[changeIdx] * sizeof(Real), cudaMemcpyHostToDevice);
}