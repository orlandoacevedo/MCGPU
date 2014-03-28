/*
	New version of GPUSimBox
	Serves as a wrapper for SimBox

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

#include "ParallelBox.cuh"

using namespace std;

//Constructor & Destructor
ParallelBox::ParallelBox(IOUtilities ioUtil): Box(ioUtil)
{
	//call Box constructor first
	
	copyDataToDevice();
}

ParallelBox::~ParallelBox()
{
	//free device memory
}

void ParallelBox::copyDataToDevice()
{
	transferMoleculesH = (Molecule*) malloc(moleculeCount * sizeof(Molecule));
	nbrMolsH = (int*) malloc(moleculeCount * sizeof(int));
	molBatchH = (int*) malloc(moleculeCount * sizeof(int));
	
	cudaMalloc(&(environmentD), sizeof(Environment));
	cudaMalloc(&(atomsD), atomCount * sizeof(Atom));
	cudaMalloc(&(moleculesD), moleculeCount * sizeof(Molecule));
	cudaMalloc(&(nbrMolsD), moleculeCount * sizeof(int));
	cudaMalloc(&(molBatchD), moleculeCount * sizeof(int));
	
	//sets up device molecules for transfer copies host molecules exactly except
	//for *atoms, which is translated to GPU pointers calculated here
	Atom *a = atomsD;
	//upper bound on number of atoms in any molecule
	maxMolSize = 0;
	for (int i = 0; i < moleculeCount; i++)
	{
		transferMoleculesH[i].atoms = a;
		transferMoleculesH[i].numOfAtoms = moleculesH[i].numOfAtoms;
		a += moleculesH[i].numOfAtoms;
		
		if (moleculesH[i].numOfAtoms > maxMolSize)
		{
			maxMolSize = moleculesH[i].numOfAtoms;
		}
	}
	
	energyCount = moleculeCount * maxMolSize * maxMolSize;
	cudaMalloc(&(energiesD), energyCount * sizeof(double));
	
	//initialize energies
	cudaMemset(energiesD, 0, sizeof(double));
	
	//copy data to device
	cudaMemcpy(environmentD, environment, sizeof(Environment), cudaMemcpyHostToDevice);
	cudaMemcpy(atomsD, atoms, atomCount * sizeof(Atom), cudaMemcpyHostToDevice);
	cudaMemcpy(moleculesD, transferMoleculesH, moleculeCount * sizeof(Molecule), cudaMemcpyHostToDevice);
}

void ParallelBox::writeChangeToDevice(int changeIdx)
{
	//create temp Molecule
	Molecule *changedMol = (Molecule*) malloc(sizeof(Molecule));
	
	//copy changed Molecule into temp Molecule
	//ready to be copied over to device, except that it still contains host pointer in .atoms
	memcpy(changedMol, ptrs->moleculesH + changeIdx, sizeof(Molecule));
	
	//changedMol.atoms will now contain a pointer to Atoms on device
	//this pointer never meant to be followed from host
	changedMol->atoms = ptrs->molTrans[changeIdx].atoms;
	
	//copy changed molecule to device
	cudaMemcpy(ptrs->moleculesD + changeIdx, changedMol, sizeof(Molecule), cudaMemcpyHostToDevice);
	
	//copy changed atoms to device
	Atom *destAtoms = ptrs->molTrans[changeIdx].atoms;
	cudaMemcpy(destAtoms, ptrs->moleculesH[changeIdx].atoms, ptrs->moleculesH[changeIdx].numOfAtoms * sizeof(Atom), cudaMemcpyHostToDevice);
}