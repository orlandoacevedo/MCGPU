/*
	New Simulation to replace linearSim and parallelSim

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

#include "Simulation.h"

using namespace std;

//Constructor & Destructor
Simulation::Simulation(){}
Simulation::~Simulation(){}

//Utility
void Simulation::run()
{
	box=initbox;
    steps=initsteps;
 	currentEnergy=0;
    oldEnergy=0;
    accepted=0;
    rejected=0;
	
	ptrs = (*SimPointers) malloc(sizeof(SimPointers));
	
	ptrs->innerbox = box->getSimBox();
	ptrs->envH = innerbox->getEnviro();
	
	ptrs->atomsH = innerbox->getAtoms();
	ptrs->moleculesH = innerbox->getMolecules();
	
	ptrs->numA = envH->numOfAtoms;
	ptrs->numM = envH->numOfMolecules;
	ptrs->numEnergies = numA * numA;//This is an upper bound. May be able to be tightened.
	
	ptrs->molTrans = (*Molecule) malloc(numM * sizeof(Molecule));
	ptrs->energiesH = (*double) malloc(ptrs->numEnergies * sizeof(double));
	
	cudaMalloc(&(ptrs->envD), sizeof(Environment));
	cudaMalloc(&(ptrs->atomsD), ptrs->numA * sizeof(Atom));
	cudaMalloc(&(ptrs->moleculesD), ptrs->numM * sizeof(Molecule));
	cudaMalloc(&(ptrs->energiesD), ptrs->numEnergies * sizeof(double));
	
	int i;
	
	//initialize energies
	for (i = 0; i < numEnergies; i++)
	{
		ptrs->energiesH[i] = 0;
	}
	
	//sets up device molecules for transfer copies host molecules exactly except
	//for *atoms, which is translated to GPU pointers calculated here
	Atom *a = ptrs->atomsD;
	//upper bound on number of atoms in any molecule
	ptrs->maxMolSize = 0;
	for (i = 0; i < numM; i++)
	{
		ptrs->molTrans[i].atoms = a;
		ptrs->molTrans[i].numOfAtoms = ptrs->moleculesH[i].numOfAtoms;
		a += ptrs->moleculesH[i].numOfAtoms;
		
		if (ptrs->moleculesH[i].numOfAtoms > ptrs->maxMolSize)
		{
			ptrs->maxMolSize = ptrs->moleculesH[i].numOfAtoms;
		}
	}
	
	//copy data to device
	cudaMemcpy(ptrs->envD, ptrs->envH, sizeof(Environment), cudaMemcpyHostToDevice);
	cudaMemcpy(ptrs->atomsD, ptrs->atomsH, ptrs->numA * sizeof(Atom), cudaMemcpyHostToDevice);
	cudaMemcpy(ptrs->moleculesD, ptrs->molTrans, ptrs->numM * sizeof(Molecule), cudaMemcpyHostToDevice);
	cudaMemcpy(ptrs->energiesD, ptrs->energiesH, ptrs->numEnergies * sizeof(double), cudaMemcpyHostToDevice);
}