/*!\file
  \Class for parallel Simulation, including Energy calculate and points to molecules,only save all states
  \author David(Xiao Zhang).
 
  This file contains implement of SimBox that are used to handle environments and common function
  for box.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "parallelSim.cuh"

#define NO 0
#define YES 1

#define MAX_WARP 32
#define MOL_BLOCK 256
#define BATCH_BLOCK 512
#define AGG_BLOCK 512

ParallelSim::ParallelSim(GPUSimBox *initbox,int initsteps)
{
	box = initbox;
	steps = initsteps;
	currentEnergy = 0;
	oldEnergy = 0;
	accepted = 0;
	rejected = 0;
	
	ptrs = (SimPointers*) malloc(sizeof(SimPointers));
	
	ptrs->innerbox = box->getSimBox();
	ptrs->envH = ptrs->innerbox->getEnviro();
	
	ptrs->atomsH = ptrs->innerbox->getAtoms();
	ptrs->moleculesH = ptrs->innerbox->getMolecules();
	
	ptrs->numA = ptrs->envH->numOfAtoms;
	ptrs->numM = ptrs->envH->numOfMolecules;
	
	ptrs->molTrans = (Molecule*) malloc(ptrs->numM * sizeof(Molecule));
	ptrs->nbrMolsH = (int*) malloc(ptrs->numM * sizeof(int));
	ptrs->molBatchH = (int*) malloc(ptrs->numM * sizeof(int));
	
	cudaMalloc(&(ptrs->envD), sizeof(Environment));
	cudaMalloc(&(ptrs->atomsD), ptrs->numA * sizeof(Atom));
	cudaMalloc(&(ptrs->moleculesD), ptrs->numM * sizeof(Molecule));
	cudaMalloc(&(ptrs->nbrMolsD), ptrs->numM * sizeof(int));
	cudaMalloc(&(ptrs->molBatchD), ptrs->numM * sizeof(int));
	
	//sets up device molecules for transfer copies host molecules exactly except
	//for *atoms, which is translated to GPU pointers calculated here
	Atom *a = ptrs->atomsD;
	//upper bound on number of atoms in any molecule
	ptrs->maxMolSize = 0;
	for (int i = 0; i < ptrs->numM; i++)
	{
		ptrs->molTrans[i].atoms = a;
		ptrs->molTrans[i].numOfAtoms = ptrs->moleculesH[i].numOfAtoms;
		a += ptrs->moleculesH[i].numOfAtoms;
		
		if (ptrs->moleculesH[i].numOfAtoms > ptrs->maxMolSize)
		{
			ptrs->maxMolSize = ptrs->moleculesH[i].numOfAtoms;
		}
	}
	
	ptrs->numEnergies = ptrs->numM * ptrs->maxMolSize * ptrs->maxMolSize;
	cudaMalloc(&(ptrs->energiesD), ptrs->numEnergies * sizeof(double));
	
	//initialize energies
	cudaMemset(ptrs->energiesD, 0, sizeof(double));
	
	//copy data to device
	cudaMemcpy(ptrs->envD, ptrs->envH, sizeof(Environment), cudaMemcpyHostToDevice);
	cudaMemcpy(ptrs->atomsD, ptrs->atomsH, ptrs->numA * sizeof(Atom), cudaMemcpyHostToDevice);
	cudaMemcpy(ptrs->moleculesD, ptrs->molTrans, ptrs->numM * sizeof(Molecule), cudaMemcpyHostToDevice);
}

ParallelSim::~ParallelSim()
{
    /*if (energySum_host!=NULL)
    {
        free(energySum_host);
        energySum_host=NULL;
    }
  
    if (energySum_device!=NULL)
    {
        cudaFree(energySum_device);
        energySum_device=NULL;
    }*/
}

void ParallelSim::writeChangeToDevice(int changeIdx)
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

double ParallelSim::calcSystemEnergy()
{
	double totalEnergy = 0;
	
	//for each molecule
	for (int mol = 0; mol < ptrs->numM; mol++)
	{
		totalEnergy += calcMolecularEnergyContribution(mol, mol);
	}

    return totalEnergy;
}

double ParallelSim::calcMolecularEnergyContribution(int molIdx, int startIdx)
{
	return calcBatchEnergy(createMolBatch(molIdx, startIdx), molIdx);
}

int ParallelSim::createMolBatch(int curentMol, int startIdx)
{
	//initialize neighbor molecule slots to NO
	cudaMemset(ptrs->nbrMolsD, NO, ptrs->numM * sizeof(int));
	
	checkMoleculeDistances<<<ptrs->numM / MOL_BLOCK + 1, MOL_BLOCK>>>(ptrs->moleculesD, curentMol, startIdx, ptrs->numM, ptrs->envD, ptrs->nbrMolsD);
	
	cudaMemcpy(ptrs->nbrMolsH, ptrs->nbrMolsD, ptrs->numM * sizeof(int), cudaMemcpyDeviceToHost);
	
	memset(ptrs->molBatchH, -1, ptrs->numM * sizeof(int));
	
	int batchSize = 0;
	
	for (int i = startIdx; i < ptrs->numM; i++)
	{
		if (ptrs->nbrMolsH[i] == YES)
		{
			ptrs->molBatchH[batchSize++] = i;
		}
	}
	
	return batchSize;
}

double ParallelSim::calcBatchEnergy(int numMols, int molIdx)
{
	if (numMols > 0)
	{
		//initialize energies to 0
		cudaMemset(ptrs->energiesD, 0, sizeof(double));
		
		cudaMemcpy(ptrs->molBatchD, ptrs->molBatchH, ptrs->numM * sizeof(int), cudaMemcpyHostToDevice);
		
		calcInterAtomicEnergy<<<ptrs->numEnergies / BATCH_BLOCK + 1, BATCH_BLOCK>>>
		(ptrs->moleculesD, molIdx, ptrs->envD, ptrs->energiesD, ptrs->numEnergies, ptrs->molBatchD, ptrs->maxMolSize);
		
		return getEnergyFromDevice();
	}
	else
	{
		return 0;
	}
}

double ParallelSim::getEnergyFromDevice()
{
	double totalEnergy = 0;
	
	//a batch size of 3 seems to offer the best tradeoff
	int batchSize = 3, blockSize = AGG_BLOCK;
	int numBaseThreads = ptrs->numEnergies / (batchSize);
	for (int i = 1; i < ptrs->numEnergies; i *= batchSize)
	{
		if (blockSize > MAX_WARP && numBaseThreads / i + 1 < blockSize)
		{
			blockSize /= 2;
		}
		aggregateEnergies<<<numBaseThreads / (i * blockSize) + 1, blockSize>>>
		(ptrs->energiesD, ptrs->numEnergies, i, batchSize);
	}
	
	cudaMemcpy(&totalEnergy, ptrs->energiesD, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemset(ptrs->energiesD, 0, sizeof(double));
	
	return totalEnergy;
}

double ParallelSim::makePeriodicH(double x, double box)
{
    
    while(x < -0.5 * box)
    {
        x += box;
    }

    while(x > 0.5 * box)
    {
        x -= box;
    }

    return x;

}

__global__ void checkMoleculeDistances(Molecule *molecules, int currentMol, int startIdx, int numM, Environment *enviro, int *inCutoff)
{
	int otherMol = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (otherMol < numM  && otherMol >= startIdx && otherMol != currentMol)
	{
		Atom atom1 = molecules[currentMol].atoms[enviro->primaryAtomIndex];
		Atom atom2 = molecules[otherMol].atoms[enviro->primaryAtomIndex];
			
		//calculate difference in coordinates
		double deltaX = makePeriodic(atom1.x - atom2.x, enviro->x);
		double deltaY = makePeriodic(atom1.y - atom2.y, enviro->y);
		double deltaZ = makePeriodic(atom1.z - atom2.z, enviro->z);
	  
		double r2 = (deltaX * deltaX) +
					(deltaY * deltaY) + 
					(deltaZ * deltaZ);

		if (r2 < enviro->cutoff * enviro->cutoff)
		{
			inCutoff[otherMol] = YES;
		}
	}
}

__global__ void calcInterAtomicEnergy(Molecule *molecules, int currentMol, Environment *enviro, double *energies, int numEnergies, int *molBatch, int maxMolSize)
{
	int energyIdx = blockIdx.x * blockDim.x + threadIdx.x, segmentSize = maxMolSize * maxMolSize;
	
	if (energyIdx < numEnergies and molBatch[energyIdx / segmentSize] != -1)
	{
		Molecule mol1 = molecules[currentMol], mol2 = molecules[molBatch[energyIdx / segmentSize]];
		int x = (energyIdx % segmentSize) / maxMolSize, y = (energyIdx % segmentSize) % maxMolSize;
		
		if (x < mol1.numOfAtoms && y < mol2.numOfAtoms)
		{
			Atom atom1 = mol1.atoms[x], atom2 = mol2.atoms[y];
		
			if (atom1.sigma >= 0 && atom1.epsilon >= 0 && atom2.sigma >= 0 && atom2.epsilon >= 0)
			{
				double totalEnergy = 0;
			  
				//calculate distance between atoms
				double deltaX = makePeriodic(atom1.x - atom2.x, enviro->x);
				double deltaY = makePeriodic(atom1.y - atom2.y, enviro->y);
				double deltaZ = makePeriodic(atom1.z - atom2.z, enviro->z);
				
				double r2 = (deltaX * deltaX) +
					 (deltaY * deltaY) + 
					 (deltaZ * deltaZ);
				
				totalEnergy += calc_lj(atom1, atom2, r2);
				totalEnergy += calcCharge(atom1.charge, atom2.charge, sqrt(r2));
				
				energies[energyIdx] = totalEnergy;
			}
		}
	}
}

__global__ void aggregateEnergies(double *energies, int numEnergies, int interval, int batchSize)
{
	int idx = batchSize * interval * (blockIdx.x * blockDim.x + threadIdx.x), i;
	
	for (i = 1; i < batchSize; i++)
	{
		if (idx + i * interval < numEnergies)
		{
			energies[idx] += energies[idx + i * interval];
			energies[idx + i * interval] = 0;
		}
	}
}

__device__ double calc_lj(Atom atom1, Atom atom2, double r2)
{
    //store LJ constants locally
    double sigma = calcBlending(atom1.sigma, atom2.sigma);
    double epsilon = calcBlending(atom1.epsilon, atom2.epsilon);
    
    if (r2 == 0.0)
    {
        return 0.0;
    }
    else
    {
    	//calculate terms
    	const double sig2OverR2 = (sigma*sigma) / r2;
		const double sig6OverR6 = (sig2OverR2*sig2OverR2*sig2OverR2);
    	const double sig12OverR12 = (sig6OverR6*sig6OverR6);
    	const double energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);
        return energy;
    }
}

__device__ double calcCharge(double charge1, double charge2, double r)
{  
    if (r == 0.0)
    {
        return 0.0;
    }
    else
    {
    	// conversion factor below for units in kcal/mol
    	const double e = 332.06;
        return (charge1 * charge2 * e) / r;
    }
}

__device__ double makePeriodic(double x, double box)
{
    
    while(x < -0.5 * box)
    {
        x += box;
    }

    while(x > 0.5 * box)
    {
        x -= box;
    }

    return x;

}

__device__ double calcBlending(double d1, double d2)
{
    return sqrt(d1 * d2);
}

__device__ int getXFromIndex(int idx)
{
    int c = -2 * idx;
    int discriminant = 1 - 4 * c;
    int qv = (-1 + sqrtf(discriminant)) / 2;
    return qv + 1;
}

__device__ int getYFromIndex(int x, int idx)
{
    return idx - (x * x - x) / 2;
}

void ParallelSim::runParallel(int steps)
{	
    double temperature = ptrs->envH->temperature;
    double kT = kBoltz * temperature;
    double newEnergyCont, oldEnergyCont;
	
    if (oldEnergy == 0)
	{
		oldEnergy = calcSystemEnergy();
	}
	
    for(int move = 0; move < steps; move++)
    {
        int changeIdx = ptrs->innerbox->chooseMolecule();
		
		oldEnergyCont = calcMolecularEnergyContribution(changeIdx);
		
		ptrs->innerbox->changeMolecule(changeIdx);
		writeChangeToDevice(changeIdx);
		
		newEnergyCont = calcMolecularEnergyContribution(changeIdx);

        bool accept = false;

        if(newEnergyCont < oldEnergyCont)
        {
            accept = true;
        }
        else
        {
            double x = exp(-(newEnergyCont - oldEnergyCont) / kT);

            if(x >= randomFloat(0.0, 1.0))
            {
                accept = true;
            }
            else
            {
                accept = false;
            }
        }
		
        if(accept)
        {
            accepted++;
            oldEnergy += newEnergyCont - oldEnergyCont;
        }
        else
        {
            rejected++;
            //restore previous configuration
            ptrs->innerbox->Rollback(changeIdx);
			writeChangeToDevice(changeIdx);
        }
    }
    currentEnergy=oldEnergy;
}