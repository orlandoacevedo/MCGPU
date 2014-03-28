/*
	Contains calculations for ParallelBox
	Same functions as SerialCalcs with function qualifiers and CUDA code

	Author: Nathan Coleman
*/

#include "ParallelCalcs.h"
#include "ParallelCalcs.cuh"

using namespace std;

Real ParallelCalcs::calcSystemEnergy(Box box)
{
	Real totalEnergy = 0;
	
	//for each molecule
	for (int mol = 0; mol < box->moleculeCount; mol++)
	{
		totalEnergy += calcMolecularEnergyContribution(box, mol, mol);
	}

    return totalEnergy;
}

Real ParallelCalcs::calcMolecularEnergyContribution(Box box, int molIdx, int startIdx)
{
	return calcBatchEnergy(createMolBatch(box, molIdx, startIdx), molIdx);
}

int ParallelCalcs::createMolBatch(Box box, int currentMol, int startIdx)
{
	//initialize neighbor molecule slots to NO
	cudaMemset(box->nbrMolsD, NO, box->moleculeCount * sizeof(int));
	
	checkMoleculeDistances<<<box->moleculeCount / MOL_BLOCK + 1, MOL_BLOCK>>>(box->moleculesD, currentMol, startIdx, box->moleculeCount, box->environmentD, box->nbrMolsD);
	
	cudaMemcpy(box->nbrMolsH, box->nbrMolsD, box->moleculeCount * sizeof(int), cudaMemcpyDeviceToHost);
	
	memset(box->molBatchH, -1, box->moleculeCount * sizeof(int));
	
	int batchSize = 0;
	
	for (int i = startIdx; i < box->moleculeCount; i++)
	{
		if (box->nbrMolsH[i] == YES)
		{
			box->molBatchH[batchSize++] = i;
		}
	}
	
	return batchSize;
}

Real ParallelCalcs::calcBatchEnergy(Box box, int numMols, int molIdx)
{
	if (numMols > 0)
	{
		//initialize energies to 0
		cudaMemset(box->energiesD, 0, sizeof(Real));
		
		cudaMemcpy(box->molBatchD, box->molBatchH, box->moleculeCount * sizeof(int), cudaMemcpyHostToDevice);
		
		calcInterAtomicEnergy<<<box->energyCount / BATCH_BLOCK + 1, BATCH_BLOCK>>>
		(box->moleculesD, molIdx, box->environmentD, box->energiesD, box->energyCount, box->molBatchD, box->maxMolSize);
		
		return getEnergyFromDevice(box);
	}
	else
	{
		return 0;
	}
}

Real ParallelCalcs::getEnergyFromDevice(Box box)
{
	Real totalEnergy = 0;
	
	//a batch size of 3 seems to offer the best tradeoff
	int batchSize = 3, blockSize = AGG_BLOCK;
	int numBaseThreads = box->energyCount / (batchSize);
	for (int i = 1; i < box->energyCount; i *= batchSize)
	{
		if (blockSize > MAX_WARP && numBaseThreads / i + 1 < blockSize)
		{
			blockSize /= 2;
		}
		aggregateEnergies<<<numBaseThreads / (i * blockSize) + 1, blockSize>>>
		(box->energiesD, box->energyCount, i, batchSize);
	}
	
	cudaMemcpy(&totalEnergy, box->energiesD, sizeof(Real), cudaMemcpyDeviceToHost);
	cudaMemset(box->energiesD, 0, sizeof(Real));
	
	return totalEnergy;
}

__global__ void checkMoleculeDistances(Molecule *molecules, int currentMol, int startIdx, int moleculeCount, Environment *enviro, int *inCutoff)
{
	int otherMol = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (otherMol < moleculeCount  && otherMol >= startIdx && otherMol != currentMol)
	{
		Atom atom1 = molecules[currentMol].atoms[enviro->primaryAtomIndex];
		Atom atom2 = molecules[otherMol].atoms[enviro->primaryAtomIndex];
			
		//calculate difference in coordinates
		Real deltaX = makePeriodic(atom1.x - atom2.x, enviro->x);
		Real deltaY = makePeriodic(atom1.y - atom2.y, enviro->y);
		Real deltaZ = makePeriodic(atom1.z - atom2.z, enviro->z);
	  
		Real r2 = (deltaX * deltaX) +
					(deltaY * deltaY) + 
					(deltaZ * deltaZ);

		if (r2 < enviro->cutoff * enviro->cutoff)
		{
			inCutoff[otherMol] = YES;
		}
	}
}

__global__ void calcInterAtomicEnergy(Molecule *molecules, int currentMol, Environment *enviro, Real *energies, int energyCount, int *molBatch, int maxMolSize)
{
	int energyIdx = blockIdx.x * blockDim.x + threadIdx.x, segmentSize = maxMolSize * maxMolSize;
	
	if (energyIdx < energyCount and molBatch[energyIdx / segmentSize] != -1)
	{
		Molecule mol1 = molecules[currentMol], mol2 = molecules[molBatch[energyIdx / segmentSize]];
		int x = (energyIdx % segmentSize) / maxMolSize, y = (energyIdx % segmentSize) % maxMolSize;
		
		if (x < mol1.numOfAtoms && y < mol2.numOfAtoms)
		{
			Atom atom1 = mol1.atoms[x], atom2 = mol2.atoms[y];
		
			if (atom1.sigma >= 0 && atom1.epsilon >= 0 && atom2.sigma >= 0 && atom2.epsilon >= 0)
			{
				Real totalEnergy = 0;
			  
				//calculate distance between atoms
				Real deltaX = makePeriodic(atom1.x - atom2.x, enviro->x);
				Real deltaY = makePeriodic(atom1.y - atom2.y, enviro->y);
				Real deltaZ = makePeriodic(atom1.z - atom2.z, enviro->z);
				
				Real r2 = (deltaX * deltaX) +
					 (deltaY * deltaY) + 
					 (deltaZ * deltaZ);
				
				totalEnergy += calc_lj(atom1, atom2, r2);
				totalEnergy += calcCharge(atom1.charge, atom2.charge, sqrt(r2));
				
				energies[energyIdx] = totalEnergy;
			}
		}
	}
}

__global__ void aggregateEnergies(Real *energies, int energyCount, int interval, int batchSize)
{
	int idx = batchSize * interval * (blockIdx.x * blockDim.x + threadIdx.x), i;
	
	for (i = 1; i < batchSize; i++)
	{
		if (idx + i * interval < energyCount)
		{
			energies[idx] += energies[idx + i * interval];
			energies[idx + i * interval] = 0;
		}
	}
}

__device__ Real calc_lj(Atom atom1, Atom atom2, Real r2)
{
    //store LJ constants locally
    Real sigma = calcBlending(atom1.sigma, atom2.sigma);
    Real epsilon = calcBlending(atom1.epsilon, atom2.epsilon);
    
    if (r2 == 0.0)
    {
        return 0.0;
    }
    else
    {
    	//calculate terms
    	const Real sig2OverR2 = (sigma*sigma) / r2;
		const Real sig6OverR6 = (sig2OverR2*sig2OverR2*sig2OverR2);
    	const Real sig12OverR12 = (sig6OverR6*sig6OverR6);
    	const Real energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);
        return energy;
    }
}

__device__ Real calcCharge(Real charge1, Real charge2, Real r)
{  
    if (r == 0.0)
    {
        return 0.0;
    }
    else
    {
    	// conversion factor below for units in kcal/mol
    	const Real e = 332.06;
        return (charge1 * charge2 * e) / r;
    }
}

__device__ Real makePeriodic(Real x, Real box)
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

__device__ Real calcBlending(Real d1, Real d2)
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