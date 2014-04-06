/*
	Contains calculations for ParallelBox
	Same functions as SerialCalcs with function qualifiers and CUDA code

	Author: Nathan Coleman
*/

#include "ParallelCalcs.h"
#include "ParallelCalcs.cuh"
#include "ParallelBox.cuh"
#include <string>
#include "Metropolis/Utilities/FileUtilities.h"
#include "Metropolis/Box.h"

#define NO 0
#define YES 1

#define MAX_WARP 32
#define MOL_BLOCK 256
#define BATCH_BLOCK 512
#define AGG_BLOCK 512

using namespace std;

Box* ParallelCalcs::createBox(string configpath, long* steps)
{
	ParallelBox* box = new ParallelBox();
	if (!loadBoxData(configpath, box, steps))
	{
		std::cerr << "Error: Cannot create ParallelBox from config: " << configpath << std::endl;
		return NULL;
	}
	box->copyDataToDevice();
	return (Box*) box;
}

Real ParallelCalcs::calcSystemEnergy(Box *box)
{
	Real totalEnergy = 0;
	
	//for each molecule
	for (int mol = 0; mol < box->moleculeCount; mol++)
	{
		totalEnergy += calcMolecularEnergyContribution(box, mol, mol);
	}

    return totalEnergy;
}

Real ParallelCalcs::calcMolecularEnergyContribution(Box *box, int molIdx, int startIdx)
{
	ParallelBox *pBox = (ParallelBox*) box;
	
	if (pBox == NULL)
	{
		return 0;
	}
	
	return calcBatchEnergy(pBox, createMolBatch(pBox, molIdx, startIdx), molIdx);
}

int ParallelCalcs::createMolBatch(ParallelBox *box, int currentMol, int startIdx)
{
	//initialize neighbor molecule slots to NO
	cudaMemset(box->nbrMolsD, NO, box->moleculeCount * sizeof(int));
	
	checkMoleculeDistances<<<box->moleculeCount / MOL_BLOCK + 1, MOL_BLOCK>>>(box->moleculesD, box->atomsD, currentMol, startIdx, box->environmentD, box->nbrMolsD);
	
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

Real ParallelCalcs::calcBatchEnergy(ParallelBox *box, int numMols, int molIdx)
{
	if (numMols > 0)
	{
		int validEnergies = numMols * box->maxMolSize * box->maxMolSize;
		//initialize energies to 0
		cudaMemset(box->energiesD, 0, sizeof(Real));
		
		cudaMemcpy(box->molBatchD, box->molBatchH, box->moleculeCount * sizeof(int), cudaMemcpyHostToDevice);
		
		calcInterAtomicEnergy<<<box->energyCount / BATCH_BLOCK + 1, BATCH_BLOCK>>>
		(box->moleculesD, box->atomsD, molIdx, box->environmentD, box->energiesD, validEnergies, box->molBatchD, box->maxMolSize);
		
		return getEnergyFromDevice(box, validEnergies);
	}
	else
	{
		return 0;
	}
}

Real ParallelCalcs::getEnergyFromDevice(ParallelBox *box, int validEnergies)
{
	Real totalEnergy = 0;
	
	//a batch size of 3 seems to offer the best tradeoff
	int batchSize = 3, blockSize = AGG_BLOCK;
	int numBaseThreads = validEnergies / (batchSize);
	for (int i = 1; i < validEnergies; i *= batchSize)
	{
		if (blockSize > MAX_WARP && numBaseThreads / i + 1 < blockSize)
		{
			blockSize /= 2;
		}
		aggregateEnergies<<<numBaseThreads / (i * blockSize) + 1, blockSize>>>
		(box->energiesD, validEnergies, i, batchSize);
	}
	
	cudaMemcpy(&totalEnergy, box->energiesD, sizeof(Real), cudaMemcpyDeviceToHost);
	cudaMemset(box->energiesD, 0, sizeof(Real));
	
	return totalEnergy;
}

__global__ void ParallelCalcs::checkMoleculeDistances(MoleculeData *molecules, AtomData *atoms, int currentMol, int startIdx, Environment *enviro, int *inCutoff)
{
	int otherMol = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (otherMol < molecules->moleculeCount && otherMol >= startIdx && otherMol != currentMol)
	{
		int atom1 = molecules->atomsIdx[currentMol] + enviro->primaryAtomIndex;
		int atom2 = molecules->atomsIdx[otherMol] + enviro->primaryAtomIndex;
			
		//calculate difference in coordinates
		Real deltaX = makePeriodic(atoms->x[atom1] - atoms->x[atom2], enviro->x);
		Real deltaY = makePeriodic(atoms->y[atom1] - atoms->y[atom2], enviro->y);
		Real deltaZ = makePeriodic(atoms->z[atom1] - atoms->z[atom2], enviro->z);
	  
		Real r2 = (deltaX * deltaX) +
					(deltaY * deltaY) + 
					(deltaZ * deltaZ);

		if (r2 < enviro->cutoff * enviro->cutoff)
		{
			inCutoff[otherMol] = YES;
		}
	}
}

__global__ void ParallelCalcs::calcInterAtomicEnergy(MoleculeData *molecules, AtomData *atoms, int currentMol, Environment *enviro, Real *energies, int energyCount, int *molBatch, int maxMolSize)
{
	int energyIdx = blockIdx.x * blockDim.x + threadIdx.x, segmentSize = maxMolSize * maxMolSize;
	
	if (energyIdx < energyCount and molBatch[energyIdx / segmentSize] != -1)
	{
		int otherMol = molBatch[energyIdx / segmentSize];
		int x = (energyIdx % segmentSize) / maxMolSize, y = (energyIdx % segmentSize) % maxMolSize;
		
		if (x < molecules->numOfAtoms[currentMol] && y < molecules->numOfAtoms[otherMol])
		{
			int atom1 = molecules->atomsIdx[currentMol] + x;
			int atom2 = molecules->atomsIdx[otherMol] + y;
		
			if (atoms->sigma[atom1] >= 0 && atoms->epsilon[atom1] >= 0 && atoms->sigma[atom2] >= 0 && atoms->epsilon[atom2] >= 0)
			{
				Real totalEnergy = 0;
			  
				//calculate distance between atoms
				Real deltaX = makePeriodic(atoms->x[atom1] - atoms->x[atom2], enviro->x);
				Real deltaY = makePeriodic(atoms->y[atom1] - atoms->y[atom2], enviro->y);
				Real deltaZ = makePeriodic(atoms->z[atom1] - atoms->z[atom2], enviro->z);
				
				Real r2 = (deltaX * deltaX) +
					 (deltaY * deltaY) + 
					 (deltaZ * deltaZ);
				
				totalEnergy += calc_lj(atoms, atom1, atom2, r2);
				totalEnergy += calcCharge(atoms->charge[atom1], atoms->charge[atom2], sqrt(r2));
				
				energies[energyIdx] = totalEnergy;
			}
		}
	}
}

__global__ void ParallelCalcs::aggregateEnergies(Real *energies, int energyCount, int interval, int batchSize)
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

__device__ Real ParallelCalcs::calc_lj(AtomData *atoms, int atom1, int atom2, Real r2)
{
    //store LJ constants locally
    Real sigma = calcBlending(atoms->sigma[atom1], atoms->sigma[atom2]);
    Real epsilon = calcBlending(atoms->epsilon[atom1], atoms->epsilon[atom2]);
    
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

__device__ Real ParallelCalcs::calcCharge(Real charge1, Real charge2, Real r)
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

__device__ Real ParallelCalcs::makePeriodic(Real x, Real box)
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

__device__ Real ParallelCalcs::calcBlending(Real d1, Real d2)
{
    return sqrt(d1 * d2);
}

__device__ int ParallelCalcs::getXFromIndex(int idx)
{
    int c = -2 * idx;
    int discriminant = 1 - 4 * c;
    int qv = (-1 + sqrtf(discriminant)) / 2;
    return qv + 1;
}

__device__ int ParallelCalcs::getYFromIndex(int x, int idx)
{
    return idx - (x * x - x) / 2;
}