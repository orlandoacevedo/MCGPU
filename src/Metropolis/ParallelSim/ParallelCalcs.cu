/*
	Contains the methods required to calculate energies in parallel.

	Created: February 21, 2014
	
	-> February 26, by Albert Wallace
	-> March 28, by Joshua Mosby
	-> April 21, by Nathan Coleman
*/

#include "ParallelCalcs.h"
#include "ParallelCalcs.cuh"
#include "ParallelBox.cuh"
#include <string>
#include "Metropolis/Utilities/FileUtilities.h"
#include "Metropolis/Box.h"
#include "Metropolis/SimulationArgs.h"
#include <thrust/reduce.h>
#include <thrust/count.h>
#include <thrust/remove.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>

#define NO -1

#define MAX_WARP 32
#define MOL_BLOCK 256
#define BATCH_BLOCK 512
#define AGG_BLOCK 512

using namespace std;

Box* ParallelCalcs::createBox(string inputPath, InputFileType inputType, long* startStep, long* steps)
{
	ParallelBox* box = new ParallelBox();
	if (!loadBoxData(inputPath, inputType, box, startStep, steps))
	{
		if (inputType != InputFile::Unknown)
		{
			std::cerr << "Error: Could not build from file: " << inputPath << std::endl;
			return NULL;
		}
		else
		{
			std::cerr << "Error: Can not build environment with unknown file: " << inputPath << std::endl;
			return NULL;
		}
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
		//use startIdx parameter to prevent double-calculating energies (Ex mols 3->5 and mols 5->3)
		totalEnergy += calcMolecularEnergyContribution(box, mol, mol + 1);
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

struct isThisTrue {
	__device__ bool operator()(const int &x) {
	  return x != NO;
	}
};

int ParallelCalcs::createMolBatch(ParallelBox *box, int currentMol, int startIdx)
{
	//initialize neighbor molecule slots to NO
	cudaMemset(box->nbrMolsD, NO, box->moleculeCount * sizeof(int));
	
	//check molecule distances in parallel, conditionally replacing NO with index value in box->nbrMolsD
	checkMoleculeDistances<<<box->moleculeCount / MOL_BLOCK + 1, MOL_BLOCK>>>(box->moleculesD, box->atomsD, currentMol, startIdx, box->environmentD, box->nbrMolsD);
	
	thrust::device_ptr<int> neighborMoleculesOnDevice = thrust::device_pointer_cast(&box->nbrMolsD[0]);
	thrust::device_ptr<int> moleculesInBatchOnDevice = thrust::device_pointer_cast(&box->molBatchD[0]);
	
	//copy over neighbor molecules that don't have NO as their index value
	thrust::device_ptr<int> lastElementFound = thrust::copy_if(neighborMoleculesOnDevice, neighborMoleculesOnDevice + box->moleculeCount, moleculesInBatchOnDevice, isThisTrue());
	
	return lastElementFound - moleculesInBatchOnDevice;
}

Real ParallelCalcs::calcBatchEnergy(ParallelBox *box, int numMols, int molIdx)
{
	if (numMols <= 0) return 0;
	
	//There will only be as many energy segments filled in as there are molecules in the batch.
	int validEnergies = numMols * box->maxMolSize * box->maxMolSize;
	
	//calculate interatomic energies between changed molecule and all molecules in batch
	calcInterAtomicEnergy<<<validEnergies / BATCH_BLOCK + 1, BATCH_BLOCK>>>
	(box->moleculesD, box->atomsD, molIdx, box->environmentD, box->energiesD, validEnergies, box->molBatchD, box->maxMolSize);
	
	//Using Thrust here for a sum reduction on all of the individual energy contributions in box->energiesD.
	thrust::device_ptr<Real> energiesOnDevice = thrust::device_pointer_cast(&box->energiesD[0]);
	return thrust::reduce(energiesOnDevice, energiesOnDevice + validEnergies, (Real) 0, thrust::plus<Real>());
}

__global__ void ParallelCalcs::checkMoleculeDistances(MoleculeData *molecules, AtomData *atoms, int currentMol, int startIdx, Environment *enviro, int *inCutoff)
{
	int otherMol = blockIdx.x * blockDim.x + threadIdx.x;
	
	//checks validity of molecule pair
	if (otherMol < molecules->moleculeCount && otherMol >= startIdx && otherMol != currentMol)
	{
		bool included = false;    	
		for (int i = 0; i < molecules->totalPrimaryIndexSize; i++)
		{
		    int currentMoleculeIndexCount = molecules->primaryIndexes[i];
		    int currentTypeIndex = i+1;
		    int potentialCurrentMoleculeType = molecules->primaryIndexes[currentTypeIndex];
			
		    if (potentialCurrentMoleculeType == molecules->type[currentMol])
		    {
	    	        int *currentMolPrimaryIndexArray = molecules->primaryIndexes + currentTypeIndex + 1;
		        int currentMolPrimaryIndexArrayLength = currentMoleculeIndexCount - 1;		

			for (int k = 0; k < molecules->totalPrimaryIndexSize; k++)
			{
		    	    int otherMoleculeIndexCount = molecules->primaryIndexes[k];
    			    int otherTypeIndex = k+1;
			    int potentialOtherMoleculeType = molecules->primaryIndexes[otherTypeIndex];
			
			    if (potentialOtherMoleculeType == molecules->type[otherMol])
			    {
				int *otherMolPrimaryIndexArray = molecules->primaryIndexes + otherTypeIndex + 1;
				int otherMolPrimaryIndexArrayLength = otherMoleculeIndexCount - 1;
					
				for (int m = 0; m < currentMolPrimaryIndexArrayLength; m++)
				{
				    for (int n = 0; n < otherMolPrimaryIndexArrayLength; n++)
				    {
					//find primary atom indices for this pair of molecules
					int atom1 = molecules->atomsIdx[currentMol] + *(currentMolPrimaryIndexArray + m);
					int atom2 = molecules->atomsIdx[otherMol] + *(otherMolPrimaryIndexArray + n);
			
					//calculate periodic difference in coordinates
					Real deltaX = makePeriodic(atoms->x[atom1] - atoms->x[atom2], enviro->x);
					Real deltaY = makePeriodic(atoms->y[atom1] - atoms->y[atom2], enviro->y);
					Real deltaZ = makePeriodic(atoms->z[atom1] - atoms->z[atom2], enviro->z);
		
					Real r2 = (deltaX * deltaX) +
							    (deltaY * deltaY) + 
							    (deltaZ * deltaZ);
		
					//if within curoff, write index to inCutoff
					if (r2 < enviro->cutoff * enviro->cutoff)
					{
					    inCutoff[otherMol] = otherMol;
					    included = true;
					    break;
					}	
				    }
				    if (included)
					break;
				}
			    }
			    if (included)
				break;
			    else
			    	k += otherMoleculeIndexCount;
			 }
		    }
		    if (included)
			break;
		    else
			i += currentMoleculeIndexCount;
		}
	
		/*//find primary atom indices for this pair of molecules
		for (int i = 0; i < molecules->totalPrimaryIndexSize; i++)
		{
		    printf("checkMoleculeDistances:totalPrimaryIndexSize: %d Array: ",molecules->totalPrimaryIndexSize);
		    printf("%d: ", i);
		    printf("%d", molecules->primaryIndexes[i]);
		} 
		printf("\n");
		//int atom1 = molecules->atomsIdx[currentMol] + enviro->primaryAtomIndex;
		//int atom2 = molecules->atomsIdx[otherMol] + enviro->primaryAtomIndex;


		int atom1 = molecules->atomsIdx[currentMol] + enviro->primaryAtomIndex;
		int atom2 = molecules->atomsIdx[otherMol] + enviro->primaryAtomIndex;
			
		//calculate periodic difference in coordinates
		Real deltaX = makePeriodic(atoms->x[atom1] - atoms->x[atom2], enviro->x);
		Real deltaY = makePeriodic(atoms->y[atom1] - atoms->y[atom2], enviro->y);
		Real deltaZ = makePeriodic(atoms->z[atom1] - atoms->z[atom2], enviro->z);
		
		Real r2 = (deltaX * deltaX) +
					(deltaY * deltaY) + 
					(deltaZ * deltaZ);
		
		//if within curoff, write index to inCutoff
		if (r2 < enviro->cutoff * enviro->cutoff)
		{
			inCutoff[otherMol] = otherMol;
		}*/
	}
}

__global__ void ParallelCalcs::calcInterAtomicEnergy(MoleculeData *molecules, AtomData *atoms, int currentMol, Environment *enviro, Real *energies, int energyCount, int *molBatch, int maxMolSize)
{
	int energyIdx = blockIdx.x * blockDim.x + threadIdx.x;
	int segmentSize = maxMolSize * maxMolSize;
	
	//check validity of thread
	if (energyIdx < energyCount and molBatch[energyIdx / segmentSize] != NO)
	{
		//get other molecule index
		int otherMol = molBatch[energyIdx / segmentSize];
		
		//get atom pair for this thread
		int x = (energyIdx % segmentSize) / maxMolSize;
		int y = (energyIdx % segmentSize) % maxMolSize;
		
		//check validity of atom pair
		if (x < molecules->numOfAtoms[currentMol] && y < molecules->numOfAtoms[otherMol])
		{
			//get atom indices
			int atom1 = molecules->atomsIdx[currentMol] + x;
			int atom2 = molecules->atomsIdx[otherMol] + y;
			
			//check validity of atoms (ensure they are not dummy atoms)
			if (atoms->sigma[atom1] >= 0 && atoms->epsilon[atom1] >= 0 && atoms->sigma[atom2] >= 0 && atoms->epsilon[atom2] >= 0)
			{
				Real totalEnergy = 0;
			  
				//calculate periodic distance between atoms
				Real deltaX = makePeriodic(atoms->x[atom1] - atoms->x[atom2], enviro->x);
				Real deltaY = makePeriodic(atoms->y[atom1] - atoms->y[atom2], enviro->y);
				Real deltaZ = makePeriodic(atoms->z[atom1] - atoms->z[atom2], enviro->z);
				
				Real r2 = (deltaX * deltaX) +
					 (deltaY * deltaY) + 
					 (deltaZ * deltaZ);
				
				//calculate interatomic energies
				totalEnergy += calc_lj(atoms, atom1, atom2, r2);
				totalEnergy += calcCharge(atoms->charge[atom1], atoms->charge[atom2], sqrt(r2));
				
				//store energy
				energies[energyIdx] = totalEnergy;
			}
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

__device__ Real ParallelCalcs::makePeriodic(Real x, Real boxDim)
{
    
    while(x < -0.5 * boxDim)
    {
        x += boxDim;
    }

    while(x > 0.5 * boxDim)
    {
        x -= boxDim;
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
