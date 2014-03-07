/*
	Contains calculations for ParallelBox
	Same functions as SerialCalcs with function qualifiers and CUDA code

	Author: Nathan Coleman
*/

#include <stdio.h>
#include <math.h>
#include "Metropolis/DataTypes.h"
#include "ParallelCalcs.h"
#include "ParallelCalcs.cuh"

using namespace std;

Real calcMolecularEnergyContribution(int molIdx, int startIdx)
{
	Real totalEnergy = 0;
	
	//initialize energies to 0
	for (int i = 0; i < ptrs->numEnergies; i++)
	{
		ptrs->energiesH[i] = 0;
	}
	
	cudaMemcpy(ptrs->energiesD, ptrs->energiesH, ptrs->numEnergies * sizeof(Real), cudaMemcpyHostToDevice);
	
	//calculate intermolecular energies (cutoff check for each molecule)
	//using startIdx this way has the potential to waste a significant
	//amount of GPU resources, look into other methods later.
	calcInterMolecularEnergy<<<ptrs->numM / BLOCK_SIZE + 1, BLOCK_SIZE>>>
	(ptrs->moleculesD, molIdx, ptrs->numM, startIdx, ptrs->envD, ptrs->energiesD, ptrs->maxMolSize * ptrs->maxMolSize);
	
	//calculate intramolecular energies for changed molecule
	int numAinM = ptrs->moleculesD[molIdx].numOfAtoms;
	int numIntraEnergies = numAinM * (numAinM - 1) / 2;
	calcIntraMolecularEnergy<<<numIntraEnergies / BLOCK_SIZE + 1, BLOCK_SIZE>>>
	(ptrs->moleculesD, molIdx, numIntraEnergies, ptrs->envD, ptrs->energiesD, ptrs->maxMolSize * ptrs->maxMolSize);
						
	cudaMemcpy(ptrs->energiesH, ptrs->energiesD, ptrs->numEnergies * sizeof(Real), cudaMemcpyDeviceToHost);
	
	for (i = 0; i < numEnergies; i++)
	{
		totalEnergy += energiesH[i];
	}
	
	return totalEnergy;
}

Real ParallelSim::calcSystemEnergy()
{
	Real totalEnergy = 0;

	//for each molecule
	for (int mol = 0; mol < ptrs->numM; mol++)
	{
		totalEnergy += calcMolecularEnergyContribution(mol, mol);
	}
	
    return totalEnergy;
}

__global__ void calcInterMolecularEnergy(Molecule *molecules, int currentMol, int numM, int startIdx, Environment *environment, Real *energies, int segmentSize)
{
	int otherMol = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (otherMol >= startIdx && otherMol < numM && otherMol != currentMol)
	{
		Atom atom1 = molecules[currentMol].atoms[environment->primaryAtomIndex];
		Atom atom2 = molecules[otherMol].atoms[environment->primaryAtomIndex];
			
		//calculate distance between atoms
		Real deltaX = makePeriodic(atom1.x - atom2.x, environment->x);
		Real deltaY = makePeriodic(atom1.y - atom2.y, environment->y);
		Real deltaZ = makePeriodic(atom1.z - atom2.z, environment->z);
		
		Real r2 = (deltaX * deltaX) +
					(deltaY * deltaY) + 
					(deltaZ * deltaZ);

		if (r2 < environment->cutoff * environment->cutoff)
		{
			//calculate intermolecular energies
			calcInterAtomicEnergy
			<<<molecules[currentMol].numOfAtoms, molecules[otherMol].numOfAtoms>>>
			(molecules, currentMol, otherMol, environment, energies, segmentSize);
		}
	}
}

__global__ void calcInterAtomicEnergy(Molecule *molecules, int currentMol, int otherMol, Environment *environment, Real *energies, int segmentSize)
{
	Atom atom1 = molecules[currentMol].atoms[blockIdx.x],
		 atom2 = molecules[otherMol].atoms[threadIdx.x];
	int energyIdx = otherMol * segmentSize + blockIdx.x * blockDim.x + threadIdx.x;
	
	if (!(currentMol == otherMol && threadIdx.x == blockIdx.x) && atom1.sigma >= 0 && atom1.epsilon >= 0 && atom2.sigma >= 0 && atom2.epsilon >= 0)
	{
		Real totalEnergy = 0;
		
		//calculate difference in coordinates
		Real deltaX = atom1.x - atom2.x;
		Real deltaY = atom1.y - atom2.y;
		Real deltaZ = atom1.z - atom2.z;
	  
		//calculate distance between atoms
		deltaX = makePeriodic(deltaX, environment->x);
		deltaY = makePeriodic(deltaY, environment->y);
		deltaZ = makePeriodic(deltaZ, environment->z);
		
		Real r2 = (deltaX * deltaX) +
			 (deltaY * deltaY) + 
			 (deltaZ * deltaZ);
		
		totalEnergy += calc_lj(atom1, atom2, r2);
		totalEnergy += calcCharge(atom1.charge, atom2.charge, sqrt(r2));
		
		energies[energyIdx] = totalEnergy;
	}
}

__global__ void calcIntraMolecularEnergy(Molecule *molecules, int currentMol, int numE, Environment *environment, Real *energies, int segmentSize)
{
	Molecule cMol = molecules[currentMol];
	int energyIdx = blockIdx.x * blockDim.x + threadIdx.x,
		x = getXFromIndex(energyIdx);
	Atom atom1 = cMol.atoms[x], atom2 = cMol.atoms[getYFromIndex(x, energyIdx)];
	
	if (energyIdx < numE)
	{
		energyIdx += currentMol * segmentSize;
		
		Real totalEnergy = 0;
			
		//calculate difference in coordinates
		Real deltaX = atom1.x - atom2.x;
		Real deltaY = atom1.y - atom2.y;
		Real deltaZ = atom1.z - atom2.z;
	  
		//calculate distance between atoms
		deltaX = makePeriodic(deltaX, environment->x);
		deltaY = makePeriodic(deltaY, environment->y);
		deltaZ = makePeriodic(deltaZ, environment->z);
		
		Real r2 = (deltaX * deltaX) +
			 (deltaY * deltaY) + 
			 (deltaZ * deltaZ);
			
		//gets the fValue if in the same molecule
		Real fvalue = 1;
		
		totalEnergy += calc_lj(atom1, atom2, r2) * fvalue;
		totalEnergy += calcCharge(atom1.charge, atom2.charge, sqrt(r2)) * fvalue;
		
		energies[energyIdx] = totalEnergy;
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