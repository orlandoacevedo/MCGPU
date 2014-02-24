#include "ParallelCalcs.h"
#include "ParallelCalcs.cuh"
#include <stdio.h>

using namespace std;

void calcCont(Molecule *mol){}

int main()
{
	printf("hello\n");
	return 0;
}

__global__ void calcInterMolecularEnergy(Molecule *molecules, int currentMol, int numM, Environment *enviro, double *energies, int segmentSize)
{
	int otherMol = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (otherMol < numM && otherMol != currentMol)
	{
		Atom atom1 = molecules[currentMol].atoms[enviro->primaryAtomIndex];
		Atom atom2 = molecules[otherMol].atoms[enviro->primaryAtomIndex];
			
		//calculate distance between atoms
		double deltaX = makePeriodic(atom1.x - atom2.x, enviro->x);
		double deltaY = makePeriodic(atom1.y - atom2.y, enviro->y);
		double deltaZ = makePeriodic(atom1.z - atom2.z, enviro->z);
		
		double r2 = (deltaX * deltaX) +
					(deltaY * deltaY) + 
					(deltaZ * deltaZ);

		if (r2 < enviro->cutoff * enviro->cutoff)
		{
			//calculate intermolecular energies
			calcInterAtomicEnergy
			<<<molecules[currentMol].numAtoms, molecules[otherMol].numAtoms>>>
			(molecules, currentMol, otherMol, enviro, energies, segmentSize);
		}
	}
}

__global__ void calcInterAtomicEnergy(Molecule *molecules, int currentMol, int otherMol, Environment *enviro, double *energies, int segmentSize)
{
	Atom atom1 = molecules[currentMol].atoms[blockIdx.x],
		 atom2 = molecules[otherMol].atoms[threadIdx.x];
	//energyIdx will need to take into account different mol sizes (in atoms)
	int energyIdx = otherMol * segmentSize + blockIdx.x * blockDim.x + threadIdx.x;
	
	if (!(currentMol == otherMol && threadIdx.x == blockIdx.x) && atom1.sigma >= 0 && atom1.epsilon >= 0 && atom2.sigma >= 0 && atom2.epsilon >= 0)
	{
		double totalEnergy = 0;
		
		//calculate difference in coordinates
		double deltaX = atom1.x - atom2.x;
		double deltaY = atom1.y - atom2.y;
		double deltaZ = atom1.z - atom2.z;
	  
		//calculate distance between atoms
		deltaX = makePeriodic(deltaX, enviro->x);
		deltaY = makePeriodic(deltaY, enviro->y);
		deltaZ = makePeriodic(deltaZ, enviro->z);
		
		double r2 = (deltaX * deltaX) +
			 (deltaY * deltaY) + 
			 (deltaZ * deltaZ);
			
		//gets the fValue if in the same molecule
		double fvalue = 1.0;
		
		totalEnergy += calc_lj(atom1, atom2, r2) * fvalue;
		totalEnergy += calcCharge(atom1.charge, atom2.charge, sqrt(r2)) * fvalue;
		
		energies[energyIdx] = totalEnergy;
	}
}

__global__ void calcIntraMolecularEnergy(Molecule *molecules, int currentMol, int numE, Environment *enviro, double *energies, int segmentSize)
{
	Molecule cMol = molecules[currentMol];
	int energyIdx = blockIdx.x * blockDim.x + threadIdx.x,
		x = getXFromIndex(energyIdx);
	Atom atom1 = cMol.atoms[x], atom2 = cMol.atoms[getYFromIndex(x, energyIdx)];
	
	if (energyIdx < numE)
	{
		energyIdx += currentMol * segmentSize;
		
		double totalEnergy = 0;
			
		//calculate difference in coordinates
		double deltaX = atom1.x - atom2.x;
		double deltaY = atom1.y - atom2.y;
		double deltaZ = atom1.z - atom2.z;
	  
		//calculate distance between atoms
		deltaX = makePeriodic(deltaX, enviro->x);
		deltaY = makePeriodic(deltaY, enviro->y);
		deltaZ = makePeriodic(deltaZ, enviro->z);
		
		double r2 = (deltaX * deltaX) +
			 (deltaY * deltaY) + 
			 (deltaZ * deltaZ);
			
		//gets the fValue if in the same molecule
		double fvalue = 1.0;
		
		totalEnergy += calc_lj(atom1, atom2, r2) * fvalue;
		totalEnergy += calcCharge(atom1.charge, atom2.charge, sqrt(r2)) * fvalue;
		
		energies[energyIdx] = totalEnergy;
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
    	const double sig2OverR2 = pow(sigma, 2) / r2;
		const double sig6OverR6 = pow(sig2OverR2, 3);
    	const double sig12OverR12 = pow(sig6OverR6, 2);
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