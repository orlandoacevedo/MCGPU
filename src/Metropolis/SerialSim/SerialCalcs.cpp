/*
	Contains the methods required to calculate energies serially.

	Created: February 21, 2014
	
	-> February 26, by Albert Wallace
	-> March 28, by Joshua Mosby
	-> April 21, by Nathan Coleman
*/

#include <math.h>
#include <string>
#include "Metropolis/DataTypes.h"
#include "Metropolis/SimulationArgs.h"
#include "Metropolis/Utilities/FileUtilities.h"
#include "SerialCalcs.h"

using namespace std;

Box* SerialCalcs::createBox(std::string inputPath, InputFileType inputType, long* startStep, long* steps)
{
	SerialBox* box = new SerialBox();
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
	return (Box*) box;
}

Real SerialCalcs::calcSystemEnergy(Molecule *molecules, Environment *enviro)
{
    Real totalEnergy = 0;

	//for each molecule
	for (int mol = 0; mol < enviro->numOfMolecules; mol++)
	{
		totalEnergy += calcMolecularEnergyContribution(molecules, enviro, mol, mol);
	}
	
    return totalEnergy;
}

Real SerialCalcs::calcMolecularEnergyContribution(Molecule *molecules, Environment *environment, int currentMol, int startIdx)
{
    Real totalEnergy = 0;
	
	//for every other molecule
	#pragma omp parallel for //num_threads(4) <- this is set in Simulation.cpp
	for (int otherMol = startIdx; otherMol < environment->numOfMolecules; otherMol++)
	{
		if (otherMol != currentMol)
		{
			Atom atom1 = molecules[currentMol].atoms[environment->primaryAtomIndex];
			Atom atom2 = molecules[otherMol].atoms[environment->primaryAtomIndex];
			
			//square cutoff value for easy comparison
			Real cutoffSQ = environment->cutoff * environment->cutoff;
				
			//calculate difference in coordinates
			Real deltaX = makePeriodic(atom1.x - atom2.x, environment->x);
			Real deltaY = makePeriodic(atom1.y - atom2.y, environment->y);
			Real deltaZ = makePeriodic(atom1.z - atom2.z, environment->z);
			
			Real r2 = (deltaX * deltaX) +
						(deltaY * deltaY) + 
						(deltaZ * deltaZ);

			if (r2 < cutoffSQ)
			{
				Real tempEnergy = calcInterMolecularEnergy(molecules, currentMol, otherMol, environment);
				//this addition needs to be atomic since multiple threads will be modifying totalEnergy
				#pragma omp atomic
				totalEnergy += tempEnergy;
			}
		}
	}
	return totalEnergy;
}

Real SerialCalcs::calcInterMolecularEnergy(Molecule *molecules, int mol1, int mol2, Environment *enviro)
{
	Real totalEnergy = 0;
	
	for (int i = 0; i < molecules[mol1].numOfAtoms; i++)
	{
		Atom atom1 = molecules[mol1].atoms[i];
		
		for (int j = 0; j < molecules[mol2].numOfAtoms; j++)
		{
			Atom atom2 = molecules[mol2].atoms[j];
		
			if (atom1.sigma >= 0 && atom1.epsilon >= 0 && atom2.sigma >= 0 && atom2.epsilon >= 0)
			{
				//calculate difference in coordinates
				Real deltaX = makePeriodic(atom1.x - atom2.x, enviro->x);
				Real deltaY = makePeriodic(atom1.y - atom2.y, enviro->y);
				Real deltaZ = makePeriodic(atom1.z - atom2.z, enviro->z);
				
				Real r2 = (deltaX * deltaX) +
					 (deltaY * deltaY) + 
					 (deltaZ * deltaZ);
				
				totalEnergy += calc_lj(atom1, atom2, r2);
				totalEnergy += calcCharge(atom1.charge, atom2.charge, sqrt(r2));
			}
		}
		
	}
	return totalEnergy;
}

Real SerialCalcs::calc_lj(Atom atom1, Atom atom2, Real r2)
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
    	const Real sig2OverR2 = pow(sigma, 2) / r2;
		const Real sig6OverR6 = pow(sig2OverR2, 3);
    	const Real sig12OverR12 = pow(sig6OverR6, 2);
    	const Real energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);
        return energy;
    }

}

Real SerialCalcs::calcCharge(Real charge1, Real charge2, Real r)
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

Real SerialCalcs::makePeriodic(Real x, Real boxDim)
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

Real SerialCalcs::calcBlending(Real d1, Real d2)
{
    return sqrt(d1 * d2);
}
