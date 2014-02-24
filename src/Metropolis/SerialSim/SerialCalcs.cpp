/*
	Contains calculations for SerialBox

	Author: Nathan Coleman
*/

#include <math.h>
#include "SerialCalcs.h"

using namespace std;

void calcContribution(Molecule *mol){}

double calcMolecularEnergyContribution(Molecule *molecules, Environment *enviro, int currentMol, int startIdx = 0)
{
    double totalEnergy = 0;
	
	//for every other molecule
	for (int otherMol = startIdx; otherMol < enviro->numOfMolecules; otherMol++)
	{
		Atom atom1 = molecules[currentMol].atoms[enviro->primaryAtomIndex];
		Atom atom2 = molecules[otherMol].atoms[enviro->primaryAtomIndex];
		
		//square cutoff value for easy comparison
		double cutoffSQ = enviro->cutoff * enviro->cutoff;
			
		//calculate difference in coordinates
		double deltaX = atom1.x - atom2.x;
		double deltaY = atom1.y - atom2.y;
		double deltaZ = atom1.z - atom2.z;
	  
		//calculate distance between atoms
		deltaX = box->makePeriodic(deltaX, enviro->x);
		deltaY = box->makePeriodic(deltaY, enviro->y);
		deltaZ = box->makePeriodic(deltaZ, enviro->z);
		
		double r2 = (deltaX * deltaX) +
					(deltaY * deltaY) + 
					(deltaZ * deltaZ);

		if (r2 < cutoffSQ)
		{
			totalEnergy += calcInterMolecularEnergy(molecules, currentMol, otherMol, enviro);
		}
	}
	return totalEnergy;
}

double calcSystemEnergy(Molecule *molecules, Environment *enviro)
{
    double totalEnergy = 0;

	//for each molecule
	for (int mol = 0; mol < enviro->numOfMolecules; mol++)
	{
		totalEnergy += calcMolecularEnergyContribution(molecules, enviro, mol, mol);
	}
	
    return totalEnergy;
}

double calcInterMolecularEnergy(Molecule *molecules, int currentMol, int numM, Environment *environment, double *energies, int segmentSize)
{
	double totalEnergy = 0;
	
	for (int i = 0; i < molecules[mol1].numOfAtoms; i++)
	{
		Atom atom1 = molecules[mol1].atoms[i];
	
		for (int j = 0; j < molecules[mol2].numOfAtoms; j++)
		{
			Atom atom2 = molecules[mol2].atoms[j];
		
			if (atom1.id <= atom2.id && atom1.sigma >= 0 && atom1.epsilon >= 0 && atom2.sigma >= 0 && atom2.epsilon >= 0)
			{
				//calculate difference in coordinates
				double deltaX = atom1.x - atom2.x;
				double deltaY = atom1.y - atom2.y;
				double deltaZ = atom1.z - atom2.z;
			  
				//calculate distance between atoms
				deltaX = box->makePeriodic(deltaX, enviro->x);
				deltaY = box->makePeriodic(deltaY, enviro->y);
				deltaZ = box->makePeriodic(deltaZ, enviro->z);
				
				double r2 = (deltaX * deltaX) +
					 (deltaY * deltaY) + 
					 (deltaZ * deltaZ);
					
				//gets the fValue if in the same molecule
				double fValue = 1.0;
				if(mol1 == mol2)
				{
					int ** hopTab1 = box->tables[mol1 % box->molecTypenum].hopTable;
					fValue = box->getFValue(i, j, hopTab1);
				}
				
				totalEnergy += calcLJ(atom1, atom2, r2) * fvalue;
				totalEnergy += calcCharge(atom1.charge, atom2.charge, sqrt(r2)) * fValue;
			}
		}
		
	}
	return totalEnergy;
}

double calcLJ(Atom atom1, Atom atom2, double r2)
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

double calcCharge(double charge1, double charge2, double r)
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

double makePeriodic(double x, double box)
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

double calcBlending(double d1, double d2)
{
	return sqrt(d1 * d2);
}

int getXFromIndex(int index)
{
	int c = -2 * index;
    int discriminant = 1 - 4 * c;
    int qv = (-1 + sqrtf(discriminant)) / 2;
    return qv + 1;
}

int getYFromIndex(int x, int index)
{
	return index - (x * x - x) / 2;
}