/*
	Contains the methods required to calculate energies serially.

	Created: February 21, 2014
	
	-> February 26, by Albert Wallace
	-> March 28, by Joshua Mosby
	-> April 21, by Nathan Coleman
	-> February 25, 2015 by Jared Brown
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

Real SerialCalcs::calcSystemEnergy(Molecule *molecules, Environment *enviro, Real &subLJ, Real &subCharge)
{
    Real totalEnergy = 0;

	//for each molecule
	for (int mol = 0; mol < enviro->numOfMolecules; mol++)
	{
		totalEnergy += calcMolecularEnergyContribution(molecules, enviro, subLJ, subCharge, mol, mol);
	}
	
    return totalEnergy;
}

Real SerialCalcs::calcMolecularEnergyContribution(Molecule *molecules, Environment *enviro, Real &subLJ, Real &subCharge, int currentMol, int startIdx)
{
    Real totalEnergy = 0;

	//for every other molecule
	#pragma omp parallel for //num_threads(4) <- this is set in Simulation.cpp
	for (int otherMol = startIdx; otherMol < enviro->numOfMolecules; otherMol++)
	{
	    bool included = false;

		if (otherMol != currentMol)
		{
			//grab the primary indexes for the current and other molecules	
			std::vector<int> currentMolPrimaryIndexArray = (*(*(enviro->primaryAtomIndexArray))[molecules[currentMol].type]);
			std::vector<int> otherMolPrimaryIndexArray;
			if (molecules[currentMol].type == molecules[otherMol].type)
			{
				otherMolPrimaryIndexArray = currentMolPrimaryIndexArray;
			}
			else 
			{
				otherMolPrimaryIndexArray = (*(*(enviro->primaryAtomIndexArray))[molecules[otherMol].type]);
			}

			//If any of otherMol's primary index atoms are within any of currentMol's
			//cutoffs, the inter-molecular energy between the two molcules will be calculated
			for (int i = 0; i < currentMolPrimaryIndexArray.size(); i++)
			{
			    for (int j = 0; j < otherMolPrimaryIndexArray.size(); j++)
    			{
					int primaryIndex1 = currentMolPrimaryIndexArray[i];
					int primaryIndex2 = otherMolPrimaryIndexArray[j];
					Atom atom1 = molecules[currentMol].atoms[primaryIndex1];
					Atom atom2 = molecules[otherMol].atoms[primaryIndex2];
				
					//square cutoff value for easy comparison
					Real cutoffSQ = enviro->cutoff * enviro->cutoff;
			
					//calculate squared distance between atoms 
					Real r2 = calcAtomDist(atom1, atom2, enviro);

					if (r2 < cutoffSQ)
					{
						Real lj_energy = 0, charge_energy = 0;
						Real tempEnergy = calcInterMolecularEnergy(molecules, currentMol, otherMol, enviro, lj_energy, charge_energy);
						
						if (tempEnergy > 1000000)
							cout << "currentMol: " << currentMol << ", otherMol: " << otherMol << ", tempEnergy: " << tempEnergy << endl;	
						//this addition needs to be atomic since multiple threads will be modifying totalEnergy
						#pragma omp critical
						{
							totalEnergy += tempEnergy;
							subLJ += lj_energy;
							subCharge += charge_energy;
						}
						
						//Molecule has been included. Skipping rest of primary indexes for otherMol
						included = true;
						break;
					}
			    }
			    if (included)
			    {
					break;
			    }
			}
		}
	}

	return totalEnergy;
}


/* ------ Neighbor-List System Energy Calculation Functions ------ */

Real SerialCalcs::calcSystemEnergy_NLC(NeighborList *nl, Molecule *molecules, Environment *enviro, Real &subLJ, Real &subCharge)
{
	Real totalEnergy = 0.0;			/* Total system energy */
	
  /* Calculate pair interaction-----------------------------------------------*/
	int vectorCells[3];			  	/* Vector cells */
	
	// Scan inner cells
	for (vectorCells[0] = 0; vectorCells[0] < nl->numCells[0]; (vectorCells[0])++)
		for (vectorCells[1] = 0; vectorCells[1] < nl->numCells[1]; (vectorCells[1])++)
			for (vectorCells[2] = 0; vectorCells[2] < nl->numCells[2]; (vectorCells[2])++)
			{

				// Calculate a scalar cell index
				int c = vectorCells[0] * nl->numCellsYZ 
					+ vectorCells[1] * nl->numCells[2] 
					+ vectorCells[2];
				// Skip this cell if empty
				if (nl->head[c] == EMPTY) continue;
				
				// For each molecule i in cell c
				int i = nl->head[c];
				while (i != EMPTY)
				{
					// Add the molecular energy contribution to the system total
					std::vector<int> neighbors;
					getNeighbors_NLC(nl, molecules, enviro, i, neighbors, true);
					totalEnergy += calcMolecularEnergyContribution_NLC(molecules, enviro, subLJ, subCharge, i, neighbors);
					
					i = nl->linkedCellList[i];
				} /* Endwhile i not empty */
			} /* Endfor central cell, c */
	
	return totalEnergy;
}

void SerialCalcs::getNeighbors_NLC(NeighborList *nl, Molecule *molecules, Environment *enviro, int currentMol, std::vector<int>& neighbors, bool isSysCalc) {
	// Find the vector cell for the currentMol (based on 1st primary index)
	// NOTE: for multiple solvents each with multiple primary indexes, this is slightly off for
	//	special cases when a molecule has primary indexes across cell boundaries
	std::vector<int> currentMolPrimaryIndexArray = (*(*(enviro->primaryAtomIndexArray))[molecules[currentMol].type]);
	int primaryIndex = currentMolPrimaryIndexArray[0]; // Use first primary index to determine cell placement
		
	int vectorCells[3];			
	vectorCells[0] = molecules[currentMol].atoms[primaryIndex].x / nl->lengthCell[0]; 
	vectorCells[1] = molecules[currentMol].atoms[primaryIndex].y / nl->lengthCell[1];
	vectorCells[2] = molecules[currentMol].atoms[primaryIndex].z / nl->lengthCell[2];
	
	// Scan the neighbor cells (including itself) of cell c
	int neighborCells[3];			/* Neighbor cells */
	for (neighborCells[0] = vectorCells[0]-1; neighborCells[0] <= vectorCells[0]+1; (neighborCells[0])++)
		for (neighborCells[1] = vectorCells[1]-1; neighborCells[1] <= vectorCells[1]+1; (neighborCells[1])++)
			for (neighborCells[2] = vectorCells[2]-1; neighborCells[2] <= vectorCells[2]+1; (neighborCells[2])++)
			{
				// Periodic boundary condition by shifting coordinates
				Real rshift[3];	  		/* Shift coordinates for periodicity */
				for (int a = 0; a < 3; a++)
				{
					if (neighborCells[a] < 0)
					{
						rshift[a] = -nl->region[a];
					}
					else if (neighborCells[a] >= nl->numCells[a])
					{
						rshift[a] = nl->region[a];
					}
					else
					{
						rshift[a] = 0.0;
					}
				}
		
				// Calculate the scalar cell index of the neighbor cell
				int c1 = ((neighborCells[0] + nl->numCells[0]) % nl->numCells[0]) * nl->numCellsYZ
					+((neighborCells[1] + nl->numCells[1]) % nl->numCells[1]) * nl->numCells[2]
					+((neighborCells[2] + nl->numCells[2]) % nl->numCells[2]);
				
				// Skip this neighbor cell if empty
				if (nl->head[c1] == EMPTY) continue;
				
				int otherMol = nl->head[c1];
				
				while (otherMol != EMPTY)
				{
					bool included = false;

					// For other molecules (and if total system calc, then avoid double counting of pairs)
					if (((currentMol != otherMol) && !isSysCalc) || (currentMol < otherMol))
					{	
						std::vector<int> otherMolPrimaryIndexArray;
						if (molecules[currentMol].type == molecules[otherMol].type)
						{
							otherMolPrimaryIndexArray = currentMolPrimaryIndexArray;
						}
						else 
						{
							otherMolPrimaryIndexArray = (*(*(enviro->primaryAtomIndexArray))[molecules[otherMol].type]);
						}

						for (int i1 = 0; i1 < currentMolPrimaryIndexArray.size(); i1++)
						{
							for (int i2 = 0; i2 < otherMolPrimaryIndexArray.size(); i2++)
							{
								int primaryIndex1 = currentMolPrimaryIndexArray[i1];
								int primaryIndex2 = otherMolPrimaryIndexArray[i2];
								Atom atom1 = molecules[currentMol].atoms[primaryIndex1];
								Atom atom2 = molecules[otherMol].atoms[primaryIndex2];
							
								Real dr[3];		  /* Pair vector dr = atom[currentMol]-atom[otherMol] */
								dr[0] = atom1.x - (atom2.x + rshift[0]);
								dr[1] = atom1.y - (atom2.y + rshift[1]);
								dr[2] = atom1.z - (atom2.z + rshift[2]);
								Real rr = (dr[0] * dr[0]) + (dr[1] * dr[1]) + (dr[2] * dr[2]);			
								
								// Calculate energy for entire molecule interaction if rij < Cutoff for atom index
								if (rr < nl->rrCut)
								{	
									neighbors.push_back(otherMol);
									
									included = true;
									break;
								} /* Endif rr < rrCut */
									
							}
							if (included)
							{
								break;
							}
						}
					} /* Endif i<j */
			
					otherMol = nl->linkedCellList[otherMol];
				} /* Endwhile otherMol not empty */
			} /* Endfor neighbor cells, c1 */
}

Real SerialCalcs::calcMolecularEnergyContribution_NLC(Molecule *molecules, Environment *enviro, Real &subLJ, Real &subCharge, int currentMol, std::vector<int> neighbors) {
	Real totalEnergy = 0.0;			/* Total nonbonded energy x fudge factor */
	Real fValue = 1.0;				/* Holds 1,4-fudge factor value */
	
	#pragma omp parallel for
	for (int index = 0; index < neighbors.size(); index++)
	{
		Real lj_energy = 0, charge_energy = 0;
		Real tempEnergy = calcInterMolecularEnergy(molecules, currentMol, neighbors[index], enviro, lj_energy, charge_energy) * fValue;
		
		if (tempEnergy > 1000000)
		{
			cout << "OtherMol: " << index << ", tempEnergy: " << tempEnergy << endl;
		}							
		#pragma omp critical
		{
			totalEnergy += tempEnergy;
			subLJ += lj_energy;
			subCharge += charge_energy;
		}
	}
	
	return totalEnergy;
}


/* ------ Utility Calculation Functions ------ */

Real SerialCalcs::calcInterMolecularEnergy(Molecule *molecules, int mol1, int mol2, Environment *enviro, Real &subLJ, Real &subCharge)
{
	Real totalEnergy = 0, lj_energy = 0, charge_energy = 0;
	
	for (int i = 0; i < molecules[mol1].numOfAtoms; i++)
	{
		Atom atom1 = molecules[mol1].atoms[i];
		
		for (int j = 0; j < molecules[mol2].numOfAtoms; j++)
		{
			Atom atom2 = molecules[mol2].atoms[j];
		
			if (atom1.sigma >= 0 && atom1.epsilon >= 0 && atom2.sigma >= 0 && atom2.epsilon >= 0)
			{
				//calculate squared distance between atoms 
				Real r2 = calcAtomDist(atom1, atom2, enviro);
				
				lj_energy = calc_lj(atom1, atom2, r2);
				subLJ += lj_energy;
				
				charge_energy = calcCharge(atom1.charge, atom2.charge, sqrt(r2));
				subCharge += charge_energy;
				
				totalEnergy += lj_energy + charge_energy;
			}
		}
		
	}
	
	return totalEnergy;
}

Real SerialCalcs::calcIntraMolecularEnergy(Molecule *molecules, Environment *enviro, Real &subLJ, Real &subCharge)
{
    Real totalEnergy = 0, lj_energy = 0, charge_energy = 0, nonbonded_energy = 0;
	
    // NOTE: Current code only handles one or two solvent type systems. ****
	int NMOL = 2;
	
	for (int molIndex = 0; molIndex < NMOL - 1; molIndex++)
	{
		//calculate all interatomic energies
		Atom atom1, atom2;
		for (int atomIn1_i = 0; atomIn1_i < molecules[molIndex].numOfAtoms; atomIn1_i++)
		{	
			atom1 = molecules[molIndex].atoms[atomIn1_i];
						
			for (int atomIn2_i = atomIn1_i; atomIn2_i < molecules[molIndex].numOfAtoms; atomIn2_i++)
			{
				atom2 = molecules[molIndex].atoms[atomIn2_i];
					
				if (atom1.sigma < 0 || atom1.epsilon < 0 || atom2.sigma < 0 || atom2.epsilon < 0) continue;
					
				//calculate squared distance between atoms 
				Real r2 = calcAtomDist(atom1, atom2, enviro);
					
				if (r2 == 0.0) 
					continue;
						
				//calculate LJ energies
				lj_energy = calc_lj(atom1, atom2, r2);
				subLJ += lj_energy;
					
				//calculate Coulombic energies
				charge_energy = calcCharge(atom1.charge, atom2.charge, sqrt(r2));
				subCharge += charge_energy;
					
				//gets the fValue in the same molecule
				Real fValue = 0.0;
					
				int hops = 0;
				for (int k = 0; k < molecules[molIndex].numOfHops; k++)
				{
					Hop currentHop = molecules[molIndex].hops[k];
					if (currentHop.atom1 == atomIn1_i && currentHop.atom2 == atomIn2_i)
					{
						hops = currentHop.hop;
					}
				}
					
				if (hops == 3)
					fValue = 0.5;
				else if (hops > 3)
					fValue = 1.0;
							
				Real subtotal = (lj_energy + charge_energy) * fValue;
				totalEnergy += subtotal;
			} /* EndFor atomIn2_i */
		} /* EndFor atomIn1_i */
	}

	totalEnergy *= (enviro->numOfMolecules / NMOL);

    return totalEnergy;
}

Real SerialCalcs::calcEnergy_LRC(Molecule *molecules, Environment *enviro)
{
	Real Ecut = 0.0;		// Holds LJ long-range cutoff energy correction 
	
	Real Vnew = enviro->x * enviro->y * enviro->z;	// Volume of box in Ang^3
	Real RC3 = 1.00 / pow(enviro->cutoff, 3);		// 1 / cutoff^3
	Real RC9 = pow(RC3, 3);							// 1 / cutoff^9
	
	// Note: currently only supports at most TWO solvents (needs to be updated for more)
	int a = 0, b = 1;
	Real NMOL1 = enviro->numOfMolecules / 2;	// Number of molecules of solvent1
	Real NMOL2 = enviro->numOfMolecules / 2;	// Number of molecules of solvent2
	int NATOM1 = molecules[a].numOfAtoms;			// Number of atoms in solvent1
	int NATOM2 = molecules[b].numOfAtoms;			// Number of atoms in solvent2
	int NATMX = NATOM1;
	if (NATMX < NATOM2)		// NATMX = MAX(NAT0M1, NAT0M2)
	{
		NATMX = NATOM2;
	}
	
	Real sig2, sig6, sig12;
	// get LJ-values for solvent1 and store in A6, A12
	Real SigmaA[NATOM1], EpsilonA[NATOM1];
	Real A6[NATOM1], A12[NATOM1];
	for(int i = 0; i < NATOM1; i++)
	{
		if (molecules[a].atoms[i].sigma < 0 || molecules[a].atoms[i].epsilon < 0)
		{
			SigmaA[i] = 0.0;
			EpsilonA[i] = 0.0;
		}
		else
		{
			SigmaA[i] = molecules[a].atoms[i].sigma;
			EpsilonA[i] = molecules[a].atoms[i].epsilon;
		}
		
		sig2 = pow(SigmaA[i], 2);
        sig6 = pow(sig2, 3);
    	sig12 = pow(sig6, 2);
		A6[i] = sqrt(4 * EpsilonA[i] * sig6);
		A12[i] = sqrt(4 * EpsilonA[i] * sig12);
	}
	
	// get LJ-values for solvent2 and store in B6, B12
	Real SigmaB[NATOM2], EpsilonB[NATOM2];
	Real B6[NATOM2], B12[NATOM2];
	for(int i = 0; i < NATOM2; i++)
	{
		if (molecules[b].atoms[i].sigma < 0 || molecules[b].atoms[i].epsilon < 0)
		{
			SigmaB[i] = 0.0;
			EpsilonB[i] = 0.0;
		}
		else
		{
			SigmaB[i] = molecules[b].atoms[i].sigma;
			EpsilonB[i] = molecules[b].atoms[i].epsilon;
		}
		
		sig2 = pow(SigmaB[i], 2);
        sig6 = pow(sig2, 3);
    	sig12 = pow(sig6, 2);
		B6[i] = sqrt(4 * EpsilonB[i] * sig6);
		B12[i] = sqrt(4 * EpsilonB[i] * sig12);
	}
	
	// loop over all atoms in a pair
	for(int i = 0; i < NATOM1; i++)
	{
		for(int j = 0; j < NATOM2; j++)
		{
			Ecut += (2*PI*NMOL1*NMOL1/(3.0*Vnew)) * (A12[i]*A12[j]*RC9/3.0 - A6[i]*A6[j]*RC3);
			Ecut += (2*PI*NMOL2*NMOL2/(3.0*Vnew)) * (B12[i]*B12[j]*RC9/3.0 - B6[i]*B6[j]*RC3);
			Ecut += (4*PI*NMOL1*NMOL2/(3.0*Vnew)) * (A12[i]*B12[j]*RC9/3.0 - A6[i]*B6[j]*RC3);
		}
	}

	//std::cout << "Energy_LRC = " << Ecut << std::endl;
	return Ecut;
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
    
    if(x < -0.5 * boxDim)
    {
        x += boxDim;
    }

    else if(x > 0.5 * boxDim)
    {
        x -= boxDim;
    }

    return x;

} 

Real SerialCalcs::calcBlending(Real d1, Real d2)
{
    return sqrt(d1 * d2);
}

Real SerialCalcs::calcAtomDist(Atom atom1, Atom atom2, Environment *enviro)
{
	//calculate difference in coordinates
	Real deltaX = makePeriodic(atom1.x - atom2.x, enviro->x);
	Real deltaY = makePeriodic(atom1.y - atom2.y, enviro->y);
	Real deltaZ = makePeriodic(atom1.z - atom2.z, enviro->z);
				
	//calculate squared distance (r2 value) and return
	return (deltaX * deltaX) + (deltaY * deltaY) + (deltaZ * deltaZ);
}
