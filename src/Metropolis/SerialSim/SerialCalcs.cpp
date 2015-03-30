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

Real SerialCalcs::calcSystemEnergy(Molecule *molecules, Environment *enviro)
{
    Real totalEnergy = 0;

	//for each molecule
	for (int mol = 0; mol < enviro->numOfMolecules; mol++)
	{
		totalEnergy += calcMolecularEnergyContribution(molecules, enviro, mol, mol);
	}
	
    return totalEnergy + Energy_LRC(molecules, enviro);
}

Real SerialCalcs::calcMolecularEnergyContribution(Molecule *molecules, Environment *environment, int currentMol, int startIdx)
{
    Real totalEnergy = 0;

	//for every other molecule
	#pragma omp parallel for //num_threads(4) <- this is set in Simulation.cpp
	for (int otherMol = startIdx; otherMol < environment->numOfMolecules; otherMol++)
	{
	    bool included = false;

		if (otherMol != currentMol)
		{
			Atom atom1 = molecules[currentMol].atoms[environment->primaryAtomIndex];
			Atom atom2 = molecules[otherMol].atoms[environment->primaryAtomIndex];
			
			std::vector<int> currentMolPrimaryIndexArray = (*(*(environment->primaryAtomIndexArray))[molecules[currentMol].type]);
			std::vector<int> otherMolPrimaryIndexArray;
			if (molecules[currentMol].type == molecules[otherMol].type)
			{
				otherMolPrimaryIndexArray = currentMolPrimaryIndexArray;
			}
			else 
			{
				otherMolPrimaryIndexArray = (*(*(environment->primaryAtomIndexArray))[molecules[otherMol].type]);
			}

			for (int i = 0; i < currentMolPrimaryIndexArray.size(); i++)
			{
			    for (int j = 0; j < otherMolPrimaryIndexArray.size(); j++)
    			{
					int primaryIndex1 = currentMolPrimaryIndexArray[i];
					int primaryIndex2 = otherMolPrimaryIndexArray[j];
					Atom atom1 = molecules[currentMol].atoms[primaryIndex1];
					Atom atom2 = molecules[otherMol].atoms[primaryIndex2];
				
					//square cutoff value for easy comparison
					Real cutoffSQ = environment->cutoff * environment->cutoff;
					
			
					//calculate squared distance between atoms 
					Real r2 = calcAtomDist(atom1, atom2, environment);

					if (r2 < cutoffSQ)
					{
						Real tempEnergy = calcInterMolecularEnergy(molecules, currentMol, otherMol, environment);
						//this addition needs to be atomic since multiple threads will be modifying totalEnergy
						#pragma omp atomic
						totalEnergy += tempEnergy;
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

/**
	Calculates the nonbonded energy for intermolecular molecule pairs using a linked-cell
	neighbor list. The function then calls a separate function to the calculate the
	intramolecular nonbonded interactions for every molecule and sums it to the total
	energy.
*/
Real SerialCalcs::calcEnergy_NLC(Molecule *molecules, Environment *enviro)
{
	// Variables for linked-cell neighbor list	
	int lc[3];            	/* Number of cells in the x|y|z direction */
	Real rc[3];         	/* Length of a cell in the x|y|z direction */
	int head[NCLMAX];    	/* Headers for the linked cell lists */
	int mc[3];			  	/* Vector cell */
	int lscl[NMAX];       	/* Linked cell lists */
	int mc1[3];				/* Neighbor cells */
	Real rshift[3];	  		/* Shift coordinates for periodicity */
	const Real Region[3] = {enviro->x, enviro->y, enviro->z};  /* MD box lengths */
	int c1;				  	/* Used for scalar cell index */
	Real dr[3];		  		/* Pair vector dr = atom[i]-atom[j] */
	Real rrCut = enviro->cutoff * enviro->cutoff;	/* Cutoff squared */
	Real rr;			  	/* Distance between atoms */
	Real lj_energy;			/* Holds current Lennard-Jones energy */
	Real charge_energy;		/* Holds current coulombic charge energy */
	Real fValue = 1.0;		/* Holds 1,4-fudge factor value */
	Real totalEnergy = 0.0;	/* Total nonbonded energy x fudge factor */
			
	// Compute the # of cells for linked cell lists
	for (int k=0; k<3; k++)
	{
		lc[k] = Region[k] / enviro->cutoff; 
		rc[k] = Region[k] / lc[k];
	}
		
  /* Make a linked-cell list, lscl--------------------------------------------*/
	int lcyz = lc[1]*lc[2];
	int lcxyz = lc[0]*lcyz;
		
	// Reset the headers, head
	for (int c = 0; c < lcxyz; c++) 
	{
		head[c] = EMPTY;
	}

	// Scan cutoff index atom in each molecule to construct headers, head, & linked lists, lscl
	for (int i = 0; i < enviro->numOfMolecules; i++)
	{
		mc[0] = molecules[i].atoms[enviro->primaryAtomIndex].x / rc[0]; 
		mc[1] = molecules[i].atoms[enviro->primaryAtomIndex].y / rc[1];
		mc[2] = molecules[i].atoms[enviro->primaryAtomIndex].z / rc[2];
		
		// Translate the vector cell index, mc, to a scalar cell index
		int c = mc[0]*lcyz + mc[1]*lc[2] + mc[2];

		// Link to the previous occupant (or EMPTY if you're the 1st)
		lscl[i] = head[c];

		// The last one goes to the header
		head[c] = i;
	} /* Endfor molecule i */

  /* Calculate pair interaction-----------------------------------------------*/
		
	// Scan inner cells
	for (mc[0] = 0; mc[0] < lc[0]; (mc[0])++)
	{
		for (mc[1] = 0; mc[1] < lc[1]; (mc[1])++)
		{
			for (mc[2] = 0; mc[2] < lc[2]; (mc[2])++)
			{

				// Calculate a scalar cell index
				int c = mc[0]*lcyz + mc[1]*lc[2] + mc[2];
				// Skip this cell if empty
				if (head[c] == EMPTY) continue;

				// Scan the neighbor cells (including itself) of cell c
				for (mc1[0] = mc[0]-1; mc1[0] <= mc[0]+1; (mc1[0])++)
					for (mc1[1] = mc[1]-1; mc1[1] <= mc[1]+1; (mc1[1])++)
						for (mc1[2] = mc[2]-1; mc1[2] <= mc [2]+1; (mc1[2])++)
						{
							// Periodic boundary condition by shifting coordinates
							for (int a = 0; a < 3; a++)
							{
								if (mc1[a] < 0)
								{
									rshift[a] = -Region[a];
								}
								else if (mc1[a] >= lc[a])
								{
									rshift[a] = Region[a];
								}
								else
								{
									rshift[a] = 0.0;
								}
							}
							// Calculate the scalar cell index of the neighbor cell
							c1 = ((mc1[0] + lc[0]) % lc[0]) * lcyz
							    +((mc1[1] + lc[1]) % lc[1]) * lc[2]
							    +((mc1[2] + lc[2]) % lc[2]);
							// Skip this neighbor cell if empty
							if (head[c1] == EMPTY)
							{
								continue;
							}

							// Scan atom i in cell c
							int i = head[c];
							while (i != EMPTY)
							{

								// Scan atom j in cell c1
								int j = head[c1];
								while (j != EMPTY)
								{

									// Avoid double counting of pairs
									if (i < j)
									{
										// Pair vector dr = atom[i]-atom[j]
										rr = 0.0;
										dr[0] = molecules[i].atoms[enviro->primaryAtomIndex].x - (molecules[j].atoms[enviro->primaryAtomIndex].x + rshift[0]);
										dr[1] = molecules[i].atoms[enviro->primaryAtomIndex].y - (molecules[j].atoms[enviro->primaryAtomIndex].y + rshift[1]);
										dr[2] = molecules[i].atoms[enviro->primaryAtomIndex].z - (molecules[j].atoms[enviro->primaryAtomIndex].z + rshift[2]);
										rr = (dr[0] * dr[0]) + (dr[1] * dr[1]) + (dr[2] * dr[2]);			
										
										// Calculate energy for entire molecule interaction if rij < Cutoff for atom index
										if (rr < rrCut)
										{	
											totalEnergy += calcInterMolecularEnergy(molecules, i, j, enviro) * fValue;
											
										} /* Endif rr < rrCut */
									} /* Endif i<j */
									
									j = lscl[j];
								} /* Endwhile j not empty */

								i = lscl[i];
							} /* Endwhile i not empty */
						} /* Endfor neighbor cells, c1 */
			} /* Endfor central cell, c */
		}
	}
	
	// *** Note: this function returns values similar to calcSystemEnergy before the calcIntramolEnergy_NLC is called
	return totalEnergy + calcIntramolEnergy_NLC(enviro, molecules) + Energy_LRC(molecules, enviro);
	//return totalEnergy;
}

/**
	Calculates the nonbonded energy for intramolecular nonbonded interactions for every 
	solvent molecule and sums it to the total energy. Uses getFValue().
*/
Real SerialCalcs::calcIntramolEnergy_NLC(Environment *enviro, Molecule *molecules)
{
    //setup storage
    Real totalEnergy = 0.0;
    Real *energySum_device;
    // Molecule to be computed. Currently code only handles single solvent type systems.
    // will need to update to handle more than one solvent type (i.e., co-solvents)
	int mol1_i = 0;

    //determine number of energy calculations
    int N =(int) ( pow( (float) molecules[mol1_i].numOfAtoms,2)-molecules[mol1_i].numOfAtoms)/2;	 
    size_t energySumSize = N * sizeof(Real);
	Real* energySum = (Real*) malloc(energySumSize);

    //calculate all energies
    Real lj_energy, charge_energy, fValue, nonbonded_energy;
    Atom atom1, atom2;

	for (int atomIn1_i = 0; atomIn1_i < molecules[mol1_i].numOfAtoms; atomIn1_i++)
	{	
		atom1 = molecules[mol1_i].atoms[atomIn1_i];
					
		for (int atomIn2_i = atomIn1_i; atomIn2_i < molecules[mol1_i].numOfAtoms; atomIn2_i++)
		{
			atom2 = molecules[mol1_i].atoms[atomIn2_i];
						
				if (atom1.sigma < 0 || atom1.epsilon < 0 || atom2.sigma < 0 || atom2.epsilon < 0)
				{
					continue;
				}
					  
				//calculate squared distance between atoms 
				Real r2 = calcAtomDist(atom1, atom2, enviro);
										  
				if (r2 == 0.0)
				{
					continue;
				}
					
				//calculate LJ energies
				lj_energy = calc_lj(atom1, atom2, r2);
						
				//calculate Coulombic energies
				charge_energy = calcCharge(atom1.charge, atom2.charge, sqrt(r2));
						
				//gets the fValue in the same molecule
				fValue = 0.0;
				
				int hops = 0;
				for (int k = 0; k < molecules[mol1_i].numOfHops; k++)
				{
					Hop currentHop = molecules[mol1_i].hops[k];
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
	
	// Multiply single solvent molecule energy by number of solvent molecules in the system
	totalEnergy *= enviro->numOfMolecules;
	
    free(energySum);
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
				//calculate squared distance between atoms 
				Real r2 = calcAtomDist(atom1, atom2, enviro);
				
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

Real SerialCalcs::calcAtomDist(Atom atom1, Atom atom2, Environment *enviro)
{
	//calculate difference in coordinates
	Real deltaX = makePeriodic(atom1.x - atom2.x, enviro->x);
	Real deltaY = makePeriodic(atom1.y - atom2.y, enviro->y);
	Real deltaZ = makePeriodic(atom1.z - atom2.z, enviro->z);
				
	//calculate squared distance (r2 value) and return
	return (deltaX * deltaX) + (deltaY * deltaY) + (deltaZ * deltaZ);
}


Real SerialCalcs::Energy_LRC(Molecule *molec, Environment *enviro)
{
	Real Ecut = 0.0;		// Holds LJ long-range cutoff energy correction 
	
	Real Vnew = enviro->x * enviro->y * enviro->z;	// Volume of box in Ang^3
	Real RC3 = 1.00 / pow(enviro->cutoff, 3);		// 1 / cutoff^3
	Real RC9 = pow(RC3, 3);							// 1 / cutoff^9
	
	int a = 0, b = 1;
	Real NMOL1 = enviro->numOfMolecules / 2;	// Number of molecules of solvent1
	Real NMOL2 = enviro->numOfMolecules / 2;	// Number of molecules of solvent1
	int NATOM1 = molec[a].numOfAtoms;		// Number of atoms in solvent1
	int NATOM2 = molec[b].numOfAtoms;		// Number of atoms in solvent2
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
		if (molec[a].atoms[i].sigma < 0 || molec[a].atoms[i].epsilon < 0)
		{
			SigmaA[i] = 0.0;
			EpsilonA[i] = 0.0;
		}
		else
		{
			SigmaA[i] = molec[a].atoms[i].sigma;
			EpsilonA[i] = molec[a].atoms[i].epsilon;
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
		if (molec[b].atoms[i].sigma < 0 || molec[b].atoms[i].epsilon < 0)
		{
			SigmaB[i] = 0.0;
			EpsilonB[i] = 0.0;
		}
		else
		{
			SigmaB[i] = molec[b].atoms[i].sigma;
			EpsilonB[i] = molec[b].atoms[i].epsilon;
		}
		
		sig2 = pow(SigmaB[i], 2);
        sig6 = pow(sig2, 3);
    	sig12 = pow(sig6, 2);
		B6[i] = sqrt(4 * EpsilonB[i] * sig6);
		B12[i] = sqrt(4 * EpsilonB[i] * sig12);
	}
	
	// loop over all atoms in a pair
	for(int i = 0; i < NATMX; i++)
	{
		for(int j = 0; j < NATMX; j++)
		{
			Ecut += (2*PI*NMOL1*NMOL1/(3.0*Vnew)) * (A12[i]*A12[j]*RC9/3.0 - A6[i]*A6[j]*RC3);
			Ecut += (2*PI*NMOL2*NMOL2/(3.0*Vnew)) * (B12[i]*B12[j]*RC9/3.0 - B6[i]*B6[j]*RC3);
			Ecut += (4*PI*NMOL1*NMOL2/(3.0*Vnew)) * (A12[i]*B12[j]*RC9/3.0 - A6[i]*B6[j]*RC3);
		}
	}

	std::cout << "Energy_LRC = " << Ecut << std::endl;
	return Ecut;
}
