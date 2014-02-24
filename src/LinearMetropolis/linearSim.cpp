/*!\file
  \Class for simulation Box, including Enviroments and points to molecules, only save all states
  \author David(Xiao Zhang) and Orlando Acevedo.
 
  This file contains implement of SimBox that are used to handle enviroments and common function
  for box.
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <time.h>
#include "linearSim.h"

/**
	Initializes data to be used in the simulation
	@param initbox - Initial box to be used in the simulation
	@param initsteps - Steps to be taken in the simulation - determines how long the simulation will run 
*/
LinearSim::LinearSim(SimBox *initbox,int initsteps)
{
	box=initbox;
	steps=initsteps;
	currentEnergy=0;
	oldEnergy=0;
	accepted=0;
	rejected=0;
}

/**

	Calculates the energy between each atom pair.
	If the atoms are in the same molecule it applies a fudge factor
	@param enviro - pointer to the environment
	@param molecules - pointer to the molecules in the simulation
	@return - returns the new sum of energy in the system
	
*/
double LinearSim::calcSystemEnergy(Molecule *molecules, Environment *enviro)
{
    double totalEnergy = 0;

	//for each molecule
	for (int mol = 0; mol < enviro->numOfMolecules; mol++)
	{
		totalEnergy += calcMolecularEnergyContribution(molecules, enviro, mol, mol);
	}
	
    return totalEnergy;
}

/**
	Calculates the energy contributions of a single molecule to the system.
*/
double LinearSim::calcMolecularEnergyContribution(Molecule *molecules, Environment *enviro, int currentMol, int startIdx = 0)
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

/**
	Calculates energy between two molecules, atom by atom.
*/
double LinearSim::calcInterMolecularEnergy(Molecule *molecules, int mol1, int mol2, Environment *enviro)
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
				double fvalue = 1.0;
				if(mol1 == mol2)
				{
					int ** hopTab1 = box->tables[mol1 % box->molecTypenum].hopTable;
					fValue = box->getFValue(i, j, hopTab1);
				}
				
				totalEnergy += calc_lj(atom1, atom2, r2) * fvalue;
				totalEnergy += calcCharge(atom1.charge, atom2.charge, sqrt(r2)) * fValue;
			}
		}
		
	}
	return totalEnergy;
}

/**
	Calculates the Lennard-Jones energies between atoms
	@param atom1 - first atom
	@param atom2 - second atom
	@param enviro - the environment in the simulation
*/
double LinearSim::calc_lj(Atom atom1, Atom atom2, double r2)
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

/**
	Calculates the charge relative from one atom to another
*/
double LinearSim::calcCharge(double charge1, double charge2, double r)
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

/**
	Calculates the nonbonded energy for intermolecular molecule pairs using a linked-cell
	neighbor list. The function then calls a separate function to the calculate the
	intramolecular nonbonded interactions for every molecule and sums it to the total
	energy.
*/
double LinearSim::calcEnergy_NLC(Molecule *molecules, Environment *enviro)
{
	// Variables for linked-cell neighbor list	
	int lc[3];            /* Number of cells in the x|y|z direction */
	double rc[3];         /* Length of a cell in the x|y|z direction */
	int head[NCLMAX];     /* Headers for the linked cell lists */
	int mc[3];			  /* Vector cell */
	int lscl[NMAX];       /* Linked cell lists */
	int mc1[3];			  /* Neighbor cells */
	double rshift[3];	  /* Shift coordinates for periodicity */
	const double Region[3] = {enviro->x, enviro->y, enviro->z};  /* MD box lengths */
	int c1;				  /* Used for scalar cell index */
	double dr[3];		  /* Pair vector dr = atom[i]-atom[j] */
	double rrCut = enviro->cutoff * enviro->cutoff;	/* Cutoff squared */
	double rr;			  /* Distance between atoms */
	double lj_energy;		/* Holds current Lennard-Jones energy */
	double charge_energy;	/* Holds current coulombic charge energy */
	double fValue = 1.0;		/* Holds 1,4-fudge factor value */
	double totalEnergy = 0.0;	/* Total nonbonded energy x fudge factor */
			
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
											Atom xAtom, yAtom;
											
				        					for (int atomIn1_i = 0; atomIn1_i < molecules[i].numOfAtoms; atomIn1_i++)
				        					{		
												xAtom = molecules[i].atoms[atomIn1_i];
								
												for (int atomIn2_i = 0; atomIn2_i < molecules[j].numOfAtoms; atomIn2_i++)
												{
													yAtom = molecules[j].atoms[atomIn2_i];
									
													if (xAtom.sigma < 0 || xAtom.epsilon < 0 || yAtom.sigma < 0 || yAtom.epsilon < 0)
													{
														continue;
													}
										
													if(xAtom.id > yAtom.id)
													{
														continue;
													}
				        												
													//store LJ constants locally and define terms in kcal/mol
													const double e = 332.06;
													double sigma = calcBlending(xAtom.sigma, yAtom.sigma);
													double epsilon = calcBlending(xAtom.epsilon, yAtom.epsilon);
									
													//calculate difference in coordinates
													double deltaX = xAtom.x - yAtom.x;
													double deltaY = xAtom.y - yAtom.y;
													double deltaZ = xAtom.z - yAtom.z;
												
													//calculate distance between atoms
													deltaX = box->makePeriodic(deltaX, enviro->x);
													deltaY = box->makePeriodic(deltaY, enviro->y);
													deltaZ = box->makePeriodic(deltaZ, enviro->z);
										
													double r2 = (deltaX * deltaX) +
											 		 	 		(deltaY * deltaY) + 
											 		 			(deltaZ * deltaZ);
														  
													if (r2 == 0.0)
													{
														continue;
													}							

													//calculate LJ energies
													double sig2OverR2 = (sigma * sigma) / r2;
													double sig6OverR6 = sig2OverR2 * sig2OverR2 * sig2OverR2;
													double sig12OverR12 = sig6OverR6 * sig6OverR6;
													lj_energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);

													//calculate Coulombic energies
													double r = sqrt(r2);
													charge_energy = (xAtom.charge * yAtom.charge * e) / r;
										
													double subtotal = (lj_energy + charge_energy) * fValue;
													totalEnergy += subtotal;							
												} /* Endfor atomIn2_i */
											} /* Endfor atomIn_1 */
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
	return totalEnergy + calcIntramolEnergy_NLC(enviro, molecules);
}

/**
	Calculates the nonbonded energy for intramolecular nonbonded interactions for every 
	solvent molecule and sums it to the total energy. Uses getFValue().
*/
double LinearSim::calcIntramolEnergy_NLC(Environment *enviro, Molecule *molecules)
{
    //setup storage
    double totalEnergy = 0.0;
    double *energySum_device;
    // Molecule to be computed. Currently code only handles single solvent type systems.
    // will need to update to handle more than one solvent type (i.e., co-solvents)
	int mol1_i = 0;

    //determine number of energy calculations
    int N =(int) ( pow( (float) molecules[mol1_i].numOfAtoms,2)-molecules[mol1_i].numOfAtoms)/2;	 
    size_t energySumSize = N * sizeof(double);
	double* energySum = (double*) malloc(energySumSize);

    //calulate all energies
    double lj_energy,charge_energy, fValue, nonbonded_energy;
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

				//store LJ constants locally and define terms in kcal/mol
				const double e = 332.06;
				double sigma = calcBlending(atom1.sigma, atom2.sigma);
				double epsilon = calcBlending(atom1.epsilon, atom2.epsilon);
					
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
										  
				if (r2 == 0.0)
				{
					continue;
				}
					
				//calculate LJ energies
				double sig2OverR2 = (sigma * sigma) / r2;
				double sig6OverR6 = sig2OverR2 * sig2OverR2 * sig2OverR2;
				double sig12OverR12 = sig6OverR6 * sig6OverR6;
				double lj_energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);
						
				//calculate Coulombic energies
				double r = sqrt(r2);
				double charge_energy = (atom1.charge * atom2.charge * e) / r;
						
				//gets the fValue in the same molecule
				fValue = 1.0;
				int ** hopTab1 = box->tables[mol1_i % box->molecTypenum].hopTable;
				fValue = box->getFValue(atomIn1_i,atomIn2_i,hopTab1);
						
				double subtotal = (lj_energy + charge_energy) * fValue;
				totalEnergy += subtotal;

		} /* EndFor atomIn2_i */
	} /* EndFor atomIn1_i */
	
	// Multiply single solvent molecule energy by number of solvent molecules in the system
	totalEnergy *= enviro->numOfMolecules;
	
    free(energySum);
    return totalEnergy;
}

double LinearSim::calcBlending(double d1, double d2)
{
    return sqrt(d1 * d2);
}

double LinearSim::Energy_LRC(Molecule *molec, Environment *enviro)
{
	double Ecut = 0.0;		/* Holds LJ long-range cutoff energy correction */
	double NMOL = enviro->numOfMolecules;	/* Number of molecules in simulation */
	// Need to update to accept more than 1 molecule from z-matrix
	int a = 0;
	int NATM = molec[a].numOfAtoms;	/* Number of atoms in molecule */
	// Volume below will have to accept new values in NPT simulations
	double Vnew = enviro->x * enviro->y * enviro->z;	/* Volume of box in Ang^3 */
	
	// Setup arrays to hold all sigma and epsilon atom values
	double Sigma[NATM], Epsilon[NATM];
	for(int i = 0; i < NATM; i++)
	{
		if (molec[a].atoms[i].sigma < 0 || molec[a].atoms[i].epsilon < 0)
		{
			Sigma[i] = 0.0;
			Epsilon[i] = 0.0;
		}
		else
		{
			Sigma[i] = molec[a].atoms[i].sigma;
			Epsilon[i] = molec[a].atoms[i].epsilon;
		}
	}
		
	// (4 * epsilon * sigma^6 or 12)^0.5
	double A6[NATM], A12[NATM];
	double sig2, sig6, sig12;
	for(int i = 0; i < NATM; i++)
	{
		sig2 = pow(Sigma[i], 2);
        sig6 = pow(sig2, 3);
    	sig12 = pow(sig6, 2);
		A6[i] = sqrt(4 * Epsilon[i] * sig6);
		A12[i] = sqrt(4 * Epsilon[i] * sig12);
	}	
	
	double RC3 = 1.00 / pow(enviro->cutoff, 3);		/* 1 / cutoff^3 */
	double RC9 = pow(RC3, 3);		/* 1 / cutoff^9 */
	// Loop over all atoms in a pair
	for(int i = 0; i < NATM; i++)
	{
		for(int j = 0; j < NATM; j++)
		{
			Ecut += (2*PI*NMOL*NMOL/(3.0*Vnew)) * (A12[i]*A12[j]*RC9/3.0 - A6[i]*A6[j]*RC3);
		}
	}
	
	return Ecut;
}

/**
	Starts the linear simulation.
*/
void LinearSim::runLinear(int steps)
{
    Molecule *molecules=box->getMolecules();
 	Environment *enviro=box->getEnviro();
 	  
    int numberOfAtoms = enviro->numOfAtoms;
    double maxTranslation = enviro->maxTranslation;
    double maxRotation = enviro->maxRotation;
    double temperature = enviro->temperature;
    double kT = kBoltz * temperature;

    int atomTotal = 0;
    int aIndex = 0;
    int mIndex = 0;
    double newEnergyCont, oldEnergyCont;

    if (oldEnergy == 0)
	{
		oldEnergy = calcSystemEnergy(molecules, enviro);
	}
	 
    for(int move = 0; move < steps; move++)
    {
            
        int changeIdx = box->chooseMolecule();

	oldEnergyCont = calcMolecularEnergyContribution(molecules, enviro, changeIdx);
	
	box->changeMolecule(changeIdx);
	
	newEnergyCont = calcMolecularEnergyContribution(molecules, enviro, changeIdx);

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
            box->Rollback(changeIdx);
        }
    }
    currentEnergy=oldEnergy;
}
