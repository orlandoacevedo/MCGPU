/*!\file
  \Class for parallel Simulation, including Energy calculate and points to molocoles,only save all states
  \author David(Xiao Zhang).
 
  This file contains implement of SimBox that are used to handle enviroments and common function
  for box.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <time.h>
#include "parallelSim.cuh"

__global__ void calcEnergy(Atom *atoms, Environment *enviro, double *energySum, DeviceMolecule *dev_molecules, Hop *hops);

ParallelSim::ParallelSim(GPUSimBox *initbox,int initsteps)
{
	box=initbox;
    steps=initsteps;
 	currentEnergy=0;
    oldEnergy=0;
    accepted=0;
    rejected=0;
  
    Environment *enviro_host=box->getSimBox()->getEnviro();

    //calculate CUDA thread mgmt
    N =(int) ( pow( (float) enviro_host->numOfAtoms,2)-enviro_host->numOfAtoms)/2;
    blocks = N / THREADS_PER_BLOCK + (N % THREADS_PER_BLOCK == 0 ? 0 : 1); 
  
    size_t energySumSize = N * sizeof(double);
    
   //allocate memory on the device
   energySum_host = (double *) malloc(energySumSize);
   cudaMalloc((void **) &energySum_device, energySumSize);
  
}
ParallelSim::~ParallelSim()
{
    if (energySum_host!=NULL)
    {
        free(energySum_host);
        energySum_host=NULL;
    }
  
    if (energySum_device!=NULL)
    {
        cudaFree(energySum_device);
        energySum_device=NULL;
    }
}


double ParallelSim::calcEnergyWrapper(Molecule *molecules, Environment *enviro)
{
    
    Atom *atoms = (Atom *) malloc(sizeof(Atom) * enviro->numOfAtoms);
    int atomIndex = 0;
    double energy=0;
    
    for(int i = 0; i < enviro->numOfMolecules; i++)
    {
        Molecule currentMolecule = molecules[i];
        for(int j = 0; j < currentMolecule.numOfAtoms; j++)
        {
            atoms[atomIndex] = currentMolecule.atoms[j];
            atomIndex++;
        }
    }
    
    energy=calcEnergyWrapper(this->getGPUSimBox());
    free(atoms);
    
    return energy; 
}

//double ParallelSim::calcEnergyWrapper(Atom *atoms, Environment *enviro, Molecule *molecules){
double ParallelSim::calcEnergyWrapper(GPUSimBox *box)
{
    SimBox *innerbox=box->getSimBox();
    Atom *atoms=innerbox->getAtoms();
    Environment *enviro=innerbox->getEnviro();
    Molecule *molecules=innerbox->getMolecules();
    
    
    //setup CUDA storage
    double totalEnergy = 0.0;
    
    //calculate CUDA thread mgmt
    int N =(int) ( pow( (float) enviro->numOfAtoms,2)-enviro->numOfAtoms)/2;
    int blocks = N / THREADS_PER_BLOCK + (N % THREADS_PER_BLOCK == 0 ? 0 : 1); 

    //The number of bytes of shared memory per block of
    //size_t sharedSize = sizeof(double) * THREADS_PER_BLOCK;
    
    //size_t atomSize = enviro->numOfAtoms * sizeof(Atom);
    size_t energySumSize = N * sizeof(double);
    
    //copy data to the device
    box->CopyBoxtoDevice(innerbox);
    Atom *atoms_device=box->getdevAtom();
    Environment *enviro_device=box->getdevEnvironment();
    DeviceMolecule *molec_d=box->getdevDeviceMolecule();
    Hop *hops_d=box->getdevHop();    

    energySum_device =getdevEnergySum();
    //energySum_host=gethostEnergySum();
    energySum_host = (double *) malloc(energySumSize);
    memset(energySum_host,0,energySumSize);

    calcEnergy <<<blocks, THREADS_PER_BLOCK>>>(atoms_device, enviro_device, energySum_device, molec_d, hops_d);
    
    cudaErrorCheck(cudaMemcpy(energySum_host, energySum_device, energySumSize, cudaMemcpyDeviceToHost));

    for(int i = 0; i < N; i++)
    {

        //get atom IDs for each calculation
        int c = -2 * i;
        int discriminant = 1 - 4 * c;
        int qv = (-1 + sqrtf(discriminant)) / 2;
        int atomXid = qv + 1;
        
        int atomYid =  i - (atomXid * atomXid - atomXid) / 2;
        
        //check for stray calculations that returned invalid results
        if (isnan(energySum_host[i]) != 0 || isinf(energySum_host[i]) != 0)
        {
            energySum_host[i] = calcEnergyOnHost(atoms[atomXid], atoms[atomYid], enviro, molecules);
        }
           
        //sum up energies 
        totalEnergy += energySum_host[i];
    }

    return totalEnergy;
}

//calc_lj() //method not needed; functional definition enclosed in calcEnergyOnHost,
// which is the only place the method is used in the corresponding linearSim method, calcEnergy
// !!maybe could also be used in calcEnergy_NLC, as it does have to calc the LJ constant

double ParallelSim::calcEnergyOnHost(Atom atom1, Atom atom2, Environment *enviro, Molecule *molecules)
{
    //define terms in kcal/mol
    const double e = 332.06;
    
    //**************************************************************
    //begin chunk of code that duplicates calc_lj() in linearSim
    //**************************************************************
    double sigma = sqrt(atom1.sigma * atom2.sigma);
    double epsilon = sqrt(atom1.epsilon * atom2.epsilon);
    
    //calculate distance between atoms
    double deltaX = atom1.x - atom2.x;
    double deltaY = atom1.y - atom2.y;
    double deltaZ = atom1.z - atom2.z;
  
    deltaX = box->getSimBox()->makePeriodic(deltaX, enviro->x);
    deltaY = box->getSimBox()->makePeriodic(deltaY, enviro->y);
    deltaZ = box->getSimBox()->makePeriodic(deltaZ, enviro->z);

    double r2 = (deltaX * deltaX) +
                      (deltaY * deltaY) + 
                      (deltaZ * deltaZ);

    //check if atoms overlap
    if (r2 == 0.0)
    {
        //lj_energy = 0.0;
        //charge_energy = 0.0;
        return 0.0;
    }
    else
    {
    	//calculate LJ energies
    	double sig2OverR2 = pow(sigma, 2) / r2;
        double sig6OverR6 = pow(sig2OverR2, 3);
    	double sig12OverR12 = pow(sig6OverR6, 2);
    	double lj_energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);
    //**************************************************************
	//end chunk of code that duplicates calc_lj() in linearSim
	//**************************************************************
	
    	//calculate Coulombic energies
    	double r = sqrt(r2);
    	double charge_energy = (atom2.charge * atom1.charge * e) / r;
    	
    	// Check for 1,3-nonbonded interaction
    	double fValue = 1.0;
		if (molecules != NULL)
        {
        	fValue = getFValueHost(atom1, atom2, molecules, enviro);
        }

    	return fValue * (lj_energy + charge_energy);
	}
}


void ParallelSim::runParallel(int steps)
{
	SimBox   *innerbox=box->getSimBox();
    Molecule *molecules=innerbox->getMolecules();
 	Environment *enviro=innerbox->getEnviro();
 	  
    //int numberOfAtoms = enviro->numOfAtoms;
    //double maxTranslation = enviro->maxTranslation;
    //double maxRotation = enviro->maxRotation;
    double temperature = enviro->temperature;
    double kT = kBoltz * temperature;

    //int atomTotal = 0;
    //int aIndex = 0;
    //int mIndex = 0;
    double newEnergy;
	 
    for(int move = 0; move < steps; move++)
    {
        if (oldEnergy==0)
        {
            oldEnergy = calcEnergyWrapper(molecules, enviro);
        }

        int changeno=innerbox->ChangeMolecule();
		//box->CopyBoxtoDevice(innerbox);
				
        newEnergy = calcEnergyWrapper(this->getGPUSimBox());

        bool accept = false;

        if(newEnergy < oldEnergy)
        {
            accept = true;            
        }
        else
        {
            double x = exp(-(newEnergy - oldEnergy) / kT);
            

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
            oldEnergy=newEnergy;
            
            //box->CopyBoxtoHost(innerbox);
        }
        else{
            rejected++;
            //restore previous configuration
            innerbox->Rollback(changeno);
            oldEnergy=oldEnergy;//meanless, just show we assigned the same energy values.
        }
    }
    currentEnergy=oldEnergy;
}

/**
	Calculates the nonbonded energy for intermolecular molecule pairs using a linked-cell
	neighbor list. The function then calls a separate function to the calculate the
	intramolecular nonbonded interactions for every molecule and sums it to the total
	energy.
*/
/*
double ParallelSim::calcEnergy_NLC(Molecule *molecules, Environment *enviro)
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
		head[c] = EMPTY;

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
		for (mc[1] = 0; mc[1] < lc[1]; (mc[1])++)
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
									rshift[a] = -Region[a];
								else if (mc1[a] >= lc[a])
									rshift[a] = Region[a];
								else
									rshift[a] = 0.0;
							}
							// Calculate the scalar cell index of the neighbor cell
							c1 = ((mc1[0] + lc[0]) % lc[0]) * lcyz
							    +((mc1[1] + lc[1]) % lc[1]) * lc[2]
							    +((mc1[2] + lc[2]) % lc[2]);
							// Skip this neighbor cell if empty
							if (head[c1] == EMPTY) continue;

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
									
													if (xAtom.sigma < 0 || xAtom.epsilon < 0 || yAtom.sigma < 0 || yAtom.epsilon < 0) continue;
										
													if(xAtom.id > yAtom.id) continue;										
				        												
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
														  
													if (r2 == 0.0) continue;								

													//calculate LJ energies //maybe do a calc_lj() call here when function is implemented!
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
	
	return totalEnergy + calcIntramolEnergy_NLC(enviro, molecules);
}
*/