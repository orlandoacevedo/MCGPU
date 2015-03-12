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
Real ParallelCalcs::calcIntramolEnergy_NLC(Environment *enviro, Molecule *molecules)
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

Real ParallelCalcs::calcSystemEnergy(Box *box){ 
	            Molecule *molecules = box->getMolecules();
		        Environment *enviro = box->getEnvironment();
	            const Real Region[3] = {enviro->x, enviro->y, enviro->z};
	            int head[NCLMAX];
	            int lscl[NMAX];
	            int lc[3];            	/* Number of cells in the x|y|z direction */
	            int mc[3];			  	/* Vector cell */
	            Real rc[3];

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

	            Real *d_totalEnergy;
	            Real totalEnergy;
	            Real oldEnergy;
	            Molecule *d_molecules;
	            Environment *d_enviro;
	            int *d_head;
	            int *d_lscl;

	            cudaMalloc(&d_molecules, enviro->numOfMolecules*sizeof(Molecule));
	            cudaMalloc(&d_enviro, sizeof(Environment));
	            cudaMalloc(&d_head, sizeof(int)*NCLMAX);
	            cudaMalloc(&d_lscl, sizeof(int)*NMAX);
	            cudaMalloc(&d_totalEnergy, sizeof(Real));

	            cudaMemcpy(d_enviro, enviro, sizeof(Environment), cudaMemcpyHostToDevice);
	            cudaMemcpy(d_molecules, molecules, enviro->numOfMolecules*sizeof(Molecule), cudaMemcpyHostToDevice);
	            cudaMemcpy(d_head, head, sizeof(int)*NCLMAX, cudaMemcpyHostToDevice);
	            cudaMemcpy(d_lscl, lscl, sizeof(int)*NMAX, cudaMemcpyHostToDevice);
                dim3 dimGrid(lc[0], lc[1], lc[2]);
                dim3 dimBlock(3, 3, 3);
	            calcEnergy_NLC<<<dimGrid, dimBlock>>>(d_molecules, d_enviro, d_head, d_lscl, d_totalEnergy);
				cudaMemcpy(&totalEnergy, d_totalEnergy, sizeof(Real), cudaMemcpyDeviceToHost);

	            oldEnergy = totalEnergy + calcIntramolEnergy_NLC(enviro, molecules);
	            cudaFree(d_totalEnergy);
	            cudaFree(d_molecules);
	            cudaFree(d_enviro);
	            cudaFree(d_head);
	            cudaFree(d_lscl);
	            return oldEnergy;
}

__global__ void ParallelCalcs::calcEnergy_NLC(Molecule *molecules, Environment *enviro, int *head, int *lscl, Real *totalEnergy)
{
	// Variables for linked-cell neighbor list	
	int lc[3];            	/* Number of cells in the x|y|z direction */
	Real rc[3];         	/* Length of a cell in the x|y|z direction */
	//int head[NCLMAX];    	/* Headers for the linked cell lists */
	int mc[3];			  	/* Vector cell */
	//int lscl[NMAX];       	/* Linked cell lists */
	int mc1[3];				/* Neighbor cells */
	Real rshift[3];	  		/* Shift coordinates for periodicity */
	Real Region[3] = {enviro->x, enviro->y, enviro->z};  /* MD box lengths */
	int c1;				  	/* Used for scalar cell index */
	Real dr[3];		  		/* Pair vector dr = atom[i]-atom[j] */
	Real rrCut = enviro->cutoff * enviro->cutoff;	/* Cutoff squared */
	Real rr;			  	/* Distance between atoms */
	//Real lj_energy;			/* Holds current Lennard-Jones energy */
	//Real charge_energy;		/* Holds current coulombic charge energy */
	Real fValue = 1.0;		/* Holds 1,4-fudge factor value */
	//Real totalEnergy = 0.0;	/* Total nonbonded energy x fudge factor */
	
    mc[0] = blockIdx.x;
    mc[1] = blockIdx.y;
    mc[2] = blockIdx.z;
    
    mc1[0] = mc[0] + (threadIdx.x - 1);
    mc1[1] = mc[1] + (threadIdx.y - 1);
    mc1[2] = mc[2] + (threadIdx.z - 1);
	// Compute the # of cells for linked cell lists
	for (int k=0; k<3; k++)
	{
		lc[k] = Region[k] / enviro->cutoff; 
		rc[k] = Region[k] / lc[k];
	}
		
  /* Make a linked-cell list, lscl--------------------------------------------*/
	int lcyz = lc[1]*lc[2];
	//int lcxyz = lc[0]*lcyz;

  /* Calculate pair interaction-----------------------------------------------*/
		
	// Scan inner cells
    // Calculate a scalar cell index
    if(mc[0] < lc[0]&&mc[1] < lc[1]&&mc[2]<lc[2]){
        int c = mc[0]*lcyz + mc[1]*lc[2] + mc[2];
				// Skip this cell if empty		
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
											(*totalEnergy) += calcInterMolecularEnergy(molecules, i, j, enviro) * fValue;
											
										} /* Endif rr < rrCut */
									} /* Endif i<j */
									
									j = lscl[j];
								} /* Endwhile j not empty */

								i = lscl[i];
							} /* Endwhile i not empty */
    }
		
}

__device__ Real ParallelCalcs::calcInterMolecularEnergy(Molecule *molecules, int mol1, int mol2, Environment *enviro)
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
		//find primary atom indices for this pair of molecules
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
		}
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

__device__ __host__ Real ParallelCalcs::calc_lj(Atom atom1, Atom atom2, Real r2)
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

__device__ __host__ Real ParallelCalcs::calcAtomDist(Atom atom1, Atom atom2, Environment *enviro)
{
	//calculate difference in coordinates
	Real deltaX = makePeriodic(atom1.x - atom2.x, enviro->x);
	Real deltaY = makePeriodic(atom1.y - atom2.y, enviro->y);
	Real deltaZ = makePeriodic(atom1.z - atom2.z, enviro->z);
				
	//calculate squared distance (r2 value) and return
	return (deltaX * deltaX) + (deltaY * deltaY) + (deltaZ * deltaZ);
}

__device__ __host__ Real ParallelCalcs::calcCharge(Real charge1, Real charge2, Real r)
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

__device__ __host__ Real ParallelCalcs::makePeriodic(Real x, Real boxDim)
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

__device__ __host__ Real ParallelCalcs::calcBlending(Real d1, Real d2)
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
