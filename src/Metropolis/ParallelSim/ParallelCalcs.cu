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
#define THREADS_PER_BLOCK 192

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
    box->neighborList = new NeighborList(box->molecules, box->environment);
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

Real ParallelCalcs::calcSystemEnergy_NLC(Box *box){ 
	NeighborList *nl = box->neighborList;
	Molecule *molecules = box->getMolecules();
	Environment *enviro = box->getEnvironment();
	
	int maxPairs = enviro->numOfMolecules * enviro->numOfMolecules;
	int pair_i[maxPairs];
	int pair_j[maxPairs];
	int iterater_i = 0;

	int vectorCells[3];
	for (vectorCells[0] = 0; vectorCells[0] < nl->numCells[0]; (vectorCells[0])++)
	{
		for (vectorCells[1] = 0; vectorCells[1] < nl->numCells[1]; (vectorCells[1])++)
		{
			for (vectorCells[2] = 0; vectorCells[2] < nl->numCells[2]; (vectorCells[2])++)
			{
				// Calculate a scalar cell index
				int c = vectorCells[0]*nl->numCellsYZ + vectorCells[1]*nl->numCells[2] + vectorCells[2];
				// Skip this cell if empty
				if (nl->head[c] == EMPTY) continue;
				// Scan the neighbor cells (including itself) of cell c
				int neighborCells[3];
				for (neighborCells[0] = vectorCells[0]-1; neighborCells[0] <= vectorCells[0]+1; (neighborCells[0])++)
					for (neighborCells[1] = vectorCells[1]-1; neighborCells[1] <= vectorCells[1]+1; (neighborCells[1])++)
						for (neighborCells[2] = vectorCells[2]-1; neighborCells[2] <= vectorCells[2]+1; (neighborCells[2])++)
						{
							// Periodic boundary condition by shifting coordinates
							Real rshift[3];
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
							if (nl->head[c1] == EMPTY)
							{
								continue;
							}
							// Scan atom i in cell c
							int i = nl->head[c];
							while (i != EMPTY)
							{
								// Scan atom j in cell c1
								int j = nl->head[c1];
								while (j != EMPTY)
								{
									bool included = false;
									// Avoid double counting of pairs
									if (i < j)
									{	
										std::vector<int> currentMolPrimaryIndexArray = (*(*(enviro->primaryAtomIndexArray))[molecules[i].type]);
										std::vector<int> otherMolPrimaryIndexArray;
										if (molecules[i].type == molecules[j].type)
										{
											otherMolPrimaryIndexArray = currentMolPrimaryIndexArray;
										}
										else 
										{
											otherMolPrimaryIndexArray = (*(*(enviro->primaryAtomIndexArray))[molecules[j].type]);
										}
										for (int i1 = 0; i1 < currentMolPrimaryIndexArray.size(); i1++)
										{
											for (int i2 = 0; i2 < otherMolPrimaryIndexArray.size(); i2++)
											{
												int primaryIndex1 = currentMolPrimaryIndexArray[i1];
												int primaryIndex2 = otherMolPrimaryIndexArray[i2];
												Atom atom1 = molecules[i].atoms[primaryIndex1];
												Atom atom2 = molecules[j].atoms[primaryIndex2];
												Real dr[3];		  /* Pair vector dr = atom[i]-atom[j] */
												dr[0] = atom1.x - (atom2.x + rshift[0]);
												dr[1] = atom1.y - (atom2.y + rshift[1]);
												dr[2] = atom1.z - (atom2.z + rshift[2]);
												Real rr = (dr[0] * dr[0]) + (dr[1] * dr[1]) + (dr[2] * dr[2]);			
												
												// Calculate energy for entire molecule interaction if rij < Cutoff for atom index
												if (rr < nl->rrCut)
												{	
													//totalEnergy += calcInterMolecularEnergy(molecules, i, j, enviro, subLJ, subCharge) * fValue;
													pair_i[iterater_i] = i;
													pair_j[iterater_i] = j;
													iterater_i++;
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
									
									j = nl->linkedCellList[j];
								} /* Endwhile j not empty */
								
								i = nl->linkedCellList[i];
							} /* Endwhile i not empty */
						} /* Endfor neighbor cells, c1 */
			} /* Endfor central cell, c */
		}
	}
	int *d_pair_i;
	int *d_pair_j;
	cudaMalloc((void **)&d_pair_i, sizeof(int)*maxPairs);
	cudaMalloc((void **)&d_pair_j, sizeof(int)*maxPairs);
	cudaMemcpy(d_pair_i, pair_i, sizeof(int)*maxPairs, cudaMemcpyHostToDevice);
	cudaMemcpy(d_pair_j, pair_j, sizeof(int)*maxPairs, cudaMemcpyHostToDevice);
	thrust::device_vector<Real> part_energy(maxPairs, 0);//this will store the result
	ParallelBox *pBox = (ParallelBox*) box;
	if (pBox == NULL)
	{
		return 0;
	}
	MoleculeData *d_molecules = pBox->moleculesD;
	AtomData *d_atoms = pBox->atomsD;
	Environment *d_enviro = pBox->environmentD;
	Real *raw_ptr = thrust::raw_pointer_cast(&part_energy[0]);
	int blocksPerGrid = 53;
	calcEnergy_NLC<<<blocksPerGrid, THREADS_PER_BLOCK>>>(d_pair_i, d_pair_j, raw_ptr, d_molecules, d_atoms, d_enviro, iterater_i);
	Real total_energy = thrust::reduce(part_energy.begin(), part_energy.end());
	cudaFree(d_pair_i);
	cudaFree(d_pair_j);
	return total_energy;// + calcIntramolEnergy_NLC(enviro, pBox->moleculesD, pBox->atomsD);
}

__global__ void ParallelCalcs::calcEnergy_NLC(int* d_pair_i, int* d_pair_j, Real *part_energy, MoleculeData *molecules, AtomData *atoms, Environment *enviro, int limit)
{
 
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i < limit){
		part_energy[i] = part_energy[i] + calcInterMolecularEnergy(molecules, atoms, d_pair_i[i], d_pair_j[i], enviro) * 1.0;
	}
}

__device__ Real ParallelCalcs::calcInterMolecularEnergy(MoleculeData *molecules, AtomData *atoms, int mol1, int mol2, Environment *enviro)
{
	Real totalEnergy = 0;
	for (int i = 0; i < molecules->numOfAtoms[mol1]; i++)
	{
		for (int j = 0; j < molecules->numOfAtoms[mol2]; j++)
		{
			int atom1 = molecules->atomsIdx[mol1] + i;
			int atom2 = molecules->atomsIdx[mol2] + j;
			if (atoms->sigma[atom1] >= 0 && atoms->epsilon[atom1] >= 0 && atoms->sigma[atom2] >= 0 && atoms->epsilon[atom2] >= 0)
			{
				//calculate squared distance between atoms 
				Real r2 = calcAtomDist(atoms, atom1, atom2, enviro);
				totalEnergy += calc_lj(atoms, atom1, atom2, r2);
				totalEnergy += calcCharge(atoms->charge[atom1], atoms->charge[atom2], sqrt(r2));
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
	calcInterMolecularEnergy<<<validEnergies / BATCH_BLOCK + 1, BATCH_BLOCK>>>
	(box->moleculesD, box->atomsD, molIdx, box->environmentD, box->energiesD, validEnergies, box->molBatchD, box->maxMolSize);
	
	//Using Thrust here for a sum reduction on all of the individual energy contributions in box->energiesD.
	thrust::device_ptr<Real> energiesOnDevice = thrust::device_pointer_cast(&box->energiesD[0]);
	double reduction = thrust::reduce(energiesOnDevice, energiesOnDevice + validEnergies, (Real) 0, thrust::plus<Real>());
	cudaMemset(box->energiesD, 0, validEnergies * sizeof(Real));
	return reduction;
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
	
	}
}

Real ParallelCalcs::calcMolecularEnergyContribution_NLC(Box *box, int currentMol) {
	NeighborList *nl = box->neighborList; 
	Molecule *molecules = box->getMolecules();
	Environment *enviro = box->getEnvironment();
	thrust::device_vector<Real> part_energy(19200, 0);
	Real *raw_ptr = thrust::raw_pointer_cast(&part_energy[0]);

	ParallelBox *pBox = (ParallelBox*) box;
	if (pBox == NULL)
	{
		return 0;
	}
	MoleculeData *d_molecules = pBox->moleculesD;
	AtomData *d_atoms = pBox->atomsD;
	Environment *d_enviro = pbox->environmentD;
	int counter = 0;//count how many iterations 
	
	// Find the vector cell for the currentMol (based on 1st primary index)
	std::vector<int> currentMolPrimaryIndexArray = (*(*(enviro->primaryAtomIndexArray))[molecules[currentMol].type]);
	int primaryIndex = currentMolPrimaryIndexArray[0]; // Use first primary index to determine cell placement
		
	int vectorCells[3];			
	vectorCells[0] = molecules[currentMol].atoms[primaryIndex].x / nl->lengthCell[0]; 
	vectorCells[1] = molecules[currentMol].atoms[primaryIndex].y / nl->lengthCell[1];
	vectorCells[2] = molecules[currentMol].atoms[primaryIndex].z / nl->lengthCell[2];
	
	// Scan the neighbor cells (including itself) of cell c
	int neighborCells[3];			/* Neighbor cells */
	//#pragma omp parallel for collapse(3)
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

				// Scan atom otherMol in cell c1
				int otherMol = nl->head[c1];
				while (otherMol != EMPTY)
				{
					bool included = false;

					// For other molecules (and if total system calc, then avoid double counting of pairs)
					if (currentMol != otherMol)
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
									
									calcInterMolecularEnergy_NLC<<<1, THREADS_PER_BLOCK>>>(d_molecules, d_atoms, d_enviro, currentMol, otherMol, counter, raw_ptr, pBox->maxMolSize);
									counter++;
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
	
	Real total_energy = thrust::reduce(part_energy.begin(), part_energy.end());
	
	return total_energy;
}

__global__ void ParallelCalcs::calcInterMolecularEnergy_NLC(MoleculeData *molecules, AtomData *atoms, Environment *enviro, int currentMol, int otherMol, int counter, Real *energies, int maxMolSize)
{
	int energyIdx = threadIdx.x;

	//check validity of thread

	//get atom pair for this thread
	int x = energyIdx / maxMolSize;
	int y = energyIdx % maxMolSize;

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
			energies[energyIdx + blockDim.x*counter] = totalEnergy;
		}
	}

}

__global__ void ParallelCalcs::calcInterMolecularEnergy(MoleculeData *molecules, AtomData *atoms, int currentMol, Environment *enviro, Real *energies, int energyCount, int *molBatch, int maxMolSize)
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

__host__ __device__ Real ParallelCalcs::calc_lj(AtomData *atoms, int atom1, int atom2, Real r2)
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

__device__ Real ParallelCalcs::calcAtomDist(AtomData *atoms, int atomIdx1, int atomIdx2, Environment *enviro)
{
	//calculate difference in coordinates
	Real deltaX = makePeriodic(atoms->x[atomIdx1] - atoms->x[atomIdx2], enviro->x);
	Real deltaY = makePeriodic(atoms->y[atomIdx1] - atoms->y[atomIdx2], enviro->y);
	Real deltaZ = makePeriodic(atoms->z[atomIdx1] - atoms->z[atomIdx2], enviro->z);
				
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
