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

#define BLOCK_SIZE 512

ParallelSim::ParallelSim(GPUSimBox *initbox,int initsteps)
{
	box = initbox;
    steps = initsteps;
 	currentEnergy = 0;
    oldEnergy = 0;
    accepted = 0;
    rejected = 0;
	
	ptrs = (SimPointers*) malloc(sizeof(SimPointers));
	
	ptrs->innerbox = box->getSimBox();
	ptrs->envH = ptrs->innerbox->getEnviro();
	
	ptrs->atomsH = ptrs->innerbox->getAtoms();
	ptrs->moleculesH = ptrs->innerbox->getMolecules();
	
	ptrs->numA = ptrs->envH->numOfAtoms;
	ptrs->numM = ptrs->envH->numOfMolecules;
	ptrs->numEnergies = ptrs->numA * ptrs->numA;//This is an upper bound. May be able to be tightened.
	
	ptrs->molTrans = (Molecule*) malloc(ptrs->numM * sizeof(Molecule));
	ptrs->energiesH = (double*) malloc(ptrs->numEnergies * sizeof(double));
	
	cudaMalloc(&(ptrs->envD), sizeof(Environment));
	cudaMalloc(&(ptrs->atomsD), ptrs->numA * sizeof(Atom));
	cudaMalloc(&(ptrs->moleculesD), ptrs->numM * sizeof(Molecule));
	cudaMalloc(&(ptrs->energiesD), ptrs->numEnergies * sizeof(double));
	
	int i;
	
	//initialize energies
	for (i = 0; i < ptrs->numEnergies; i++)
	{
		ptrs->energiesH[i] = 0;
	}
	
	//sets up device molecules for transfer copies host molecules exactly except
	//for *atoms, which is translated to GPU pointers calculated here
	Atom *a = ptrs->atomsD;
	//upper bound on number of atoms in any molecule
	ptrs->maxMolSize = 0;
	for (i = 0; i < ptrs->numM; i++)
	{
		ptrs->molTrans[i].atoms = a;
		ptrs->molTrans[i].numOfAtoms = ptrs->moleculesH[i].numOfAtoms;
		a += ptrs->moleculesH[i].numOfAtoms;
		
		if (ptrs->moleculesH[i].numOfAtoms > ptrs->maxMolSize)
		{
			ptrs->maxMolSize = ptrs->moleculesH[i].numOfAtoms;
		}
	}
	
	//copy data to device
	cudaMemcpy(ptrs->envD, ptrs->envH, sizeof(Environment), cudaMemcpyHostToDevice);
	cudaMemcpy(ptrs->atomsD, ptrs->atomsH, ptrs->numA * sizeof(Atom), cudaMemcpyHostToDevice);
	cudaMemcpy(ptrs->moleculesD, ptrs->molTrans, ptrs->numM * sizeof(Molecule), cudaMemcpyHostToDevice);
	cudaMemcpy(ptrs->energiesD, ptrs->energiesH, ptrs->numEnergies * sizeof(double), cudaMemcpyHostToDevice);
}

ParallelSim::~ParallelSim()
{
    /*if (energySum_host!=NULL)
    {
        free(energySum_host);
        energySum_host=NULL;
    }
  
    if (energySum_device!=NULL)
    {
        cudaFree(energySum_device);
        energySum_device=NULL;
    }*/
}

void ParallelSim::writeChangeToDevice(int changeIdx)
{
	//create temp Molecule
	Molecule *changedMol = (Molecule*) malloc(sizeof(Molecule));
	
	//copy changed Molecule into temp Molecule
	//ready to be copied over to device, except that it still contains host pointer in .atoms
	memcpy(changedMol, ptrs->moleculesH + changeIdx, sizeof(Molecule));
	
	//get device pointer to device Atoms from device Molecule, and store in temp Molecule
	//temp Molecule.atoms will now contain a pointer to Atoms on device
	//this pointer never meant to be followed from host
	cudaMemcpy(&(changedMol->atoms), &((ptrs->moleculesD + changeIdx)->atoms), sizeof(Atom*), cudaMemcpyDeviceToHost);
	//copy changed molecule to device
	cudaMemcpy(ptrs->moleculesD + changeIdx, changedMol, sizeof(Molecule), cudaMemcpyHostToDevice);
	//copy changed atoms to device
	cudaMemcpy(ptrs->moleculesD[changeIdx].atoms, changedMol->atoms, changedMol->numOfAtoms * sizeof(Atom), cudaMemcpyHostToDevice);
}

double ParallelSim::calcSystemEnergy()
{
	double totalEnergy = 0;

	//for each molecule
	for (int mol = 0; mol < ptrs->numM; mol++)
	{
		totalEnergy += calcMolecularEnergyContribution(mol, mol);
	}
	
    return totalEnergy;
}

double ParallelSim::calcMolecularEnergyContribution(int molIdx, int startIdx)
{
	double totalEnergy = 0;
	int i;
	
	//initialize energies to 0
	for (i = 0; i < ptrs->numEnergies; i++)
	{
		ptrs->energiesH[i] = 0;
	}
	
	cudaMemcpy(ptrs->energiesD, ptrs->energiesH, ptrs->numEnergies * sizeof(double), cudaMemcpyHostToDevice);
	
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
						
	cudaMemcpy(ptrs->energiesH, ptrs->energiesD, ptrs->numEnergies * sizeof(double), cudaMemcpyDeviceToHost);
	
	for (i = 0; i < ptrs->numEnergies; i++)
	{
		totalEnergy += ptrs->energiesH[i];
	}
	
	return totalEnergy;
}

__global__ void calcInterMolecularEnergy(Molecule *molecules, int currentMol, int numM, int startIdx, Environment *enviro, double *energies, int segmentSize)
{
	int otherMol = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (otherMol >= startIdx && otherMol < numM && otherMol != currentMol)
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
			<<<molecules[currentMol].numOfAtoms, molecules[otherMol].numOfAtoms>>>
			(molecules, currentMol, otherMol, enviro, energies, segmentSize);
		}
	}
}

__global__ void calcInterAtomicEnergy(Molecule *molecules, int currentMol, int otherMol, Environment *enviro, double *energies, int segmentSize)
{
	Atom atom1 = molecules[currentMol].atoms[blockIdx.x],
		 atom2 = molecules[otherMol].atoms[threadIdx.x];
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
		
		totalEnergy += calc_lj(atom1, atom2, r2);
		totalEnergy += calcCharge(atom1.charge, atom2.charge, sqrt(r2));
		
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
		
		totalEnergy += calc_lj(atom1, atom2, r2);
		totalEnergy += calcCharge(atom1.charge, atom2.charge, sqrt(r2));
		
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
    	const double sig2OverR2 = (sigma*sigma) / r2;
		const double sig6OverR6 = (sig2OverR2*sig2OverR2*sig2OverR2);
    	const double sig12OverR12 = (sig6OverR6*sig6OverR6);
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

/*
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

double ParallelSim::getFValueHost(Atom atom1, Atom atom2, Molecule *molecules, Environment *enviro)
{
    Molecule *m1 = getMoleculeFromAtomIDHost(atom1, molecules, *enviro);
    Molecule *m2 = getMoleculeFromAtomIDHost(atom2, molecules, *enviro);
    Molecule molec = molecules[0];
    for(int i = 0; i < enviro->numOfMolecules; i++)
    {
        if(molecules[i].id == m1->id)
        {
            molec = molecules[i];
            break;
        }
    }

    if(m1->id != m2->id)
    {
        return 1.0;
    }
    else
    {
        int hops = hopGE3Host(atom1.id, atom2.id, *m1);
        if (hops == 3)
        {
            return 0.5;
        }
        else if (hops > 3)
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
    }
}

int ParallelSim::hopGE3Host(int atom1, int atom2, Molecule molecule)
{
    for(int x=0; x< molecule.numOfHops; x++)
    {
        Hop myHop = molecule.hops[x];
        if((myHop.atom1==atom1 && myHop.atom2==atom2) || (myHop.atom1 == atom2 && myHop.atom2 == atom1) )
        {
            return myHop.hop;
        }
     }
     return 0;
}

Molecule* ParallelSim::getMoleculeFromAtomIDHost(Atom a1, Molecule *molecules, Environment enviro)
{
    int atomId = a1.id;
    int currentIndex = enviro.numOfMolecules - 1;
    Molecule molec = molecules[currentIndex];
    int molecId = molec.atoms[0].id;
    while(atomId < molecId && currentIndex > 0)
    {
        currentIndex -= 1;
        molec = molecules[currentIndex];
        molecId = molec.atoms[0].id;
    }
    return &molecules[currentIndex];

}
*/

void ParallelSim::runParallel(int steps)
{	
    double temperature = ptrs->envH->temperature;
    double kT = kBoltz * temperature;
    double newEnergyCont, oldEnergyCont;

    if (oldEnergy == 0)
	{
		oldEnergy = calcSystemEnergy();
	}
	 
    for(int move = 0; move < steps; move++)
    {
            
        int changeIdx = ptrs->innerbox->chooseMolecule();

		oldEnergyCont = calcMolecularEnergyContribution(changeIdx);
		
		ptrs->innerbox->changeMolecule(changeIdx);
		writeChangeToDevice(changeIdx);
		
		newEnergyCont = calcMolecularEnergyContribution(changeIdx);

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
            ptrs->innerbox->Rollback(changeIdx);
        }
    }
    currentEnergy=oldEnergy;
}