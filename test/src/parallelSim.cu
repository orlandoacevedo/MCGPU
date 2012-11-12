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

double ParallelSim::wrapBox(double x, double box){
    while(x >  box){
        x -= box;
    }
    while(x < 0){
        x += box;
    }

    return x;
}


double ParallelSim::calcEnergyWrapper(Molecule *molecules, Environment *enviro){
    
    Atom *atoms = (Atom *) malloc(sizeof(Atom) * enviro->numOfAtoms);
    int atomIndex = 0;
    double energy=0;
    
    for(int i = 0; i < enviro->numOfMolecules; i++){
        Molecule currentMolecule = molecules[i];
        for(int j = 0; j < currentMolecule.numOfAtoms; j++){
            atoms[atomIndex] = currentMolecule.atoms[j];
            atomIndex++;
        }
    }
    
    energy=calcEnergyWrapper(this->getGPUSimBox());
    free(atoms);
    
    return energy; 
}

//double ParallelSim::calcEnergyWrapper(Atom *atoms, Environment *enviro, Molecule *molecules){
double ParallelSim::calcEnergyWrapper(GPUSimBox *box){
    SimBox *innerbox=box->getSimBox();
    Atom *atoms=innerbox->getAtom();
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

    for(int i = 0; i < N; i++){

        //get atom IDs for each calculation
        int c = -2 * i;
        int discriminant = 1 - 4 * c;
        int qv = (-1 + sqrtf(discriminant)) / 2;
        int atomXid = qv + 1;
        
        int atomYid =  i - (atomXid * atomXid - atomXid) / 2;
        
        //check for stray calculations that returned invalid results
        if (isnan(energySum_host[i]) != 0 || isinf(energySum_host[i]) != 0){
            energySum_host[i] = calcEnergyOnHost(atoms[atomXid], atoms[atomYid], enviro, molecules);
        }
           
        //sum up energies 
        totalEnergy += energySum_host[i];
    }

    return totalEnergy;
}

double ParallelSim::calcEnergyOnHost(Atom atom1, Atom atom2, Environment *enviro, Molecule *molecules){
    //define terms in kcal/mol
    const double e = 332.06;
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
    if (r2 == 0.0){
        //lj_energy = 0.0;
        //charge_energy = 0.0;
        return 0.0;
    }
    else{
    	//calculate LJ energies
    	double sig2OverR2 = pow(sigma, 2) / r2;
        double sig6OverR6 = pow(sig2OverR2, 3);
    	double sig12OverR12 = pow(sig6OverR6, 2);
    	double lj_energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);
    	
    	//calculate Coulombic energies
    	double r = sqrt(r2);
    	double charge_energy = (atom2.charge * atom1.charge * e) / r;
    	
    	// Check for 1,3-nonbonded interaction
    	double fValue = 1.0;
		if (molecules != NULL)
        	fValue = getFValueHost(atom1, atom2, molecules, enviro);

    	return fValue * (lj_energy + charge_energy);
	}
}

double ParallelSim::getFValueHost(Atom atom1, Atom atom2, Molecule *molecules, Environment *enviro){
    Molecule *m1 = getMoleculeFromAtomIDHost(atom1, molecules, *enviro);
    Molecule *m2 = getMoleculeFromAtomIDHost(atom2, molecules, *enviro);
    Molecule molec = molecules[0];
    for(int i = 0; i < enviro->numOfMolecules; i++){
        if(molecules[i].id == m1->id){
            molec = molecules[i];
            break;
        }
    }

    if(m1->id != m2->id)
        return 1.0;
	else{
        int hops = hopGE3Host(atom1.id, atom2.id, *m1);
        if (hops == 3)
            return 0.5;
        else if (hops > 3)
            return 1.0;
        else
            return 0.0;
     }
}

int ParallelSim::hopGE3Host(int atom1, int atom2, Molecule molecule){
    for(int x=0; x< molecule.numOfHops; x++){
		      Hop myHop = molecule.hops[x];
				if((myHop.atom1==atom1 && myHop.atom2==atom2) ||
                        (myHop.atom1 == atom2 && myHop.atom2 == atom1) )
				    return myHop.hop;
	 }
	 return 0;
}

Molecule* ParallelSim::getMoleculeFromAtomIDHost(Atom a1, Molecule *molecules, Environment enviro){
    int atomId = a1.id;
    int currentIndex = enviro.numOfMolecules - 1;
    int molecId = molecules[currentIndex].id;
    while(atomId < molecId && currentIndex > 0){
        currentIndex -= 1;
        molecId = molecules[currentIndex].id;
    }
    return &molecules[currentIndex];

}

void ParallelSim::runParallel(int steps){
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
	 
    for(int move = 0; move < steps; move++){
        if (oldEnergy==0)
            oldEnergy = calcEnergyWrapper(molecules, enviro);

        int changeno=innerbox->ChangeMolecule();
				//box->CopyBoxtoDevice(innerbox);
				
        newEnergy = calcEnergyWrapper(this->getGPUSimBox());

        bool accept = false;

        if(newEnergy < oldEnergy){
            accept = true;            
        }
        else{
            double x = exp(-(newEnergy - oldEnergy) / kT);
            

            if(x >= randomFloat(0.0, 1.0)){
                accept = true;
            }
            else{
                accept = false;
            }
        }

        if(accept){
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

