/*!\file
  \Utility functions  for parallel Simulation, including Energy calculate and points to molocoles,only save all states
  \author David(Xiao Zhang).
 
  This file contains implement of functions running in GPU, but all allocation and free should be implemented outside
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <time.h>
#include "parallelSim.cuh"

__device__ int getXFromIndex(int idx){
    int c = -2 * idx;
    int discriminant = 1 - 4 * c;
    int qv = (-1 + sqrtf(discriminant)) / 2;
    return qv + 1;
}

__device__ int getYFromIndex(int x, int idx){
    return idx - (x * x - x) / 2;
}

__device__ double calcBlending(double d1, double d2){
    return sqrt(d1 * d2);
}

__device__ double makePeriodic(double x, double box){
    while(x < -0.5 * box){
        x += box;
    }

    while(x > 0.5 * box){
        x -= box;
    }

    return x;
}

double wrapBox(double x, double box){
    while(x >  box){
        x -= box;
    }
    while(x < 0){
        x += box;
    }

    return x;
}

__device__ double calc_lj(Atom atom1, Atom atom2, Environment enviro){
    //store LJ constants locally
    double sigma = calcBlending(atom1.sigma, atom2.sigma);
    double epsilon = calcBlending(atom1.epsilon, atom2.epsilon);
    
    //calculate difference in coordinates
    double deltaX = atom1.x - atom2.x;
    double deltaY = atom1.y - atom2.y;
    double deltaZ = atom1.z - atom2.z;

    //calculate distance between atoms
    deltaX = makePeriodic(deltaX, enviro.x);
    deltaY = makePeriodic(deltaY, enviro.y);
    deltaZ = makePeriodic(deltaZ, enviro.z);

    const double r2 = (deltaX * deltaX) +
                      (deltaY * deltaY) + 
                      (deltaZ * deltaZ);

    if (r2 == 0){
        return 0.0;
    }
    else{
    	//calculate terms
    	const double sig2OverR2 = pow(sigma, 2) / r2;
    	const double sig6OverR6 = pow(sig2OverR2, 3);
    	const double sig12OverR12 = pow(sig6OverR6, 2);
    	const double energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);
        return energy;
    }    
}

/*
double calcEnergyWrapper(Molecule *molecules, Environment *enviro){
    
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

*/
/*
//double ParallelSim::calcEnergyWrapper(Atom *atoms, Environment *enviro, Molecule *molecules){
double ParallelSim::calcEnergyWrapper(GPUSimBox *box){
    SimBox *innerbox=box->getSimBox();
    Atom *atoms=innerbox->getAom();
    Environment *enviro=innerbox->getEnviro();
    Molecule *molecules=innerbox->getMolecules();
    
    
    //setup CUDA storage
    double totalEnergy = 0.0;
    Atom *atoms_device;
//    double *energySum_device =getdevEnergySum();
//    double *energySum_host=gethostEnergySum();
    Environment *enviro_device;

    //calculate CUDA thread mgmt
    int N =(int) ( pow( (float) enviro->numOfAtoms,2)-enviro->numOfAtoms)/2;
    int blocks = N / THREADS_PER_BLOCK + (N % THREADS_PER_BLOCK == 0 ? 0 : 1); 

    //The number of bytes of shared memory per block of
    //size_t sharedSize = sizeof(double) * THREADS_PER_BLOCK;
    
    size_t atomSize = enviro->numOfAtoms * sizeof(Atom);
    size_t energySumSize = N * sizeof(double);
    
    //allocate memory on the device
    energySum_host = (double *) malloc(energySumSize);
    cudaMalloc((void **) &atoms_device, atomSize);
    cudaMalloc((void **) &energySum_device, energySumSize);
    cudaMalloc((void **) &enviro_device, sizeof(Environment));

    //copy data to the device
    cudaErrorCheck(cudaMemcpy(atoms_device, atoms, atomSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(enviro_device, enviro, sizeof(Environment), cudaMemcpyHostToDevice));
    memset(energySum_host,0,energySumSize);
    cudaErrorCheck(cudaMemcpy(energySum_device, energySum_host,  energySumSize, cudaMemcpyHostToDevice));
   

    if (molecules != NULL){
        int bondCount = 0;
        int angleCount = 0;
        int dihedralCount = 0;
        int hopCount = 0;
        for (int i = 0; i < enviro->numOfMolecules; i++){
            bondCount += molecules[i].numOfBonds;
            angleCount += molecules[i].numOfAngles;
            dihedralCount += molecules[i].numOfDihedrals;
            hopCount += molecules[i].numOfHops;
        }

        size_t dMolecSize = sizeof(DeviceMolecule) * enviro->numOfMolecules;
        size_t bondSize = sizeof(Bond) * bondCount;
        size_t angleSize = sizeof(Angle) * angleCount;
        size_t dihedralSize = sizeof(Dihedral) * dihedralCount;
        size_t hopSize = sizeof(Hop) * hopCount;
        
        DeviceMolecule *molec_d;
        Bond *bonds_d;
        Angle *angles_d;
        Dihedral *dihedrals_d;
        Hop *hops_d;
        
        cudaMalloc((void **) &molec_d, dMolecSize);
        cudaMalloc((void **) &bonds_d, bondSize);
        cudaMalloc((void **) &angles_d, angleSize);
        cudaMalloc((void **) &dihedrals_d, dihedralSize);
        cudaMalloc((void **) &hops_d, hopSize);


//        moleculeDeepCopyToDevice(molec_d, molecules, enviro->numOfMolecules, atoms_device, bonds_d, angles_d, dihedrals_d, hops_d);

//        calcEnergy <<<blocks, THREADS_PER_BLOCK>>>(atoms_device, enviro_device, energySum_device, molec_d, hops_d);

        calcEnergy <<<blocks, THREADS_PER_BLOCK>>>(this, atoms_device, enviro_device, energySum_device, molec_d, hops_d);
    cudaFree(molec_d);
    cudaFree(bonds_d);
    cudaFree(angles_d);
    cudaFree(dihedrals_d);
    cudaFree(hops_d);    
    }
    else{
        calcEnergy <<<blocks, THREADS_PER_BLOCK>>>(this,atoms_device, enviro_device, energySum_device, NULL, NULL);
    }
    
    cout<<"energySum_host"<<energySum_host<<"energySum_device"<<energySum_device<<"size"<<energySumSize<<endl;
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

    //cleanup
    cudaFree(atoms_device);
    cudaFree(energySum_device);
    free(energySum_host);

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
*/

/*
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

*/

__device__ double calcCharge(Atom atom1, Atom atom2, Environment *enviro){
    // conversion factor below for units in kcal/mol
    const double e = 332.06;
 
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
    
    if (r2 == 0.0){
        return 0.0;
    }
    else{
    	double r = sqrt(r2);
        return (atom1.charge * atom2.charge * e) / r;
    }
}


__device__ int getMoleculeFromAtomID(Atom a1, DeviceMolecule *dev_molecules, Environment enviro){
    int atomId = a1.id;
    int currentIndex = enviro.numOfMolecules - 1;
    int molecId = dev_molecules[currentIndex].id;
    while(atomId < molecId && currentIndex > 0){
        currentIndex -= 1;
        molecId = dev_molecules[currentIndex].id;
    }
    return molecId;

}

__device__ int hopGE3(int atom1, int atom2, DeviceMolecule dev_molecule,Hop *molecule_hops){
    Hop myHop;
    for(int x=0; x< dev_molecule.numOfHops; x++){
	    myHop = molecule_hops[x];
	    if((myHop.atom1==atom1 && myHop.atom2==atom2) || (myHop.atom1==atom2 && myHop.atom2==atom1))
	        return myHop.hop;
	 }
	 return 0;
}

__device__ double getFValue(Atom atom1, Atom atom2, DeviceMolecule *dev_molecules, Environment *enviro, Hop *dev_hops){
    int m1 = getMoleculeFromAtomID(atom1, dev_molecules, *enviro);
    int m2 = getMoleculeFromAtomID(atom2, dev_molecules, *enviro);
    
    if(m1 != m2){
        return 1.0;
    }
    else{
        int moleculeIndex = 0;
        for (int i = 0; i < enviro->numOfMolecules; i++){
            if (dev_molecules[i].id == m1)
                moleculeIndex = i;
        }
/*         size_t molecHopSize = sizeof(Hop) * dev_molecules[moleculeIndex].numOfHops;
        Hop *molecHops = (Hop *)malloc(molecHopSize);
        int hopStart = dev_molecules[moleculeIndex].hopStart;
        for (int i = 0; i < dev_molecules[moleculeIndex].numOfHops; i++){
            molecHops[i] = hops[hopStart + i];
        }*/
        
        int hopStart = dev_molecules[moleculeIndex].hopStart;
        int hopChain = hopGE3(atom1.id, atom2.id, dev_molecules[moleculeIndex], &dev_hops[hopStart]);
//        free(molecHops);
        if (hopChain == 3)
            return 0.5;
        else if (hopChain > 3)
            return 1.0;
        else
            return 0.0;
    } 
}



__global__ void calcEnergy(Atom *dev_atoms, Environment *dev_enviro, double *dev_energySum, DeviceMolecule *dev_molecules, Hop *dev_hops){
 
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    double lj_energy, charge_energy;

    int N =(int) ( pow( (float) dev_enviro->numOfAtoms,2)-dev_enviro->numOfAtoms)/2;

    if(idx < N ){
        //calculate the x and y positions in the Atom array
        int xAtom_pos, yAtom_pos;
        xAtom_pos = getXFromIndex(idx);
        yAtom_pos = getYFromIndex(xAtom_pos, idx);

        Atom xAtom, yAtom;
        xAtom = dev_atoms[xAtom_pos];
        yAtom = dev_atoms[yAtom_pos];

        if(xAtom.sigma < 0 || xAtom.epsilon < 0 || yAtom.sigma < 0 || yAtom.epsilon < 0){
            dev_energySum[idx] = 0.0;
        }
        else{
            lj_energy = calc_lj(xAtom,yAtom,*dev_enviro);
            charge_energy = calcCharge(xAtom, yAtom, dev_enviro);
            double fValue = 1.0;
            if (dev_molecules != NULL){
               fValue = getFValue(xAtom, yAtom, dev_molecules, dev_enviro, dev_hops);
            }
            
            dev_energySum[idx] = fValue * (lj_energy + charge_energy);
        }
    }
}