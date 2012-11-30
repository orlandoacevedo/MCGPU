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
    int molecId = dev_molecules[currentIndex].atomStart;
    while(atomId < molecId && currentIndex > 0){
        currentIndex -= 1;
        molecId = dev_molecules[currentIndex].atomStart;
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
//            lj_energy=0;charge_energy=0;
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