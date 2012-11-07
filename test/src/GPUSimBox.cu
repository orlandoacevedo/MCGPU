/*!\file
  \Class for simulation Box used for GPU, including Enviroments and points to molocoles,only save all states
  \author David(Xiao Zhang).
 
  This file contains implement of SimBox that are used to handle enviroments and common function
  for box.
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "GPUSimBox.cuh"

#define THREADS_PER_BLOCK 128
#define PI 3.14159265

using namespace std;

void cudaAssert(const cudaError err, const char *file, const int line)
{ 
    if( cudaSuccess != err) {                                                
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        
                file, line, cudaGetErrorString(err) );
    } 
}

GPUSimBox::GPUSimBox(Config_Scan configScan)
{
   molec_d = NULL;
   bonds_d = NULL;
   angles_d = NULL;
   dihedrals_d = NULL;
   hops_d = NULL;
   atoms_device = NULL;
   enviro_device = NULL; 
   	
   atomSize = 0;
   dMolecSize= 0;	
   bondSize= 0;
   angleSize= 0;
   dihedralSize= 0;
   hopSize= 0;	
   
   innerbox = new SimBox(configScan);
   initGPUSimBox(innerbox);
}

int GPUSimBox::initGPUSimBox(SimBox *hostbox)
{

		Environment *enviro=hostbox->getEnviro();
		Molecule *molecules=hostbox->getMolecules();
		
	  atomSize = enviro->numOfAtoms * sizeof(Atom);
    dMolecSize = enviro->numOfMolecules *sizeof(struct DeviceMolecule) ;
    
    //allocate memory on the device
    struct DeviceMolecule *dMole_h=(struct DeviceMolecule *)malloc(dMolecSize);
    
    cudaMalloc((void **) &atoms_device, atomSize);
    cudaMalloc((void **) &enviro_device, sizeof(Environment));
    cudaMalloc((void **) &molec_d, dMolecSize);


    if (molecules != NULL){
        int atomCount = 0;
        int bondCount = 0;
        int angleCount = 0;
        int dihedralCount = 0;
        int hopCount = 0;
        for (int i = 0; i < enviro->numOfMolecules; i++){
            dMole_h[i].id=i;
            
            dMole_h[i].atomStart=atomCount;
            atomCount += molecules[i].numOfAtoms;
            dMole_h[i].numOfAtoms=molecules[i].numOfAtoms;
            
            dMole_h[i].bondStart = bondCount;
            bondCount += molecules[i].numOfBonds;
            dMole_h[i].numOfBonds=molecules[i].numOfBonds;
            
						dMole_h[i].angleStart = angleCount;
            angleCount += molecules[i].numOfAngles;
            dMole_h[i].numOfAngles=molecules[i].numOfAngles;
            
						dMole_h[i].angleStart = angleCount;
            dihedralCount += molecules[i].numOfDihedrals;
            dMole_h[i].numOfDihedrals=molecules[i].numOfDihedrals;

						dMole_h[i].angleStart = angleCount;
            hopCount += molecules[i].numOfHops;
            dMole_h[i].numOfHops=molecules[i].numOfHops;

        }

    //copy data to the device    
    cudaErrorCheck(cudaMemcpy(enviro_device, enviro, sizeof(Environment), cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(molec_d, dMole_h, dMolecSize, cudaMemcpyHostToDevice));

        bondSize = sizeof(Bond) * bondCount;
        angleSize = sizeof(Angle) * angleCount;
        dihedralSize = sizeof(Dihedral) * dihedralCount;
        hopSize = sizeof(Hop) * hopCount;
                
        cudaMalloc((void **) &molec_d, dMolecSize);
        cudaMalloc((void **) &bonds_d, bondSize);
        cudaMalloc((void **) &angles_d, angleSize);
        cudaMalloc((void **) &dihedrals_d, dihedralSize);
        cudaMalloc((void **) &hops_d, hopSize);
    }
        
   return 0;
}

GPUSimBox::~GPUSimBox()
{
		cudaFREE(bonds_d);
		cudaFREE(angles_d);
		cudaFREE(dihedrals_d);
		cudaFREE(hops_d);
		cudaFREE(atoms_device);
		cudaFREE(enviro_device);
		cudaFREE(molec_d);
		
		delete innerbox;
}

int GPUSimBox::CopyBoxtoHost(SimBox hostbox)
{
		Environment *enviro=hostbox.getEnviro();
		Molecule *molecules=hostbox.getMolecules();

    cudaErrorCheck(cudaMemcpy(molecules[0].atoms,atoms_device,  atomSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(molecules[0].bonds, bonds_d, bondSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(molecules[0].angles, angles_d, angleSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy( molecules[0].dihedrals, dihedrals_d, dihedralSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(molecules[0].hops, hops_d, hopSize, cudaMemcpyHostToDevice));
    
    return 0;

}
int GPUSimBox::CopyBoxtoDevice(SimBox hostbox)
{
		Environment *enviro=hostbox.getEnviro();
		Molecule *molecules=hostbox.getMolecules();

    cudaErrorCheck(cudaMemcpy(atoms_device, molecules[0].atoms, atomSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(bonds_d, molecules[0].bonds, bondSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(angles_d, molecules[0].angles, angleSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(dihedrals_d, molecules[0].dihedrals, dihedralSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(hops_d, molecules[0].hops, hopSize, cudaMemcpyHostToDevice));
    
    return 0;

}