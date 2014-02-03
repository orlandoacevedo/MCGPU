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
    
    cudaMalloc((void **) &enviro_device, sizeof(Environment));
    cudaMalloc((void **) &molec_d, dMolecSize);
    cudaMalloc((void **) &atoms_device, atomSize);
    
    int atomCount = 0;
    int bondCount = 0;
    int angleCount = 0;
    int dihedralCount = 0;
    int hopCount = 0;
   
    if (molecules != NULL)
    {

        for (int i = 0; i < enviro->numOfMolecules; i++)
        {
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
            
						dMole_h[i].dihedralStart = dihedralCount;
            dihedralCount += molecules[i].numOfDihedrals;
            dMole_h[i].numOfDihedrals=molecules[i].numOfDihedrals;

						dMole_h[i].hopStart = hopCount;
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
                
        //cudaMalloc((void **) &molec_d, dMolecSize);
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

int GPUSimBox::CopyBoxtoHost(SimBox *hostbox)
{
	Environment *enviro=hostbox->getEnviro();
	Molecule *molecules=hostbox->getMolecules();

    cudaErrorCheck(cudaMemcpy(molecules[0].atoms,atoms_device,  atomSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(molecules[0].bonds, bonds_d, bondSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(molecules[0].angles, angles_d, angleSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy( molecules[0].dihedrals, dihedrals_d, dihedralSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(molecules[0].hops, hops_d, hopSize, cudaMemcpyHostToDevice));
    
    return 0;

}
int GPUSimBox::CopyBoxtoDevice(SimBox *hostbox)
{
	Environment *enviro=hostbox->getEnviro();
	Molecule *molecules=hostbox->getMolecules();

    cudaErrorCheck(cudaMemcpy(atoms_device, molecules[0].atoms, atomSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(bonds_d, molecules[0].bonds, bondSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(angles_d, molecules[0].angles, angleSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(dihedrals_d, molecules[0].dihedrals, dihedralSize, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(hops_d, molecules[0].hops, hopSize, cudaMemcpyHostToDevice));
    
    return 0;

}

/*
//***************************
//These methods have been transplanted from parallelSim.cu
// to match the format of SimBox.cpp. [Feb 03, 2014]
// Be sure to check for proper fit & finish of all variables, etc.
//***************************
*/

/*
--Allows the sim to have a given item wrap back around when it hits the border of a given box.
@param x: the coordinate being checked; may be out-of-bounds, and if so, will be fixed
@param box: the maximum boundary for the given item, based on the boundaries of the box along that axis
@return: returns the new value of "x", corrected for the maximum size of the box.
[end comments]
*/
double GPUSimBox::wrapBox(double x, double box)
{
    while(x >  box)
    {
        x -= box;
    }
    while(x < 0)
    {
        x += box;
    }

    return x;
}

double GPUSimBox::getFValueHost(Atom atom1, Atom atom2, Molecule *molecules, Environment *enviro)
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

int GPUSimBox::hopGE3Host(int atom1, int atom2, Molecule molecule)
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

Molecule* GPUSimBox::getMoleculeFromAtomIDHost(Atom a1, Molecule *molecules, Environment enviro)
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


int GPUSimBox::getXFromIndex(int idx)
{
    int c = -2 * idx;
    int discriminant = 1 - 4 * c;
    int qv = (-1 + sqrtf(discriminant)) / 2;
    return qv + 1;
}

int GPUSimBox::getYFromIndex(int x, int idx)
{
    return idx - (x * x - x) / 2;
}

double GPUSimBox::makePeriodic(double x, double box)
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
