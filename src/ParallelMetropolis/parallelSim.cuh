/*!\file
  \Class for parallel simulation caculate, change state of SimBox.
  \author David(Xiao Zhang).
 
  This file contains functions that are used to handle enviroments and common function
  for box.
 */
#ifndef PARALLELSIM_H
#define PARALLELSIM_H 

#include "GPUSimBox.cuh"
#include "Utilities/Opls_Scan.h"
#include "Utilities/Config_Scan.h"
#include "Utilities/metroUtil.h"
#include "Utilities/Zmatrix_Scan.h"
#include "Utilities/State_Scan.h"

#define THREADS_PER_BLOCK 128
#define PI 3.14159265
// Linked-cell neighbor list constants
#define NMAX 100000  /* Maximum number of atoms which can be simulated */
#define NCLMAX 10000 /* Maximum number of linked-list cells */
#define EMPTY -1

// boltzman constant in kcal mol-1 K-1
const double kBoltz = 0.00198717;

class ParallelSim
{
    private:
        GPUSimBox *box;
        int steps;
        float currentEnergy;
        float oldEnergy;
        int accepted;
        int rejected;
        int N;  //number of need compute
        int blocks; //number of per node
        size_t energySumSize;
        double *energySum_device;
        double *energySum_host;

    public:
        ParallelSim(GPUSimBox *initbox,int initsteps); 	
        ~ParallelSim(); 	
        float getcurrentEnergy(){return currentEnergy;}; 	
        int getaccepted() {return accepted;};
        int getrejected() {return rejected;};
        GPUSimBox * getGPUSimBox() {return box;};
        GPUSimBox * getdevGPUSimBox() {return box;};
        double *getdevEnergySum() { return energySum_device;};
        double *gethostEnergySum() { return energySum_host;};
        double calcEnergyWrapper(Molecule *molecules, Environment *enviro);
        double calcEnergyWrapper(GPUSimBox *box);
        double calcEnergyOnHost(Atom atom1, Atom atom2, Environment *enviro, Molecule *molecules);
        void runParallel(int steps);
        //double calcEnergy_NLC(Molecule *molecules, Environment *enviro);

    public:
        double wrapBox(double x, double box);
/*	__device__ int getXFromIndex(int idx);
	__device__ int getYFromIndex(int x, int idx);
	__device__ double makePeriodic(double x, double box);
	__device__ double calc_lj(Atom atom1, Atom atom2, Environment enviro);
	__device__ double calcCharge(Atom atom1, Atom atom2, Environment *enviro);
	__device__ double calcBlending(double d1, double d2);
	__device__ int getMoleculeFromAtomID(Atom a1, DeviceMolecule *dev_molecules, Environment enviro);
	__device__ double getFValue(Atom atom1, Atom atom2, DeviceMolecule *dev_molecules, Environment *enviro, Hop *hops);
	__device__ int hopGE3(int atom1, int atom2, DeviceMolecule dev_molecule, Hop *molecule_hops);
	*/
	Molecule* getMoleculeFromAtomIDHost(Atom a1, Molecule *molecules, Environment enviro);
 	int hopGE3Host(int atom1, int atom2, Molecule molecule);
    double getFValueHost(Atom atom1, Atom atom2, Molecule *molecules, Environment *enviro);
 	
};
 	
#endif