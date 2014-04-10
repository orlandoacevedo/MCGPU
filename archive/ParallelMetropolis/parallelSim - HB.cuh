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

struct SimPointers
{
	SimBox *innerbox;
    Environment *envH, *envD;
    Atom *atomsH, *atomsD;
    Molecule *moleculesH, *moleculesD, *molTrans;
	int numA, numM, numEnergies, maxMolSize;
	int *molBatchH, *molBatchD;
	double *energiesH, *energiesD;
};

__global__ void calcInterAtomicEnergy(Molecule *molecules, int curentMol, Environment *enviro, double *energies, int numEnergies, int *molBatch, int maxMolSize);
__global__ void aggregateEnergies(double *energies, int numEnergies, int interval, int batchSize);
__device__ double calc_lj(Atom atom1, Atom atom2, double r2);
__device__ double calcCharge(double charge1, double charge2, double r);
__device__ double makePeriodic(double x, double box);
__device__ double calcBlending(double d1, double d2);
__device__ int getXFromIndex(int idx);
__device__ int getYFromIndex(int x, int idx);
		
class ParallelSim
{
    private:
        GPUSimBox *box;
		SimPointers *ptrs;
        int steps;
        float currentEnergy;
        float oldEnergy;
        int accepted;
        int rejected;

    public:
        ParallelSim(GPUSimBox *initbox,int initsteps); 	
        ~ParallelSim(); 	
        float getcurrentEnergy(){return currentEnergy;}; 	
        int getaccepted() {return accepted;};
        int getrejected() {return rejected;};
        GPUSimBox * getGPUSimBox() {return box;};
        GPUSimBox * getdevGPUSimBox() {return box;};
		void writeChangeToDevice(int changeIdx);
		double calcSystemEnergy();
		double calcMolecularEnergyContribution(int molIdx, int startIdx = 0);
		double calcBatchEnergy(int numMols, int molIdx);
		double getEnergyFromDevice();
		double makePeriodicH(double x, double box);
        void runParallel(int steps);

    public:
        double wrapBox(double x, double box);
    
	Molecule* getMoleculeFromAtomIDHost(Atom a1, Molecule *molecules, Environment enviro);
 	int hopGE3Host(int atom1, int atom2, Molecule molecule);
    double getFValueHost(Atom atom1, Atom atom2, Molecule *molecules, Environment *enviro);
 	
};
 	
#endif