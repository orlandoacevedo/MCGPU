/*
	Contains calculations for ParallelBox
	Same functions as SerialCalcs with function qualifiers and CUDA code

	Author: Nathan Coleman
*/

#ifndef PARALLELCALCS_CUH
#define PARALLELCALCS_CUH

#include "Metropolis/Utilities/StructLibrary.h"

__global__ void calcInterMolecularEnergy(Molecule *molecules, int currentMol, int numM, Environment *environment, double *energies, int segmentSize);
__global__ void calcInterAtomicEnergy(Molecule *molecules, int currentMol, int otherMol, Environment *environment, double *energies, int segmentSize);
__global__ void calcIntraMolecularEnergy(Molecule *molecules, int currentMol, int numE, Environment *environment, double *energies, int segmentSize);
__device__ double calc_lj(Atom atom1, Atom atom2, double r2);
__device__ double calcCharge(double charge1, double charge2, double r);
__device__ double makePeriodic(double x, double box);
__device__ double calcBlending(double d1, double d2);
__device__ int getXFromIndex(int idx);
__device__ int getYFromIndex(int x, int idx);

#endif