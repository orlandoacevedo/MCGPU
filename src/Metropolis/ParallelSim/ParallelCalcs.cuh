/*
	Contains calculations for ParallelBox
	Same functions as SerialCalcs with function qualifiers and CUDA code

	Author: Nathan Coleman
*/

#ifndef PARALLELCALCS_CUH
#define PARALLELCALCS_CUH

#include "Metropolis/Utilities/StructLibrary.h"
#include "Metropolis/DataTypes.h"

__global__ void calcInterMolecularEnergy(Molecule *molecules, int currentMol, int numM, Environment *environment, Real *energies, int segmentSize);
__global__ void calcInterAtomicEnergy(Molecule *molecules, int currentMol, int otherMol, Environment *environment, Real *energies, int segmentSize);
__global__ void calcIntraMolecularEnergy(Molecule *molecules, int currentMol, int numE, Environment *environment, Real *energies, int segmentSize);
__device__ Real calc_lj(Atom atom1, Atom atom2, Real r2);
__device__ Real calcCharge(Real charge1, Real charge2, Real r);
__device__ Real makePeriodic(Real x, Real box);
__device__ Real calcBlending(Real d1, Real d2);
__device__ int getXFromIndex(int idx);
__device__ int getYFromIndex(int x, int idx);

#endif