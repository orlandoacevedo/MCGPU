/*
	Contains the methods required to calculate energies in parallel.

	Created: February 21, 2014
	
	-> February 26, by Albert Wallace
	-> March 28, by Joshua Mosby
	-> April 21, by Nathan Coleman
*/

#ifndef PARALLELCALCS_CUH
#define PARALLELCALCS_CUH

#include "Metropolis/Utilities/StructLibrary.h"
#include "Metropolis/Utilities/Coalesced_Structs.h"
#include "Metropolis/DataTypes.h"
#include "ParallelBox.cuh"

namespace ParallelCalcs
{
	
	/// Creates a batch of molecule IDs within the cutoff distance of
	///   the chosen molecule, and returns the batch size.
	/// @param box A ParallelBox containing the molecule data.
	/// @param currentMol The index of the current changed molecule.
	/// @param startIdx The optional starting index for other molecules.
	///   Used for system energy calculation.
	/// @return Returns batch size (number of molecules within cutoff).
	int createMolBatch(ParallelBox *box, int currentMol, int startIdx);
	
	/// Given a box with a filled-in molecule batch, calculate the
	///   inter-molecular energy contribution of a given molecule,
	///   with every molecule specified in the batch.
	/// @param box A ParallelBox containing the molecule data.
	/// @param numMols The number of molecules within the cutoff.
	/// @param molIdx The index of the current changed molecule.
	/// @return Returns total molecular energy contribution, without
	///   intramolecular energy.
	Real calcBatchEnergy(ParallelBox *box, int numMols, int molIdx);
	
	/// Given an array of calculated energies on the device, get
	///   the total energy and return it on the host.
	/// @param box  A ParallelBox containing the molecule data.
	/// @param validEnergies The number of energies over which to
	///   aggregate.
	/// @return Returns the total energy as calculated from the
	///   device energies array.
	Real getEnergyFromDevice(ParallelBox *box, int validEnergies);
	
	/// Each thread of this kernel checks the distance between
	///   the chosen molecule (constant across threads) and one
	///   other molecule. If within cutoff, store index of other
	///   molecule in the 'inCutoff' array. This array will be
	///   compacted to form the molecule batch for the energy
	///   calculations.
	/// @param molecules Device pointer to MoleculeData struct.
	/// @param atoms Device pointer to AtomData struct.
	/// @param currentMol The index of the current changed molecule.
	/// @param startIdx The optional starting index for other molecules.
	///   Used for system energy calculation.
	/// @param enviro Device pointer to Environment struct.
	/// @param inCutoff Device pointer to valid neighbor molecules
	///   array.
	Real calcIntramolEnergy_NLC(Environment *enviro, Molecule *molecules);
    __device__ Real calcAtomDist(Atom atom1, Atom atom2, Environment *enviro);
    Real calcAtomDist(Atom atom1, Atom atom2, Environment *enviro);
    __device__ Real calc_lj(Atom atom1, Atom atom2, Real r2);
    Real calc_lj(Atom atom1, Atom atom2, Real r2);
    __device__ Real calcInterMolecularEnergy(Molecule *molecules, int mol1, int mol2, Environment *enviro);
    __global__ void calcEnergy_NLC(Molecule *molecules, Environment *enviro, int *head, int *lscl, Real *totalEnergy);
	__global__ void checkMoleculeDistances(MoleculeData *molecules, AtomData *atoms, int currentMol, int startIdx, Environment *enviro, int *inCutoff);
	
	/// Each thread in this kernel calculates the inter-atomic
	///   energy between one pair of atoms in the molecule pair
	///   (the chosen molecule, a neighbor molecule from the batch).
	/// @param molecules Device pointer to MoleculeData struct.
	/// @param atoms Device pointer to AtomData struct.
	/// @param currentMol The index of the current changed molecule.
	/// @param enviro Device pointer to Environment struct.
	/// @param energies Device pointer to energies array.
	/// @param numEnergies The maximum index to be used in the
	///   energies array.
	/// @param molBatch Device pointer to list of valid neighbor
	///   molecule indexes.
	/// @param maxMolSize The size (in Atoms) of the largest molecule.
	///   Used for energy segmentation size calculation.
	__global__ void calcInterAtomicEnergy(MoleculeData *molecules, AtomData *atoms, int curentMol, Environment *enviro, Real *energies, int numEnergies, int *molBatch, int maxMolSize);
	
	/// This kernel performs parallel energy aggregation,
	///   and also performs the service of resetting energies
	///   after they have been aggregated. For example, after
	///   any given aggregation pass, the sum of all elements
	///   in the energies array will not have changed, but there
	///   will be fewer, larger energies left in the array.
	/// @param energies Device pointer to energies array.
	/// @param numEnergies The maximum index to be used in the
	///   energies array.
	/// @param interval The distance between each consecutive
	///   energy to be aggregated.
	/// @param batchSize The number of energies to aggregate
	///   per thread.
	/// @note After energy aggregation, all energies that were
	///   found will have been reset to 0.
	__global__ void aggregateEnergies(Real *energies, int numEnergies, int interval, int batchSize);
	
	/// Calculates the LJ energy between two atoms.
	/// @param atoms Device pointer to AtomData struct.
	/// @param atom1 The index of the first atom.
	/// @param atom2 The index of the second atom.
	/// @param r2 The distance between the two atoms, squared.
	/// @return Returns the LJ energy between the two specified
	///   atoms.
	__device__ Real calc_lj(AtomData *atoms, int atom1, int atom2, Real r2);
	
	/// Calculates the charge energy between two atoms.
	/// @param charge1 The charge of atom 1.
	/// @param charge2 The charge of atom 2.
	/// @param r The distance between the two atoms.
	/// @returns Returns the charge energy between two atoms.
	__device__ Real calcCharge(Real charge1, Real charge2, Real r);
	Real calcCharge(Real charge1, Real charge2, Real r);
	
	/// Makes a distance periodic within a specified range.
	/// @param x The distance to be made periodic.
	/// @param boxDim The magnitude of the periodic range.
	/// @return Returns the periodic distance.
	__device__ Real makePeriodic(Real x, Real boxDim);
	
	/// Calculates the geometric mean of two values.
	/// @param d1 The first value.
	/// @param d2 The second value.
	/// @return Returns the geometric mean of the two supplied
	///   values.
	__device__ Real calcBlending(Real d1, Real d2);
	
	/// This method, combined with getYFromIndex(),
	///   provides functionality to generate a non-overlapping
	///   sequence of unique index pairs. For example,
	///   for (int i = 0; i < 6; i++) {
	///      int x = getXFromIndex(i);
	///      int y = getYFromIndex(x, i);
	///      cout << i << ' ' << x << ',' << y << endl; }
	///   Will print the following:
	///   0 1,0
	///   1 2,0
	///   2 2,1
	///   3 3,0
	///   4 3,1
	///   5 3,2
	///   This is used to generate a unique pair of indices
	///   from a unique global pair id.
	/// @param idx The global pair id (index in result array).
	/// @return Returns the X portion of the unique pair.
	__device__ int getXFromIndex(int idx);
	
	/// See documentation for getXFromIndex().
	/// @param x The result of getXFromIndex(idx).
	/// @param idx The global pair id (index in result array).
	__device__ int getYFromIndex(int x, int idx);
}

#endif