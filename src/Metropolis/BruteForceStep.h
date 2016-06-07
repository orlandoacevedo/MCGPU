/**
 * BruteForceStep.h
 *
 * A subclass of SimulationStep that uses a "brute force" strategy for energy
 * calculations
 */

#ifndef METROPOLIS_BRUTEFORCE_H
#define METROPOLIS_BRUTEFORCE_H

#include "SimulationStep.h"

class BruteForceStep: public SimulationStep {
 public:
  /** Construct a BruteForceStep object from a SimBox pointer */
  BruteForceStep(SimBox* box): SimulationStep(box) {}

  /**
   * Determines the energy contribution of a particular molecule.
   *
   * @param currMol The index of the molecule to calculate the contribution of
   * @param startMol The index of the molecule to begin searching from to
   *        determine interaction energies.
   * @return The total energy of the box (discounts initial lj /
   * charge energy)
   */
  virtual Real calcMolecularEnergyContribution(int currMol, int startMol);
};


/**
 * BruteForceCalcs namespace
 *
 * Contains logic for caclulations consumed by the BruteForceStep class.
 * Although logically related to the BruteForceStep class, must be separated
 * out to a namespace to acurrately run on the GPU.
 */
namespace BruteForceCalcs {
  /**
   * Determines the energy contribution of a particular molecule.
   *
   * @param currMol The index of the molecule to calculate the contribution of
   * @param startMol The index of the molecule to begin searching from to
   * determine interaction energies.
   * @return The total energy of the box (discounts initial lj &
   * charge energy)
   */
  Real calcMolecularEnergyContribution(int currMol, int startMol);

  /**
   * Determines whether or not two molecule's primaryIndexes are within the
   * cutoff range of one another.
   *
   * @param p1Start The index of the first primary index from molecule 1.
   * @param p1End index of the last primary index from molecule 1,  + 1.
   * @param p2Start index of the first primary index from molecule 2.
   * @param p2End index of the last primary index from molecule 2,  + 1.
   * @param atomCoords The coordinates of the atoms to check.
   * @return true if the molecules are in range, false otherwise.
   */
  #pragma acc routine vector
  Real calcMoleculeInteractionEnergy (int m1, int m2, int** molData,
                                      Real** aData, Real** aCoords,
                                      Real* bSize);
}

#endif
