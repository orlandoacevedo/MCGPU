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
  #pragma acc routine seq
  bool moleculesInRange(int p1Start, int p1End, int p2Start, int p2End,
                        Real** atomCoords, Real* bSize, int* primaryIndexes,
                        Real cutoff);

  /**
   * Calculates the square of the distance between two atoms.
   *
   * @param a1 The atom index of the first atom.
   * @param a2 The atom index of the second atom.
   * @return The square of the distance between the atoms.
   */
  #pragma acc routine seq
  Real calcAtomDistSquared(int a1, int a2, Real** aCoords,
                           Real* bSize);

  /**
   * Given a distance, makes the distance periodic to mitigate distances
   * greater than half the length of the box.
   *
   * @param x The measurement to make periodic.
   * @param dimension The dimension the measurement must be made periodic in.
   * @return The measurement, scaled to be periodic.
   */
  #pragma acc routine seq
  Real makePeriodic(Real x, int dimension, Real* bSize);

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

  /**
   * Calculates the Lennard - Jones potential between two atoms.
   *
   * @param a1  The atom index of the first atom.
   * @param a2  The atom index of the second atom.
   * @param r2  The distance between the atoms, squared.
   * @return The Lennard - Jones potential from the two atoms' interaction.
   */
  #pragma acc routine seq
  Real calcLJEnergy(int a1, int a2, Real r2, Real** aData);

  /**
   * Calculates the Coloumb potential between two atoms.
   *
   * @param a1 The atom index of the first atom.
   * @param a2 The atom index of the second atom.
   * @param r  The distance between the atoms.
   * @return The Coloumb potential from the two atoms' interaction.
   */
  #pragma acc routine seq
  Real calcChargeEnergy(int a1, int a2, Real r, Real** aData);

  /**
   * Calculates the geometric mean of two real numbers.
   *
   * @param a The first real number.
   * @param b The second real number.
   * @return sqrt(|a*b|), the geometric mean of the two numbers.
   */
  #pragma acc routine seq
  Real calcBlending (Real a, Real b);
}

#endif
