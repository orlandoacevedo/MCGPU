/**
 * Base class for an iteration of the simulation.
 *
 * Contains logic common to progressing the simulation, regardless of the
 * particular strategy employed by the configuration.
 */

#ifndef METROPOLIS_SIMULATIONSTEP_H
#define METROPOLIS_SIMULATIONSTEP_H

#include "SimBox.h"
#include "Metropolis/Utilities/MathLibrary.h"

#define VERBOSE true
#define ENABLE_INTRA false
#define ENABLE_BOND 1
#define ENABLE_ANGLE 1
#define ENABLE_DIHEDRAL 0

#define ENABLE_TUNING true
#define RATIO_MARGIN 0.0001
#define TARGET_RATIO 0.4

class SimulationStep {
 public:
  /** Construct a new SimulationStep object from a SimBox pointer */
  SimulationStep(SimBox *box);

  /**
   * Returns the index of a random molecule within the simulation box.
   * @return A random integer from 0 to (numMolecules - 1)
   */
  int chooseMolecule(SimBox *box);

  /**
   * Determines the total energy of a particular molecule in the system.
   *
   * Combines the intramolecular energy with the intermolecular energy.
   * @param currMol The index of the molecule to calculate the contribution of
   * @param startMol The index of the molecule to begin searching from to
   * determine interaction energies (for intermolecular forces)
   * @return the sum of the intermolecular and intramolecular energies
   */
  Real calcMoleculeEnergy(int currMol, int startMol);

  /**
   * Determines the energy contribution of a particular molecule.
   *
   * @param currMol The index of the molecule to calculate the contribution of.
   * @param startMol The index of the molecule to begin searching from to
   *        determine interaction energies.
   * @return The total energy of the box (discounts initial lj / charge energy)
   */
  virtual Real calcMolecularEnergyContribution(int currMol, int startMol) = 0;


  /**
   * Calculates the energy contribution from intramolecular forces within the
   *     given molecule.
   *
   * @param molIdx The index of the molecule to calculate the intramolecular
   *      energy contribution of.
   * @return The intramolecular energy contribution of the molecule.
   */
  Real calcIntraMolecularEnergy(int molIdx);


  /**
   * Randomly moves the molecule within the box. This performs a translation
   * and a rotation, and saves the old position.
   *
   * @param molIdx The index of the molecule to change.
   */
  virtual void changeMolecule(int molIdx, SimBox *box);


  /**
   * Move a molecule back to its original poisition. Performs translation and
   * rotation.
   *
   * @param molIdx The index of the molecule to change.
   */
  virtual void rollback(int molIdx, SimBox *box);


  /**
   * Determines the total energy of the box
   *
   * @param subLJ Initial Lennard - Jones energy.
   * @param subCharge Initial Coulomb energy.
   * @return The total energy of the box.
   */
  virtual Real calcSystemEnergy(Real &subLJ, Real &subCharge, int numMolecules);
};

/**
 * SimCalcs namespace
 *
 * Contains logic for caclulations consumed by the SimulationStep class.
 * Although logically related to the SimulationStep class, must be separated
 * out to a namespace to acurrately run on the GPU.
 */
namespace SimCalcs {
  extern SimBox* sb;
  extern int on_gpu;

  /**
   * Calculates the energy contribution from intramolecular forces within the
   *     given molecule.
   *
   * @param molIdx The index of the molecule to calculate the intramolecular
   *      energy contribution of.
   * @return The intramolecular energy contribution of the molecule.
   */
  Real calcIntraMolecularEnergy(int molIdx);

  /**
   * Calculates the energy from various flexible angles within the molecule.
   *
   * @param molIdx The index of the molecule to calculate the energy of.
   * @return The energy produced by angles being different from their
   *     equilibrium measurement.
   */
  Real angleEnergy(int molIdx);

  /**
   * Expands or contracts the given angle in the given molecule.
   *
   * @param molIdx The index of the molecule that the angle is in.
   * @param angleIdx The index of teh angle within the molecule. For example,
   *     the first angle in molecule 12 has index 0.
   * @param expandDeg The amount to expand or contract the angle. Positive
   *     values correspond to expansion, negative to contraction. This is
   *     measured in degrees.
   */
  void expandAngle(int molIdx, int angleIdx, Real expandDeg);

  /**
   * Calculates the energy from various flexible bonds within the molecule.
   *
   * @param molIdx The index of the molecule to calculate the energy of.
   * @return The energy produced by bonds being different from their
   *     equilibrium measurement.
   */
   Real bondEnergy(int molIdx);

  /**
   * Stretches or compresses the given bond in the given molecule.
   *
   * @param molIdx The index of the molecule that the bond is in.
   * @param bondIdx The index of the bond within the molecule. For example, the
   *     first bond in molecule 12 has index 0.
   * @param stretchDist The amount to stretch or compress the bond. Positive
   *     values correspond to stretching, negative to compression.
   */
  void stretchBond(int molIdx, int bondIdx, Real stretchDist);

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

  /**
   * Performs all the movements of a molecule for a simulation step. This
   * includes intermolecular movements (rotation and translation) as well
   * as a random selection of intramolecular movements (bond stretching,
   * angle expanding, and dihedral movements.
   *
   * @param molIdx The index of the molecule to change.
   */
  void changeMolecule(int molIdx);


  bool acceptVolumeMove(Real deltaE, Real oldVolume, Real pressure);

  void resizeBox(Real factor);

  Real calcDirectEwaldSum(int radius);

  /**
   * Given a molecule to move, randomly rotates and translates the molecule.
   * Saves the old position in case of a rollback.
   *
   * @param molIdx The index of the molecule to change.
   */
  void intermolecularMove(int molIdx);

  /**
   * Given a molecule, perform a random series of bond stretches, angle
   * expansions, an dihedral movements on a molecule. The number of
   * translations will be a random number between 1 and the user-defined
   * maximum.
   *
   * @param molIdx The index of the molecule to peform movements on
   */
  void intramolecularMove(int molIdx);

  /**
   * Keeps a molecule in the simulation box's boundaries, based on the location
   * of its first primary index atom.
   *
   * @param molIdx The index of the molecule to keep inside the box.
   */
  #pragma acc routine seq
  void keepMoleculeInBox(int molIdx, Real** aCoords, int** molData,
                         int* pIdxes, Real* bsize);


  /**
   * Given an atom to translate, and the directions to tranlate it, moves it.
   *
   * @param aIdx The index of the atom to change.
   * @param dX The amount to translate in the x direction.
   * @param dY The amount to translate in the y direction.
   * @param dZ The amount to tranlsate in the z direction.
   */
  #pragma acc routine seq
  void translateAtom(int aIdx, Real dX, Real dY, Real dZ, Real** aCoords);

  /**
   * Given an atom to rotate, its pivot point, and the amounts to rotate,
   * rotates the atom around the pivot (rotations in degrees).
   *
   * @param aIdx The index of the atom to change.
   * @param pivotIdx The index of the atom about which to rotate.
   * @param rotX The amount of rotation around the x-axis that is done.
   * @param rotY The amount of rotation around the y-axis that is done.
   * @param rotZ The amount of rotation around the z-axis that is done.
   */
  #pragma acc routine seq
  void rotateAtom(int aIdx, int pivotIdx, Real rotX, Real rotY, Real rotZ,
                  Real** aCoords);

  /**
   * Given an atom and an amount to rotate, rotates about the x-axis.
   *
   * @param aIdx The index of the atom to rotate.
   * @param angleDeg The angle to rotate it (in degrees).
   */
  #pragma acc routine seq
  void rotateX(int aIdx, Real angleDeg, Real** aCoords);

  /**
   * Given an atom and an amount to rotate, rotates it about the y-axis.
   *
   * @param aIdx The index of the atom to rotate.
   * @param angleDeg The angle to rotate it (in degrees).
   */
  #pragma acc routine seq
  void rotateY(int aIdx, Real angleDeg, Real** aCoords);

  /**
   * Given an atom and an amount to rotate, rotates it about the z-axis.
   *
   * @param aIdx The index of the atom to rotate.
   * @param angleDeg The angle to rotate it (in degrees).
   */
  #pragma acc routine seq
  void rotateZ(int aIdx, Real angleDeg, Real** aCoords);

  /**
   * Roll back a molecule to its original poisition. Performs translation and
   * rotation.
   *
   * @param molIdx The index of the molecule to change.
   */
  void rollback(int molIdx);

  /**
   * Rolls back the size of the most recently changed angle. Does not change
   *     any angle coordinates.
   */
  void rollbackAngles(int molIdx);

  /**
   * Rolls back the length of the most recently changed bond. Does not change
   *     any atom coordinates.
   */
  void rollbackBonds(int molIdx);

  /**
   * Saves the bonds of a particular molecule in the event of a rollback
   */
  void saveBonds(int molIdx);

  /**
   * Saves the angles of a particular molecule in the event of a rollback
   */
  void saveAngles(int molIdx);

  /**
   * Unions the sets of two different atoms connected by a bond. Note that atoms
   *     are indexed within the molecule, so the first atom in the 12th molecule
   *     has index 0.
   *
   * @param atom1 The index (within the molecule) of the first atom to union.
   * @param atom2 The index (within the molecule) second atom in the union.
   */
  void unionAtoms(int atom1, int atom2);

  /**
   * Finds the representative element of an atom in the disjoint set array.
   * Note that the atom is indexed within the molecule, so the first atom in the
   * 12th molecule has index 0. This function also performs path compression.
   *
   * @param atomIdx The atom's index within the molecule.
   * @return The representative element of the atom index.
   */
  int find (int atomIdx);

  /**
   * Determines if a move should be accepted based on the before and after
   * energy values.
   *
   * @param oldEnergy the total molecular energy prior to the move
   * @param newEnergy the total molecular energy proceeding the move
   */
  bool acceptMove(Real deltaE);

  /** Set the current SimBox instance for this namespace */
  void setSB(SimBox* sb);
}

#endif
