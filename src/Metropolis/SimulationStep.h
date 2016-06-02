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
   * Determines the energy contribution of a particular molecule.
   *
   * @param currMol The index of the molecule to calculate the contribution of.
   * @param startMol The index of the molecule to begin searching from to
   *        determine interaction energies.
   * @return The total energy of the box (discounts initial lj / charge energy)
   */
  virtual Real calcMolecularEnergyContribution(int currMol, int startMol) = 0;


  /**
   * Randomly moves the molecule within the box. This performs a translation
   * and a rotation, and saves the old position.
   *
   * @param molIdx The index of the molecule to change.
   */
  void changeMolecule(int molIdx, SimBox *box);


  /**
   * Move a molecule back to its original poisition. Performs translation and
   * rotation.
   *
   * @param molIdx The index of the molecule to change.
   */
  void rollback(int molIdx, SimBox *box);


  /**
   * Determines the total energy of the box
   *
   * @param subLJ Initial Lennard - Jones energy.
   * @param subCharge Initial Coulomb energy.
   * @return The total energy of the box.
   */
  Real calcSystemEnergy(Real &subLJ, Real &subCharge, int numMolecules);
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
   * Given a molecule to change, randomly moves the molecule within the box.
   * This performs a translation and a rotation, and saves the old position.
   *
   * @param molIdx The index of the molecule to change.
   */
  void changeMolecule(int molIdx);

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

  /** Set the current SimBox instance for this namespace */
  void setSB(SimBox* sb);
}

#endif
