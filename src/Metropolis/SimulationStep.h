/**
 * Base class for an iteration of the simulation.
 *
 * Contains logic common to progressing the simulation, regardless of the
 * particular strategy employed by the configuration.
 */

#include "SimBox.h"
#include "Metropolis/Utilities/MathLibrary.h"

class SimulationStep {
 public:
  /**
   * Returns the index of a random molecule within the simulation box.
   * @return A random integer from 0 to (numMolecules - 1)
   */
  int chooseMolecule(SimBox *box);


  /**
   * Determines the energy contribution of a particular molecule.
   *
   * @param subLJ Initial Lennard - Jones energy. Final L-J energy passed out
   *        by reference.
   * @param subCharge Initial Coulomb energy. Final Coulomb energy passed out
   *        by reference.
   * @param currMol The index of the molecule to calculate the contribution of.
   * @param startMol The index of the molecule to begin searching from to
   *        determine interaction energies.
   * @return The total energy of the box (discounts initial lj / charge energy)
   */
  Real calcMolecularEnergyContribution(Real &subLJ, Real &subCharge,
                                       int currMol, int startMol,
                                       SimBox *box);


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
   * @param subLJ Initial Lennard - Jones energy. Final L-J energy passed out by reference.
   * @param subCharge Initial Coulomb energy. Final Coulomb energy passed out by reference.
   * @return The total energy of the box.
   */
  Real calcSystemEnergy(Real &subLJ, Real &subCharge, SimBox *box);
};
