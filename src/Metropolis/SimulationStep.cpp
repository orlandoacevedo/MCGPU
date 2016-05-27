/**
 * Base class for an iteration of the simulation.
 */

#include "SimBox.h"
#include "SimulationStep.h"


/** Returns the index of a random molecule within the simulation box */
int SimulationStep::chooseMolecule(SimBox *box) {
  return (int) randomReal(0, box->numMolecules);
}


/** Determines the energy contribution of a particular molecule */
Real SimulationStep::calcMolecularEnergyContribution(Real &subLJ,
                                                     Real &subCharge,
                                                     int currMol,
                                                     int startMol,
                                                     SimBox *box) {
  // TODO: Replace with actual implementation
  return box->calcMolecularEnergyContribution(subLJ, subCharge, currMol,
                                              startMol);
}


/** Perturb a given molecule */
void SimulationStep::changeMolecule(int molIdx, SimBox *box) {
  // TODO: Replace with actual implementation
  box->changeMolecule(molIdx);
}


/** Move a molecule back to its original position */
void SimulationStep::rollback(int molIdx, SimBox *box) {
  // TODO: Replace with actual implementation
  box->rollback(molIdx);
}


/** Determines the total energy of the box */
Real SimulationStep::calcSystemEnergy(Real &subLJ, Real &subCharge,
                                      SimBox *box) {
  // TODO: Replace witha actual implementation
  return box->calcSystemEnergy(subLJ, subCharge);
}

