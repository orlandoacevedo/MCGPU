/**
 * VerletStep.h
 *
 * A subclass of SimulationStep that uses a "neighbor linked cell" for energy
 * calculations
 *
 *
 * In a neighbor linked cell implementation, the simulation box is divided into
 * smaller cells, each of which has dimension >= the cutoff radius. When looking
 * at molecules in range of a given molecule, we only examine the current cell
 * and all other adjacent cells.
 *
 * As molecules move from one cell to another, the NLC structure is automatically
 * updated.
 *
 * (Future Work: Rework from linked list to array for more parallelism).
 */

#ifndef METROPOLIS_VERLET_STEP_H
#define METROPOLIS_VERLET_STEP_H

#include "SimulationStep.h"

struct Verlet {
  Real** startingLocation;
  int**  neighbors;
  int*   nNeighbors;
  Real   shellRadius;
  Real   maxMove;
};

class VerletStep: public SimulationStep {
 public:
  explicit VerletStep(SimBox* box): SimulationStep(box),
                                             verlet(NULL) {}
  virtual ~VerletStep();
  virtual Real calcSystemEnergy(Real &subLJ, Real &subCharge,
                                int numMolecules);
  virtual Real calcMolecularEnergyContribution(int currMol, int startMol);
  virtual void changeMolecule(int molIdx, SimBox *box);
  virtual void rollback(int molIdx, SimBox *box);

 private:
  Verlet* verlet;
};

namespace VerletCalcs {

  Real calcMolecularEnergyContribution(int currMol, int startMol,
                                       Verlet *verlet);

  Real calcMoleculeInteractionEnergy (int m1, int m2, int** molData,
                                      Real** aData, Real** aCoords,
                                      Real* bSize);

  Verlet* createVerlet();

  void updateVerlet(int molIdx, Verlet *verlet);

  void freeVerlet(Verlet * verlet);

  void fillVerlet(Verlet * verlet);
}

#endif
