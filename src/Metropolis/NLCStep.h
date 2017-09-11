/**
 * NLCStep.h
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

#ifndef METROPOLIS_NLC_STEP_H
#define METROPOLIS_NLC_STEP_H

#include "SimulationStep.h"

struct NLC {
  NLC_Node* nlc_heap;

  NLC_Node** neighbors;

  NLC_Node**** neighborCells;

  NLC_Node* prevNode;

  int* numCells;
  Real* cellWidth;
  int* prevCell;
};

class NLCStep: public SimulationStep {
 public:
  explicit NLCStep(SimBox* box): SimulationStep(box),
                                             nlc(NULL) {}
  virtual ~NLCStep();
  virtual Real calcSystemEnergy(Real &subLJ, Real &subCharge,
                                int numMolecules);
  virtual Real calcMolecularEnergyContribution(int currMol, int startMol);
  virtual void changeMolecule(int molIdx, SimBox *box);
  virtual void rollback(int molIdx, SimBox *box);

 private:
  NLC* nlc;
};

namespace NLCCalcs {

  Real calcMolecularEnergyContribution(int currMol, int startMol,
                                       NLC *nlc);

  Real calcMoleculeInteractionEnergy (int m1, int m2, int** molData,
                                      Real** aData, Real** aCoords,
                                      Real* bSize);

  NLC* createNLC();

  void updateNLC(int molIdx, NLC *nlc);

  void freeNLC(NLC * nlc);

  int findNeighbors(int startMol, NLC* nlc);

  int getCell(Real loc, int dim, NLC* nlc);

  int wrapCell(int idx, int dim, NLC* nlc);
}

#endif
