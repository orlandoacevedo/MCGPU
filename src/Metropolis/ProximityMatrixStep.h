/**
 * ProximityMatrixStep.h
 *
 * A subclass of SimulationStep that uses a "proximity matrix" for energy calculations
 *
 * We define a "proximity matrix" analogously to an adjacency matrix in graph theory:
 * There is one row and one column for each molecule, and the entry in row i column j
 * is 1 if moleculesInRange() indicates that molecules i and j are in range of each
 * other (and 0 if those molecules are not in range).
 *
 * This matrix is symmetric, since if entry (i,j) indicates that molecules i and j are
 * in range, entry (j,i) should also indicate the same.
 *
 * For GPU efficiency, the matrix is linearized as a one-dimensional char array, where
 * entry (i,j) is stored in proximityMatrix[i*numMolecules + j].
 *
 * This is not an especially space-efficient data structure, but it can be updated
 * quickly (in parallel) on the GPU.
 *
 * (Future Work: This could be stored as a bit matrix to save space...)
 */

#ifndef METROPOLIS_PROXIMITYMATRIX_H
#define METROPOLIS_PROXIMITYMATRIX_H

#include "SimulationStep.h"

class ProximityMatrixStep: public SimulationStep {
 public:
  explicit ProximityMatrixStep(SimBox* box): SimulationStep(box), proximityMatrix(NULL) {}
  virtual ~ProximityMatrixStep();
  virtual Real calcSystemEnergy(Real &subLJ, Real &subCharge, int numMolecules);
  virtual Real calcMolecularEnergyContribution(int currMol, int startMol);
  virtual void changeMolecule(int molIdx, SimBox *box);
  virtual void rollback(int molIdx, SimBox *box);
private:
  char *proximityMatrix;
};

namespace ProximityMatrixCalcs {

  Real calcMolecularEnergyContribution(int currMol, int startMol, char *proximityMatrix);

  #pragma acc routine seq
  bool moleculesInRange(int p1Start, int p1End, int p2Start, int p2End,
                        Real** atomCoords, Real* bSize, int* primaryIndexes,
                        Real cutoff);

  #pragma acc routine seq
  Real calcAtomDistSquared(int a1, int a2, Real** aCoords,
                           Real* bSize);

  #pragma acc routine seq
  Real makePeriodic(Real x, int dimension, Real* bSize);

  #pragma acc routine vector
  Real calcMoleculeInteractionEnergy (int m1, int m2, int** molData,
                                      Real** aData, Real** aCoords,
                                      Real* bSize);

  #pragma acc routine seq
  Real calcLJEnergy(int a1, int a2, Real r2, Real** aData);

  #pragma acc routine seq
  Real calcChargeEnergy(int a1, int a2, Real r, Real** aData);

  #pragma acc routine seq
  Real calcBlending (Real a, Real b);

  char *createProximityMatrix();

  void updateProximityMatrix(char *matrix, int i);

  void freeProximityMatrix(char *matrix);
}

#endif
