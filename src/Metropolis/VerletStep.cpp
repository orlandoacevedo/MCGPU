#include "BruteForceStep.h"
#include "VerletStep.h"
#include "SimulationStep.h"
#include "GPUCopy.h"

#ifdef _OPENACC
#include <openacc.h>
#endif


VerletStep::~VerletStep() {
  VerletCalcs::freeVerlet(this->verlet);
  this->verlet = NULL;
}

Real VerletStep::calcMolecularEnergyContribution(int currMol, int startMol) {
  return VerletCalcs::calcMolecularEnergyContribution( currMol,
      startMol, this->verlet);
}

Real VerletStep::calcSystemEnergy(Real &subLJ, Real &subCharge, int numMolecules) {
  Real result = SimulationStep::calcSystemEnergy(subLJ, subCharge,
                                                 numMolecules);
  this->verlet = VerletCalcs::createVerlet();
  return result;
}

void VerletStep::changeMolecule(int molIdx, SimBox *box) {
  SimulationStep::changeMolecule(molIdx, box);
  VerletCalcs::updateVerlet(molIdx, this->verlet);
}

void VerletStep::rollback(int molIdx, SimBox *box) {
  SimulationStep::rollback(molIdx, box);
  VerletCalcs::updateVerlet(molIdx, this->verlet);
}

// ----- VerletCalcs Definitions -----

Real VerletCalcs::calcMolecularEnergyContribution(
    int currMol, int startMol, Verlet* verlet) {

  Real total = 0;

  int **molData = GPUCopy::moleculeDataPtr();
  Real **atomCoords = GPUCopy::atomCoordinatesPtr();
  Real *bSize = GPUCopy::sizePtr();
  int *pIdxes = GPUCopy::primaryIndexesPtr();
  Real **aData = GPUCopy::atomDataPtr();
  Real cutoff = SimCalcs::sb->cutoff;
  const long numMolecules = SimCalcs::sb->numMolecules;

  const int p1Start = SimCalcs::sb->moleculeData[MOL_PIDX_START][currMol];
  const int p1End = (SimCalcs::sb->moleculeData[MOL_PIDX_COUNT][currMol]
                     + p1Start);

  if (verlet == NULL) {
    for (int otherMol = startMol; otherMol < numMolecules; otherMol++) {
      if (otherMol != currMol) {
        int p2Start = molData[MOL_PIDX_START][otherMol];
        int p2End = molData[MOL_PIDX_COUNT][otherMol] + p2Start;
        if (SimCalcs::moleculesInRange(p1Start, p1End, p2Start, p2End,
                                       atomCoords, bSize, pIdxes, cutoff)) {
          total += calcMoleculeInteractionEnergy(currMol, otherMol, molData,
                                                 aData, atomCoords, bSize);
        }
      }
    }
  } else {
    for (int i = 0; i < verlet->nNeighbors[startMol]; i++) {
      int otherMol = verlet->neighbors[startMol][i];
      if (otherMol != currMol) {
        int p2Start = molData[MOL_PIDX_START][otherMol];
        int p2End = molData[MOL_PIDX_COUNT][otherMol] + p2Start;
        if (SimCalcs::moleculesInRange(p1Start, p1End, p2Start, p2End,
                                       atomCoords, bSize, pIdxes, cutoff)) {
          total += calcMoleculeInteractionEnergy(currMol, otherMol, molData,
                                                 aData, atomCoords, bSize);
        }
      }
    }
  }
  return total;
}

// TODO: Duplicate; abstract out when PGCC supports it
Real VerletCalcs::calcMoleculeInteractionEnergy (int m1, int m2,
                                                          int** molData,
                                                          Real** aData,
                                                          Real** aCoords,
                                                          Real* bSize) {
  Real energySum = 0;

  const int m1Start = molData[MOL_START][m1];
  const int m1End = molData[MOL_LEN][m1] + m1Start;

  const int m2Start = molData[MOL_START][m2];
  const int m2End = molData[MOL_LEN][m2] + m2Start;

  for (int i = m1Start; i < m1End; i++) {
    for (int j = m2Start; j < m2End; j++) {
      if (aData[ATOM_SIGMA][i] >= 0 && aData[ATOM_SIGMA][j] >= 0
          && aData[ATOM_EPSILON][i] >= 0 && aData[ATOM_EPSILON][j] >= 0) {

        const Real r2 = SimCalcs::calcAtomDistSquared(i, j, aCoords, bSize);
        if (r2 == 0.0) {
          energySum += 0.0;
        } else {
          energySum += SimCalcs::calcLJEnergy(i, j, r2, aData);
          energySum += SimCalcs::calcChargeEnergy(i, j, sqrt(r2), aData);
        }
      }
    }
  }

  return (energySum);
}

Verlet* VerletCalcs::createVerlet() {

  const Real cutoff = SimCalcs::sb->cutoff;
  const int numMolecules = SimCalcs::sb->numMolecules;

  Real* bSize = GPUCopy::sizePtr();
  int* pIdxes = GPUCopy::primaryIndexesPtr();
  int** molData = GPUCopy::moleculeDataPtr();
  Real** atomCoords = GPUCopy::atomCoordinatesPtr();


  Verlet* simulation_Verlet = new Verlet;
  // Shell width may need to be tuned. This value is taken from the
  // Spring 2016 Senior Design Group.
  simulation_Verlet->shellRadius = 1.2875 * cutoff;
  simulation_Verlet->maxMove = (simulation_Verlet->shellRadius - cutoff) / 2.0;
  simulation_Verlet->startingLocation = new Real*[NUM_DIMENSIONS];
  simulation_Verlet->neighbors = new int*[numMolecules];
  simulation_Verlet->nNeighbors = new int[numMolecules];

  for (int dim = 0; dim < NUM_DIMENSIONS; dim++) {
    simulation_Verlet->startingLocation[dim] = new Real[numMolecules];
  }

  for (int i = 0; i < numMolecules; i++) {
    simulation_Verlet->neighbors[i] = new int[numMolecules];
  }

  fillVerlet(simulation_Verlet);

  return simulation_Verlet;
}

void VerletCalcs::fillVerlet(Verlet * verlet) {



  int numMolecules = SimCalcs::sb->numMolecules;

  Real* bSize = GPUCopy::sizePtr();
  int* pIdxes = GPUCopy::primaryIndexesPtr();
  int** molData = GPUCopy::moleculeDataPtr();
  Real** atomCoords = GPUCopy::atomCoordinatesPtr();

  for (int dim = 0; dim < NUM_DIMENSIONS; dim++) {
    for (int mol = 0; mol < numMolecules; mol++) {
      int atom = pIdxes[molData[MOL_PIDX_START][mol]];
      verlet->startingLocation[dim][mol] = atomCoords[dim][atom];
    }
  }


  Real sr2 = verlet->shellRadius * verlet->shellRadius;

  for (int i = 0; i < numMolecules; i++) {
    int count = 0;
    int a1 = pIdxes[molData[MOL_PIDX_START][i]];
    for (int j = 0; j < numMolecules; j++) {
      if (i == j) continue;
      int a2 = pIdxes[molData[MOL_PIDX_START][j]];
      Real r2 = SimCalcs::calcAtomDistSquared(a1, a2, atomCoords, bSize);
      if (sr2 <= r2) {
        verlet->neighbors[i][count++] = j;
      }
    }
    verlet->nNeighbors[i] = count;
  }

}

void VerletCalcs::updateVerlet(int molIdx, Verlet * v) {


  int* pIdxes = GPUCopy::primaryIndexesPtr();
  int** molData = GPUCopy::moleculeDataPtr();
  Real** atomCoords = GPUCopy::atomCoordinatesPtr();

  int atom = pIdxes[molData[MOL_PIDX_START][molIdx]];
  Real dx = atomCoords[X_COORD][atom] - v->startingLocation[X_COORD][molIdx];
  Real dy = atomCoords[Y_COORD][atom] - v->startingLocation[Y_COORD][molIdx];
  Real dz = atomCoords[Z_COORD][atom] - v->startingLocation[Z_COORD][molIdx];
  Real dist = dx * dx + dy * dy + dz * dz;
  Real max = v->maxMove * v->maxMove;
  if (dist > max) {
    fillVerlet(v);
  }
}

void VerletCalcs::freeVerlet(Verlet *verlet) {

  int numMolecules = SimCalcs::sb->numMolecules;

  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    free(verlet->startingLocation[i]);
  }
  for (int i = 0; i < numMolecules; i++) {
    free(verlet->neighbors[i]);
  }
  free(verlet->neighbors);
  free(verlet->nNeighbors);
  free(verlet->startingLocation);
  free(verlet);
}
