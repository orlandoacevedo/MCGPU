#include "BruteForceStep.h"
#include "ProximityMatrixStep.h"
#include "SimulationStep.h"
#include "GPUCopy.h"
#include <openacc.h>


ProximityMatrixStep::~ProximityMatrixStep() {
  ProximityMatrixCalcs::freeProximityMatrix(this->proximityMatrix);
  this->proximityMatrix = NULL;
}

Real ProximityMatrixStep::calcMolecularEnergyContribution(int currMol,
                                                     int startMol) {
  return ProximityMatrixCalcs::calcMolecularEnergyContribution(currMol, startMol, this->proximityMatrix);
}

Real ProximityMatrixStep::calcSystemEnergy(Real &subLJ, Real &subCharge, int numMolecules) {
  Real result = SimulationStep::calcSystemEnergy(subLJ, subCharge, numMolecules);
  this->proximityMatrix = ProximityMatrixCalcs::createProximityMatrix();
  return result;
}

void ProximityMatrixStep::changeMolecule(int molIdx, SimBox *box) {
  SimulationStep::changeMolecule(molIdx, box);
  ProximityMatrixCalcs::updateProximityMatrix(this->proximityMatrix, molIdx);
}

void ProximityMatrixStep::rollback(int molIdx, SimBox *box) {
  SimulationStep::rollback(molIdx, box);
  ProximityMatrixCalcs::updateProximityMatrix(this->proximityMatrix, molIdx);
}

// ----- ProximityMatrixCalcs Definitions -----

Real ProximityMatrixCalcs::calcMolecularEnergyContribution(int currMol,
                                                      int startMol, char *proximityMatrix) {
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

  if (proximityMatrix == NULL) {
    #pragma acc parallel loop gang deviceptr(molData, atomCoords, bSize, \
        pIdxes, aData) if (SimCalcs::on_gpu) vector_length(64)
    for (int otherMol = startMol; otherMol < numMolecules; otherMol++) {
      if (otherMol != currMol) {
        int p2Start = molData[MOL_PIDX_START][otherMol];
        int p2End = molData[MOL_PIDX_COUNT][otherMol] + p2Start;
        if (SimCalcs::moleculesInRange(p1Start, p1End, p2Start, p2End, atomCoords, bSize,
                             pIdxes, cutoff)) {
          total += calcMoleculeInteractionEnergy(currMol, otherMol, molData,
                                                 aData, atomCoords, bSize);
        }
      }
    }
  } else {
    #pragma acc parallel loop gang deviceptr(molData, atomCoords, bSize, \
        pIdxes, aData, proximityMatrix) if (SimCalcs::on_gpu) vector_length(64)
    for (int otherMol = startMol; otherMol < numMolecules; otherMol++) {
      if (otherMol != currMol) {
        //int p2Start = molData[MOL_PIDX_START][otherMol];
        //int p2End = molData[MOL_PIDX_COUNT][otherMol] + p2Start;
        if (proximityMatrix[currMol*numMolecules + otherMol]) {
          total += calcMoleculeInteractionEnergy(currMol, otherMol, molData,
                                                 aData, atomCoords, bSize);
        }
      }
    }
  }

  return total;
}

// FIXME BEGIN FUNCTIONS DUPLICATED FROM BRUTEFORCESTEP...

Real ProximityMatrixCalcs::calcMoleculeInteractionEnergy (int m1, int m2,
                                                     int** molData,
                                                     Real** aData,
                                                     Real** aCoords,
                                                     Real* bSize) {
  Real energySum = 0;

  const int m1Start = molData[MOL_START][m1];
  const int m1End = molData[MOL_LEN][m1] + m1Start;

  const int m2Start = molData[MOL_START][m2];
  const int m2End = molData[MOL_LEN][m2] + m2Start;

  #pragma acc loop vector collapse(2) reduction(+:energySum)
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

// FIXME ...END FUNCTIONS DUPLICATED FROM BRUTEFORCESTEP

char *ProximityMatrixCalcs::createProximityMatrix() {
  const long numMolecules = SimCalcs::sb->numMolecules;
  const Real cutoff = SimCalcs::sb->cutoff;

  int** molData = GPUCopy::moleculeDataPtr();
  Real** atomCoords = GPUCopy::atomCoordinatesPtr();
  Real* bSize = GPUCopy::sizePtr();
  int* pIdxes = GPUCopy::primaryIndexesPtr();
  Real** aData = GPUCopy::atomDataPtr();

  char *matrix;
  if (SimCalcs::on_gpu) {
    matrix = (char *)acc_malloc(numMolecules * numMolecules * sizeof(char));
  } else {
    matrix = (char *)malloc(numMolecules * numMolecules * sizeof(char));
  }
  assert(matrix != NULL);
#pragma acc parallel loop deviceptr(molData, atomCoords, bSize, pIdxes, aData, matrix) if(SimCalcs::on_gpu)
  for (int i = 0; i < numMolecules; i++) {
    const int p1Start = molData[MOL_PIDX_START][i];
    const int p1End   = molData[MOL_PIDX_COUNT][i] + p1Start;
#pragma acc loop seq
    for (int j = 0; j < numMolecules; j++) {
      const int p2Start = molData[MOL_PIDX_START][j];
      const int p2End = molData[MOL_PIDX_COUNT][j] + p2Start;
      matrix[i*numMolecules + j] = (j != i && SimCalcs::moleculesInRange(p1Start, p1End, p2Start, p2End, atomCoords, bSize, pIdxes, cutoff));
    }
  }
  return matrix;
}

void ProximityMatrixCalcs::updateProximityMatrix(char *matrix, int i) {
  const long numMolecules = SimCalcs::sb->numMolecules;
  const Real cutoff = SimCalcs::sb->cutoff;

  int** molData = GPUCopy::moleculeDataPtr();
  Real** atomCoords = GPUCopy::atomCoordinatesPtr();
  Real* bSize = GPUCopy::sizePtr();
  int* pIdxes = GPUCopy::primaryIndexesPtr();
  Real** aData = GPUCopy::atomDataPtr();

#pragma acc parallel loop deviceptr(molData, atomCoords, bSize, pIdxes, aData, matrix) if(SimCalcs::on_gpu)
  for (int j = 0; j < numMolecules; j++) {
    const int p1Start = molData[MOL_PIDX_START][i];
    const int p1End   = molData[MOL_PIDX_COUNT][i] + p1Start;
    const int p2Start = molData[MOL_PIDX_START][j];
    const int p2End = molData[MOL_PIDX_COUNT][j] + p2Start;

    const char entry = (j != i && SimCalcs::moleculesInRange(p1Start, p1End, p2Start, p2End, atomCoords, bSize, pIdxes, cutoff));
    matrix[i*numMolecules + j] = entry;
    matrix[j*numMolecules + i] = entry;
  }
}

void ProximityMatrixCalcs::freeProximityMatrix(char *matrix) {
  if (SimCalcs::on_gpu)
    acc_free(matrix);
  else
    free(matrix);
}
