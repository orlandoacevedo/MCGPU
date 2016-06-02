#include "BruteForceStep.h"
#include "SimulationStep.h"
#include "GPUCopy.h"


Real BruteForceStep::calcMolecularEnergyContribution(int currMol,
                                                     int startMol) {
  return BruteForceCalcs::calcMolecularEnergyContribution(currMol, startMol);
}


// ----- BruteForceCalcs Definitions -----


Real BruteForceCalcs::calcMolecularEnergyContribution(int currMol,
                                                      int startMol) {
  Real total = 0;

  int** molData = GPUCopy::moleculeDataPtr();
  Real** atomCoords = GPUCopy::atomCoordinatesPtr();
  Real* bSize = GPUCopy::sizePtr();
  int* pIdxes = GPUCopy::primaryIndexesPtr();
  Real** aData = GPUCopy::atomDataPtr();
  Real cutoff = SimCalcs::sb->cutoff;
  const long numMolecules = SimCalcs::sb->numMolecules;

  const int p1Start = SimCalcs::sb->moleculeData[MOL_PIDX_START][currMol];
  const int p1End = (SimCalcs::sb->moleculeData[MOL_PIDX_COUNT][currMol]
                     + p1Start);

  #pragma acc parallel loop gang deviceptr(molData, atomCoords, bSize, \
      pIdxes, aData) if (SimCalcs::on_gpu) vector_length(64)
  for (int otherMol = startMol; otherMol < numMolecules; otherMol++) {
    if (otherMol != currMol) {
      int p2Start = molData[MOL_PIDX_START][otherMol];
      int p2End = molData[MOL_PIDX_COUNT][otherMol] + p2Start;
      if (moleculesInRange(p1Start, p1End, p2Start, p2End, atomCoords, bSize,
                           pIdxes, cutoff)) {
        total += calcMoleculeInteractionEnergy(currMol, otherMol, molData,
                                               aData, atomCoords, bSize);
      }
    }
  }

  return total;
}

bool BruteForceCalcs::moleculesInRange(int p1Start, int p1End, int p2Start,
                                       int p2End, Real** atomCoords,
                                       Real* bSize, int* primaryIndexes,
                                       Real cutoff) {
  bool out = false;
  for (int p1Idx = p1Start; p1Idx < p1End; p1Idx++) {
    int p1 = primaryIndexes[p1Idx];
    for (int p2Idx = p2Start; p2Idx < p2End; p2Idx++) {
      int p2 = primaryIndexes[p2Idx];
      out |= (calcAtomDistSquared(p1, p2, atomCoords, bSize) <= cutoff * cutoff);
    }
  }
  return out;
}

Real BruteForceCalcs::calcAtomDistSquared(int a1, int a2, Real** aCoords,
                                          Real* bSize) {
  Real dx = makePeriodic(aCoords[X_COORD][a2] - aCoords[X_COORD][a1],
                         X_COORD, bSize);
  Real dy = makePeriodic(aCoords[Y_COORD][a2] - aCoords[Y_COORD][a1],
                         Y_COORD, bSize);
  Real dz = makePeriodic(aCoords[Z_COORD][a2] - aCoords[Z_COORD][a1],
                         Z_COORD, bSize);

  return dx * dx + dy * dy + dz * dz;
}

Real BruteForceCalcs::makePeriodic(Real x, int dimension, Real* bSize) {
  Real dimLength = bSize[dimension];

  int lt = (x < -0.5 * dimLength); // 1 or 0
  x += lt * dimLength;
  int gt = (x > 0.5 * dimLength);  // 1 or 0
  x -= gt * dimLength;
  return x;
}

Real BruteForceCalcs::calcMoleculeInteractionEnergy (int m1, int m2,
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

        const Real r2 = calcAtomDistSquared(i, j, aCoords, bSize);
        if (r2 == 0.0) {
          energySum += 0.0;
        } else {
          energySum += calcLJEnergy(i, j, r2, aData);
          energySum += calcChargeEnergy(i, j, sqrt(r2), aData);
        }
      }
    }
  }

  return (energySum);
}

Real BruteForceCalcs::calcLJEnergy(int a1, int a2, Real r2, Real** aData) {
  if (r2 == 0.0) {
    return 0.0;
  } else {
    const Real sigma = calcBlending(aData[ATOM_SIGMA][a1],
        aData[ATOM_SIGMA][a2]);
    const Real epsilon = calcBlending(aData[ATOM_EPSILON][a1],
        aData[ATOM_EPSILON][a2]);

    const Real s2r2 = pow(sigma, 2) / r2;
    const Real s6r6 = pow(s2r2, 3);
    const Real s12r12 = pow(s6r6, 2);
    return 4.0 * epsilon * (s12r12 - s6r6);
  }
}

Real BruteForceCalcs::calcChargeEnergy(int a1, int a2, Real r, Real** aData) {
  if (r == 0.0) {
    return 0.0;
  } else {
    const Real e = 332.06;
    return (aData[ATOM_CHARGE][a1] * aData[ATOM_CHARGE][a2] * e) / r;
  }
}

Real BruteForceCalcs::calcBlending (Real a, Real b) {
  if (a * b >= 0) {
    return sqrt(a*b);
  } else {
    return sqrt(-1*a*b);
  }
}

