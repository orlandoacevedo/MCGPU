#include "SimCalcs.h"

SimBox* sb;
int on_gpu;

bool SimCalcs::moleculesInRange(int p1Start, int p1End, int p2Start, int p2End,
    Real** atomCoords, Real* bSize, int* primaryIndexes, Real cutoff) {

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

Real SimCalcs::calcAtomDistSquared(int a1, int a2, Real** aCoords,
    Real* bSize) {

  Real dx = makePeriodic(aCoords[X_COORD][a2] - aCoords[X_COORD][a1],
      X_COORD, bSize);
  Real dy = makePeriodic(aCoords[Y_COORD][a2] - aCoords[Y_COORD][a1],
      Y_COORD, bSize);
  Real dz = makePeriodic(aCoords[Z_COORD][a2] - aCoords[Z_COORD][a1],
      Z_COORD, bSize);

  return dx * dx + dy * dy + dz * dz;
}

Real SimCalcs::makePeriodic(Real x, int dimension, Real* bSize) {
  Real dimLength = bSize[dimension];

  int lt = (x < -0.5 * dimLength); // 1 or 0
  x += lt * dimLength;
  int gt = (x > 0.5 * dimLength);  // 1 or 0
  x -= gt * dimLength;
  return x;
}

Real SimCalcs::calcMolecularEnergyContribution(int currMol, int startMol) {
  Real total = 0;

  int** molData = GPUCopy::moleculeDataPtr();
  Real** atomCoords = GPUCopy::atomCoordinatesPtr();
  Real* bSize = GPUCopy::sizePtr();
  int* pIdxes = GPUCopy::primaryIndexesPtr();
  Real** aData = GPUCopy::atomDataPtr();
  Real cutoff = sb->cutoff;
  const long numMolecules = sb->numMolecules;

  const int p1Start = sb->moleculeData[MOL_PIDX_START][currMol];
  const int p1End   = sb->moleculeData[MOL_PIDX_COUNT][currMol] + p1Start;

  #pragma acc parallel loop gang deviceptr(molData, atomCoords, bSize, pIdxes, aData) if (on_gpu) \
              vector_length(64)
  for (int otherMol = startMol; otherMol < numMolecules; otherMol++) {
    if (otherMol != currMol) {
      int p2Start = molData[MOL_PIDX_START][otherMol];
      int p2End = molData[MOL_PIDX_COUNT][otherMol] + p2Start;
      if (moleculesInRange(p1Start, p1End, p2Start, p2End, atomCoords, bSize, pIdxes, cutoff)) {
        total += calcMoleculeInteractionEnergy(currMol, otherMol, molData, aData, atomCoords, bSize);
      }
    }
  }

  return total;
}

Real SimCalcs::calcMoleculeInteractionEnergy (int m1, int m2, int** molData,
    Real** aData, Real** aCoords, Real* bSize) {

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

Real SimCalcs::calcLJEnergy(int a1, int a2, Real r2, Real** aData) {

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

Real SimCalcs::calcChargeEnergy(int a1, int a2, Real r, Real** aData) {
  if (r == 0.0) {
    return 0.0;
  } else {
    const Real e = 332.06;
    return (aData[ATOM_CHARGE][a1] * aData[ATOM_CHARGE][a2] * e) / r;
  }
}

Real SimCalcs::calcBlending (Real a, Real b) {
  if (a * b >= 0) {
    return sqrt(a*b);
  } else {
    return sqrt(-1*a*b);
  }
}

void SimCalcs::rotateAtom(int aIdx, int pivotIdx, Real rotX, Real rotY, Real rotZ, Real** aCoords) {
  Real pX = aCoords[X_COORD][pivotIdx];
  Real pY = aCoords[Y_COORD][pivotIdx];
  Real pZ = aCoords[Z_COORD][pivotIdx];

  translateAtom(aIdx, -pX, -pY, -pZ, aCoords);
  rotateX(aIdx, rotX, aCoords);
  rotateY(aIdx, rotY, aCoords);
  rotateZ(aIdx, rotZ, aCoords);
  translateAtom(aIdx, pX, pY, pZ, aCoords);
}


void SimCalcs::rotateX(int aIdx, Real angleDeg, Real** aCoords) {
  Real angleRad = angleDeg * 3.14159265358979 / 180.0;
  Real oldY = aCoords[Y_COORD][aIdx];
  Real oldZ = aCoords[Z_COORD][aIdx];
  aCoords[Y_COORD][aIdx] = oldY * cos(angleRad) + oldZ * sin(angleRad);
  aCoords[Z_COORD][aIdx] = oldZ * cos(angleRad) - oldY * sin(angleRad);
}

void SimCalcs::rotateY(int aIdx, Real angleDeg, Real** aCoords) {
  Real angleRad = angleDeg * 3.14159265358979 / 180.0;
  Real oldZ = aCoords[Z_COORD][aIdx];
  Real oldX = aCoords[X_COORD][aIdx];
  aCoords[Z_COORD][aIdx] = oldZ * cos(angleRad) + oldX * sin(angleRad);
  aCoords[X_COORD][aIdx] = oldX * cos(angleRad) - oldZ * sin(angleRad);
}

void SimCalcs::rotateZ(int aIdx, Real angleDeg, Real** aCoords) {
  Real angleRad = angleDeg * 3.14159265358979 / 180.0;
  Real oldX = aCoords[X_COORD][aIdx];
  Real oldY = aCoords[Y_COORD][aIdx];
  aCoords[X_COORD][aIdx] = oldX * cos(angleRad) + oldY * sin(angleRad);
  aCoords[Y_COORD][aIdx] = oldY * cos(angleRad) - oldX * sin(angleRad);
}

void SimCalcs::changeMolecule(int molIdx) {
  Real maxT = sb->maxTranslate;
  Real maxR = sb->maxRotate;

  int molStart = sb->moleculeData[MOL_START][molIdx];
  int molLen = sb->moleculeData[MOL_LEN][molIdx];

  int vertexIdx = (int) randomReal(0, molLen);

  const Real deltaX = randomReal(-maxT, maxT);
  const Real deltaY = randomReal(-maxT, maxT);
  const Real deltaZ = randomReal(-maxT, maxT);

  const Real rotX = randomReal(-maxR, maxR);
  const Real rotY = randomReal(-maxR, maxR);
  const Real rotZ = randomReal(-maxR, maxR);

  Real ** rBCoords = GPUCopy::rollBackCoordinatesPtr();
  Real ** aCoords = GPUCopy::atomCoordinatesPtr();
  Real * bSize = GPUCopy::sizePtr();
  int* pIdxes = GPUCopy::primaryIndexesPtr();
  int** molData = GPUCopy::moleculeDataPtr();

  // do the move here
  #pragma acc parallel loop deviceptr(aCoords, rBCoords) if(on_gpu)
  for (int i = 0; i < molLen; i++) {
    for (int j = 0; j < NUM_DIMENSIONS; j++) {
      rBCoords[j][i] = aCoords[j][molStart + i];
    }
    if (i == vertexIdx)
      continue;
    rotateAtom(molStart + i, molStart + vertexIdx, rotX, rotY, rotZ, aCoords);
    translateAtom(molStart + i, deltaX, deltaY, deltaZ, aCoords);
  }

  #pragma acc parallel loop deviceptr(aCoords, molData, pIdxes, bSize) if(on_gpu)
  for (int i = 0; i < 1; i++) {
    aCoords[0][molStart + vertexIdx] += deltaX;
    aCoords[1][molStart + vertexIdx] += deltaY;
    aCoords[2][molStart + vertexIdx] += deltaZ;
    keepMoleculeInBox(molIdx, aCoords, molData, pIdxes, bSize);
  }
}

void SimCalcs::translateAtom(int aIdx, Real dX, Real dY, Real dZ, Real** aCoords) {
  aCoords[X_COORD][aIdx] += dX;
  aCoords[Y_COORD][aIdx] += dY;
  aCoords[Z_COORD][aIdx] += dZ;
}

void SimCalcs::keepMoleculeInBox(int molIdx, Real** aCoords, int** molData, int* pIdxes, Real* bSize) {
  int start = molData[MOL_START][molIdx];
  int end = start + molData[MOL_LEN][molIdx];
  int pIdx = pIdxes[molData[MOL_PIDX_START][molIdx]];

  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    if (aCoords[i][pIdx] < 0) {
      #pragma acc loop independent 
      for (int j = start; j < end; j++) {
        aCoords[i][j] += bSize[i];
      }
    } else if (aCoords[i][pIdx] > bSize[i]) {
      #pragma acc loop independent
      for (int j = start; j < end; j++) {
        aCoords[i][j] -= bSize[i];
      }
    }
  }
}

void SimCalcs::rollback(int molIdx) {

  int molStart = sb->moleculeData[MOL_START][molIdx];
  int molLen = sb->moleculeData[MOL_LEN][molIdx];

  Real** aCoords = GPUCopy::atomCoordinatesPtr();
  Real** rBCoords = GPUCopy::rollBackCoordinatesPtr();

  #pragma acc parallel loop deviceptr(aCoords, rBCoords) if(on_gpu)
  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    #pragma acc loop independent 
    for (int j = 0; j < molLen; j++) {
      aCoords[i][molStart + j] = rBCoords[i][j];
    }
  }
}

void SimCalcs::setSB(SimBox* sb_in) {
  sb = sb_in;
  on_gpu = GPUCopy::onGpu();
}
