/**
 * Base class for an iteration of the simulation.
 */

#include "SimBox.h"
#include "GPUCopy.h"
#include "SimulationStep.h"

/** Construct a new SimulationStep from a SimBox pointer */
SimulationStep::SimulationStep(SimBox *box) {
  SimCalcs::setSB(box);
}

Real SimulationStep::calcMoleculeEnergy(int currMol, int startMol) {
  return calcMolecularEnergyContribution(currMol, startMol) +
    calcIntraMolecularEnergy(currMol);
}

Real SimulationStep::calcIntraMolecularEnergy(int molIdx) {
  return SimCalcs::calcIntraMolecularEnergy(molIdx);
}

/** Returns the index of a random molecule within the simulation box */
int SimulationStep::chooseMolecule(SimBox *box) {
  return (int) randomReal(0, box->numMolecules);
}

/** Perturb a given molecule */
void SimulationStep::changeMolecule(int molIdx, SimBox *box) {
  SimCalcs::changeMolecule(molIdx);
}

/** Move a molecule back to its original position */
void SimulationStep::rollback(int molIdx, SimBox *box) {
  SimCalcs::rollback(molIdx);
}

/** Determines the total energy of the box */
Real SimulationStep::calcSystemEnergy(Real &subLJ, Real &subCharge,
                                      int numMolecules) {
  Real total = subLJ + subCharge;
  for (int mol = 0; mol < numMolecules; mol++) {
    total += calcMolecularEnergyContribution(mol, mol) +
      calcIntraMolecularEnergy(mol);
  }

  return total;
}


// ----- SimCalcs Definitions -----


SimBox* SimCalcs::sb;
int SimCalcs::on_gpu;

Real SimCalcs::calcIntraMolecularEnergy(int molIdx) {
  int molStart = sb->moleculeData[MOL_START][molIdx];
  int molEnd = molStart + sb->moleculeData[MOL_LEN][molIdx];
  int molType = sb->moleculeData[MOL_TYPE][molIdx];
  Real out = 0.0;
  out += angleEnergy(molIdx);
  out += bondEnergy(molIdx);

  // Calculate intramolecular LJ and Coulomb energy if necessary
  for (int i = molStart; i < molEnd; i++) {
    for (int j = i + 1; j < molEnd; j++) {
      Real fudgeFactor = 1.0;
      for (int k = 0; ; k++) {
        int val = sb->excludeAtoms[molType][i - molStart][k];
        if (val == -1) {
          break;
        } else if (val == j - molStart) {
          fudgeFactor = 0.0;
          break;
        }
      }
      if (fudgeFactor > 0.0) {
        for (int k = 0; ; k++) {
          int val = sb->fudgeAtoms[molType][i - molStart][k];
          if (val == -1) {
            break;
          } else if (val == j - molStart) {
            fudgeFactor = 0.5;
            break;
          }
        }
      }
      if (fudgeFactor > 0.0) {
        Real r2 = calcAtomDistSquared(i, j, sb->atomCoordinates, sb->size);
        Real r = sqrt(r2);
        Real energy = calcLJEnergy(i, j, r2, sb->atomData);
        energy += calcChargeEnergy(i, j, r, sb->atomData);
        out += fudgeFactor * energy;
      }
    }
  }
  return out;
}

Real SimCalcs::angleEnergy(int molIdx) {
  Real out = 0;
  int angleStart = sb->moleculeData[MOL_ANGLE_START][molIdx];
  int angleEnd = angleStart + sb->moleculeData[MOL_ANGLE_COUNT][molIdx];
  for (int i = angleStart; i < angleEnd; i++) {
    if ((bool) sb->angleData[ANGLE_VARIABLE][i]) {
      Real diff = sb->angleData[ANGLE_EQANGLE][i] - sb->angleSizes[i];
      out += sb->angleData[ANGLE_KANGLE][i] * diff * diff;
    }
  }
  return out;
}

void SimCalcs::expandAngle(int molIdx, int angleIdx, Real expandDeg) {
  int bondStart = sb->moleculeData[MOL_BOND_START][molIdx];
  int bondEnd = bondStart + sb->moleculeData[MOL_BOND_COUNT][molIdx];
  int angleStart = sb->moleculeData[MOL_ANGLE_START][molIdx];
  int startIdx = sb->moleculeData[MOL_START][molIdx];
  int molSize = sb->moleculeData[MOL_LEN][molIdx];
  int end1 = (int)sb->angleData[ANGLE_A1_IDX][angleStart + angleIdx];
  int end2 = (int)sb->angleData[ANGLE_A2_IDX][angleStart + angleIdx];
  int mid = (int)sb->angleData[ANGLE_MID_IDX][angleStart + angleIdx];


  // Create a disjoint set of the atoms in the molecule
  for (int i = 0; i < molSize; i++) {
    sb->unionFindParent[i] = i;
  }

  // Union atoms connected by a bond
  for (int i = bondStart; i < bondEnd; i++) {
    int a1 = (int)sb->bondData[BOND_A1_IDX][i];
    int a2 = (int)sb->bondData[BOND_A2_IDX][i];
    if (a1 == mid || a2 == mid)
      continue;
    unionAtoms(a1 - startIdx, a2 - startIdx);
  }

  int group1 = find(end1 - startIdx);
  int group2 = find(end2 - startIdx);
  if (group1 == group2) {
    // std::cout << "ERROR: EXPANDING ANGLE IN A RING!" << std::endl;
    return;
  }
  Real DEG2RAD = 3.14159256358979323846264 / 180.0;
  Real end1Mid[NUM_DIMENSIONS];
  Real end2Mid[NUM_DIMENSIONS];
  Real normal[NUM_DIMENSIONS];
  Real mvector[NUM_DIMENSIONS];
  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    end1Mid[i] = sb->atomCoordinates[i][mid] - sb->atomCoordinates[i][end1];
    end2Mid[i] = sb->atomCoordinates[i][mid] - sb->atomCoordinates[i][end2];
    mvector[i] = sb->atomCoordinates[i][mid];
  }
  normal[0] = end1Mid[1] * end2Mid[2] - end2Mid[1] * end1Mid[2];
  normal[1] = end2Mid[0] * end1Mid[2] - end1Mid[0] * end2Mid[2];
  normal[2] = end1Mid[0] * end2Mid[1] - end2Mid[0] * end1Mid[1];
  Real normLen = 0.0;
  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    normLen += normal[i] * normal[i];
  }
  normLen = sqrt(normLen);
  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    normal[i] = normal[i] / normLen;
  }


  for (int i = startIdx; i < startIdx + molSize; i++) {
    Real theta;
    Real point[NUM_DIMENSIONS];
    Real dot = 0.0;
    Real cross[NUM_DIMENSIONS];
    if (find(i - startIdx) == group1) {
      theta = expandDeg * -DEG2RAD;
    } else if (find(i - startIdx) == group2) {
      theta = expandDeg * DEG2RAD;
    } else {
      continue;
    }

    for (int j = 0; j < NUM_DIMENSIONS; j++) {
      point[j] = sb->atomCoordinates[j][i] - mvector[j];
      dot += point[j] * normal[j];
    }

    cross[0] = normal[1] * point[2] - point[1] * normal[2];
    cross[1] = point[0] * normal[2] - normal[0] * point[2];
    cross[2] = normal[0] * point[1] - point[0] * normal[1];

    for (int j = 0; j < NUM_DIMENSIONS; j++) {
      point[j] = (normal[j] * dot * (1 - cos(theta)) + point[j] * cos(theta) +
                  cross[j] * sin(theta));
      sb->atomCoordinates[j][i] = point[j] + mvector[j];
    }
  }

  sb->angleSizes[angleStart + angleIdx] += expandDeg;
}

Real SimCalcs::bondEnergy(int molIdx) {
  Real out = 0;
  int bondStart = sb->moleculeData[MOL_BOND_START][molIdx];
  int bondEnd = bondStart + sb->moleculeData[MOL_BOND_COUNT][molIdx];
  for (int i = bondStart; i < bondEnd; i++) {
    if ((bool) sb->bondData[BOND_VARIABLE][i]) {
      Real diff = sb->bondData[BOND_EQDIST][i] - sb->bondLengths[i];
      out += sb->bondData[BOND_KBOND][i] * diff * diff;
    }
  }
  return out;
}

void SimCalcs::stretchBond(int molIdx, int bondIdx, Real stretchDist) {
  int bondStart = sb->moleculeData[MOL_BOND_START][molIdx];
  int bondEnd = bondStart + sb->moleculeData[MOL_BOND_COUNT][molIdx];
  int startIdx = sb->moleculeData[MOL_START][molIdx];
  int molSize = sb->moleculeData[MOL_LEN][molIdx];
  int end1 = (int)sb->bondData[BOND_A1_IDX][bondStart + bondIdx];
  int end2 = (int)sb->bondData[BOND_A2_IDX][bondStart + bondIdx];

  for (int i = 0; i < molSize; i++) {
    sb->unionFindParent[i] = i;
  }

  // Split the molecule atoms into two disjoint sets around the bond
  for (int i = bondStart; i < bondEnd; i++) {
    if (i == bondIdx + bondStart)
      continue;
    int a1 = (int)sb->bondData[BOND_A1_IDX][i] - startIdx;
    int a2 = (int)sb->bondData[BOND_A2_IDX][i] - startIdx;
    unionAtoms(a1, a2);
  }
  int side1 = find(end1 - startIdx);
  int side2 = find(end2 - startIdx);
  if (side1 == side2) {
    // std::cerr << "ERROR: EXPANDING BOND IN A RING!" << std::endl;
    return;
  }

  // Move each atom the appropriate distance for the bond stretch
  Real v[NUM_DIMENSIONS];
  Real denon = 0.0;
  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    v[i] = sb->atomCoordinates[i][side2] - sb->atomCoordinates[i][side1];
    denon += v[i] * v[i];
  }
  denon = sqrt(denon);
  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    v[i] = v[i] / denon / 2.0;
  }
  for (int i = 0; i < molSize; i++) {
    if (find(i) == side2) {
      for (int j = 0; j < NUM_DIMENSIONS; j++) {
        sb->atomCoordinates[j][i + startIdx] += v[j];
      }
    } else {
      for (int j = 0; j < NUM_DIMENSIONS; j++) {
        sb->atomCoordinates[j][i + startIdx] -= v[j];
      }
    }
  }

  // Record the actual bond stretch
  sb->bondLengths[bondStart + bondIdx] += stretchDist;
}

bool SimCalcs::moleculesInRange(int p1Start, int p1End, int p2Start, int p2End,
                                Real** atomCoords, Real* bSize,
                                int* primaryIndexes, Real cutoff) {
  bool out = false;
  for (int p1Idx = p1Start; p1Idx < p1End; p1Idx++) {
    int p1 = primaryIndexes[p1Idx];
    for (int p2Idx = p2Start; p2Idx < p2End; p2Idx++) {
      int p2 = primaryIndexes[p2Idx];
      out |= (calcAtomDistSquared(p1, p2, atomCoords, bSize) <=
              cutoff * cutoff);
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

Real SimCalcs::calcLJEnergy(int a1, int a2, Real r2, Real** aData) {
  if (r2 == 0.0) {
    return 0.0;
  } else {
    const Real sigma = SimCalcs::calcBlending(aData[ATOM_SIGMA][a1],
        aData[ATOM_SIGMA][a2]);
    const Real epsilon = SimCalcs::calcBlending(aData[ATOM_EPSILON][a1],
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

Real SimCalcs::makePeriodic(Real x, int dimension, Real* bSize) {
  Real dimLength = bSize[dimension];

  int lt = (x < -0.5 * dimLength); // 1 or 0
  x += lt * dimLength;
  int gt = (x > 0.5 * dimLength);  // 1 or 0
  x -= gt * dimLength;
  return x;
}

void SimCalcs::rotateAtom(int aIdx, int pivotIdx, Real rotX, Real rotY,
                          Real rotZ, Real** aCoords) {
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
  // Intermolecular moves first, to save proper rollback positions
  intermolecularMove(molIdx);
  intramolecularMove(molIdx);
}

void SimCalcs::intermolecularMove(int molIdx) {
  Real maxT = sb->maxTranslate;
  Real maxR = sb->maxRotate;

  int molStart = sb->moleculeData[MOL_START][molIdx];
  int molLen = sb->moleculeData[MOL_LEN][molIdx];

  int vertexIdx = (int)randomReal(0, molLen);

  const Real deltaX = randomReal(-maxT, maxT);
  const Real deltaY = randomReal(-maxT, maxT);
  const Real deltaZ = randomReal(-maxT, maxT);

  const Real rotX = randomReal(-maxR, maxR);
  const Real rotY = randomReal(-maxR, maxR);
  const Real rotZ = randomReal(-maxR, maxR);

  Real** rBCoords = GPUCopy::rollBackCoordinatesPtr();
  Real** aCoords = GPUCopy::atomCoordinatesPtr();
  Real* bSize = GPUCopy::sizePtr();
  int* pIdxes = GPUCopy::primaryIndexesPtr();
  int** molData = GPUCopy::moleculeDataPtr();

  // Do the move here
  #pragma acc parallel loop deviceptr(aCoords, rBCoords) \
      if (on_gpu)
  for (int i = 0; i < molLen; i++) {
    for (int j = 0; j < NUM_DIMENSIONS; j++) {
      rBCoords[j][i] = aCoords[j][molStart + i];
    }
    if (i == vertexIdx)
      continue;
    rotateAtom(molStart + i, molStart + vertexIdx, rotX, rotY, rotZ, aCoords);
    translateAtom(molStart + i, deltaX, deltaY, deltaZ, aCoords);
  }

  #pragma acc parallel loop deviceptr(aCoords, molData, pIdxes, bSize) \
      if (on_gpu)
  for (int i = 0; i < 1; i++) {
    aCoords[0][molStart + vertexIdx] += deltaX;
    aCoords[1][molStart + vertexIdx] += deltaY;
    aCoords[2][molStart + vertexIdx] += deltaZ;
    keepMoleculeInBox(molIdx, aCoords, molData, pIdxes, bSize);
  }
}

void SimCalcs::intramolecularMove(int molIdx) {
  // Save the molecule data for rolling back
  // TODO (blm): Put these in the GPU with GPUCopy
  saveBonds(molIdx);
  saveAngles(molIdx);

  // TODO (blm): allow max to be configurable
  int numMoves = (int)round(randomReal(1, sb->maxIntraMoves)); 
  int numBonds = sb->moleculeData[MOL_BOND_COUNT][molIdx];
  int numAngles = sb->moleculeData[MOL_ANGLE_COUNT][molIdx];
  Real bondDelta = sb->maxBondDelta, angleDelta = sb->maxAngleDelta;
  for (int i = 0; i < numMoves; i++) {
    Real moveType = randomReal(0, 1);
    int selectedBond, selectedAngle;
    if (moveType > 0.5) {
        if (numBonds == 0) break;
        selectedBond = (int)randomReal(0, numBonds);
        // TODO (blm): Make bond delta more accurate
        stretchBond(molIdx, selectedBond, randomReal(-bondDelta, bondDelta));
        // TODO (blm): Do an MC test to self-correct bond delta
    } else {
        if (numAngles == 0) break;
        selectedAngle = (int)randomReal(0, numAngles);
        // TODO (blm): Make bond delta more accurate
        expandAngle(molIdx, selectedAngle, randomReal(-angleDelta, angleDelta));
        // TODO (blm): Do an MC test to self-correct angle delta
    } // TODO: Add dihedrals here
  }
}

void SimCalcs::saveBonds(int molIdx) {
  int start = sb->moleculeData[MOL_BOND_START][molIdx];
  int end = sb->moleculeData[MOL_BOND_COUNT][molIdx];

  for (int i = start; i < end; i++) {
    sb->rollBackBondLengths[i] = sb->bondLengths[i];
  }
}

void SimCalcs::saveAngles(int molIdx) {
  int start = sb->moleculeData[MOL_ANGLE_START][molIdx];
  int end = sb->moleculeData[MOL_ANGLE_COUNT][molIdx];

  for (int i = start; i < end; i++) {
    sb->rollBackAngleSizes[i] = sb->angleSizes[i];
  }
}

void SimCalcs::translateAtom(int aIdx, Real dX, Real dY, Real dZ,
                             Real** aCoords) {
  aCoords[X_COORD][aIdx] += dX;
  aCoords[Y_COORD][aIdx] += dY;
  aCoords[Z_COORD][aIdx] += dZ;
}

void SimCalcs::keepMoleculeInBox(int molIdx, Real** aCoords, int** molData,
                                 int* pIdxes, Real* bSize) {
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

  #pragma acc parallel loop deviceptr(aCoords, rBCoords) if (on_gpu)
  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    #pragma acc loop independent 
    for (int j = 0; j < molLen; j++) {
      aCoords[i][molStart + j] = rBCoords[i][j];
    }
  }

  rollbackAngles(molIdx);
  rollbackBonds(molIdx);
}

void SimCalcs::rollbackBonds(int molIdx) {
  int start = sb->moleculeData[MOL_BOND_START][molIdx];
  int end = sb->moleculeData[MOL_BOND_COUNT][molIdx];

  for (int i = start; i < end; i++) {
    sb->angleSizes[i] = sb->rollBackAngleSizes[i];
  }
}

void SimCalcs::rollbackAngles(int molIdx) {
  int start = sb->moleculeData[MOL_ANGLE_START][molIdx];
  int end = sb->moleculeData[MOL_ANGLE_COUNT][molIdx];

  for (int i = start; i < end; i++) {
    sb->angleSizes[i] = sb->rollBackAngleSizes[i];
  }
}

void SimCalcs::unionAtoms(int atom1, int atom2) {
  int a1Parent = find(atom1);
  int a2Parent = find(atom2);
  if (a1Parent != a2Parent) {
    sb->unionFindParent[a1Parent] = a2Parent;
  }
}

int SimCalcs::find(int atomIdx) {
  if (sb->unionFindParent[atomIdx] == atomIdx) {
    return atomIdx;
  } else {
    sb->unionFindParent[atomIdx] = find(sb->unionFindParent[atomIdx]);
    return sb->unionFindParent[atomIdx];
  }
}

void SimCalcs::setSB(SimBox* sb_in) {
  sb = sb_in;
  on_gpu = GPUCopy::onGpu();
}
