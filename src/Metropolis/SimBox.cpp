#ifdef _OPENACC
#include <openacc.h>
#endif

#include "SimBox.h"
#include "GPUCopy.h"


// ----- Experimental -----

Real SimBox::calcAtomDistSquared(int a1, int a2, Real** aCoords, Real* bSize) {
  Real dx = makePeriodic(aCoords[X_COORD][a2] - aCoords[X_COORD][a1], X_COORD, bSize);
  Real dy = makePeriodic(aCoords[Y_COORD][a2] - aCoords[Y_COORD][a1], Y_COORD, bSize);
  Real dz = makePeriodic(aCoords[Z_COORD][a2] - aCoords[Z_COORD][a1], Z_COORD, bSize);

  return dx * dx + dy * dy + dz * dz;
}

Real SimBox::makePeriodic(Real x, int dimension, Real* bSize) {
  Real dimLength = bSize[dimension];
  if (x < -0.5 * dimLength) {
    x += dimLength;
  } else if (x > 0.5 * dimLength) {
    x -= dimLength;
  }
  return x;
}

Real SimBox::calcLJEnergy(int a1, int a2, const Real& r2, Real** aData) {

  if (r2 == 0.0) {
    return 0.0;
  } else {

    const Real sigma = calcBlending(aData[ATOM_SIGMA][a1], aData[ATOM_SIGMA][a2]);
    const Real epsilon = calcBlending(aData[ATOM_EPSILON][a1], aData[ATOM_EPSILON][a2]);

    const Real s2r2 = pow(sigma, 2) / r2;
    const Real s6r6 = pow(s2r2, 3);
    const Real s12r12 = pow(s6r6, 2);
    return 4.0 * epsilon * (s12r12 - s6r6);
  }
}

Real SimBox::calcChargeEnergy(int a1, int a2, const Real& r, Real** aData) {
  if (r == 0.0) {
    return 0.0;
  } else {
    const Real e = 332.06;
    return (aData[ATOM_CHARGE][a1] * aData[ATOM_CHARGE][a2] * e) / r;
  }
}

Real SimBox::calcBlending (const Real &a, const Real &b) {
  if (a * b >= 0)
    return sqrt(a*b);
  else
    return sqrt(-1*a*b);
}

int SimBox::findNeighbors(int molIdx) {
  int outIdx = 0;

  int pIdx = primaryIndexes[moleculeData[MOL_PIDX_START][molIdx]];

  int base[3];
  for (int i = 0; i < 3; i++) {
    base[i] = getCell(atomCoordinates[i][pIdx], i);
  }

  for (int i = -1; i <= 1; i++) {
    if (i + 1 >= numCells[0]) break;
    int c0 = wrapCell(base[0] + i, 0);
    for (int j = -1; j <= 1; j++) {
      if (j + 1 >= numCells[1]) break;
      int c1 = wrapCell(base[1] + j, 1);
      for (int k = -1; k <= 1; k++) {
        if (k + 1 >= numCells[2]) break;
        int c2 = wrapCell(base[2] + k, 2);
        neighbors[outIdx] = neighborCells[c0][c1][c2];
        outIdx++;
      }
    }
  }

  return outIdx;
}

int SimBox::wrapCell(int idx, int dimension) {
  return (idx + numCells[dimension]) % numCells[dimension];
}

int SimBox::getCell(Real loc, int dimension) {
  int idx = (int) (loc / cellWidth[dimension]);
  int numC = numCells[dimension];
  return (idx + numC) % numC;
}

void SimBox::updateNLC(int molIdx) {
  bool update = false;
  int newCell[3];
  int pIdx = primaryIndexes[moleculeData[MOL_PIDX_START][molIdx]];
  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    newCell[i] = getCell(atomCoordinates[i][pIdx], i);
    if (newCell[i] != prevCell[i])
      update = true;
  }

  if (update) {
    if (prevNode->index == molIdx) {
      neighborCells[prevCell[0]][prevCell[1]][prevCell[2]] = (prevNode->next);
      prevNode->next = neighborCells[newCell[0]][newCell[1]][newCell[2]];
      neighborCells[newCell[0]][newCell[1]][newCell[2]] = prevNode;
    } else {
      NLC_Node* removeNode = prevNode->next;
      prevNode->next = (prevNode->next->next);
      removeNode->next = (neighborCells[newCell[0]][newCell[1]][newCell[2]]);
      neighborCells[newCell[0]][newCell[1]][newCell[2]] = removeNode;
    }
  }
}

Real SimBox::calcIntraMolecularEnergy(int molIdx) {
  int molStart = moleculeData[MOL_START][molIdx];
  int molEnd = molStart + moleculeData[MOL_LEN][molIdx];
  int molType = moleculeData[MOL_TYPE][molIdx];
  Real out = 0.0;
  out += angleEnergy(molIdx);
  out += bondEnergy(molIdx);

  // Calculate intramolecular LJ and Coulomb energy if necessary
  for (int i = molStart; i < molEnd; i++) {
    for (int j = i + 1; j < molEnd; j++) {
      Real fudgeFactor = 1.0;
      for (int k = 0; ; k++) {
        int val = excludeAtoms[molType][i - molStart][k];
        if (val == -1) {
          break;
        } else if (val == j - molStart) {
          fudgeFactor = 0.0;
          break;
        }
      }
      if (fudgeFactor > 0.0) {
        for (int k = 0; ; k++) {
          int val = fudgeAtoms[molType][i - molStart][k];
          if (val == -1) {
            break;
          } else if (val == j - molStart) {
            fudgeFactor = 0.5;
            break;
          }
        }
      }
      if (fudgeFactor > 0.0) {
        Real r2 = calcAtomDistSquared(i, j, atomCoordinates, size);
        Real r = sqrt(r2);
        Real energy = calcLJEnergy(i, j, r2, atomData);
        energy += calcChargeEnergy(i, j, r, atomData);
        out += fudgeFactor * energy;
      }
    }
  }
  return out;
}

void SimBox::expandAngle(int molIdx, int angleIdx, Real expandDeg) {
  int bondStart = moleculeData[MOL_BOND_START][molIdx];
  int bondEnd = bondStart + moleculeData[MOL_BOND_COUNT][molIdx];
  int angleStart = moleculeData[MOL_ANGLE_START][molIdx];
  changedAngle = angleStart + angleIdx;
  prevAngle = angleSizes[angleStart + angleIdx];
  int startIdx = moleculeData[MOL_START][molIdx];
  int molSize = moleculeData[MOL_LEN][molIdx];
  int end1 = (int) angleData[ANGLE_A1_IDX][angleStart + angleIdx];
  int end2 = (int) angleData[ANGLE_A2_IDX][angleStart + angleIdx];
  int mid = (int) angleData[ANGLE_MID_IDX][angleStart + angleIdx];

  for (int i = 0; i < molSize; i++) {
    unionFindParent[i] = i;
  }

  for (int i = bondStart; i < bondEnd; i++) {
    int a1 = (int) bondData[BOND_A1_IDX][i];
    int a2 = (int) bondData[BOND_A2_IDX][i];
    if (a1 == mid || a2 == mid)
      continue;
    unionAtoms(a1 - startIdx, a2 - startIdx);
  }
  int group1 = find(end1 - startIdx);
  int group2 = find(end2 - startIdx);
  if (group1 == group2) {
    std::cout << "ERROR: EXPANDING ANGLE IN A RING!" << std::endl;
    return;
  }
  Real DEG2RAD = 3.14159256358979323846264 / 180.0;
  Real end1Mid[NUM_DIMENSIONS];
  Real end2Mid[NUM_DIMENSIONS];
  Real normal[NUM_DIMENSIONS];
  Real mvector[NUM_DIMENSIONS];
  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    end1Mid[i] = atomCoordinates[i][mid] - atomCoordinates[i][end1];
    end2Mid[i] = atomCoordinates[i][mid] - atomCoordinates[i][end2];
    mvector[i] = atomCoordinates[i][mid];
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
      point[j] = atomCoordinates[j][i] - mvector[i];
      dot += point[j] * normal[j];
    }

    cross[0] = normal[1] * point[2] - point[1] * normal[2];
    cross[1] = point[0] * normal[2] - normal[0] * point[2];
    cross[2] = normal[0] * point[1] - point[0] * normal[1];

    for (int j = 0; j < NUM_DIMENSIONS; j++) {
      point[j] = normal[j] * dot * (1 - cos(theta)) + point[j] * cos(theta) + cross[j] * sin(theta);
      atomCoordinates[j][i] = point[j] + mvector[j];
    }
  }

  angleSizes[angleStart + angleIdx] += expandDeg;
}

void SimBox::stretchBond(int molIdx, int bondIdx, Real stretchDist) {
  int bondStart = moleculeData[MOL_BOND_START][molIdx];
  int bondEnd = bondStart + moleculeData[MOL_BOND_COUNT][molIdx];
  changedBond = bondStart + bondIdx;
  prevBond = bondLengths[bondStart + bondIdx];
  int startIdx = moleculeData[MOL_START][molIdx];
  int molSize = moleculeData[MOL_LEN][molIdx];
  int end1 = (int) bondData[BOND_A1_IDX][bondStart + bondIdx];
  int end2 = (int) bondData[BOND_A2_IDX][bondStart + bondIdx];

  for (int i = 0; i < molSize; i++) {
    unionFindParent[i] = i;
  }

  for (int i = bondStart; i < bondEnd; i++) {
    if (i == bondIdx + bondStart)
      continue;
    int a1 = (int) bondData[BOND_A1_IDX][i] - startIdx;
    int a2 = (int) bondData[BOND_A2_IDX][i] - startIdx;
    unionAtoms(a1, a2);
  }
  int side1 = find(end1 - startIdx);
  int side2 = find(end2 - startIdx);
  if (side1 == side2) {
    std::cout << "ERROR: EXPANDING BOND IN A RING!" << std::endl;
    return;
  }
  Real v[NUM_DIMENSIONS];
  Real denon;
  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    v[i] = atomCoordinates[i][side2] - atomCoordinates[i][side1];
    denon += v[i] * v[i];
  }
  denon = sqrt(denon);
  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    v[i] = v[i] / denon / 2.0;
  }
  for (int i = 0; i < molSize; i++) {
    if(find(i) == side2) {
      for (int j = 0; j < NUM_DIMENSIONS; j++) {
        atomCoordinates[j][i + startIdx] += v[i];
      }
    } else {
      for (int j = 0; j < NUM_DIMENSIONS; j++) {
        atomCoordinates[j][i + startIdx] -= v[i];
      }
    }
  }
  bondLengths[bondStart + bondIdx] += stretchDist;
}

Real SimBox::angleEnergy(int molIdx) {
  Real out = 0;
  int angleStart = moleculeData[MOL_ANGLE_START][molIdx];
  int angleEnd = angleStart + moleculeData[MOL_ANGLE_COUNT][molIdx];
  for (int i = angleStart; i < angleEnd; i++) {
    if ((bool) angleData[ANGLE_VARIABLE][i]) {
      Real diff = angleData[ANGLE_EQANGLE][i] - angleSizes[i];
      out += angleData[ANGLE_KANGLE][i] * diff * diff;
    }
  }
  return out;
}

void SimBox::rollbackBond() {
  bondLengths[changedBond] = prevBond;
}

void SimBox::rollbackAngle() {
  angleSizes[changedAngle] = prevAngle;
}

Real SimBox::bondEnergy(int molIdx) {
  Real out = 0;
  int bondStart = moleculeData[MOL_BOND_START][molIdx];
  int bondEnd = bondStart + moleculeData[MOL_BOND_COUNT][molIdx];
  for (int i = bondStart; i < bondEnd; i++) {
    if ((bool) bondData[BOND_VARIABLE][i]) {
      Real diff = bondData[BOND_EQDIST][i] - bondLengths[i];
      out += bondData[BOND_KBOND][i] * diff * diff;
    }
  }
  return out;
}

void SimBox::unionAtoms(int atom1, int atom2) {
  int a1Parent = find(atom1);
  int a2Parent = find(atom2);
  if (a1Parent != a2Parent) {
    unionFindParent[a1Parent] = a2Parent;
  }
}

int SimBox::find(int atomIdx) {
  if (unionFindParent[atomIdx] == atomIdx) {
    return atomIdx;
  } else {
    unionFindParent[atomIdx] = find(unionFindParent[atomIdx]);
    return unionFindParent[atomIdx];
  }
}
