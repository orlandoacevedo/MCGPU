#ifdef _OPENACC
#include <openacc.h>
#endif

#include "SimBox.h"
#include "GPUCopy.h"
#include "SimCalcs.h"

void SimBox::keepMoleculeInBox(int molIdx, Real** aCoords, int** molData, int* pIdxes, Real* bSize) {
  int start = molData[MOL_START][molIdx];
  int end = start + molData[MOL_LEN][molIdx];
  int pIdx = pIdxes[molData[MOL_PIDX_START][molIdx]];

  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    if (aCoords[i][pIdx] < 0) {
      #pragma acc loop independent
      for (int j = start; j < end; j++) {
        aCoords[i][j] += bSize[i];
      }
    } else if (aCoords[i][pIdx] > size[i]) {
      #pragma acc loop independent
      for (int j = start; j < end; j++) {
        aCoords[i][j] -= bSize[i];
      }
    }
  }
}

void SimBox::changeMolecule(int molIdx) {
  SimCalcs::changeMolecule(molIdx);
  /*
  Real maxT = maxTranslate;
  Real maxR = maxRotate;

  int molStart = moleculeData[MOL_START][molIdx];
  int molLen = moleculeData[MOL_LEN][molIdx];

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
  #pragma acc parallel loop deviceptr(rBCoords, aCoords)
  for (int i = 0; i < molLen; i++) {
    #pragma acc loop independent
    for (int j = 0; j < NUM_DIMENSIONS; j++) {
      rBCoords[j][i] = aCoords[j][molStart + i];
    }
    if (i == vertexIdx)
      continue;
    rotateAtom(molStart + i, molStart + vertexIdx, rotX, rotY, rotZ, aCoords);
    translateAtom(molStart + i, deltaX, deltaY, deltaZ, aCoords);
  }

  #pragma acc parallel loop deviceptr(aCoords, molData, pIdxes, bSize)
  for (int i = 0; i < 1; i++) {
    aCoords[0][molStart + vertexIdx] += deltaX;
    aCoords[1][molStart + vertexIdx] += deltaY;
    aCoords[2][molStart + vertexIdx] += deltaZ;
    keepMoleculeInBox(molIdx, aCoords, molData, pIdxes, bSize);
  }

  if (useNLC) {
    updateNLC(molIdx);
  } */
}

void SimBox::rollback(int molIdx) {
  SimCalcs::rollback(molIdx);
  /*
  int molStart = moleculeData[MOL_START][molIdx];
  int molLen = moleculeData[MOL_LEN][molIdx];

  Real** aCoords = GPUCopy::atomCoordinatesPtr();
  Real** rBCoords = GPUCopy::rollBackCoordinatesPtr();

  #pragma acc parallel loop deviceptr(aCoords, rBCoords)
  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    #pragma acc loop independent
    for (int j = 0; j < molLen; j++) {
      aCoords[i][molStart + j] = rBCoords[i][j];
    }
  }

  if (useNLC) {
    updateNLC(molIdx);
  } */
}

void SimBox::translateAtom(int aIdx, Real dX, Real dY, Real dZ, Real** aCoords) {
  aCoords[X_COORD][aIdx] += dX;
  aCoords[Y_COORD][aIdx] += dY;
  aCoords[Z_COORD][aIdx] += dZ;
}

void SimBox::rotateAtom(int aIdx, int pivotIdx, Real rotX, Real rotY, Real rotZ, Real** aCoords) {
  Real pX = aCoords[X_COORD][pivotIdx];
  Real pY = aCoords[Y_COORD][pivotIdx];
  Real pZ = aCoords[Z_COORD][pivotIdx];

  translateAtom(aIdx, -pX, -pY, -pZ, aCoords);
  rotateX(aIdx, rotX, aCoords);
  rotateY(aIdx, rotY, aCoords);
  rotateZ(aIdx, rotZ, aCoords);
  translateAtom(aIdx, pX, pY, pZ, aCoords);
}

void SimBox::rotateX(int aIdx, Real angleDeg, Real** aCoords) {
//Real angleRad = degreesToRadians(angleDeg);
  Real angleRad = angleDeg * 3.14159265358979 / 180.0;
  Real oldY = aCoords[Y_COORD][aIdx];
  Real oldZ = aCoords[Z_COORD][aIdx];
  aCoords[Y_COORD][aIdx] = oldY * cos(angleRad) + oldZ * sin(angleRad);
  aCoords[Z_COORD][aIdx] = oldZ * cos(angleRad) - oldY * sin(angleRad);
}

void SimBox::rotateY(int aIdx, Real angleDeg, Real** aCoords) {
  //Real angleRad = degreesToRadians(angleDeg);
  Real angleRad = angleDeg * 3.14159265358979 / 180.0;
  Real oldZ = aCoords[Z_COORD][aIdx];
  Real oldX = aCoords[X_COORD][aIdx];
  aCoords[Z_COORD][aIdx] = oldZ * cos(angleRad) + oldX * sin(angleRad);
  aCoords[X_COORD][aIdx] = oldX * cos(angleRad) - oldZ * sin(angleRad);
}

void SimBox::rotateZ(int aIdx, Real angleDeg, Real** aCoords) {
  //Real angleRad = degreesToRadians(angleDeg);
  Real angleRad = angleDeg * 3.14159265358979 / 180.0;
  Real oldX = aCoords[X_COORD][aIdx];
  Real oldY = aCoords[Y_COORD][aIdx];
  aCoords[X_COORD][aIdx] = oldX * cos(angleRad) + oldY * sin(angleRad);
  aCoords[Y_COORD][aIdx] = oldY * cos(angleRad) - oldX * sin(angleRad);
}

int SimBox::chooseMolecule() const {
  return (int) randomReal(0, numMolecules);
}

Real SimBox::calcSystemEnergy (Real &subLJ, Real &subCharge) {
  Real total = subLJ + subCharge;

  for (int mol = 0; mol < numMolecules; mol++) {
    total += calcMolecularEnergyContribution(subLJ, subCharge, mol, mol);
  }

  return total;
}

Real SimBox::calcMolecularEnergyContribution(Real &subLJ, Real &subCharge, int currMol, int startMol) {

  return SimCalcs::calcMolecularEnergyContribution(currMol, startMol);

 /*
  Real total = 0;
  const int p1Start = moleculeData[MOL_PIDX_START][currMol];
  const int p1End   = moleculeData[MOL_PIDX_COUNT][currMol] + p1Start;
  const int p1 = primaryIndexes[p1Start];

  if (useNLC) {
      int p2Start, p2End;
    for (int i = 0; i < NUM_DIMENSIONS; i++) {
      prevCell[i] = getCell(atomCoordinates[i][p1], i);
    }
    int neighborSize = findNeighbors(currMol);
    for (int cellHead = 0; cellHead < neighborSize; cellHead++) {
      NLC_Node* head = neighbors[cellHead];

      if (head->index == currMol)
        prevNode = head;
      while(head->next != NULL) {
        int idx = head->index;
        if (head->next->index == currMol)
          prevNode = head;
        if (idx >= startMol && idx != currMol) {
          p2Start = moleculeData[MOL_PIDX_START][idx];
          p2End = moleculeData[MOL_PIDX_COUNT][idx] + p2Start;
          if (moleculesInRange(p1Start, p1End, p2Start, p2End)) {
            total += calcMoleculeInteractionEnergy(currMol, idx);
          }
        }
        head = head->next;
      }
    }

  } else {


    int** molData = GPUCopy::moleculeDataPtr();
    Real** aCoords = GPUCopy::atomCoordinatesPtr();
    int* pIdxes = GPUCopy::primaryIndexesPtr();
    Real** aData = GPUCopy::atomDataPtr();
    Real* bSize = GPUCopy::sizePtr();

  //  #pragma acc parallel loop deviceptr(molData, aCoords, pIdxes, aData) reduction(+:total)
    for (int otherMol = startMol; otherMol < numMolecules; otherMol++) {
      if (otherMol != currMol) {
        int p2Start = moleculeData[MOL_PIDX_START][otherMol];
        int p2End = moleculeData[MOL_PIDX_COUNT][otherMol] + p2Start;
        if (moleculesInRange(p1Start, p1End, p2Start, p2End)) {
          total += calcMoleculeInteractionEnergy(currMol, otherMol);
        }
      }
    }
  }

  return total; */
}

bool SimBox::moleculesInRange(int p1Start, int p1End, int p2Start, int p2End) {

  int p1;
  int p2;

  bool out = false;

  for (int p1Idx = p1Start; p1Idx < p1End; p1Idx++) {
    p1 = primaryIndexes[p1Idx];
    for (int p2Idx = p2Start; p2Idx < p2End; p2Idx++) {
      p2 = primaryIndexes[p2Idx];
      if (calcAtomDistSquared(p1, p2, atomCoordinates, size) <= cutoff * cutoff) {
        out = true;
      }
    }
  }

  return out;
}

Real SimBox::calcMoleculeInteractionEnergy (int m1, int m2) {
  Real energySum = 0;

  const int m1Start = moleculeData[MOL_START][m1];
  const int m1End = moleculeData[MOL_LEN][m1] + m1Start;

  const int m2Start = moleculeData[MOL_START][m2];
  const int m2End = moleculeData[MOL_LEN][m2] + m2Start;

  Real** aData = GPUCopy::atomDataPtr();
  Real** aCoords = GPUCopy::atomCoordinatesPtr();
  Real* bSize = GPUCopy::sizePtr();

  #pragma acc parallel loop deviceptr(aData, aCoords, bSize) reduction(+:energySum)
  for (int i = m1Start; i < m1End; i++) {
    #pragma acc loop independent
    for (int j = m2Start; j < m2End; j++) {
      if (aData[ATOM_SIGMA][i] >= 0 && aData[ATOM_SIGMA][j] >= 0
          && aData[ATOM_EPSILON][i] >= 0 && aData[ATOM_EPSILON][j] >= 0) {

        const Real r2 = calcAtomDistSquared(i, j, aCoords, bSize);
        if (r2 == 0.0) {
          energySum += 0.0;
        } else {

          const Real sigma = calcBlending(aData[ATOM_SIGMA][i], aData[ATOM_SIGMA][j]);
          const Real epsilon = calcBlending(aData[ATOM_EPSILON][i], aData[ATOM_EPSILON][j]);

          const Real s2r2 = pow(sigma, 2) / r2;
          const Real s6r6 = pow(s2r2, 3);
          const Real s12r12 = pow(s6r6, 2);
          energySum +=  4.0 * epsilon * (s12r12 - s6r6);
        }

        if (r2 == 0.0) {
          energySum += 0.0;
        } else {
          const Real e = 332.06;
          energySum += (aData[ATOM_CHARGE][i] * aData[ATOM_CHARGE][j] * e) / sqrt(r2);
        }

      }
    }
  }

  return (energySum);
}

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
      if (fudgeFactor == 1.0) {
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
