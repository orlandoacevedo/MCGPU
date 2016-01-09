#include "SimBox.h"

void SimBox::keepMoleculeInBox(int molIdx) {
  int start = moleculeData[MOL_START][molIdx];
  int end = start + moleculeData[MOL_LEN][molIdx];
  int pIdx = primaryIndexes[moleculeData[MOL_PIDX_START][molIdx]];

  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    if (atomCoordinates[i][pIdx] < 0) {
      for (int j = start; j < end; j++) {
        atomCoordinates[i][j] += size[i];
      }
    } else if (atomCoordinates[i][pIdx] > size[i]) {
      for (int j = start; j < end; j++) {
        atomCoordinates[i][j] -= size[i];
      }
    }
  }
}

bool SimBox::addMolecules(Molecule* molecules, int nMolecules) {
  int nAtoms = 0, largestMolecule = 0;

  for (int i = 0; i < nMolecules; i++) {
    nAtoms += molecules[i].numOfAtoms;
    if (nAtoms > largestMolecule) {
      largestMolecule = nAtoms;
    }
  }

  numAtoms = nAtoms;
  numMolecules = nMolecules;

  rollBackCoordinates = new Real*[NUM_DIMENSIONS];
  atomCoordinates = new Real*[NUM_DIMENSIONS];
  atomData = new Real*[ATOM_DATA_SIZE];
  moleculeData = new int*[MOL_DATA_SIZE];

  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    atomCoordinates[i] = new Real[numAtoms];
    rollBackCoordinates[i] = new Real[largestMolecule];
  }
  for (int i = 0; i < ATOM_DATA_SIZE; i++) {
    atomData[i] = new Real[numAtoms];
  }
  for (int i = 0; i < MOL_DATA_SIZE; i++) {
    moleculeData[i] = new int[numMolecules];
  }

  int atomIdx = 0;
  for (int i = 0; i < numMolecules; i++) {

    moleculeData[MOL_START][i] = atomIdx;
    moleculeData[MOL_LEN][i] = molecules[i].numOfAtoms;
    moleculeData[MOL_TYPE][i] = molecules[i].type;

    for (int j = 0; j < molecules[i].numOfAtoms; j++) {
      Atom a = molecules[i].atoms[j];
      atomData[ATOM_SIGMA][atomIdx] = a.sigma;
      atomData[ATOM_EPSILON][atomIdx] = a.epsilon;
      atomData[ATOM_CHARGE][atomIdx] = a.charge;
      atomCoordinates[X_COORD][atomIdx] = a.x;
      atomCoordinates[Y_COORD][atomIdx] = a.y;
      atomCoordinates[Z_COORD][atomIdx] = a.z;
      atomIdx++;
    }
  }
  return true;
}


void SimBox::addPrimaryIndexes(std::vector<std::vector<int>* >* in) {
  int numPIdxes = 0;
  for (int i = 0; i < numMolecules; i++) {
    numPIdxes += in->at(moleculeData[MOL_TYPE][i])->size();
  }
  primaryIndexes = new int[numPIdxes];
  int idx = 0;
  for (int i = 0; i < numMolecules; i++) {
    vector<int>* v = in->at(moleculeData[MOL_TYPE][i]);
    moleculeData[MOL_PIDX_START][i] = idx;
    moleculeData[MOL_PIDX_COUNT][i] = v->size();
    for (int j = 0; j < v->size(); j++) {
      primaryIndexes[idx++] = v->at(j) + moleculeData[MOL_START][i];
    }
  }
}

void SimBox::changeMolecule(int molIdx, int vIdx, Real dX, Real dY, Real dZ, Real rX, Real rY, Real rZ) {
  int molStart = moleculeData[MOL_START][molIdx];
  int molLen = moleculeData[MOL_LEN][molIdx];

  for (int i = 0; i < molLen; i++) {
    for (int j = 0; j < NUM_DIMENSIONS; j++) {
      rollBackCoordinates[j][i] = atomCoordinates[j][molStart+i];
    }
    if (i == vIdx)
      continue;
    rotateAtom(molStart + i, molStart + vIdx, rX, rY, rZ);
    translateAtom(molStart + i, dX, dY, dZ);
  }
  translateAtom(molStart + vIdx, dX, dY, dZ);
  keepMoleculeInBox(molIdx);

  if (useNLC) {
    updateNLC(molIdx);
  }
}

void SimBox::changeMolecule(int molIdx) {
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

  // do the move here
  for (int i = 0; i < molLen; i++) {
    for (int j = 0; j < NUM_DIMENSIONS; j++) {
      rollBackCoordinates[j][i] = atomCoordinates[j][molStart + i];
    }
    if (i == vertexIdx)
      continue;
    rotateAtom(molStart + i, molStart + vertexIdx, rotX, rotY, rotZ);
    translateAtom(molStart + i, deltaX, deltaY, deltaZ);
  }
  translateAtom(molStart + vertexIdx, deltaX, deltaY, deltaZ);

  keepMoleculeInBox(molIdx);

  if (useNLC) {
    updateNLC(molIdx);
  }

}

void SimBox::rollback(int molIdx) {

  int molStart = moleculeData[MOL_START][molIdx];
  int molLen = moleculeData[MOL_LEN][molIdx];

  for (int i = 0; i < molLen; i++) {
    for (int j = 0; j < NUM_DIMENSIONS; j++) {
      atomCoordinates[j][molStart + i] = rollBackCoordinates[j][i];
    }
  }

  if (useNLC) {
    updateNLC(molIdx);
  }
}

void SimBox::translateAtom(int aIdx, Real dX, Real dY, Real dZ) {
  atomCoordinates[X_COORD][aIdx] += dX;
  atomCoordinates[Y_COORD][aIdx] += dY;
  atomCoordinates[Z_COORD][aIdx] += dZ;
}

void SimBox::rotateAtom(int aIdx, int pivotIdx, Real rotX, Real rotY, Real rotZ) {
  Real pX = atomCoordinates[X_COORD][pivotIdx];
  Real pY = atomCoordinates[Y_COORD][pivotIdx];
  Real pZ = atomCoordinates[Z_COORD][pivotIdx];

  translateAtom(aIdx, -pX, -pY, -pZ);
  rotateX(aIdx, rotX);
  rotateY(aIdx, rotY);
  rotateZ(aIdx, rotZ);
  translateAtom(aIdx, pX, pY, pZ);
}

void SimBox::rotateX(int aIdx, Real angleDeg) {
  Real angleRad = degreesToRadians(angleDeg);
  Real oldY = atomCoordinates[Y_COORD][aIdx];
  Real oldZ = atomCoordinates[Z_COORD][aIdx];
  atomCoordinates[Y_COORD][aIdx] = oldY * cos(angleRad) + oldZ * sin(angleRad);
  atomCoordinates[Z_COORD][aIdx] = oldZ * cos(angleRad) - oldY * sin(angleRad);
}

void SimBox::rotateY(int aIdx, Real angleDeg) {
  Real angleRad = degreesToRadians(angleDeg);
  Real oldZ = atomCoordinates[Z_COORD][aIdx];
  Real oldX = atomCoordinates[X_COORD][aIdx];
  atomCoordinates[Z_COORD][aIdx] = oldZ * cos(angleRad) + oldX * sin(angleRad);
  atomCoordinates[X_COORD][aIdx] = oldX * cos(angleRad) - oldZ * sin(angleRad);
}

void SimBox::rotateZ(int aIdx, Real angleDeg) {
  Real angleRad = degreesToRadians(angleDeg);
  Real oldX = atomCoordinates[X_COORD][aIdx];
  Real oldY = atomCoordinates[Y_COORD][aIdx];
  atomCoordinates[X_COORD][aIdx] = oldX * cos(angleRad) + oldY * sin(angleRad);
  atomCoordinates[Y_COORD][aIdx] = oldY * cos(angleRad) - oldX * sin(angleRad);
}

void SimBox::buildBox(Box* box) {
  size = new Real[NUM_DIMENSIONS];
  size[X_COORD] = box->environment->x;
  size[Y_COORD] = box->environment->y;
  size[Z_COORD] = box->environment->z;
  cutoff = box->environment->cutoff;
  temperature = box->environment->temp;
  maxTranslate = box->environment->maxTranslation;
  maxRotate = box->environment->maxRotation;
  numAtoms = box->environment->numOfAtoms;
  numMolecules = box->environment->numOfMolecules;
  addMolecules(box->molecules, numMolecules);
  addPrimaryIndexes(box->environment->primaryAtomIndexArray);
  if (useNLC) {
    fillNLC();
  }
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
  Real total = 0;

  const int p1Start = moleculeData[MOL_PIDX_START][currMol];
  const int p1End   = moleculeData[MOL_PIDX_COUNT][currMol] + p1Start;
  const int p1 = primaryIndexes[p1Start];
  int p2Start, p2End;

  if (useNLC) {
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
            total += calcMoleculeInteractionEnergy(subLJ, subCharge, currMol, idx);
          }
        }
        head = head->next;
      }
    }

  } else {

    for (int otherMol = startMol; otherMol < numMolecules; otherMol++) {
      if (otherMol != currMol) {
        p2Start = moleculeData[MOL_PIDX_START][otherMol];
        p2End = moleculeData[MOL_PIDX_COUNT][otherMol] + p2Start;
        if (moleculesInRange(p1Start, p1End, p2Start, p2End)) {
          total += calcMoleculeInteractionEnergy(subLJ, subCharge, currMol, otherMol);
        }
      }
    }
  }

  return total;
}

bool SimBox::moleculesInRange(refInt p1Start, refInt p1End, refInt p2Start, refInt p2End) {

  int p1;
  int p2;

  for (int p1Idx = p1Start; p1Idx < p1End; p1Idx++) {
    p1 = primaryIndexes[p1Idx];
    for (int p2Idx = p2Start; p2Idx < p2End; p2Idx++) {
      p2 = primaryIndexes[p2Idx];
    //  std::cout << p1 << " and " << p2 << ": dist is " << calcAtomDistSquared(box, p1, p2);
      if (calcAtomDistSquared(p1, p2) <= cutoff * cutoff) {
        return true;
      }
    }
  }

  return false;
}

Real SimBox::calcMoleculeInteractionEnergy (Real &subLJ, Real &subCharge, refInt m1, refInt m2) {
  Real ljSum = 0, chargeSum = 0;

  const int m1Start = moleculeData[MOL_START][m1];
  const int m1End = moleculeData[MOL_LEN][m1] + m1Start;

  const int m2Start = moleculeData[MOL_START][m2];
  const int m2End = moleculeData[MOL_LEN][m2] + m2Start;

  for (int i = m1Start; i < m1End; i++) {
    for (int j = m2Start; j < m2End; j++) {
      if (atomData[ATOM_SIGMA][i] >= 0 && atomData[ATOM_SIGMA][j] >= 0
          && atomData[ATOM_EPSILON][i] >= 0 && atomData[ATOM_EPSILON][j] >= 0) {

        const Real r2 = calcAtomDistSquared(i, j);
        ljSum += calcLJEnergy(i, j, r2);
        chargeSum += calcChargeEnergy(i, j, sqrt(r2));
      }
    }
  }

  subLJ += ljSum;
  subCharge += chargeSum;

  return (ljSum + chargeSum);
}

Real SimBox::calcAtomDistSquared(refInt a1, refInt a2) {
  Real dx = makePeriodic(atomCoordinates[X_COORD][a2] - atomCoordinates[X_COORD][a1], X_COORD);
  Real dy = makePeriodic(atomCoordinates[Y_COORD][a2] - atomCoordinates[Y_COORD][a1], Y_COORD);
  Real dz = makePeriodic(atomCoordinates[Z_COORD][a2] - atomCoordinates[Z_COORD][a1], Z_COORD);

  return dx * dx + dy * dy + dz * dz;
}

Real SimBox::makePeriodic(Real x, int dimension) {
  Real dimLength = size[dimension];
  if (x < -0.5 * dimLength) {
    x += dimLength;
  } else if (x > 0.5 * dimLength) {
    x -= dimLength;
  }
  return x;
}

Real SimBox::calcLJEnergy(refInt a1, refInt a2, const Real& r2) {

  if (r2 == 0.0) {
    return 0.0;
  } else {

    const Real sigma = calcBlending(atomData[ATOM_SIGMA][a1], atomData[ATOM_SIGMA][a2]);
    const Real epsilon = calcBlending(atomData[ATOM_EPSILON][a1], atomData[ATOM_EPSILON][a2]);

    const Real s2r2 = pow(sigma, 2) / r2;
    const Real s6r6 = pow(s2r2, 3);
    const Real s12r12 = pow(s6r6, 2);
    return 4.0 * epsilon * (s12r12 - s6r6);
  }
}

Real SimBox::calcChargeEnergy(refInt a1, refInt a2, const Real& r) {
  if (r == 0.0) {
    return 0.0;
  } else {
    const Real e = 332.06;
    return (atomData[ATOM_CHARGE][a1] * atomData[ATOM_CHARGE][a2] * e) / r;
  }
}

Real SimBox::calcBlending (const Real &a, const Real &b) {
  return sqrt(abs(a*b));
}

void SimBox::fillNLC() {
  neighbors = new NLC_Node*[27];
  numCells = new int[NUM_DIMENSIONS];
  cellWidth = new Real[NUM_DIMENSIONS];
  prevCell = new int[NUM_DIMENSIONS];
  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    numCells[i] = (int) (size[i] / cutoff);
    if (numCells[i] == 0)
      numCells[i] = 1;
    cellWidth[i] = size[i] / numCells[i];
  }
  neighborCells = new NLC_Node***[numCells[0]];
  for (int i = 0; i < numCells[0]; i++) {
    neighborCells[i] = new NLC_Node**[numCells[1]];
    for (int j = 0; j < numCells[1]; j++) {
      neighborCells[i][j] = new NLC_Node*[numCells[2]];
      for (int k = 0; k < numCells[2]; k++) {
        neighborCells[i][j][k] = new NLC_Node();
        neighborCells[i][j][k]->index = -1;
        neighborCells[i][j][k]->next = NULL;
      }
    }
  }
  nlc_heap = new NLC_Node[numMolecules];
  for (int i = 0; i < numMolecules; i++) {
    int pIdx = primaryIndexes[moleculeData[MOL_PIDX_START][i]];
    int cloc[3];
    for (int j = 0; j < NUM_DIMENSIONS; j++) {
      cloc[j] = getCell(atomCoordinates[j][pIdx], j);
    }
    nlc_heap[i].next = neighborCells[cloc[0]][cloc[1]][cloc[2]];
    nlc_heap[i].index = i;
    neighborCells[cloc[0]][cloc[1]][cloc[2]] = &nlc_heap[i];
  }
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
