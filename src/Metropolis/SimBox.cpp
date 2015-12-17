#include "SimBox.h"


void SimBox::rollBack(refInt molIdx, Real &dx, Real &dy, Real &dz) {
  int startAtomIdx = moleculeData[MOL_START][molIdx];
  int endAtomIdx = moleculeData[MOL_LEN][molIdx] + startAtomIdx;
  for (int idx = startAtomIdx; idx < endAtomIdx; idx++) {
    atomCoordinates[X_COORD][idx] = keepInBox(atomCoordinates[X_COORD][idx]-dx, X_COORD);
    atomCoordinates[Y_COORD][idx] = keepInBox(atomCoordinates[Y_COORD][idx]-dy, Y_COORD);
    atomCoordinates[Z_COORD][idx] = keepInBox(atomCoordinates[Z_COORD][idx]-dz, Z_COORD);
  }
}

void SimBox::moveMolecule(refInt molIdx, Real &dx, Real &dy, Real &dz) {
  dx = rand() / ((Real) RAND_MAX);
  dx = dx * 2 * maxTranslate - maxTranslate;
  dy = rand() / ((Real) RAND_MAX);
  dy = dy * 2 * maxTranslate - maxTranslate;
  dz = rand() / ((Real) RAND_MAX);
  dz = dy * 2 * maxTranslate - maxTranslate;

  int startAtomIdx = moleculeData[MOL_START][molIdx];
  int endAtomIdx = moleculeData[MOL_LEN][molIdx] + startAtomIdx;

  for (int idx = startAtomIdx; idx < endAtomIdx; idx++) {
    atomCoordinates[X_COORD][idx] = keepInBox(atomCoordinates[X_COORD][idx]+dx, X_COORD);
    atomCoordinates[Y_COORD][idx] = keepInBox(atomCoordinates[Y_COORD][idx]+dy, Y_COORD);
    atomCoordinates[Z_COORD][idx] = keepInBox(atomCoordinates[Z_COORD][idx]+dz, Z_COORD);
  }
}

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

void SimBox::printCoordinates() {
  for (int i = 0; i < numMolecules; i++) {
    int start = moleculeData[MOL_START][i];
    int end = start + moleculeData[MOL_LEN][i];
    std::cout << "MOL # " << (i+1) << std::endl;
    for (int j = start; j < end; j++) {
      for (int k = 0; k < NUM_DIMENSIONS; k++) {
        std::cout << atomCoordinates[k][j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
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

}

void SimBox::rollback(int molIdx) {

  int molStart = moleculeData[MOL_START][molIdx];
  int molLen = moleculeData[MOL_LEN][molIdx];

  for (int i = 0; i < molLen; i++) {
    for (int j = 0; j < NUM_DIMENSIONS; j++) {
      atomCoordinates[j][molStart + i] = rollBackCoordinates[j][i];
    }
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

bool SimBox::addAtoms(std::vector<Atom> atoms) {

  int numAtoms = atoms.size();

  atomCoordinates = new Real*[NUM_DIMENSIONS];
  atomData = new Real*[ATOM_DATA_SIZE];

  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    atomCoordinates[i] = new Real[numAtoms];
  }
  for (int i = 0; i < ATOM_DATA_SIZE; i++) {
    atomData[i] = new Real[numAtoms];
  }

  for (int i = 0; i < numAtoms; i++) {
    atomCoordinates[X_COORD][i] = atoms.at(i).x;
    atomCoordinates[Y_COORD][i] = atoms.at(i).y;
    atomCoordinates[Z_COORD][i] = atoms.at(i).z;
    atomData[ATOM_SIGMA][i] = atoms.at(i).sigma;
    atomData[ATOM_CHARGE][i] = atoms.at(i).charge;
    atomData[ATOM_EPSILON][i] = atoms.at(i).epsilon;
  }

  return true;
}

Real SimBox::keepInBox(const Real &value, const int dimension) {
  Real out = value;
  Real comparision = size[dimension];
  while (out < 0) {
    out += comparision;
  }
  while (out > comparision) {
    out -= comparision;
  }
  return out;
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
}

int SimBox::chooseMolecule() const {
  return (int) randomReal(0, numMolecules);
}

Real SimBox::calcSystemEnergy (Real &subLJ, Real &subCharge) {
  Real total = 0;

  for (int mol = 0; mol < numMolecules; mol++) {
    total += calcMolecularEnergyContribution(subLJ, subCharge, mol, mol);
  }

  return total;
}

Real SimBox::calcMolecularEnergyContribution(Real &subLJ, Real &subCharge, refInt currMol, refInt startMol) {
  Real total = 0;

  //std::cout << "Calculating energy for " << currMol << std::endl;


  const int p1Start = moleculeData[MOL_PIDX_START][currMol];
  const int p1End   = moleculeData[MOL_PIDX_COUNT][currMol] + p1Start;
  int p2Start, p2End;

  for (int otherMol = startMol; otherMol < numMolecules; otherMol++) {
    if (otherMol != currMol) {
      //std::cout << "Molecule # " << otherMol << " is ";
      p2Start = moleculeData[MOL_PIDX_START][otherMol];
      p2End = moleculeData[MOL_PIDX_COUNT][otherMol] + p2Start;
      if (moleculesInRange(p1Start, p1End, p2Start, p2End)) {
      //  std::cout << "in range." << std::endl;
        total += calcMoleculeInteractionEnergy(subLJ, subCharge, currMol, otherMol);
      } else {
      //  std::cout << "not in range." << std::endl;
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

  /*std::cout << "(" << box.atomCoordinates[box.X_COORD][a1] << ", " <<
    box.atomCoordinates[box.Y_COORD][a1] << ", " <<
    box.atomCoordinates[box.Z_COORD][a1] << ") and (" <<
    box.atomCoordinates[box.X_COORD][a2] << ", " <<
    box.atomCoordinates[box.Y_COORD][a2] << ", " <<
    box.atomCoordinates[box.Z_COORD][a2] << ")" << std::endl; */

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
