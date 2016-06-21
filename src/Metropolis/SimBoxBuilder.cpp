#include "SimBoxBuilder.h"


SimBoxBuilder::SimBoxBuilder(bool useNLC, SBScanner* sbData_in) {
  sb = new SimBox();
  sb->useNLC = useNLC;
  sbData = sbData_in;
}

SimBox* SimBoxBuilder::build(Box* box) {
  // TODO (blm): Make this command line/config file configurable
  sb->maxIntraMoves = 15;
  initEnvironment(box->environment);
  addMolecules(box->molecules, box->environment->primaryAtomIndexArray->size());
  addPrimaryIndexes(box->environment->primaryAtomIndexArray);
  if (sb->useNLC) {
    fillNLC();
  }
  return sb;
}

void SimBoxBuilder::initEnvironment(Environment* environment) {
  sb->size = new Real[NUM_DIMENSIONS];
  sb->size[X_COORD] = environment->x;
  sb->size[Y_COORD] = environment->y;
  sb->size[Z_COORD] = environment->z;
  sb->cutoff = environment->cutoff;
  sb->temperature = environment->temp;
  sb->kT = sb->temperature * kBoltz;
  sb->maxTranslate = environment->maxTranslation;
  sb->maxRotate = environment->maxRotation;
  sb->numAtoms = environment->numOfAtoms;
  sb->numMolecules = environment->numOfMolecules;
  sb->maxBondDelta = environment->maxBondDelta;
  sb->maxAngleDelta = environment->maxAngleDelta;

  sb->stepNum = 0;
  sb->numBondMoves = 0;
  sb->numAcceptedBondMoves = 0;
  sb->numAngleMoves = 0;
  sb->numAcceptedAngleMoves = 0;
}

void SimBoxBuilder::addMolecules(Molecule* molecules, int numTypes) {
  int largestMolecule = 0, nAtoms = 0;
  int mostBonds = 0, nBonds = 0;
  int mostAngles = 0, nAngles = 0;

  sb->excludeAtoms = new int**[numTypes];
  sb->fudgeAtoms =  new int**[numTypes];

  for (int i = 0; i < numTypes; i++) {
    sb->excludeAtoms[i] = NULL;
    sb->fudgeAtoms[i] = NULL;
  }

  for (int i = 0; i < sb->numMolecules; i++) {
    nAtoms += molecules[i].numOfAtoms;
    nBonds += molecules[i].numOfBonds;
    nAngles += molecules[i].numOfAngles;
    if (molecules[i].numOfAtoms > largestMolecule) {
      largestMolecule = molecules[i].numOfAtoms;
    }
    if (molecules[i].numOfBonds > mostBonds) {
      mostBonds = molecules[i].numOfBonds;
    }
    if (molecules[i].numOfAngles > mostAngles) {
      mostAngles = molecules[i].numOfAngles;
    }
  }

  sb->numAtoms = nAtoms;
  sb->numBonds = nBonds;
  sb->numAngles = nAngles;

  sb->rollBackCoordinates = new Real*[NUM_DIMENSIONS];
  sb->atomCoordinates = new Real*[NUM_DIMENSIONS];
  sb->atomData = new Real*[ATOM_DATA_SIZE];
  sb->moleculeData = new int*[MOL_DATA_SIZE];
  sb->bondData = new Real*[BOND_DATA_SIZE];
  sb->angleData = new Real*[ANGLE_DATA_SIZE];
  sb->bondLengths = new Real[nBonds];
  sb->rollBackBondLengths = new Real[nBonds];
  sb->unionFindParent = new int[largestMolecule];
  sb->largestMol = largestMolecule;
  sb->angleSizes = new Real[nAngles];
  sb->rollBackAngleSizes = new Real[nAngles];

  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    sb->atomCoordinates[i] = new Real[sb->numAtoms];
    sb->rollBackCoordinates[i] = new Real[largestMolecule];
  }

  for (int i = 0; i < ATOM_DATA_SIZE; i++) {
    sb->atomData[i] = new Real[sb->numAtoms];
  }

  for (int i = 0; i < MOL_DATA_SIZE; i++) {
    sb->moleculeData[i] = new int[sb->numMolecules];
  }

  for (int i = 0; i < BOND_DATA_SIZE; i++) {
    sb->bondData[i] = new Real[sb->numBonds];
  }

  for (int i = 0; i < ANGLE_DATA_SIZE; i++) {
    sb->angleData[i] = new Real[sb->numAngles];
  }
  int atomIdx = 0, bondIdx = 0, angleIdx = 0;

  std::map<int, int> idToIdx;
  std::map<int, std::string*> idToName;

  for (int i = 0; i < sb->numMolecules; i++) {
    sb->moleculeData[MOL_START][i] = atomIdx;
    sb->moleculeData[MOL_LEN][i] = molecules[i].numOfAtoms;
    sb->moleculeData[MOL_TYPE][i] = molecules[i].type;
    sb->moleculeData[MOL_BOND_START][i] = bondIdx;
    sb->moleculeData[MOL_BOND_COUNT][i] = molecules[i].numOfBonds;
    sb->moleculeData[MOL_ANGLE_START][i] = angleIdx;
    sb->moleculeData[MOL_ANGLE_COUNT][i] = molecules[i].numOfAngles;

    // Store all the atom data for the molecule
    for (int j = 0; j < molecules[i].numOfAtoms; j++) {
      Atom a = molecules[i].atoms[j];
      idToIdx[a.id] = atomIdx;
      idToName[a.id] = a.name;
      sb->atomData[ATOM_SIGMA][atomIdx] = a.sigma;
      sb->atomData[ATOM_EPSILON][atomIdx] = a.epsilon;
      sb->atomData[ATOM_CHARGE][atomIdx] = a.charge;
      sb->atomCoordinates[X_COORD][atomIdx] = a.x;
      sb->atomCoordinates[Y_COORD][atomIdx] = a.y;
      sb->atomCoordinates[Z_COORD][atomIdx] = a.z;
      atomIdx++;
    }

    // Store all the bond data for the molecule
    for (int j = 0; j < molecules[i].numOfBonds; j++) {
      Bond b = molecules[i].bonds[j];
      sb->bondData[BOND_A1_IDX][bondIdx] = idToIdx[b.atom1];
      sb->bondData[BOND_A2_IDX][bondIdx] = idToIdx[b.atom2];
      std::string name1 = *(idToName[b.atom1]);
      std::string name2 = *(idToName[b.atom2]);
      sb->bondData[BOND_KBOND][bondIdx] = sbData->getKBond(name1, name2);
      sb->bondData[BOND_EQDIST][bondIdx] = sbData->getEqBondDist(name1, name2);
      sb->bondLengths[bondIdx] = b.distance;
      sb->bondData[BOND_VARIABLE][bondIdx] = b.variable;
      bondIdx++;
    }

    // Store all the angle data for the molecule
    for (int j = 0; j < molecules[i].numOfAngles; j++) {
      Angle a = molecules[i].angles[j];
      sb->angleData[ANGLE_A1_IDX][angleIdx] = idToIdx[a.atom1];
      sb->angleData[ANGLE_A2_IDX][angleIdx] = idToIdx[a.atom2];
      sb->angleData[ANGLE_MID_IDX][angleIdx] = idToIdx[a.commonAtom];
      std::string a1Name = *(idToName[a.atom1]);
      std::string a2Name = *(idToName[a.atom2]);
      std::string midName = *(idToName[a.commonAtom]);
      sb->angleData[ANGLE_KANGLE][angleIdx] = sbData->getKAngle(a1Name, midName, a2Name);
      sb->angleData[ANGLE_EQANGLE][angleIdx] = sbData->getEqAngle(a1Name, midName, a2Name);
      sb->angleSizes[angleIdx] = a.value;
      sb->angleData[ANGLE_VARIABLE][angleIdx] = a.variable;
      angleIdx++;
    }

    int type = molecules[i].type;
    if (sb->excludeAtoms[type] == NULL) {
      int numOfAtoms = sb->moleculeData[MOL_LEN][i];
      sb->excludeAtoms[type] = new int*[numOfAtoms];
      sb->fudgeAtoms[type] = new int*[numOfAtoms];
      int startIdx = sb->moleculeData[MOL_START][i];
      int *excludeCount = new int[numOfAtoms];
      int *fudgeCount = new int[numOfAtoms];
      for (int j = 0; j < numOfAtoms; j++) {
        excludeCount[j] = 0;
        fudgeCount[j] = 0;
      }

      // Count the number of atoms that will be excluded
      for (int j = 0; j < molecules[i].numOfBonds; j++) {
        int idx1 = idToIdx[molecules[i].bonds[j].atom1] - startIdx;
        int idx2 = idToIdx[molecules[i].bonds[j].atom2] - startIdx;
        if (idx1 >= 0 && idx1 < numOfAtoms && idx2 >= 0 && idx2 < numOfAtoms) {
          excludeCount[idx1]++;
          excludeCount[idx2]++;
        }
      }
      for (int j = 0; j < molecules[i].numOfAngles; j++) {
        int idx1 = idToIdx[molecules[i].angles[j].atom1] - startIdx;
        int idx2 = idToIdx[molecules[i].angles[j].atom2] - startIdx;
        if (idx1 >= 0 && idx1 < numOfAtoms && idx2 >= 0 && idx2 < numOfAtoms) {
          excludeCount[idx1]++;
          excludeCount[idx2]++;
        }
      }

      // Count the number of atoms that will be fudged
      for (int j = 0; j < molecules[i].numOfHops; j++) {
        int idx1 = idToIdx[molecules[i].hops[j].atom1] - startIdx;
        int idx2 = idToIdx[molecules[i].hops[j].atom2] - startIdx;
        int hopDist = molecules[i].hops[j].hop;
        if (idx1 >= 0 && idx1 < numOfAtoms && idx2 >= 0 && idx2 < numOfAtoms && hopDist == 3) {
          fudgeCount[idx1]++;
          fudgeCount[idx2]++;
        }
      }

      // Build the exclusion and fudge matrix
      for (int j = 0; j < numOfAtoms; j++) {
        sb->excludeAtoms[type][j] = new int[excludeCount[j] + 1];
        sb->fudgeAtoms[type][j] = new int[fudgeCount[j] + 1];
        excludeCount[j] = -1;
        fudgeCount[j] = -1;
      }

      // Exclude two atoms if they're joined by a bond
      for (int j = 0; j < molecules[i].numOfBonds; j++) {
        int idx1 = idToIdx[molecules[i].bonds[j].atom1] - startIdx;
        int idx2 = idToIdx[molecules[i].bonds[j].atom2] - startIdx;
        if (idx1 >= 0 && idx1 < numOfAtoms && idx2 >= 0 && idx2 < numOfAtoms) {
          sb->excludeAtoms[type][idx1][++excludeCount[idx1]] = idx2;
          sb->excludeAtoms[type][idx2][++excludeCount[idx2]] = idx1;
        }
      }

      // Exclude two atoms if they're joined by an angle
      for (int j = 0; j < molecules[i].numOfAngles; j++) {
        int idx1 = idToIdx[molecules[i].angles[j].atom1] - startIdx;
        int idx2 = idToIdx[molecules[i].angles[j].atom2] - startIdx;
        if (idx1 >= 0 && idx1 < numOfAtoms && idx2 >= 0 && idx2 < numOfAtoms) {
          sb->excludeAtoms[type][idx1][++excludeCount[idx1]] = idx2;
          sb->excludeAtoms[type][idx2][++excludeCount[idx2]] = idx1;
        }
      }

      // Fudge two atoms if they're separated by exactly 3 hops
      for (int j = 0; j < molecules[i].numOfHops; j++) {
        int idx1 = idToIdx[molecules[i].hops[j].atom1] - startIdx;
        int idx2 = idToIdx[molecules[i].hops[j].atom2] - startIdx;
        int hopDist = molecules[i].hops[j].hop;
        if (idx1 >= 0 && idx1 < numOfAtoms && idx2 >= 0 && idx2 < numOfAtoms && hopDist == 3) {
          sb->fudgeAtoms[type][idx1][++fudgeCount[idx1]] = idx2;
          sb->fudgeAtoms[type][idx2][++fudgeCount[idx2]] = idx1;
        }
      }

      // End the exclude and fudge list with a -1 flag
      for (int j = 0; j < numOfAtoms; j++) {
        sb->excludeAtoms[type][j][++excludeCount[j]] = -1;
        sb->fudgeAtoms[type][j][++fudgeCount[j]] = -1;
      }

      delete[] excludeCount;
      delete[] fudgeCount;
    }
  }
}

void SimBoxBuilder::addPrimaryIndexes(std::vector< std::vector<int>* >* in) {
  int numPIdxes = 0;
  for (int i = 0; i < sb->numMolecules; i++) {
    numPIdxes += in->at(sb->moleculeData[MOL_TYPE][i])->size();
  }
  sb->primaryIndexes = new int[numPIdxes];
  sb->numPIdxes = numPIdxes;

  int idx = 0;
  for (int i = 0; i < sb->numMolecules; i++) {
    vector<int>* v = in->at(sb->moleculeData[MOL_TYPE][i]);
    sb->moleculeData[MOL_PIDX_START][i] = idx;
    sb->moleculeData[MOL_PIDX_COUNT][i] = v->size();
    for (int j = 0; j < v->size(); j++) {
      sb->primaryIndexes[idx++] = v->at(j) + sb->moleculeData[MOL_START][i];
    }
  }
}

void SimBoxBuilder::fillNLC() {
  sb->neighbors = new NLC_Node*[27];
  sb->numCells = new int[NUM_DIMENSIONS];
  sb->cellWidth = new Real[NUM_DIMENSIONS];
  sb->prevCell = new int[NUM_DIMENSIONS];
  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    sb->numCells[i] = (int) (sb->size[i] / sb->cutoff);
    if (sb->numCells[i] == 0) {
      sb->numCells[i] = 1;
    }
    sb->cellWidth[i] = sb->size[i] / sb->numCells[i];
  }

  sb->neighborCells = new NLC_Node***[sb->numCells[0]];
  for (int i = 0; i < sb->numCells[0]; i++) {
    sb->neighborCells[i] = new NLC_Node**[sb->numCells[1]];
    for (int j = 0; j < sb->numCells[1]; j++) {
      sb->neighborCells[i][j] = new NLC_Node*[sb->numCells[2]];
      for (int k = 0; k < sb->numCells[2]; k++) {
        sb->neighborCells[i][j][k] = new NLC_Node();
        sb->neighborCells[i][j][k]->index = -1;
        sb->neighborCells[i][j][k]->next = NULL;
      }
    }
  }
  sb->nlc_heap = new NLC_Node[sb->numMolecules];
  for (int i = 0; i < sb->numMolecules; i++) {
    int pIdx = sb->primaryIndexes[sb->moleculeData[MOL_PIDX_START][i]];
    int cloc[3];
    for (int j = 0; j < NUM_DIMENSIONS; j++) {
      cloc[j] = sb->getCell(sb->atomCoordinates[j][pIdx], j);
    }
    sb->nlc_heap[i].next = sb->neighborCells[cloc[0]][cloc[1]][cloc[2]];
    sb->nlc_heap[i].index = i;
    sb->neighborCells[cloc[0]][cloc[1]][cloc[2]] = &(sb->nlc_heap[i]);
  }
}
