#include "SimBoxBuilder.h"


SimBoxBuilder::SimBoxBuilder(bool useNLC, SBScanner* sbData_in) {
  sb = new SimBox();
  sb->useNLC = useNLC;
  sbData = sbData_in;
}

SimBox* SimBoxBuilder::build(Box* box) {
  initEnvironment(box->environment);
  addMolecules(box->molecules);
  addPrimaryIndexes(box->environment->primaryAtomIndexArray);
  if (sb->useNLC) {
    fillNLC();
  }
  return sb;
}

void SimBoxBuilder::initEnvironment(Environment* environment) {
  sb->size = new Real[SimBox::NUM_DIMENSIONS];
  sb->size[SimBox::X_COORD] = environment->x;
  sb->size[SimBox::Y_COORD] = environment->y;
  sb->size[SimBox::Z_COORD] = environment->z;
  sb->cutoff = environment->cutoff;
  sb->temperature = environment->temp;
  sb->maxTranslate = environment->maxTranslation;
  sb->maxRotate = environment->maxRotation;
  sb->numAtoms = environment->numOfAtoms;
  sb->numMolecules = environment->numOfMolecules;
}

void SimBoxBuilder::addMolecules(Molecule* molecules) {
  int largestMolecule = 0, nAtoms = 0;
  int mostBonds = 0, nBonds = 0;
  int mostAngles = 0, nAngles = 0;

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

  sb->rollBackCoordinates = new Real*[SimBox::NUM_DIMENSIONS];
  sb->atomCoordinates = new Real*[SimBox::NUM_DIMENSIONS];
  sb->atomData = new Real*[SimBox::ATOM_DATA_SIZE];
  sb->moleculeData = new int*[SimBox::MOL_DATA_SIZE];
  sb->bondData = new Real*[SimBox::BOND_DATA_SIZE];
  sb->angleData = new Real*[SimBox::ANGLE_DATA_SIZE];
  sb->bondLengths = new Real[nBonds];
  sb->unionFindParent = new int[largestMolecule];
  sb->angleSizes = new Real[nAngles];

  for (int i = 0; i < SimBox::NUM_DIMENSIONS; i++) {
    sb->atomCoordinates[i] = new Real[sb->numAtoms];
    sb->rollBackCoordinates[i] = new Real[largestMolecule];
  }

  for (int i = 0; i < SimBox::ATOM_DATA_SIZE; i++) {
    sb->atomData[i] = new Real[sb->numAtoms];
  }

  for (int i = 0; i < SimBox::MOL_DATA_SIZE; i++) {
    sb->moleculeData[i] = new int[sb->numMolecules];
  }

  for (int i = 0; i < SimBox::BOND_DATA_SIZE; i++) {
    sb->bondData[i] = new Real[sb->numBonds];
  }

  for (int i = 0; i < SimBox::ANGLE_DATA_SIZE; i++) {
    sb->angleData[i] = new Real[sb->numAngles];
  }
  int atomIdx = 0, bondIdx = 0, angleIdx = 0;

  std::map<int, int> idToIdx;
  std::map<int, std::string*> idToName;

  for (int i = 0; i < sb->numMolecules; i++) {

    sb->moleculeData[SimBox::MOL_START][i] = atomIdx;
    sb->moleculeData[SimBox::MOL_LEN][i] = molecules[i].numOfAtoms;
    sb->moleculeData[SimBox::MOL_TYPE][i] = molecules[i].type;
    sb->moleculeData[SimBox::MOL_BOND_START][i] = bondIdx;
    sb->moleculeData[SimBox::MOL_BOND_COUNT][i] = molecules[i].numOfBonds;
    sb->moleculeData[SimBox::MOL_ANGLE_START][i] = angleIdx;
    sb->moleculeData[SimBox::MOL_ANGLE_COUNT][i] = molecules[i].numOfAngles;

    for (int j = 0; j < molecules[i].numOfAtoms; j++) {
      Atom a = molecules[i].atoms[j];
      idToIdx[a.id] = atomIdx;
      idToName[a.id] = a.name;
      sb->atomData[SimBox::ATOM_SIGMA][atomIdx] = a.sigma;
      sb->atomData[SimBox::ATOM_EPSILON][atomIdx] = a.epsilon;
      sb->atomData[SimBox::ATOM_CHARGE][atomIdx] = a.charge;
      sb->atomCoordinates[SimBox::X_COORD][atomIdx] = a.x;
      sb->atomCoordinates[SimBox::Y_COORD][atomIdx] = a.y;
      sb->atomCoordinates[SimBox::Z_COORD][atomIdx] = a.z;
      atomIdx++;
    }

    for (int j = 0; j < molecules[i].numOfBonds; j++) {
      Bond b = molecules[i].bonds[j];
      sb->bondData[SimBox::BOND_A1_IDX][bondIdx] = idToIdx[b.atom1];
      sb->bondData[SimBox::BOND_A2_IDX][bondIdx] = idToIdx[b.atom2];
      std::string name1 = *(idToName[b.atom1]);
      std::string name2 = *(idToName[b.atom2]);
      sb->bondData[SimBox::BOND_KBOND][bondIdx] = sbData->getKBond(name1, name2);
      sb->bondData[SimBox::BOND_EQDIST][bondIdx] = sbData->getEqBondDist(name1, name2);
      sb->bondLengths[bondIdx] = b.distance;
      sb->bondData[SimBox::BOND_VARIABLE][bondIdx] = b.variable;
      bondIdx++;
    }

    for (int j = 0; j < molecules[i].numOfAngles; j++) {
      Angle a = molecules[i].angles[j];
      sb->angleData[SimBox::ANGLE_A1_IDX][angleIdx] = idToIdx[a.atom1];
      sb->angleData[SimBox::ANGLE_A2_IDX][angleIdx] = idToIdx[a.atom2];
      sb->angleData[SimBox::ANGLE_MID_IDX][angleIdx] = idToIdx[a.commonAtom];
      if (a.commonAtom != -1) {
        std::string a1Name = *(idToName[a.atom1]);
        std::string a2Name = *(idToName[a.atom2]);
        std::string midName = *(idToName[a.commonAtom]);
        sb->angleData[SimBox::ANGLE_KANGLE][angleIdx] = sbData->getKAngle(a1Name, midName, a2Name);
        sb->angleData[SimBox::ANGLE_EQANGLE][angleIdx] = sbData->getEqAngle(a1Name, midName, a2Name);
      }
      sb->angleSizes[angleIdx] = a.value;
      sb->angleData[SimBox::ANGLE_VARIABLE][angleIdx] = a.variable;
      angleIdx++;
    }
  }
}

void SimBoxBuilder::addPrimaryIndexes(std::vector< std::vector<int>* >* in) {
  int numPIdxes = 0;
  for (int i = 0; i < sb->numMolecules; i++) {
    numPIdxes += in->at(sb->moleculeData[SimBox::MOL_TYPE][i])->size();
  }
  sb->primaryIndexes = new int[numPIdxes];

  int idx = 0;
  for (int i = 0; i < sb->numMolecules; i++) {
    vector<int>* v = in->at(sb->moleculeData[SimBox::MOL_TYPE][i]);
    sb->moleculeData[SimBox::MOL_PIDX_START][i] = idx;
    sb->moleculeData[SimBox::MOL_PIDX_COUNT][i] = v->size();
    for (int j = 0; j < v->size(); j++) {
      sb->primaryIndexes[idx++] = v->at(j) + sb->moleculeData[SimBox::MOL_START][i];
    }
  }
}

void SimBoxBuilder::fillNLC() {
  sb->neighbors = new NLC_Node*[27];
  sb->numCells = new int[SimBox::NUM_DIMENSIONS];
  sb->cellWidth = new Real[SimBox::NUM_DIMENSIONS];
  sb->prevCell = new int[SimBox::NUM_DIMENSIONS];
  for (int i = 0; i < SimBox::NUM_DIMENSIONS; i++) {
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
    int pIdx = sb->primaryIndexes[sb->moleculeData[SimBox::MOL_PIDX_START][i]];
    int cloc[3];
    for (int j = 0; j < SimBox::NUM_DIMENSIONS; j++) {
      cloc[j] = sb->getCell(sb->atomCoordinates[j][pIdx], j);
    }
    sb->nlc_heap[i].next = sb->neighborCells[cloc[0]][cloc[1]][cloc[2]];
    sb->nlc_heap[i].index = i;
    sb->neighborCells[cloc[0]][cloc[1]][cloc[2]] = &(sb->nlc_heap[i]);
  }
}
