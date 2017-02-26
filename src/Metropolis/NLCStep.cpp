#include "BruteForceStep.h"
#include "NLCStep.h"
#include "SimulationStep.h"
#include "GPUCopy.h"

#ifdef _OPENACC
#include <openacc.h>
#endif


NLCStep::~NLCStep() {
  NLCCalcs::freeNLC(this->nlc);
  this->nlc = NULL;
}

Real NLCStep::calcMolecularEnergyContribution(int currMol, int startMol) {
  return NLCCalcs::calcMolecularEnergyContribution( currMol,
      startMol, this->nlc);
}

Real NLCStep::calcSystemEnergy(Real &subLJ, Real &subCharge, int numMolecules) {
  Real result = SimulationStep::calcSystemEnergy(subLJ, subCharge,
                                                 numMolecules);
  this->nlc = NLCCalcs::createNLC();
  return result;
}

void NLCStep::changeMolecule(int molIdx, SimBox *box) {
  SimulationStep::changeMolecule(molIdx, box);
  NLCCalcs::updateNLC(molIdx, this->nlc);
}

void NLCStep::rollback(int molIdx, SimBox *box) {
  SimulationStep::rollback(molIdx, box);
  NLCCalcs::updateNLC(molIdx, this->nlc);
}

// ----- NLCCalcs Definitions -----

Real NLCCalcs::calcMolecularEnergyContribution(
    int currMol, int startMol, NLC* nlc) {
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

  if (nlc == NULL) {
    for (int otherMol = startMol; otherMol < numMolecules; otherMol++) {
      if (otherMol != currMol) {
        int p2Start = molData[MOL_PIDX_START][otherMol];
        int p2End = molData[MOL_PIDX_COUNT][otherMol] + p2Start;
        if (SimCalcs::moleculesInRange(p1Start, p1End, p2Start, p2End,
                                       atomCoords, bSize, pIdxes, cutoff)) {
          total += calcMoleculeInteractionEnergy(currMol, otherMol, molData,
                                                 aData, atomCoords, bSize);
        }
      }
    }
  } else {
    int nNeighbors = NLCCalcs::findNeighbors(startMol, nlc);
    for (int i = 0; i < nNeighbors; i++) {
      NLC_Node* node = nlc->neighbors[i];
      while(node->next != NULL) {
        if (node->next->index == currMol) nlc->prevNode = node;
        int otherMol = node->index;
        if (otherMol != currMol) {
          int p2Start = molData[MOL_PIDX_START][otherMol];
          int p2End = molData[MOL_PIDX_COUNT][otherMol] + p2Start;
          if (SimCalcs::moleculesInRange(p1Start, p1End, p2Start, p2End,
                                         atomCoords, bSize, pIdxes, cutoff)) {
            total += calcMoleculeInteractionEnergy(currMol, otherMol, molData,
                                                   aData, atomCoords, bSize);
          }
        }
        node = node->next;
      }
    }
  }

  return total;
}

// TODO: Duplicate; abstract out when PGCC supports it
Real NLCCalcs::calcMoleculeInteractionEnergy (int m1, int m2,
                                                          int** molData,
                                                          Real** aData,
                                                          Real** aCoords,
                                                          Real* bSize) {
  Real energySum = 0;

  const int m1Start = molData[MOL_START][m1];
  const int m1End = molData[MOL_LEN][m1] + m1Start;

  const int m2Start = molData[MOL_START][m2];
  const int m2End = molData[MOL_LEN][m2] + m2Start;

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

NLC* NLCCalcs::createNLC() {

  const Real cutoff = SimCalcs::sb->cutoff;
  const int numMolecules = SimCalcs::sb->numMolecules;

  Real* bSize = GPUCopy::sizePtr();
  int* pIdxes = GPUCopy::primaryIndexesPtr();
  int** molData = GPUCopy::moleculeDataPtr();
  Real** atomCoords = GPUCopy::atomCoordinatesPtr();


  NLC* simulation_NLC = new NLC;
  simulation_NLC->neighbors = new NLC_Node*[27];
  simulation_NLC->numCells = new int[NUM_DIMENSIONS];
  simulation_NLC->cellWidth = new Real[NUM_DIMENSIONS];
  simulation_NLC->prevCell = new int[NUM_DIMENSIONS];

  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    simulation_NLC->numCells[i] = (int) (bSize[i] / cutoff);
    if (simulation_NLC->numCells[i] == 0) {
      simulation_NLC->numCells[i] = 1;
    }
    simulation_NLC->cellWidth[i] = bSize[i] / simulation_NLC->numCells[i];
  }

  simulation_NLC->neighborCells = new NLC_Node***[simulation_NLC->numCells[0]];
  for (int i = 0; i < simulation_NLC->numCells[0]; i++) {
    int nCells = simulation_NLC->numCells[1];
    simulation_NLC->neighborCells[i] = new NLC_Node**[nCells];

    for (int j = 0; j < simulation_NLC->numCells[1]; j++) {
      int nCells = simulation_NLC->numCells[2];
      simulation_NLC->neighborCells[i][j] = new NLC_Node*[nCells];
      for (int k = 0; k < nCells; k++) {
        simulation_NLC->neighborCells[i][j][k] = new NLC_Node();
        simulation_NLC->neighborCells[i][j][k]->index = -1;
        simulation_NLC->neighborCells[i][j][k]->next = NULL;
      }
    }
  }

  simulation_NLC->nlc_heap = new NLC_Node[numMolecules];
  for (int i = 0; i < numMolecules; i++) {
    int pIdx = pIdxes[molData[MOL_PIDX_START][i]];
    int cloc[3];
    for (int j = 0; j < NUM_DIMENSIONS; j++) {
      cloc[j] = NLCCalcs::getCell(atomCoords[j][pIdx], j, simulation_NLC);
    }

    simulation_NLC->nlc_heap[i].next = simulation_NLC->neighborCells[cloc[0]][cloc[1]][cloc[2]];
    simulation_NLC->nlc_heap[i].index = i;
    simulation_NLC->neighborCells[cloc[0]][cloc[1]][cloc[2]] = &(simulation_NLC->nlc_heap[i]);
  }

  return simulation_NLC;
}

void NLCCalcs::updateNLC(int molIdx, NLC * nlc) {

  int* pIdxes = GPUCopy::primaryIndexesPtr();
  int** molData = GPUCopy::moleculeDataPtr();
  Real** atomCoords = GPUCopy::atomCoordinatesPtr();

  bool update = false;
  int newCell[3];
  int pIdx = pIdxes[molData[MOL_PIDX_START][molIdx]];
  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    newCell[i] = getCell(atomCoords[i][pIdx], i, nlc);
    if (newCell[i] != nlc->prevCell[i]) {
      update = true;
    }
  }

  if (update) {
    if (nlc->prevNode->index == molIdx) {
      nlc->neighborCells[nlc->prevCell[0]][nlc->prevCell[1]][nlc->prevCell[2]] = (nlc->prevNode->next);
      nlc->prevNode->next = nlc->neighborCells[newCell[0]][newCell[1]][newCell[2]];
      nlc->neighborCells[newCell[0]][newCell[1]][newCell[2]] = nlc->prevNode;
    } else {
      NLC_Node* removeNode = nlc->prevNode->next;
      nlc->prevNode->next = (nlc->prevNode->next->next);
      removeNode->next = (nlc->neighborCells[newCell[0]][newCell[1]][newCell[2]]);
      nlc->neighborCells[newCell[0]][newCell[1]][newCell[2]] = removeNode;
    }
  }
}

int NLCCalcs::findNeighbors(int molIdx, NLC * nlc) {

  int* pIdxes = GPUCopy::primaryIndexesPtr();
  int** molData = GPUCopy::moleculeDataPtr();
  Real** atomCoords = GPUCopy::atomCoordinatesPtr();

  int outIdx = 0;
  int pIdx = pIdxes[molData[MOL_PIDX_START][molIdx]];

  int base[3];
  for (int i = 0; i < 3; i++) {
    base[i] = getCell(atomCoords[i][pIdx], i, nlc);
    nlc->prevCell[i] = base[i];
  }

  for (int i = -1; i <= 1; i++) {
    if (i + 1 >= nlc->numCells[0]) break;
    int c0 = wrapCell(base[0] + i, 0, nlc);
    for (int j = -1; j <= 1; j++) {
      if (j + 1 >= nlc->numCells[1]) break;
      int c1 = wrapCell(base[1] + j, 1, nlc);
      for (int k = -1; k <= 1; k++) {
        if (k + 1 >= nlc->numCells[2]) break;
        int c2 = wrapCell(base[2] + k, 2, nlc);
        nlc->neighbors[outIdx] = nlc->neighborCells[c0][c1][c2];
        outIdx++;
      }
    }
  }

  return outIdx;
}

int NLCCalcs::wrapCell(int idx, int dimension, NLC * nlc) {
  return (idx + nlc->numCells[dimension]) % nlc->numCells[dimension];
}

int NLCCalcs::getCell(Real loc, int dimension, NLC * nlc) {
  int idx = (int) (loc / nlc->cellWidth[dimension]);
  int numC = nlc->numCells[dimension];
  return (idx + numC) % numC;
}

void NLCCalcs::freeNLC(NLC *nlc) {
  free(nlc);
}
