#ifdef _OPENACC
#include <openacc.h>
#endif

#include "SimBox.h"
#include "GPUCopy.h"


// ----- Experimental -----

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
