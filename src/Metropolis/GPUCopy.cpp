#include "GPUCopy.h"
#include <openacc.h>

bool parallel = false;

Real**  h_atomData = NULL;
Real**  d_atomData = NULL;

Real** h_rollBackCoordinates = NULL;
Real** d_rollBackCoordinates = NULL;

Real** h_atomCoordinates = NULL;
Real** d_atomCoordinates = NULL;

int* h_primaryIndexes = NULL;
int* d_primaryIndexes = NULL;

int** h_moleculeData = NULL;
int** d_moleculeData = NULL;

Real*   h_size = NULL;
Real*   d_size = NULL;

Real** GPUCopy::atomDataPtr() {
  if (parallel) {
    return d_atomData;
  } else {
    return h_atomData;
  }
}

Real** GPUCopy::rollBackCoordinatesPtr() {
  if (parallel) {
    return d_rollBackCoordinates;
  } else {
    return h_rollBackCoordinates;
  }
}

Real** GPUCopy::atomCoordinatesPtr() {
  if (parallel) {
    return d_atomCoordinates;
  } else {
    return h_atomCoordinates;
  }
}

int* GPUCopy::primaryIndexesPtr() {
  if (parallel) {
    return d_primaryIndexes;
  } else {
    return h_primaryIndexes;
  }
}

int** GPUCopy::moleculeDataPtr() {
  if (parallel) {
    return d_moleculeData;
  } else {
    return h_moleculeData;
  }
}

Real* GPUCopy::sizePtr() {
  if (parallel) {
    return d_size;
  } else {
    return h_size;
  }
}

void GPUCopy::copyIn(SimBox * sb) {
  d_moleculeData = (int**)acc_malloc(MOL_DATA_SIZE * sizeof(int *));
  assert(d_moleculeData != NULL);
  for (int row = 0; row < MOL_DATA_SIZE; row++) {
    int *h_moleculeData_row = sb->moleculeData[row];
    int *d_moleculeData_row = (int *)acc_copyin(h_moleculeData_row, sb->numMolecules * sizeof(int));
    assert(d_moleculeData_row != NULL);
    #pragma acc parallel deviceptr(d_moleculeData)
    d_moleculeData[row] = d_moleculeData_row;
  }
  h_moleculeData = sb->moleculeData;

  d_atomData = (Real**)acc_malloc(ATOM_DATA_SIZE * sizeof(Real *));
  assert(d_atomData != NULL);
  for (int row = 0; row < ATOM_DATA_SIZE; row++) {
    Real *h_atomData_row = sb->atomData[row];
    Real *d_atomData_row = (Real *)acc_copyin(h_atomData_row, sb->numAtoms * sizeof(Real));
    assert(d_atomData_row != NULL);
    #pragma acc parallel deviceptr(d_atomData)
    d_atomData[row] = d_atomData_row;
  }
  h_atomData = sb->atomData;

  d_atomCoordinates = (Real**)acc_malloc(NUM_DIMENSIONS * sizeof(Real *));
  assert(d_atomCoordinates != NULL);
  for (int row = 0; row < NUM_DIMENSIONS; row++) {
    Real *h_atomCoordinates_row = sb->atomCoordinates[row];
    Real *d_atomCoordinates_row = (Real *)acc_copyin(h_atomCoordinates_row, sb->numAtoms * sizeof(Real));
    assert(d_atomCoordinates_row != NULL);
    #pragma acc parallel deviceptr(d_atomCoordinates)
    d_atomCoordinates[row] = d_atomCoordinates_row;
  }
  h_atomCoordinates = sb->atomCoordinates;

  d_rollBackCoordinates = (Real**)acc_malloc(NUM_DIMENSIONS * sizeof(Real *));
  assert(d_rollBackCoordinates != NULL);
  for (int row = 0; row < NUM_DIMENSIONS; row++) {
    Real *h_rollBackCoordinates_row = sb->rollBackCoordinates[row];
    Real *d_rollBackCoordinates_row = (Real *)acc_copyin(h_rollBackCoordinates_row, sb->largestMol * sizeof(Real));
    assert(d_rollBackCoordinates_row != NULL);
    #pragma acc parallel deviceptr(d_rollBackCoordinates)
    d_rollBackCoordinates[row] = d_rollBackCoordinates_row;
  }
  h_rollBackCoordinates = sb->rollBackCoordinates;

  d_primaryIndexes = (int *)acc_copyin(sb->primaryIndexes, sb->numPIdxes * sizeof(int));
  h_primaryIndexes = sb->primaryIndexes;

  d_size = (Real *)acc_copyin(sb->size, NUM_DIMENSIONS * sizeof(Real));
  h_size = sb->size;
}

void GPUCopy::copyOut(SimBox* sb) {
  for (int row = 0; row < MOL_DATA_SIZE; row++) {
    int *h_moleculeData_row = h_moleculeData[row];
    acc_copyout(h_moleculeData_row, sb->numMolecules * sizeof(int));
  }

  for (int row = 0; row < ATOM_DATA_SIZE; row++) {
    Real *h_atomData_row = h_atomData[row];
    acc_copyout(h_atomData_row, sb->numAtoms * sizeof(Real));
  }

  for (int row = 0; row < NUM_DIMENSIONS; row++) {
    Real *h_atomCoordinates_row = h_atomCoordinates[row];
    acc_copyout(h_atomCoordinates_row, sb->numAtoms * sizeof(Real));
  }

  for (int row = 0; row < NUM_DIMENSIONS; row++) {
    Real *h_rollBackCoordinates_row = h_rollBackCoordinates[row];
    acc_copyout(h_rollBackCoordinates_row, sb->largestMol * sizeof(Real));
  }

  acc_copyout(h_primaryIndexes, sb->numPIdxes);
  acc_copyout(h_size, NUM_DIMENSIONS);
}

void GPUCopy::setParallel(bool in) {
  parallel = in;
}
