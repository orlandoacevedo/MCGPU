#ifdef _OPENACC
#include <openacc.h>
#endif

#include "GPUCopy.h"

bool parallel = false;

Real** h_atomData = NULL;
Real** d_atomData = NULL;

Real** h_rollBackCoordinates = NULL;
Real** d_rollBackCoordinates = NULL;

Real** h_atomCoordinates = NULL;
Real** d_atomCoordinates = NULL;

int* h_primaryIndexes = NULL;
int* d_primaryIndexes = NULL;

int** h_moleculeData = NULL;
int** d_moleculeData = NULL;

Real* h_size = NULL;
Real* d_size = NULL;

// Angle values

Real** h_angleData = NULL;
Real** d_angleData = NULL;

Real* h_angleSizes = NULL;
Real* d_angleSizes = NULL;

Real* h_rollBackAngleSizes = NULL;
Real* d_rollBackAngleSizes = NULL;

// Bond values

Real** h_bondData = NULL;
Real** d_bondData = NULL;

Real* h_bondLengths = NULL;
Real* d_bondLengths = NULL;

Real* h_rollBackBondLengths = NULL;
Real* d_rollBackBondLengths = NULL;

void GPUCopy::setParallel(bool in) { parallel = in; }

int GPUCopy::onGpu() { return parallel; }

Real** GPUCopy::atomDataPtr() { return parallel ? d_atomData : h_atomData; }

Real** GPUCopy::rollBackCoordinatesPtr() {
  return parallel ? d_rollBackCoordinates : h_rollBackCoordinates;
}

Real** GPUCopy::atomCoordinatesPtr() {
  return parallel ? d_atomCoordinates : h_atomCoordinates;
}

int* GPUCopy::primaryIndexesPtr() {
  return parallel ? d_primaryIndexes : h_primaryIndexes;
}

int** GPUCopy::moleculeDataPtr() {
  return parallel ? d_moleculeData : h_moleculeData;
}

Real** GPUCopy::bondDataPtr() {
  return parallel ? d_bondData : h_bondData;
}

Real* GPUCopy::bondLengthsPtr() {
  return parallel ? d_bondLengths : h_bondLengths;
}

Real* GPUCopy::rollBackBondsPtr() {
  return parallel ? d_rollBackBondLengths : h_rollBackBondLengths;
}

Real** GPUCopy::angleDataPtr() {
  return parallel ? d_angleData : h_angleData;
}

Real* GPUCopy::angleSizesPtr() {
  return parallel ? d_angleSizes : h_angleSizes;
}

Real* GPUCopy::rollBackAnglesPtr() {
  return parallel ? d_rollBackAngleSizes : h_rollBackAngleSizes;
}

Real* GPUCopy::sizePtr() { return parallel ? d_size : h_size; }

void GPUCopy::copyIn(SimBox *sb) {
  h_moleculeData = sb->moleculeData;
  h_atomData = sb->atomData;
  h_atomCoordinates = sb->atomCoordinates;
  h_rollBackCoordinates = sb->rollBackCoordinates;
  h_size = sb->size;
  h_primaryIndexes = sb->primaryIndexes;
  h_bondData = sb->bondData;
  h_bondLengths = sb->bondLengths;
  h_rollBackBondLengths = sb->rollBackBondLengths;
  h_angleData = sb->angleData;
  h_angleSizes = sb->angleSizes;
  h_rollBackAngleSizes = sb->rollBackAngleSizes;
  if (!parallel) { return; }

#ifdef _OPENACC
  d_moleculeData = (int**)acc_malloc(MOL_DATA_SIZE * sizeof(int*));
  assert(d_moleculeData != NULL);
  for (int row = 0; row < MOL_DATA_SIZE; row++) {
    int *h_moleculeData_row = sb->moleculeData[row];
    int *d_moleculeData_row = (int*)acc_copyin(h_moleculeData_row,
                                               sb->numMolecules * sizeof(int));
    assert(d_moleculeData_row != NULL);
    #pragma acc parallel deviceptr(d_moleculeData)
    d_moleculeData[row] = d_moleculeData_row;
  }

  d_atomData = (Real**)acc_malloc(ATOM_DATA_SIZE * sizeof(Real*));
  assert(d_atomData != NULL);
  for (int row = 0; row < ATOM_DATA_SIZE; row++) {
    Real* h_atomData_row = sb->atomData[row];
    Real* d_atomData_row = (Real*)acc_copyin(h_atomData_row, sb->numAtoms * sizeof(Real));
    assert(d_atomData_row != NULL);
    #pragma acc parallel deviceptr(d_atomData)
    d_atomData[row] = d_atomData_row;
  }

  d_atomCoordinates = (Real**)acc_malloc(NUM_DIMENSIONS * sizeof(Real *));
  assert(d_atomCoordinates != NULL);
  for (int row = 0; row < NUM_DIMENSIONS; row++) {
    Real* h_atomCoordinates_row = sb->atomCoordinates[row];
    Real* d_atomCoordinates_row = (Real *)acc_copyin(h_atomCoordinates_row, sb->numAtoms * sizeof(Real));
    assert(d_atomCoordinates_row != NULL);
    #pragma acc parallel deviceptr(d_atomCoordinates)
    d_atomCoordinates[row] = d_atomCoordinates_row;
  }

  d_rollBackCoordinates = (Real**)acc_malloc(NUM_DIMENSIONS * sizeof(Real *));
  assert(d_rollBackCoordinates != NULL);
  for (int row = 0; row < NUM_DIMENSIONS; row++) {
    Real* h_rollBackCoordinates_row = sb->rollBackCoordinates[row];
    Real* d_rollBackCoordinates_row = (Real *)acc_copyin(h_rollBackCoordinates_row, sb->largestMol * sizeof(Real));
    assert(d_rollBackCoordinates_row != NULL);
    #pragma acc parallel deviceptr(d_rollBackCoordinates)
    d_rollBackCoordinates[row] = d_rollBackCoordinates_row;
  }

  d_primaryIndexes = (int*)acc_copyin(sb->primaryIndexes, sb->numPIdxes * sizeof(int));

  d_size = (Real*)acc_copyin(sb->size, NUM_DIMENSIONS * sizeof(Real));

  d_angleData = (Real**)acc_malloc(ANGLE_DATA_SIZE * sizeof(Real*));
  assert(d_angleData != NULL);
  for (int row = 0; row < ANGLE_DATA_SIZE; row++) {
    Real *h_angleDataRow = sb->angleData[row];
    Real *d_angleDataRow = (Real*)acc_copyin(h_angleDataRow, sb->numAngles * sizeof(Real));
    assert(d_angleDataRow != NULL);
    #pragma acc parallel deviceptr(d_angleData)
    d_angleData[row] = d_angleDataRow;
  }

  d_angleSizes = (Real*)acc_copyin(sb->angleSizes,
                                   sb->numAngles * sizeof(Real));
  d_rollBackAngleSizes = (Real*)acc_copyin(sb->rollBackAngleSizes,
                                           sb->numAngles * sizeof(Real));

  d_bondData = (Real**)acc_malloc(BOND_DATA_SIZE * sizeof(Real*));
  assert(d_bondData != NULL);
  for (int row = 0; row < BOND_DATA_SIZE; row++) {
    Real *h_bondDataRow = sb->bondData[row];
    Real *d_bondDataRow = (Real*)acc_copyin(h_bondDataRow, sb->numBonds * sizeof(Real));
    assert(d_bondDataRow != NULL);
    #pragma acc parallel deviceptr(d_bondData)
    d_bondData[row] = d_bondDataRow;
  }

  d_bondLengths = (Real*)acc_copyin(sb->bondLengths, 
                                    sb->numBonds * sizeof(Real));
  d_rollBackBondLengths = (Real*)acc_copyin(sb->rollBackBondLengths, 
                                            sb->numBonds * sizeof(Real));
#endif
}

void GPUCopy::copyOut(SimBox* sb) {
  if (!parallel) return;

#ifdef _OPENACC
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

  acc_copyout(h_angleSizes, sb->numAngles);
  acc_copyout(h_rollBackAngleSizes, sb->numAngles);
  acc_copyout(h_bondLengths, sb->numBonds);
  acc_copyout(h_rollBackBondLengths, sb->numBonds);
#endif
}

