#ifndef GPUCOPY_H
#define GPUCOPY_H

#include "SimBox.h"

namespace GPUCopy {
  void copyIn(SimBox* sb);
  void copyOut(SimBox* sb);
  void setParallel(bool in);
  Real** atomDataPtr();
  Real** rollBackCoordinatesPtr();
  Real** atomCoordinatesPtr();
  int* primaryIndexesPtr();
  int** moleculeDataPtr();
  Real* sizePtr();
  int onGpu();
}


#endif
