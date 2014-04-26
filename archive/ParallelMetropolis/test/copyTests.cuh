#ifndef COPYTEST_H
#define COPYTEST_H

#include <assert.h>
#include "../src/metroCudaUtil.cuh"

/**
  Tests the moleculeDeepCopyToDevice() function
  Tests the moleculeDeepCopyToHost() function.
*/
void testCopyMolecules();

#endif // COPYTEST_H
