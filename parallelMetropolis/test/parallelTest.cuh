#ifndef PARALLELTEST_H
#define PARALLELTEST_H

#include <float.h>
#include <assert.h>
#include <cuda.h>
#include "baseTests.h"
#include <vector>
#include "../src/metroCudaUtil.cuh"
#include "../../Utilities/src/geometricUtil.h"
#include "copyTests.cuh"
#include <sys/time.h>

/**
Tests the getXFromIndex function. The global Kernel is used to setup the tests.
*/
void setupGetXFromIndex();

/**
Tests the getYFromIndex function. The global Kernel is used to setup the tests.
*/
void setupGetYFromIndex();

/**
Tests the makePeriodic function. The global Kernel is used to setup the tests.
*/
void setupMakePeriodic();

/**
Tests the wrapBox Function. The global Kernel is used to setup the tests.
*/
void setupWrapBox();

/**
Tests the keepMoleculeInBox() function.
*/
void testKeepMoleculeInBox();

/**
Tests the calc_lj Function. The global Kernel is used to setup the tests.
*/
void setupCalc_lj();

/**
  Tests the function that generates points against the known bounds for the
  generated points
*/
void testGeneratePoints();

/**
  Tests the energy calculated using parallel method vs the energy caclulated
  by the linear method.
*/
void testCalcEnergy();

/**
  Tests the calc energy function that accepts an array of molecules instead
  of an array of atoms.
*/
void testCalcEnergyWithMolecules();

/**
  Tests getMoleculeFromID by calling testGetMoleculeFromID()
*/
void testGetMoleculeFromIDWrapper();

/**
  Tests calcBlending by calling testCalcBlending()
*/
void testCalcBlendingWrapper();

/**
  Tests getFValue by calling testGetFValue()
*/
void testGetFValueWrapper();

/**
  Tests calcCharge() by calling testCalcCharge()
*/
void testCalcChargeWrapper();

/**
  Tests getDistance() device function via testGetDistance() global function.
*/
void testGetDistanceWrapper();

/**
  @param pivot - the atom about which the rotation will take place.
  @param toRotate - the atom to rotate
  @param rotation - the amount to rotate about each axis
  @return - the position of the atom after the max rotation.
*/
Atom findMaxRotation(Atom pivot, Atom toRotate, double rotation);

/**
  Tests that the function that rotates about a certain atom in a molecule is correct
*/
void testRotateMolecule();

/**
Wrapper function that calls a global function to test the calc charge device function.
*/
void testCalcChargeWrapper();

#endif //PARALLELTEST_H
