#ifndef ZMATRIXTEST_H
#define ZMATRIXTEST_H

#include "../src/Opls_Scan.h"
#include "../src/Zmatrix_Scan.h"
#include "../src/State_Scan.h"
#include "stateTest.h"
#include "configurationTest.h"
#include "geometricTest.h"
#include <assert.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

/**
  creates a copy of the molecule represented in the mesh.z z-matrix file 
  @param scanner - an Opls_Scan scanner object used to asaign the charge, sigma, and epsilon
  @returns - a Molecule object with the data from mesh.z matrix file
*/
Molecule createMeshZMolecules(Opls_Scan *scanner);

/**
  sets the x y z positions to a molecule based on the relationships defined in mesh.z z-matrix file
  @param meshMolec - the molecule to assign the positions to 
*/
void addPositionsToMeshZ(Molecule *meshMolec);

/**
  creates a copy of the molecule represented in the t3pdim.z z-matrix file 
  @param scanner - an Opls_Scan scanner object used to asaign the charge, sigma, and epsilon
  @returns - a vector with the 2 bonded molecules
*/
vector<Molecule> createT3pdimMolecules(Opls_Scan *scanner);

/**
  sets the x y z positions to a molecule based on the relationships defined in t3pdim.z z-matrix file
  @param t3pdimMolec - the vector of bonded molecules to assign the positions to. 
*/
void addPositionsToT3pdim(vector<Molecule> t3pdimMolec);

/**
  Compares all the values of two Molecules using asserts.
  @param molec1 - the first molecule to compare 
  @param molec2 - the second molecule to compare
*/
void compareTestMolecules(Molecule molec1, Molecule molec2);

/**
  test the Zmatrix Scanner with the mesh.z file
  @param opls - Opls_Scan object
*/
void testZmatrixScanner(Opls_Scan *opls);

/**
  test the Zmatrix Scanner with the test.z file to check for the the id assignments
  when resusing the buildMolecule() function
  @param opls - Opls_Scan object
*/
void testZmatrixScanner_multipleSingle(Opls_Scan *opls);

/**
  test the Zmatrix Scanner with the t3pdim.z file 
  @param opls - Opls_Scan object
*/
void testZmatrixScanner_BondedMolecules(Opls_Scan *opls);

/**
  test all the Zmatrix 
  @param opls - Opls_Scan object
*/
void testZmatrix(Opls_Scan *opls);

#endif //ZMATRIXTEST_H
