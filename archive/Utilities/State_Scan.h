/*!\file
  \brief Structures and functions used to read and write state files.
  \author Alexander Luchs, Riley Spahn, Seth Wooten
 
 */

#ifndef STATE_SCAN_H
#define STATE_SCAN_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <string.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "metroUtil.h"

using namespace std;


/**
  @param enviro - the environment state
  @param molecules - array of molecules to be printed out
  @param numOfMolecules - the number of molecules to be written out
  @param fileName - the name of the file to be written
*/
void printState(Environment *enviro, Molecule *molecules, int numOfMolecules, string filename);

/**
  @param filename - the name of the state file
  @return - the environment recorded in the state file.
*/
Environment readInEnvironment(string filename);

/**
  @pararm filename - the name of the state file
  @return - an array of molecules
*/
vector<Molecule> readInMolecules(string filename);

/**
  input line is of the format:
  "id x y z sigma epsilon"
  @param line - the line from the state file that contains atom information.
  @return - atom containing information from the line
*/
Atom getAtomFromLine(string line);

/**
  input line is of the format:
  "x y z numOfAtoms"
  @param line - the line to be parsed
  @return - Environment read from the line
*/
Environment getEnvironmentFromLine(string line);

/**
  expected input line:
  "atom1 atom2 distance [0|1]"
  0 represents not variable
  1 represents a variable bond
  @param line - the line to be read.
  @return - bond representing the information read from the line.
*/
Bond getBondFromLine(string line);

/**
  expected input line:
  "atom1 atom2 value [0|1]"
  0 represents not variable
  1 represents variable angle

  @param line - line containing information about the angle
  @return - returns a bond 
*/
Angle getAngleFromLine(string line);

/**
  expected input line:
  "atom1 atom2 value [0|1]"
  0 represents not variable
  1 represents variable dihedral 

  @param line - line containing information about the dihedral 
  @return - returns the dihedral represented on the line
*/
Dihedral getDihedralFromLine(string line);

/**
  expected input line:
  "atom1 atom2 hopDistance"

  @param line - line containing information about the hop
  @return - a hop represented by the information on the line.
*/
Hop getHopFromLine(string line);

#endif // STATE_SCAN_H
