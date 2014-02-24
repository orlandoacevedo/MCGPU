/*
	Contains calculations for SerialBox

	Author: Nathan Coleman
*/

#ifndef SERIALCALCS_H
#define SERIALCALCS_H

#include "Metropolis/Utilities/StructLibrary.h"

void calcContribution(Molecule *mol);

void calcInterMolecularEnergy(Molecule *molecules, int currentMol, int numM, Environment *environment, double *energies, int segmentSize);
void calcInterAtomicEnergy(Molecule *molecules, int currentMol, int otherMol, Environment *environment, double *energies, int segmentSize);
void calcIntraMolecularEnergy(Molecule *molecules, int currentMol, int numE, Environment *environment, double *energies, int segmentSize);
double calcLJ(Atom atom1, Atom atom2, double r2);
double calcCharge(double charge1, double charge2, double r);
double makePeriodic(double x, double box);
double calcBlending(double d1, double d2);
int getXFromIndex(int index);
int getYFromIndex(int x, int index);

#endif