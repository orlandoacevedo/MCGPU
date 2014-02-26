/*
	Contains calculations for SerialBox

	Author: Nathan Coleman
*/

#ifndef SERIALCALCS_H
#define SERIALCALCS_H

#include "Metropolis/Utilities/StructLibrary.h"

double calcBlending(double d1, double d2);
double calcCharge(double charge1, double charge2, double r);
double calcInterMolecularEnergy(Molecule *molecules, int mol1, int mol2, Environment *environment);
double calc_lj(Atom atom1, Atom atom2, double r2);
double calcMolecularEnergyContribution(Molecule *molecules, Environment *environment, int currentMol, int startIdx = 0);
double calcSystemEnergy(Molecule *molecules, Environment *environment);
double getFValue(int atom1, int atom2, int **table1);
double makePeriodic(double x, double box);

#endif