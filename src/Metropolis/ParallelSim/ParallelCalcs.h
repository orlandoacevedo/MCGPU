/*
	Contains calculations for ParallelBox

	Author: Nathan Coleman
*/

#ifndef PARALLELCALCS_H
#define PARALLELCALCS_H

#include "Metropolis/Box.h"

double calcMolecularEnergyContribution(Box *box, int molIdx, int startIdx = 0);
double ParallelSim::calcSystemEnergy();

#endif