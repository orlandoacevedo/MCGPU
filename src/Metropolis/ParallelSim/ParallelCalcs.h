/*
	Contains calculations for ParallelBox

	Author: Nathan Coleman
*/

#ifndef PARALLELCALCS_H
#define PARALLELCALCS_H

#include "Metropolis/Box.h"
#include "Metropolis/DataTypes.h"

Real calcMolecularEnergyContribution(Box *box, int molIdx, int startIdx = 0);
Real calcSystemEnergy();

#endif