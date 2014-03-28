/*
	Contains calculations for ParallelBox

	Author: Nathan Coleman
*/

#ifndef PARALLELCALCS_H
#define PARALLELCALCS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Metropolis/Box.h"
#include "Metropolis/DataTypes.h"

Real calcSystemEnergy();
Real calcMolecularEnergyContribution(int molIdx, int startIdx = 0);
int createMolBatch(int curentMol, int startIdx);
Real calcBatchEnergy(int numMols, int molIdx);
Real getEnergyFromDevice();

#endif