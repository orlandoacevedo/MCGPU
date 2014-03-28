/*
	Contains calculations for SerialBox

	Author: Nathan Coleman
*/

#ifndef SERIALCALCS_H
#define SERIALCALCS_H

#include "Metropolis/DataTypes.h"
#include "Metropolis/Utilities/StructLibrary.h"

namespace SerialCalcs
{

Real calcBlending(Real d1, Real d2);
Real calcCharge(Real charge1, Real charge2, Real r);
Real calcInterMolecularEnergy(Molecule *molecules, int mol1, int mol2, Environment *environment);
Real calc_lj(Atom atom1, Atom atom2, Real r2);
Real calcMolecularEnergyContribution(Molecule *molecules, Environment *environment, int currentMol, int startIdx = 0);
Real calcSystemEnergy(Molecule *molecules, Environment *environment);
Real getFValue(int atom1, int atom2, int **table1);
Real makePeriodic(Real x, Real box);
}

#endif