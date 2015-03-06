#ifndef SERIALCALCS_H
#define SERIALCALCS_H

#include <string>
#include "Metropolis/Box.h"
#include "SerialBox.h"
#include "Metropolis/DataTypes.h"
#include "Metropolis/SimulationArgs.h"
#include "Metropolis/Utilities/StructLibrary.h"

// Linked-cell neighbor list constants
#define NMAX 100000  /* Maximum number of atoms which can be simulated */
#define NCLMAX 10000 /* Maximum number of linked-list cells */
#define EMPTY -1

namespace SerialCalcs
{
    __global__ void calcEnergy_NLC(Molecule *molecules, Environment *enviro, int *head, int *lscl, Real *totalEnergy);
}

#endif
