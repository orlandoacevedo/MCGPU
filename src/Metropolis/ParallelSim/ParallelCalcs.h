/*
	Declares ParallellCalcs methods accessed from C++ compiled code (no Cuda).

	Created: February 21, 2014
	
	-> February 26, by Albert Wallace
	-> March 28, by Joshua Mosby
	-> April 21, by Nathan Coleman
*/

#ifndef PARALLELCALCS_H
#define PARALLELCALCS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include "Metropolis/Box.h"
#include "Metropolis/DataTypes.h"
#include "Metropolis/SimulationArgs.h"

#define NMAX 100000  /* Maximum number of atoms which can be simulated */
#define NCLMAX 10000 /* Maximum number of linked-list cells */
#define EMPTY -1

namespace ParallelCalcs
{
	/// Factory method for creating a Box from a configuration file.
	/// @param configpath The path to the configuration file.
	/// @param steps The number of steps desired in the simulation,
	/// @return Returns a pointer to the filled-in Box.
	/// @note This functionality should ideally reside in ParallelBox,
	///   but it was placed here due to time constraints.
	///   TODO for future group.
	Box* createBox(std::string inputPath, InputFileType inputType, long* startStep, long* steps);
	
	/// Calculates the system energy using consecutive calls to
	///   calcMolecularEnergyContribution.
	/// @param box A ParallelBox cast as a Box, passed from Simulation.
	/// @return Returns total system energy.
	Real calcSystemEnergy(Box *box);
	
	/// Parallel version of the Neighbor List function
	/// Calculates the system energy using consecutive calls to
	///   calcMolecularEnergyContribution.
	/// @param box A ParallelBox cast as a Box, passed from Simulation.
	/// @return Returns total system energy.
	Real calcSystemEnergy_NLC(Box *box);
	

	/// Calculates the inter-molecular energy contribution of a given molecule,
	///   without intramolecular energy, using a batch method.
	/// @param box A ParallelBox cast as a Box, passed from Simulation.
	/// @param molIdx the index of the current changed molecule.
	/// @param startIdx The optional starting index for other molecules.
	///   Used for system energy calculation.
	/// @return Returns total molecular energy contribution, without
	///   intramolecular energy.
	Real calcMolecularEnergyContribution(Box *box, int molIdx, int startIdx = 0);
    Real calcMolecularEnergyContribution_NLC(Box *box, int currentMol);
}

#endif
