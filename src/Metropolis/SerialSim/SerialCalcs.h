/*
	Contains the methods required to calculate energies serially.

	Created: February 21, 2014

	-> February 26, by Albert Wallace
	-> March 28, by Joshua Mosby
	-> April 21, by Nathan Coleman
*/

#ifndef SERIALCALCS_H
#define SERIALCALCS_H

#include <string>
#include "Metropolis/Box.h"
#include "SerialBox.h"
#include "NeighborList.h"
#include "Metropolis/DataTypes.h"
#include "Metropolis/SimulationArgs.h"
#include "Metropolis/Utilities/StructLibrary.h"
#include "Metropolis/Utilities/FileUtilities.h"

namespace SerialCalcs {
	
	/**
	 * Factory method for creating a Box from a configuration file.
	 * 	@param configpath The path to the configuration file.
	 * 	@param steps The number of steps desired in the simulation,
	 * 	@return Returns a pointer to the filled-in Box.
	 * 	@note This functionality should ideally reside in SerialBox,
	 *   		but it was placed here due to time constraints.
	 *   		TODO for future group.
   */
	Box* createBox(SimulationArgs& simArgs, long* startStep, long* steps,
                   SBScanner* sbScanner);

	/**
	 * Calculates the long-range correction energy value for molecules outside the cutoff.
	 * 	@note *** This needs to be updated if the volume ever changes or to support more than 2 solvents
	 * 	@param box A Box containing the molecule data.
	 * 	@return Returns r2, the distance between the two atoms squared.
   */
	Real calcEnergy_LRC(Box* box);
}

#endif
