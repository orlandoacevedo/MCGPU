/// @file Application.h
///
/// Contains the declaration for the method to run the application.
///
/// @author Tavis Maclellan
/// @date Created 2/23/2014
/// @date Updated 2/26/2014

#ifndef METROSIM_APPLICATION_H
#define METROSIM_APPLICATION_H

#include <stdlib.h>
#include <stdio.h>

#include "CommandParsing.h"
#include "Metropolis/Simulation.h"
#include "Metropolis/SimulationArgs.h"


namespace metrosim
{

	/// Run the application.
	///
	/// @param argc The count of command line arguments passed to the program.
	/// @param argv The array of command line argument strings passed to the
	///   application.
	/// @note The command line argument list will always have at least one entry
	///   because the name of the program as it was entered via the command
	///   line will always be the first entry in the arguments list (i.e. 
	///   @e argv[0] will be the invoked program name).
	int run(int argc, char** argv);
}

#endif