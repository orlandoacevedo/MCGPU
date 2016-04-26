/// @file Application.cpp
///
/// Contains the main entry point for the program and the method to start
/// running the simulation.
///
/// @author Tavis Maclellan
/// @date Created 2/23/2014
/// @date Updated 2/27/2014

// #include <stdlib.h>
// #include <stdio.h>

#include "Application.h"
// #include "CommandParsing.h"
// #include "Metropolis/Simulation.h"
// #include "Metropolis/SimulationArgs.h"
#include "Metropolis/Utilities/DeviceQuery.h"
#include <iostream>
#include <fstream>


int metrosim::run(int argc, char** argv)
{
	SimulationArgs args = SimulationArgs();
	//args.useNeighborList = false;
	if (!getCommands(argc, argv, &args))
	{
		exit(EXIT_FAILURE);
	}

	DeviceContext context = DeviceContext();
	if (args.simulationMode == SimulationMode::Parallel)
	{
		if (!openDeviceContext(&context, MIN_MAJOR_VER, MIN_MINOR_VER, args.deviceIndex))
		{
			exit(EXIT_FAILURE);
		}
	}

	if (args.simulationMode == SimulationMode::Parallel)
	{
		fprintf(stdout, "Beginning simulation using GPU...\n");
	}
	else
	{
		fprintf(stdout, "Beginning simulation using CPU...\n");
	}

	Simulation sim = Simulation(args);
	sim.run();

	fprintf(stdout, "Finishing simulation...\n\n");

	if (args.simulationMode == SimulationMode::Parallel)
	{
		if (!closeDeviceContext(&context))
		{
			exit(EXIT_FAILURE);
		}
	}

	exit(EXIT_SUCCESS);
}
