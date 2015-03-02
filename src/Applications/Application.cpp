/// @file Application.cpp
///
/// Contains the main entry point for the program and the method to start
/// running the simulation.
///
/// @author Tavis Maclellan
/// @date Created 2/23/2014
/// @date Updated 2/27/2014

#include <stdlib.h>
#include <stdio.h>

#include "Application.h"
#include "CommandParsing.h"
#include "Metropolis/Simulation.h"
#include "Metropolis/SimulationArgs.h"
#include "Metropolis/Utilities/DeviceQuery.h"
#include <iostream>
#include <fstream>


int metrosim::run(int argc, char** argv)
{
	SimulationArgs args = SimulationArgs();
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

	std::streambuf* cout_sbuf;
	if (!args.verboseOutput) {
		std::cout << "Silent run, integration test started..." << endl; // save original sbuf
		std::streambuf* cout_sbuf = std::cout.rdbuf();
		std::ofstream fout("/dev/null");
		std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
	}
	

	Simulation sim = Simulation(args);
	sim.run();
	
	if (!args.verboseOutput) {
		 std::cout.rdbuf(cout_sbuf); // restore the original stream buffer
	}
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
