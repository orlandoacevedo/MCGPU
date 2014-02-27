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


int metrosim::run(int argc, char** argv)
{
	SimulationArgs args = SimulationArgs();
	if (!getCommands(argc, argv, &args))
	{
		exit(EXIT_FAILURE);
	}

	fprintf(stdout, "Running simulation...\n\n");

	Simulation sim = Simulation(args);
	sim.run();

	exit(EXIT_SUCCESS);
}

int main(int argc, char** argv)
{
	metrosim::run(argc, argv);
	
}