/// @file Application.cpp
///
/// Contains the main entry point for the program and the method to start
/// running the simulation.
///
/// @author Tavis Maclellan
/// @date Created 2/23/2014
/// @date Updated 2/26/2014

#include <stdlib.h>
#include <stdio.h>
#include "Application.h"
#include "Metropolis/SimulationArgs.h"
#include "CommandParsing.h"
#include "Metropolis/Simulation.h"

int metrosim::run(int argc, char** argv)
{
	SimulationArgs args;
	if (!getCommands(argc, argv, &args))
	{
		exit(EXIT_FAILURE);
	}

	fprintf(stdout, "Running simulation...\n\n");

	Simulation* sim = new Simulation(args);
	sim->run();
	delete sim;

	exit(EXIT_SUCCESS);
}

int main(int argc, char** argv)
{
	metrosim::run(argc, argv);
}