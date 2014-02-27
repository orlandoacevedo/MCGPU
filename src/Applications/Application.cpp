/// @file Application.cpp
///
/// Contains the main entry point for the program and the method to start
/// running the simulation.
///
/// @author Tavis Maclellan
/// @date Created 2/23/2014
/// @date Updated 2/26/2014

#ifndef METROSIM_APPLICATION_CPP
#define METROSIM_APPLICATION_CPP

#include <stdlib.h>
#include <stdio.h>


#include "Application.h"
//#include "Metropolis/SimulationArgs.h" // AlbertIncludes
#include "../Metropolis/SimulationArgs.h" // AlbertIncludes
#include "CommandParsing.cpp"
//#include "Metropolis/Simulation.h" // AlbertIncludes
#include "../Metropolis/Simulation.cpp" // AlbertIncludes



int metrosim::run(int argc, char** argv)
{
	SimulationArgs args = SimulationArgs();
	if (!getCommands(argc, argv, &args))
	{
		exit(EXIT_FAILURE);
	}

	fprintf(stdout, "Running simulation...\n\n");

	Simulation::Simulation* sim = new Simulation::Simulation(args);
	sim->run();
	delete sim;

	exit(EXIT_SUCCESS);
}

int main(int argc, char** argv)
{
	metrosim::run(argc, argv);
	
}

#endif