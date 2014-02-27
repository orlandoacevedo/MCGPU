/*
	New Simulation to replace linearSim and parallelSim

	Author: Nathan Coleman
	Last Changed: February 21, 2014
	
	-> February 26, by Albert Wallace
*/
#ifndef SIMULATION_CPP
#define SIMULATION_CPP


#include "Simulation.h"
//#include "Metropolis/SimulationArgs.h" // AlbertIncludes
#include "Box.h"
#include "SimulationArgs.h" // AlbertIncludes
#include "SerialSim/SerialBox.cpp" // AlbertIncludes
//#include "Metropolis/SerialSim/SerialBox.h" //AlbertIncludes
//#include "Metropolis/SerialSim/SerialCalcs.h"
//#include "Metropolis/ParallelSim/ParallelBox.cuh"
//#include "Metropolis/ParallelSim/ParallelCalcs.h"


//Constructor & Destructor
Simulation::Simulation(SimulationArgs args)
{
	IOUtilities configScan(args.configPath);
	configScan.readInConfig();
	box = (Box*) (new SerialBox(configScan));
}

Simulation::~Simulation()
{
	if( box == NULL )
	{
		delete box;
		box = NULL;
	}
}

//Utility
void Simulation::run()
{
	Environment *boxEnviro = box->getEnvironment();
	printf("X: %f\n",boxEnviro->x);
	printf("Y: %f\n",boxEnviro->y);
	printf("Z: %f\n",boxEnviro->z);
}

#endif