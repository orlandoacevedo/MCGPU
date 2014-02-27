/*
	New Simulation to replace linearSim and parallelSim

	Author: Nathan Coleman
	Last Changed: February 21, 2014
	
	-> February 26, by Albert Wallace
*/
#include <stdio.h>

#include "Simulation.h"
#include "SimulationArgs.h"
#include "Box.h"
#include "SerialSim/SerialBox.h"
#include "Utilities/IOUtilities.h"


//Constructor & Destructor
Simulation::Simulation(SimulationArgs args)
{
	IOUtilities configScan = IOUtilities(args.configPath);
	configScan.readInConfig();
	box = new SerialBox(configScan);
}

Simulation::~Simulation()
{
	if(box != NULL)
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