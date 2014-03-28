/*
	New Simulation to replace linearSim and parallelSim

	Author: Nathan Coleman
	Last Changed: February 21, 2014
	
	-> February 26, by Albert Wallace
*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include "SimulationArgs.h"
#include "Box.h"

const double kBoltz = 0.00198717;

class Simulation
{
	public:
		Simulation(SimulationArgs simArgs);
		~Simulation();
		void run();
		
	private:
		Box *box;
		SimulationArgs args;
		int simSteps;
};

#endif