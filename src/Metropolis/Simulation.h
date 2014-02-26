/*
	New Simulation to replace linearSim and parallelSim

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include "Metropolis/SimulationArgs.h"
#include "Box.h"
#include "Metropolis/SerialSim/SerialBox.h"
// #include "Metropolis/SerialSim/SerialCalcs.h"
// #include "Metropolis/ParallelSim/ParallelBox.cuh"
// #include "Metropolis/ParallelSim/ParallelCalcs.h"

const double kBoltz = 0.00198717;

class Simulation
{
	private:
		Box *box;
	public:
		Simulation(SimulationArgs args);
		~Simulation();
		void run();
};

#endif