/*
	New Simulation to replace linearSim and parallelSim

	Author: Nathan Coleman
	Last Changed: February 21, 2014
	
	-> February 26, by Albert Wallace
*/

#ifndef SIMULATION_H
#define SIMULATION_H

//#include "Metropolis/SimulationArgs.h" // AlbertIncludes
#include "SimulationArgs.h" // AlbertIncludes
#include "Box.h"
#include "SerialSim/SerialBox.cpp" // AlbertIncludes
//#include "Metropolis/SerialSim/SerialBox.h" // AlbertIncludes
// #include "Metropolis/SerialSim/SerialCalcs.h"
// #include "Metropolis/ParallelSim/ParallelBox.cuh"
// #include "Metropolis/ParallelSim/ParallelCalcs.h"

const double kBoltz = 0.00198717;

class Simulation
{
	public:
		Simulation(SimulationArgs args);
		~Simulation();
		void run();
	private:
		Box *box;
};

#endif