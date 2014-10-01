/*
	Driver for the simulation. Takes in a SimulationArgs object and creates the
	the necessary Box type, state file output path, etc.

	Author: Nathan Coleman
	Created: February 21, 2014
	
	-> February 26, by Albert Wallace
	-> March 28, by Joshua Mosby
	-> April 21, by Nathan Coleman
*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include "SimulationArgs.h"
#include "Box.h"

#define OUT_INTERVAL 100

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
		long simSteps;
		long stepStart;

		int writePDB(Environment sourceEnvironment, Molecule * sourceMoleculeCollection);
		void saveState(const std::string& simName, int simStep);
		const std::string currentDateTime();
};

#endif
