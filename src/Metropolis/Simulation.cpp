/*
	New Simulation to replace linearSim and parallelSim

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

#include "Simulation.h"
#include "Metropolis/SimulationArgs.h"
#include "Box.h"
#include "Metropolis/SerialSim/SerialBox.h"
#include "Metropolis/SerialSim/SerialCalcs.h"
#include "Metropolis/ParallelSim/ParallelBox.cuh"
#include "Metropolis/ParallelSim/ParallelCalcs.h"

using namespace std;

//Constructor & Destructor
Simulation::Simulation(SimulationArguments args){}
Simulation::~Simulation(){}

//Utility
void Simulation::run(){}