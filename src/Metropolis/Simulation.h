/*
	New Simulation to replace linearSim and parallelSim

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include "Utilities/Opls_Scan.h"
#include "Utilities/Config_Scan.h"
#include "Utilities/metroUtil.h"
#include "Utilities/Zmatrix_Scan.h"
#include "Utilities/State_Scan.h"
#include "Metropolis/SimulationArgs.h"
#include "Metropolis/Box.h"
#include "Metropolis/SerialSim/SerialBox.h"
#include "Metropolis/ParallelSim/ParallelBox.cuh"

#define THREADS_PER_BLOCK 128
#define PI 3.14159265
// Linked-cell neighbor list constants
#define NMAX 100000  /* Maximum number of atoms which can be simulated */
#define NCLMAX 10000 /* Maximum number of linked-list cells */
#define EMPTY -1

const double kBoltz = 0.00198717;

class Simulation
{
	Simulation(SimulationArguments args);
	~Simulation();
	void run();
};

#endif