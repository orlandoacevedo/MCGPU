/**
 * @file Application.cpp
 *
 *  Contains the main entry point for the program and the method to start
 *  running the simulation.
 *
 * @author Tavis Maclellan
 * @date Created 2/23/2014
 * @date Updated 1/7/2016
 */

#include "Application.h"
#include <iostream>
#include <fstream>
#include <openacc.h>


int metrosim::run(int argc, char** argv) {
	SimulationArgs args = SimulationArgs();

	if (!getCommands(argc, argv, &args)) {
		exit(EXIT_FAILURE);
	}

	int numDevices = acc_get_num_devices(acc_device_nvidia);

	// Ensure args.simulationMode is either Serial or Parallel
	if (args.simulationMode == SimulationMode::Default) {
		if (numDevices < 1) {
			fprintf(stdout, "No CUDA devices found; defaulting to CPU execution.\n");
			args.simulationMode = SimulationMode::Serial;
		} else {
			fprintf(stdout, "%d CUDA device(s) found; running on GPU.\n", numDevices);
			args.simulationMode = SimulationMode::Parallel;
		}
	}

	bool parallelMode = (args.simulationMode == SimulationMode::Parallel);

	if (parallelMode) {
		if (numDevices == 0) {
			fprintf(stdout, "ERROR: Cannot find suitable GPU!\n");
			exit(EXIT_FAILURE);
		}
		acc_init(acc_device_nvidia);
		acc_set_device_type(acc_device_nvidia);
	}

	if (parallelMode) {
		fprintf(stdout, "Beginning simulation using GPU...\n");
	} else {
		fprintf(stdout, "Beginning simulation using CPU...\n");
	}

	Simulation sim = Simulation(args);
	sim.run();

	fprintf(stdout, "Finishing simulation...\n\n");

	if (parallelMode) {
		acc_shutdown(acc_device_nvidia);
	}

	exit(EXIT_SUCCESS);
}
