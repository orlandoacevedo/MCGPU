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

#ifdef _OPENACC
#include <openacc.h>
#endif


int metrosim::run(int argc, char** argv) {
	SimulationArgs args = SimulationArgs();

	if (!getCommands(argc, argv, &args)) {
		exit(EXIT_FAILURE);
	}

#ifdef _OPENACC
	int numDevices = acc_get_num_devices(acc_device_nvidia);

	// Ensure args.simulationMode is either Serial or Parallel
	if (args.simulationMode == SimulationMode::Default) {
		if (numDevices < 1) {
			fprintf(stdout, "No GPU devices found; defaulting to CPU "
                            "execution.\n");
			args.simulationMode = SimulationMode::Serial;
		} else {
			fprintf(stdout, "%d GPU device(s) found; running on GPU.\n",
                    numDevices);
      // FIXME (blm)
      fprintf(stdout, "WARNING: Not all features support GPU offloading at "
                      "this time. Results may be inaccurate.\n");
			args.simulationMode = SimulationMode::Parallel;
		}
	}

    if (args.simulationMode == SimulationMode::Parallel) {
        if (numDevices == 0) {
            fprintf(stdout, "ERROR: Cannot find suitable GPU!\n");
            exit(EXIT_FAILURE);
        }

        // Implicitly sets the device
        acc_init(acc_device_nvidia);
    }

#else
    // Without OpenACC, only serial calculations are supported
    if (args.simulationMode == SimulationMode::Parallel) {
        fprintf(stdout, "Must compile with OpenACC capability to run in "
                "parallel mode.\n");
        exit(EXIT_FAILURE);
    } else {
        args.simulationMode = SimulationMode::Serial;
		fprintf(stdout, "Beginning simulation using CPU...\n");
    }
#endif

	Simulation sim = Simulation(args);
	sim.run();

	fprintf(stdout, "Finishing simulation...\n\n");

// Shutdown the device if we were using the GPU
#ifdef _OPENACC
	if (args.simulationMode == SimulationMode::Parallel) {
		acc_shutdown(acc_device_nvidia);
	}
#endif

	exit(EXIT_SUCCESS);
}
