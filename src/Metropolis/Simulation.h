/**
 * Simulation.h
 *
 * Driver for the simulation. Takes in a SimulationArgs object and creates the
 * the necessary Box type, state file output path, etc.
 */

#ifndef SIMULATION_H
#define SIMULATION_H

#include "SimulationArgs.h"
#include "Box.h"
#include "Utilities/Logger.h"
#include "SimBox.h"
#include "Utilities/FileUtilities.h"

#define OUT_INTERVAL 100

class Simulation
{
  public:
    /** Construct the Simulation with configuration arguments */
    Simulation(SimulationArgs simArgs);

    /** Destruct the simulation */
    ~Simulation();

    /** Execute the simulation */
    void run();

  private:
    /** The periodic box the simulation run in */
    Box *box;

    /** The aguments that configure the simualtion */
    SimulationArgs args;

    /** The number of steps in the simulation */
    long simSteps;

    /** The starting step number for the simulation */
    long stepStart;

    /** The object that logs simulation events */
    Logger log;

    /** Pointer to the SB Scanner used to create the box */
    SBScanner *sbScanner;

    /** Writes the final state of the application to a PDB file */
    int writePDB(Environment sourceEnvironment,
                 Molecule *sourceMoleculeCollection, SimBox* sb);

    /** Saves the state of the simulation to a file */
    void saveState(const std::string& simName, int simStep, const SimBox* sb);

    /** The current date */
    const std::string currentDateTime();
};

#endif
