/**
 * SimulationArgs.h
 * 
 * Contains the definitions for the data structures and data types used to
 * communicate command-line information to the simulation.
 */

#ifndef METROSIM_SIMULATION_ARGS_H
#define METROSIM_SIMULATION_ARGS_H

#include <string>

/** Contains SimulationModeType enumeration */
namespace SimulationMode {
  /**
   * Specifies whether to run the simulation on the CPU or the GPU.
   *
   * @note If default is specified, then the program will utilize device
   *   querying to decide whether to run in serial or parallel.
   */
  enum Type {
    Default,
    Serial,
    Parallel
  };
}

/** Allows easy access to the SimulationMode::Type enumeration. */
typedef SimulationMode::Type SimulationModeType;

/** Enumeration for the type of input file */
namespace InputFile {
  enum Type {
    Unknown,
    Configuration,
    State
  };
}

typedef InputFile::Type InputFileType;

/** Enumeration for the simulation strategy */
namespace Strategy {
  enum Type {
    Default,
    BruteForce,
    ProximityMatrix,
    NLC,
    Unknown
  };

  Type fromString(std::string type);
}
typedef Strategy::Type SimulationStrategy;

/**
 * A list of commands and arguments that define the settings for the
 * simulation set by the user.
 * 
 * @remarks This is the data structure that will pass information read
 *   from the command line to the simulation.
 */
struct SimulationArgs {
  /**
   * A relative file path to the input file specified by the user.
   * @note This file path may lead to different types of files, such
   *     as configuration files or state files.
   */
  std::string filePath;

  /** The type of input file given to the program by the user. */
  InputFileType fileType;

  /** An optional name given to the current simulation run by the user. */
  std::string simulationName;

  /** Determines whether to run on the CPU or the GPU. */
  SimulationModeType simulationMode;

  /** Determines the particular strategy used for energy calculations */
  SimulationStrategy strategy;

  /**
   * If executing in parallel, the index of the graphics card being
   * used to run the simulation. Set to DEVICE_ANY if running 
   * using the CPU.
   */
  int deviceIndex;

  /** The number of simulation steps to execute in the current run. */
  int stepCount;

  /** If true, prints to cout (silenced for integration tests) */
  bool verboseOutput;

  /** Enables or disables use of linked-list neigborlist in calcs */
  bool useNeighborList;

  /** The number of simultion steps between updating the neighborlist */
  int neighborListInterval;

  /**
   * The number of simulation steps between status updates printed to
   * the console. A value of 0 means that status updates are only
   * printed at the beginning and the end of the simulation.
   */
  int statusInterval;

  /**
   * The number of simulation steps between state file saves.
   * @note This interval must be a non-negative number, and specifying
   *    an interval of 0 means a single state file should be saved
   *    at the very end of the simulation (which is the default
   *    behavior with no interval specified).
   */
  int stateInterval;

  /** Path to directory for PDB output data */
  std::string pdbOutputPath;

  /** Path to directory for simulation state data */
  std::string stateOutputPath;
};

#endif
