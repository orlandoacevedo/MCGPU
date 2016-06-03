/** @file CommandParsing.h
 *
 * Contains definitions for data structures and functions that parse command
 * line input into the application.
 *
 * @author Tavis Maclellan
 * @date Created 2/23/2014
 * @date Updated 1/7/2016
 */

#ifndef METROSIM_COMMAND_PARSING_H
#define METROSIM_COMMAND_PARSING_H

#include <string>
#include "Metropolis/SimulationArgs.h"

#ifndef APP_NAME
#define APP_NAME "metrosim"
#endif

#define DEFAULT_STATUS_INTERVAL 1000
#define DEFAULT_NEIGHBORLIST_INTERVAL 100

/**
 * Contains the intermediate values and flags read in from the command
 * line.
 *
 * @note This data structure will have to be updated when new commands
 * are added.
 */
struct CommandParameters {

  /**
   * The number of simulation steps between status updates.
   * @note This interval must be a non-negative number, and specifying
   *     an interval of 0 means status updates will occur only
   *     at the beginning and end of the simulation.
   */
  int statusInterval;

  /**
   * The number of simulation steps between state file saves.
   *
   * @note This interval must be a non-negative number, and specifying
   *    an interval of 0 means a single state file should be saved
   *    at the very end of the simulation (which is the default
   *    behavior with no interval specified).
   */
  int stateInterval;

  /**
   * The number of simulation steps to execute in the current run.
   * This must be a valid integer number greater than zero.
   */
  int stepCount;

  /** Declares whether the help option was specified. */
  bool helpFlag;

  /** Declares whether the version option was specified. */
  bool versionFlag;

  /** The number of non-option arguments given by the user. */
  int argCount;

  /** The list of non-option arguments given by the user. */
  char** argList;

  /** The optional name of the simulation run. */
  std::string simulationName;

  /** Declares whether the status interval option was specified. */
  bool statusFlag;

  /** Declare whether the state interval option was specified */
  bool stateFlag;

  /** Declares whether the steps option was specified.  */
  bool stepFlag;

  /** Declares whether the serial execution option was specified. */
  bool serialFlag;

  /** Declares whether the parallel execution option was specified. */
  bool parallelFlag;

  /**
   * Declares whether or not to display the printed run
   * information to standard cout.
   */
  bool verboseOutputFlag;

  /** Declares whether the use of the Neighbor-list is enabled. */
  bool neighborListFlag;

  /**
   * The number of simulation steps between updating the neighborlist.
   * @note This interval must be a non-negative number, and specifying
   *     an interval of 0 means the initial neighborlist will not be updated
   */
  int neighborListInterval;

  /** The simulation strategy specified by the user */
  std::string simStrategy;

  /** Default constructor */
  CommandParameters() : statusInterval(DEFAULT_STATUS_INTERVAL),
              stateInterval(0),
              stepCount(0),
              argCount(0),
              argList(NULL),
              helpFlag(false),
              versionFlag(false),
              statusFlag(false),
              stepFlag(false),
              serialFlag(false),
              parallelFlag(false),
              verboseOutputFlag(false),
              neighborListFlag(false),
              neighborListInterval(DEFAULT_NEIGHBORLIST_INTERVAL)   {}
};

/**
 * Goes through each argument specified from the command line and checks
 * for usage errors, then fills a list of simulation settings and values
 * that will be used by the application.
 *
 * @param[in] argc The number of entries in the arguments list.
 * @param[in] argv The list of command line arguments given to the program.
 * @param[out] args The list of settings that define the behavior and
 *     values for the simulation to use.
 * @returns True if the simulation should be executed; false if the
 *     simulation should not be executed.
 */
bool getCommands(int argc, char** argv, SimulationArgs* args);

/**
 * Analyzes the syntax of the command line arguments and fills a parameter
 * list with flags and raw values.
 *
 * @param[in] argc The number of entries in the arguments list.
 * @param[in] argv The list of command line arguments given to the program.
 * @param[out] params The list of command line parameters that contain
 *     intermediate values ready to be parsed.
 * @returns True if no errors occurred; false if a syntax error occurred.
 */
bool readCommandLine(int argc, char** argv, CommandParameters* params);

/**
 * Analyzes the semantics of the command line parameters and fills the
 * simulation arguments list with settings and values.
 *
 * @param[in] params The list of command line parameters that contain
 *     intermediate values ready to be parsed.
 * @param[out] args The list of settings that define the behavior and
 *     values for the simulation to use.
 * @returns True if no errors occurred and the simulation shoud execute;
 *     false if a semantics error occurred or the simulation was halted.
 * @remarks Currently if the user specifies the --help or -h options via
 *     the command line input to the program, the simulation will not
 *     execute, and instead usage information will be printed to the user
 *     and the program will terminate. This complies with GNU conventions.
 */
bool parseCommandLine(CommandParameters* params, SimulationArgs* args);

/**
 * Analyzes a file name and determines its type.
 *
 * @param[in] filename The name of the file, as specified at the command line.
 * @param[out] name The filename, converted to std::string, if applicable.
 * @param[out] type The type of the file, either InputFile::Configuration,
 *     InputFile::State, or InputFile::Unknown
 * @returns True if the filename is valid, false otherwise.
 */
bool parseInputFile(char* filename, std::string& name, InputFileType& type);

/**
 * Outputs the help documentation to the standard output stream and
 * displays how to use the application.
 */
void printHelpScreen();

/**
 * Outputs the current program version and build information to the
 * standard output stream.
 */
void printVersionInformation();

#endif
