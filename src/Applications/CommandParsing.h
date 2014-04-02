/// @file CommandParsing.h
///
/// Contains definitions for data structures and functions that parse command
/// line input into the application.
///
/// @author Tavis Maclellan
/// @date Created 2/23/2014
/// @date Updated 2/26/2014

#ifndef METROSIM_COMMAND_PARSING_H
#define METROSIM_COMMAND_PARSING_H

#include <stdlib.h>
#include "Metropolis/SimulationArgs.h"


#ifndef APP_NAME
#define APP_NAME "metrosim"		///< The executable name of the application.
#endif

#define DEFAULT_STATUS_INTERVAL 100

namespace metrosim
{

	/// Contains the intermediate values and flags read in from the command
	/// line.
	///
	/// @note This data structure will have to be updated when new commands
	/// are added.
	struct CommandParameters
	{
		/// The number of simulation steps between status updates.
		/// @note This interval must be a non-negative number, and specifying
		///     an interval of 0 means status updates will occur only
		///     at the beginning and end of the simulation.
		int statusInterval;

		/// Declares whether the help option was specified.
		bool helpFlag;

		/// Declares whether the query device option was specified.
		bool listDevicesFlag;

		/// The number of non-option arguments given by the user.
		int argCount;

		/// The list of non-option arguments given by the user.
		char** argList;

		/// Declares whether the status option was specified.
		bool statusFlag;

		/// Declares whether the serial execution option was specified.
		bool serialFlag;

		/// Declares whether the parallel execution option was specified.
		bool parallelFlag;

		/// Declares whether the single-precision option was specified.
		bool floatFlag;

		/// Declares whether the double-precision option was specified.
		bool doubleFlag;

		/// Default constructor
		CommandParameters() :	statusInterval(DEFAULT_STATUS_INTERVAL),
								argCount(0),
								argList(NULL),
								helpFlag(false),
								listDevicesFlag(false),
								statusFlag(false), 
								serialFlag(false),
								parallelFlag(false), 
								floatFlag(false), 
								doubleFlag(false) {}
	};

	/// Goes through each argument specified from the command line and checks
	/// for usage errors, then fills a list of simulation settings and values
	/// that will be used by the application.
	///
	/// @param[in] argc The number of entries in the arguments list.
	/// @param[in] argv The list of command line arguments given to the program.
	/// @param[out] args The list of settings that define the behavior and
	///     values for the simulation to use.
	/// @returns True if the simulation should be executed; false if the
	///     simulation should not be executed.
	bool getCommands(int argc, char** argv, SimulationArgs* args);

	/// Analyzes the syntax of the command line arguments and fills a parameter
	/// list with flags and raw values.
	///
	/// @param[in] argc The number of entries in the arguments list.
	/// @param[in] argv The list of command line arguments given to the program.
	/// @param[out] params The list of command line parameters that contain
	///     intermediate values ready to be parsed.
	/// @returns True if no errors occurred; false if a syntax error occurred.
	bool readCommandLine(int argc, char** argv, CommandParameters* params);

	/// Analyzes the semantics of the command line parameters and fills the
	/// simulation arguments list with settings and values.
	///
	/// @param[in] params The list of command line parameters that contain
	///     intermediate values ready to be parsed.
	/// @param[out] args The list of settings that define the behavior and
	///     values for the simulation to use.
	/// @returns True if no errors occurred and the simulation shoud execute; 
	///     false if a semantics error occurred or the simulation was halted.
	/// @remarks Currently if the user specifies the --help or -h options via
	///     the command line input to the program, the simulation will not
	///     execute, and instead usage information will be printed to the user
	///     and the program will terminate. This complies with GNU conventions.
	bool parseCommandLine(CommandParameters* params, SimulationArgs* args);

	/// Attempts to parse a config filepath from a command line argument.
	///
	/// @param[out] dest The destination location that will store the
	///     parsed config filepath.
	/// @param[in] arg The command line argument to parse.
	/// @returns True if no errors occurred; false if an error occurred while
	///     parsing the config filepath.
	bool parseConfigPath(char** dest, char* arg);

	/// Attempts to parse a status interval from a command line argument.
	///
	/// @param[out] dest The destination location that will store the parsed
	///     status interval.
	/// @param[in] arg The command line argument to parse.
	/// @returns True if no errors occurred; false if an error occurred while
	///     parsing the status interval.
	bool parseStatusInterval(int* dest, char* arg);

	/// Outputs the help documentation to the standard output stream and
	/// displays how to use the application.
	void printHelpScreen();
}

#endif