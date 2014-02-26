/// @file CommandArgs.h
///
/// Contains the definitions for the data structures and data types used to communicate
/// command-line information to the simulation.
///
/// @author Tavis Maclellan
/// @date Created 2/23/2014
/// @date Updated 2/26/2014

#ifndef METROSIM_COMMAND_ARGS_H
#define METROSIM_COMMAND_ARGS_H

/// Application namespace for metro simulation
namespace metrosim
{

	/// Contains SimulationModeType enum
	namespace SimulationMode
	{
		/// Specifies whether to run the simulation on the CPU or the GPU.
		///
		/// @note If default is specified, then the program will utilize device
		///   querying to decide whether to run in serial or parallel.
		enum Type
		{
			/// Run the simulation using a default mode.
			Default,

			/// Run the simulation serially on the CPU.
			Serial,

			/// Run the simulation in parallel on the GPU. The simulation
			/// will not run if a valid CUDA device is not found.
			Parallel
		};
	}

	/// Allows easy access to the SimulationMode::Type enumeration.
	typedef SimulationMode::Type SimulationModeType;

	/// Contains PrecisionModeType enum
	namespace PrecisionMode
	{

		/// Specifies the floating-point precision to use for all calculations
		/// used during the simulation.
		enum Type
		{
			/// Use the default precision defined by the simulation.
			Default,

			/// Use single-precision (32 bits) for calculations and data.
			Single,

			/// Use double-precision (64 bits) for calculations and data.
			Double
		};
	}

	/// Allows easy access to the PrecisionMode::Type enumeration.
	typedef PrecisionMode::Type PrecisionModeType;

	/// A list of commands and arguments that define the settings for the
	/// simulation set by the user.
	///
	/// @remarks This is the data structure that will pass information read
	///   from the command line to the simulation.
	struct CommandArguments
	{
		/// A relative filepath to a config file that contains simulation data and settings.
		char* configPath;

		/// The number of simulation steps between status updates printed to
		/// the console. A value of 0 means that status updates are only
		/// printed at the beginning and the end of the simulation.
		int statusInterval;

		/// The simulation mode that determines whether to run on the
		/// CPU or the GPU.
		SimulationModeType simulationMode;

		/// The precision to use for floating-point calculations and data.
		PrecisionModeType precisionMode;
	};
}

#endif