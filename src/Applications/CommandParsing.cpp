/// @file CommandParsing.cpp
///
/// Contains definitions for methods that read and parse command line
/// arguments and output a list of Simulation settings and values.
///
/// @author Tavis Maclellan
/// @date Created 2/23/2014
/// @date Updated 2/26/2014 -> further shoehorning done by Albert Wallace on 27 February


#include <getopt.h>
#include <iostream>
#include <string>

#include "CommandParsing.h"
#include "Metropolis/SimulationArgs.h"
#include "Metropolis/Utilities/DeviceQuery.h"
#include "Metropolis/Utilities/Parsing.h"

using std::string;

#define LONG_NAME 400


	//bool metrosim::getCommands(int argc, char** argv, SimulationArgs* args) //RBAl
	bool getCommands(int argc, char** argv, SimulationArgs* args)
	{
		CommandParameters params = CommandParameters();
		if (!readCommandLine(argc, argv, &params))
		{
			return false;	/* Incorrect syntax for command line input */
		}

		if (!parseCommandLine(&params, args))
		{
			return false;	/* Incorrect semantics of command line input */
		}

		return true;
	}

	//bool metrosim::readCommandLine(int argc, char** argv, CommandParameters* params) //RBAl
	bool readCommandLine(int argc, char** argv, CommandParameters* params)
	{
		// The getopt variables
		int getopt_ret, long_index;

		// Disable default getopt error messages. We will handle errors manually
		opterr = 0;

		// The short options recognized by the program
		const char* short_options = ":i:I:n:d:sphQV";

		// The long options recognized by the program
		struct option long_options[] = 
		{
			{"status-interval", 	required_argument, 	0, 	'i'},
			{"state-interval",		required_argument,	0,	'I'},
			{"steps",				required_argument,	0,	'n'},
			{"serial", 				no_argument, 		0, 	's'},
			{"parallel", 			no_argument, 		0, 	'p'},
			{"help", 				no_argument, 		0, 	'h'},
			{"list-devices",		no_argument,		0,	'Q'},
			{"device",				required_argument,	0,	'd'},
			{"version",				no_argument,		0,	'V'},
			{"name",				required_argument,	0,	LONG_NAME},
			{0, 0, 0, 0} 
		};

		// Iterate over all command-line arguments and match any option entries.
		while ((getopt_ret = getopt_long(argc, argv, short_options, long_options, &long_index)) != -1)
		{ 	
 			// Switch the returned option value
			switch (getopt_ret)
			{
				case 0:		/* long options */
					if (long_options[long_index].flag != 0)
						break;
					/* handle other long options without return values */
					break;
				case 'V':
					printVersionInformation();
					return false;
					break;
				case LONG_NAME:
					if (!fromString<string>(optarg, params->simulationName))
					{
						std::cerr << APP_NAME << ": ";
						std::cerr << " --name: Invalid simulation name" << std::endl;
						return false;
					}
					break;
				case 'i':	/* status interval */
					params->statusFlag = true;
					if (!fromString<int>(optarg, params->statusInterval))
					{
						std::cerr << APP_NAME << ": ";
						std::cerr << " --status_interval (-i): Invalid status interval" << std::endl;
						return false;
					}
					if (params->statusInterval < 0)
					{
						std::cerr << APP_NAME << ": ";
						std::cerr << " --status_interval (-i): Status Interval must be non-negative" << std::endl;
						return false;
					}
					break;
				case 'I':	/* state interval */
					params->stateFlag = true;
					if (!fromString<int>(optarg, params->stateInterval))
					{
						std::cerr << APP_NAME << ": ";
						std::cerr << " --state_interval (-I): Invalid state interval" << std::endl;
						return false;
					}
					break;
				case 'n':	/* simulation steps */
					params->stepFlag = true;
					if (!fromString<int>(optarg, params->stepCount))
					{
						std::cerr << APP_NAME << ": ";
						std::cerr << " --steps (-n): Invalid step count" << std::endl;
						return false;
					}
					if (params->stepCount <= 0)
					{
						std::cerr << APP_NAME << ": ";
						std::cerr << " --steps (-n): Step count must be greater than zero" << std::endl;
						return false;
					}
					break;
				case 's':	/* run serial */
					params->serialFlag = true;
					break;
				case 'p':	/* run parallel */
					params->parallelFlag = true;
					break;
				case 'h':	/* print help */
					printHelpScreen();
					return false;
					break;
				case 'Q':	/* list all devices */
					printDeviceInformation();
					return false;
					break;
				case 'd':
					params->deviceFlag = true;
					if (!fromString<int>(optarg, params->deviceIndex))
					{
						std::cerr << APP_NAME << ": ";
						std::cerr << " --device (-d): Device index must be a valid integer" << std::endl;
						return false;
					}
					if (params->deviceIndex < 0)
					{
						std::cerr << APP_NAME << ": ";
						std::cerr << " --device (-d): Device index cannot be a negative index" << std::endl;
						return false;
					}
					break;
				case ':':	/* missing argument */
					if (sizeof(argv[optind-1]) > 2 && argv[optind-1][0] == '-' && argv[optind-1][1] == '-')
					{
						std::cerr << APP_NAME << ": A required argument is missing for ";
						std::cerr << argv[optind-1] << std::endl;
					}
					else
					{
						std::cerr << APP_NAME << ": A required argument is missing for ";
						std::cerr << "-" << (char) optopt << std::endl;
					}
					return false;
					break;
				case '?':	/* unknown option */
					if (optopt)
					{
						std::cerr << APP_NAME << ": Unknown option -" << (char) optopt << std::endl;
					}
					else if ((optind - 1) < argc)
					{
						std::cerr << APP_NAME << ": Unknown option " << argv[optind-1] << std::endl;
					}
					else
					{
						std::cerr << APP_NAME << ": Unknown option not recognized" << std::endl;
					}
					return false;
					break;
				default:	/* unknown error */
					std::cerr << APP_NAME << ": Fatal error occurred while parsing commnd-line" << std::endl;
					return false;
					break;
			}
		}

		// Get the count and list of non-option arguments specified by the user.
		params->argCount = argc - optind;
		params->argList = &(argv[optind]);

		return true;
	}

	//bool metrosim::parseCommandLine(CommandParameters* params, SimulationArgs* args) //RBAl
	bool parseCommandLine(CommandParameters* params, SimulationArgs* args)
	{
		if (params->argCount < 1)
		{
			std::cerr << APP_NAME << ": Input file not specified" << std::endl;
			return false;
		}
		
		if (params->argCount > 1)
		{
			std::cerr << APP_NAME << ": Too many input files specified" << std::endl;
			return false;
		}

		if (params->argList == NULL)
		{
			std::cerr << APP_NAME << ": Must specify a valid file path" << std::endl;
			return false;
		}

		if (params->serialFlag && params->parallelFlag) /* conflicting flags */
		{
			std::cerr << APP_NAME << ": Cannot specify both GPU and CPU modes" << std::endl;
			return false;
		}

		// Assign the relevant information that will be used in the simulation
		// to the command arguments container.

		if (params->parallelFlag)	/* parallel flag set */
		{
			args->simulationMode = SimulationMode::Parallel;
		}
		else if (params->serialFlag) /* serial flag set */
		{
			args->simulationMode = SimulationMode::Serial;
		}

		if (!parseInputFile(params->argList[0], args->filePath, args->fileType))
		{
			std::cerr << APP_NAME << ": Must specify a config or state file" << std::endl;
			return false;
		}

		args->statusInterval = params->statusInterval;
		args->stateInterval = params->stateInterval;
		args->stepCount = params->stepCount;
		args->simulationName = params->simulationName;

		if (!params->parallelFlag && params->deviceFlag)
		{
			std::cerr << APP_NAME << ": Cannot use graphics device in serial simulation" << std::endl;
			return false;
		}

		// Set the user selected device index.
		if (params->deviceFlag)
			args->deviceIndex = params->deviceIndex;
		else
			args->deviceIndex = DEVICE_ANY;

		return true;
	}

	bool parseInputFile(char* filename, string& name, InputFileType& type)
	{
		if (!filename)
			return false;
		if (!fromString<string>(filename, name))
			return false;
		if (name.empty())
		{
			std::cerr << APP_NAME << ": Input file must be a valid path" << std::endl;
			return false;
		}

		std::string extension = getExtension(name);
		if (extension == "config")
		{
			type = InputFile::Configuration;
			return true;
		}
		else if (extension == "state")
		{
			type = InputFile::State;
			return true;
		}

		type = InputFile::Unknown;
		return false;
	}

	//void metrosim::printHelpScreen() //RBAl
	void printHelpScreen()
	{
		using std::endl;
		using std::cout;

		cout << endl;
		cout << "Usage: " << APP_NAME << " <inputfile> [options]" << endl << endl;

		cout << "Input File Types\n"
			  "=================\n";
		cout << "The file extension of the input file will determine what phase of the\n"
			  "simulation to execute:\n";
		cout << "\t.config\t: Create a fresh simulation using the given set of\n"
			  "\t\t  configuration settings and parameters.\n";
		cout << "\t.state\t: Resume a previous simulation that was saved.\n";
		cout << "\n";

		cout << "GPU Operation Flags\n"
			  "====================\n";
		cout << "This tool can execute a given simulation on both the CPU and the GPU.\n"
			  "The following flags can query the current availability of GPU devices\n"
			  "and control how the simulation will execute on the GPU.\n\n";
		cout << "--parallel\t(-p)\n";
		cout << "\tRun the simulation in parallel by executing steps on a CUDA\n"
			  "\tcapable device. The program will terminate if no available devices\n"
			  "\tare found or if the tool's minimum specifications are not met.\n"
			  "\tIf you specify this flag you cannot also specify the --serial\n"
			  "\tflag.\n\n";
		cout << "--serial\t(-s)\n";
		cout << "\tRun the simulation in serial on the host CPU. If you specify this\n"
			  "\tflag you cannot also specify the --parallel flag.\n\n";
		cout << "--list-devices\t(-Q)\n";
		cout << "\tQueries the host machine for available CUDA capable devices. The\n"
			  "\tentries will contain the following information about each\n"
			  "\tavailable device:\n\n";
		cout << "\tDevice Index\t: The index of the device on the host system\n";
		cout << "\tDevice Name\t: The name of the available device.\n";
		cout << "\tClock Rate\t: The clock rate of the device in kilohertz\n";
		cout << "\tCompute Capability : The major and minor version of the device.\n";
		cout << "\tProcessor Count\t: The number of multiprocessors on the device.\n";
		cout << "\tGlobal Memory\t: The total amount of global memory on the device.\n";
		cout << "\tConstant Memory\t: The total amount of constant memory on the device\n";
		cout << "\tWarp Size\t: The number of threads in a warp.\n";
		cout << "\n";
		cout << "--device <index>\t(-d)\n";
		cout << "\tWhen running the simulation on the GPU, this option specifies which\n"
				"\tavailable graphics device to use if multiple are installed on the\n"
				"\thost system. If this option is not specified then the simulation\n"
				"\twill attempt to find the best device. Additionally, this flag must\n"
				"\tbe accompanied with the --parallel flag.\n\n";

		cout << "Simulation Run Parameters\n"
			  "==========================\n";
		cout << "These options define simulation settings and parameters that specify\n"
				"how to the run will execute.\n\n";
		cout << "--name <title>\n";
		cout << "\tSpecifies the name of the simulation that will be run. This name\n"
				"\twill be used as the basename for saved state files, as well as\n"
				"\tthe title of the simulation results.\n\n";
		cout << "--steps <count>\t\t(-n)\n";
		cout << "\tSpecifies how many simulation steps to execute in the Monte\n"
				"\tCarlo Metropolis algorithm. This value must a valid integer\n"
				"\tthat is greater than zero.\n\n";
		cout << "--status-interval <interval>\t(-i)\n";
		cout << "\tSpecifies the number of simulation steps between status updates.\n"
				"\tThese status updates will periodically be printed out that list\n"
				"\tthe current step number and the current total system energy.\n"
				"\tThis interval must be a valid non-negative integer value. An\n"
				"\tinterval of 0 (zero) means that status updates should only\n"
				"\tprint out at the beginning and the end of the simulation.\n\n";
		cout <<"--state-interval <interval>\t(-I)\n";
		cout << "\tSpecifies the number of simulation steps between state file\n"
				"\tsnapshots of the current simulation run. The program will\n"
				"\tperiodically save the simulation to a state file (*.state)\n"
				"\tthat can be resumed later by the user. This provides checkpoints\n"
				"\tjust in case a long running simulation is interrupted.\n\n";
		cout << "\tSpecifying a state interval of zero, or not specifying the\n"
				"\toption at all, means that only a final state file should be\n"
				"\twritten at the end of the simulation run. An interval greater\n"
				"\tthan zero means that multiple state files shouls be written\n"
				"\tat the specified intervals. Finally, a state interval less than\n"
				"\tzero means that no state files should be written at all.\n\n";
		cout << "\tThe filename for the saved state files will be the simulation\n"
				"\tname with the current step number appended at the end:\n\n";
		cout << "\t\t<simulation-name>_<step-num>.state\n\n";

		cout << "Generic Tool Options\n"
			  "=====================\n";
		cout << "\n";
		cout << "--help\t\t(-h)\n";
		cout << "\tPrint this help information to the standard output stream.\n\n";
		cout << "--version\t(-V)\n";
		cout << "\tPrint version information of this tool to the standard output\n"
			  "\tstream.\n\n";
	}

	void printVersionInformation()
	{
		std::cout << APP_NAME << ": MCGPU Monte Carlo Simulator" << std::endl;
		#ifdef DEBUG
		std::cout << "Current Build: Debug" << std::endl;
		#else
		std::cout << "Current Build: Release" << std::endl;
		#endif
		#ifdef DOUBLE_PRECISION
		std::cout << "Floating-point Precision: Double" << std::endl;
		#else
		std::cout << "Floating-point Precision: Single" << std::endl;
		#endif
		std::cout << "Minimum Required CUDA Version: " << MIN_MAJOR_VER << "." << MIN_MINOR_VER << std::endl;
	}