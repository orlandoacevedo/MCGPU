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

#define LONG_VERSION 300
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
		const char* short_options = ":i:I:n:sphQ";

		// The long options recognized by the program
		struct option long_options[] = 
		{
			{"status_interval", 	required_argument, 	0, 	'i'},
			{"state_interval",		required_argument,	0,	'I'},
			{"steps",				required_argument,	0,	'n'},
			{"serial", 				no_argument, 		0, 	's'},
			{"parallel", 			no_argument, 		0, 	'p'},
			{"help", 				no_argument, 		0, 	'h'},
			{"list_devices",		no_argument,		0,	'Q'},
			{"version",				no_argument,		0,	LONG_VERSION},
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
				case LONG_VERSION:
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
		std::cout << "Usage: " << APP_NAME << ": ";
		std::cout << "<config_file> [-i interval] [-sp] [-fd] [-hQ]" << std::endl << std::endl;
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