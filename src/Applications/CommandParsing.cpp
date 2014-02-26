/// @file CommandParsing.cpp
///
/// Contains definitions for methods that read and parse command line
/// arguments and output a list of Simulation settings and values.
///
/// @author Tavis Maclellan
/// @date Created 2/23/2014
/// @date Updated 2/26/2014

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include "CommandParsing.h"
#include "Metropolis/SimulationArgs.h"

namespace metrosim
{

	bool getCommands(int argc, char** argv, SimulationArgs* args)
	{
		CommandParameters params;
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

	bool readCommandLine(int argc, char** argv, CommandParameters* params)
	{
		// The getopt variables
		int getopt_ret, long_index;

		// Disable default getopt error messages. We will handle errors manually
		opterr = 0;

		// The short options recognized by the program
		const char* short_options = ":i:spfdh";

		// The long options recognized by the program
		struct option long_options[] = {
			{"status_interval", 	required_argument, 	0, 	'i'},
			{"serial", 				no_argument, 		0, 	's'},
			{"parallel", 			no_argument, 		0, 	'p'},
			{"single_precision",  	no_argument, 		0, 	'f'},
			{"double_precision", 	no_argument, 		0, 	'd'},
			{"help", 				no_argument, 		0, 	'h'},
			{0, 0, 0, 0} };

		// Iterate over all command-line arguments and match any option entries.
		while ((getopt_ret = getopt_long(argc, argv, short_options, long_options, &long_index)) != -1)
		{ 	
 			// Switch the returned option value
			switch (getopt_ret)
			{
				case 0:		/* long option */
					if (long_options[long_index].flag != 0)
						break;
					/* Handle long options that are not automatically handled */
					break;
				case 'i':	/* status interval */
					params->statusFlag = true;
					if (!parseStatusInterval(&(params->statusInterval), optarg))
						return false;
					break;
				case 's':	/* run serial */
					params->serialFlag = true;
					break;
				case 'p':	/* run parallel */
					params->parallelFlag = true;
					break;
				case 'f':	/* use floats */
					params->floatFlag = true;
					break;
				case 'd':	/* use doubles */
					params->doubleFlag = true;
					break;
				case 'h':	/* print help */
					params->helpFlag = true;
					break;
				case ':':	/* missing argument */
					if (sizeof(argv[optind-1]) > 2 && argv[optind-1][0] == '-' && argv[optind-1][1] == '-')
					{
						fprintf(stderr, "%s: A required argument is missing for %s\n", APP_NAME, argv[optind-1]);
					}
					else
					{
						fprintf(stderr, "%s: A required argument is missing for -%c\n", APP_NAME, optopt);
					}
					return false;
					break;
				case '?':	/* unknown option */
					if (optopt)
					{
						fprintf(stderr, "%s: Unknown option -%c\n", APP_NAME, optopt);
					}
					else if ((optind - 1) < argc)
					{
						fprintf(stderr, "%s: Unknown option %s\n", APP_NAME, argv[optind-1]);
					}
					else
					{
						fprintf(stderr, "%s: Unknown option not recognized\n", APP_NAME);
					}
					return false;
					break;
				default:	/* unknown error */
					fprintf(stderr, "%s: Fatal error occurred while parsing command-line\n", APP_NAME);
					return false;
					break;
			}
		}

		// Get the count and list of non-option arguments specified by the user.
		params->argCount = argc - optind;
		params->argList = &(argv[optind]);

		return true;
	}

	bool parseCommandLine(CommandParameters* params, SimulationArgs* args)
	{
		if (params->helpFlag)		/* print help screen and exit */
		{
			printHelpScreen();
			return false;
		}

		if (params->argCount < 1)
		{
			fprintf(stderr, "%s: input config file not specified.\n", APP_NAME);
			return false;
		}
		
		if (params->argCount > 1)
		{
			fprintf(stderr, "%s: too many input files specified.\n", APP_NAME);
			return false;
		}

		if (params->argList == NULL)
		{
			fprintf(stderr, "%s: must specify a valid config file path.\n", APP_NAME);
			return false;
		}

		if (params->serialFlag && params->parallelFlag) /* conflicting flags */
		{
			fprintf(stderr, "%s: cannot specify both GPU and CPU modes.\n", APP_NAME);
			return false;
		}

		if (params->floatFlag && params->doubleFlag)	/* conflicting flags */
		{
			fprintf(stderr, "%s: cannot specify both float and double precisions.\n", APP_NAME);
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

		if (params->floatFlag)		/* float flag set */
		{
			args->precisionMode = PrecisionMode::Single;
		}
		else if (params->doubleFlag) /* double flag set */
		{
			args->precisionMode = PrecisionMode::Double;
		}

		if (!parseConfigPath(&(args->configPath), params->argList[0]))
		{
			return false;
		}

		args->statusInterval = params->statusInterval;

		return true;
	}

	bool parseConfigPath(char** dest, char* arg)
	{
		if (arg == NULL || arg[0] == '\0')
		{
			fprintf(stderr, "%s: config filepath must be a valid filename.\n", APP_NAME);
			return false;
		}

		*dest = arg;
		return true;
	}

	bool parseStatusInterval(int* dest, char* arg)
	{
		if (arg == NULL) /* argument is invalid */
		{
			fprintf(stderr, "%s: a valid status interval must be specified.\n", APP_NAME);
			return false;
		}

		char* temp;

		// Parse the string to an integer value
		int interval = strtol(arg, &temp, 10);

		if (temp == arg || *temp != '\0') /* could not parse argument */
		{
			fprintf(stderr, "%s: status interval must be a valid non-negative integer.\n", APP_NAME);
			return false;
		}
		else if (interval < 0)	/* less than zero */
		{
			fprintf(stderr, "%s: status interval must be greater than or equal to zero.\n", APP_NAME);
			return false;
		}
		else if (errno)			/* integer out of range */
		{
			fprintf(stderr, "%s: status interval out of range.\n", APP_NAME);
			return false;
		}

		*dest = interval;
		return true;
	}

	void printHelpScreen()
	{
		const char* usage = "Usage  : %s -c config_file [-i stat_interval] [-sp] [-fd] [-h]\n";
		fprintf(stdout, usage, APP_NAME);
	}
}