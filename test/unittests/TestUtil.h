#ifndef TEST_UTIL_H
#define TEST_UTIL_H

#include<string>
#include<sstream>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<unistd.h>

struct ConfigFileData {
	double sizeX;
	double sizeY;
	double sizeZ;
	double tempKelvin;
	double maxTranslation;
	int numSteps;
	int numMolecules;
	std::string oplsaa_par_path;
	std::string working_path;
	std::string z_matrix_file;
	double cutoffDistance;
	double maxRotation;
	int randomSeed;

	/**
	 * Constructor for ConfigFileData.
	 * @param sizeX_in The length of the box on the x-axis in angstroms.
	 * @param sizeY_in The length of the box on the y-axis in angstroms.
	 * @param sizeZ_in The length of the box on the z-axis in angstroms.
	 * @param tempKelvin_in The temperature of the box, in kelvin.
	 * @param maxTranslation_in The maximum distance a molecule can translate in one step (measured in angstroms).
	 * @param numSteps_in The number of steps to run the simulation for.
	 * @param numMolecules_in The number of molecules in the box.
	 * @param oplsaa_par_path_in The path to the oplsaa.par constants file.
	 * @param working_path_in The path to the z-matrix, config, and state files.
	 * @param z_matrix_file_in The path to and name of the z-matrix file.
	 * @param cutoffDistance_in The cutoffDistance to use (in angstroms).
	 * @param maxRotation_in The maximum amount a molecule can rotate in one step (measured in degrees).
	 * @param randomSeed_in The random seed to use.
	 */
	ConfigFileData(double sizeX_in, double sizeY_in, double sizeZ_in, double tempKelvin_in, double maxTranslation_in,
				   int numSteps_in, int numMolecules_in, std::string oplsaa_par_path_in, std::string z_matrix_file_in,
				   std::string working_path_in, double cutoffDistance_in, double maxRotation_in, int randomSeed_in) {

		sizeX = sizeX_in;
		sizeY = sizeY_in;
		sizeZ = sizeZ_in;
		tempKelvin = tempKelvin_in;
		maxTranslation = maxTranslation_in;
		numSteps = numSteps_in;
		numMolecules = numMolecules_in;
		oplsaa_par_path = oplsaa_par_path_in;
		working_path = working_path_in;
		z_matrix_file = z_matrix_file_in;
		cutoffDistance = cutoffDistance_in;
		maxRotation = maxRotation_in;
		randomSeed = randomSeed_in;
	}
};

/**
 * Generates a configuration file to be used in testing.
 * @param MCGPU_path The path to MCGPU's root directory.
 * @param fileName The name of the configuration file to generate.
 * @param primaryAtomIndexString The entry for the Primary Atom Index line of the config file.
 * @param settings The specifics of the test type's config file.
 */
void createConfigFile(std::string MCGPU, std::string fileName, std::string primaryAtomIndexString, ConfigFileData settings);

/**
 * Returns the path to MCGPU's root directory.
 */
std::string getMCGPU_path ();


/**
 * Generates a string that is used as a command to run the necessary test.
 * @param MCGPU The path to MCGPU's root.
 * @param configFile The name of the config file to use for the test.
 * @param outputName The name of the simulation, or the file to pipe cerr to (if expecting an error).
 * @param series true if the simulation is to be run in series, false otherwise.
 * @param neighborlist true if the simulation is to be run with a neighborlist, false otherwise.
 * @param errorExpected true if passing behavior for the test throws an error.
 * @param working_path The path in which the config file is located.
 */
std::string buildCommand(std::string MCGPU, std::string configFile, std::string outputName, bool series, bool neighborlist, bool errorExpected, std::string working_path);

/**
 * Examines a result file to determine the final energy. Results file is
 * assumed to be in the current working directory.
 * @param resultsFile The name of the results file.
 * @return The final energy.
 */
double getEnergyResult(std::string resultsFile);

/**
 * Examines a file that cerr was piped to and returns what function call the error occurred in.
 * @param MCGPU The path to MCGPU's root.
 * @param errorFile The name of the file cerr was piped to.
 * @return The function the error occurred in.
 */
std::string getErrorResult(std::string MCGPU, std::string errorFile);

/** Simple wrapper to obtain the current working directory */
std::string get_cwd();

#endif
