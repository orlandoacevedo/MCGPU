#include "Applications/Application.h"
#include "gtest/gtest.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

/**
 * Generates a configuration file to be used in testing.
 * @param MCGPU_path The path to MCGPU's root directory.
 * @param fileName The name of the configuration file to generate.
 * @param primaryAtomIndexString The entry for the Primary Atom Index line of the config file.
 */
void createConfigFile(std::string MCGPU, std::string fileName, std::string primaryAtomIndexString) {
	ofstream configFile;
	std::string configFilePath (std::string (MCGPU + "test/unittests/MultipleSolvents/t3pdimoneTests/" + fileName));
	configFile.open(configFilePath.c_str());
	std::stringstream cfg;
	cfg << ""
		<< "#size of periodic box (x, y, z in angstroms)\n"
		<< "26.15\n"
		<< "26.15\n"
		<< "26.15\n"
		<< "#temperature in Kelvin\n"
		<< "298.15\n"
		<< "#max translation\n"
		<< ".12\n"
		<< "#number of steps\n"
		<< "100000\n"
		<< "#number of molecules\n"
		<< "256\n"
		<< "#path to opla.par file\n"
		<< MCGPU << "resources/bossFiles/oplsaa.par\n"
		<< "#path to z matrix file\n"
		<< MCGPU << "test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimone.z\n"
		<< "#path to state input\n"
		<< MCGPU << "test/unittests/MultipleSolvents/t3pdimoneTests\n"
        << "#path to state output\n"
        << MCGPU << "test/unittests/MultipleSolvents/t3pdimoneTests\n"
        << "#pdb output path\n"
        << MCGPU << "test/unittests/MultipleSolvents/t3pdimoneTests\n"
        << "#cutoff distance in angstroms\n"
        << "11.0\n"
        << "#max rotation\n"
        << "12.0\n"
        << "#Random Seed Input\n"
        << "12345\n"
        << "#Primary Atom Index\n"
        << primaryAtomIndexString;
        configFile << cfg.str();
        configFile.close();
}

/**
 * Returns the path to MCGPU's root directory.
 */
std::string getMCGPU_path () {
	string directory = get_current_dir_name();
	std::string mc ("MCGPU");
	std::size_t found = directory.find(mc);
	
	if(found != std::string::npos) {
		directory = directory.substr(0, found + 6);
	}
	std::string MCGPU = directory + "/";
	return MCGPU;
}

/**
 * Generates a string that is used as a command to run the necessary test.
 * @param MCGPU The path to MCGPU's root.
 * @param configFile The name of the config file to use for the test.
 * @param outputName The name of the simulation, or the file to pipe cerr to (if expecting an error).
 * @param series true if the simulation is to be run in series, false otherwise.
 * @param neighborlist true if the simulation is to be run with a neighborlist, false otherwise.
 * @param errorExpected true if passing behavior for the test throws an error.
 */
std::string buildCommand(std::string MCGPU, std::string configFile, std::string outputName, bool series, bool neighborlist, bool errorExpected) {
	//Setting up standard build command.
	std::stringstream ss;
	ss << MCGPU << "/bin/metrosim " << MCGPU << "test/unittests/MultipleSolvents/t3pdimoneTests/" << configFile << " ";
	
	if(series) {
		ss << "-s --threads 12 ";	//If testing in series, give the flag and specify a the number of threads.
	} else {
		ss << "-p ";				//If testing in parallel, give the corresponding flag.
	}
	
	if(neighborlist) {
		ss << "-l 100 ";				//Add the neighborlist flag if applicable.
	}
	
	if(errorExpected) {
		ss << "-i 10000 > " << MCGPU << "bin/" << outputName << " 2>&1 ";	//If we expect an error, pipe cerr to a textfile where it can be read.
	} else {
		ss << "--name " << outputName << " -i 10000 ";							//If we do not expect an error, simply give the name for the results file.
	}
	std::string output = ss.str();
	std::cout << "RUNNING: " << output << std::endl;
	return output;
}

/**
 * Examines a result file to determine the final energy.
 * @param MCGPU The path to MCGPU's root.
 * @param resultsFile The name of the results file.
 * @return The final energy.
 */
double getEnergyResult(std::string MCGPU, std::string resultsFile) {
	std::ifstream infile(std::string(MCGPU + "bin/" + resultsFile).c_str());
	std::size_t found;
	
	for(std::string line; getline(infile, line);) {
		std::string str2 ("Final-Energy");
		std::string result;
		found = line.find(str2);
		if(found != std::string::npos) {
			result = line.substr(15);
			return strtod(result.c_str(), NULL);
		}
	}
	
	return -1;
}

/**
 * Examines a file that cerr was piped to and returns what function call the error occurred in.
 * @param MCGPU The path to MCGPU's root.
 * @param errorFile The name of the file cerr was piped to.
 * @return The function the error occurred in.
 */
std::string getErrorResult(std::string MCGPU, std::string errorFile) {
    std::ifstream infile(std::string(MCGPU + "bin/" + errorFile).c_str());
	std::size_t found;
	
    for(std::string line; getline(infile, line);) {
        std::string str2 ("Error");
        found = line.find(str2);
        if (found != std::string::npos) {
            return line.substr(7,13);
        }
    }
	return "ERROR: COULD NOT PARSE ERROR FILE!";
}

//Test t3pdimone with 1 primary index on CPU
TEST (t3pdimoneTest, OnePrimaryIndex)
{
    std::string MCGPU = getMCGPU_path();
	createConfigFile(MCGPU, "t3pdimone1MPI.config", "1");		
    system(buildCommand(MCGPU, "t3pdimone1MPI.config", "t3pdimone1MPI", true, false, false).c_str());
    double expected = -1810;
    double energyResult = getEnergyResult(MCGPU, "t3pdimone1MPI.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test t3pdimone with 1 primary index on GPU
TEST (t3pdimoneTest, PrimaryIndexGPU)
{
	std::string MCGPU = getMCGPU_path();
	system(buildCommand(MCGPU, "t3pdimone1MPI.config", "t3pdimone1MPI-GPU", false, false, false).c_str());
	double expected = -1810;
    double energyResult = getEnergyResult(MCGPU, "t3pdimone1MPI-GPU.results");
	EXPECT_NEAR(expected, energyResult, 100);
}

//Test t3pdimone with 1 primary index on CPU
//Using neighborlist
TEST (t3pdimoneTest, NeighborListFunction1MPI)
{
    std::string MCGPU = getMCGPU_path();
    system(buildCommand(MCGPU, "t3pdimone1MPI.config", "t3pdimone1MPI-NL", true, true, false).c_str());
    double expected = -1810;
    double energyResult = getEnergyResult(MCGPU, "t3pdimone1MPI-NL.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test t3pdimone with 1 primary index on GPU
//Using neighborlist
TEST (t3pdimoneTest, NeighborListFunction1MPI_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildCommand(MCGPU, "t3pdimone1MPI.config", "t3pdimone1MPI-NL-GPU", false, true, false).c_str());
    double expected = -1810;
    double energyResult = getEnergyResult(MCGPU, "t3pdimone1MPI-NL-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test t3pdimone with 1 primary [1,2] index on CPU
TEST (t3pdimoneTest, TwoPrimaryIndex)
{
    std::string MCGPU = getMCGPU_path();
	createConfigFile(MCGPU, "t3pdimone2MPI.config", "[1,2]");
    system(buildCommand(MCGPU, "t3pdimone2MPI.config", "t3pdimone2MPI", true, false, false).c_str());
    double expected = -1770;
    double energyResult = getEnergyResult(MCGPU, "t3pdimone2MPI.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test t3pdimone with 1 primary [1,2] index on GPU
TEST (t3pdimoneTest, TwoPrimaryIndexGPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildCommand(MCGPU, "t3pdimone2MPI.config", "t3pdimone2MPI-GPU", false, false, false).c_str());
    double expected = -1770;
    double energyResult = getEnergyResult(MCGPU, "t3pdimone2MPI-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test t3pdimone with multiple solvents 1,2 primary indexes on CPU
TEST (t3pdimoneTest, MultipleSolventDefinition)
{
	std::string MCGPU = getMCGPU_path();
	createConfigFile(MCGPU, "t3pdimoneMulSolvent.config", "1,2");
    system(buildCommand(MCGPU, "t3pdimoneMulSolvent.config", "t3pdimoneMulSolvent.txt", true, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "t3pdimoneMulSolvent.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}
       
//Test t3pdimone with multiple solvents 1,2 primary indexes on GPU
TEST (t3pdimoneTest, MultipleSolventDefinition_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildCommand(MCGPU, "t3pdimoneMulSolvent.config", "t3pdimoneMulSolvent-GPU.txt", false, false, true).c_str());
	std::string errorResult = getErrorResult(MCGPU, "t3pdimoneMulSolvent-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test t3pdimone with multiple solvents ([1,2],[3,4]) primary indexes on CPU
TEST (t3pdimoneTest, MultipleSolventDefinitionMPI)
{
    std::string MCGPU = getMCGPU_path();
	createConfigFile(MCGPU, "t3pdimoneMulSolventMPI.config", "[1,2],[3,4]");
    system(buildCommand(MCGPU, "t3pdimoneMulSolventMPI.config", "t3pdimoneMulSolventMPI.txt", true, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "t3pdimoneMulSolventMPI.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test t3pdimone with multiple solvents ([1,2],[3,4]) primary indexes on GPU
TEST (t3pdimoneTest, MultipleSolventDefinitionMPI_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildCommand(MCGPU, "t3pdimoneMulSolventMPI.config", "t3pdimoneMulSolventMPI-GPU.txt", false, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "t3pdimoneMulSolventMPI-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test t3pdimone with multiple solvents (1,[1,2]) primary indexes on CPU
TEST (t3pdimoneTest, SingleMultipleIndexes)
{
    std::string MCGPU = getMCGPU_path();
	createConfigFile(MCGPU, "t3pdimoneSingleMultipleIndexes.config", "1,[1,2]");
    system(buildCommand(MCGPU, "t3pdimoneSingleMultipleIndexes.config", "t3pdimoneSingleMultipleIndexes.txt", true, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "t3pdimoneSingleMultipleIndexes.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test t3pdimone with multiple solvents (1,[1,2]) primary indexes on GPU
TEST (t3pdimoneTest, SingleMultipleIndexes_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildCommand(MCGPU, "t3pdimoneSingleMultipleIndexes.config", "t3pdimoneSingleMultipleIndexes-GPU.txt", false, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "t3pdimoneSingleMultipleIndexes-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test t3pdimone with multiple solvents ([1,2],1) primary indexes on CPU
TEST (t3pdimoneTest, SingleMultipleIndexes2)
{
    std::string MCGPU = getMCGPU_path();
	createConfigFile(MCGPU, "t3pdimoneSingleMultipleIndexes2.config", "[1,2],1");
    system(buildCommand(MCGPU, "t3pdimoneSingleMultipleIndexes2.config", "t3pdimoneSingleMultipleIndexes2.txt", true, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "t3pdimoneSingleMultipleIndexes2.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test t3pdimone with multiple solvents ([1,2],1) primary indexes on GPU
TEST (t3pdimoneTest, SingleMultipleIndexes2_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildCommand(MCGPU, "t3pdimoneSingleMultipleIndexes2.config", "t3pdimoneSingleMultipleIndexes2-GPU.txt", false, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "t3pdimoneSingleMultipleIndexes2-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
} 
