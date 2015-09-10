#include "Applications/Application.h"
#include "TestUtil.h"
#include "gtest/gtest.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>


/**
 * Generates a string that is used as a command to run the necessary test.
 * @param MCGPU The path to MCGPU's root.
 * @param configFile The name of the config file to use for the test.
 * @param outputName The name of the simulation, or the file to pipe cerr to (if expecting an error).
 * @param series true if the simulation is to be run in series, false otherwise.
 * @param neighborlist true if the simulation is to be run with a neighborlist, false otherwise.
 * @param errorExpected true if passing behavior for the test throws an error.
 */
std::string buildMeoh500Command(std::string MCGPU, std::string configFile, std::string outputName, bool series, bool neighborlist, bool errorExpected) {
	return buildCommand(MCGPU, configFile, outputName, series, neighborlist, errorExpected, "test/unittests/MultipleSolvents/Methanol500Tests");
}

/**
 * Generates a configuration file to be used in testing.
 * @param MCGPU_path The path to MCGPU's root directory.
 * @param fileName The name of the configuration file to generate.
 * @param primaryAtomIndexString The entry for the Primary Atom Index line of the config file.
 */
void createMeoh500ConfigFile(std::string MCGPU, std::string fileName, std::string primaryAtomIndexString) {

	ConfigFileData settings = ConfigFileData(32.91, 32.91, 32.91, 298.15, .12, 100000, 500, "resources/bossFiles/oplsaa.par",
	"test/unittests/MultipleSolvents/MethanolTests/meoh.z", "test/unittests/MultipleSolvents/Methanol500Tests", 11.0,
	12.0, 12345);

	createConfigFile(MCGPU, fileName, primaryAtomIndexString, settings);
}

//Test methanol with 500 molecules and a  primary index of 1 on CPU
TEST (Meoh500Test, OnePrimaryIndex)
{
	std::string MCGPU = getMCGPU_path();
	createMeoh500ConfigFile(MCGPU, "meoh5001MPI.config", "1");
    system(buildMeoh500Command(MCGPU, "meoh5001MPI.config", "meoh5001MPI", true, false, false).c_str());
    double expected = -3454;
    double energyResult = getEnergyResult(MCGPU, "meoh5001MPI.results");
    EXPECT_NEAR(expected, energyResult, 100);
}



//Test methanol with 500 molecules and a  primary index of 1 on GPU
TEST (Meoh500Test, OnePrimaryIndex_GPU)
{
	std::string MCGPU = getMCGPU_path();
    system(buildMeoh500Command(MCGPU, "meoh5001MPI.config", "meoh5001MPI-GPU", false, false, false).c_str());
    double expected = -3454;
    double energyResult = getEnergyResult(MCGPU, "meoh5001MPI-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test methanol with 500 molecules and a  primary index of 1 on CPU
//Using neighborlist
TEST (Meoh500Test, NeighborListFunction1MPI)
{
    std::string MCGPU = getMCGPU_path();
    system(buildMeoh500Command(MCGPU, "meoh5001MPI.config", "meoh5001MPI-NL", true, true, false).c_str());
	double expected = -3402;
    double energyResult = getEnergyResult(MCGPU, "meoh5001MPI-NL.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test methanol with 500 molecules and a  primary index of 1 on GPU
//Using neighborlist
TEST (Meoh500Test, NeighborListFunction1MPI_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildMeoh500Command(MCGPU, "meoh5001MPI.config", "meoh5001MPI-NL-GPU", false, true, false).c_str());
    double expected = -3402;
    double energyResult = getEnergyResult(MCGPU, "meoh5001MPI-NL-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test methanol with 500 molecules and two  primary indexes [1,2] on CPU
TEST (Meoh500Test, TwoPrimaryIndex)
{
	std::string MCGPU = getMCGPU_path();
	createMeoh500ConfigFile(MCGPU, "meoh5002MPI.config", "[1,2]");
    system(buildMeoh500Command(MCGPU, "meoh5002MPI.config", "meoh5002MPI", true, false, false).c_str());
    double expected = -3511;
    double energyResult = getEnergyResult(MCGPU, "meoh5002MPI.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test methanol with 500 molecules and two  primary indexes [1,2] on GPU
TEST (Meoh500Test, TwoPrimaryIndex_GPU)
{
	std::string MCGPU = getMCGPU_path();
    system(buildMeoh500Command(MCGPU, "meoh5002MPI.config", "meoh5002MPI-GPU", false, false, false).c_str());
    double expected = -3511;
    double energyResult = getEnergyResult(MCGPU, "meoh5002MPI-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test methanol with 500 molecules and primary index of [1,2] on CPU
//Using neighborlist
TEST (Meoh500Test, NeighborListFunction2MPI)
{
    std::string MCGPU = getMCGPU_path();
    system(buildMeoh500Command(MCGPU, "meoh5002MPI.config", "meoh5002MPI-NL", true, true, false).c_str());
    double expected = -3433;
    double energyResult = getEnergyResult(MCGPU, "meoh5002MPI-NL.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test methanol with 500 molecules and primary index of [1,2] on GPU
//Using neighborlist
TEST (Meoh500Test, NeighborListFunction2MPI_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildMeoh500Command(MCGPU, "meoh5002MPI.config", "meoh5002MPI-NL-GPU", false, true, false).c_str());
	double expected = -3433;
    double energyResult = getEnergyResult(MCGPU, "meoh5002MPI-NL-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test methanol with 500 molecules and multiple solvent primary indexes 1,2  on CPU
TEST (Meoh500Test, MultipleSolventDefinition)
{
    std::string MCGPU = getMCGPU_path();
	createMeoh500ConfigFile(MCGPU, "meoh500MulSolvent.config", "1,2");
    system(buildMeoh500Command(MCGPU, "meoh500MulSolvent.config", "meoh500MulSolvent.txt", true, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "meoh500MulSolvent.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test methanol with 500 molecules and multiple solvent primary indexes 1,2  on GPU
TEST (Meoh500Test, MultipleSolventDefinition_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildMeoh500Command(MCGPU, "meoh500MulSolvent.config", "meoh500MulSolvent-GPU.txt", false, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "meoh500MulSolvent-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test methanol with 500 molecules and multiple solvent primary indexes [1,2],[3,4]  on CPU
TEST (Meoh500Test, MultipleSolventDefinitionMPI)
{
    std::string MCGPU = getMCGPU_path();
	createMeoh500ConfigFile(MCGPU, "meoh500MulSolventMPI.config", "[1,2],[3,4]");
    system(buildMeoh500Command(MCGPU, "meoh500MulSolventMPI.config", "meoh500MulSolventMPI.txt", true, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "meoh500MulSolventMPI.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test methanol with 500 molecules and multiple solvent primary indexes [1,2],[3,4] on GPU
TEST (Meoh500Test, MultipleSolventDefinitionMPI_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildMeoh500Command(MCGPU, "meoh500MulSolventMPI.config", "meoh500MulSolventMPI-GPU.txt", false, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "meoh500MulSolventMPI-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test methanol with 500 molecules and multiple solvent primary indexes 1, [1,2] on CPU
TEST (Meoh500Test, SingleMultipleIndexes)
{
    std::string MCGPU = getMCGPU_path();
	createMeoh500ConfigFile(MCGPU, "meoh500SingleMultipleIndexes.config", "1, [1,2]");
    system(buildMeoh500Command(MCGPU, "meoh500SingleMultipleIndexes.config", "meoh500SingleMultipleIndexes.txt", true, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "meoh500SingleMultipleIndexes.txt");
	EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test methanol with 500 molecules and multiple solvent primary indexes 1, [1,2] on GPU
TEST (Meoh500Test, SingleMultipleIndexes_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildMeoh500Command(MCGPU, "meoh500SingleMultipleIndexes.config", "meoh500SingleMultipleIndexes-GPU.txt", false, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "meoh500SingleMultipleIndexes-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test methanol with 500 molecules and multiple solvent primary indexes [1,2],1 on CPU
TEST (Meoh500Test, SingleMultipleIndexes2)
{
    std::string MCGPU = getMCGPU_path();
	createMeoh500ConfigFile(MCGPU, "meoh500SingleMultipleIndexes2.config", "[1,2],1");
    system(buildMeoh500Command(MCGPU, "meoh500SingleMultipleIndexes2.config", "meoh500SingleMultipleIndexes2.txt", true, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "meoh500SingleMultipleIndexes2.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test methanol with 500 molecules and multiple solvent primary indexes [1,2],1 on GPU
TEST (Meoh500Test, SingleMultipleIndexes2_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildMeoh500Command(MCGPU, "meoh500SingleMultipleIndexes2.config", "meoh500SingleMultipleIndexes2-GPU.txt", false, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "meoh500SingleMultipleIndexes2-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}
