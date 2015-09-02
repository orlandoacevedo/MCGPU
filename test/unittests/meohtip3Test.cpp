#include "Applications/Application.h"
#include "gtest/gtest.h"
#include "TestUtil.h"
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
void createMeohTip3ConfigFile(std::string MCGPU, std::string fileName, std::string primaryAtomIndexString) {
	ConfigFileData settings = ConfigFileData(23.7856, 23.7856, 23.7856, 298.15, .12, 100000, 256, "resources/bossFiles/oplsaa.par", 
	"test/unittests/MultipleSolvents/meohtip3Tests/meohtip3.z", "test/unittests/MultipleSolvents/meohtip3Tests", 11.0, 
	12.0, 12345);
	createConfigFile(MCGPU, fileName, primaryAtomIndexString, settings);
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
std::string buildMeohTip3Command(std::string MCGPU, std::string configFile, std::string outputName, bool series, bool neighborlist, bool errorExpected) {
	return buildCommand(MCGPU, configFile, outputName, series, neighborlist, errorExpected, "test/unittests/MultipleSolvents/meohtip3Tests");
}

//Test meohtip3 with 1 primary index on CPU
TEST (meohtip3Test, OnePrimaryIndex)
{
    std::string MCGPU = getMCGPU_path();
	createMeohTip3ConfigFile(MCGPU, "meohtip31MPI.config", "1");
	system(buildMeohTip3Command(MCGPU, "meohtip31MPI.config", "meohtip31MPI.txt", true, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "meohtip31MPI.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test meohtip3 with 1 primary index on GPU
TEST (meohtip3Test, OnePrimaryIndexGPU)
{
    std::string MCGPU = getMCGPU_path();
	system(buildMeohTip3Command(MCGPU, "meohtip31MPI.config", "meohtip31MPI-GPU.txt", false, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "meohtip31MPI-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test meohtip3 with 2 primary indexes [1,2] on CPU
TEST (meohtip3Test, TwoPrimaryIndex)
{
    std::string MCGPU = getMCGPU_path();
	createMeohTip3ConfigFile(MCGPU, "meohtip32MPI.config", "[1,2]");
	system(buildMeohTip3Command(MCGPU, "meohtip32MPI.config", "meohtip32MPI.txt", true, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "meohtip32MPI.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test meohtip3 with 2 primary indexes [1,2] on GPU
TEST (meohtip3Test, TwoPrimaryIndexGPU)
{
    std::string MCGPU = getMCGPU_path();
	system(buildMeohTip3Command(MCGPU, "meohtip32MPI.config", "meohtip32MPI-GPU.txt", false, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "meohtip32MPI-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//test meohtip3 with multiple solvents 1,2 primary indexes on CPU
TEST (meohtip3Test, MultpleSolventDefinition)
{
    std::string MCGPU = getMCGPU_path();
	createMeohTip3ConfigFile(MCGPU, "meohtip3MulSolvent.config", "1,2");
	system(buildMeohTip3Command(MCGPU, "meohtip3MulSolvent.config", "meohtip3MulSolvent", true, false, false).c_str());
    double expected = -1980;
	double energyResult = getEnergyResult(MCGPU, "meohtip3MulSolvent.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//test meohtip3 with multiple solvents 1,2 primary indexes on GPU
TEST (meohtip3Test, MultpleSolventDefinition_GPU)
 {
    std::string MCGPU = getMCGPU_path();
	system(buildMeohTip3Command(MCGPU, "meohtip3MulSolvent.config", "meohtip3MulSolvent-GPU", false, false, false).c_str());
    double expected = -1980;
	double energyResult = getEnergyResult(MCGPU, "meohtip3MulSolvent-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test meohtip3 with multiple solvents 1,2 primary indexes on CPU
//Using neighborlist
TEST (meohtip3Test, NeighborListFunctionMulSolvent) 
 {
    std::string MCGPU = getMCGPU_path();
	system(buildMeohTip3Command(MCGPU, "meohtip3MulSolvent.config", "meohtip3MulSolvent-NL", true, true, false).c_str());
    double expected = -2400;
	double energyResult = getEnergyResult(MCGPU, "meohtip3MulSolvent-NL.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test meohtip3 with multiple solvents 1,2 primary indexes on GPU
//Using neighborlist
TEST (meohtip3Test,NeighborListFunctionMulSolvent_GPU)
 {
    std::string MCGPU = getMCGPU_path();
	system(buildMeohTip3Command(MCGPU, "meohtip3MulSolvent.config", "meohtip3MulSolvent-NL-GPU", false, true, false).c_str());
    double expected = -2400;
	double energyResult = getEnergyResult(MCGPU, "meohtip3MulSolvent-NL-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);	 
}


//test meohtip3 with multiple solvents [1,2],[3,4] primary indexes on CPU
TEST (meohtip3Test, MultpleSolventDefinitionMPI)
 {
	std::string MCGPU = getMCGPU_path();
	createMeohTip3ConfigFile(MCGPU, "meohtip3MulSolventMPI.config", "[1,2],[3,4]");
	system(buildMeohTip3Command(MCGPU, "meohtip3MulSolventMPI.config", "meohtip3MulSolventMPI", true, false, false).c_str());
    double expected = -2040;
	double energyResult = getEnergyResult(MCGPU, "meohtip3MulSolventMPI.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//test meohtip3 with multiple solvents [1,2],[3,4] primary indexes on GPU
TEST (meohtip3Test, MultpleSolventDefinitionMPI_GPU)
 {
	std::string MCGPU = getMCGPU_path();
	system(buildMeohTip3Command(MCGPU, "meohtip3MulSolventMPI.config", "meohtip3MulSolventMPI-GPU", false, false, false).c_str());
    double expected = -2040;
	double energyResult = getEnergyResult(MCGPU, "meohtip3MulSolventMPI-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//test meohtip3 with multiple solvents (1,[1,2]) primary indexes on CPU
TEST (meohtip3Test, SingleMultipleIndexes)
 {
	std::string MCGPU = getMCGPU_path();
	createMeohTip3ConfigFile(MCGPU, "meohtip3SingleMultipleIndexes.config", "1,[1,2]");
	system(buildMeohTip3Command(MCGPU, "meohtip3SingleMultipleIndexes.config", "meohtip3SingleMultipleIndexes", true, false, false).c_str());
    double expected = -1990;
	double energyResult = getEnergyResult(MCGPU, "meohtip3SingleMultipleIndexes.results");
    EXPECT_NEAR(expected, energyResult, 100);	 
}

//test meohtip3 with multiple solvents (1,[1,2]) primary indexes on GPU
TEST (meohtip3Test, SingleMultipleIndexes_GPU)
 {
	std::string MCGPU = getMCGPU_path();
	system(buildMeohTip3Command(MCGPU, "meohtip3SingleMultipleIndexes.config", "meohtip3SingleMultipleIndexes-GPU", false, false, false).c_str());
    double expected = -1990;
	double energyResult = getEnergyResult(MCGPU, "meohtip3SingleMultipleIndexes-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);         
}

//test meohtip3 with multiple solvents ([1,2],1) primary indexes on CPU
TEST (meohtip3Test, SingleMultipleIndexes2)
 {
	std::string MCGPU = getMCGPU_path();
	createMeohTip3ConfigFile(MCGPU, "meohtip3SingleMultipleIndexes2.config", "[1,2],1");
	system(buildMeohTip3Command(MCGPU, "meohtip3SingleMultipleIndexes2.config", "meohtip3SingleMultipleIndexes2", true, false, false).c_str());
    double expected = -2010;
	double energyResult = getEnergyResult(MCGPU, "meohtip3SingleMultipleIndexes2.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//test meohtip3 with multiple solvents ([1,2], 1) primary indexes on GPU
TEST (meohtip3Test, SingleMultipleIndexes2_GPU)
 {
	std::string MCGPU = getMCGPU_path();
	system(buildMeohTip3Command(MCGPU, "meohtip3SingleMultipleIndexes2.config", "meohtip3SingleMultipleIndexes2-GPU", false, false, false).c_str());
    double expected = -2010;
	double energyResult = getEnergyResult(MCGPU, "meohtip3SingleMultipleIndexes2-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}
