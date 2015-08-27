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
 * Generates a configuration file to be used in testing.
 * @param MCGPU_path The path to MCGPU's root directory.
 * @param fileName The name of the configuration file to generate.
 * @param primaryAtomIndexString The entry for the Primary Atom Index line of the config file.
 */
void createT3pDimConfigFile(std::string MCGPU, std::string fileName, std::string primaryAtomIndexString) {
	ConfigFileData settings = ConfigFileData(26.15, 26.15, 26.15, 298.15, .12, 100000, 256, "resources/bossFiles/oplsaa.par", 
	"test/unittests/MultipleSolvents/t3pdimTests/t3pdim.z", "test/unittests/MultipleSolvents/t3pdimTests", 11.0, 
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
std::string buildT3pDimCommand(std::string MCGPU, std::string configFile, std::string outputName, bool series, bool neighborlist, bool errorExpected) {
	return buildCommand(MCGPU, configFile, outputName, series, neighborlist, errorExpected, "test/unittests/MultipleSolvents/t3pdimTests");
}

//Test t3pdim with 1 primary index on CPU
//Should result in an error since t3pdim has two molecules
TEST (t3pdimTest, OnePrimaryIndex) 
{
	std::string MCGPU = getMCGPU_path();
	createT3pDimConfigFile(MCGPU, "t3pdim1MPI.config", "1");
    system(buildT3pDimCommand(MCGPU, "t3pdim1MPI.config", "t3pdim1MPI.txt", true, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "t3pdim1MPI.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test t3pdim with 1 primary index on GPU
//Should result in an error since t3pdim has two molecules
TEST (t3pdimTest, OnePrimaryIndexGPU)
{
	std::string MCGPU = getMCGPU_path();
	system(buildT3pDimCommand(MCGPU, "t3pdim1MPI.config", "t3pdim1MPI-GPU.txt", false, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "t3pdim1MPI-GPU.txt");
	EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

 
//Test t3pdim with [1,2] primary index on CPU
//Should result in an error since t3pdim has two molecules
TEST (t3pdimTest, TwoPrimaryIndex)
{
	std::string MCGPU = getMCGPU_path();
	createT3pDimConfigFile(MCGPU, "t3pdim2MPI.config", "[1,2]");
    system(buildT3pDimCommand(MCGPU, "t3pdim2MPI.config", "t3pdim2MPI.txt", true, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "t3pdim2MPI.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test t3pdim with [1,2] primary index on GPU
//Should result in an error since t3pdim has two molecules
TEST (t3pdimTest, TwoPrimaryIndexGPU)
{
	std::string MCGPU = getMCGPU_path();
    system(buildT3pDimCommand(MCGPU, "t3pdim2MPI.config", "t3pdim2MPI-GPU.txt", false, false, true).c_str());
    std::string errorResult = getErrorResult(MCGPU, "t3pdim2MPI-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}


//Test t3pdim with multiple solvents with indexes 1,2 on CPU
TEST (t3pdimTest, MultipleSolventDefinition)
{
	std::string MCGPU = getMCGPU_path();
	createT3pDimConfigFile(MCGPU, "t3pdimMulSolvents.config", "1,2");
    system(buildT3pDimCommand(MCGPU, "t3pdimMulSolvents.config", "t3pdimMulSolvents", true, false, false).c_str());
    double expected = -1832;
    double energyResult = getEnergyResult(MCGPU, "t3pdimMulSolvents.results");
    EXPECT_NEAR(expected, energyResult, 100);
}


//Test t3pdim with multiple solvents with indexes 1,2 on GPU
TEST (t3pdimTest, MultipleSolventDefinitionGPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildT3pDimCommand(MCGPU, "t3pdimMulSolvents.config", "t3pdimMulSolvents-GPU", false, false, false).c_str());
    double expected = -1832;
    double energyResult = -getEnergyResult(MCGPU, "t3pdimMulSolvents-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test t3pdim with multiple solvents with indexes 1,2 on CPU
//Using neighborlist
TEST (t3pdimTest, NeighborListFunctionMulSolvents)
{
    std::string MCGPU = getMCGPU_path();
    system(buildT3pDimCommand(MCGPU, "t3pdimMulSolvents.config", "t3pdimMulSolvents-NL", true, true, false).c_str());
    double expected = -1727;
    double energyResult = -getEnergyResult(MCGPU, "t3pdimMulSolvents-NL.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test t3pdim with multiple solvents with indexes 1,2 on GPU
//Using neighborlist
TEST (t3pdimTest, NeighborListFunctionMulSolvents_GPU)
{
	std::string MCGPU = getMCGPU_path();
	system(buildT3pDimCommand(MCGPU, "t3pdimMulSolvents.config", "t3pdimMulSolvents-NL-GPU", false, true, false).c_str());
    double expected = -1727;
    double energyResult = getEnergyResult(MCGPU, "t3pdimMulSolvents-NL-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//test t3pdim for multiple solvents with [1,2],[3,4] primary indexes on CPU
TEST (t3pdimTest, MultipleSolventDefinitionMPI)
{
    std::string MCGPU = getMCGPU_path();
	createT3pDimConfigFile(MCGPU, "t3pdimMulSolventsMPI.config", "[1,2],[3,4]");
	system(buildT3pDimCommand(MCGPU, "t3pdimMulSolventsMPI.config", "t3pdimMulSolventsMPI", true, false, false).c_str());
    double expected = -1790;
    double energyResult = getEnergyResult(MCGPU, "t3pdimMulSolventsMPI.results");
    EXPECT_NEAR(expected, energyResult, 100);
}


//Test t3pdim with multiple solvents with indexes [1,2],[3,4] on GPU
TEST (t3pdimTest, MultipleSolventDefinitionMPI_GPU)
{
    std::string MCGPU = getMCGPU_path();
	system(buildT3pDimCommand(MCGPU, "t3pdimMulSolventsMPI.config", "t3pdimMulSolventsMPI-GPU", false, false, false).c_str());
    double expected = -1790;
    double energyResult = getEnergyResult(MCGPU, "t3pdimMulSolventsMPI-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//test t3pdim for multiple solvents with (1,[1,2]) primary indexes on CPU
TEST (t3pdimTest, SingleMultipleIndexes)
{
	std::string MCGPU = getMCGPU_path();
	createT3pDimConfigFile(MCGPU, "t3pdimSingleMultipleIndexes.config", "1,[1,2]");
	system(buildT3pDimCommand(MCGPU, "t3pdimSingleMultipleIndexes.config", "t3pdimSingleMultipleIndexes", true, false, false).c_str());
    double expected = -1770;
    double energyResult = getEnergyResult(MCGPU, "t3pdimSingleMultipleIndexes.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//test t3pdim for multiple solvents with (1,[1,2]) primary indexes on GPU
TEST (t3pdimTest, SingleMultipleIndexes_GPU)
{
    std::string MCGPU = getMCGPU_path();
	system(buildT3pDimCommand(MCGPU, "t3pdimSingleMultipleIndexes.config", "t3pdimSingleMultipleIndexes-GPU", false, false, false).c_str());
    double expected = -1770;
    double energyResult = getEnergyResult(MCGPU, "t3pdimSingleMultipleIndexes-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//test t3pdim for multiple solvents with ([1,2],1) primary indexes on CPU
TEST (t3pdimTest, SingleMultipleIndexes2)
{
	std::string MCGPU = getMCGPU_path();
	createT3pDimConfigFile(MCGPU, "t3pdimSingleMultipleIndexes2.config", "[1,2],1");
	system(buildT3pDimCommand(MCGPU, "t3pdimSingleMultipleIndexes2.config", "t3pdimSingleMultipleIndexes2", true, false, false).c_str());
    double expected = -1800;
    double energyResult = getEnergyResult(MCGPU, "t3pdimSingleMultipleIndexes2.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//test t3pdim for multiple solvents with ([1,2],1) primary indexes on GPU
TEST (t3pdimTest, SingleMultipleIndexes2_GPU)
{
	std::string MCGPU = getMCGPU_path();
	system(buildT3pDimCommand(MCGPU, "t3pdimSingleMultipleIndexes2.config", "t3pdimSingleMultipleIndexes-GPU", false, false, false).c_str());
    double expected = -1800;
    double energyResult = getEnergyResult(MCGPU, "t3pdimSingleMultipleIndexes2-GPU.results");
    EXPECT_NEAR(expected, energyResult,100 );
}
