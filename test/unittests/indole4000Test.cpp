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
void createIndole4000ConfigFile(std::string MCGPU, std::string fileName, std::string primaryAtomIndexString) {
	ConfigFileData settings = ConfigFileData(87.17, 87.17, 87.17, 298.15, .06, 100000, 4000, "resources/bossFiles/oplsaa.par", 
	"test/unittests/MultipleSolvents/indole4000Tests/indole.z", "test/unittests/MultipleSolvents/indole4000Tests", 12.0, 
	6.0, 12345);
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
std::string buildIndole4000Command(std::string MCGPU, std::string configFile, std::string outputName, bool series, bool neighborlist, bool errorExpected) {
	return buildCommand(MCGPU, configFile, outputName, series, neighborlist, errorExpected, "test/unittests/MultipleSolvents/indole4000Tests");
}


//Test indole with 4000 molecules and a  primary index of 1 on CPU
TEST (Indole4000Test, OnePrimaryIndex)
{
	std::string MCGPU = getMCGPU_path();
	createIndole4000ConfigFile(MCGPU, "indole4000-1MPI.config", "1");
    system(buildIndole4000Command(MCGPU, "indole4000-1MPI.config", "indole4000-1MPI", true, false, false));
    double expected = 588107;
    double energyResult = getEnergyResult(MCGPU, "indole4000-1MPI.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test indole with 4000 molecules and a  primary index of 1 on GPU
TEST (Indole4000Test, OnePrimaryIndex_GPU)
{
    std::string MCGPU = getMCGPU_path();
	system(buildIndole4000Command(MCGPU, "indole4000-1MPI.config", "indole4000-1MPI-GPU", false, false, false));
    double expected = 588107;
    double energyResult = getEnergyResult(MCGPU, "indole4000-1MPI-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}


//Test indole with 4000 molecules and a  primary index of 1 on CPU
//Uses neighbor list function
TEST (Indole4000Test, NeighborListFunction1MPI)
{
    std::string MCGPU = getMCGPU_path();
	system(buildIndole4000Command(MCGPU, "indole4000-1MPI.config", "indole4000-1MPI-NL", true, true, false));
    double expected = 530000;
    double energyResult = getEnergyResult(MCGPU, "indole4000-1MPI-NL.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test indole with 4000 molecules and a  primary index of 1 on GPU
//Uses neighbor list function
TEST (Indole4000Test, NeighborListFunction1MPI_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildIndole4000Command(MCGPU, "indole4000-1MPI.config", "indole4000-1MPI-NL-GPU", false, true, false));
    double expected = 530000;
    double energyResult = getEnergyResult(MCGPU, "indole4000-1MPI-NL-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test indole with 4000 molecules and a primary index of [1,2] on CPU
TEST (Indole4000Test, TwoPrimaryIndex)
{
    std::string MCGPU = getMCGPU_path();
    createIndole4000ConfigFile(MCGPU, "indole4000-2MPI.config", "[1,2]");
    system(buildIndole4000Command(MCGPU, "indole4000-2MPI.config", "indole4000-2MPI", true, false, false));
    double expected = 557412;
    double energyResult = getEnergyResult(MCGPU, "indole4000-2MPI.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test indole with 4000 molecules and a primary index of [1,2] on GPU
TEST (Indole4000Test, TwoPrimaryIndex_GPU)
{
    std::string MCGPU = getMCGPU_path();
	system(buildIndole4000Command(MCGPU, "indole4000-2MPI.config", "indole4000-2MPI-GPU", false, false, false);
    double expected = 557412;
    double energyResult = getEnergyResult(MCGPU, "indole4000-2MPI-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}


//Test indole with 4000 molecules and a primary index of [1,2] on CPU
//Uses neighborlist function
TEST (Indole4000Test, NeighborListFunction2MPI)
{
    std::string MCGPU = getMCGPU_path();
    system(buildIndole4000Command(MCGPU, "indole4000-2MPI.config", "indole4000-2MPI-NL", true, true, false));
    double expected = 547878;
    double energyResult = getEnergyResult(MCGPU, "indole4000-2MPI-NL.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test indole with 4000 molecules and a primary index of [1,2] on GPU
//Uses neighborlist function
TEST (Indole4000Test, NeighborListFunction2MPI_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildIndole4000Command(MCGPU, "indole4000-2MPI.config", "indole4000-2MPI-NL-GPU", false, true, false));
    double expected = 547878;
    double energyResult = getEnergyResult(MCGPU, "indole4000-2MPI-NL-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test indole with 4000 molecules and multiple solvent primary indexes of 1,2 on CPU
TEST (Indole4000Test, MultipleSolventsDefinition)
{
    std::string MCGPU = getMCGPU_path();
	createIndole4000ConfigFile(MCGPU, "indole4000MulSolvent.config", "1,2");
    system(buildIndole4000Command(MCGPU, "indole4000MulSolvent.config", "indole4000MulSolvent.txt", true, false, true));
	std::string errorResult = getErrorResult(MCGPU, "indole4000MulSolvent.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}
        

//Test indole with 4000 molecules and multiple solvent primary indexes of 1,2 on GPU
TEST (Indole4000Test, MultipleSolventsDefinition_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildIndole4000Command(MCGPU, "indole4000MulSolvent.config", "indole4000MulSolvent-GPU.txt", false, false, true));
	std::string errorResult = getErrorResult(MCGPU, "indole4000MulSolvent-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test indole with 4000 molecules and multiple solvent primary indexes of [1,2],[3,4] on CPU
TEST (Indole4000Test, MultipleSolventsDefinitionMPI)
{
    std::string MCGPU = getMCGPU_path();
	createIndole4000ConfigFile(MCGPU, "indole4000MulSolvent-MPI.config", "[1,2],[3,4]");
    system(buildIndole4000Command(MCGPU, "indole4000MulSolvent-MPI.config", "indole4000MulSolvent-MPI.txt", true, false, true));
	std::string errorResult = getErrorResult(MCGPU, "indole4000MulSolvent.txt-MPI");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test indole with 4000 molecules and multiple solvent primary indexes of [1,2],[3,4] on GPU
TEST (Indole4000Test, MultipleSolventsDefinitionMPI_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildIndole4000Command(MCGPU, "indole4000MulSolvent-MPI.config", "indole4000MulSolvent-MPI-GPU.txt", false, false, true));
	std::string errorResult = getErrorResult(MCGPU, "indole4000MulSolvent-MPI-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test indole with 4000 molecules and multiple solvent primary indexes of 1,[1,2] on CPU
TEST (Indole4000Test, SingleMultipleIndexes)
{
    std::string MCGPU = getMCGPU_path();
	createIndole4000ConfigFile(MCGPU, "indole4000SingleMultipleIndexes.config", "1,[1,2]");
    system(buildIndole4000Command(MCGPU, "indole4000SingleMultipleIndexes.config", "indole4000SingleMultipleIndexes.txt", true, false, true));
	std::string errorResult = getErrorResult(MCGPU, "indole4000SingleMultipleIndexes.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test indole with 4000 molecules and multiple solvent primary indexes of 1,[1,2] on GPU
TEST (Indole4000Test, SingleMultipleIndexes_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildIndole4000Command(MCGPU, "indole4000SingleMultipleIndexes.config", "indole4000SingleMultipleIndexes-GPU.txt", false, false, true));
	std::string errorResult = getErrorResult(MCGPU, "indole4000SingleMultipleIndexes-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test indole with 4000 molecules and multiple solvent primary indexes of [1,2],1 on CPU
TEST (Indole4000Test, SingleMultipleIndexes2)
{
    std::string MCGPU = getMCGPU_path();
	createIndole4000ConfigFile(MCGPU, "indole4000SingleMultipleIndexes.config", "[1,2],1");
    system(buildIndole4000Command(MCGPU, "indole4000SingleMultipleIndexes2.config", "indole4000SingleMultipleIndexes2.txt", true, false, true));
	std::string errorResult = getErrorResult(MCGPU, "indole4000SingleMultipleIndexes2.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test indole with 4000 molecules and multiple solvent primary indexes of [1,2],1 on GPU
TEST (Indole4000Test, SingleMultipleIndexes2_GPU)
{
    std::string MCGPU = getMCGPU_path();
    system(buildIndole4000Command(MCGPU, "indole4000SingleMultipleIndexes2.config", "indole4000SingleMultipleIndexes2-GPU.txt", true, false, true));
	std::string errorResult = getErrorResult(MCGPU, "indole4000SingleMultipleIndexes2-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}
