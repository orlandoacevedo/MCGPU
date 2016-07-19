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
 * Generates a string that is used as a command to run the necessary test.
 * @param MCGPU The path to MCGPU's root.
 * @param configFile The name of the config file to use for the test.
 * @param outputName The name of the simulation, or the file to pipe cerr to
 *   (if expecting an error).
 * @param series true if the simulation is to be run in series, false otherwise
 * @param neighborlist true if the simulation is to be run with a neighborlist,
 *   false otherwise.
 * @param errorExpected true if passing behavior for the test throws an error.
 */
std::string buildT3pDimoneCommand(std::string MCGPU, std::string configFile,
                                  std::string outputName, bool series,
                                  bool neighborlist, bool errorExpected) {
  return buildCommand(MCGPU, configFile, outputName, series, neighborlist,
                      errorExpected,
                      "test/unittests/MultipleSolvents/t3pdimoneTests");
}
/**
 * Generates a configuration file to be used in testing.
 * @param MCGPU_path The path to MCGPU's root directory.
 * @param fileName The name of the configuration file to generate.
 * @param primaryAtomIndexString The entry for the Primary Atom Index line of the config file.
 */
void createT3pDimoneConfigFile(std::string MCGPU, std::string fileName,
                               std::string primaryAtomIndexString) {
  ConfigFileData settings = ConfigFileData(
      26.15, 26.15, 26.15, 298.15, .12, 100000, 256,
      "resources/bossFiles/oplsaa.par",
      "test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimone.z",
      "test/unittests/MultipleSolvents/t3pdimoneTests",
      11.0, 12.0, 12345);
  createConfigFile(MCGPU, fileName, primaryAtomIndexString, settings);
}

// Test t3pdimone with 1 primary index on CPU
TEST(t3pdimoneTest, OnePrimaryIndex)
{
  std::string MCGPU = getMCGPU_path();
  createT3pDimoneConfigFile(MCGPU, "t3pdimone1MPI.config", "1");
  system(buildT3pDimoneCommand(MCGPU, "t3pdimone1MPI.config", "t3pdimone1MPI",
                               true, false, false).c_str());
  double expected = -1810;
  double energyResult = getEnergyResult("t3pdimone1MPI.results");
  EXPECT_NEAR(expected, energyResult, 100);
}

//Test t3pdimone with 1 primary index on GPU
TEST(t3pdimoneTest, PrimaryIndexGPU)
{
  std::string MCGPU = getMCGPU_path();
  system(buildT3pDimoneCommand(MCGPU, "t3pdimone1MPI.config",
                               "t3pdimone1MPI-GPU", false, false,
                               false).c_str());
  double expected = -1810;
  double energyResult = getEnergyResult("t3pdimone1MPI-GPU.results");
  EXPECT_NEAR(expected, energyResult, 100);
}

// Test t3pdimone with 1 primary index on CPU
// Using neighborlist
TEST(t3pdimoneTest, NeighborListFunction1MPI)
{
  std::string MCGPU = getMCGPU_path();
  system(buildT3pDimoneCommand(MCGPU, "t3pdimone1MPI.config",
                               "t3pdimone1MPI-NL", true, true, false).c_str());
  double expected = -1810;
  double energyResult = getEnergyResult("t3pdimone1MPI-NL.results");
  EXPECT_NEAR(expected, energyResult, 100);
}

// Test t3pdimone with 1 primary index on GPU
// Using neighborlist
#ifdef _OPENACC
TEST(t3pdimoneTest, NeighborListFunction1MPI_GPU)
{
  std::string MCGPU = getMCGPU_path();
  system(buildT3pDimoneCommand(MCGPU, "t3pdimone1MPI.config",
                               "t3pdimone1MPI-NL-GPU", false, true,
                               false).c_str());
  double expected = -1810;
  double energyResult = getEnergyResult("t3pdimone1MPI-NL-GPU.results");
  EXPECT_NEAR(expected, energyResult, 100);
}
#endif

// Test t3pdimone with 1 primary [1,2] index on CPU
TEST(t3pdimoneTest, TwoPrimaryIndex)
{
  std::string MCGPU = getMCGPU_path();
  createT3pDimoneConfigFile(MCGPU, "t3pdimone2MPI.config", "[1,2]");
  system(buildT3pDimoneCommand(MCGPU, "t3pdimone2MPI.config", "t3pdimone2MPI",
                               true, false, false).c_str());
  double expected = -1770;
  double energyResult = getEnergyResult("t3pdimone2MPI.results");
  EXPECT_NEAR(expected, energyResult, 100);
}

// Test t3pdimone with 1 primary [1,2] index on GPU
TEST(t3pdimoneTest, TwoPrimaryIndexGPU)
{
  std::string MCGPU = getMCGPU_path();
  system(buildT3pDimoneCommand(MCGPU, "t3pdimone2MPI.config",
                               "t3pdimone2MPI-GPU", false, false,
                               false).c_str());
  double expected = -1770;
  double energyResult = getEnergyResult("t3pdimone2MPI-GPU.results");
  EXPECT_NEAR(expected, energyResult, 100);
}

// Test t3pdimone with multiple solvents 1,2 primary indexes on CPU
TEST(t3pdimoneTest, MultipleSolventDefinition)
{
  std::string MCGPU = getMCGPU_path();
  createT3pDimoneConfigFile(MCGPU, "t3pdimoneMulSolvent.config", "1,2");
  system(buildT3pDimoneCommand(MCGPU, "t3pdimoneMulSolvent.config",
                               "t3pdimoneMulSolvent.txt", true, false,
                               true).c_str());
  std::string errorResult = getErrorResult(MCGPU, "t3pdimoneMulSolvent.txt");
  EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

// Test t3pdimone with multiple solvents 1,2 primary indexes on GPU
#ifdef _OPENACC
TEST(t3pdimoneTest, MultipleSolventDefinition_GPU)
{
  std::string MCGPU = getMCGPU_path();
  system(buildT3pDimoneCommand(MCGPU, "t3pdimoneMulSolvent.config",
                               "t3pdimoneMulSolvent-GPU.txt", false, false,
                               true).c_str());
  std::string errorResult = getErrorResult(MCGPU,
                                           "t3pdimoneMulSolvent-GPU.txt");
  EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}
#endif

// Test t3pdimone with multiple solvents ([1,2],[3,4]) primary indexes on CPU
TEST(t3pdimoneTest, MultipleSolventDefinitionMPI)
{
  std::string MCGPU = getMCGPU_path();
  createT3pDimoneConfigFile(MCGPU, "t3pdimoneMulSolventMPI.config",
                            "[1,2],[3,4]");
  system(buildT3pDimoneCommand(MCGPU, "t3pdimoneMulSolventMPI.config",
                               "t3pdimoneMulSolventMPI.txt", true, false,
                               true).c_str());
  std::string errorResult = getErrorResult(MCGPU,
                                           "t3pdimoneMulSolventMPI.txt");
  EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

// Test t3pdimone with multiple solvents ([1,2],[3,4]) primary indexes on GPU
#ifdef _OPENACC
TEST(t3pdimoneTest, MultipleSolventDefinitionMPI_GPU)
{
  std::string MCGPU = getMCGPU_path();
  system(buildT3pDimoneCommand(MCGPU, "t3pdimoneMulSolventMPI.config",
                               "t3pdimoneMulSolventMPI-GPU.txt", false, false,
                               true).c_str());
  std::string errorResult = getErrorResult(MCGPU,
                                           "t3pdimoneMulSolventMPI-GPU.txt");
  EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}
#endif

// Test t3pdimone with multiple solvents (1,[1,2]) primary indexes on CPU
TEST(t3pdimoneTest, SingleMultipleIndexes)
{
  std::string MCGPU = getMCGPU_path();
  createT3pDimoneConfigFile(MCGPU, "t3pdimoneSingleMultipleIndexes.config",
                            "1,[1,2]");
  system(buildT3pDimoneCommand(MCGPU, "t3pdimoneSingleMultipleIndexes.config",
                               "t3pdimoneSingleMultipleIndexes.txt", true,
                               false, true).c_str());
  std::string errorResult = getErrorResult(
      MCGPU, "t3pdimoneSingleMultipleIndexes.txt");
  EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

// Test t3pdimone with multiple solvents (1,[1,2]) primary indexes on GPU
#ifdef _OPENACC
TEST(t3pdimoneTest, SingleMultipleIndexes_GPU)
{
  std::string MCGPU = getMCGPU_path();
  system(buildT3pDimoneCommand(MCGPU, "t3pdimoneSingleMultipleIndexes.config",
                               "t3pdimoneSingleMultipleIndexes-GPU.txt", false,
                               false, true).c_str());
  std::string errorResult = getErrorResult(
      MCGPU, "t3pdimoneSingleMultipleIndexes-GPU.txt");
  EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}
#endif

// Test t3pdimone with multiple solvents ([1,2],1) primary indexes on CPU
TEST(t3pdimoneTest, SingleMultipleIndexes2)
{
  std::string MCGPU = getMCGPU_path();
  createT3pDimoneConfigFile(MCGPU, "t3pdimoneSingleMultipleIndexes2.config",
                            "[1,2],1");
  system(buildT3pDimoneCommand(MCGPU, "t3pdimoneSingleMultipleIndexes2.config",
                               "t3pdimoneSingleMultipleIndexes2.txt", true,
                               false, true).c_str());
  std::string errorResult = getErrorResult(
      MCGPU, "t3pdimoneSingleMultipleIndexes2.txt");
  EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test t3pdimone with multiple solvents ([1,2],1) primary indexes on GPU
#ifdef _OPENACC
TEST(t3pdimoneTest, SingleMultipleIndexes2_GPU)
{
  std::string MCGPU = getMCGPU_path();
  system(buildT3pDimoneCommand(MCGPU, "t3pdimoneSingleMultipleIndexes2.config",
                               "t3pdimoneSingleMultipleIndexes2-GPU.txt",
                               false, false, true).c_str());
  std::string errorResult = getErrorResult(
      MCGPU, "t3pdimoneSingleMultipleIndexes2-GPU.txt");
  EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}
#endif
