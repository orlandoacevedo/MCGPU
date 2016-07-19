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
void createWaterConfigFile(std::string MCGPU, std::string fileName, std::string primaryAtomIndexString) {
  ConfigFileData settings = ConfigFileData(
      55.0, 55.0, 55.0, 298.15, .15, 1000, 5120,
      "resources/bossFiles/oplsaa.par",
      "test/unittests/Integration/WaterTest/watt4p.z",
      "test/unittests/Integration/WaterTest", 25.0, 15.0, 12345);
  createConfigFile(MCGPU, fileName, primaryAtomIndexString, settings);
}

double waterFrontToEndIntegrationTest(bool parallel) {
  std::string MCGPU = getMCGPU_path();
  createWaterConfigFile(MCGPU, "WaterTest.config", "1");

  std::stringstream ss;
  ss << MCGPU << "/bin/metrosim "
     << MCGPU << "/test/unittests/Integration/WaterTest/WaterTest.config "
     << (parallel ? "-p" : "-s") << " --name waterCPU > /dev/null";

  // Launch MCGPU application, expect output files in current directory
  system(ss.str().c_str());
  double energyResult = getEnergyResult("waterCPU.results");

  // Delete the files used for the test
  remove(std::string(MCGPU + "test/unittests/Integration/WaterTest.config").c_str());
  remove(std::string(MCGPU + "/waterCPU.pdb").c_str());
  remove(std::string(MCGPU + "/waterCPU.results").c_str());
  remove(std::string(MCGPU + "/waterCPU_1000.state").c_str());

  return energyResult;
}

// Descr: evident
// Implementation details: See gtest/samples for GTest syntax and usage
TEST(WaterTest, FrontToEndIntegrationTest)
{
  double result = waterFrontToEndIntegrationTest(false);
  double expected = -14000, margin = 250;
  EXPECT_NEAR(expected, result, margin);
}

// Same water test, but will run on the GPU
TEST(WaterTest, FrontToEndIntegrationTest_GPU)
{
  double result = waterFrontToEndIntegrationTest(true);
  double expected = -14000, margin = 250;
  EXPECT_NEAR(expected, result, margin);
}

