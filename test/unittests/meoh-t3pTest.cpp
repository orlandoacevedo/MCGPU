#include "Applications/Application.h"
#include "gtest/gtest.h"
#include "TestUtil.h"

void createMeohT3PConfigFile(std::string MCGPU, std::string fileName,
                            double boxSize, int molecules) {
  ConfigFileData settings = ConfigFileData(
      boxSize, boxSize, boxSize, 298.15, 0.12, 100000, molecules,
      "resources/bossFiles/oplsaa.par",
      "test/unittests/MultipleSolvents/meoh-t3pTests/meoh-t3p.z",
      "test/unittests/MultipleSolvents/meoh-t3pTests", 11.0, 12.0, 12345);
  createConfigFile(MCGPU, fileName, "[1],[1]", settings);
}

std::string buildMeohT3PCommand(std::string MCGPU, std::string configFile,
                               std::string outputName, bool series,
                               bool neighborlist, bool errorExpected) {
  return buildCommand(MCGPU, configFile, outputName, series, neighborlist,
                      errorExpected,
                      "test/unittests/MultipleSolvents/meoh-t3pTests");
}

// Test meoh-t3p with 256 molecules
TEST(meoht3pTest, mol256) {
  std::string MCGPU = getMCGPU_path();
  std::string testName = "meoh-t3p256";
  createMeohT3PConfigFile(MCGPU, testName + ".config", 23.7856, 256);
  system(buildMeohT3PCommand(MCGPU, testName + ".config", testName,
                             true, false, false).c_str());
  double expected = -2007.0;
  double energyResult = getEnergyResult(MCGPU, testName + ".results");
  EXPECT_NEAR(expected, energyResult, 100);
}

// Test meoh-t3p with 256 molecules on the GPU
TEST(meoht3pTest, mol256_GPU) {
  std::string MCGPU = getMCGPU_path();
  std::string testName = "meoh-t3p256";
  createMeohT3PConfigFile(MCGPU, testName + ".config", 23.7856, 256);
  system(buildMeohT3PCommand(MCGPU, testName + ".config", testName,
                             false, false, false).c_str());
  double expected = -2007.0;
  double energyResult = getEnergyResult(MCGPU, testName + ".results");
  EXPECT_NEAR(expected, energyResult, 100);
}

// Test meoh-t3p with 500 molecules
TEST(meoht3pTest, mol500) {
  std::string MCGPU = getMCGPU_path();
  std::string testName = "meoh-t3p500";
  createMeohT3PConfigFile(MCGPU, testName + ".config", 32.91, 500);
  system(buildMeohT3PCommand(MCGPU, testName + ".config", testName,
                             true, false, false).c_str());
  double expected = -3148.0;
  double energyResult = getEnergyResult(MCGPU, testName + ".results");
  EXPECT_NEAR(expected, energyResult, 100);
}

// Test meoh-t3p with 500 molecules on the GPU
TEST(meoht3pTest, mol500_GPU) {
  std::string MCGPU = getMCGPU_path();
  std::string testName = "meoh-t3p500";
  createMeohT3PConfigFile(MCGPU, testName + ".config", 32.91, 500);
  system(buildMeohT3PCommand(MCGPU, testName + ".config", testName,
                             false, false, false).c_str());
  double expected = -3148.0;
  double energyResult = getEnergyResult(MCGPU, testName + ".results");
  EXPECT_NEAR(expected, energyResult, 100);
}

// Test meoh-t3p with 8788 molecules
TEST(meoht3pTest, mol8788) {
  std::string MCGPU = getMCGPU_path();
  std::string testName = "meoh-t3p8788";
  createMeohT3PConfigFile(MCGPU, testName + ".config", 85.57, 8788);
  system(buildMeohT3PCommand(MCGPU, testName + ".config", testName,
                             true, false, false).c_str());
  double expected = -13219.0;
  double energyResult = getEnergyResult(MCGPU, testName + ".results");
  EXPECT_NEAR(expected, energyResult, 100);
}

// Test meoh-t3p with 8788 molecules on the GPU
TEST(meoht3pTest, mol8788_GPU) {
  std::string MCGPU = getMCGPU_path();
  std::string testName = "meoh-t3p8788";
  createMeohT3PConfigFile(MCGPU, testName + ".config", 85.57, 8788);
  system(buildMeohT3PCommand(MCGPU, testName + ".config", testName,
                             false, false, false).c_str());
  double expected = -13219.0;
  double energyResult = getEnergyResult(MCGPU, testName + ".results");
  EXPECT_NEAR(expected, energyResult, 100);
}
