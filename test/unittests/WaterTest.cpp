
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
	ConfigFileData settings = ConfigFileData(55.0, 55.0, 55.0, 298.15, .15, 1000, 5120, "resources/bossFiles/oplsaa.par", 
	"test/unittests/Integration/WaterTest/watt4p.z", "test/unittests/Integration/WaterTest", 25.0, 
	15.0, 12345);
	createConfigFile(MCGPU, fileName, primaryAtomIndexString, settings);
}

// Descr: evident
// Implementation details: See gtest/samples for GTest syntax and usage
TEST(WaterTest, FrontToEndIntegrationTest)
{
	std::string MCGPU = getMCGPU_path();
	createWaterConfigFile(MCGPU, "WaterTest.config", "1");
	
    std::stringstream ss;    
	ss << MCGPU << "/bin/metrosim "
       << " " // don't forget a space between the path and the arguments
       << MCGPU << "/test/unittests/Integration/WaterTest/WaterTest.config -s --name waterCPU -k";
	
    // Launch MCGPU application in serial, expect output files in /MCGPU/ directory   
	system(ss.str().c_str());
	
	double expected = -14000;
	double energyResult = getEnergyResult(MCGPU, "waterCPU.results");
	
	remove(  std::string(MCGPU + "test/unittests/Integration/WaterTest.config").c_str());
	remove(  std::string(MCGPU + "/waterCPU.pdb").c_str()  );
	remove(  std::string(MCGPU + "/waterCPU.results").c_str()  );
	remove(  std::string(MCGPU + "/waterCPU_1000.state").c_str()  );
	
	EXPECT_NEAR(expected, energyResult, 250 );
}


