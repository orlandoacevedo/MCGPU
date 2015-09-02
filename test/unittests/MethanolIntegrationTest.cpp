
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
void createMethanolConfigFile(std::string MCGPU, std::string fileName, std::string primaryAtomIndexString) {
	ConfigFileData settings = ConfigFileData(26.15, 26.15, 26.15, 298.15, .12, 100000, 256, "resources/bossFiles/oplsaa.par", 
	"test/unittests/Integration/MethanolTest/meoh.z", "test/unittests/Integration/MethanolTest", 11.0, 
	12.0, 12345);
	createConfigFile(MCGPU, fileName, primaryAtomIndexString, settings);
}

// Descr: evident
// Implementation details: See gtest/samples for GTest syntax and usage
TEST(MethanolTest, FrontToEndIntegrationTest)
{
	std::string MCGPU = getMCGPU_path();
	createMethanolConfigFile(MCGPU, "MethanolTest.config", "1");
	
    std::stringstream ss;    
	ss << MCGPU << "/bin/metrosim "
       << " " // don't forget a space between the path and the arguments
       << MCGPU << "/test/unittests/Integration/MethanolTest/MethanolTest.config -s --name methanolCPU -k";
	
    // Launch MCGPU application in serial, expect output files in /MCGPU/ directory   
	system(ss.str().c_str());
	
	double expected = -1900;
	double energyResult = getEnergyResult(MCGPU, "methanolCPU.results");
	
	// Clean up
	remove(  std::string(MCGPU + "/test/unittests/Integration/MethanolTest/MethanolTest.config").c_str());
	remove(  std::string(MCGPU + "/bin/methanolCPU.pdb").c_str()  );
	remove(  std::string(MCGPU + "/bin/methanolCPU.results").c_str()  );
	remove(  std::string(MCGPU + "/bin/methanolCPU_100000.state").c_str()  );
	
	EXPECT_NEAR(expected, energyResult, 100 );
}


