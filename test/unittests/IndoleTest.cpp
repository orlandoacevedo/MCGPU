
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
void createIndoleConfigFile(std::string MCGPU, std::string fileName, std::string primaryAtomIndexString) {
	ConfigFileData settings = ConfigFileData(35.36, 35.36, 35.36, 298.15, .12, 40000, 267, "resources/bossFiles/oplsaa.par", 
	"test/unittests/Integration/IndoleTest/indole.z", "test/unittests/Integration/IndoleTest", 12.0, 
	12.0, 12345);
	createConfigFile(MCGPU, fileName, primaryAtomIndexString, settings);
}

// Descr: evident
// Implementation details: See gtest/samples for GTest syntax and usage
TEST(IndoleTest, FrontToEndIntegrationTest)
{
	std::string MCGPU = getMCGPU_path();
	createIndoleConfigFile(MCGPU, "IndoleTest.config", "1");
    std::stringstream ss;    
	ss << MCGPU << "/bin/metrosim "
       << " " // don't forget a space between the path and the arguments
       << MCGPU << "/test/unittests/Integration/IndoleTest/IndoleTest.config -s --name indoleCPU -k";
	
	
    // Launch MCGPU application in serial, expect output files in /MCGPU/ directory   
	system(ss.str().c_str());
	
	double expected = 59800;
	double energyResult = getEnergyResult(MCGPU, "indoleCPU.results");
	// Clean up
	remove(  std::string(MCGPU + "test/unittests/Integration/IndoleTest/IndoleTest.config").c_str());
	remove(  std::string(MCGPU + "/indoleCPU.pdb").c_str()  );
	remove(  std::string(MCGPU + "/indoleCPU.results").c_str()  );
	remove(  std::string(MCGPU + "/indoleCPU_40000.state").c_str()  );
	
	EXPECT_NEAR(expected, energyResult, 100 );
}


