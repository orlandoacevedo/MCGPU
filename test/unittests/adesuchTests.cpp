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
void createAdesuchConfigFile(std::string MCGPU, std::string fileName, std::string primaryAtomIndexString) {
	ConfigFileData settings = ConfigFileData(26.15, 26.15, 26.15, 298.15, .12, 100000, 256, "resources/bossFiles/oplsaa.par", 
	"test/unittests/MultipleSolvents/adesuchTests/adesuch.z", "test/unittests/MultipleSolvents/adesuchTests", 11.0, 
	12.0, 12345);
	createConfigFile(MCGPU, fileName, primaryAtomIndexString, settings);
}

//Test adesuch with 1 primary index on CPU
//Should result in an error since t3pdim has two molecules
TEST (adesuchTest, OnePrimaryIndex)
{
        std::string MCGPU = getMCGPU_path();
		createAdesuchConfigFile(MCGPU, "adesuch1MPI.config", "1");
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/adesuchTests/adesuch1MPI.config -s --verbose --threads 12 -i 10000 >"
                << " "
                << MCGPU << "/bin/adesuch1MPI.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult = getErrorResult(MCGPU, "adesuch1MPI.txt");
        EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}
