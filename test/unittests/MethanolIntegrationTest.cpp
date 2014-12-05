
#include "Applications/Application.h"
#include "gtest/gtest.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

// Descr: evident
// Implementation details: See gtest/samples for GTest syntax and usage
TEST(MethanolTest, FrontToEndIntegrationTest)
{
	string directory = get_current_dir_name();
	
	std::string mc ("MCGPU");
	std::size_t found = directory.find(mc);
	
	if (found != std::string::npos) {
		directory = directory.substr(0,found+6);
	}
	std::string MCGPU = directory;
	
	
	// Create the test config file,
	// 		because it requires full hardcoded filepath names that will change on each user's computer
	
	ofstream configFile;
	std::string configFilePath (std::string( MCGPU + "/test/unittests/Integration/MethanolTest/MethanolTest.config"));
	configFile.open( configFilePath.c_str() );
	std::stringstream cfg;
	cfg << ""
		<< "#size of periodic box (x, y, z in angstroms)\n"
		<< "26.15\n"
		<< "26.15\n"
		<< "26.15\n"
		<< "#temperature in Kelvin\n"
		<< "298.15\n"
		<< "#max translation\n"
		<< ".12\n"
		<< "#number of steps\n"
		<< "100000\n"
		<< "#number of molecues\n"
		<< "256\n"
		<< "#path to opls.par file\n"
		<< MCGPU << "/resources/bossFiles/oplsaa.par\n"
		<< "#path to z matrix file\n"
		<< MCGPU << "/test/unittests/Integration/MethanolTest/meoh.z\n"
		<< "#path to state input\n"
		<< "\n"
		<< "#path to state output\n"
		<< MCGPU << "/test/unittests/Integration/MethanolTest\n"
		<< "#pdb output path\n"
		<< "methanoltest.pdb\n"
		<< "#cutoff distance in angstroms\n"
		<< "11.0\n"
		<< "#max rotation\n"
		<< "12.0\n"
		<< "#Random Seed Input\n"
		<< "12345\n"
		<< "#Primary Atom Index\n"
		<< "1";
		
	configFile << cfg.str();
	configFile.close();
	
    std::stringstream ss;    
	ss << MCGPU << "/bin/metrosim "
       << " " // don't forget a space between the path and the arguments
       << MCGPU << "/test/unittests/Integration/MethanolTest/MethanolTest.config -s --name methanolCPU -k";
	
	
    // Launch MCGPU application in serial, expect output files in /MCGPU/ directory   
	system(ss.str().c_str());
	
	double expected = -1800;
	double energyResult = -1;
	
	std::ifstream infile( std::string(MCGPU + "/methanolCPU.results").c_str() );
	std::cout << std::string(MCGPU + "/bin/methanolCPU.results").c_str() << endl;
	
	
	for( std::string line; getline( infile, line ); )
	{
		std::string str2 ("Final-Energy");
		std::string result;
		found = line.find(str2);
		
		if (found != std::string::npos) {
			result = line.substr(15);
			energyResult  = strtod(result.c_str(), NULL);
			break;
		}
	}
	
	// Clean up
	remove(configFilePath.c_str());
	remove(  std::string(MCGPU + "/methanolCPU.pdb").c_str()  );
	remove(  std::string(MCGPU + "/methanolCPU.results").c_str()  );
	remove(  std::string(MCGPU + "/methanolCPU_100000.state").c_str()  );
	
	EXPECT_NEAR(expected, energyResult, 100 );
}


