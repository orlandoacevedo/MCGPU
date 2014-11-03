
#include "Metropolis/Utilities/FileUtilities.h"
#include "gtest/gtest.h"

// Descr: evident
// Implementation details: See gtest/samples for GTest syntax and usage
TEST(IOTests, ConfigScan)
{
  
    string configPath = "bin/configurationTest.txt"; 
    string oplsPath = "path/to/opls/file";
    string zMatrixPath = "path/to/zMatrix/file";
    string stateInputPath = "path/to/state/input";
    string stateOutputPath = "path/to/State/output";
    string pdbOutputPath = "path/to/pdb/output";
    ConfigScanner cs;
    cs.readInConfig(configPath);
	
	cout << cs.getConfigPath() << endl;
	cout << cs.getOplsusaparPath() << endl;
    //test configuration path
	EXPECT_EQ(0, configPath.compare(cs.getConfigPath())  );
    //test opls path
	EXPECT_EQ(0, oplsPath.compare(cs.getOplsusaparPath())  );
    //test zmatrix path
	EXPECT_EQ(0, zMatrixPath.compare(cs.getZmatrixPath())  );
    //test state input path
	EXPECT_EQ(0, stateInputPath.compare(cs.getStatePath())  );
    //test state ouput path
	EXPECT_EQ(0, stateOutputPath.compare(cs.getStateOutputPath())  );
    //test pdb output
	EXPECT_EQ(0, pdbOutputPath.compare(cs.getPdbOutputPath())  );
	
    //test box dimensions
    Environment* enviro = cs.getEnviro();
    double x = 10.342;
    double y = 1234.45;
    double z = 100.3;
    double temperature = 293;
    double maxTranslation = .5;
    int numberOfSteps = 10000;
    int numberOfMolecules = 500;
    
    
    //EXPECT_EQ should work, and tester even shows values are the same but fails anyway
    // Best luck, m8.
	EXPECT_NEAR(x, enviro->x, .0001 );
	EXPECT_NEAR(y, enviro->y, .0001);
	EXPECT_NEAR(z, enviro->z, .0001 );
	
	EXPECT_NEAR(maxTranslation, enviro->maxTranslation, .0001 );
	EXPECT_NEAR(numberOfMolecules, enviro->numOfMolecules, .0001 );
	EXPECT_NEAR(temperature, enviro->temp, .0001 );
	EXPECT_NEAR(numberOfSteps, cs.getSteps(), .0001 );

    
}

