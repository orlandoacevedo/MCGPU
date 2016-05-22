
#include "Metropolis/Utilities/FileUtilities.h"
#include "gtest/gtest.h"

/**
 * Test the configuration scanner
 *
 * The ConfigScanner reads .config files and sets up the simulation based on
 * those values. Test that those values are properly read and processed.
 */
TEST(IOTests, ConfigScan)
{
	const double ERROR_MARGIN = 0.0001;

	// Find the MCGPU project path
    size_t size = 4096;
    char *path = (char*) malloc(size);
    path = getcwd(path, size);
    std::string directory(path);
    std::string mc ("MCGPU");
    std::size_t found = directory.find(mc);

    if (found != std::string::npos) {
        directory = directory.substr(0,found+6);
    }
    std::string MCGPU = directory;

    string configPath = MCGPU;
    configPath.append("stuff/configurationTest.txt");
    string oplsPath = "path/to/opla.par/file";
    string zMatrixPath = "path/to/zMatrix/file";
    string stateInputPath = "path/to/state/input";
    string stateOutputPath = "path/to/State/output";
    string pdbOutputPath = "path/to/pdb/output";
	string simName = "MyTestSimulation";
    ConfigScanner cs;
    cs.readInConfig(configPath);

    // Test various path elements
	EXPECT_EQ(configPath, cs.getConfigPath());
	EXPECT_EQ(oplsPath, cs.getOplsusaparPath());
	EXPECT_EQ(zMatrixPath, cs.getZmatrixPath());
	EXPECT_EQ(stateInputPath, cs.getStatePath());
	EXPECT_EQ(stateOutputPath, cs.getStateOutputPath());
	EXPECT_EQ(pdbOutputPath, cs.getPdbOutputPath());
	EXPECT_EQ(simName, cs.getSimulationName());

    // Test box dimensions
    double x = 10.342;
    double y = 1234.45;
    double z = 100.3;
    double temperature = 293;
    double maxTranslation = .5;
    int numberOfSteps = 10000;
    int numberOfMolecules = 500;

    Environment* enviro = cs.getEnviro();
	EXPECT_NEAR(x, enviro->x, ERROR_MARGIN);
	EXPECT_NEAR(y, enviro->y, ERROR_MARGIN);
	EXPECT_NEAR(z, enviro->z, ERROR_MARGIN);
	EXPECT_NEAR(maxTranslation, enviro->maxTranslation, ERROR_MARGIN);
	EXPECT_NEAR(numberOfMolecules, enviro->numOfMolecules, ERROR_MARGIN);
	EXPECT_NEAR(temperature, enviro->temp, ERROR_MARGIN);
	EXPECT_NEAR(numberOfSteps, cs.getSteps(), ERROR_MARGIN);
}
