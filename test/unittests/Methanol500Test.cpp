#include "Applications/Application.h"
#include "gtest/gtest.h"
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
void createMeoh500ConfigFile(std::string MCGPU, std::string fileName, std::string primaryAtomIndexString) {
	ofstream configFile;
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/Methanol500Tests/" + fileName));
        configFile.open(configFilePath.c_str());
        std::stringstream cfg;
        cfg << ""
            << "#size of periodic box (x, y, z in angstroms)\n"
            << "32.91\n"
            << "32.91\n"
            << "32.91\n"
            << "#temperature in Kelvin\n"
            << "298.15\n"
            << "#max translation\n"
            << ".12\n"
            << "#number of steps\n"
            << "100000\n"
            << "#number of molecules\n"
            << "500\n"
            << "#path to opls.par file\n"
            << MCGPU << "resources/bossFiles/oplsaa.par\n"
			<< "#path to z matrix file\n"
            << MCGPU << "test/unittests/MultipleSolvents/MethanolTests/meoh.z\n"
            << "#path to state input\n"
            << MCGPU << "test/unittests/MultipleSolvents/Methanol500Tests\n"
            << "#path to state output\n"
            << MCGPU << "test/unittests/MultipleSolvents/Methanol500Tests\n"
            << "#pdb output path\n"
            << MCGPU << "test/unittests/MultipleSolvents/Methanol500Tests\n"
            << "#cutoff distance in angstroms\n"
            << "11.0\n"
            << "#max rotation\n"
            << "12.0\n"
            << "#Random Seed Input\n"
            << "12345\n"
            << "#Primary Atom Index\n"
            << primaryAtomIndexString;
        configFile << cfg.str();
        configFile.close();
}

std::string getMCGPU_Meoh500_path() {
	std::string directory = get_current_dir_name();
	std::string mc ("MCGPU");
	std::size_t found = directory.find(mc);
	
	if(found != std::string::npos) {
		directory = directory.substr(0, found + 6);
	}
	
	return (directory + "/");
}

/**
 * Generates a string that is used as a command to run the necessary test.
 * @param MCGPU The path to MCGPU's root.
 * @param configFile The name of the config file to use for the test.
 * @param outputName The name of the simulation, or the file to pipe cerr to (if expecting an error).
 * @param series true if the simulation is to be run in series, false otherwise.
 * @param neighborlist true if the simulation is to be run with a neighborlist, false otherwise.
 * @param errorExpected true if passing behavior for the test throws an error.
 */
std::string buildMeoh500Command(std::string MCGPU, std::string configFile, std::string outputName, bool series, bool neighborlist, bool errorExpected) {
	//Setting up standard build command.
	std::stringstream ss;
	ss << MCGPU << "/bin/metrosim " << MCGPU << "test/unittests/MultipleSolvents/Methanol500Tests/" << configFile << " ";
	
	if(series) {
		ss << "-s --threads 12 ";	//If testing in series, give the flag and specify a the number of threads.
	} else {
		ss << "-p ";				//If testing in parallel, give the corresponding flag.
	}
	
	if(neighborlist) {
		ss << "-l 100 ";				//Add the neighborlist flag if applicable.
	}
	
	if(errorExpected) {
		ss << "-i 10000 > " << MCGPU << "bin/" << outputName << " 2>&1 ";	//If we expect an error, pipe cerr to a textfile where it can be read.
	} else {
		ss << "--name " << outputName << " -i 10000 ";							//If we do not expect an error, simply give the name for the results file.
	}
	std::string output = ss.str();
	std::cout << "RUNNING: " << output << std::endl;
	return output;
}

/**
 * Examines a result file to determine the final energy.
 * @param MCGPU The path to MCGPU's root.
 * @param resultsFile The name of the results file.
 * @return The final energy.
 */
double getMeoh500EnergyResult(std::string MCGPU, std::string resultsFile) {
	std::ifstream infile(std::string(MCGPU + "bin/" + resultsFile).c_str());
	std::size_t found;
	
	for(std::string line; getline(infile, line);) {
		std::string str2 ("Final-Energy");
		std::string result;
		found = line.find(str2);
		if(found != std::string::npos) {
			result = line.substr(15);
			return strtod(result.c_str(), NULL);
		}
	}
	
	return -1;
}

/**
 * Examines a file that cerr was piped to and returns what function call the error occurred in.
 * @param MCGPU The path to MCGPU's root.
 * @param errorFile The name of the file cerr was piped to.
 * @return The function the error occurred in.
 */
std::string getMeoh500ErrorResult(std::string MCGPU, std::string errorFile) {
    std::ifstream infile(std::string(MCGPU + "bin/" + errorFile).c_str());
	std::size_t found;
	
    for(std::string line; getline(infile, line);) {
        std::string str2 ("Error");
        found = line.find(str2);
        if (found != std::string::npos) {
            return line.substr(7,13);
        }
    }
	return "ERROR: COULD NOT PARSE ERROR FILE!";
}

//Test methanol with 500 molecules and a  primary index of 1 on CPU
TEST (Meoh500Test, OnePrimaryIndex)
{   
	std::string MCGPU = getMCGPU_Meoh500_path();
	createMeoh500ConfigFile(MCGPU, "meoh5001MPI.config", "1");
    system(buildMeoh500Command(MCGPU, "meoh5001MPI.config", "meoh5001MPI", true, false, false).c_str());
    double expected = -3454;
    double energyResult = getMeoh500EnergyResult(MCGPU, "meoh5001MPI.results");
    EXPECT_NEAR(expected, energyResult, 100);
}



//Test methanol with 500 molecules and a  primary index of 1 on GPU
TEST (Meoh500Test, OnePrimaryIndex_GPU)
{   
	std::string MCGPU = getMCGPU_Meoh500_path();
    system(buildMeoh500Command(MCGPU, "meoh5001MPI.config", "meoh5001MPI-GPU", false, false, false).c_str());
    double expected = -3454;
    double energyResult = getMeoh500EnergyResult(MCGPU, "meoh5001MPI-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test methanol with 500 molecules and a  primary index of 1 on CPU
//Using neighborlist
TEST (Meoh500Test, NeighborListFunction1MPI)
{
    std::string MCGPU = getMCGPU_Meoh500_path();
    system(buildMeoh500Command(MCGPU, "meoh5001MPI.config", "meoh5001MPI-NL", true, true, false).c_str());
	double expected = -3402;
    double energyResult = getMeoh500EnergyResult(MCGPU, "meoh5001MPI-NL.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test methanol with 500 molecules and a  primary index of 1 on GPU
//Using neighborlist
TEST (Meoh500Test, NeighborListFunction1MPI_GPU)
{
    std::string MCGPU = getMCGPU_Meoh500_path();
    system(buildMeoh500Command(MCGPU, "meoh5001MPI.config", "meoh5001MPI-NL-GPU", false, true, false).c_str());
    double expected = -3402;
    double energyResult = getMeoh500EnergyResult(MCGPU, "meoh5001MPI-NL-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test methanol with 500 molecules and two  primary indexes [1,2] on CPU
TEST (Meoh500Test, TwoPrimaryIndex)
{
	std::string MCGPU = getMCGPU_Meoh500_path();
	createMeoh500ConfigFile(MCGPU, "meoh5002MPI.config", "[1,2]");
    system(buildMeoh500Command(MCGPU, "meoh5002MPI.config", "meoh5002MPI", true, false, false).c_str());
    double expected = -3511;
    double energyResult = getMeoh500EnergyResult(MCGPU, "meoh5002MPI.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test methanol with 500 molecules and two  primary indexes [1,2] on GPU
TEST (Meoh500Test, TwoPrimaryIndex_GPU)
{
	std::string MCGPU = getMCGPU_Meoh500_path();
    system(buildMeoh500Command(MCGPU, "meoh5002MPI.config", "meoh5002MPI-GPU", false, false, false).c_str());
    double expected = -3511;
    double energyResult = getMeoh500EnergyResult(MCGPU, "meoh5002MPI-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test methanol with 500 molecules and primary index of [1,2] on CPU
//Using neighborlist
TEST (Meoh500Test, NeighborListFunction2MPI)
{
    std::string MCGPU = getMCGPU_Meoh500_path();
    system(buildMeoh500Command(MCGPU, "meoh5002MPI.config", "meoh5002MPI-NL", true, true, false).c_str());
    double expected = -3433;
    double energyResult = getMeoh500EnergyResult(MCGPU, "meoh5002MPI-NL.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test methanol with 500 molecules and primary index of [1,2] on GPU
//Using neighborlist
TEST (Meoh500Test, NeighborListFunction2MPI_GPU)
{
    std::string MCGPU = getMCGPU_Meoh500_path();
    system(buildMeoh500Command(MCGPU, "meoh5002MPI.config", "meoh5002MPI-NL-GPU", false, true, false).c_str()); 
	double expected = -3433;
    double energyResult = getMeoh500EnergyResult(MCGPU, "meoh5002MPI-NL-GPU.results");
    EXPECT_NEAR(expected, energyResult, 100);
}

//Test methanol with 500 molecules and multiple solvent primary indexes 1,2  on CPU
TEST (Meoh500Test, MultipleSolventDefinition)
{
    std::string MCGPU = getMCGPU_Meoh500_path();	
	createMeoh500ConfigFile(MCGPU, "meoh500MulSolvent.config", "1,2");
    system(buildMeoh500Command(MCGPU, "meoh500MulSolvent.config", "meoh500MulSolvent.txt", true, false, true).c_str());
    std::string errorResult = getMeoh500ErrorResult(MCGPU, "meoh500MulSolvent.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}        

//Test methanol with 500 molecules and multiple solvent primary indexes 1,2  on GPU
TEST (Meoh500Test, MultipleSolventDefinition_GPU)
{
    std::string MCGPU = getMCGPU_Meoh500_path();
    system(buildMeoh500Command(MCGPU, "meoh500MulSolvent.config", "meoh500MulSolvent-GPU", false, false, true).c_str());
    std::string errorResult = getMeoh500ErrorResult(MCGPU, "meoh500MulSolvent-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test methanol with 500 molecules and multiple solvent primary indexes [1,2],[3,4]  on CPU
TEST (Meoh500Test, MultipleSolventDefinitionMPI)
{
    std::string MCGPU = getMCGPU_Meoh500_path();
	createMeoh500ConfigFile(MCGPU, "meoh500MulSolventMPI.config", "[1,2],[3,4]");
    system(buildMeoh500Command(MCGPU, "meoh500MulSolventMPI.config", "meoh500MulSolventMPI.txt", true, false, true).c_str());
    std::string errorResult = getMeoh500ErrorResult(MCGPU, "meoh500MulSolventMPI.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test methanol with 500 molecules and multiple solvent primary indexes [1,2],[3,4] on GPU
TEST (Meoh500Test, MultipleSolventDefinitionMPI_GPU)
{
    std::string MCGPU = getMCGPU_Meoh500_path();
    system(buildMeoh500Command(MCGPU, "meoh500MulSolventMPI.config", "meoh500MulSolventMPI-GPU.txt", false, false, true).c_str());
    std::string errorResult = getMeoh500ErrorResult(MCGPU, "meoh500MulSolventMPI-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test methanol with 500 molecules and multiple solvent primary indexes 1, [1,2] on CPU
TEST (Meoh500Test, SingleMultipleIndexes)
{
    std::string MCGPU = getMCGPU_Meoh500_path();
	createMeoh500ConfigFile(MCGPU, "meoh500SingleMultipleIndexes.config", "1, [1,2]");
    system(buildMeoh500Command(MCGPU, "meoh500SingleMultipleIndexes.config", "meoh500SingleMultipleIndexes.txt", true, false, true).c_str());
    std::string errorResult = getMeoh500ErrorResult(MCGPU, "meoh500SingleMultipleIndexes.txt");
	EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test methanol with 500 molecules and multiple solvent primary indexes 1, [1,2] on GPU
TEST (Meoh500Test, SingleMultipleIndexes_GPU)
{
    std::string MCGPU = getMCGPU_Meoh500_path();
    system(buildMeoh500Command(MCGPU, "meoh500SingleMultipleIndexes.config", "meoh500SingleMultipleIndexes-GPU.txt", false, false, true).c_str());
    std::string errorResult = getMeoh500ErrorResult(MCGPU, "meoh500SingleMultipleIndexes-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test methanol with 500 molecules and multiple solvent primary indexes [1,2],1 on CPU
TEST (Meoh500Test, SingleMultipleIndexes2)
{
    std::string MCGPU = getMCGPU_Meoh500_path();
	createMeoh500ConfigFile(MCGPU, "meoh500SingleMultipleIndexes2.config", "[1,2],1");
    system(buildMeoh500Command(MCGPU, "meoh500SingleMultipleIndexes2.config", "meoh500SingleMultipleIndexes2.txt", true, false, true).c_str());
    std::string errorResult = getMeoh500ErrorResult(MCGPU, "meoh500SingleMultipleIndexes2.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test methanol with 500 molecules and multiple solvent primary indexes [1,2],1 on GPU
TEST (Meoh500Test, SingleMultipleIndexes2_GPU)
{
    std::string MCGPU = getMCGPU_Meoh500_path();
    system(buildMeoh500Command(MCGPU, "meoh500SingleMultipleIndexes2.config", "meoh500SingleMultipleIndexes2-GPU.txt", false, false, true).c_str());
    std::string errorResult = getMeoh500ErrorResult(MCGPU, "meoh500SingleMultipleIndexes2-GPU.txt");
    EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}
