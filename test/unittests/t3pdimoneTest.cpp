#include "Applications/Application.h"
#include "gtest/gtest.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

void createConfigFile(std::string MCGPU_path, std::string fileName, std::string primaryAtomIndexString) {
	ofstream configFile;
	std::string configFilePath (std::string (MCGPU + "test/unittests/MultipleSolvents/t3pdimoneTests/" + fileName));
	configFile.open(configFilePath.c_str());
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
		<< "#number of molecules\n"
		<< "256\n"
		<< "#path to opla.par file\n"
		<< MCGPU << "resources/bossFiles/oplsaa.par\n"
		<< "#path to z matrix file\n"
		<< MCGPU << "test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimone.z\n"
		<< "#path to state input\n"
		<< MCGPU << "test/unittests/MultipleSolvents/t3pdimoneTests\n"
        << "#path to state output\n"
        << MCGPU << "test/unittests/MultipleSolvents/t3pdimoneTests\n"
        << "#pdb output path\n"
        << MCGPU << "test/unittests/MultipleSolvents/t3pdimoneTests\n"
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

std::string getMCGPU_path () {
	string directory = get_current_dir_name();
	std::string mc ("MCGPU");
	std::size_t found = directory.find(mc);
	
	if(found != std::string::npos) {
		directory = directory.substr(0, found + 6);
	}
	std::string MCGPU = directory + "/";
	return MCGPU;
}

//Test t3pdimone with 1 primary index on CPU
TEST (t3pdimoneTest, OnePrimaryIndex)
{
        std::string MCGPU = getMCGPU_path();

        //create test config file 
        //hardcode since file path will change on each user

		createConfigFile(MCGPU, "t3pdimone1MPI.config", "1");		
		
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimone1MPI.config -s --threads 12 --name t3pdimone1MPI -i 10000";

        //launch MCGPU in serial
        system(ss.str().c_str());
        double expected = -1810;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/t3pdimone1MPI.results").c_str());

        for(std::string line; getline(infile, line);) {
                std::string str2 ("Final-Energy");
                std::string result;
                found = line.find(str2);
                if (found != std::string::npos) {
                        result = line.substr(15);
                        energyResult = strtod(result.c_str(), NULL);
                        break;
                }
        }

        EXPECT_NEAR(expected, energyResult, 100);
}

//Test t3pdimone with 1 primary index on GPU
TEST (t3pdimoneTest, PrimaryIndexGPU)
{
	std::string MCGPU = getMCGPU_path();

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimone1MPI.config -p --name t3pdimone1MPI-GPU -i 10000";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        double expected = -1810;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/t3pdimone1MPI-GPU.results").c_str());

        for(std::string line; getline(infile, line);) {
                std::string str2 ("Final-Energy");
                std::string result;
                found = line.find(str2);

                if (found != std::string::npos) {
                        result = line.substr(15);
                        energyResult = strtod(result.c_str(), NULL);
                        break;
                }
        }

        EXPECT_NEAR(expected, energyResult, 100);
}

//Test t3pdimone with 1 primary index on CPU
//Using neighborlist
TEST (t3pdimoneTest, NeighborListFunction1MPI)
{
        std::string MCGPU = getMCGPU_path();

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimone1MPI.config -s --threads 12 --name t3pdimone1MPI-NL -i 10000 -l";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        double expected = -1810;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/t3pdimone1MPI-NL.results").c_str());

        for(std::string line; getline(infile, line);) {
                std::string str2 ("Final-Energy");
                std::string result;
                found = line.find(str2);

                if (found != std::string::npos) {
                        result = line.substr(15);
                        energyResult = strtod(result.c_str(), NULL);
                        break;
                }
        }

        EXPECT_NEAR(expected, energyResult, 100);
}

//Test t3pdimone with 1 primary index on GPU
//Using neighborlist
TEST (t3pdimoneTest, NeighborListFunction1MPI_GPU)
{
        std::string MCGPU = getMCGPU_path();

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimone1MPI.config -p --threads 12 --name t3pdimone1MPI-NL-GPU -i 10000 -l";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        double expected = -1810;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/t3pdimone1MPI-NL-GPU.results").c_str());

        for(std::string line; getline(infile, line);) {
                std::string str2 ("Final-Energy");
                std::string result;
                found = line.find(str2);

                if (found != std::string::npos) {
                        result = line.substr(15);
                        energyResult = strtod(result.c_str(), NULL);
                        break;
                }
        }

        EXPECT_NEAR(expected, energyResult, 100);
}

//Test t3pdimone with 1 primary [1,2] index on CPU
TEST (t3pdimoneTest, TwoPrimaryIndex)
{
        std::string MCGPU = getMCGPU_path();

        //create test config file 
        //hardcode since file path will change on each user
		createConfigFile(MCGPU, "t3pdimone2MPI.config", "[1,2]");

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimone2MPI.config -s --threads 12 --name t3pdimone2MPI -i 10000";

        //launch MCGPU in serial
        system(ss.str().c_str());
        double expected = -1770;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/t3pdimone2MPI.results").c_str());

        for(std::string line; getline(infile, line);) {
                std::string str2 ("Final-Energy");
                std::string result;
                found = line.find(str2);

                if (found != std::string::npos) {
                        result = line.substr(15);
                        energyResult = strtod(result.c_str(), NULL);
                        break;
                }
        }

        EXPECT_NEAR(expected, energyResult, 100);
}

//Test t3pdimone with 1 primary [1,2] index on GPU
TEST (t3pdimoneTest, TwoPrimaryIndexGPU)
{
        std::string MCGPU = getMCGPU_path();

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimone2MPI.config -p --threads 12 --name t3pdimone2MPI-GPU -i 10000";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        double expected = -1770;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/t3pdimone2MPI-GPU.results").c_str());

        for(std::string line; getline(infile, line);) {
                std::string str2 ("Final-Energy");
                std::string result;
                found = line.find(str2);

                if (found != std::string::npos) {
                        result = line.substr(15);
                        energyResult = strtod(result.c_str(), NULL);
                        break;
                }
        }

        EXPECT_NEAR(expected, energyResult, 100);
}

//Test t3pdimone with multiple solvents 1,2 primary indexes on CPU
TEST (t3pdimoneTest, MultipleSolventDefinition)
{

        std::string MCGPU = getMCGPU_path();

        //create test config file 
        //hardcode since file path will change on each user

		createConfigFile(MCGPU, "t3pdimoneMulSolvent.config", "1,2");
	
        //Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimoneMulSolvent.config -s -i 10000 >"
                << " "
                << MCGPU << "/bin/t3pdimoneMulSolvent.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/t3pdimoneMulSolvent.txt").c_str());

        for(std::string line; getline(infile, line);) {
                std::string str2 ("Error");
                found = line.find(str2);
                if (found != std::string::npos) {
                       errorResult = line.substr(7,13);
                        break;
                }
        }

        EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}
       
//Test t3pdimone with multiple solvents 1,2 primary indexes on GPU
TEST (t3pdimoneTest, MultipleSolventDefinition_GPU)
{
        std::string MCGPU = getMCGPU_path();

        //Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimoneMulSolvent.config -p -i 10000 >"
                << " "
                << MCGPU << "/bin/t3pdimoneMulSolvent-GPU.txt 2>&1";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/t3pdimoneMulSolvent-GPU.txt").c_str());

        for(std::string line; getline(infile, line);) {
                std::string str2 ("Error");
                found = line.find(str2);
                if (found != std::string::npos) {
                       errorResult = line.substr(7,13);
                        break;
                }
        }

        EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test t3pdimone with multiple solvents ([1,2],[3,4]) primary indexes on CPU
TEST (t3pdimoneTest, MultipleSolventDefinitionMPI)
{
        std::string MCGPU = getMCGPU_path();

        //create test config file 
        //hardcode since file path will change on each user

		createConfigFile(MCGPU, "t3pdimoneMulSolventMPI.config", "[1,2],[3,4]");

        //Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimoneMulSolventMPI.config -s -i 10000 >"
                << " "
                << MCGPU << "/bin/t3pdimoneMulSolventMPI.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/t3pdimoneMulSolventMPI.txt").c_str());

        for(std::string line; getline(infile, line);) {
                std::string str2 ("Error");
                found = line.find(str2);
                if (found != std::string::npos) {
                       errorResult = line.substr(7,13);
                        break;
                }
        }

        EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test t3pdimone with multiple solvents ([1,2],[3,4]) primary indexes on GPU
TEST (t3pdimoneTest, MultipleSolventDefinitionMPI_GPU)
{
        std::string MCGPU = getMCGPU_path();

        //Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimoneMulSolventMPI.config -p -i 10000 >"
                << " "
                << MCGPU << "/bin/t3pdimoneMulSolventMPI-GPU.txt 2>&1";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/t3pdimoneMulSolventMPI-GPU.txt").c_str());

        for(std::string line; getline(infile, line);) {
                std::string str2 ("Error");
                found = line.find(str2);
                if (found != std::string::npos) {
                       errorResult = line.substr(7,13);
                        break;
                }
        }

        EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test t3pdimone with multiple solvents (1,[1,2]) primary indexes on CPU
TEST (t3pdimoneTest, SingleMultipleIndexes)
{
        std::string MCGPU = getMCGPU_path();

        //create test config file 
        //hardcode since file path will change on each user

		createConfigFile(MCGPU, "t3pdimoneSingleMultipleIndexes.config", "1,[1,2]");

        //Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimoneSingleMultipleIndexes.config -s -i 10000 >"
                << " "
                << MCGPU << "/bin/t3pdimoneSingleMultipleIndexes.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/t3pdimoneSingleMultipleIndexes.txt").c_str());

        for(std::string line; getline(infile, line);) {
                std::string str2 ("Error");
                found = line.find(str2);
                if (found != std::string::npos) {
                       errorResult = line.substr(7,13);
                        break;
                }
        }

        EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test t3pdimone with multiple solvents (1,[1,2]) primary indexes on GPU
TEST (t3pdimoneTest, SingleMultipleIndexes_GPU)
{
        std::string MCGPU = getMCGPU_path();

        //Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimoneSingleMultipleIndexes.config -p -i 10000 >"
                << " "
                << MCGPU << "/bin/t3pdimoneSingleMultipleIndexes-GPU.txt 2>&1";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/t3pdimoneSingleMultipleIndexes-GPU.txt").c_str());

        for(std::string line; getline(infile, line);) {
                std::string str2 ("Error");
                found = line.find(str2);
                if (found != std::string::npos) {
                       errorResult = line.substr(7,13);
                        break;
                }
        }

        EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test t3pdimone with multiple solvents ([1,2],1) primary indexes on CPU
TEST (t3pdimoneTest, SingleMultipleIndexes2)
{
        std::string MCGPU = getMCGPU_path();

        //create test config file 
        //hardcode since file path will change on each user

		createConfigFile(MCGPU, "t3pdimoneSingleMultipleIndexes2.config", "[1,2],1");

        //Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimoneSingleMultipleIndexes2.config -s -i 10000 >"
                << " "
                << MCGPU << "/bin/t3pdimoneSingleMultipleIndexes2.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/t3pdimoneSingleMultipleIndexes2.txt").c_str());

        for(std::string line; getline(infile, line);) {
                std::string str2 ("Error");
                found = line.find(str2);
                if (found != std::string::npos) {
                       errorResult = line.substr(7,13);
                        break;
                }
        }

        EXPECT_STREQ("loadBoxData()", errorResult.c_str());
}

//Test t3pdimone with multiple solvents ([1,2],1) primary indexes on GPU
TEST (t3pdimoneTest, SingleMultipleIndexes2_GPU)
{
        std::string MCGPU = getMCGPU_path();

        //Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimoneTests/t3pdimoneSingleMultipleIndexes2.config -p -i 10000 >"
                << " "
                << MCGPU << "/bin/t3pdimoneSingleMultipleIndexes2-GPU.txt 2>&1";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/t3pdimoneSingleMultipleIndexes2-GPU.txt").c_str());

        for(std::string line; getline(infile, line);) {
                std::string str2 ("Error");
                found = line.find(str2);
                if (found != std::string::npos) {
                       errorResult = line.substr(7,13);
                        break;
                }
        }

        EXPECT_STREQ("loadBoxData()", errorResult.c_str());
} 
