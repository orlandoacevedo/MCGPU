#include "Applications/Application.h"
#include "TestUtil.h"
#include "gtest/gtest.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

//Test t3pdim with 1 primary index on CPU
//Should result in an error since t3pdim has two molecules
TEST (t3pdimTest, OnePrimaryIndex) 
{
		std::size_t found;
		std::string MCGPU = getMCGPU_path();

        //create test config file 
        //hardcode since file path will change on each user

        ofstream configFile;
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/t3pdimTests/t3pdim1MPI.config"));
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
                << "#path to opls.par file\n"
                << MCGPU << "resources/bossFiles/oplsaa.par\n"
                << "#path to z matrix file\n"
	        << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests/t3pdim.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#cutoff distance in angrstoms\n"
                << "11.0\n"
                << "#max rotation\n"
                << "12.0\n"
                << "#Random Seed Input\n"
                << "12345\n"
                << "#Primary Atom Index\n"
                << "1";
        configFile << cfg.str();
        configFile.close();

        //Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimTests/t3pdim1MPI.config -s --threads 12 -i 10000 >"
                << " "
                << MCGPU << "/bin/t3pdim1MPI.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/t3pdim1MPI.txt").c_str());

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

//Test t3pdim with 1 primary index on GPU
//Should result in an error since t3pdim has two molecules
TEST (t3pdimTest, OnePrimaryIndexGPU)
{

        string directory = get_current_dir_name();
        std::string mc ("MCGPU");
        std::size_t found = directory.find(mc);

        if(found != std::string::npos) {
                directory = directory.substr(0,found+6);

        }

        std::string MCGPU = directory;

        //Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimTests/t3pdim1MPI.config -p --threads 12 -i 10000 >"
                << " "
                << MCGPU << "/bin/t3pdim1MPI-GPU.txt 2>&1";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/t3pdim1MPI-GPU.txt").c_str());

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

 
//Test t3pdim with [1,2] primary index on CPU
//Should result in an error since t3pdim has two molecules
TEST (t3pdimTest, TwoPrimaryIndex)
{

        string directory = get_current_dir_name();
        std::string mc ("MCGPU");
        std::size_t found = directory.find(mc);

        if(found != std::string::npos) {
                directory = directory.substr(0,found+6);

        }

        std::string MCGPU = directory;

        //create test config file 
        //hardcode since file path will change on each user

        ofstream configFile;
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/t3pdimTests/t3pdim2MPI.config"));
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
                << "#path to opls.par file\n"
                << MCGPU << "resources/bossFiles/oplsaa.par\n"
                << "#path to z matrix file\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests/t3pdim.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#cutoff distance in angrstoms\n"
                << "11.0\n"
                << "#max rotation\n"
                << "12.0\n"
                << "#Random Seed Input\n"
                << "12345\n"
                << "#Primary Atom Index\n"
                << "[1,2]";
        configFile << cfg.str();
        configFile.close();

        
        //Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimTests/t3pdim2MPI.config -s --threads 12 -i 10000 >"
                << " "
                << MCGPU << "/bin/t3pdim2MPI.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/t3pdim2MPI.txt").c_str());

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

//Test t3pdim with [1,2] primary index on GPU
//Should result in an error since t3pdim has two molecules
TEST (t3pdimTest, TwoPrimaryIndexGPU)
{

        string directory = get_current_dir_name();
        std::string mc ("MCGPU");
        std::size_t found = directory.find(mc);

        if(found != std::string::npos) {
                directory = directory.substr(0,found+6);

        }

        std::string MCGPU = directory;

        //Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimTests/t3pdim2MPI.config -p --threads 12 -i 10000 >"
                << " "
                << MCGPU << "/bin/t3pdim2MPI-GPU.txt 2>&1";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/t3pdim2MPI-GPU.txt").c_str());

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


//Test t3pdim with multiple solvents with indexes 1,2 on CPU
TEST (t3pdimTest, MultipleSolventDefinition)
{

        string directory = get_current_dir_name();
        std::string mc ("MCGPU");
        std::size_t found = directory.find(mc);

        if(found != std::string::npos) {
                directory = directory.substr(0,found+6);

        }

        std::string MCGPU = directory;

        //create test config file 
        //hardcode since file path will change on each user

        ofstream configFile;
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/t3pdimTests/t3pdimMulSolvents.config"));
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
                << "#path to opls.par file\n"
                << MCGPU << "resources/bossFiles/oplsaa.par\n"
                << "#path to z matrix file\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests/t3pdim.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#cutoff distance in angrstoms\n"
                << "11.0\n"
                << "#max rotation\n"
                << "12.0\n"
                << "#Random Seed Input\n"
                << "12345\n"
                << "#Primary Atom Index\n"
                << "1,2";
        configFile << cfg.str();
        configFile.close();

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimTests/t3pdimMulSolvents.config -s --threads 12 --name t3pdimMulSolvents -i 10000";

        //launch MCGPU in serial
        system(ss.str().c_str());
        double expected = -1832;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/t3pdimMulSolvents.results").c_str());

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


//Test t3pdim with multiple solvents with indexes 1,2 on GPU
TEST (t3pdimTest, MultipleSolventDefinitionGPU)
{

        string directory = get_current_dir_name();
        std::string mc ("MCGPU");
        std::size_t found = directory.find(mc);

        if(found != std::string::npos) {
                directory = directory.substr(0,found+6);

        }

        std::string MCGPU = directory;

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimTests/t3pdimMulSolvents.config -p --threads 12 --name t3pdimMulSolvents-GPU -i 10000";

       
        system(ss.str().c_str());
        double expected = -1832;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/t3pdimMulSolvents-GPU.results").c_str());

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

//Test t3pdim with multiple solvents with indexes 1,2 on CPU
//Using neighborlist
TEST (t3pdimTest, NeighborListFunctionMulSolvents)
{

        string directory = get_current_dir_name();
        std::string mc ("MCGPU");
        std::size_t found = directory.find(mc);

        if(found != std::string::npos) {
                directory = directory.substr(0,found+6);

        }

        std::string MCGPU = directory;

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimTests/t3pdimMulSolvents.config -s --threads 12 --name t3pdimMulSolvents-NL -i 10000 -l";


        system(ss.str().c_str());
        double expected = -1727;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/t3pdimMulSolvents-NL.results").c_str());

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

//Test t3pdim with multiple solvents with indexes 1,2 on GPU
//Using neighborlist
TEST (t3pdimTest, NeighborListFunctionMulSolvents_GPU)
{

        string directory = get_current_dir_name();
        std::string mc ("MCGPU");
        std::size_t found = directory.find(mc);
        
        if(found != std::string::npos) {
                directory = directory.substr(0,found+6);

        }

        std::string MCGPU = directory;

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimTests/t3pdimMulSolvents.config -p --threads 12 --name t3pdimMulSolvents-NL-GPU -i 10000 -l";


        system(ss.str().c_str());
        double expected = -1727;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/t3pdimMulSolvents-NL-GPU.results").c_str());

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

//test t3pdim for multiple solvents with [1,2],[3,4] primary indexes on CPU
TEST (t3pdimTest, MultipleSolventDefinitionMPI)
{

        string directory = get_current_dir_name();
        std::string mc ("MCGPU");
        std::size_t found = directory.find(mc);

        if(found != std::string::npos) {
                directory = directory.substr(0,found+6);

        }

        std::string MCGPU = directory;

        //create test config file 
        //hardcode since file path will change on each user

        ofstream configFile;
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/t3pdimTests/t3pdimMulSolventsMPI.config"));
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
                << "#path to opls.par file\n"
                << MCGPU << "resources/bossFiles/oplsaa.par\n"
                << "#path to z matrix file\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests/t3pdim.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#cutoff distance in angstroms\n"
                << "11.0\n"
                << "#max rotation\n"
                << "12.0\n"
                << "#Random Seed Input\n"
                << "12345\n"
                << "#Primary Atom Index\n"
                << "[1,2],[3,4]";
        configFile << cfg.str();
        configFile.close();

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimTests/t3pdimMulSolventsMPI.config -s --threads 12 --name t3pdimMulSolventsMPI -i 10000";

        //launch MCGPU in serial
        system(ss.str().c_str());
        double expected = -1790;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/t3pdimMulSolventsMPI.results").c_str());

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


//Test t3pdim with multiple solvents with indexes [1,2],[3,4] on GPU
TEST (t3pdimTest, MultipleSolventDefinitionMPI_GPU)
{

        string directory = get_current_dir_name();
        std::string mc ("MCGPU");
        std::size_t found = directory.find(mc);

        if(found != std::string::npos) {
                directory = directory.substr(0,found+6);

        }

        std::string MCGPU = directory;

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimTests/t3pdimMulSolventsMPI.config -p --name t3pdimMulSolventsMPI-GPU -i 10000";


        system(ss.str().c_str());
        double expected = -1790;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/t3pdimMulSolventsMPI-GPU.results").c_str());

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

//test t3pdim for multiple solvents with (1,[1,2]) primary indexes on CPU
TEST (t3pdimTest, SingleMultipleIndexes)
{

        string directory = get_current_dir_name();
        std::string mc ("MCGPU");
        std::size_t found = directory.find(mc);

        if(found != std::string::npos) {
                directory = directory.substr(0,found+6);

        }

        std::string MCGPU = directory;

        //create test config file 
        //hardcode since file path will change on each user

        ofstream configFile;
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/t3pdimTests/t3pdimSingleMultipleIndexes.config"));
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
                << "#path to opls.par file\n"
                << MCGPU << "resources/bossFiles/oplsaa.par\n"
                << "#path to z matrix file\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests/t3pdim.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#pdb output path\n"
		<< MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#cutoff distance in angstroms\n"
                << "11.0\n"
                << "#max rotation\n"
                << "12.0\n"
                << "#Random Seed Input\n"
                << "12345\n"
                << "#Primary Atom Index\n"
                << "1,[1,2]";
        configFile << cfg.str();
        configFile.close();

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimTests/t3pdimSingleMultipleIndexes.config -s --threads 12 --name t3pdimSingleMultipleIndexes -i 10000";

        //launch MCGPU in serial
        system(ss.str().c_str());
        double expected = -1770;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/t3pdimSingleMultipleIndexes.results").c_str());

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

//test t3pdim for multiple solvents with (1,[1,2]) primary indexes on GPU
TEST (t3pdimTest, SingleMultipleIndexes_GPU)
{

        string directory = get_current_dir_name();
        std::string mc ("MCGPU");
        std::size_t found = directory.find(mc);

        if(found != std::string::npos) {
                directory = directory.substr(0,found+6);

        }

        std::string MCGPU = directory;

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimTests/t3pdimSingleMultipleIndexes.config -p --threads 12 --name t3pdimSingleMultipleIndexes-GPU -i 10000";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        double expected = -1770;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/t3pdimSingleMultipleIndexes-GPU.results").c_str());

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

//test t3pdim for multiple solvents with ([1,2], 1) primary indexes on CPU
TEST (t3pdimTest, SingleMultipleIndexes2)
{

        string directory = get_current_dir_name();
        std::string mc ("MCGPU");
        std::size_t found = directory.find(mc);

        if(found != std::string::npos) {
                directory = directory.substr(0,found+6);

        }

        std::string MCGPU = directory;

        //create test config file 
        //hardcode since file path will change on each user

        ofstream configFile;
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/t3pdimTests/t3pdimSingleMultipleIndexes2.config"));
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
                << "#path to opls.par file\n"
                << MCGPU << "resources/bossFiles/oplsaa.par\n"
                << "#path to z matrix file\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests/t3pdim.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/t3pdimTests\n"
                << "#cutoff distance in angrstoms\n"
                << "11.0\n"
                << "#max rotation\n"
                << "12.0\n"
                << "#Random Seed Input\n"
                << "12345\n"
                << "#Primary Atom Index\n"
                << "[1,2],1";
        configFile << cfg.str();
        configFile.close();

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimTests/t3pdimSingleMultipleIndexes2.config -s --threads 12 --name t3pdimSingleMultipleIndexes2 -i 10000";

        //launch MCGPU in serial
        system(ss.str().c_str());
        double expected = -1800;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/t3pdimSingleMultipleIndexes2.results").c_str());

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

//test t3pdim for multiple solvents with ([1,2], 1) primary indexes on GPU
TEST (t3pdimTest, SingleMultipleIndexes2_GPU)
{

        string directory = get_current_dir_name();
        std::string mc ("MCGPU");
        std::size_t found = directory.find(mc);

        if(found != std::string::npos) {
                directory = directory.substr(0,found+6);

        }

        std::string MCGPU = directory;

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/t3pdimTests/t3pdimSingleMultipleIndexes2.config -p --threads 12 --name t3pdimSingleMultipleIndexes2-GPU -i 10000";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        double expected = -1800;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/t3pdimSingleMultipleIndexes2-GPU.results").c_str());

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

        EXPECT_NEAR(expected, energyResult,100 );
}
