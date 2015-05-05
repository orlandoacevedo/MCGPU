#include "Applications/Application.h"
#include "gtest/gtest.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

//Test meohtip3 with 1 primary index on CPU
TEST (meohtip3Test, OnePrimaryIndex)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip31MPI.config"));
        configFile.open(configFilePath.c_str());
        std::stringstream cfg;
        cfg << ""
                << "#size of periodic box (x, y, z in angstroms)\n"
                << "23.7856\n"
                << "23.7856\n"
                << "23.7856\n"
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
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests/meohtip3.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
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
                << MCGPU << "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip31MPI.config -s -i 10000 >"
                << " "
                << MCGPU << "/bin/meohtip31MPI.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/meohtip31MPI.txt").c_str());

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

//Test meohtip3 with 1 primary index on GPU
TEST (meohtip3Test, OnePrimaryIndexGPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip31MPI.config -p -i 10000 >"
                << " "
                << MCGPU << "/bin/meohtip31MPI-GPU.txt 2>&1";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/meohtip31MPI-GPU.txt").c_str());

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

//Test meohtip3 with 2 primary indexes [1,2] on CPU
TEST (meohtip3Test, TwoPrimaryIndex)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip32MPI.config"));
        configFile.open(configFilePath.c_str());
        std::stringstream cfg;
        cfg << ""
                << "#size of periodic box (x, y, z in angstroms)\n"
                << "23.7856\n"
                << "23.7856\n"
                << "23.7856\n"
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
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests/meohtip3.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
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
                << MCGPU << "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip32MPI.config -s -i 10000 >"
                << " "
                << MCGPU << "/bin/meohtip32MPI.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/meohtip32MPI.txt").c_str());

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

//Test meohtip3 with 2 primary indexes [1,2] on GPU
TEST (meohtip3Test, TwoPrimaryIndexGPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip32MPI.config -p -i 10000 >"
                << " "
                << MCGPU << "/bin/meohtip32MPI-GPU.txt 2>&1";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/meohtip32MPI-GPU.txt").c_str());

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

//test meohtip3 with multiple solvents 1,2 primary indexes on CPU
TEST (meohtip3Test, MultpleSolventDefinition)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip3MulSolvent.config"));
        configFile.open(configFilePath.c_str());
        std::stringstream cfg;
        cfg << ""
                << "#size of periodic box (x, y, z in angstroms)\n"
                << "23.7856\n"
                << "23.7856\n"
                << "23.7856\n"
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
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests/meohtip3.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
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
                << MCGPU << "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip3MulSolvent.config -s --threads 12 --name meohtip3MulSolvent -i 10000";

        //launch MCGPU in serial
        system(ss.str().c_str());
        double expected = -1980;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meohtip3MulSolvent.results").c_str());

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

//test meohtip3 with multiple solvents 1,2 primary indexes on GPU
TEST (meohtip3Test, MultpleSolventDefinition_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip3MulSolvent.config -p --threads 12 --name meohtip3MulSolvent-GPU -i 10000";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        double expected = -1980;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meohtip3MulSolvent-GPU.results").c_str());

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

//Test meohtip3 with multiple solvents 1,2 primary indexes on CPU
//Using neighborlist
TEST (meohtip3Test, NeighborListFunctionMulSolvent) 
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
                << MCGPU << "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip3MulSolvent.config -s --threads 12 --name meohtip3MulSolvent-NL -i 10000 -l";

        //launch MCGPU in serial
        system(ss.str().c_str());
        double expected = -2400;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meohtip3MulSolvent-NL.results").c_str());

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

//Test meohtip3 with multiple solvents 1,2 primary indexes on GPU
//Using neighborlist
TEST (meohtip3Test,NeighborListFunctionMulSolvent_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip3MulSolvent.config -p --threads 12 --name meohtip3MulSolvent-NL-GPU -i 10000 -l";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        double expected = -2400;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meohtip3MulSolvent-NL-GPU.results").c_str());

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


//test meohtip3 with multiple solvents [1,2],[3,4] primary indexes on CPU
TEST (meohtip3Test, MultpleSolventDefinitionMPI)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip3MulSolventMPI.config"));
        configFile.open(configFilePath.c_str());
        std::stringstream cfg;
        cfg << ""
                << "#size of periodic box (x, y, z in angstroms)\n"
                << "23.7856\n"
                << "23.7856\n"
                << "23.7856\n"
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
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests/meohtip3.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
		<< "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
                << "#cutoff distance in angrstoms\n"
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
                << MCGPU << "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip3MulSolventMPI.config -s --threads 12 --name meohtip3MulSolventMPI -i 10000";

        //launch MCGPU in serial
        system(ss.str().c_str());
        double expected = -2040;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meohtip3MulSolventMPI.results").c_str());

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

//test meohtip3 with multiple solvents [1,2],[3,4] primary indexes on GPU
TEST (meohtip3Test, MultpleSolventDefinitionMPI_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip3MulSolventMPI.config -p --threads 12 --name meohtip3MulSolventMPI-GPU -i 10000";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        double expected = -2040;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meohtip3MulSolventMPI-GPU.results").c_str());

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

//test meohtip3 with multiple solvents (1,[1,2]) primary indexes on CPU
TEST (meohtip3Test, SingleMultipleIndexes)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip3SingleMultipleIndexes.config"));
        configFile.open(configFilePath.c_str());
        std::stringstream cfg;
        cfg << ""
                << "#size of periodic box (x, y, z in angstroms)\n"
                << "23.7856\n"
                << "23.7856\n"
                << "23.7856\n"
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
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests/meohtip3.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
                << "#cutoff distance in angrstoms\n"
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
                << MCGPU << "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip3SingleMultipleIndexes.config -s --threads 12 --name meohtip3SingleMultipleIndexes -i 10000";

        //launch MCGPU in serial
        system(ss.str().c_str());
        double expected = -1990;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meohtip3SingleMultipleIndexes.results").c_str());

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

//test meohtip3 with multiple solvents (1,[1,2]) primary indexes on GPU
TEST (meohtip3Test, SingleMultipleIndexes_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip3SingleMultipleIndexes.config -s --name meohtip3SingleMultipleIndexes-GPU -i 10000";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        double expected = -1990;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meohtip3SingleMultipleIndexes-GPU.results").c_str());

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

//test meohtip3 with multiple solvents ([1,2], 1) primary indexes on CPU
TEST (meohtip3Test, SingleMultipleIndexes2)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip3SingleMultipleIndexes2.config"));
        configFile.open(configFilePath.c_str());
        std::stringstream cfg;
        cfg << ""
                << "#size of periodic box (x, y, z in angstroms)\n"
                << "23.7856\n"
                << "23.7856\n"
                << "23.7856\n"
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
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests/meohtip3.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/meohtip3Tests\n"
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
                << MCGPU << "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip3SingleMultipleIndexes2.config -s --name meohtip3SingleMultipleIndexes2 -i 10000";

        //launch MCGPU in serial
        system(ss.str().c_str());
        double expected = -2010;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meohtip3SingleMultipleIndexes2.results").c_str());

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

//test meohtip3 with multiple solvents ([1,2], 1) primary indexes on GPU
TEST (meohtip3Test, SingleMultipleIndexes2_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/meohtip3Tests/meohtip3SingleMultipleIndexes2.config -p --name meohtip3SingleMultipleIndexes2-GPU -i 10000";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        double expected = -2010;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meohtip3SingleMultipleIndexes2-GPU.results").c_str());

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
