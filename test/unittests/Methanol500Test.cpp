#include "Applications/Application.h"
#include "gtest/gtest.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

//Test methanol with 500 molecules and a  primary index of 1 on CPU
TEST (Meoh500Test, OnePrimaryIndex)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/Methanol500Tests/meoh5001MPI.config"));
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
                << "1";
        configFile << cfg.str();
        configFile.close();

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/Methanol500Tests/meoh5001MPI.config -s --threads 12 --name meoh5001MPI -i 10000";

        //launch MCGPU in serial
         system(ss.str().c_str());
        double expected = -3454;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meoh5001MPI.results").c_str());

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

//Test methanol with 500 molecules and a  primary index of 1 on GPU
TEST (Meoh500Test, OnePrimaryIndex_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/Methanol500Tests/meoh5001MPI.config -p --threads 12 --name meoh5001MPI-GPU -i 10000";

        //launch MCGPU in parallel
         system(ss.str().c_str());
        double expected = -3454;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meoh5001MPI-GPU.results").c_str());

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

//Test methanol with 500 molecules and a  primary index of 1 on CPU
//Using neighborlist
TEST (Meoh500Test, NeighborListFunction1MPI)
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
                << MCGPU << "/test/unittests/MultipleSolvents/Methanol500Tests/meoh5001MPI.config -s --threads 12 --name meoh5001MPI-NL -i 10000 -l";

        //launch MCGPU in serial
         system(ss.str().c_str());
        double expected = -3402;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meoh5001MPI-NL.results").c_str());

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

//Test methanol with 500 molecules and a  primary index of 1 on GPU
//Using neighborlist
TEST (Meoh500Test, NeighborListFunction1MPI_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/Methanol500Tests/meoh5001MPI.config -p --threads 12 --name meoh5001MPI-NL-GPU -i 10000 -l";

        //launch MCGPU in parallel
         system(ss.str().c_str());
        double expected = -3402;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meoh5001MPI-NL-GPU.results").c_str());

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

//Test methanol with 500 molecules and two  primary indexes [1,2] on CPU
TEST (Meoh500Test, TwoPrimaryIndex)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/Methanol500Tests/meoh5002MPI.config"));
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

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/Methanol500Tests/meoh5002MPI.config -s --threads 12 --name meoh5002MPI -i 10000";

        //launch MCGPU in serial
         system(ss.str().c_str());
        double expected = -3511;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meoh5002MPI.results").c_str());

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

//Test methanol with 500 molecules and two  primary indexes [1,2] on GPU
TEST (Meoh500Test, TwoPrimaryIndex_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/Methanol500Tests/meoh5002MPI.config -p --threads 12 --name meoh5002MPI-GPU -i 10000";

        //launch MCGPU in parallel
         system(ss.str().c_str());
        double expected = -3511;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meoh5002MPI-GPU.results").c_str());

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

//Test methanol with 500 molecules and primary index of [1,2] on CPU
//Using neighborlist
TEST (Meoh500Test, NeighborListFunction2MPI)
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
                << MCGPU << "/test/unittests/MultipleSolvents/Methanol500Tests/meoh5002MPI.config -s --threads 12 --name meoh5002MPI-NL -i 10000 -l";

        //launch MCGPU in serial
         system(ss.str().c_str());
        double expected = -3433;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meoh5002MPI-NL.results").c_str());

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

//Test methanol with 500 molecules and primary index of [1,2] on GPU
//Using neighborlist
TEST (Meoh500Test, NeighborListFunction2MPI_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/Methanol500Tests/meoh5002MPI.config -p --threads 12 --name meoh5002MPI-NL-GPU -i 10000 -l";

        //launch MCGPU in parallel
         system(ss.str().c_str());
        double expected = -3433;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/meoh5002MPI-NL-GPU.results").c_str());

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

//Test methanol with 500 molecules and multiple solvent primary indexes 1,2  on CPU
TEST (Meoh500Test, MultipleSolventDefinition)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/Methanol500Tests/meoh500MulSolvent.config"));
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
                << "1,2";
        configFile << cfg.str();
        configFile.close();

	//Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/Methanol500Tests/meoh500MulSolvent.config -s -i 10000 >"
                << " "
                << MCGPU << "/bin/meoh500MulSolvent.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/meoh500MulSolvent.txt").c_str());

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

//Test methanol with 500 molecules and multiple solvent primary indexes 1,2  on GPU
TEST (Meoh500Test, MultipleSolventDefinition_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/Methanol500Tests/meoh500MulSolvent.config -p -i 10000 >"
                << " "
                << MCGPU << "/bin/meoh500MulSolvent-GPU.txt 2>&1";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/meoh500MulSolvent-GPU.txt").c_str());

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

//Test methanol with 500 molecules and multiple solvent primary indexes [1,2],[3,4]  on CPU
TEST (Meoh500Test, MultipleSolventDefinitionMPI)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/Methanol500Tests/meoh500MulSolventMPI.config"));
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
                << "[1,2],[3,4]";
        configFile << cfg.str();
        configFile.close();

        //Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/Methanol500Tests/meoh500MulSolventMPI.config -s -i 10000 >"
                << " "
                << MCGPU << "/bin/meoh500MulSolventMPI.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/meoh500MulSolventMPI.txt").c_str());

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

//Test methanol with 500 molecules and multiple solvent primary indexes [1,2],[3,4] on GPU
TEST (Meoh500Test, MultipleSolventDefinitionMPI_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/Methanol500Tests/meoh500MulSolventMPI.config -s -i 10000 >"
                << " "
                << MCGPU << "/bin/meoh500MulSolventMPI.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/meoh500MulSolventMPI.txt").c_str());

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

//Test methanol with 500 molecules and multiple solvent primary indexes 1, [1,2] on CPU
TEST (Meoh500Test, SingleMultipleIndexes)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/Methanol500Tests/meoh500SingleMultipleIndexes.config"));
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
                << "1, [1,2]";
        configFile << cfg.str();
        configFile.close();

        //Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/Methanol500Tests/meoh500SingleMultipleIndexes.config -s -i 10000 >"
                << " "
                << MCGPU << "/bin/meoh500SingleMultipleIndexes.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/meoh500SingleMultipleIndexes.txt").c_str());

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

//Test methanol with 500 molecules and multiple solvent primary indexes 1, [1,2] on GPU
TEST (Meoh500Test, SingleMultipleIndexes_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/Methanol500Tests/meoh500SingleMultipleIndexes.config -p -i 10000 >"
                << " "
                << MCGPU << "/bin/meoh500SingleMultipleIndexes-GPU.txt 2>&1";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/meoh500SingleMultipleIndexes-GPU.txt").c_str());

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

//Test methanol with 500 molecules and multiple solvent primary indexes [1,2],1 on CPU
TEST (Meoh500Test, SingleMultipleIndexes2)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/Methanol500Tests/meoh500SingleMultipleIndexes2.config"));
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
                << "[1,2],1";
        configFile << cfg.str();
        configFile.close();

        //Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/Methanol500Tests/meoh500SingleMultipleIndexes2.config -s -i 10000 >"
                << " "
                << MCGPU << "/bin/meoh500SingleMultipleIndexes2.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/meoh500SingleMultipleIndexes2.txt").c_str());

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

//Test methanol with 500 molecules and multiple solvent primary indexes [1,2],1 on GPU
TEST (Meoh500Test, SingleMultipleIndexes2_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/Methanol500Tests/meoh500SingleMultipleIndexes2.config -p -i 10000 >"
                << " "
                << MCGPU << "/bin/meoh500SingleMultipleIndexes2-GPU.txt 2>&1";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/meoh500SingleMultipleIndexes2-GPU.txt").c_str());

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
