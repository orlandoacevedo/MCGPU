#include "Applications/Application.h"
#include "gtest/gtest.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

//Test indole with 4000 molecules and a  primary index of 1 on CPU
TEST (Indole4000Test, OnePrimaryIndex)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/indole4000Tests/indole4000-1MPI.config"));
        configFile.open(configFilePath.c_str());
        std::stringstream cfg;
        cfg << ""
                << "#size of periodic box (x, y, z in angstroms)\n"
                << "87.17\n"
                << "87.17\n"
                << "87.17\n"
                << "#temperature in Kelvin\n"
                << "298.15\n"
                << "#max translation\n"
                << ".12\n"
                << "#number of steps\n"
                << "100000\n"
                << "#number of molecules\n"
                << "4000\n"
                << "#path to opls.par file\n"
                << MCGPU << "resources/bossFiles/oplsaa.par\n"
                << "#path to z matrix file\n"
		<< MCGPU << "test/unittests/MultipleSolvents/indole4000Tests/indole.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
                << "#cutoff distance in angstroms\n"
                << "11.0\n"
                << "#max rotation\n"
                << "6.0\n"
                << "#Random Seed Input\n"
                << "12345\n"
                << "#Primary Atom Index\n"
                << "1";
        configFile << cfg.str();
        configFile.close();

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/indole4000Tests/indole4000-1MPI.config -s --threads 2 --name indole4000-1MPI -i 10000";

        //launch MCGPU in serial
         system(ss.str().c_str());
        double expected = 480016;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/indole4000-1MPI.results").c_str());

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

//Test indole with 4000 molecules and a  primary index of 1 on GPU
TEST (Indole4000Test, OnePrimaryIndex_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/indole4000Tests/indole4000-1MPI.config -p --threads 12 --name indole4000-1MPI-GPU -i 10000";

        //launch MCGPU in parallel
         system(ss.str().c_str());
        double expected = 480016;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/indole4000-1MPI-GPU.results").c_str());

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


//Test indole with 4000 molecules and a  primary index of 1 on CPU
//Uses neighbor list function
TEST (Indole4000Test, NeighborListFunction1MPI)
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
                << MCGPU << "/test/unittests/MultipleSolvents/indole4000Tests/indole4000-1MPI.config -s --threads 12 --name indole4000-1MPI-NL -i 10000 -l";

        //launch MCGPU in parallel
         system(ss.str().c_str());
        double expected = 500239;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/indole4000-1MPI-NL.results").c_str());

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

//Test indole with 4000 molecules and a  primary index of 1 on GPU
//Uses neighbor list function
TEST (Indole4000Test, NeighborListFunction1MPI_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/indole4000Tests/indole4000-1MPI.config -p --threads 12 --name indole4000-1MPI-NL-GPU -i 10000 -l";

        //launch MCGPU in parallel
         system(ss.str().c_str());
        double expected = 500239;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/indole4000-1MPI-NL-GPU.results").c_str());

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

//Test indole with 4000 molecules and a primary index of [1,2] on CPU
TEST (Indole4000Test, TwoPrimaryIndex)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/indole4000Tests/indole4000-2MPI.config"));
        configFile.open(configFilePath.c_str());
        std::stringstream cfg;
        cfg << ""
                << "#size of periodic box (x, y, z in angstroms)\n"
                << "87.17\n"
                << "87.17\n"
                << "87.17\n"
                << "#temperature in Kelvin\n"
                << "298.15\n"
                << "#max translation\n"
                << ".12\n"
                << "#number of steps\n"
                << "100000\n"
                << "#number of molecules\n"
                << "4000\n"
                << "#path to opls.par file\n"
                << MCGPU << "resources/bossFiles/oplsaa.par\n"
                << "#path to z matrix file\n"
                 << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests/indole.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
                << "#path to state output\n"
		<< MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
                << "#cutoff distance in angstroms\n"
                << "11.0\n"
                << "#max rotation\n"
                << "6.0\n"
                << "#Random Seed Input\n"
                << "12345\n"
                << "#Primary Atom Index\n"
                << "[1,2]";
        configFile << cfg.str();
        configFile.close();

        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/indole4000Tests/indole4000-2MPI.config -s --threads 12 --name indole4000-2MPI -i 10000";

        //launch MCGPU in serial
         system(ss.str().c_str());
        double expected = 449963;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/indole4000-2MPI.results").c_str());

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

//Test indole with 4000 molecules and a primary index of [1,2] on GPU
TEST (Indole4000Test, TwoPrimaryIndex_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/indole4000Tests/indole4000-2MPI.config -p --threads 12 --name indole4000-2MPI-GPU -i 10000";

        //launch MCGPU in parallel
         system(ss.str().c_str());
        double expected = 449963;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/indole4000-2MPI-GPU.results").c_str());

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


//Test indole with 4000 molecules and a primary index of [1,2] on CPU
//Uses neighborlist function
TEST (Indole4000Test, NeighborListFunction2MPI)
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
                << MCGPU << "/test/unittests/MultipleSolvents/indole4000Tests/indole4000-2MPI.config -s --threads 12 --name indole4000-2MPI-NL -i 10000 -l";

        //launch MCGPU in serial
         system(ss.str().c_str());
        double expected = 499427;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/indole4000-2MPI-NL.results").c_str());

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

//Test indole with 4000 molecules and a primary index of [1,2] on GPU
//Uses neighborlist function
TEST (Indole4000Test, NeighborListFunction2MPI_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/indole4000Tests/indole4000-2MPI.config -p --threads 12 --name indole4000-2MPI-NL-GPU -i 10000 -l";

        //launch MCGPU in parallel
         system(ss.str().c_str());
        double expected = 499427;
        double energyResult = -1;

        std::ifstream infile(std::string(MCGPU + "bin/indole4000-2MPI-NL-GPU.results").c_str());

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


//Test indole with 4000 molecules and multiple solvent primary indexes of 1,2 on CPU
TEST (Indole4000Test, MultipleSolventsDefinition)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/indole4000Tests/indole4000MulSolvent.config"));
        configFile.open(configFilePath.c_str());
        std::stringstream cfg;
        cfg << ""
                << "#size of periodic box (x, y, z in angstroms)\n"
                << "87.17\n"
                << "87.17\n"
                << "87.17\n"
                << "#temperature in Kelvin\n"
                << "298.15\n"
                << "#max translation\n"
                << ".12\n"
                << "#number of steps\n"
                << "100000\n"
                << "#number of molecules\n"
                << "4000\n"
                << "#path to opls.par file\n"
                << MCGPU << "resources/bossFiles/oplsaa.par\n"
                << "#path to z matrix file\n"
                 << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests/indole.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
                << "#cutoff distance in angrstoms\n"
                << "11.0\n"
                << "#max rotation\n"
                << "6.0\n"
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
                << MCGPU << "/test/unittests/MultipleSolvents/indole4000Tests/indole4000MulSolvent.config -s -i 10000 >"
                << " "
                << MCGPU << "/bin/indole4000MulSolvent.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/indole4000MulSolvent.txt").c_str());

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
        

//Test indole with 4000 molecules and multiple solvent primary indexes of 1,2 on GPU
TEST (Indole4000Test, MultipleSolventsDefinition_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/indole4000Tests/indole4000MulSolvent.config -p -i 10000 >"
                << " "
                << MCGPU << "/bin/indole4000MulSolvent-GPU.txt 2>&1";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/indole4000MulSolvent-GPU.txt").c_str());

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

//Test indole with 4000 molecules and multiple solvent primary indexes of [1,2],[3,4] on CPU
TEST (Indole4000Test, MultipleSolventsDefinitionMPI)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/indole4000Tests/indole4000MulSolvent-MPI.config"));
        configFile.open(configFilePath.c_str());
        std::stringstream cfg;
        cfg << ""
                << "#size of periodic box (x, y, z in angstroms)\n"
                << "87.17\n"
                << "87.17\n"
                << "87.17\n"
                << "#temperature in Kelvin\n"
                << "298.15\n"
                << "#max translation\n"
                << ".12\n"
                << "#number of steps\n"
                << "100000\n"
                << "#number of molecules\n"
                << "4000\n"
                << "#path to opls.par file\n"
                << MCGPU << "resources/bossFiles/oplsaa.par\n"
                << "#path to z matrix file\n"
                 << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests/indole.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
		<< "#cutoff distance in angrstoms\n"
                << "11.0\n"
                << "#max rotation\n"
                << "6.0\n"
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
                << MCGPU << "/test/unittests/MultipleSolvents/indole4000Tests/indole4000MulSolvent-MPI.config -s -i 10000 >"
                << " "
                << MCGPU << "/bin/indole4000MulSolvent-MPI.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/indole4000MulSolvent-MPI.txt").c_str());

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

//Test indole with 4000 molecules and multiple solvent primary indexes of [1,2],[3,4] on GPU
TEST (Indole4000Test, MultipleSolventsDefinitionMPI_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/indole4000Tests/indole4000MulSolvent-MPI.config -p -i 10000 >"
                << " "
                << MCGPU << "/bin/indole4000MulSolvent-MPI-GPU.txt 2>&1";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/indole4000MulSolvent-MPI-GPU.txt").c_str());

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

//Test indole with 4000 molecules and multiple solvent primary indexes of 1,[1,2] on CPU
TEST (Indole4000Test, SingleMultipleIndexes)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/indole4000Tests/indole4000SingleMultipleIndexes.config"));
        configFile.open(configFilePath.c_str());
        std::stringstream cfg;
        cfg << ""
                << "#size of periodic box (x, y, z in angstroms)\n"
                << "87.17\n"
                << "87.17\n"
                << "87.17\n"
                << "#temperature in Kelvin\n"
                << "298.15\n"
                << "#max translation\n"
                << ".12\n"
                << "#number of steps\n"
                << "100000\n"
                << "#number of molecules\n"
                << "4000\n"
                << "#path to opls.par file\n"
                << MCGPU << "resources/bossFiles/oplsaa.par\n"
                << "#path to z matrix file\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests/indole.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
		<< "#cutoff distance in angrstoms\n"
                << "11.0\n"
                << "#max rotation\n"
                << "6.0\n"
                << "#Random Seed Input\n"
                << "12345\n"
                << "#Primary Atom Index\n"
                << "1,[1,2]";
        configFile << cfg.str();
        configFile.close();

        //Since this case results in an error it does not print to a '.results' file.
        //So the commandline error output is pipelined to a text file
        //The textfile is used to compare the expected and actual output
        std::stringstream ss;
                ss << MCGPU << "/bin/metrosim "
                << " "
                << MCGPU << "/test/unittests/MultipleSolvents/indole4000Tests/indole4000SingleMultipleIndexes.config -s -i 10000 >"
                << " "
                << MCGPU << "/bin/indole4000SingleMultipleIndexes.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/indole4000SingleMultipleIndexes.txt").c_str());

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

//Test indole with 4000 molecules and multiple solvent primary indexes of 1,[1,2] on GPU
TEST (Indole4000Test, SingleMultipleIndexes_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/indole4000Tests/indole4000SingleMultipleIndexes.config -p -i 10000 >"
                << " "
                << MCGPU << "/bin/indole4000SingleMultipleIndexes-GPU.txt 2>&1";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/indole4000SingleMultipleIndexes-GPU.txt").c_str());

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

//Test indole with 4000 molecules and multiple solvent primary indexes of [1,2],1 on CPU
TEST (Indole4000Test, SingleMultipleIndexes2)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/indole4000Tests/indole4000SingleMultipleIndexes2.config"));
        configFile.open(configFilePath.c_str());
        std::stringstream cfg;
        cfg << ""
                << "#size of periodic box (x, y, z in angstroms)\n"
                << "87.17\n"
                << "87.17\n"
                << "87.17\n"
                << "#temperature in Kelvin\n"
                << "298.15\n"
                << "#max translation\n"
                << ".12\n"
                << "#number of steps\n"
                << "100000\n"
                << "#number of molecules\n"
                << "4000\n"
                << "#path to opls.par file\n"
                << MCGPU << "resources/bossFiles/oplsaa.par\n"
                << "#path to z matrix file\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests/indole.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/indole4000Tests\n"
                << "#cutoff distance in angrstoms\n"
		<< "11.0\n"
                << "#max rotation\n"
                << "6.0\n"
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
                << MCGPU << "/test/unittests/MultipleSolvents/indole4000Tests/indole4000SingleMultipleIndexes2.config -s -i 10000 >"
                << " "
                << MCGPU << "/bin/indole4000SingleMultipleIndexes2.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/indole4000SingleMultipleIndexes2.txt").c_str());

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

//Test indole with 4000 molecules and multiple solvent primary indexes of [1,2],1 on GPU
TEST (Indole4000Test, SingleMultipleIndexes2_GPU)
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
                << MCGPU << "/test/unittests/MultipleSolvents/indole4000Tests/indole4000SingleMultipleIndexes2.config -p -i 10000 >"
                << " "
                << MCGPU << "/bin/indole4000SingleMultipleIndexes2-GPU.txt 2>&1";

        //launch MCGPU in parallel
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/indole4000SingleMultipleIndexes2-GPU.txt").c_str());

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
