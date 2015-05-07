#include "Applications/Application.h"
#include "gtest/gtest.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

//Test adesuch with 1 primary index on CPU
//Should result in an error since t3pdim has two molecules
TEST (adesuchTest, OnePrimaryIndex)
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
        std::string configFilePath (std::string (MCGPU + "/test/unittests/MultipleSolvents/adesuchTests/adesuch1MPI.config"));
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
                << MCGPU << "test/unittests/MultipleSolvents/adesuchTests/adesuch.z\n"
                << "#path to state input\n"
                << MCGPU << "test/unittests/MultipleSolvents/adesuchTests\n"
                << "#path to state output\n"
                << MCGPU << "test/unittests/MultipleSolvents/adesuchTests\n"
                << "#pdb output path\n"
                << MCGPU << "test/unittests/MultipleSolvents/adesuchTests\n"
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
                << MCGPU << "/test/unittests/MultipleSolvents/adesuchTests/adesuch1MPI.config -s --verbose --threads 12 -i 10000 >"
                << " "
                << MCGPU << "/bin/adesuch1MPI.txt 2>&1";

        //launch MCGPU in serial
        system(ss.str().c_str());
        std::string errorResult;
        std::ifstream infile(std::string(MCGPU + "bin/adesuch1MPI.txt").c_str());

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
