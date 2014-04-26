/*!\file
  \brief Class used to read and write configuration files.

  \author Alexander Luchs, Riley Spahn, Seth Wooten, and Orlando Acevedo
  Class used to read and write configuration files.
 
 */

#ifndef CONFIG_SCAN_H
#define CONFIG_SCAN_H

#include <iostream>
#include <fstream>
#include "metroUtil.h"
#include <cstdlib>

using namespace std;

/**
  Class used to read and write configuration files.


*/
class Config_Scan
{
    private:
        /**
          The environment used in the simulation.
        */
        Environment enviro;
        /**
          The path to the configuration file to be read.
        */
        string configpath;
        /**
          The number of step to run the simulation.
        */
        long numOfSteps;
        /**
          The path to the opls file.
        */
        string oplsuaparPath;
        /**
          The path to the Z-matrix file to be used in the simulation.
        */
        string zmatrixPath;
        /**
          The path to the state file to be used.
        */
        string statePath;
        /**
          The path to write the state output files.
        */
        string stateOutputPath;
        /**
          The path where to write the pdb output files.
        */
        string pdbOutputPath;
        /**
          The nonbonded cutoff distance.
        */
        long cutoff;
		
		void throwScanError(string message);

    public:

        /**
          Config file scanner instantiated with a path.
          @param configPath - path to parameters.cfg
        */
        Config_Scan(string configPath);

        /**
          Reads in the config file located at config path
          given in the constructor.
        */
        void readInConfig();

        /**
          @return - returns the environment variable in this Config_Scan 
        */
        Environment * getEnviro();

        /**
          @return - returns the path to the configuration file.
        */
        string getConfigPath();

        /**
          @return getSteps - returns the number of steps to run the simulation.
        */
        long getSteps();
        
        /**
          @return - returns the path to the opls file for the simulation.
        */
        string getOplsusaparPath();

        /**
          @return - returns the path to the z-matrix file to be used in the
          simulation
        */
        string getZmatrixPath();

        /**
          @return - returns the path to the state input file for the simulation.
        */
        string getStatePath();

        /**
          @return - returns the path to the state output file for the simulation.
        */
        string getStateOutputPath();

        /**
          @return - returns the path to the pdb output file(s) for the simulation.
        */
        string getPdbOutputPath();
        /**
          @return getSteps - returns the nonbonded cutoff in the simulation.
        */
        long getcutoff();

        /**
          @return getSteps - returns the nonbonded cutoff in the simulation.
        */        
        unsigned int getrandomseed() {return enviro.randomseed;};

};

#endif //CONFIG_SCAN_H
