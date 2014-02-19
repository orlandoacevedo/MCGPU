/*Intended goal: support read, parse, and extract operations on configuration files to properly initialize 
*  a simulation environment.
*
*Created 19 February 2014. N. Coleman, A. Wallace
*/

#ifndef IOUTILITIES_H
#define IOUTILITIES_H


#include <iostream>
#include <fstream>
#include "../Utilities/metroUtil.h"
#include <cstdlib>

using namespace std;

// this should supplant the need for all the getters in Config_Scan.cpp and .h -Albert
// [we made this decision knowing that it's unsafe; if someone modified these variables willy-nilly, it will be a problem]
Struct UtilitiesInfo
{
	Environment currentEnvironment; //The current working environment for the 
    string configPath; //The path to the main configuration file read in for the simulation
    unsigned int numOfSteps; //The number of steps to run the simulation
    string oplsuaparPath; //The path to the opls files containing additional geometry data, to be used (eventually) during simulation
    string zmatrixPath; //The path to the Z-matrix files to be used during simulation
    string statePath; //The path to the state information file to be used in the simulation
    string stateOutputPath; //The path where we write the state output files after simulation
    string pdbOutputPath; //The path where we write the pdb output files after simulation
    unsigned int cutoff; //The nonbonded cutoff distance.
};    

class IOUtilities
{	
	public:
		Config_Scan(string configPath);
		void readInConfig(); //this should ideally be merged into the constructor
		int ReadStateFile(char const* StateFile);
	 	//int ReadStateFile(string StateFile) { return ReadStateFile(StateFile.c_str());}; //preference: do not use this
	 	int WriteStateFile(char const* StateFile); 	
	 	//int WriteStateFile(string StateFile) { return WriteStateFile(StateFile.c_str());};
	 	int writePDB(char const* pdbFile);
	 	
	 private:
		void throwScanError(string message);
};
#endif

