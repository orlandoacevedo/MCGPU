/*Intended goal: support read, parse, and extract operations on configuration files to properly initialize 
*  a simulation environment.
*
*Created 19 February 2014. N. Coleman, A. Wallace
*/
//Changes on:
//	Sun, 23 Feb 2014. 1530PM to 1558PM, 1611 to 1655PM, 1757 to 2031PM

#ifndef IOUTILITIES_H
#define IOUTILITIES_H


#include <iostream>
#include <fstream>
#include <sstream>
//#include "../../Utilities/metroUtil.h"
#include "StructLibrary.h" //this is our goal for inclusion, but for now...

#include <cstdlib>

using namespace std;

// this should supplant the need for all the getters in Config_Scan.cpp and .h -Albert
// [we made this decision knowing that it's unsafe; if someone modified these variables willy-nilly, it will be a problem]
struct UtilitiesInfo
{
	Environment currentEnvironment; //The current working environment for the simulation
    string configPath; //The path to the main configuration file read in for the simulation
    unsigned int numOfSteps; //The number of steps to run the simulation
    string oplsuaparPath; //The path to the opls files containing additional geometry data, to be used (eventually) during simulation
    string zmatrixPath; //The path to the Z-matrix files to be used during simulation
    string statePath; //The path to the state information file to be used in the simulation
    string stateOutputPath; //The path where we write the state output files after simulation
    string pdbOutputPath; //The path where we write the pdb output files after simulation
    unsigned int cutoff; //The nonbonded cutoff distance.
    
    UtilitiesInfo() //constructor for the structure
    {
    memset(&currentEnvironment,0,sizeof(Environment)); //just to be safe, zero out this area. (done in original code)
    }
};    

class IOUtilities
{	
	public:
		IOUtilities(string configPath); //this should be the constructor, which does very little on its own
		void readInConfig(); //this should ideally be merged into the constructor, but run separately, everything is okay
		int ReadStateFile(char const* StateFile, Environment * destinationEnvironment, Molecule * destinationMoleculeCollection);
	 	//int ReadStateFile(string StateFile) { return ReadStateFile(StateFile.c_str());}; //preference: do not use this
	 	int WriteStateFile(char const* StateFile, Environment * sourceEnvironment, Molecule * sourceMoleculeCollection); 	
	 	//int WriteStateFile(string StateFile) { return WriteStateFile(StateFile.c_str());}; //preference: do not use this
	 	int writePDB(char const* pdbFile); //this is the old version from SimBox, probably not useful or usable
	 	int writePDB(char const* pdbFile, Environment sourceEnvironment, Molecule * sourceMoleculeCollection); //this is a new version, possibly required for this implementation with IOUtilities
	 	UtilitiesInfo filePathsEtc;
	 	
	 	//unsigned int getrandomseed() {return filePathsEtc.currentEnvironment.randomseed;}; //this random seed isn't so random; it's found within the configuration file.
	 			//as you can guess, the constructor for this class stores that very info in the Environment struct. For sake of safety, we'll leave it there.
	 	
	 private:
		void throwScanError(string message);
		bool readInConfigAlreadyDone;
};
#endif

