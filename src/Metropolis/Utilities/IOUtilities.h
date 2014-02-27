/*Intended goal: support read, parse, and extract operations on configuration files to properly initialize 
*  a simulation environment.
* [1] This file/class will, ideally, replace Config_Scan.cpp & will augment MetroUtil.cpp. (Started, 19 Feb. Almost done, 26 Feb. -Albert)
* [2] This file/class will, ideally, replace Opls_Scan and Zmatrix_Scan. (Started, 26 Feb. -Albert)
*Created 19 February 2014. Authors: Nathan Coleman, Albert Wallace
*/
/*
--Changes made on:
		->Sun, 23 Feb 2014 (Albert)
		->Wed, 26 Feb (Albert)
*/
/*Based on work from earlier sessions by Alexander Luchs, Riley Spahn, Seth Wooten, and Orlando Acevedo*/

#ifndef IOUTILITIES_H
#define IOUTILITIES_H

#include <string>
#include "StructLibrary.h" //this is our goal for inclusion, but for now...

//using namespace std;

//_________________________________________________________________________________________________________________
//  UtilitiesInfo structure
//_________________________________________________________________________________________________________________
// this should supplant the need for all the getters in Config_Scan.cpp and .h -Albert
// [we made this decision knowing that it's unsafe; if someone modified these variables willy-nilly, it will be a problem]
struct UtilitiesInfo
{
	Environment * currentEnvironment; //The current working environment for the simulation
    std::string configPath; //The path to the main configuration file read in for the simulation
    unsigned int numOfSteps; //The number of steps to run the simulation
    std::string oplsuaparPath; //The path to the opls files containing additional geometry data, to be used (eventually) during simulation
    std::string zmatrixPath; //The path to the Z-matrix files to be used during simulation
    std::string statePath; //The path to the state information file to be used in the simulation
    std::string stateOutputPath; //The path where we write the state output files after simulation
    std::string pdbOutputPath; //The path where we write the pdb output files after simulation
    unsigned int cutoff; //The nonbonded cutoff distance.
    
    UtilitiesInfo() //constructor for the structure
    {
    		// memset(&currentEnvironment,0,sizeof(Environment)); //just to be safe, zero out this area. (done in original code) 
    				//commented out in favor of letting Environment's constructor do this dynamically
    				//to revert, ensure you remove the asterisk from "Environment * currentEnvironment" above
    currentEnvironment = new Environment(); //dynamically allocate space for the Environment
    configPath = "";
    numOfSteps = 0;
    oplsuaparPath = "";
    zmatrixPath = "";
    statePath = "";
    stateOutputPath = "";
    pdbOutputPath = "";
    cutoff = 0;
    }
};    

//_________________________________________________________________________________________________________________
//Actual IOUtilities class definition; handles most of the setup duties when it comes to parsing configuration files and the like.
//_________________________________________________________________________________________________________________
class IOUtilities
{	
	public:
		IOUtilities(std::string configPath); //this should be the constructor, which does very little on its own
		void readInConfig(); //this should ideally be merged into the constructor, but run separately, everything is okay
		int ReadStateFile(char const* StateFile, Environment * destinationEnvironment, Molecule * destinationMoleculeCollection);
	 	int WriteStateFile(char const* StateFile, Environment * sourceEnvironment, Molecule * sourceMoleculeCollection); 	
	 	int writePDB(char const* pdbFile); //this is the old version from SimBox, probably not useful or usable
	 	int writePDB(char const* pdbFile, Environment sourceEnvironment, Molecule * sourceMoleculeCollection); //this is a new version, possibly required for this implementation with IOUtilities
	 	UtilitiesInfo * filePathsEtc;
	 	
	 	//unsigned int getrandomseed() {return filePathsEtc.currentEnvironment.randomseed;}; //this random seed isn't so random; it's found within the configuration file.
	 			//as you can guess, the constructor for this class stores that very info in the Environment struct. For sake of safety, we'll leave it there.
	 	
	 private:
		void throwScanError(std::string message);
		bool readInConfigAlreadyDone;
};

void writeToLog(std::string text,int stamp);
void writeToLog(std::stringstream& ss, int stamp);

#endif

