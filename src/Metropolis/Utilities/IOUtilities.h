/*Intended goal: support read, parse, and extract operations on configuration files to properly initialize 
*  a simulation environment.
* [1] This file/class will, ideally, replace Config_Scan.cpp & will augment MetroUtil.cpp. (Started, 19 Feb. Almost done, 26 Feb. -Albert)
* [2] This file/class will, ideally, replace Opls_Scan and Zmatrix_Scan. (Started, 26 Feb. -Albert)
*Created 19 February 2014. Authors: Albert Wallace
*/
/*
--Changes made on:
		->Sun, 23 Feb 2014 (Albert)
		->Wed, 26 Feb (Albert)
*/
/*Based on work from earlier sessions by Alexander Luchs, Riley Spahn, Seth Wooten, and Orlando Acevedo*/

#ifndef IOUTILITIES_H
#define IOUTILITIES_H

//_________________________________________________________________________________________________________________
//  INCLUDE statements
//_________________________________________________________________________________________________________________
#include <string>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>

//#include "metroUtil.h"
#include "StructLibrary.h"

//_________________________________________________________________________________________________________________
//  Specific namespace/using requirements
//_________________________________________________________________________________________________________________
		//to avoid "using namespace std;"
using std::vector;
using std::string;
using std::map;

//_________________________________________________________________________________________________________________
//  Custom structures
//_________________________________________________________________________________________________________________

/**
  Structure used to represent the 4 Fourier Coeficients
*/
struct Fourier
{
    double vValues[4];
};

/**
	Structure to supplant the need for all the getters in Config_Scan.cpp and .h
	Stores file paths, as well as the environment being set up using the files at those paths.
	That way, we have one near package to access them all.
*/
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
			
			void deleteOpls_Scan();
			
			//From OPLS_Scan
			/**
			Scans in the opls File calls sub-function addLineToTable
			@param filename - the name/path of the opls file
			@return - success code 
					  0: successful
					  1: error
			*/
			int scanInOpls(string filename); 
		
			/**
			Parses out a line from the opls file and gets (sigma, epsilon, charge)
			adds an entry to the map at the hash number position, 
			calls sub-function checkFormat()
			@param line - a line from the opls file
			@param numOflines - number of lines previously read in, used for error output
			*/
			void addLineToTable(string line, int numOfLines);
		
			/**
			  Checks the format of the line being read in
			  returns false if the format of the line is invalid
			  @param line -  a line from the opls file
			@return - Format code
					  -1: Invalid Format
					  1: Normal OPLS format
					  2: Fourier Coefficent format
			*/
			int checkFormat(string line);
				
			/**
			  Logs all the Errors found in the OPLS file to the output log.
			*/
			void logErrors();
		
			/**
			Returns an Atom struct based on the hashNum (1st col) in Z matrix file
			The Atom struct has -1 for x,y,z and has the hashNum for an id. 
			sigma, epsilon, charges
			@param hashNum -  the hash number (1st col) in Z matrix file
			@return - the atom with that is the value to the hasNum key.
			*/
			Atom getAtom(string hashNum);
		
			/**
			Returns the sigma value based on the hashNum (1st col) in Z matrix file
			@param hashNum -  the hash number (1st col) in Z matrix file
			@return - the sigma value of the atom that is the value associated with the hasNum key.
			*/
			double getSigma(string hashNum);
		
			/**
			Returns the epsilon value based on the hashNum (1st col) in Z matrix file
			@param hashNum - the hash number (1st col) in Z matrix file
			@return - the epsilon value of the atom that is the value associated with the hashNum key.
			*/
			double getEpsilon(string hashNum);
		
			/**
			Returns the charge value based on the hashNum (1st col) in Z matrix file
			@param hashNum -  the hash number (1st col) in Z matrix file
			@return - the charge value of the atom that is the value associated with the hashNum key.
			*/
			double getCharge(string hashNum);
		
			/**
			Returns the V values value based on the hashNum (1st col) in Z matrix file
			@param hashNum -  the hash number (1st col) in Z matrix file
			@return - TODO
			*/
			Fourier getFourier(string hashNum);

	 	
	 private:
		void throwScanError(std::string message);
		bool readInConfigAlreadyDone;
		
		//From OPLS_Scan
		 	/**
			 HashTable that holds all opls references
		   */
			map<string,Atom> oplsTable;
			/**
				HashTable that holds all the Fourier Coefficents
				  stored
			*/
			map<string,Fourier> fourierTable;
			/**
			  the path to the OPLS file.
			*/
			//		//stored in the UtilitiesInfo struct! see also: filePathsEtc.
			 /**
			 Vector used to keep track of errors when lines in 
			  the OPLS file that don't match the correct format.
			*/
			 vector<int> errLinesOPLS;
			 /**
			  Vector used to keep track of errors when the hashes in 
				the OPLS file that exist more than once. 
			*/
			 vector<string> errHashes;
			  /**
			  Vector used to keep track of errors when the hashes in the fourier coefficent
				section of the OPLS file that exist more than once. 
			*/
			 vector<string> errHashesFourier;
};

void writeToLog(std::string text,int stamp);
void writeToLog(std::stringstream& ss, int stamp);

#endif

