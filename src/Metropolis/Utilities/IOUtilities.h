/*Intended goal: support read, parse, and extract operations on configuration files to properly initialize 
*  a simulation environment.
* [1] This file/class will, ideally, replace Config_Scan.cpp & will augment MetroUtil.cpp. (Started, 19 Feb. Beta completion, 28 Feb. -Albert)
* [2] This file/class will, ideally, replace Opls_Scan and Zmatrix_Scan. (Started, 26 Feb. -Albert)
*Created 19 February 2014. Authors: Albert Wallace
*/
/*
--Changes made on:
		->Sun, 23 Feb 2014 (Albert)
		->Wed, 26 Feb (Albert)
		->Thu, 27 Feb (Albert, then Tavis)
		->Fri, 28 Feb (Albert)
		->Mon, 03 Mar; Wed, 05 Mar; Thur, 06 Mar, Fri, 07 Mar; Monday, 10 Mar (Albert)
*/
/*Based on work from earlier sessions by Alexander Luchs, Riley Spahn, Seth Wooten, and Orlando Acevedo*/

#ifndef IOUTILITIES_H
#define IOUTILITIES_H

//_________________________________________________________________________________________________________________
//  INCLUDE statements
//_________________________________________________________________________________________________________________
#include <assert.h>
#include <errno.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <vector>

#include "StructLibrary.h"
//#include "../../Utilities/geometricUtil.h"
#include "MathLibrary.h"

//_________________________________________________________________________________________________________________
//  Specific namespace/using requirements
//_________________________________________________________________________________________________________________
		//to avoid "using namespace std;"
using std::vector;
using std::string;
using std::map;  

//_________________________________________________________________________________________________________________
//Actual IOUtilities class definition; handles most of the setup duties when it comes to parsing configuration files and the like.
//_________________________________________________________________________________________________________________
class IOUtilities
{	
	public:
			IOUtilities(std::string configPath); //this should be the constructor, which does very little on its own.
				///As of right now, only calls the readInConfig() method, and does nothing more. So only file paths are acquired,
				///   but nothing else is set up for further execution
			bool readInConfig(); //this represents the first of the chain of calls to configuration methods, called from
						// the constructor. (Does *not* call the second in the chain, or in other words does not continue environment setup.)
			int ReadStateFile(char const* StateFile, Environment * destinationEnvironment, Molecule * destinationMoleculeCollection);
			int WriteStateFile(char const* StateFile, Environment * sourceEnvironment, Molecule * sourceMoleculeCollection); 	
			int writePDB(char const* pdbFile); //this is the old version from SimBox, probably not useful or usable
			int writePDB(char const* pdbFile, Environment sourceEnvironment, Molecule * sourceMoleculeCollection); //this is a new version, possibly required for this implementation with IOUtilities

			
			
			//From OPLS_Scan......
			/*
			What used to be the old destructor for the OPLS_Scan object.
			*/
			void destructOpls_Scan();
			
			/**
			Scans in the opls File calls sub-function addLineToTable
			@return - success code 
					  0: successful
					  1: error
			*/
			int scanInOpls(); 
		
			/**
			Parses out a line from the opls file and gets (sigma, epsilon, charge)
			adds an entry to the map at the hash number position, 
			calls sub-function checkFormat()
			@param line - a line from the opls file
			@param numOflines - number of lines previously read in, used for error output
			*/
			void addLineToTable(std::string line, int numOfLines);
		
			/**
			  Checks the format of the line being read in
			  returns false if the format of the line is invalid
			  @param line -  a line from the opls file
			  @return - Format code
					  -1: Invalid Format
					  1: Normal OPLS format
					  2: Fourier Coefficent format
			*/
			int OPLScheckFormat(string line);
				
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
	 	
	
		void throwScanError(std::string message);
		bool readInConfigAlreadyDone;
		bool criticalErrorEncountered;
		
		//From OPLS_Scan
			/**
			OPLS scanning is required to assign sigma, epsilon and charge values.
		  	*/
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
		//From Zmatrix_Scan
			void destructZmatrix_Scan();
		
			/**
			  Scans in the z-matrix File calls sub-function parseLine
			  @param filename - the name/path of the z-matrix file
			  @return - success code
						0: Valid z-matrix path
						1: Invalid z-matrix path
			*/
			int scanInZmatrix(); 
		
			/**
			  Parses out a line from the zmatrix file and gets the atom from the OPLS hash
			  Creates structs for bond, angle, dihedral, if applicable
			  calls sub-function checkFormat()
			  @param line - a line from the zmatrix file
			  @param numOflines - number of lines previously read in, used for error output
			*/
			void parseLine(string line, int numOfLines);
		
			/**
			  Checks the format of the line being read in
			  returns false if the format of the line is invalid
			  @param line -  a line from the zmatrix file
			  @param stringCount - number of strings in a line you're looking at
			  @return - Format code
						1: Base atom row listing
						2: TERZ line
						3: Geometry Variations
						4: Variable Bonds
						5: Additional Bonds
						6: Harmonic Constraints
						7: Variable Angles
						8: Additional Angles
						9: Variable Dihedrals
						10: Additional Dihedrals
						11: Domain Definitions
						-2: Final Blank Line
						-1: Invalid line format
			*/
			int ZMcheckFormat(string line);

			/**
			  Handles the additional stuff listed at the bottom of the Z-matrix file
			  @param line -   a line from the zmatrix file
			  @param cmdFormat- the int representing the format for which the line applies :see checkFormat
			*/
			void handleZAdditions(string line, int cmdFormat);

			/**
			  Creates a vector containg the Hop distances of a molecule
			  for all hops that have a distance greater than 3.
			  @param molec - the molecule to check its bonds for valid hops
			  @return - returns a vector of hops needed for energy calculations
			*/
			vector<Hop> calculateHops(Molecule molec);

			/**
			  Checks if the int item is contained in the vector
			  @param vect - the vector to search through
			  @param item - the item to search for
			  @return - true if the int is in the vector and false otherwise.
			*/
			bool contains(vector<int> &vect, int item);

			/**
			  Returns the distance of the shortest path amongst bonds between two atoms
			  @param atom1 - the id of the starting atom
			  @param atom2 - the if of the ending atom
			  @param size - number of atoms in molecule
			  @param graph - a 2d array representing all the bonds in the molecule
			  @return - returns the number of bonds seperating atom1 and atom2.
			*/
			int findHopDistance(int atom1,int atom2,int size, int **graph);

			/**
			  Creates a 2d array representing all the bonds in the molecule
			  the array is filled with 1 or 0 for a bond or no bond between atoms
			  array is n*n where n is the number of atoms.		   
			  @param graph - a 2d array representing all the bonds in the molecule
			  @param molec - the molecule to check its bonds for valid hops
			*/
			void buildAdjacencyMatrix(int **&graph, Molecule molec);

		  
			/**
			  Creates a molecule(s)  based on a starting unique ID and the pattern specified
			  by the Z-matrix in the scan functions
			  @param startingID - first ID for the molecule being built
			  @return - vector of unique molecules that are copies of the molecules from the
			  Z-matrix file.
			*/
			vector<Molecule> buildMolecule(int startingID);
			//The path to the Z-matrix file
		  			//Stored in the UtilitiesInfo structure. See also: filePathsEtc.
	
		  /**
			Vector that holds example molecules.
		  */
		  vector<Molecule> moleculePattern_ZM;
		  /**
			Vector that holds the atoms contained in the molecule in the Z-matrix.
		  */
		  vector<Atom> atomVector_ZM;
		  /**
			Vector that holds the bonds contained in the molecule in the Z-matrix.
		  */
		  vector<Bond> bondVector_ZM;
		  /**
			Vector that holds the angles contained in the molecule in the Z-matrix.
		  */
		  vector<Angle> angleVector_ZM;
		  /**
			Vector that holds the dihedrals contained in the molecule in the Z-matrix.
		  */
		  vector<Dihedral> dihedralVector_ZM;
		  /**
			Vector of dummy atoms held seperately from the normal atoms.
		  */
			//vector<unsigned long> dummies_ZM; //potentially unused
		  /**
			Global variable that determines if there are multiple molecules in a Z-matrix.
		  */
		  bool startNewMolecule_ZM;
		  /**
			Global variable used to store the format id of the last line scanned.
			  needed to read the bottom sections of the Z-matrix file.
		  */
		  int previousFormat_ZM;
		  
		  
		  /******************
		  *Path variables and the temporary environment
		  ***************************************/
		  	Environment * currentEnvironment; //The current working environment for the simulation
    		std::string configPath; //The path to the main configuration file read in for the simulation
    		unsigned int numOfSteps; //The number of steps to run the simulation
    		std::string oplsuaparPath; //The path to the opls files containing additional geometry data, to be used (eventually) during simulation
    		std::string zmatrixPath; //The path to the Z-matrix files to be used during simulation
    		std::string statePath; //The path to the state information file to be used in the simulation
    		std::string stateOutputPath; //The path where we write the state output files after simulation
    		std::string pdbOutputPath; //The path where we write the pdb output files after simulation
    		unsigned int cutoff; //The nonbonded cutoff distance.
		  
};

void writeToLog(std::string text,int stamp);
void writeToLog(std::stringstream& ss, int stamp);

#endif

