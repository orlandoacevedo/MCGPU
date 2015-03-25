#ifndef FILE_UTILITIES_H
#define FILE_UTILITIES_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <map>
#include <sstream>
#include "Metropolis/SimulationArgs.h"
#include "StructLibrary.h"
#include "MathLibrary.h"
#include "Metropolis/Box.h"

#define DEFAULT 0
#define START 1
#define END 2
#define OPLS 3
#define Z_MATRIX 4
#define GEOM 5

#define TRIMMED_CHARS " \n\r\t"
#define COMMENT_DELIM ';'

class ConfigScanner
{
    private:
        /**
          The environment used in the simulation.
        */
        Environment enviro;
        /**
          The path to the configuration file to be read.
        */
        string configPath;
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
        
        bool isSafeToContinue; //used to help with logic flow during detection of Zmatrix and State file paths in file
        bool useStatefileSetup; //if true, always try the state file for pre-sim setup; else, try ZMatrix
        bool useZMatrixSetup; //if set to true, use the ZMatrix file unless state file is available
		
	void throwScanError(string message);
	void parsePrimaryIndexDefinitions(string definitions);
    public:
        /// Default constructor that initializes the scanner to default
        /// values but does not read in any files yet.
        ConfigScanner();

        /**
          Reads in the config file located at config path
          given in the constructor.
        */
        bool readInConfig(string configpath);

        /**
          @return - returns the environment variable in this Config_Scan 
        */
        Environment* getEnviro();

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
          @return - returns the random seed used in the simulation.
        */        
        unsigned int getrandomseed() {return enviro.randomseed;}
        
        /**
        	@return - returns whether or not we should use the state file for setup, with preference over the zmatrix file
        */
        bool doSetupFromStateFile() {return useStatefileSetup;}

		/**
			@return - return whether or not we should use the ZMatrix file for setup, in the event the state file isn't found
		*/
		bool doSetupFromZMatrixFile() {return useZMatrixSetup;}
        

};


class OplsScanner
{
   private:
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
    string fileName;
	 /**
     Vector used to keep track of errors when lines in 
	  the OPLS file that don't match the correct format.
    */
	 vector<int> errLines;
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
   public:
   /**
        Constrctor for the OplsScanner object.
        @param fileName - the path to the OPLS file
     */
    OplsScanner(); // constructor
    ~OplsScanner();
		
		/**
		Scans in the opls File calls sub-function addLineToTable
		@param filename - the name/path of the opls file
        @return - success code 
                  0: successful
                  1: error
		*/
      bool readInOpls(string filename); 
		
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
		Real getSigma(string hashNum);
		
		/**
		Returns the epsilon value based on the hashNum (1st col) in Z matrix file
		@param hashNum - the hash number (1st col) in Z matrix file
        @return - the epsilon value of the atom that is the value associated with the hashNum key.
		*/
		Real getEpsilon(string hashNum);
		
		/**
		Returns the charge value based on the hashNum (1st col) in Z matrix file
		@param hashNum -  the hash number (1st col) in Z matrix file
        @return - the charge value of the atom that is the value associated with the hashNum key.
		*/
		Real getCharge(string hashNum);
		
		/**
		Returns the V values value based on the hashNum (1st col) in Z matrix file
		@param hashNum -  the hash number (1st col) in Z matrix file
        @return - TODO
		*/
		Fourier getFourier(string hashNum);

};


class ZmatrixScanner
{
   private:
      /**
         The path to the Z-matrix file
      */
      string fileName;
      /**
        Opls_Scan object used to assign sigma, epsilon and charge values.
      */
      OplsScanner* oplsScanner;
      /**
        Vector that holds example molecules.
      */
      vector<Molecule> moleculePattern;
      /**
        Vector that holds the atoms contained in the molecule in the Z-matrix.
      */
      vector<Atom> atomVector;
      /**
        Vector that holds the bonds contained in the molecule in the Z-matrix.
      */
      vector<Bond> bondVector;
      /**
        Vector that holds the angles contained in the molecule in the Z-matrix.
      */
      vector<Angle> angleVector;
      /**
        Vector that holds the dihedrals contained in the molecule in the Z-matrix.
      */
      vector<Dihedral> dihedralVector;
      /**
        Vector of dummy atoms held seperately from the normal atoms.
      */
      vector<unsigned long> dummies;
      /**
        Global variable that determines if there are multiple molecules in a Z-matrix.
      */
      bool startNewMolecule;
      /**
        Global variable used to store the format id of the last line scanned.
		  needed to read the bottom sections of the Z-matrix file.
      */
      int previousFormat;

   public:
      ZmatrixScanner(); // constructor
      ~ZmatrixScanner();
		
		/**
          Scans in the z-matrix File calls sub-function parseLine
          @param filename - the name/path of the z-matrix file
          @return - success code
                    0: Valid z-matrix path
                    1: Invalid z-matrix path
		*/
        bool readInZmatrix(string filename, OplsScanner* oplsScanner); 
		
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
        int checkFormat(string line);

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
};


class StateScanner
{
  private:
    std::string universal_filename;
    void parsePrimaryIndexDefinitions(Environment* enviro, string definitions);

  public:
    StateScanner(std::string filename);
    ~StateScanner();

    /*
      @param filename - the name of the state file
      @return - the environment recorded in the state file.
    */
    Environment* readInEnvironment();

    /**
      @pararm filename - the name of the state file
      @return - an array of molecules
    */
    vector<Molecule> readInMolecules();

    /**
      @return - the starting step number in the state file
    */
    long readInStepNumber();

    /**
      expected input line:
      "atom1 atom2 value [0|1]"
      0 represents not variable
      1 represents variable angle

      @param line - line containing information about the angle
      @return - returns a bond 
    */
    Angle getAngleFromLine(string line);

    /**
      input line is of the format:
      "id x y z sigma epsilon"
      @param line - the line from the state file that contains atom information.
      @return - atom containing information from the line
    */
    Atom getAtomFromLine(string line);

    /**
      expected input line:
      "atom1 atom2 distance [0|1]"
      0 represents not variable
      1 represents a variable bond
      @param line - the line to be read.
      @return - bond representing the information read from the line.
    */
    Bond getBondFromLine(string line);

    /**
      expected input line:
      "atom1 atom2 value [0|1]"
      0 represents not variable
      1 represents variable dihedral 

      @param line - line containing information about the dihedral 
      @return - returns the dihedral represented on the line
    */
    Dihedral getDihedralFromLine(string line);

    /**
      input line is of the format:
      "x y z numOfAtoms"
      @param line - the line to be parsed
      @return - Environment read from the line
    */
    Environment* getEnvironmentFromLine(string line);

    /**
      expected input line:
      "atom1 atom2 hopDistance"

      @param line - line containing information about the hop
      @return - a hop represented by the information on the line.
    */
    Hop getHopFromLine(string line);
    
    /**
    	Once the state file has been read in,
    	passes the Environment for pre-box setup.
    	[If you attempt to get the Environment before
    	the state file is successfully read in, passes NULL].
    */
    
    
    /**
    	Once the state file has been read in,
    	passes the Molecule for pre-box setup.
    	[If you attempt to get the Molecule before
    	the state file is successfully read in, passes NULL].
    	@param: [None]
    	@return: a singular molecule for the environment
    */

    /**
      @param environment - the environment state
      @param molecules - array of molecules to be printed out
      @param numOfMolecules - the number of molecules to be written out
      @param fileName - the name of the file to be written
    */
    void outputState(Environment *environment, Molecule *molecules, int numOfMolecules, int step, string filename);
};

/**
*	Properly sets up the box for use in the simulation.
*	Builds the environment based on either the configuration file or the state file,
*	as determined elsewhere.
*
* @param: inputPath: the path to either the config file OR the state file, depending on logic path chosen
* @param: inputType: Whether we are running from the config file + zmatrix combo, or from the state file
* @param: box: the generic box into which all the environment information
*		is to be stored; can be modified as necessary for serial or parallel use later.
* @param: startStep: 
* @param: steps: 
*
* @returns: returns TRUE if completed successfully, or FALSE if there was a show-stopping error
*		for which you should do halt the simulation
*/
bool loadBoxData(string inputPath, InputFileType inputType, Box* box, long* startStep, long* steps);


/*************************
*Builds data for the simulation box/environment from the Zmatrix and config files, if
* so desired. Note: this doesn't fully enable the box to be used for the simulation; the geometry
* information must be set up. There's a call to generatefcc box at the end of the method, 
* at which point the box will actually be "generated"/ready.
*
* @param: enviro: the environment stored in the configuration file
* @param: molecVec: the array of molecules, partly created from the Zmatrix file
* @param: box: the box in which all the data is to be stored and the environment set up
*
* @return: returns TRUE if completed successfully, or FALSE if there was a show-stopping error
*		for which you should do halt the simulation
*/
bool buildBoxData(Environment* enviro, vector<Molecule>& molecVec, Box* box);

/*************************
*	Uses data from a given state file to reconstruct a box/environment, as the box looked
* at the time of the statefile's capture during a prior simulation run.
*	Includes geometry information, etc, so the box will be ready for use in the current 
* simulation immediately after the data has been copied over to the box.
*
* @param: enviro: the Environment data from the statefile
* @param: molecVec: the vector array of molecules, copied from the statefile
* @param: box: the box for which data will be built, into which the data will go
*
* @return: returns TRUE if completed successfully, or FALSE if there was a show-stopping error
*		for which you should do halt the simulation
*/
bool fillBoxData(Environment* enviro, vector<Molecule>& moleVec, Box* box);

/*************************
*	Once the box has been filled with environment data, takes the data and 
*	sets up the geometry and other aspects necessary to successfully run a simulation.
* Data must be loaded by either buildBoxData or fillBoxData, depending on if you are
* loading a clean simulation, or loading from the state file.
*	
* @param: box: the box to be used for simulation, with all data pre-loaded	
*
*	@return: returns TRUE if completed successfully, or FALSE if there was a show-stopping error
*		for which you should do halt the simulation
*/
bool generatefccBox(Box* box);

/*************************
*This method allows for writing to a given log file, with some measure of automation.
* These functions are namespace-less and class-less.
* Said string contains almost all of the text you wish to be written to the file.
* Opening and closing of the file will be done on the fly, and a closed file should be guaranteed once this method reaches its end.
* If stringstream is needed, you may call the overloaded version below, which relies on this version of the method.
*@param: text: the text to be written to the output file. Type: string.
*@param: stamp: under which category we log this text: start of simulation [START], end of simulation [END],
*		error while handling the OPLS files [OPLS] Zmatrix files [Z_MATRIX] or geometric config files [GEO], or generic [DEFAULT].
*
*@returns: [none]
*/
void writeToLog(string text,int stamp);

/**
*This method allows for writing to a given log file, with some measure of automation
* This overloaded version allows for a stringstream, instead of a normal string, to be input.
* Said string contains almost all of the text you wish to be written to the file.
* Opening and closing of the file will be done on the fly, and should be guaranteed once this method reaches its end.
*@param: ss: the text to be written to the output file. Type: Reference to a Stringsteam
*@param: stamp: under which category we log this text: start of simulation [START], end of simulation [END],
*		error while handling the OPLS files [OPLS] Zmatrix files [Z_MATRIX] or geometric config files [GEO], or generic [DEFAULT].
*
*@returns: [none]
*/
void writeToLog(stringstream& ss, int stamp);

static inline std::string& rtrim(std::string& s);
static inline std::string& ltrim(std::string& s);
static inline std::string& trim(std::string& s);

#endif
