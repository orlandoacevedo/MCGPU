/*!\file
  \brief Structures and functions used to read and write OPLS files.
  \author Alexander Luchs, Riley Spahn, Seth Wooten
 
 */
#ifndef OPLS_SCAN_H
#define OPLS_SCAN_H

// writing on a text file
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include "metroUtil.h"
using namespace std;

/**
  Structure used to represent the 4 Fourier Coeficients
*/
struct Fourier
{
    double vValues[4];
};

/**
  Class used to read and write Opls files.
*/
class Opls_Scan
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
        Constrctor for the Opls_Scan object.
        @param fileName - the path to the OPLS file
     */
    Opls_Scan(string filename); // constructor
    ~Opls_Scan();
		
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

};

#endif //OPLS_SCAN_H
