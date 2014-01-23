/*!\file
  \brief Structures and functions used to read and write Z-matrix files.
  \author Alexander Luchs, Riley Spahn, Seth Wooten
 
 */
#ifndef ZMATRIX_SCAN_H
#define ZMATRIX_SCAN_H

// writing on a text file
#include <vector>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <queue>
#include "metroUtil.h"
#include "Opls_Scan.h"
#include "geometricUtil.h"
using namespace std;


class Zmatrix_Scan
{
   private:
      /**
         The path to the Z-matrix file
      */
      string fileName;
      /**
        Opls_Scan object used to assign sigma, epsilon and charge values.
      */
      Opls_Scan *oplsScanner;
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
      Zmatrix_Scan(string filename, Opls_Scan* oplsScannerRef); // constructor
      ~Zmatrix_Scan();
		
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
#endif // ZMATRIX_SCAN_H
