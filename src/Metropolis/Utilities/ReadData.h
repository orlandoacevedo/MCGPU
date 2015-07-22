#ifndef READ_DATA_H
#define READ_DATA_H

#include "../DataTypes.h"
#include<stdlib.h>
#include<string>
#include<map>

using namespace std;

/**
 * This struct contains data about the bond between two atoms. 
 */
struct bondData {
	
	/**
	 * The force constant of this bond.
	 */
	Real kBond;
	
	/**
	 * The equilibrium bond distance of this bond.
	 */
	Real eqBondDist;
	
	/**
	 * Default constructor. Used when creating a map to bondData.
	 */
	bondData() {}
	
	/**
	 * Preferred constructor. Initializes variables.
	 * @param init_kBond The value to initialize the bond's force constant to.
	 * @param init_eqBondDist The value to initialize the bond's equilibrium distance to.
	 */
	bondData(Real init_kBond, Real init_eqBondDist) {
		kBond = init_kBond;
		eqBondDist = init_eqBondDist;
	}
};

/**
 * This struct holds data about the angle formed by three atoms.
 */
struct angleData {
	
	/**
	 * The force constant of the angle.
	 */
	Real kAngle;
	
	/**
	 * The equilibrium size of the angle formed by the three atoms.
	 * This is measured in degrees.
	 */
	Real eqAngle;
	
	/** 
	 * Default constructor for angleData. Used when creating a map to angleData.
	 */
	angleData() {}
	
	/**
	 * Preferred constructor for angleData. Initializes variables.
	 * @param init_kAngle The value to initialize the angle's force constant to.
	 * @param init_eqAngle The value to initialize the angle's size to.
	 */
	angleData(Real init_kAngle, Real init_eqAngle) {
		kAngle = init_kAngle;
		eqAngle = init_eqAngle;
	}
};

class SBScanner {
	
	private: 
		
		/**
		 * Double map from one atom to another atom to data about the bond between the two atoms.
		 */
		map<string, map<string, bondData> > bondDataMap;
		
		/**
		 * Triple map from one atom to another atom to a third atom about the angle formed by those three atoms.
		 */
		map<string, map<string, map<string, angleData> > > angleDataMap;
	
		/**
		 * Contains the path to the oplsaa.sb file.
		 */
		string fileName;
		
		/**
		 * Given a line of input that specifies the data pertaining to a bond, stores the data in the bondDataMap.
		 * @param line A line of data from the OPLSAA.sb file.
		 */
		void processBond(string line);
		
		/**
		 * Given a line of input that specifies the data pertaining to an angle, stores the data in the angleDataMap.
		 * @param line A line of data from the OPLSAA.sb file.
		 */
		void processAngle(string line);
		
		/**
		 * Given a line of input, stores the data in the pertinent map.
		 * @param line A line of data from the OPLSAA.sb file.
		 */
		void processLine(string line);
		
	public:

		/**
		 * Constructor for SBScanner;
		 */
		SBScanner();
		
		/**
		 * Destructor for SBScanner;
		 */
		~SBScanner();
		
		/**
		 * Reads in OPLS.sb and stores the bond and angle data.
		 * @param fileName The path and name of the file to read from.
		 * @return - success code
		 * 			 0: successful
		 *		     1: error
		 */
		bool readInSB(string filename);
		
		/**
		 * Given a pair of atoms, returns the kBond between them.
		 * @param atom1 One of the atoms in the bond.
		 * @param atom2 The other atom in the bond.
		 * @return The force constant of the bond.
		 */
		Real getKBond(string atom1, string atom2);
		
		/**
		 * Given a pair of atoms, returns the eqBondDist between them.
		 * @param atom1 One of the atoms in the bond.
		 * @param atom2 The other atom in the bond.
		 * @return The equilibrium bond distance.
		 */
		Real getEqBondDist(string atom1, string atom2);
		
		/**
		 * Given three atoms, returns the kAngle of the angle that they form.
		 * @param endpoint1 The atom at one end of the angle.
		 * @param middleAtom The atom in the middle of the angle.
		 * @param endpoint2 The atom at the other end of the angle.
		 * @return The force constant of the angle.
		 */
		Real getKAngle(string endpoint1, string middleAtom, string endpoint2);
		
		/**
		 * Given three atoms, returns the eqAngle of the angle that they form.
		 * @param endpoint1 The atom at one end of the angle.
		 * @param middleAtom The atom in the middle of the angle.
		 * @param endpoint2 The atom at the other end of the angle.
		 * @return The equilibrium angle formed by the three atoms (in degrees).
		 */
		Real getEqAngle(string endpoint1, string middleAtom, string endpoint2);
		
};

#endif