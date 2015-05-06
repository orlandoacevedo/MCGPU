/*
	Contains the methods required to calculate energies serially.

	Created: February 21, 2014
	
	-> February 26, by Albert Wallace
	-> March 28, by Joshua Mosby
	-> April 21, by Nathan Coleman
*/

#ifndef SERIALCALCS_H
#define SERIALCALCS_H

#include <string>
#include "Metropolis/Box.h"
#include "SerialBox.h"
#include "NeighborList.h"
#include "Metropolis/DataTypes.h"
#include "Metropolis/SimulationArgs.h"
#include "Metropolis/Utilities/StructLibrary.h"

// Linked-cell neighbor list constants


namespace SerialCalcs
{
	/// Factory method for creating a Box from a configuration file.
	/// 	@param configpath The path to the configuration file.
	/// 	@param steps The number of steps desired in the simulation,
	/// 	@return Returns a pointer to the filled-in Box.
	/// 	@note This functionality should ideally reside in SerialBox,
	///   		but it was placed here due to time constraints.
	///   		TODO for future group.
	Box* createBox(std::string inputPath, InputFileType inputType, long* startStep, long* steps);
	
	/* ------ Original System Energy Calculation Functions ------ */
	
	/// Calculates the system energy using consecutive calls to
	///   calcMolecularEnergyContribution.
	/// 	@param box A SerialBox containing the molecule data.
	/// 	@param subLJ A reference to the subtotal energy from the LJ values.
	/// 	@param subCharge A reference to the subtotal energy from the charge energies.
	/// 	@return Returns total system energy.
	Real calcSystemEnergy(Box* box, Real &subLJ, Real &subCharge);
	
	/// Calculates the inter-molecular energy contribution of a given molecule
	/// 	@param box A SerialBox containing the molecule data.
	/// 	@param subLJ A reference to the subtotal energy from the LJ values.
	/// 	@param subCharge A reference to the subtotal energy from the charge energies.
	/// 	@param currentMol the index of the current changed molecule.
	/// 	@param startIdx The optional starting index for other molecules.
	///   		Used for system energy calculation.
	/// 	@return Returns total inter-molecular energy contribution
	Real calcMolecularEnergyContribution(Box* box, Real &subLJ, Real &subCharge, 
		int currentMol, int startIdx = 0);
	
	/* ------ Neighbor-List System Energy Calculation Functions ------ */
	
	/// Calculates the system energy using neighbor linked-list cells (NLC).
	/// 	@param box A Box containing the molecule data.
	/// 	@param subLJ A reference to the subtotal energy from the LJ values.
	/// 	@param subCharge A reference to the subtotal energy from the charge energies.
	/// 	@return Returns total system energy.
	Real calcSystemEnergy_NLC(Box* box, Real &subLJ, Real &subCharge);
	
	/// Finds the neighbor molecules of a given molecule using the neighbor-list structure.
	/// 	@param box A Box containing the molecule data.
	/// 	@param currentMol the index of the current changed molecule.
	///		@param neighbors A reference to the vector where neighbors are stored.
	/// 	@param isSysCalc A boolean for full system version to avoid double-counting.
	void getNeighbors_NLC(Box* box, int currentMol, std::vector<int>& neighbors, bool isSysCalc);
	
	/// Calculates the inter-molecular energy contribution of a given molecule
	///		using the neighbor-list structure.
	/// 	@param box A Box containing the molecule data.
	/// 	@param subLJ A reference to the subtotal energy from the LJ values.
	/// 	@param subCharge A reference to the subtotal energy from the charge energies.
	/// 	@param currentMol the index of the current changed molecule.
	///		@param neighbors Vector containing neighbor molecules of currentMol.
	/// 	@return Returns total inter-molecular energy contribution
	Real calcMolecularEnergyContribution_NLC(Box* box, Real &subLJ, Real &subCharge,
		int currentMol, std::vector<int> neighbors);
	
	/* ------ Utility Calculation Functions ------ */
		
	/// Calculates the inter-molecular energy between two given molecules.
	/// 	@param molecules A pointer to the Molecule array.
	/// 	@param mol1 The index of the first molecule.
	/// 	@param mol2 The index of the second molecule.
	/// 	@param enviro A pointer to the Environment for the simulation.
	/// 	@param subLJ A reference to the subtotal energy from the LJ values.
	/// 	@param subCharge A reference to the subtotal energy from the charge energies.
	/// 	@return Returns the intermolecular energy between the two specified
	///   		molecules.
	Real calcInterMolecularEnergy(Molecule *molecules, int mol1, int mol2, Environment *environment, 
		Real &subLJ, Real &subCharge);
	
	/// Calculates the intra-molecular energy between atoms of the given molecule(s).
	/// 	@param box A Box containing the molecule data.
	/// 	@return Returns the intermolecular energy between the two specified
	///   		molecules.
	Real calcIntraMolecularEnergy(Box* box, Real &subLJ, Real &subCharge);
	
	/// Calculates the long-range correction energy value for molecules outside the cutoff.
	/// 	@note *** This needs to be updated if the volume ever changes or to support more than 2 solvents
	/// 	@param box A Box containing the molecule data.
	/// 	@return Returns r2, the distance between the two atoms squared.
	Real calcEnergy_LRC(Box* box);
	
	/// Calculates the LJ energy between two atoms.
	/// 	@param atom1 The first atom.
	/// 	@param atom2 The second atom.
	/// 	@param r2 The distance between the two atoms, squared.
	/// 	@return Returns the LJ energy between the two specified atoms.
	Real calc_lj(Atom atom1, Atom atom2, Real r2);
	
	/// Calculates the charge energy between two atoms.
	/// 	@param charge1 The charge of atom 1.
	/// 	@param charge2 The charge of atom 2.
	/// 	@param r The distance between the two atoms.
	/// 	@returns Returns the charge energy between two atoms.
	Real calcCharge(Real charge1, Real charge2, Real r);
	
	/// Makes a distance periodic within a specified range.
	/// 	@param x The distance to be made periodic.
	/// 	@param boxDim The magnitude of the periodic range.
	/// 	@return Returns the periodic distance.
	Real makePeriodic(Real x, Real boxDim);
	
	/// Calculates the geometric mean of two values.
	/// 	@param d1 The first value.
	/// 	@param d2 The second value.
	/// 	@return Returns the geometric mean of the two supplied
	///   		values.
	Real calcBlending(Real d1, Real d2);
	
	/// Calculates the squared distance between two atoms.
	/// 	@param atom1 The first atom.
	/// 	@param atom2 The second atom.
	/// 	@param enviro A pointer to the Environment for the simulation.
	/// 	@return Returns r2, the distance between the two atoms squared.
	Real calcAtomDist(Atom atom1, Atom atom2, Environment *enviro);
	
}

#endif