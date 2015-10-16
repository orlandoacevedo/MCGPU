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

	/**
	* Given an input matrix, returns the result of inverting the given matrix.
	* @param m The matrix to invert.
	* @param out The resultant, inverted matrix is written here.
	*/
	void invert_3x3_Matrix(Real m[3][3], Real out[3][3]);

	/**
	* Given a 3x3 matrix and a 3x1 matrix, multiplies the two into a resultant
	* 3x1 matrix.
	* @param iMat The 3x3 matrix to multiply.
	* @param coords The 3x1 matrix to multiply by.
	* @param out The resultant matrix to write to.
	*/
	void multiply_matrices(Real iMat[3][3], Real coords[3], Real out[3]);

	/**
	* Given the description of a point, along with a description of the rotation, returns the
	* transformed point.
	* @param vSet The vector description of the space.
	* @param origin The point about which the rotation occurs.
	* @param cosFactor Equal to dcos(theta).
	* @param sinFactor Equal to dsin(theta).
	* @param abc Holds, a, b, and c, where the point is represented by O + aV1 + bV2 + cV3.
	* @param out The resultant location of the point is written here.
	*/
	void getP_Next(Real vSet[3][3], Real origin[3], Real cosFactor, Real sinFactor, Real abc[3], Real out[3]);

	/**
	* Calculates the cross product of two 3-dimensional vectors.
	* @param v1 The first vector in the cross product.
	* @param v2 The second vector in the cross product.
	* @param out The resultant vector is written here.
	*/
	void cross_3d_Vectors(Real v1[3], Real v2[3], Real out[3]);

	/**
	* Calculates the length of a vector.
	* @param vector The vector to calculate the length of.
	* @return The length of the vector.
	*/
	Real len(Real vector[3]);

	/**
	* Given two points and/or position vectors, creates a vector between them.
	* @param from Where the vector originates.
	* @param to Where the vector terminates.
	* @param result The vector between the two positions.
	*/
	void makeVector(Real from[3], Real to[3], Real result[3]);

	/**
	* Given the three points that compose the angle, calculates the two essential
	* vectors that lie within the plane of rotation.
	* @param origin The center point of the angle.
	* @param anglePoint1 One endpoint of the angle.
	* @param anglePoint2 The other endpoint of the angle.
	* @param v1 The first resultant vector is written here.
	* @param v2 The second resultant vector is written here.
	*/
	void calc_v1_and_v2(Real origin[3], Real anglePoint1[3], Real anglePoint2[3], Real v1[3], Real v2[3]);

	/**
	* Given a set of points describing an angle, determines the appropriate subspace of R3 to use for calculations,
	* along with the inverse of the subspace.
	* @param o The center point of the angle.
	* @param p One endpoint of the angle.
	* @param q The other endpoint of the angle.
	* @param vInverseMatrix The inverse of the subspace matrix is written here.
	* @param vMatrix The spanning matrix is written here.
	*/
	void get_V_Inverse_Matrix(Real o[3], Real p[3], Real q[3], Real vInverseMatrix[3][3], Real vMatrix[3][3]);

	/**
	* Given three points, returns the angle that they form.
	* @param o The center point of the angle.
	* @param p The first endpoint of the angle.
	* @param q The second endpoint of the angle.
	*/
	Real getAngle(Real o[3], Real p[3], Real q[3]);

	/**
	* Given a description of an angle and the change in degrees, along with the set of atoms to move,
	* expands the angle and moves all relevant points accordingly.
	* @param o The center point of the angle.
	* @param p One endpoint of the angle.
	* @param q The other endpoint of the angle.
	* @param pAtoms The set of all atoms (including p) on same side of the angle as p in the molecule.
	* @param qAtoms The set of all atoms (including q) on same side of the angle as q in the molecule.
	* @param oldAngle The original measure of the angle (in degrees). -- TODO -- Remove the need for this parameter.
	* @param newAngle The new measure of the angle (in degrees).
	*/
	void expandAngle(Real o[3], Real p[3], Real q[3], Atom pAtoms[], Atom qAtoms[], Real newAngle);


	/**
	* Given a molecule, an angle, the id of a common atom, and the amount to change the angle (in degrees),
	* expands the angle symmetrically.
	* @param m Points to the molecule the target angle is a part of.
	* @param a Points to the target angle.
	* @param commonAtom The id of the atom that is at the center of the angle.
	* @param thetaChange The amount (in degrees) to change the angle. Positive measures increase the angle,
	* and negative measures decrease it.
	* @return true if successful, false if there is an error or if the angle is part of a ring.
	*/
	bool expandMoleculeAngle(Molecule* m, Angle* a, int commonAtom, Real thetaChange);

	/**
	* Helper method used to split the molecule into two groups, on either side of the angle.
	* This method recursively calls itself.
	* @param edgeList A list of all edges in (Bonds) in the molecules.
	* @param currentAtom The index of the current Atom.
	* @param offset The index of the first atom in the molecule.
	* @param atomStatus An array used to track whether each atom has been visited or not.
	* @param type The side of the bond on which the atom being visited falls (1 or 2).
	* @param edgeListSize The size of edgeList.
	*/
	void DFS(Bond edgeList[], int currentAtom, int offset, int atomStatus[], int type, int edgeListSize);

}

#endif
