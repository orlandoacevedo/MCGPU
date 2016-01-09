#ifndef SIMBOX_H
#define SIMBOX_H

#include "Box.h"
#include "Utilities/StructLibrary.h"
#include "Utilities/MathLibrary.h"

#include <stdlib.h>
#include <vector>
#include <set>

typedef double Real;
typedef unsigned int ID;
typedef const int & refInt;


class SimBox {

public:

  // BEGINNING OF PUBLIC CONSTANTS

  // DIMENSION CONSTANTS

  /**
   * Constant to represent the X-axis.
   */
  static const int X_COORD = 0;

  /**
   * Constant to represent the Y-axis.
   */
  static const int Y_COORD = 1;

  /**
   * Constant to represent the Z-axis.
   */
  static const int Z_COORD = 2;

  /**
   * Constant holds the number of spatial dimensions in the simulation.
   */
  static const int NUM_DIMENSIONS = 3;

  // MOLECULE DATA CONSTANTS

  /**
   * Indicates the row of moleculeData that holds the start index of each
   *     molecule in atomCoordinates and atomData.
   */
  static const int MOL_START = 0;

  /**
   * Indicates the row of moleculeData that holds the number of atoms of each
   *     molecule.
   */
  static const int MOL_LEN = 1;

  /**
   * Indicates the row of moleculeData that holds the start index of each
   *     molecule's primary index(es) in primaryIndexes.
   */
  static const int MOL_PIDX_START = 2;

  /**
   * Indicates the row of moleculeData that holds the number of primary indexes
   *     that each molecule has.
   */
  static const int MOL_PIDX_COUNT = 3;

  /**
   * Indicates the row of moleculeData that hold the type of each molecule.
   */
  static const int MOL_TYPE = 4;

  /**
   * Indicates the number of rows of moleculeData.
   */
  static const int MOL_DATA_SIZE = 5;

  // ATOM DATA CONSTANTS

  /**
   * Indicates the row of atomData that holds the value of sigma for each atom.
   */
  static const int ATOM_SIGMA = 0;

  /**
   * Indicates the row of atomData that holds the value of epsilon for each atom.
   */
  static const int ATOM_EPSILON = 1;

  /**
   * Indicates the row of atomData that holds the charge of each atom.
   */
  static const int ATOM_CHARGE = 2;

  /**
   * Indicates the number of rows of atomData.
   */
  static const int ATOM_DATA_SIZE = 3;

  // BOND DATA CONSTANTS

  /**
   * Indicates the row of bondData that holds the 1st atom index for each bond.
   */
  static const int BOND_A1_IDX = 0;

  /**
   * Indicates the row of bondData that holds the 2nd atom index for each bond.
   */
  static const int BOND_A2_IDX = 1;

  /**
   * Indicates the row of bondData that holds the force constant for each bond.
   */
  static const int BOND_KBOND = 2;

  /**
   * Indicates the row of bondData holding the equilibrium distance of each bond.
   */
  static const int BOND_EQDIST = 3;

  /**
   * Indicates the row of bondData that records whether each bond is variable.
   */
  static const int BOND_VARIABLE = 4;

  /**
   * Indicates the number of rows in bondData.
   */
  static const int BOND_DATA_SIZE = 5;

  // ANGLE DATA CONSTANTS

  /**
   * Indicates the row of angleData that holds the first atom's index for each
   *     angle.
   */
  static const int ANGLE_A1_IDX = 0;

  /**
   * Indicates the row of angleData that holds the middle atom's index for each
   *     angle.
   */
  static const int ANGLE_MID_IDX = 1;

  /**
   * Indicates the row of angleData that holds the third atom's index for each
   *     angle.
   */
  static const int ANGLE_A2_IDX = 2;

  /**
   * Indicates the row of angleData that holds the force constant of each angle.
   */
  static const int ANGLE_KANGLE = 3;

  /**
   * Indicates the row of angleData holding the equilibrium size of each angle.
   */
  static const int ANGLE_EQANGLE = 4;

  /**
   * Indicates the row of angleData holding whether each angle is variable.
   */
  static const int ANGLE_VARIABLE = 5;

  /**
   * Indicates the number of rows in angleData.
   */
  static const int ANGLE_DATA_SIZE = 6;

  // BEGINNING OF MEMBER VARIABLES

  // Basic environment conditions.

  /**
   * Real[3]. Holds the box's dimensions.
   */
  Real*   size;

  /**
   * Holds the box's temperature.
   */
  Real    temperature;

  /**
   * Holds the box's cutoff distance.
   */
  Real    cutoff;

  /**
   * Holds the maximum translation of a molecule.
   */
  Real    maxTranslate;

  /**
   * Holds the maximum rotation of a molecule (in degrees).
   */
  Real    maxRotate;

  /**
   * Holds the number of steps to run.
   */
  int numSteps;

  // Item sizes

  /**
   * The number of molecules in the box.
   */
  int numMolecules;

  /**
   * The number of atoms in the box.
   */
  int numAtoms;

  /**
   * The number of bonds in the box. (Currently unused).
   */
  int numBonds;

  /**
   * The number of angles in the box. (Currently unused).
   */
  int numAngles;

  /**
   * The number of dihedrals in the box. (Currently unused).
   */
  int numDihedrals;

  // Molecule information

  /**
   * int[MOL_DATA_SIZE][numMolecules]
   * Holds information about the start atom index, the number of atoms in the
   * molecule, the start of the molecule's primary indexes, and the number of
   * primary indexes in the molecule.
   */
  int** moleculeData;

  /**
   * int[#of total primary indexes]
   * Holds the primary indexes for every atom.
   */
  int* primaryIndexes;

  // Atom information

  /**
   * Real[3][numAtoms]
   * Holds the coordinates of every atom in the box. Which atom belongs to which
   *     molecule is specified in moleculeData.
   */
  Real** atomCoordinates;

  /**
   * Real[3][# of atoms in largest molecule]
   * Holds the previous coordinates of every atom in the most recently moved
   *      molecule. Used when rolling back a molecule following a rejected step.
   */
  Real** rollBackCoordinates;

  /**
   * Real[ATOM_DATA_SIZE][numAtoms]
   * Holds constant information about every atom, including sigma, epsilon, and
   *     the atom's charge.
   */
  Real**  atomData;

  // Bond information -- Currently unused.

  /**
   * Real[numBonds]
   * Holds the length of every bond in the box.
   */
  Real*   bondLengths;

  /**
   * Real [BOND_DATA_SIZE][numBonds]
   * Holds data about every bond in the box, including the endpoints, the
   *     force constant, and the equilibrium distance.
   */
  const Real**  bondData;

  // Angle information -- Currently unused.

  /**
   * Real [numAngles]
   * Holds the size of every angle in the box, in degrees.
   */
  Real*   angleSizes;

  /**
   * Real[ANGLE_DATA_SIZE][numAngles]
   * Holds data about every angle in the box, including the endpoints, central
   *     atom, force constant, and the equilibrium angle.
   */
  const Real**  angleData;

  // Dihedral information -- Currently unused.

  /**
   * Real [numDihedrals]
   * Holds the size of every dihedral in the box.
   */
  Real*   dihedralSizes;

  /**
   * Real[?][numDihedrals]
   * Will hold data about every dihedral in the box. Currently unimplemented.
   */
  const Real**  dihedralData;

  //NLC information.

  /**
   * True if NLC is being used, false otherwise.
   */
  bool useNLC;

  /**
   * int[3]
   * Holds the number of cells in each dimension.
   */
  int* numCells;

  /**
   * Real [3]
   * Holds the width of each cell in each dimension.
   */
  Real* cellWidth;

  /**
   * Points to NLC_Node[numMolecules]
   * Used only for allocation / deallocation of memory for the NLC linked cells.
   */
  NLC_Node* nlc_heap;

  /**
   * NLC_Node*[3^3]
   * Holds each linked list for cells adjacent to the target molecule.
   */
  NLC_Node** neighbors;

  /**
   * NLC_Node*[numCells[0]][numCells[1]][numCells[2]]
   * NLC heads for every cell in the box.
   */
  NLC_Node**** neighborCells;

  /**
   * Points to the node before the node just examined, or
   *     to the node just examined if that node is the head of its linked list.
   */
  NLC_Node* prevNode;

  /**
   * int[3]
   * Holds the location of the cell most recently examined.
   */
  int* prevCell;

  //BEGINNING OF FUNCTIONS

  /**
   * Roll back a molecule to its original poisition. Performs translation and
   * rotation.
   *
   * @param molIdx The index of the molecule to change.
   */
  void rollback(int molIdx);


  /**
   * Keeps a molecule in the simulation box's boundaries, based on the location
   * of its first primary index atom.
   *
   * @param molIdx The index of the molecule to keep inside the box.
   */
  void keepMoleculeInBox(int molIdx);

  /**
   * Given a distance, makes the distance periodic to mitigate distances greater
   * than half the length of the box.
   *
   * @param x The measurement to make periodic.
   * @param dimension The dimension the measurement must be made periodic in.
   * @return The measurement, scaled to be periodic.
   */
  Real makePeriodic(Real x, int dimension);

  /**
   * Returns the index of a random molecule within the simulation box.
   *
   * @return A random integer from 0 to (numMolecules - 1)
   */
  int chooseMolecule() const;

  /**
   * Given a molecule to change, randomly moves the molecule within the box.
   * This performs a translation and a rotation, and saves the old position.
   *
   * @param molIdx The index of the molecule to change.
   */
  void changeMolecule(int molIdx);

  /**
   * Given a molecule to change, and information to change, translates and
   * rotates the specified molecule, saving its old position.
   *
   * @param molIdx The index of the molecule to change.
   * @param vIdx The index of the pivot point of the molecule.
   * @param dX The amount to translate in the x direction.
   * @param dY The amount to translate in the y direction.
   * @param dZ The amount to translate in the z direction.
   * @param rX The amount to rotate around the x axis (degrees).
   * @param rY The amount to rotate around the y axis (degrees).
   * @param rZ The amount to rotate around the z axis (degrees).
   */
  void changeMolecule(int molIdx, int vIdx, Real dX, Real dY, Real dZ, Real rX, Real rY, Real rZ);

  /**
   * Given an atom to translate, and the directions to tranlate it, moves it.
   *
   * @param aIdx The index of the atom to change.
   * @param dX The amount to translate in the x direction.
   * @param dY The amount to translate in the y direction.
   * @param dZ The amount to tranlsate in the z direction.
   */
  void translateAtom(int aIdx, Real dX, Real dY, Real dZ);

  /**
   * Given an atom to rotate, its pivot point, and the amounts to rotate,
   * rotates the atom around the pivot (rotations in degrees).
   *
   * @param aIdx The index of the atom to change.
   * @param pivotIdx The index of the atom about which to rotate.
   * @param rotX The amount of rotation around the x-axis that is done.
   * @param rotY The amount of rotation around the y-axis that is done.
   * @param rotZ The amount of rotation around the z-axis that is done.
   */
  void rotateAtom(int aIdx, int pivotIdx, Real rotX, Real rotY, Real rotZ);

  /**
   * Given an atom and an amount to rotate, rotates about the x-axis.
   *
   * @param aIdx The index of the atom to rotate.
   * @param angleDeg The angle to rotate it (in degrees).
   */
  void rotateX(int aIdx, Real angleDeg);

  /**
   * Given an atom and an amount to rotate, rotates it about the y-axis.
   *
   * @param aIdx The index of the atom to rotate.
   * @param angleDeg The angle to rotate it (in degrees).
   */
  void rotateY(int aIdx, Real angleDeg);

  /**
   * Given an atom and an amount to rotate, rotates it about the z-axis.
   *
   * @param aIdx The index of the atom to rotate.
   * @param angleDeg The angle to rotate it (in degrees).
   */
  void rotateZ(int aIdx, Real angleDeg);

  /**
   * calcSystemEnergy Determines the total energy of the box
   *
   * @param subLJ Initial Lennard - Jones energy. Final L-J energy passed out by reference.
   * @param subCharge Initial Coulomb energy. Final Coulomb energy passed out by reference.
   * @return The total energy of the box.
   */
  Real calcSystemEnergy (Real &subLJ, Real &subCharge);

  /**
   * calcMolecularEnergyContribution Determines the energy contribution of a particular molecule.
   *
   * @param subLJ Initial Lennard - Jones energy. Final L-J energy passed out by reference.
   * @param subCharge Initial Coulomb energy. Final Coulomb energy passed out by reference.
   * @param currMol The index of the molecule to calculate the contribution of.
   * @param startMol The index of the molecule to begin searching from to determine interaction energies.
   * @return The total energy of the box (discounts initial lj / charge energy)
   */
  Real calcMolecularEnergyContribution(Real &subLJ, Real &subCharge, int currMol, int startMol);

  /**
   * moleculesInRange Determines whether or not two molecule's primaryIndexes are
   *                  within the cutoff range of one another.
   *
   * @param p1Start The index of the first primary index from molecule 1.
   * @param p1End index of the last primary index from molecule 1,  + 1.
   * @param p2Start index of the first primary index from molecule 2.
   * @param p2End index of the last primary index from molecule 2,  + 1.
   * @return true if the molecules are in range, false otherwise.
   */
  bool moleculesInRange(refInt p1Start, refInt p1End, refInt p2Start, refInt p2End);

  /**
   * calcMoleculeInteractionEnergy Calcs the energy caused by the interaction between
   *                               a given pair of molecules.
   *
   * @param subLJ The initial Lennard - Jones energy. Final energy is passed out by reference.
   * @param subCharge The initial Coloumb energy. Final energy is passed out by reference.
   * @param m1 The molecule index of the first molecule.
   * @param m2 The molecule index of the second molecule.
   * @return The energy from the molecular interaction.
   */
  Real calcMoleculeInteractionEnergy (Real &subLJ, Real &subCharge, refInt m1, refInt m2);

  /**
   * calcAtomDistSquared Calculates the square of the distance between two atoms.
   *
   * @param a1 The atom index of the first atom.
   * @param a2 The atom index of the second atom.
   * @return The square of the distance between the atoms.
   */
  Real calcAtomDistSquared(refInt a1, refInt a2);

  /**
   * calcLJEnergy Calculates the Lennard - Jones potential between two atoms.
   *
   * @param a1  The atom index of the first atom.
   * @param a2  The atom index of the second atom.
   * @param r2  The distance between the atoms, squared.
   * @return The Lennard - Jones potential from the two atoms' interaction.
   */
  Real calcLJEnergy(refInt a1, refInt a2, const Real& r2);

  /**
   * calcChargeEnergy Calculates the Coloumb potential between two atoms.
   *
   * @param a1 The atom index of the first atom.
   * @param a2 The atom index of the second atom.
   * @param r  The distance between the atoms.
   * @return The Coloumb potential from the two atoms' interaction.
   */
  Real calcChargeEnergy(refInt a1, refInt a2, const Real& r);

  /**
   * calcBlending Calculates the geometric mean of two real numbers.
   *
   * @param a The first real number.
   * @param b The second real number.
   * @return sqrt(|a*b|), the geometric mean of the two numbers.
   */
  Real calcBlending (const Real &a, const Real &b);

  /**
   * Given the index for a molecule, popuplate neighbors with the heads of
   * linked-lists holding the molecules in neighboring cells.
   *
   * The neighbors field variable is filled by this method, but cells beyond the
   * number of neighboring cells are unaffected.
   *
   * @param molIdx The index of the molecule to find the neighboring cells of.
   * @return The number of neighboring cells found.
   */
  int findNeighbors(int molIdx);

  /**
   * Given an index and a dimension, returns the index, wrapped around the box.
   *
   * @param idx The index of a cell to wrap.
   * @param dimension The dimension the index is in.
   * @return The index, wrapped around the size of the box.
   */
  int wrapCell(int idx, int dimension);

  /**
   * Given a location and a dimension, returns the corresponding cell index.
   *
   * @param location The index of the cell to wrap.
   * @param dimension The dimension to wrap in.
   * @return The index, wrapped around the size of the box.
   */
  int getCell(Real loc, int dimension);

  /**
   * Given a molecule, updates the molecule's position in the neighbor linked
   * cells (if it needs to be updated).
   *
   * updateNLC uses the prevNode pointer to avoid redundant traversal of the NLC
   * list. Therefore, it must be set to the node to remove's predecessor, or the
   * remove node itself if that node is the head of the linked list.
   *
   * Inputs: molIdx - The index of the moleucle to update.
   */
  void updateNLC(int molIdx);
};

typedef const SimBox & refBox;

#endif
