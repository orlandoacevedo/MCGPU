#ifndef SIMBOX_H
#define SIMBOX_H

#include "Box.h"
#include "Utilities/StructLibrary.h"
#include "Utilities/MathLibrary.h"
#include "SimBoxConstants.h"
#include "DataTypes.h"

#include <stdlib.h>
#include <vector>
#include <set>

typedef unsigned int ID;

const double kBoltz = 0.00198717;

class SimBox {
 public:
  /**
   * The current step in the simulation. Used for intramolecular delta tweaking
   */
  int stepNum;

  // ----- Basic Environment Conditions -----

  /**
   * Real[3] Holds the box's dimensions.
   */
  Real* size;

  /**
   * Holds the box's temperature.
   */
  Real temperature;

  /**
   * Holds tempeture value for Monte Carlo calculations
   */
  Real kT;

  /**
   * Holds the box's cutoff distance.
   */
  Real cutoff;

  /**
   * Holds the maximum translation of a molecule.
   */
  Real maxTranslate;

  /**
   * Holds the maximum rotation of a molecule (in degrees).
   */
  Real maxRotate;

  /**
   * Holds the absolute value of the maximum change allowed in a bond
   */
  Real maxBondDelta;

  /**
   * Holds the absolute value of the maximum change allowed in an angle
   */
  Real maxAngleDelta;

  /**
   * Holds the number of steps to run.
   */
  int numSteps;

  /** Holds the maximum number of intramolecular moves to make in a step */
  int maxIntraMoves;

  // ----- Item Sizes -----

  /**
   * The number of molecules in the box.
   */
  int numMolecules;

  /**
   * The number of atoms in the box.
   */
  int numAtoms;

  /**
   * The number of atoms in the largest molecule.
   */
  int largestMol;

  /**
   * The number of primary indexes in the box.
   */
  int numPIdxes;

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

  // ----- Molecule Information -----

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

  // ----- Atom information -----

  /**
   * Real[3][numAtoms]
   * Holds the coordinates of every atom in the box. Which atom belongs to
   * which molecule is specified in moleculeData.
   */
  Real** atomCoordinates;

  /**
   * Real[3][# of atoms in largest molecule]
   * Holds the previous coordinates of every atom in the most recently moved
   * molecule. Used when rolling back a molecule following a rejected step.
   */
  Real** rollBackCoordinates;

  /**
   * Real[ATOM_DATA_SIZE][numAtoms]
   * Holds constant information about every atom, including sigma, epsilon, and
   *     the atom's charge.
   */
  Real** atomData;

  // ----- Bond information -----

  /**
   * Real[numBonds]
   * Holds the length of every bond in the box.
   */
  Real* bondLengths;

  /**
   * Real[numBonds]
   * Holds the previous lengths of all the bonds in the simulation box
   */
  Real* rollBackBondLengths;

  /**
   * Real [BOND_DATA_SIZE][numBonds]
   * Holds data about every bond in the box, including the endpoints, the
   * force constant, and the equilibrium distance.
   */
   Real** bondData;

   /**
    * The number of intramolecular moves that have stretched a bond in this
    * simulation. Used for intrmolecular delta tuning.
    */
   int numBondMoves;

   /**
    * The number of intrmolecular bond moves that have passed an MC test after
    * occuring. Used for intrmolecular delta tuning.
    */
   int numAcceptedBondMoves;

  // ----- Angle Information -----

  /**
   * Real [numAngles]
   * Holds the size of every angle in the box, in degrees.
   */
  Real* angleSizes;

  /**
   * Real[numAngles]
   * Holds the previous anlge sizes of the angles in the simulation box
   */
  Real* rollBackAngleSizes;

  /**
   * Real[ANGLE_DATA_SIZE][numAngles]
   * Holds data about every angle in the box, including the endpoints, central
   * atom, force constant, and the equilibrium angle.
   */
  Real** angleData;

   /**
    * The number of intramolecular moves that have altered an angle in this
    * simulation. Used for intrmolecular delta tuning.
    */
  int numAngleMoves;

   /**
    * The number of intramolecular angle moves that have passed an MC test in
    * this simulation. Used for intrmolecular delta tuning.
    */
  int numAcceptedAngleMoves;

  // ----- Dihedral Information (Currently unused) -----

  /**
   * Real[numDihedrals]
   * Holds the size of every dihedral in the box.
   */
  Real* dihedralSizes;

  /**
   * Real[?][numDihedrals]
   * Will hold data about every dihedral in the box. Currently unimplemented.
   */
  const Real** dihedralData;

  // ----- NLC information -----

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
   * Used only for (de)allocation of memory for the NLC linked cells.
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
   * Points to the node before the node just examined, or to the node just
   * examined if that node is the head of its linked list.
   */
  NLC_Node* prevNode;

  /**
   * int[3]
   * Holds the location of the cell most recently examined.
   */
  int* prevCell;

  /**
   * int[atoms in largest molecule]
   * Used for performing union-find calculations to find what half of a bond or
   * angle a molecule is in.
   */
   int* unionFindParent;

  /**
   * int[# of molecule types][# atoms in each type][# of atoms to exclude]
   * Holds information for which intra molecular pairwise LJ and Coloumb
   * interactions to not calculate.
   */
  int*** excludeAtoms;

  /**
   * int[# of molecule types][# atoms in each type][# of atoms to half]
   * Holds information for which intra molecular pairwise LJ and Coloumb
   * interactions to multiply by the fudge factor.
   */
  int*** fudgeAtoms;

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
  void keepMoleculeInBox(int molIdx, Real** aCoords, int** molData,
                         int* pIdxes, Real* bsize);

  /**
   * Given a distance, makes the distance periodic to mitigate distances greater
   * than half the length of the box.
   *
   * @param x The measurement to make periodic.
   * @param dimension The dimension the measurement must be made periodic in.
   * @return The measurement, scaled to be periodic.
   */
  Real makePeriodic(Real x, int dimension, Real* bSize);

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
   * Given an atom to translate, and the directions to tranlate it, moves it.
   *
   * @param aIdx The index of the atom to change.
   * @param dX The amount to translate in the x direction.
   * @param dY The amount to translate in the y direction.
   * @param dZ The amount to tranlsate in the z direction.
   */
  void translateAtom(int aIdx, Real dX, Real dY, Real dZ, Real** aCoords);

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
  #pragma acc routine seq
  void rotateAtom(int aIdx, int pivotIdx, Real rotX, Real rotY, Real rotZ,
                  Real** aCoords);

  /**
   * Given an atom and an amount to rotate, rotates about the x-axis.
   *
   * @param aIdx The index of the atom to rotate.
   * @param angleDeg The angle to rotate it (in degrees).
   */
  void rotateX(int aIdx, Real angleDeg, Real** aCoords);

  /**
   * Given an atom and an amount to rotate, rotates it about the y-axis.
   *
   * @param aIdx The index of the atom to rotate.
   * @param angleDeg The angle to rotate it (in degrees).
   */
  void rotateY(int aIdx, Real angleDeg, Real** aCoords);

  /**
   * Given an atom and an amount to rotate, rotates it about the z-axis.
   *
   * @param aIdx The index of the atom to rotate.
   * @param angleDeg The angle to rotate it (in degrees).
   */
  void rotateZ(int aIdx, Real angleDeg, Real** aCoords);

  /**
   * calcSystemEnergy Determines the total energy of the box
   *
   * @param subLJ Initial Lennard - Jones energy. Final L-J energy passed out
   * by reference.
   * @param subCharge Initial Coulomb energy. Final Coulomb energy passed out
   * by reference.
   * @return The total energy of the box.
   */
  Real calcSystemEnergy(Real &subLJ, Real &subCharge);

  /**
   * calcMolecularEnergyContribution Determines the energy contribution of a particular molecule.
   *
   * @param subLJ Initial Lennard - Jones energy. Final L-J energy passed out
   * by reference.
   * @param subCharge Initial Coulomb energy. Final Coulomb energy passed out
   * by reference.
   * @param currMol The index of the molecule to calculate the contribution of.
   * @param startMol The index of the molecule to begin searching from to
   * determine interaction energies.
   * @return The total energy of the box (discounts initial lj / charge energy)
   */
  Real calcMolecularEnergyContribution(Real &subLJ, Real &subCharge,
                                       int currMol, int startMol);

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
  bool moleculesInRange(int p1Start, int p1End, int p2Start, int p2End);

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
  Real calcMoleculeInteractionEnergy (int m1, int m2);

  /**
   * calcAtomDistSquared Calculates the square of the distance between two atoms.
   *
   * @param a1 The atom index of the first atom.
   * @param a2 The atom index of the second atom.
   * @return The square of the distance between the atoms.
   */
  Real calcAtomDistSquared(int a1, int a2, Real** aCoords, Real* bSize);

  /**
   * calcLJEnergy Calculates the Lennard - Jones potential between two atoms.
   *
   * @param a1  The atom index of the first atom.
   * @param a2  The atom index of the second atom.
   * @param r2  The distance between the atoms, squared.
   * @return The Lennard - Jones potential from the two atoms' interaction.
   */
  Real calcLJEnergy(int a1, int a2, const Real& r2, Real** aData);

  /**
   * calcChargeEnergy Calculates the Coloumb potential between two atoms.
   *
   * @param a1 The atom index of the first atom.
   * @param a2 The atom index of the second atom.
   * @param r  The distance between the atoms.
   * @return The Coloumb potential from the two atoms' interaction.
   */
  Real calcChargeEnergy(int a1, int a2, const Real& r, Real** aData);

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
   * @param molIdx The index of the molecule to update.
   */
  void updateNLC(int molIdx);

  /**
   * Calculates the energy contribution from intramolecular forces within the
   *     given molecule.
   *
   * @param molIdx The index of the molecule to calculate the intramolecular
   *      energy contribution of.
   * @return The intramolecular energy contribution of the molecule.
   */
  Real calcIntraMolecularEnergy(int molIdx);

  /**
   * Stretches or compresses the given bond in the given molecule.
   *
   * @param molIdx The index of the molecule that the bond is in.
   * @param bondIdx The index of the bond within the molecule. For example, the
   *     first bond in molecule 12 has index 0.
   * @param stretchDist The amount to stretch or compress the bond. Positive
   *     values correspond to stretching, negative to compression.
   */
  void stretchBond(int molIdx, int bondIdx, Real stretchDist);

  /**
   * Expands or contracts the given angle in the given molecule.
   *
   * @param molIdx The index of the molecule that the angle is in.
   * @param angleIdx The index of teh angle within the molecule. For example,
   *     the first angle in molecule 12 has index 0.
   * @param expandDeg The amount to expand or contract the angle. Positive
   *     values correspond to expansion, negative to contraction. This is
   *     measured in degrees.
   */
  void expandAngle(int molIdx, int angleIdx, Real expandDeg);

  /**
   * Calculates the energy from various flexible angles within the molecule.
   *
   * @param molIdx The index of the molecule to calculate the energy of.
   * @return The energy produced by angles being different from their
   *     equilibrium measurement.
   */
  Real angleEnergy(int molIdx);

  /**
   * Rolls back the length of the most recently changed bond. Does not change
   *     any atom coordinates.
   */
  void rollbackBond();

  /**
   * Calculates the energy from various flexible bonds within the molecule.
   *
   * @param molIdx The index of the molecule to calculate the energy of.
   * @return The energy produced by bonds being different from their
   *     equilibrium measurement.
   */
   Real bondEnergy(int molIdx);

  /**
   * Rolls back the size of the most recently changed angle. Does not change
   *     any angle coordinates.
   */
  void rollbackAngle();

  /**
   * Unions the sets of two different atoms connected by a bond. Note that atoms
   *     are indexed within the molecule, so the first atom in the 12th molecule
   *     has index 0.
   *
   * @param atom1 The index (within the molecule) of the first atom to union.
   * @param atom2 The index (within the molecule) second atom in the union.
   */
  void unionAtoms(int atom1, int atom2);

  /**
   * Finds the representative element of an atom in the disjoint set array.
   * Note that the atom is indexed within the molecule, so the first atom in the
   * 12th molecule has index 0. This function also performs path compression.
   *
   * @param atomIdx The atom's index within the molecule.
   * @return The representative element of the atom index.
   */
  int find (int atomIdx);


};

typedef const SimBox & refBox;

#endif
