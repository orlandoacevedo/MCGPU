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

  static const int X_COORD = 0;
  static const int Y_COORD = 1;
  static const int Z_COORD = 2;
  static const int NUM_DIMENSIONS = 3;

  static const int MOL_START = 0;
  static const int MOL_LEN = 1;
  static const int MOL_PIDX_START = 2;
  static const int MOL_PIDX_COUNT = 3;
  static const int MOL_TYPE = 4;
  static const int MOL_DATA_SIZE = 5;

  static const int ATOM_SIGMA = 0;
  static const int ATOM_EPSILON = 1;
  static const int ATOM_CHARGE = 2;
  static const int ATOM_DATA_SIZE = 3;

  static const int BOND_A1_IDX = 0;
  static const int BOND_A2_IDX = 1;
  static const int BOND_KBOND = 2;
  static const int BOND_EQDIST = 3;
  static const int BOND_VARIABLE = 4;
  static const int BOND_DATA_SIZE = 5;

  static const int ANGLE_A1_IDX = 0;
  static const int ANGLE_MID_IDX = 1;
  static const int ANGLE_A2_IDX = 2;
  static const int ANGLE_KANGLE = 3;
  static const int ANGLE_EQANGLE = 4;
  static const int ANGLE_VARIABLE = 5;
  static const int ANGLE_DATA_SIZE = 6;

  Real*   size;
  Real    temperature;
  Real    cutoff;
  Real    maxTranslate;
  Real    maxRotate;

  int numSteps;

  int     numMolecules;

  int     numAtoms;
  int     numBonds;
  int     numAngles;
  int     numDihedrals;

  int**    moleculeData;
  int*     primaryIndexes;

  Real** rollBackCoordinates;

  Real**  atomCoordinates;
  Real**  atomData;

  Real*   bondLengths;
  const Real**  bondData;

  Real*   angleSizes;
  const Real**  angleData;

  Real*   dihedralSizes;
  const Real**  dihedralData;

  bool useNLC;

  NLC_Node** neighbors;
  NLC_Node**** neighborCells;
  NLC_Node* nlc_heap;
  NLC_Node* prevNode;
  int* prevCell;

  int* numCells;
  Real* cellWidth;

  /**
   * Roll back a molecule to its original poisition. Performs translation and
   * rotation.
   * Inputs: molIdx - The index of the molecule to change.
   */
  void rollback(int molIdx);

  /**
   * Adds the molecules from a Box to this simulation box.
   * Inputs: molecules - A dynamic array containing the box's molecules.
   *         numMolecules - The size of the molecules array.
   */
  bool addMolecules(Molecule* molecules, int numMolecules);

  /**
   * Adds the primary indexes from a box to the following simBox object.
   * Inputs: primaryAtomIndexArray - The box's array containing primary index
   *         definitions.
   */
  void addPrimaryIndexes(std::vector< std::vector<int> *>* primaryAtomIndexArray);

  /**
   * Keeps a molecule in the simulation box's boundaries, based on the location
   * of its first primary index atom.
   * Inputs: molIdx - The index of the molecule to keep inside the box.
   */
  void keepMoleculeInBox(int molIdx);

  /**
   * Given a distance, makes the distance periodic to mitigate distances greater
   * than half the length of the box.
   * Inputs: x - The measurement to make periodic.
   *         dimension - The dimension the measurement must be made periodic in.
   */
  Real makePeriodic(Real x, int dimension);

  /**
   * Prints the coordinates of every atom in every molecule in the simulation box.
   * NOTE: THIS IS FOR DEBUG ONLY.
   */
  void printCoordinates();

  /**
   * Given a box, constructs the corresponding simulation box.
   * Inputs: box - The box to construct the simulation box from.
   */
  void buildBox(Box* box);

  /**
   * Returns the index of a random molecule within the simulation box.
   */
  int chooseMolecule() const;

  /**
   * Given a molecule to change, randomly moves the molecule within the box.
   * This performs a translation and a rotation, and saves the old position.
   * Inputs: molIdx - The index of the molecule to change.
   */
  void changeMolecule(int molIdx);

  /**
   * Given a molecule to change, and information to change, translates and
   * rotates the specified molecule, saving its old position.
   * Inputs: molIdx - The index of the molecule to change.
   *         vIdx - The index of the pivot point of the molecule.
   *         dX - The amount to translate in the x direction.
   *         dY - The amount to translate in the y direction.
   *         dZ - The amount to translate in the z direction.
   *         rX - The amount to rotate around the x axis (degrees).
   *         rY - The amount to rotate around the y axis (degrees).
   *         rZ - The amount to rotate around the z axis (degrees).
   */
  void changeMolecule(int molIdx, int vIdx, Real dX, Real dY, Real dZ, Real rX, Real rY, Real rZ);

  /**
   * Given an atom to translate, and the directions to tranlate it, moves it.
   * Inputs: aIdx - The index of the atom to change.
   *         dX - The amount to translate in the x direction.
   *         dY - The amount to translate in the y direction.
   *         dZ - The amount to tranlsate in the z direction.
   */
  void translateAtom(int aIdx, Real dX, Real dY, Real dZ);

  /**
   * Given an atom to rotate, its pivot point, and the amounts to rotate,
   * rotates the atom around the pivot (rotations in degrees).
   * Inputs: aIdx - The index of the atom to change.
   *         pivotIdx - The index of the atom about which to rotate.
   *         rotX - The amount of rotation around the x-axis that is done.
   *         rotY - The amount of rotation around the y-axis that is done.
   *         rotZ - The amount of rotation around the z-axis that is done.
   */
  void rotateAtom(int aIdx, int pivotIdx, Real rotX, Real rotY, Real rotZ);

  /**
   * Given an atom and an amount to rotate, rotates about the x-axis.
   * Inputs: aIdx - The index of the atom to rotate.
   *         angleDeg - The angle to rotate it (in degrees).
   */
  void rotateX(int aIdx, Real angleDeg);

  /**
   * Given an atom and an amount to rotate, rotates it about the y-axis.
   * Inputs: aIdx - The index of the atom to rotate.
   *         angleDeg - The angle to rotate it (in degrees).
   */
  void rotateY(int aIdx, Real angleDeg);

  /**
   * Given an atom and an amount to rotate, rotates it about the z-axis.
   * Inputs: aIdx - The index of the atom to rotate.
   *         angleDeg - The angle to rotate it (in degrees).
   */
  void rotateZ(int aIdx, Real angleDeg);

  /**
   * calcSystemEnergy Determines the total energy of the box
   * Inputs:  box - The Simulation Box to calculate the energy of.
   *          subLJ - Initial Lennard - Jones energy.
   *          subCharge - Initial Coulomb energy.
   * Outputs: The total energy of the box (discounts initial lj / charge energy)
   *          subLJ - Final Lennard - Jones energy
   *          subCharge - Final Coulomb energy
   */
  Real calcSystemEnergy (Real &subLJ, Real &subCharge);

  /**
   * calcMolecularEnergyContribution Determines the energy contribution of a particular
   *                                 molecule.
   * Inputs:  box - The Simulation Box to calculate the energy of.
   *          subLJ - Initial Lennard - Jones energy.
   *          subCharge - Initial Coulomb energy.
   *          currMol - The index of the molecule to calculate the contribution of.
   *          startMol - The index of the molecule to begin searching from to determine
   *                     interaction energies.
   * Outputs: The total energy of the box (discounts initial lj / charge energy)
   *          subLJ - Final Lennard - Jones energy
   *          subCharge - Final Coulomb energy.
   */
  Real calcMolecularEnergyContribution(Real &subLJ, Real &subCharge, int currMol, int startMol);

  /**
   * moleculesInRange Determines whether or not two molecule's primaryIndexes are
   *                  within the cutoff range of one another.
   * Inputs:  box - The Simulation Box containing the molecules.
   *          p1Start - The index of the first primary index from molecule 1.
   *          p1End - index of the last primary index from molecule 1,  + 1.
   *          p2Start - index of the first primary index from molecule 2.
   *          p2End - index of the last primary index from molecule 2,  + 1.
   * Outputs: true if the molecules are in range, false otherwise.
   */
  bool moleculesInRange(refInt p1Start, refInt p1End, refInt p2Start, refInt p2End);

  /**
   * calcMoleculeInteractionEnergy Calcs the energy caused by the interaction between
   *                               a given pair of molecules.
   * Inputs:  box - The Simulation box containing the molecules.
   *          subLJ - The initial Lennard - Jones energy.
   *          subCharge - The initial Coloumb energy.
   *          m1 - The molecule index of the first molecule.
   *          m2 - The molecule index of the second molecule.
   * Outputs: The energy from the molecular interaction.
   *          subLJ - The Lennard - Jones energy total after the interaction.
   *          subCharge - The Coloumb energy total after the interaction.
   */
  Real calcMoleculeInteractionEnergy (Real &subLJ, Real &subCharge, refInt m1, refInt m2);


  /**
   * calcAtomDistSquared Calculates the square of the distance between two atoms.
   * Inputs: box - The Simulation box containing the atoms.
   *         a1 - The atom index of the first atom.
   *         a2 - The atom index of the second atom.
   * Outputs: The square of the distance between the atoms.
   */
  Real calcAtomDistSquared(refInt a1, refInt a2);

  /**
   * calcLJEnergy Calculates the Lennard - Jones potential between two atoms.
   * Inputs:  box - The Simulation box containing the atoms.
   *          a1 - The atom index of the first atom.
   *          a2 - The atom index of the second atom.
   *          r2 - The distance between the atoms, squared.
   * Outputs: The Lennard - Jones potential from the two atoms' interaction.
   */
  Real calcLJEnergy(refInt a1, refInt a2, const Real& r2);

  /**
   * calcChargeEnergy Calculates the Coloumb potential between two atoms.
   * Inputs:  box - The Simulation box containing the atoms.
   *          a1 - The atom index of the first atom.
   *          a2 - The atom index of the second atom.
   *          r - The distance between the atoms.
   * Outputs: The Coloumb potential from the two atoms' interaction.
   */
  Real calcChargeEnergy(refInt a1, refInt a2, const Real& r);

  /**
   * calcBlending Calculates the geometric mean of two real numbers.
   * Inputs:  a - The first real number.
   *          b - The second real number.
   * Outputs: sqrt(|a*b|), the geometric mean of the two numbers.
   */
  Real calcBlending (const Real &a, const Real &b);

  /**
   * Stores all of the molecules in this simulation box in an appropriate chain
   * of neighbor-list linked cells.
   */
  void fillNLC();

  /**
   * Given the index for a molecule, popuplate neighbors with the heads of
   * linked-lists holding the molecules in neighboring cells.
   */
  int findNeighbors(int molIdx);

  /**
   * Given an index and a dimension, returns the index, wrapped around the box.
   * Inputs:  idx - The index of a cell to wrap.
   *          dimension - The dimension the index is in.
   * Outputs: The index, wrapped around the size of the box.
   */
  int wrapCell(int idx, int dimension);

  /**
   * Given a location and a dimension, returns the corresponding cell index.
   * Inputs:  location - The index of the cell to wrap.
   *          dimension - The dimension to wrap in.
   * Outputs: The index, wrapped around the size of the box.
   */
  int getCell(Real loc, int dimension);

  /**
   * Given a molecule, updates the molecule's position in the neighbor linked
   * cells (if it needs to be updated).
   * Inputs: molIdx - The index of the moleucle to update.
   */
  void updateNLC(int molIdx);
};

typedef const SimBox & refBox;

#endif
