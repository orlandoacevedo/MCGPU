/**
 * Represents a simulation box, holding environment and molecule data.
 *
 * Superclass to SerialBox. This class serves as a temporary intermediate class
 * between constants files, z-matrices, config files, and state files, and the
 * SimBox class, which now performs energy calculations.
 *
 *	Author: Nathan Coleman
 *	Created: February 21, 2014
 *
 *	-> February 26, by Albert Wallace
 *	-> March 28, by Joshua Mosby
 *	-> April 21, by Nathan Coleman
 */

#ifndef BOX_H
#define BOX_H

#include <string.h>
#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
#include "Utilities/MathLibrary.h"
#include "Metropolis/DataTypes.h"
#include "Utilities/StructLibrary.h"
#include "SerialSim/NeighborList.h"

#define FREE(ptr) if(ptr!=NULL) { free(ptr);ptr=NULL;}

using namespace std;

class Box
{
	private:

		/**
		 * IN_BOX indicates that a coordinate measurement is between 0 and the size
		 *     of the box.
		 */
		static const int IN_BOX = 0;

		/**
		 * BELOW_ZERO indicates that a coordinate measurement is less than 0.
		 */
		static const int BELOW_ZERO = -1;

		/**
		 * ABOVE_BOX_DIM indicates that a coordinate measurement is greater than the
		 *     length of the box in the corresponding dimension.
		 */
		static const int ABOVE_BOX_DIM = 1;

	protected:

		/**
		 * Holds the old state of whatever molecule is being changed so that it can
		 * be rolled back if the change is rejected.
		 */
		Molecule changedMol;

    /**
		 * Given a coordinate and the length of the box in the corresponding
		 *     dimension, indicates whether the coordinate is cout of bounds.
     *
		 * @param coor The coordinate measurement to test for being in bounds.
		 * @param boxDim The length of the box in the same dimension as the
		 *     coordinate measurement.
		 * @returns IN_BOX if 0 <= coor <= boxDim, BELOW_ZERO if coor < 0, and
		 *     ABOVE_BOX_DIM if boxDim < coor.
		 */
		int isOutOfBounds(Real coor, Real boxDim);

	public:

		/**
		 * Pointer to the environment struct, as specified in StructLibrary.h.
		 *     Environment contains basic parameters, such as the box dimensions and
		 *     the temperature of the box.
		 */
		Environment *environment;

		/**
		 * A dynamic array of the atoms within the box.
		 */
		Atom *atoms;

		/**
		 * A dynamic array of the molecules within the box. The atoms in each
		 *     molecule and the atoms in atoms are aliased to one another, along
		 *     with bonds, angles, etc.
		 */
		Molecule *molecules;

		/**
		 * A dynamic array of the bonds within the box.
		 */
		Bond *bonds;

		/**
		 * A dynamic array of the angles within the box.
		 */
		Angle *angles;

		/**
		 * A dynamic array of the dihedrals within the box.
		 */
		Dihedral *dihedrals;

		/**
		 * A dynamic array of the hops within the box.
		 */
		Hop *hops;

		/**
		 * Points to the neighborList for the box.
		 */
		NeighborList *neighborList;

		/**
		 * The number of atoms in the box.
		 */
		int atomCount;

		/**
		 * The number of molecules in the box.
		 */
		int moleculeCount;

		/**
		 * The number of bonds in the box.
		 */
		int bondCount;

		/**
		 * The number of angles in the box.
		 */
		int angleCount;

		/**
		 * The number of dihedrals in the box.
		 */
		int dihedralCount;

		/**
		 * The number of hops in the box.
		 */
		int hopCount;

		/**
		 * Default (parameterless) constructor for box.
		 */
		Box();

		/**
		 * Destructor for box.
		 */
		~Box();

		/**
		 * Getter method for atoms.
		 * @return A pointer to the atoms dynamic array.
		 */
		Atom* getAtoms(){return atoms;};

		/**
		 * Getter method for atomCount
		 * @return The number of atoms in the box.
		 */
		int getAtomCount(){return atomCount;};

		/**
		 * Getter method for molecules.
		 * @return A pointer to the molecules dynamic array.
		 */
		Molecule* getMolecules(){return molecules;};

		/**
		 * Getter method for molecule count.
		 * @return The number of molecules in the box.
		 */
		int getMoleculeCount(){return moleculeCount;};

		/**
		 * Getter method for environment.
		 * @return A pointer to the box's environment.
		 */
		Environment* getEnvironment(){return environment;};

		/**
		 * Getter method for neighborList.
		 * @return A pointer to the box's neighbor list.
		 */
		NeighborList* getNeighborList() {return neighborList;};

		/**
		 * Creates the box's neighbor list.
		 */
		void createNeighborList();

    /**
		 * Chooses a random molecule to be changed for a given simulation step.
		 *
		 * @return Returns the index of the chosen molecule.
     */
		int chooseMolecule();

    /**
		 * Changes a given molecule (specifically its Atoms) in a way determined
		 * by given input tranlations and rotations.
		 *
		 * @param molIdx The index of the molecule to be changed.
		 * @param vIdx The index of the atom in the molecule about which it will
		 * 		 be rotated.
		 * @param dX The amount to translate the molecule in the x direction.
		 * @param dY The amount to translate the molecule in the y direction.
		 * @param dZ The amount to translate the molecule in the z direction.
		 * @param rX The amount to rotate the molecule about the x axis.
		 * @param rY The amount to rotate the molecule about the y axis.
		 * @param rZ The amount to rotate the molecule about the z axis.
		 * @return The index of the changed molecule.
		 * @deprecated This method was used solely for testing the new SimBox
		 */
		virtual int changeMolecule(int molIdx, int vIdx, Real dX, Real dY, Real dZ, Real rX, Real rY, Real rZ);

    /**
		 * Changes a given molecule (specifically its Atoms)
		 *   in a random way, constrained by maxTranslation
		 *   and maxRotation.
		 * @param molIdx The index of the molecule to be changed.
		 * @return Returns the index of the changed molecule.
		 * @note This method is virtual to be overridden by an subclass.
     */
		virtual int changeMolecule(int molIdx);

		/**
		 * Makes each of the molecule's positional attributes
		 *   periodic within the dimensions of the environment.
		 * @param molecule The index of the molecule to be fixed.
		 */
		void keepMoleculeInBox(int molIdx);

		/**
		 * Rolls back the previous molecule change.
		 * @param molIdx The index of the molecule to be reverted.
		 * @return Returns the index of the reverted molecule.
		 * @note This method is virtual to be overridden by an subclass.
		 */
		virtual int rollback(int molIdx);

		/**
		 * Saves the unchanged version of a molecule to be changed.
		 * @param molIdx The index of the molecule to be saved.
		 */
		void saveChangedMol(int molIdx);

		/**
		 * Copies the data of one molecule to another.
		 * @param mol_dst A pointer to the destination molecule.
		 * @param mol_src A pointer to the source molecule.
		 */
		void copyMolecule(Molecule *mol_dst, Molecule *mol_src);

		/**
		 * Makes a position periodic within a specified range.
		 * @param x The position to be made periodic.
		 * @param boxDim The magnitude of the periodic range.
		 * @return Returns the periodic position.
		 */
		Real wrapBox(Real x, Real boxDim, int position);

};

#endif
