/*
	Represents a simulation box, holding environment and molecule data.
	Superclass to SerialBox and ParallelBox.

	Author: Nathan Coleman
	Created: February 21, 2014
	
	-> February 26, by Albert Wallace
	-> March 28, by Joshua Mosby
	-> April 21, by Nathan Coleman
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
		bool first;
	protected:
		Molecule changedMol;

	public:
		Environment *environment;
		Atom *atoms;
		Molecule *molecules;
		Bond *bonds;
		Angle *angles;
		Dihedral *dihedrals;
		Hop *hops;
		NeighborList *neighborList;
		int atomCount, moleculeCount, bondCount, angleCount, dihedralCount, hopCount;
		
		Box();
		~Box();
		Atom *getAtoms(){return atoms;};
		int getAtomCount(){return atomCount;};
		Molecule *getMolecules(){return molecules;};
		int getMoleculeCount(){return moleculeCount;};
		Environment *getEnvironment(){return environment;};
		
		/// Chooses a random molecule to be changed for a given
		///   simulation step.
		/// @return Returns the index of the chosen molecule.
		int chooseMolecule();
		
		/// Changes a given molecule (specifically its Atoms)
		///   in a random way, constrained by maxTranslation
		///   and maxRotation.
		/// @param molIdx The index of the molecule to be changed.
		/// @return Returns the index of the changed molecule.
		/// @note This method is virtual to be overridden by an subclass.
		virtual int changeMolecule(int molIdx);
		
		/// Makes each of the molecule's positional attributes
		///   periodic within the dimensions of the environment.
		/// @param molecule The index of the molecule to be fixed.
		void keepMoleculeInBox(int molIdx);
		
		/// Rolls back the previous molecule change.
		/// @param molIdx The index of the molecule to be reverted.
		/// @return Returns the index of the reverted molecule.
		/// @note This method is virtual to be overridden by an subclass.
		virtual int rollback(int molIdx);
		
		/// Saves the unchanged version of a molecule to be changed.
		/// @param molIdx The index of the molecule to be saved.
		void saveChangedMol(int molIdx);
		
		/// Copies the data of one molecule to another.
		/// @param mol_dst A pointer to the destination molecule.
		/// @param mol_src A pointer to the source molecule.
		void copyMolecule(Molecule *mol_dst, Molecule *mol_src);

		/// Makes a position periodic within a specified range.
		/// @param x The position to be made periodic.
		/// @param boxDim The magnitude of the periodic range.
		/// @return Returns the periodic position.
	 	Real wrapBox(Real x, Real boxDim);
		
};

#endif
