/*
	SuperClass to the SerialBox and ParallelBox classes.

	Author: Nathan Coleman, Tavis Maclellan
	Last Changed by author: February 26, 2014
	
	-> other changes: February 26 (after authors), by Albert Wallace
*/

#ifndef BOX_H
#define BOX_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "Utilities/MathLibrary.h"
#include "Utilities/Config_Scan.h"
#include "Metropolis/DataTypes.h"
#include "Utilities/StructLibrary.h"

#define FREE(ptr) if(ptr!=NULL) { free(ptr);ptr=NULL;}

using namespace std;

class Box
{
	protected:
		Atom *atoms;
		Environment *environment;
		Molecule *molecules;
		Bond *bonds;
		Angle *angles;
		Dihedral *dihedrals;
		Hop *hops;
		Real *energies;
		int atomCount, moleculeCount, energyCount, maxMolSize;

	public:
		Box(IOUtilities ioUtil);
		~Box();
		Atom *getAtoms(){return atoms;};
		Environment *getEnvironment(){return environment;};
		Molecule *getMolecules(){return molecules;};
		Real *getEnergies(){return energies;};
		int getAtomCount(){return atomCount;};
		int getMoleculeCount(){return moleculeCount;};
		int getEnergyCount(){return energyCount;};
		int getMaxMolSize(){return maxMolSize;};
		int chooseMolecule();
		int changeMolecule(int molIdx);
		void keepMoleculeInBox(Molecule *molecule, Environment *enviro);
		int Rollback(int moleno);
		int saveChangedMole(int moleno);
		int copyMolecule(Molecule *mole_dst, Molecule *mole_src);
		
};

#endif