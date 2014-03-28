/*
	New version of SimBox
	Minimized to include only Atoms and Molecules

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

#ifndef PARALLELBOX_H
#define PARALLELBOX_H

#include "Metropolis/Box.h"
#include "Metropolis/DataTypes.h"

//DeviceMolecule struct needs to be moved to same location as other structs
class ParallelBox : Box
{
	private:
		Atom *atomsD;
		Environment *environmentD;
		Molecule *moleculesD, *transferMoleculesH;
		int *nbrMolsH, *nbrMolsD, *molBatchH, *molBatchD;
		Real *energiesD;
		
		void copyDataToDevice();	
		void writeChangeToDevice(int changeIdx);

	public:
		ParallelBox(IOUtilities ioUtil);
		~ParallelBox();
		virtual int changeMolecule(int molIdx);
		virtual int rollback(int moleno);
};

#endif