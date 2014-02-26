/*
	New version of SimBox
	Minimized to include only Atoms and Molecules

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

#ifndef SIMBOX_H
#define SIMBOX_H

#include "Utilities/Opls_Scan.h"
#include "Utilities/Config_Scan.h"
#include "Utilities/metroUtil.h"
#include "Utilities/Zmatrix_Scan.h"
#include "Utilities/State_Scan.h"
#include "Metropolis/Utilities/IOUtilities.cuh"

class SimBox
{
	private:
		Atom *atomPool;
		Environment *environment;
		Molecule changedMolecule;
		Molecule *molecules;

		//Utility
		int copyMolecule(Molecule *destination, Molecule *source);
		int saveMolecule(int moleculeIndex);

	public:
		int moleculeType;

		//Constructor & Destructor
		SimBox(Config_Scan configScan);
		int initGPUSimBox();
		~SimBox();

		//Getters
		Atom *getAtoms(){return atomPool;};
		Environment *getEnvironment(){return environment;};
		Molecule *getMolecules(){return molecules;};

		//IO functions
		//Being moved to Utilities directory
};

double randomFloat(const double start, const double end);

#endif