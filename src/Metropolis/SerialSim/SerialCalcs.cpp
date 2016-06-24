/**
 *	Contains the methods required to calculate energies serially.
 *
 *	Created: February 21, 2014
 *
 * -> February 26, by Albert Wallace
 *	-> March 28, by Joshua Mosby
 *	-> April 21, by Nathan Coleman
 *	-> February 25, 2015 by Jared Brown
 */

#include <math.h>
#include "Metropolis/DataTypes.h"
#include "Metropolis/SimulationArgs.h"
#include "Metropolis/Utilities/FileUtilities.h"
#include "SerialCalcs.h"

using namespace std;

Box* SerialCalcs::createBox(SimulationArgs& simArgs, long* startStep,
                            long* steps, SBScanner* sbScanner) {
	SerialBox* box = new SerialBox();
	if (!loadBoxData(simArgs, box, startStep, steps, sbScanner)) {
		if (simArgs.fileType != InputFile::Unknown) {
			std::cerr << "Error: Could not build from file: " << simArgs.filePath
								<< std::endl;
		} else {
			std::cerr << "Error: Can not build environment with unknown file: " 
								<< simArgs.filePath << std::endl;
		}
		return NULL;
	}

	return (Box*) box;
}

Real SerialCalcs::calcEnergy_LRC(Box* box) {
	Molecule *molecules = box->getMolecules();
	Environment *enviro = box->getEnvironment();

	Real Ecut = 0.0;		// Holds LJ long-range cutoff energy correction

	Real Vnew = enviro->x * enviro->y * enviro->z;	// Volume of box in Ang^3
	Real RC3 = 1.00 / pow(enviro->cutoff, 3);		// 1 / cutoff^3
	Real RC9 = pow(RC3, 3);							// 1 / cutoff^9

	// Note: currently only supports at most TWO solvents (needs to be updated for more)
	int a = 0, b = 1;
	Real NMOL1 = enviro->numOfMolecules / 2;	// Number of molecules of solvent1
	Real NMOL2 = enviro->numOfMolecules / 2;	// Number of molecules of solvent2
	int NATOM1 = molecules[a].numOfAtoms;			// Number of atoms in solvent1
	int NATOM2 = molecules[b].numOfAtoms;			// Number of atoms in solvent2
	int NATMX = NATOM1;
	if (NATMX < NATOM2) {		// NATMX = MAX(NAT0M1, NAT0M2)
		NATMX = NATOM2;
	}

	Real sig2, sig6, sig12;
	// get LJ-values for solvent1 and store in A6, A12
	Real SigmaA[NATOM1], EpsilonA[NATOM1];
	Real A6[NATOM1], A12[NATOM1];
	for(int i = 0; i < NATOM1; i++) {
		if (molecules[a].atoms[i].sigma < 0 || molecules[a].atoms[i].epsilon < 0) {
			SigmaA[i] = 0.0;
			EpsilonA[i] = 0.0;
		} else {
			SigmaA[i] = molecules[a].atoms[i].sigma;
			EpsilonA[i] = molecules[a].atoms[i].epsilon;
		}

		sig2 = pow(SigmaA[i], 2);
    sig6 = pow(sig2, 3);
    sig12 = pow(sig6, 2);
		A6[i] = sqrt(4 * EpsilonA[i] * sig6);
		A12[i] = sqrt(4 * EpsilonA[i] * sig12);
	}

	// get LJ-values for solvent2 and store in B6, B12
	Real SigmaB[NATOM2], EpsilonB[NATOM2];
	Real B6[NATOM2], B12[NATOM2];
	for(int i = 0; i < NATOM2; i++) {
		if (molecules[b].atoms[i].sigma < 0 || molecules[b].atoms[i].epsilon < 0) {
			SigmaB[i] = 0.0;
			EpsilonB[i] = 0.0;
		} else {
			SigmaB[i] = molecules[b].atoms[i].sigma;
			EpsilonB[i] = molecules[b].atoms[i].epsilon;
		}

		sig2 = pow(SigmaB[i], 2);
    sig6 = pow(sig2, 3);
    sig12 = pow(sig6, 2);
		B6[i] = sqrt(4 * EpsilonB[i] * sig6);
		B12[i] = sqrt(4 * EpsilonB[i] * sig12);
	}

	// loop over all atoms in a pair
	for(int i = 0; i < NATOM1; i++) {
		for(int j = 0; j < NATOM2; j++) {
			Ecut += (2*PI*NMOL1*NMOL1/(3.0*Vnew)) * (A12[i]*A12[j]*RC9/3.0 - A6[i]*A6[j]*RC3);
			Ecut += (2*PI*NMOL2*NMOL2/(3.0*Vnew)) * (B12[i]*B12[j]*RC9/3.0 - B6[i]*B6[j]*RC3);
			Ecut += (4*PI*NMOL1*NMOL2/(3.0*Vnew)) * (A12[i]*B12[j]*RC9/3.0 - A6[i]*B6[j]*RC3);
		}
	}

	//std::cout << "Energy_LRC = " << Ecut << std::endl;
	return Ecut;
}
