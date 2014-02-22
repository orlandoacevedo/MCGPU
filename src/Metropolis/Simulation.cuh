/*
	New Simulation to replace linearSim and parallelSim

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include "Utilities/Opls_Scan.h"
#include "Utilities/Config_Scan.h"
#include "Utilities/metroUtil.h"
#include "Utilities/Zmatrix_Scan.h"
#include "Utilities/State_Scan.h"
#include "SimBox.h"

#define THREADS_PER_BLOCK 128
#define PI 3.14159265
// Linked-cell neighbor list constants
#define NMAX 100000  /* Maximum number of atoms which can be simulated */
#define NCLMAX 10000 /* Maximum number of linked-list cells */
#define EMPTY -1

const double kBoltz = 0.00198717;

class Simulation
{
	private:
		int steps;
		float currentEnergy;
		float oldEnergy;
		int accepted;
		int rejected;
	
	public:
		//Constructor & Destructor
		//Take in bool for parallel or serial and filepath for config file
		Simulation(bool useGPU, char const* configPath,int interval);
		~Simulation();

		//Getters
		int getAccepted(){return accepted;};
		float getCurrentEnergy(){return currentEnergy;};
		int getRejected(){return rejected;};

		//Utility
		double calcBlending(double d1, double d2);
		double calcCharge(double charge1, double charge2m double r);
		double calcInteratomicEnergy(Molecule *molecules,int molecule1, int molecule2, Environment *environment);
		double calcIntermolecularEnergy(Molecule *molecules, int molecule1, int molecule2, Environment *environment);
		double calcIntramolecularEnergyNLC(Molecule *molecules, Environment *environment);
		double calcLennardJones(Atom atom1, Atom atom2, double r2);
		double calcMolecularEnergyContribution(Molecule *molecules, Environment *environment, int currentMolecule, int startIndex);
		double calcSystemEnergy(Molecule *molecules, Environment *environment);
		void run();
};

#endif