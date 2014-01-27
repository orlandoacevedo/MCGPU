/*!\file
  \Class for simulation caculate, change state of SimBox.
  \author David(Xiao Zhang) and Orlando Acevedo.
 
  This file contains functions that are used to handle enviroments and common function
  for box.
 */
#ifndef LINEARSIM_H
#define LINEARSIM_H 

#include "SimBox.h"
#include "../../Utilities/src/Opls_Scan.h"
#include "../../Utilities/src/Config_Scan.h"
#include "../../Utilities/src/metroUtil.h"
#include "../../Utilities/src/Zmatrix_Scan.h"
#include "../../Utilities/src/State_Scan.h"

#define THREADS_PER_BLOCK 128
#define PI 3.14159265
// Linked-cell neighbor list constants
#define NMAX 100000  /* Maximum number of atoms which can be simulated */
#define NCLMAX 10000 /* Maximum number of linked-list cells */
#define EMPTY -1

// boltzman constant in kcal mol-1 K-1
const double kBoltz = 0.00198717;

class LinearSim
{
	private:
	 	SimBox *box;
	 	int steps;
	 	float currentEnergy;
	  float oldEnergy;
	  int accepted;
	  int rejected;
	  
	public:
		float getcurrentEnergy(){return currentEnergy;}; 	
	 	int getaccepted() {return accepted;};
	 	int getrejected() {return rejected;};
	public:
	 	LinearSim(SimBox *initbox,int initsteps);
	 	double calc_lj(Atom atom1, Atom atom2, Environment enviro);
		double calcEnergyWrapper(Molecule *molecules, Environment *enviro);
		void calcEnergy(Atom *atoms, Environment *enviro, double *energySum);
		double calcCharge(Atom atom1, Atom atom2, Environment *enviro);
		double calcEnergy_NLC(Molecule *molecules, Environment *enviro);
		double calcIntramolEnergy_NLC(Environment *enviro, Molecule *molecules);
		double calcBlending(double d1, double d2);
		double Energy_LRC(Molecule *molec, Environment *enviro);
		void runLinear(int steps);
 	
};
 	
#endif