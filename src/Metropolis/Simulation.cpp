/*
	New Simulation to replace linearSim and parallelSim

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

#include "Simulation.cuh"

using namespace std;

//Constructor & Destructor
Simulation::Simulation(bool useGPU, char const* configPath, int interval){}
Simulation::~Simulation(){}

//Utility
double calcBlending(double d1, double d2){return 0.0;}
double calcCharge(double charge1, double charge2m double r){return 0.0;}
double calcInteratomicEnergy(Molecule *molecules,int molecule1, int molecule2, Environment *environment){return 0.0;}
double calcIntermolecularEnergy(Molecule *molecules, int molecule1, int molecule2, Environment *environment){return 0.0;}
double calcIntramolecularEnergyNLC(Molecule *molecules, Environment *environment){return 0.0;}
double calcLennardJones(Atom atom1, Atom atom2, double r2){return 0.0;}
double calcMolecularEnergyContribution(Molecule *molecules, Environment *environment, int currentMolecule, int startIndex){return 0.0;}
double calcSystemEnergy(Molecule *molecules, Environment *environment){return 0.0;}
void run(){}