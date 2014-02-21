/*
	New version of SimBox

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "SimBox.cuh"

using namespace std;

double randomFloat(const double start, const double end){return 0.0;}

//Constructor & Destructor
SimBox::SimBox(Config_Scan configScan){}
SimBox::~SimBox(){}

//Utility
void assignAtomPositions(double x, double y, double z, Molecule *molecule, Environment *environment){}

int changeMolecule(int moleculeIndex){return 0;}

int chooseMolecule(){return 0;}

void generateFCCBox(Molecule *molecules, Environment *environment){}

void generatePoints(Molecule *molecules, Environment *environment){}

double getFValue(Atom *atom1, Atom *atom2, Molecule *moelecules, Environment *environment){return 0.0;}

int getXFromIndex(int index){return 0;}

int getYFromIndex(int index){return 0;}

int rollBack(int moleculeIndex){return 0;}