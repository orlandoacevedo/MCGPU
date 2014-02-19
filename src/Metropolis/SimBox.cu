/*
	New version of SimBox

	Author: Nathan Coleman
	Last Changed: February 19, 2014
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "SimBox.cuh"

using namespace std;

double randomFloat(const double start, const double end){}

//Constructor & Destructor
SimBox::SimBox(Config_Scan configScan){}
SimBox::~SimBox(){}

//Utility
void assignAtomPositions(double x, double y, double z, Molecule *molecule, Environment *environment){}

int changeMolecule(int moleculeIndex){}

int chooseMolecule(){}

void generateFCCBox(Molecule *molecules, Environment *environment){}

void generatePoints(Molecule *molecules, Environment *environment){}

double getFValue(Atom *atom1, Atom *atom2, Molecule *moelecules, Environment *environment){}

int getXFromIndex(int index){}

int getYFromIndex(int index){}

int rollBack(int moleculeIndex){}