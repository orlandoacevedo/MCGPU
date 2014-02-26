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
#include "SerialBox.cuh"

using namespace std;

double randomFloat(const double start, const double end){return 0.0;}

//Constructor & Destructor
SimBox::SimBox(Config_Scan configScan)
{
	molecules = NULL;
	environment = NULL;

}

SimBox::~SimBox()
{
	FREE(molecules);
	FREE(environment);
	FREE(atomPool);
	FREE(changedMolecule);
}