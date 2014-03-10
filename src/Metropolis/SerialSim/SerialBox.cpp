/*
	New version of SimBox

	Author: Nathan Coleman
	Last Changed: February 21, 2014
	
	-> February 26, 27, by Albert Wallace
*/

#include <string.h>
#include <stdlib.h>
#include "SerialBox.h"

double randomFloat(const double start, const double end){return 0.0;}


SerialBox::SerialBox(IOUtilities configScan) : Box()
{
	environment = new Environment();
	//memcpy(environment, configScan.currentEnvironment, sizeof(Environment));
}

SerialBox::~SerialBox()
{
	FREE(atoms);
	FREE(environment);
	FREE(molecules);
	FREE(energies);
}