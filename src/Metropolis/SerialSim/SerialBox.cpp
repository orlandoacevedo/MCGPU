/*
	New version of SimBox

	Author: Nathan Coleman
	Last Changed: February 21, 2014
	
	-> February 26, 27, by Albert Wallace
*/

#include "SerialBox.h"

using namespace std;

SerialBox::SerialBox() : Box()
{
	
}

SerialBox::~SerialBox()
{
	FREE(angles);
	FREE(atoms);
	FREE(bonds);
	FREE(dihedrals);
	FREE(environment);
	FREE(hops);
	FREE(molecules);
	FREE(energies);
}