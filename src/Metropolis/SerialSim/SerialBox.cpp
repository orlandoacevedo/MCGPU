/*
	New version of SimBox

	Author: Nathan Coleman
	Last Changed: February 21, 2014
	
	-> February 26, 27, by Albert Wallace
*/

#include "SerialBox.h"

using namespace std;

Real randomReal(const Real start, const Real end)
{
	return 0.0;
}


SerialBox::SerialBox(IOUtilities ioUtil) : Box(ioUtil)
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