/*
	Represents a simulation box, holding environment and molecule data.
	Subclass of Box.

	Author: Nathan Coleman
	Created: February 21, 2014

	-> February 26, by Albert Wallace
	-> March 28, by Joshua Mosby
	-> April 21, by Nathan Coleman
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
}
