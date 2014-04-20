/// @file SerialBox.cpp
///
/// Represents a simulation box, holding environment and molecule data.
///   Subclass of Box.

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