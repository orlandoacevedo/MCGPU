/*
	Represents a simulation box, holding environment and molecule data.
	Subclass of Box.

	Author: Nathan Coleman
	Created: February 21, 2014

	-> February 26, by Albert Wallace
	-> March 28, by Joshua Mosby
	-> April 21, by Nathan Coleman
*/

#ifndef SERIALBOX_H
#define SERIALBOX_H

#include "Metropolis/Box.h"

class SerialBox : public Box
{
	public:
		SerialBox();
		~SerialBox();

		int molecTypenum;
		Table *tables;
};

#endif
