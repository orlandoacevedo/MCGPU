/// @file SerialBox.h
///
/// Represents a simulation box, holding environment and molecule data.
///   Subclass of Box.

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