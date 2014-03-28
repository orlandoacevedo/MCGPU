/*
	New version of SimBox
	Minimized to include only Atoms and Molecules

	Author: Nathan Coleman
	Last Changed: February 21, 2014
	-> February 26, by Albert Wallace
*/

#ifndef SERIALBOX_H
#define SERIALBOX_H

#include "Metropolis/Utilities/IOUtilities.h"
#include "Metropolis/Box.h"

class SerialBox : public Box
{
	public:
		SerialBox(IOUtilities ioUtil);
		~SerialBox();

		int molecTypenum;
		Table *tables;
};

#endif