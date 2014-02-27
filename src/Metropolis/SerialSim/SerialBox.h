/*
	New version of SimBox
	Minimized to include only Atoms and Molecules

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

#ifndef SERIALBOX_H
#define SERIALBOX_H

#include "Metropolis/Utilities/ConfigScanTemp.h"
#include "Metropolis/Box.h"

class SerialBox : Box
{
	public:
		//Constructor & Destructor
		SerialBox(Config_Scan configScan);
		~SerialBox();
};

#endif