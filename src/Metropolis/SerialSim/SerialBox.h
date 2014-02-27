/*
	New version of SimBox
	Minimized to include only Atoms and Molecules

	Author: Nathan Coleman
	Last Changed: February 21, 2014
	-> February 26, by Albert Wallace
*/

#ifndef SERIALBOX_H
#define SERIALBOX_H

//#include "Metropolis/Utilities/ConfigScanTemp.h" //AlbertIncludes
//#include "Metropolis/Box.h" //AlbertIncludes

//AlbertIncludes...
#include "../Utilities/IOUtilities.cpp"
#include "../Box.h"
//end AlbertIncludes

class SerialBox : Box
{
	public:
		//Constructor & Destructor
		SerialBox(IOUtilities configScan);
		~SerialBox();
};

#endif