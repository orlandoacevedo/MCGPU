/*
	New version of SimBox
	Minimized to include only Atoms and Molecules

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

#ifndef SERIALBOX_H
#define SERIALBOX_H

#include "Utilities/Opls_Scan.h"
#include "Utilities/Config_Scan.h"
#include "Utilities/metroUtil.h"
#include "Utilities/Zmatrix_Scan.h"
#include "Utilities/State_Scan.h"
#include "Metropolis/Utilities/IOUtilities.cuh"
#include "Metropolis/Box.h"
#include "SerialCalcs.h"

class SimBox : Box
{
	public:
		//Constructor & Destructor
		SerialBox(Config_Scan configScan);
		~SerialBox();
};

double randomFloat(const double start, const double end);

#endif