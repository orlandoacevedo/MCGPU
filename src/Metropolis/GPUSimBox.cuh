/*
	New version of SimBox
	Minimized to include only Atoms and Molecules

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

#ifndef GPUSIMBOX_H
#define GPUSIMBOX_H

#include "../Utilities/Opls_Scan.h"
#include "../Utilities/Config_Scan.h"
#include "../Utilities/metroUtil.h"
#include "../Utilities/Zmatrix_Scan.h"
#include "../Utilities/State_Scan.h"
#include "SimBox.h"

//DeviceMolecule struct needs to be moved to same location as other structs

class GPUSimBox : SimBox
{
	private:
		size_t atomSize;
		size_t moleculeSize;
		size_t environmentSize;

	public:
		//Constructor & Destructor
		GPUSimBox(Config_Scan configScan);
		~GPUSimBox();

		//Utility
		int copyBoxToHost();
		int copyBoxToDevice();
};

//Cuda Necessities
void cudaAssert(const cudaError err, const char *file, const int line);
#define cudaErrorCheck(call) { cudaAssert(call,__FILE__,__LINE__); }
#define cudaFREE(ptr) if(ptr!=NULL) { cudaFree(ptr);ptr=NULL;}

#endif