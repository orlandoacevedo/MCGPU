/*
	New version of SimBox

	Author: Nathan Coleman
	Last Changed: February 21, 2014
	
	-> February 26, 27, by Albert Wallace
*/

#ifndef SERIALBOX_CPP
#define SERIALBOX_CPP

//#include "Utilities/Opls_Scan.h"
//#include "Utilities/Zmatrix_Scan.h"
//#include "Metropolis/Box.h" AlbertIncludes
#include "SerialBox.h"
#include <string.h>
#include <stdlib.h>

using namespace std;

double randomFloat(const double start, const double end){return 0.0;}

//Constructor & Destructor
SerialBox::SerialBox(IOUtilities configScan):Box()
{
	//environment = (Environment*)malloc(sizeof(Environment));
	environment = new Environment();
	memcpy(environment, configScan.filePathsEtc->currentEnvironment, sizeof(Environment));
	// string oplsPath = configScan.getOplsusaparPath();
	// Opls_Scan oplsScan (oplsPath);
	// oplsScan.scanInOpls(oplsPath);
	// Zmatrix_Scan zMatrixScan (configScan.getZmatrixPath(), &oplsScan);
	// if (zMatrixScan.scanInZmatrix() == -1)
	// {
	// 	fprintf(stderr,"Could not open %s\n", configScan.getZmatrixPath());
	// 	exit(1);
	// }

}

SerialBox::~SerialBox()
{
	FREE(atoms);
	FREE(environment);
	FREE(molecules);
	FREE(energies);
}

#endif