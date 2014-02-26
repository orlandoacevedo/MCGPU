/*
	New version of SimBox

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

//#include "Utilities/Opls_Scan.h"
//#include "Utilities/Zmatrix_Scan.h"
#include "Metropolis/Box.h"
#include "SerialBox.h"

using namespace std;

double randomFloat(const double start, const double end){return 0.0;}

//Constructor & Destructor
SerialBox::SerialBox(Config_Scan configScan):Box()
{
	environment = (Environment*)malloc(sizeof(Environment));
	memcpy(environment, configScan.getEnviro(), sizeof(Environment));
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