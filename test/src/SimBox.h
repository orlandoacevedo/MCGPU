/*!\file
  \Class for simulation Box, including Enviroments and points to molocoles.
  \author David(Xiao Zhang).
 
  This file contains functions that are used to handle enviroments and common function
  for box.
 */
#ifndef SIMBOX_H
#define SIMBOX_H 


#include "../../Utilities/src/Opls_Scan.h"
#include "../../Utilities/src/Config_Scan.h"
#include "../../Utilities/src/metroUtil.h"
#include "../../Utilities/src/Zmatrix_Scan.h"
#include "../../Utilities/src/State_Scan.h"

class SimBox {
private:
 	Molecule *molecules;
 	Environment *enviro;
 	
public:
 	SimBox(char const* ConfigFile);
 	Molecule *GetMolecule();
 	Environment *GetEnviro();
 	int ReadStateFile(char const* StateFile);
 	int WriteStateFile(char const* StateFile); 	
};
 	
#endif