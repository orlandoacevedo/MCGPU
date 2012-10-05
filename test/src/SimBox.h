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
 	int getXFromIndex(int idx);
 	int getYFromIndex(int x, int idx);
 	double makePeriodic(double x, double box);
 	double wrapBox(double x, double box);
 	void keepMoleculeInBox(Molecule *molecule, Environment *enviro);
 	int hopGE3(int atom1, int atom2, Molecule *molecule);
 	Molecule* getMoleculeFromAtomID(Atom *a1, Molecule *molecules, Environment *enviro);
 	double getFValue(Atom *atom1, Atom *atom2, Molecule *molecules, Environment *enviro);
 	void generatefccBox(Molecule *molecules, Environment *enviro);
 	void generatePoints(Molecule *molecules, Environment *enviro);
 	void assignAtomPositions(double *dev_doublesX, double *dev_doublesY, double *dev_doublesZ, Molecule *molec, Environment *enviro);
 	
};
 	
#endif