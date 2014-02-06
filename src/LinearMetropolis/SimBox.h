/*!\file
  \Class for simulation Box, including Enviroments and points to molocoles.
  \author David(Xiao Zhang) and Orlando Acevedo.
 
  This file contains functions that are used to handle enviroments and common function
  for box.
 */
#ifndef SIMBOX_H
#define SIMBOX_H 


#include "Utilities/Opls_Scan.h"
#include "Utilities/Config_Scan.h"
#include "Utilities/metroUtil.h"
#include "Utilities/Zmatrix_Scan.h"
#include "Utilities/State_Scan.h"

extern double randomFloat(const double start, const double end);

class SimBox
{
	private:
	 	Molecule changedmole;
	 	Molecule *molecules;
	 	Environment *enviro;
	 	Atom * atompool;
	 	Bond * bondpool;
	 	Angle * anglepool;
	 	Dihedral * dihedralpool;
	 	Hop *      hoppool;

	 	
	public:
	 	SimBox(Config_Scan configScan);
	 	~SimBox();
	 	Molecule *getMolecules(){return molecules;};
	 	Atom *getAtoms(){return atompool;};
	 	Environment *getEnviro(){return enviro;};
	 	int ReadStateFile(char const* StateFile);
	 	int ReadStateFile(string StateFile) { return ReadStateFile(StateFile.c_str());};
	 	int WriteStateFile(char const* StateFile); 	
	 	int WriteStateFile(string StateFile) { return WriteStateFile(StateFile.c_str());};
	 	int writePDB(char const* pdbFile);
	 	int getXFromIndex(int idx);
	 	int getYFromIndex(int x, int idx);
	 	double makePeriodic(double x, double box);
	 	double wrapBox(double x, double box);
	 	void keepMoleculeInBox(Molecule *molecule, Environment *enviro);
		//int * molecTables[];
		Table * tables;
		int molecTypenum;
	 	double getFValue(int atom1, int atom2, int **table1);//Atom *atom1, Atom *atom2, Molecule *molecules, Environment *enviro);
	 	void generatefccBox(Molecule *molecules, Environment *enviro);
	 	void generatePoints(Molecule *molecules, Environment *enviro);
	 	void assignAtomPositions(double *dev_doublesX, double *dev_doublesY, double *dev_doublesZ, Molecule *molec, Environment *enviro);
	 	int chooseMolecule();
		int changeMolecule(int molIdx);
	  	int Rollback(int moleno);
	private:
		int copyMolecule(Molecule *mole_dst, Molecule *mole_src);
		int saveChangedMole(int moleno);
};
 	
#endif