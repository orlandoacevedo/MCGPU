/*!\file
  \Class for simulation Box, including Enviroments and points to molocoles.
  \author David(Xiao Zhang).
 
  This file contains functions that are used to handle enviroments and common function
  for box.
 */
#ifndef GPUSIMBOX_H
#define GPUSIMBOX_H 


#include "Utilities/Opls_Scan.h"
#include "Utilities/Config_Scan.h"
#include "Utilities/metroUtil.h"
#include "Utilities/Zmatrix_Scan.h"
#include "Utilities/State_Scan.h"
#include "LinearMetropolis/SimBox.h"

/*!
Representation of a molecule on the device. It is difficult to use the original
molecule structure because it contains a pointer to an array.  On the device we have
global arrays of:
DeviceAtoms, Atoms, Bonds, Angles, Dihedrals and Hops
Each of these arrays contains all of the atoms, dihedrals, etc. in the simulation.
The DeviceMolecule has a field that is the first index in the global array and a 
number of that element in the molecule.  This is used to index the correct elements
in the global array.
*/
struct DeviceMolecule
{
    /*!Number to uniquely identify this molecule.*/
    int id;

    /*!This device molecule's first index in the global array of atoms.*/
    int atomStart; 
    /*!The number of atoms contained in this moleclue.*/
    int numOfAtoms;

    /*!The DeviceMolecule's first index in the global array of bonds.*/
    int bondStart;
    /*!The number of bonds in the DeviceMolecule.*/
    int numOfBonds;
    
    /*!The DeviceMolecule's first index in the global array of angles.*/
    int angleStart;
    /*!The number of angle's in this DeviceMolecule.*/
    int numOfAngles;

    /*!The DeviceMolecule's first index in the global array of dihedrals.*/
    int dihedralStart;
    /*!The number of dihedrals in the DeviceMolecule.*/
    int numOfDihedrals;

    /*!The DeviceMolecule's first index in the global array of hops.*/
    int hopStart;
    /*!The number of hops in the DeviceMolecule*/
    int numOfHops;
};

void cudaAssert(const cudaError err, const char *file, const int line);

#define cudaErrorCheck(call) { cudaAssert(call,__FILE__,__LINE__); }
#define cudaFREE(ptr) if(ptr!=NULL) { cudaFree(ptr);ptr=NULL;}

class GPUSimBox
{
    private:
       SimBox * innerbox;
       DeviceMolecule *molec_d;
       Bond *bonds_d;
       Angle *angles_d;
       Dihedral *dihedrals_d;
       Hop *hops_d;
       Atom *atoms_device;
       Environment *enviro_device; 	
       
       size_t atomSize;
       size_t dMolecSize;	
       size_t bondSize;
       size_t angleSize;
       size_t dihedralSize;
       size_t hopSize;

    public:
     	GPUSimBox(Config_Scan configScan);
     	~GPUSimBox();

    SimBox *getSimBox() { return innerbox;};

    DeviceMolecule *getdevDeviceMolecule() { return molec_d;};
    Bond *getdevBond() { return  bonds_d;};
    Angle *getdevAngle() {return angles_d;};
    Dihedral * getdevDihedral() {return dihedrals_d;};
    Hop * getdevHop() {return hops_d;};
    Atom *getdevAtom() {return atoms_device;};
    Environment *getdevEnvironment() {return enviro_device;};  
   
    Atom *gethostAtom() {return innerbox->getAtoms();};
    Environment *gethostEnvironment() {return innerbox->getEnviro();};  
    Molecule *gethostMolecules() {return innerbox->getMolecules();};  
  
 	int initGPUSimBox(SimBox *hostbox);
 	int CopyBoxtoHost(SimBox *hostbox);
 	int CopyBoxtoDevice(SimBox *hostbox);

 	int getXFromIndex(int idx);
 	int getYFromIndex(int x, int idx);
 	double makePeriodic(double x, double box);
 	double wrapBox(double x, double box);
 	void keepMoleculeInBox(Molecule *molecule, Environment *enviro);
 	int hopGE3Host(int atom1, int atom2, Molecule molecule);
 	double getFValueHost(Atom atom1, Atom atom2, Molecule *molecules, Environment *enviro);
 	void generatefccBox(Molecule *molecules, Environment *enviro);
 	void generatePoints(Molecule *molecules, Environment *enviro);
 	void assignAtomPositions(double *dev_doublesX, double *dev_doublesY, double *dev_doublesZ, Molecule *molec, Environment *enviro);
  int ChangeMolecule();
 	int Rollback(int moleno);
};
 	
#endif