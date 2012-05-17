/*!\file
  \brief Functions that execute on the device or interface with the device.
  \author Alexander Luchs, Riley Spahn, Seth Wooten
 
 This file contains functions that are either executed on the device or 
 interface with the device.
 */

#ifndef METROCUDAUTIL_CUH
#define METROCUDAUTIL_CUH

#define DEBUG

#include <float.h>
#include <cstdlib>
#include <cuda.h>
#include <curand_kernel.h>
#include <math.h>
#include "../../Utilities/src/metroUtil.h"
#include "../test/baseTests.h"
#include <curand.h>
#define THREADS_PER_BLOCK 128
#define PI 3.14159265
#define cudaErrorCheck(call) { cudaAssert(call,__FILE__,__LINE__); }

/*!
Representation of a molecule on the device. It is difficult to use the original
molecule structure because it contains a pointer to an array.  On the device we have
glolbal arrays of:
DeviceAtoms, Atoms, Bonds, Angles, Dihedrals and Hops
Each of these arrays contains all of the atoms, dihedrals, etc. in the simulation.
The DeviceMolecule has a field that is the first index in the global array and a 
number of that element in the molecule.  This is used to index the correct elements
in the global array.
*/
struct DeviceMolecule{
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

/**
    Returns error code for cuda calls
    @param err - error code
    @param file - source file name
    @param line - line in source
*/
void cudaAssert(const cudaError err, const char *file, const int line);

/**
  Creates an instance of the device molecule structure.
  @param id - the molecules id
  @param atomStart - the first index of this molecules atom in the atom array.
  @param numOfAtoms - the number of atoms in this molecule.
  @param bondStart - this molecule's first index in the bond array
  @param numOfBonds - the number of bonds in this atom.
  @param angleStart - this molecule's first angle in the array of angles
  @param numOfAngles - the number of angles in this molecule
  @param dihedralStart - this molecule's first dihedral in the dihedral array.
  @param numOfDihedrals - the number of dihedrals in this molecule
  @param hopStart - this molecule's first hop in the hop array.
  @param numOfHops - the number of hops in this molecule.
  @return - instance of DeviceMolecule that is created
*/
DeviceMolecule createDeviceMolecule(int id, int atomStart, int numOfAtoms,
        int bondStart, int numOfBonds, int angleStart, int numOfAngles,
        int dihedralStart, int numOfDihedrals, int hopStart, int numOfHops); 

/**
  @param idx - the index in the 1 dimensional array of energies
  @return - the id of the X atom (atom in the calculation with the larger id)
*/
__device__ int getXFromIndex(int idx);

/**
  @param x - the id of the x atom
  @param idx - the index in the 1 dimensional array of energies
  @return - the id of the Y atom (atom in the calculation with the smaller id)
*/
__device__ int getYFromIndex(int x, int idx);

/**
  @param x - the variable to make periodic
  @param box - the size of the period
  @return - x after being made periodic
*/
__device__ double makePeriodic(double x, const double box);

/**
  @param x - the value to continue on the other side of the box
  @param box - the length of one side of the box (cube)
  @return - x after being wrapped around the box dimension
*/
double wrapBox(double x, double box);

/**
    Keeps a molecule intact within the box
    @param molecule - molecule to keep in box
    @param enviro - environment structure to keep bounds
*/
void keepMoleculeInBox(Molecule *molecule, Environment *enviro);

/**
  Calculates the Lennard-Jones energy between 2 atoms
  @param atom1 - the first atom in the pair
  @param atom2 - the second atom in the pair
  @param enviro - the environmental variables
*/
__device__ double calc_lj(Atom atom1, Atom atom2, Environment enviro); 

/**
    Global function for GPU to assign random doubles to atom positions.
    @param dev_doublesX - array of randomly generated doubles 0.0 to 1.0 to be
    used to translate in the x direction.
    @param dev_doublesY - array of randomly generated doubles 0.0 to 1.0 to be
    used to translate in the y direction.
    @param dev_doublesZ - array of randomly generated doubles 0.0 to 1.0 to be
    used to translate in the z direction.
    @param atoms - array of atoms to be positioned
    @param enviro - Environment of the simulation
*/
__global__ void assignAtomPositions(double *dev_doublesX, double *dev_doublesY, double *dev_doublesZ, Atom *atoms, Environment *enviro);

/**
  Generate random positions for atoms in the box
  nVidia CURAND reference: http://developer.download.nvidia.com/compute/cuda/5_0/toolkit/docs/CURAND_Library.pdf
  @param atoms - array of atoms to generate positions
  @param enviro - enviroment structure defining the box
*/
void generatePoints(Atom *atoms, Environment *enviro);

/**
  Generate random positions for atoms in all molecules in the box
  Does this on the CPU side unlike the atoms version.
  @param molecules - array of molecules to generate positions
  @param enviro - enviroment structure defining the box
*/
void generatePoints(Molecule *molecules, Environment *enviro);

/**
  This is a wrapper function for the calcEnergy kernel.
  @param *atoms - the array of atoms.
  @param enviro - the environmental variables.
  @param molecules - array of molecules to be used if calculating energy from molecules and not an array of atoms.
  @return - the total energy of the system.
*/
double calcEnergyWrapper(Atom *atoms, Environment *enviro, Molecule *molecules=NULL);

/**
  Wrapper function for the calcEnergy kernel.
  @param molecules - array of molecules in the system
  @param enviro - the environment of the system.
*/
double calcEnergyWrapper(Molecule *molecules, Environment *enviro);

/**
  This calculates nonbonded energy between two atoms on host for fringe
  disagreeable atoms.
  @param atom1 - first atom
  @param atom2 - second atom
  @param enviro - the environment for the system
  @param molecules - array of molecules
*/
double calcEnergyOnHost(Atom atom1, Atom atom2, Environment *enviro, Molecule *molecules=NULL);


/**
  Calculates the energy between n atoms where n is the
  the number of threads in the block. The block's sum is then stored 
  in the energySum array at its block index position.
  @param *atoms - the array of atoms
  @param enviro - the environmental variables
  @param *energySum - the array of block sums
  @param *dev_molecules - optional array of device representation of molecules
  @param *hops - optional array of hops for fudge factor
*/
__global__ void calcEnergy(Atom *atoms, Environment *enviro, double *energySum, DeviceMolecule *dev_molecules=NULL, Hop *hops=NULL);

/**
  Calculates the charge portion of the force field energy calculation between two atoms.
  @param atom1 - the first atom in the calculation
  @param atom2 - the second atom in the calculation
  @param enviro - simulation environment pointer
  @return - the charge portion of the force field.
*/
__device__ double calcCharge(Atom atom1, Atom atom2, Environment *enviro);

/**
  Returns the molecule id from the atomid
  @param a1 - the atom from which to find the molecule
  @param dev_molecules - the list of DeviceMolecules to be searchedi
  @param enviro - simulation environment
  @param return - returns the id of the molecule
*/
__device__ int getMoleculeFromAtomID(Atom a1, DeviceMolecule *dev_molecules, Environment enviro);

/**
  Returns the "fudge factor" to be used in force field calculation.
  @param atom1 - the first atom in calculation
  @param atom2 - the second atom in the calculation
  @param dev_molecules - array of all DeviceMolecules in simulation
  @param hops - array of all hops in simulation
  @return - 1.0 if the atoms are in seperate molecules
            1.0 if the bond traversal distance is greater than 3
            .5 if the bond traversal distance is equal to 3
            0.0 if the bond traversal is less than 3
*/
__device__ double getFValue(Atom atom1, Atom atom2, DeviceMolecule *dev_molecules, Environment *enviro, Hop *hops);

/**
  Return if the two atom ids are have a hop value >=3
  returns hop distance if true and 0 if false
  @param atom1 - the id of the starting atom
  @param atom2 - the id of the ending atom
  @param dev_molecule - the DeviceMolecule that contains atom1 and atom 2
  @param molecule_hops - array of hops in dev_molecule
  @return - 0 if hop is not found
            Traversal distance if found
*/
__device__ int hopGE3(int atom1, int atom2, DeviceMolecule dev_molecule, Hop *molecule_hops);

/**
  Returns the molecule id from the atomid (on host)
  @param atom - the atom from which to find the molecule
  @param molecules - the list of molecules to be searched
  @param enviro - simulation environment
  @param return - returns the id of the molecule
*/
Molecule* getMoleculeFromAtomIDHost(Atom a1, Molecule *molecules, Environment enviro);

/**
  Returns the "fudge factor" to be used in force field calculation. (on host)
  @param atom1 - the first atom in calculation
  @param atom2 - the second atom in the calculation
  @param molecules - array of molecules in the simulation
  @param enviro - simulation environment pointer
  @return - 1.0 if the atoms are in seperate molecules
            1.0 if the bond traversal distance is greater than 3
            0.5 if the bond traversal distance is equal to 3
            0.0 if the bond traversal is less than 3
*/
double getFValueHost(Atom atom1, Atom atom2, Molecule *molecules, Environment *enviro);

/**
  Return if the two atom ids are have a hop value >=3
  returns hop distance if true and 0 if false (on host)
  @param atom1 - the id of the starting atom
  @param atom2 - the id of the ending atom
  @param molecule - the molecule that contains atom1 and atom 2
  @return - 0 if hop is not found
            Traversal distance if found
*/
int hopGE3Host(int atom1, int atom2, Molecule molecule);

/**
  Returns sqrt(d1 * d2)
  @param d1 - the first double in the calculation
  @param d2 - the second double in the calculation
  @return - sqrt(d1 * d2)
*/
__device__ double calcBlending(double d1, double d2);

/**
  Rotates a molecule about a given atom a random amount
  @param molecule - the molecule to be rotated
  @param pivotAtom - the atom that the molecule is rotated about
  @param maxRotation - the maximum number of degrees for each axis of rotation
*/
void rotateMolecule(Molecule molecule, Atom pivotAtom, double maxRotation);

/**

  Deep copies an array of molecules from the host to the device.
  Assumes that the memory on the device has been allocated.
  @param molec_d - pointer on the device to contain an array of molecules.
  @param molec_h - pointer to the host array of molecules to be copied to the device.
  @param numOfMolecules - the number of molecules in the array

  ====All arrays below are allocated on the device.===
  @param atoms_d - array of atoms as long as the total number of atoms.
  @param bonds_d - array of bonds as long as the total number of bonds.
  @param angles_d - array of angles as long as the total number of angles.
  @param dihedrals_d - array of dihedrals as long as the total number of dihedrals
  @param hops_d - array of hops as long as the total number of hops.
*/
void  moleculeDeepCopyToDevice(DeviceMolecule *molec_d, Molecule *molec_h,
       int numOfMolecules, Atom *atoms_d, Bond *bonds_d, Angle *angles_d,
        Dihedral *dihedrals_d, Hop *hops_d);

/**Deep copies an array of molecules from the device to the host.
  Assumes that the memory on the host has already been correctly allocated.
  @param molec_h - pointer to the host array of molecules to be copied from the device
  @param molec_d - pointer to the host array of molecules to be copied from the device to the host.
  @param numOfMolecules - the number of molecules in the array
  @param atoms_d - pointer to the device array of atoms
  @param bonds_d - pointer to the device array of bonds
  @param angles_d - pointer to the device array of angles
  @param dihedrals_d - pointer to the device array of dihedrals
  @param hops_d - pointer to the device array of hops
*/
void moleculeDeepCopyToHost(Molecule *molec_h, DeviceMolecule *molec_d,
        int numOfMolecules,Atom *atoms_d, Bond *bonds_d, Angle *angles_d,
        Dihedral *dihedrals_d, Hop *hops_d);

#ifdef DEBUG

/**
   @param atoms - array of atoms to be used in the atom1 parameter
   @param molecules - array of molecules to be used in the molecule parameter
   @param enviros - the array of environments to be used in the enviro parameter
   @param numberOfTests - the number of tests to be run. All other arrays must 
   be of length numberOfTests
   @param answers - array where results of tests are stored
 */
__global__ void testGetMoleculeFromID(Atom *atoms, DeviceMolecule *molecules,
        Environment enviros, int numberOfTests, int *answers);

/**
  Kernel function that will be used to test test calcBlending() function.
*/
__global__ void testCalcBlending(double *d1, double *d2, double *answers, int numberOfTests);

/**
  Kernel function that will be used to test the getFValue() function.
*/
__global__ void testGetFValue(Atom *atom1List, Atom *atom2List, DeviceMolecule *molecules, Environment *enviro, double *fValues, int numberOfTests, Hop *dev_hops);

/**
  Kernel function that will be used to test the calcCharge() function.
*/
__global__ void testCalcCharge(Atom *atoms1, Atom *atoms2, double *charges, Environment *enviro);

/**
    Kernel call that will be used to test the getXFromIndexFunction
*/
__global__ void testGetXKernel(int *xValues, int totalTests);

/**
  Kernel function that will be used to test getYFromIndex device function
*/
__global__ void testGetYKernel(int *xValues, int *yValues, int numberOfTests);

/**
  Kernel function used to facilitate tests of the makePeriodic device function
*/
__global__ void testMakePeriodicKernel(double *x, double *box, int numberOfTests);

/**
  Kernel function used to facilitate tess of the wrap box device function
*/
__global__ void testWrapBoxKernel(double *x, double *box, int numberOfTests);

/**
  Kernel to test the calc_lj device function
*/
__global__ void testCalcLJ(Atom *atoms, Environment *enviro, double *energy);

/**
  Kernel function to test the calculation of bond distance between two atoms
*/
__global__ void testGetDistance(Atom *atom1List, Atom *atom2List, Molecule molecule, Environment *enviro, int *distances, int numberOfTests);

#endif //DEBUG

#endif //METROCUDAUTIL_CUH
