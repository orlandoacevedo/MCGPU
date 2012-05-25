/*!\file
  \brief Functions used for benchmarking against a linear implementation.
  \author Alexander Luchs, Riley Spahn, Seth Wooten, and Orlando Acevedo
 
  This file contains functions that are used as a linear benchmark to compare
  to our parallel implmentation.
 */
#ifndef METROLINEARUTIL_H
#define METROLINEARUTIL_H

#include <float.h>
#include <cstdlib>
#include <math.h>
#include "../../Utilities/src/metroUtil.h"
#include <curand.h>
#define THREADS_PER_BLOCK 128
#define PI 3.14159265
// Linked-cell neighbor list constants
#define NMAX 100000  /* Maximum number of atoms which can be simulated */
#define NCLMAX 10000 /* Maximum number of linked-list cells */
#define EMPTY -1

/**
  Gets the X atom in an energy calculation based on the thread index.  
  Implementation is based on the quadratic equation.  See architectural spike
  cycle binder for further information.
  @param idx - the index in the 1 dimensional array of energies
  @return - the id of the X atom
*/
int getXFromIndex(int idx);

/**
  Returns the y atom to be used in the energy calculation.
  See arch. spike cycle binder for full explanation.
  @param x - the id of the x atom.
  @param idx - the index in the 1 dimensional array of energies.
  @return - the id of the y atom to be used in the energy calculation.
*/
int getYFromIndex(int x, int idx);

/**
  Make periodic function from the original example.
  @param x - the variable to make periodic
  @param box - the size of the period
  @return - x after being made periodic
*/
double makePeriodic(double x, const double box);

/**
  If an atom moves out of the "box" then it will be wrapped to the other
  side of the box like a torus.
  @param x - the value to continue on the other side of the box
  @param box - the length of one side of the box (cube)
  @return - the position after wrapping around the box.
*/
double wrapBox(double x, double box);

/**
    Keeps a molecule intact within the box
    @param molecule - molecule to keep in box
    @param enviro - simulation environment pointer
*/
void keepMoleculeInBox(Molecule *molecule, Environment *enviro);

/**
  Calculates the energy between 2 atoms
  @param atom1 - the first atom in the pair
  @param atom2 - the second atom in the pair
  @param enviro - simulation environment pointer
  @param the energy between the two atoms.
*/
double calc_lj(Atom atom1, Atom atom2, Environment enviro); 

/**
    Global function for GPU to assign random doubles to atom positions.
    @param dev_doublesX - array of randomly generated doubles 0.0 to 1.0 (x dimension)
    @param dev_doublesY - array of randomly generated doubles 0.0 to 1.0 (y dimension)
    @param dev_doublesZ - array of randomly generated doubles 0.0 to 1.0 (z dimension)
    @param atoms - array of atoms to be positioned
    @param enviro - Environment of the simulation
*/
void assignAtomPositions(double *dev_doublesX, double *dev_doublesY, double *dev_doublesZ, Molecule *molec, Environment *enviro);

/**
  Generate random positions for atoms in the box
  @param molec - array of molecules to generate positions
  @param enviro - enviroment structure defining the box
*/
void generatePoints(Molecule *molec, Environment *enviro);

/**
  Generate positions for molecules in the box based on fcc-lattice.
  Ideal # = 4*molecules^3, i.e., 4, 32, 108, 256.
  @param molec - array of molecules to generate positions
  @param enviro - enviroment structure defining the box
*/
void generatefccBox(Molecule *molecules, Environment *enviro);

/**
  This is a wrapper function for the calcEnergy kernel.
  @param atoms - the array of atoms to calculate energies between
  @param enviro - simulation environment pointer
  @param molecules - optional array of molecules within which the atoms reside
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
  Calculates the energy between n atoms where n is the
  the number of threads in the block. The block's sum is then stored 
  in the energySum array at its block index position.
  @param atoms - the array of atoms
  @param enviro - simulation environment pointer
  @param energySum - the array of block sums
*/
void calcEnergy(Atom *atoms, Environment *enviro, double *energySum);

/**
  Calculates the charge portion of the force field energy calculation between two atoms.
  @param atom1 - the first atom in the calculation
  @param atom2 - the second atom in the calculation
  @param enviro - simulation environment pointer
  @return - the charge portion of the force field.
*/
double calcCharge(Atom atom1, Atom atom2, Environment *enviro);

/**
  Calculates the nonbonded energy portion of the force field energy calculation between two atoms.
  @param atom1 - the first atom in the calculation
  @param atom2 - the second atom in the calculation
  @param enviro - simulation environment pointer
  @return - the nonbonded energy portion of the force field.
*/
double calcNonBondEnergy(Atom atom1, Atom atom2, Environment *enviro);

/**
  Wrapper function for the calcEnergy_NLC kernel.
  @param molecules - array of molecules in the system
  @param enviro - the environment of the system.
*/
double calcEnergyWrapper_NLC(Molecule *molecules, Environment *enviro);

/**
  Calculates the energy between all atoms using a 
  linked-cell neighbor list.
  @param atoms - the array of atoms
  @param enviro - simulation environment pointer
  @param molecules - array of molecules in the system
*/
double calcEnergy_NLC(Atom *atoms, Environment *enviro, Molecule *molecules);

/**
  Returns the "fudge factor" to be used in force field calculation.
  @param atom1 - the first atom in calculation
  @param atom2 - the second atom in the calculation
  @param molecules - array of molecules in the simulation
  @param enviro - simulation environment pointer
  @return - 1.0 if the atoms are in seperate molecules
            1.0 if the bond traversal distance is greater than 3
            0.5 if the bond traversal distance is equal to 3
            0.0 if the bond traversal is less than 3
*/
double getFValue(Atom *atom1, Atom *atom2, Molecule *molecules, Environment *enviro);

/**
  Return if the two atom ids are have a hop value >=3
  returns 1 if true and 0 if false
  @param atom1 - the id of the starting atom
  @param atom2 - the id of the ending atom
  @param molecule - the molecule that contains atom1 and atom 2
  @return - 1 if true and 0 otherwise.
*/
int hopGE3(int atom1, int atom2, Molecule *molecule);

/**
  Returns sqrt(d1 * d2)
  @param d1 - the first double in the calculation
  @param d2 - the second double in the calculation
  @return - sqrt(d1 * d2)
*/
double calcBlending(double d1, double d2);

/**
  Returns the molecule id from the atomid
  @param atom - the atom from which to find the molecule
  @param molecules - the list of molecules to be searched
  @return - returns a pointer to the molecule
*/
Molecule* getMoleculeFromAtomID(Atom *a1, Molecule *molecules, Environment *enviro);

/**
  Rotates a molecule about a given atom a random amount
  @param molecule - the molecule to be rotated
  @param pivotAtom - the atom that the molecule is rotated about
  @param maxRotation - the maximum number of degrees for each axis of rotation
*/
void rotateMolecule(Molecule molecule, Atom pivotAtom, double maxRotation);

#endif
