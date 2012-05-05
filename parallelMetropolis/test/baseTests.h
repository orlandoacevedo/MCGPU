#ifndef BASETEST_H
#define BASETEST_H

#include <math.h>
#include <sys/time.h>
#include "../../Utilities/src/metroUtil.h"

/**
  Make a dimensional position periodic
  @param x - point to be made periodic
  @param box - environment dimension
  @return - altered point
*/
double make_periodic(double x, double box);

/**
  Wrap a point around into the box
  @param x - dimension to wrap
  @param box - dimension it wrap around
  @return - altered point
*/
double wrap_into_box(double x, double box);

/**
  DEPRECATED
  Original calculate energy function
  @param coords - 2D environment array
  @param n_atoms - number of atoms
  @param box_size - box dimensions (cube)
  @param sigma - sigma value for atoms in simulation
  @param epsilon - epsilon value for atoms in simulation
  @return - total energy
*/
double calculate_energy(double **coords, int n_atoms, double *box_size, double sigma, double epsilon);

/**
  Calculates the energy assuming that all atoms are the same element.
  @param atoms - array of atoms that will be used to calculate energy
  @param enviro - simulation environment pointer
  @param molecules - pointer to molecules array (optional)
  @return - total energy
*/
double calculate_energy(Atom *atoms, Environment *enviro, Molecule *molecules=NULL);

/**
  Calculate differences between two times
  @param starttime - beginning time
  @param finishtime - ending time
  @return - difference between times
*/
long timevaldiff(struct timeval *starttime, struct timeval *finishtime);

/**
    Calculate r between two atoms
    @param a1 - first atom
    @param a2 - second atom
    @param enviro - simulation environment
    @return - r value
*/
double calc_r_value(Atom a1, Atom a2, Environment enviro);

/**
  Calculate charge energy between two atoms
  @param a1 - first atom
  @param a2 - second atom
  @param enviro - simulation environment
  @return - charge energy
*/
double calc_charge(Atom a1, Atom a2, Environment enviro);

/**
  Returns the molecule id from the atomid (on host)
  @param a1 - the atom from which to find the molecule
  @param molecules - the list of molecules to be searched
  @param enviro - simulation environment
  @param return - returns the id of the molecule
*/
int getMoleculeFromIDLinear(Atom a1, Molecule *molecules, Environment enviro);

/**
  Returns the "fudge factor" to be used in force field calculation. (on host)
  @param atom1 - the first atom in calculation
  @param atom2 - the second atom in the calculation
  @return - 1.0 if the atoms are in seperate molecules
            1.0 if the bond traversal distance is greater than 3
            0.5 if the bond traversal distance is equal to 3
            0.0 if the bond traversal is less than 3
*/
double getFValueLinear(Atom atom1, Atom atom2, Molecule *molecules, Environment *enviro);

/**
  Returns hop information between two atoms in a molecule (on host)
  @param atom1 - the id of the starting atom
  @param atom2 - the id of the ending atom
  @param molecule - the molecule that contains atom1 and atom 2
  @return - Traversal distance if hop exists, 0 otherwise
*/
int hopGE3Linear(int atom1, int atom2, Molecule molecule);

#endif //BASETEST_H
