#ifndef SIMBOX_CONSTANTS_H
#define SIMBOX_CONSTANTS_H

// Represents the X - axis.
#define X_COORD 0

// Represents the Y - axis.
#define Y_COORD 1

// Represents Z - axis.
#define Z_COORD 2

// Holds the number of spatial dimenions in the simulation.
#define NUM_DIMENSIONS 3

// MOLECULE DATA CONSTANTS

// Indicates the row of moleculeData that holds the start index of each
//     molecule in atomCoordinates and atomData.
#define MOL_START 0

// Indicates the row of moleculeData that holds the number of atoms of each
//     molecule.
#define MOL_LEN  1

// Indicates the row of moleculeData that holds the start index of each
//     molecule's primary index(es) in primaryIndexes.
#define MOL_PIDX_START 2

// Indicates the row of moleculeData that holds the number of primary indexes
//     that each molecule has.
#define MOL_PIDX_COUNT 3

// Indicates the row of moleculeData that hold the type of each molecule.
#define MOL_TYPE 4


// Indicates the row of moleculeData that holds the start index of each
//     molecule's bonds in bondData and bondLengths.
#define MOL_BOND_START 5

// Indicates the row of moleculeData that holds the number of bonds in each
//     molecule.
#define MOL_BOND_COUNT 6

// Indicates the row of moleculeData that holds the start index of each
//     molecule's angles in angleData and angleSizes.
#define MOL_ANGLE_START 7

// Indicates the row of moleculeData that holds the number of angles in each
//     molecules.
#define MOL_ANGLE_COUNT 8

// Indicates the number of rows of moleculeData.
#define MOL_DATA_SIZE 9

// ATOM DATA CONSTANTS

// Indicates the row of atomData that holds the value of sigma for each atom.
#define ATOM_SIGMA 0

// Indicates the row of atomData that holds the value of epsilon for each atom.
#define ATOM_EPSILON 1

// Indicates the row of atomData that holds the charge of each atom.
#define ATOM_CHARGE 2

// Indicates the number of rows of atomData.
#define ATOM_DATA_SIZE 3

// BOND DATA CONSTANTS

// Indicates the row of bondData that holds the 1st atom index for each bond.
#define BOND_A1_IDX 0

// Indicates the row of bondData that holds the 2nd atom index for each bond.
#define BOND_A2_IDX 1

// Indicates the row of bondData that holds the force constant for each bond.
#define BOND_KBOND 2

// Indicates the row of bondData holding the equilibrium distance of each bond.
#define BOND_EQDIST 3

// Indicates the row of bondData that records whether each bond is variable.
#define BOND_VARIABLE 4

// Indicates the number of rows in bondData.
#define BOND_DATA_SIZE 5

// ANGLE DATA CONSTANTS

// Indicates the row of angleData that holds the first atom's index for each
//     angle.
#define ANGLE_A1_IDX 0

// Indicates the row of angleData that holds the middle atom's index for each
//     angle
#define ANGLE_MID_IDX 1

// Indicates the row of angleData that holds the third atom's index for each
//     angle.
#define ANGLE_A2_IDX 2

// Indicates the row of angleData that holds the force constant of each angle.
#define ANGLE_KANGLE 3

// Indicates the row of angleData holding the equilibrium size of each angle.
#define ANGLE_EQANGLE 4

// Indicates the row of angleData holding whether each angle is variable.
#define ANGLE_VARIABLE 5

// Indicates the number of rows in angleData.
#define ANGLE_DATA_SIZE 6

#endif
