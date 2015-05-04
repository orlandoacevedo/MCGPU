/*
	Neighbor linded-list class for use in NLC energy calculations.

	Author: Jared Brown
	Created: April 17, 2014
*/

#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

#include <string>
#include "Metropolis/Box.h"
#include "Metropolis/DataTypes.h"
#include "Metropolis/SimulationArgs.h"
#include "Metropolis/Utilities/StructLibrary.h"

#define NMAX 100000  /* Maximum number of atoms which can be simulated */
#define NCLMAX 10000 /* Maximum number of linked-list cells */
#define EMPTY -1

class NeighborList
{
	public:
		NeighborList(Molecule *molecules, Environment *enviro);
		~NeighborList();
		
		int numCells[3];            	/* Number of cells in the x|y|z direction */
		int numCellsYZ;					/* Total number of cells in YZ plane */
		int numCellsXYZ;				/* Total number of cells in XYZ area*/
		Real lengthCell[3];         	/* Length of a cell in the x|y|z direction */
		int head[NCLMAX];    			/* Headers for the linked cell lists */
		int linkedCellList[NMAX];       /* Linked cell lists */
		
		Real region[3];
		Real rrCut;
};

#endif
