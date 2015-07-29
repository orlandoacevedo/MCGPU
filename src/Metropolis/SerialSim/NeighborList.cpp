/*
	Neighbor linded-list class for use in NLC energy calculations.

	Author: Jared Brown
	Created: April 17, 2014
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>

#include <math.h>
//#include <string>
#include "Metropolis/DataTypes.h"
#include "Metropolis/SimulationArgs.h"
#include "Metropolis/Utilities/FileUtilities.h"
#include "NeighborList.h"

//Constructor & Destructor
NeighborList::NeighborList(Molecule *molecules, Environment *enviro)
{
	rrCut = enviro->cutoff * enviro->cutoff;
	
	// Compute the # of cells for linked cell lists
	for (int k = 0; k < 3; k++)
	{
		if (k == 0)
			region[k] = enviro->x;
		else if (k == 1)
			region[k] = enviro->y;
		else if (k == 2)
			region[k] = enviro->z;
		
		numCells[k] = region[k] / enviro->cutoff; 
		lengthCell[k] = region[k] / numCells[k];
	}
		
  	/* Make a linked-cell list --------------------------------------------*/
	numCellsYZ = numCells[1] * numCells[2];
	numCellsXYZ = numCells[0] * numCellsYZ;
		
	// Reset the headers, head
	for (int c = 0; c < numCellsXYZ; c++) 
	{
		head[c] = EMPTY;
	}

	// Scan cutoff index atom in each molecule to construct headers, head, & linked lists, linkedCellList
	for (int i = 0; i < enviro->numOfMolecules; i++)
	{
		std::vector<int> molPrimaryIndexArray = (*(*(enviro->primaryAtomIndexArray))[molecules[i].type]);
		int primaryIndex = molPrimaryIndexArray[0]; // Use first primary index to determine cell placement
		
		int vectorCells[3];			
		vectorCells[0] = molecules[i].atoms[primaryIndex].x / lengthCell[0]; 
		vectorCells[1] = molecules[i].atoms[primaryIndex].y / lengthCell[1];
		vectorCells[2] = molecules[i].atoms[primaryIndex].z / lengthCell[2];
		
		// Translate the vector cell index to a scalar cell index
		int c = vectorCells[0] * numCellsYZ 
			+ vectorCells[1] * numCells[2] 
			+ vectorCells[2];

		// Link to the previous occupant (or EMPTY if you're the 1st)
		linkedCellList[i] = head[c];

		// The last one goes to the header
		head[c] = i;
	} /* Endfor molecule i */
}

NeighborList::~NeighborList() {
	free(numCells);
	free(lengthCell);
	free(head);
	free(linkedCellList);
	free(region);
}
	
