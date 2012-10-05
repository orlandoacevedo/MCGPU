/*!\file
  \Class for simulation Box, including Enviroments and points to molocoles,only save all states
  \author David(Xiao Zhang).
 
  This file contains implement of SimBox that are used to handle enviroments and common function
  for box.
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <time.h>
#include "SimBox.h"


SimBox::SimBox(char const* ConfigFile)
{
	molecules=NULL;
	enviro=NULL;
}

Molecule * SimBox::GetMolecule()
{
	return molecules;
}

Environment * SimBox::GetEnviro()
{
	return enviro;
}

int SimBox::ReadStateFile(char const* StateFile)
{
	printf("%s\n",StateFile);
	return 0;
}

int SimBox::WriteStateFile(char const* StateFile)
{
	printf("%s\n",StateFile);
	return 0;
}

void SimBox::assignAtomPositions(double *dev_doublesX, double *dev_doublesY, double *dev_doublesZ, Molecule *molec, Environment *enviro){
    //Translates each Molecule a random X,Y,and Z direction
	 //By translating every atom in that molecule by that translation

    //for each Molecule...
	 for(int i=0; i<enviro->numOfMolecules; i++){
	     for(int a=0; a<molec[i].numOfAtoms;a++){
		      Atom myAtom  =  molec[i].atoms[a];
		      myAtom.x =  dev_doublesX[i] * enviro->x + myAtom.x;
				myAtom.y =  dev_doublesY[i] * enviro->y + myAtom.y;
				myAtom.z =  dev_doublesZ[i] * enviro->z + myAtom.z;
		  }
		   keepMoleculeInBox(&molec[i],enviro);
    }
}

void SimBox::generatePoints(Molecule *molecules, Environment *enviro){

    srand((unsigned int) time(NULL));
	 //for each Molecule assign a new XYZ
    for (int i = 0; i < enviro->numOfMolecules; i++){
        double baseX = ( (double) rand() / RAND_MAX) * enviro->x;
        double baseY = ( (double) rand() / RAND_MAX) * enviro->y;
        double baseZ = ( (double) rand() / RAND_MAX) * enviro->z;
        for (int j = 0; j < molecules[i].numOfAtoms; j++){
            molecules[i].atoms[j].x += baseX;
            molecules[i].atoms[j].y += baseY;
            molecules[i].atoms[j].z += baseZ;
        }

        keepMoleculeInBox(&(molecules[i]), enviro);
    }
}

void SimBox::generatefccBox(Molecule *molecules, Environment *enviro){
	
	double cells, dcells, cellL, halfcellL;
	
	//Determine the number of unit cells in each coordinate direction
	dcells = pow(0.25 * (double) enviro->numOfMolecules, 1.0/3.0);
	cells = (int)(dcells + 0.5);
		
	//Check if numOfMolecules is a non-fcc number of molecules
	//and increase the number of cells if necessary
	while((4 * cells * cells * cells) < enviro->numOfMolecules)
			cells++;
			
	//Determine length of unit cell
	cellL = enviro->x/ (double) cells;
	halfcellL = 0.5 * cellL;
	
	//Construct the unit cell
	for (int j = 0; j < molecules[0].numOfAtoms; j++){
	molecules[0].atoms[j].x += 0.0;
    molecules[0].atoms[j].y += 0.0;
    molecules[0].atoms[j].z += 0.0;
	}
	
	for (int j = 0; j < molecules[1].numOfAtoms; j++){
	molecules[1].atoms[j].x += halfcellL;
    molecules[1].atoms[j].y += halfcellL;
    molecules[1].atoms[j].z += 0.0;
    }
    
    for (int j = 0; j < molecules[2].numOfAtoms; j++){	
	molecules[2].atoms[j].x += 0.0;
    molecules[2].atoms[j].y += halfcellL;
    molecules[2].atoms[j].z += halfcellL;
    }
    
    for (int j = 0; j < molecules[3].numOfAtoms; j++){
    molecules[3].atoms[j].x += halfcellL;
    molecules[3].atoms[j].y += 0.0;
    molecules[3].atoms[j].z += halfcellL;
    }
    
	//Init all other molecules to initial coordinates
	//Build the lattice from the unit cell by repeatedly translating
	//the four vectors of the unit cell through a distance cellL in
	//the x, y, and z directions
	for(int i = 4; i < enviro->numOfMolecules; i++){
		for (int j = 0; j < molecules[i].numOfAtoms; j++){
			molecules[i].atoms[j].x += 0.0;
    		molecules[i].atoms[j].y += 0.0;
   	 		molecules[i].atoms[j].z += 0.0;
   	 	}		
	}
	
	int offset = 0;
	for(int z = 1; z <= cells; z++)
		for(int y = 1; y <= cells; y++)
			for(int x = 1; x <= cells; x++){
				for(int a = 0; a < 4; a++){
					int i = a + offset;
					if(i < enviro->numOfMolecules){								
						for (int j = 0; j < molecules[i].numOfAtoms; j++){
							molecules[i].atoms[j].x = molecules[a].atoms[j].x + cellL * (x-1);
							molecules[i].atoms[j].y = molecules[a].atoms[j].y + cellL * (y-1);
							molecules[i].atoms[j].z = molecules[a].atoms[j].z + cellL * (z-1);
						}
					}
				}
				offset += 4;
			}
	
	//Shift center of box to the origin
	for(int i = 0; i < enviro->numOfMolecules; i++){
		for (int j = 0; j < molecules[i].numOfAtoms; j++){
			molecules[i].atoms[j].x -= halfcellL;
			molecules[i].atoms[j].y -= halfcellL;
			molecules[i].atoms[j].z -= halfcellL;
		}
	}
}
int SimBox::getXFromIndex(int idx){
    int c = -2 * idx;
    int discriminant = 1 - 4 * c;
    int qv = (-1 + sqrtf(discriminant)) / 2;
    return qv + 1;
}

int SimBox::getYFromIndex(int x, int idx){
    return idx - (x * x - x) / 2;
}

double SimBox::makePeriodic(double x, double box){
    
    while(x < -0.5 * box){
        x += box;
    }

    while(x > 0.5 * box){
        x -= box;
    }

    return x;

}

double SimBox::wrapBox(double x, double box){

    while(x > box){
        x -= box;
    }
    while(x < 0){
        x += box;
    }

    return x;
}

void SimBox::keepMoleculeInBox(Molecule *molecule, Environment *enviro){		
		for (int j = 0; j < molecule->numOfAtoms; j++){
		//X axis
			wrapBox(molecule->atoms[j].x, enviro->x);
		//Y axis
			wrapBox(molecule->atoms[j].y, enviro->y);
		//Z axis
			wrapBox(molecule->atoms[j].z, enviro->z);
		}
}

Molecule* SimBox::getMoleculeFromAtomID(Atom *a1, Molecule *molecules, Environment *enviro){
    int atomId = a1->id;
    int currentIndex = enviro->numOfMolecules - 1;
    int molecId = molecules[currentIndex].id;
    while(atomId < molecId && currentIndex > 0){
        currentIndex -= 1;
        molecId = molecules[currentIndex].id;
    }

    return &(molecules[currentIndex]);
}

double SimBox::getFValue(Atom *atom1, Atom *atom2, Molecule *molecules, Environment *enviro){
    Molecule *m1 = getMoleculeFromAtomID(atom1, molecules, enviro);
    Molecule *m2 = getMoleculeFromAtomID(atom2, molecules, enviro);
    
    if(m1->id != m2->id)
        return 1.0;
    else{
        int hops = hopGE3(atom1->id, atom2->id, m1);
        if (hops == 3)
            return 0.5;
        else if (hops > 3)
            return 1.0;
        else
            return 0.0;
    }
}

int SimBox::hopGE3(int atom1, int atom2, Molecule *molecule){
    for(int x=0; x< molecule->numOfHops; x++){
        Hop *myHop = &(molecule->hops[x]);
        //compare atoms to each hop struct in molecule
		if((myHop->atom1==atom1 && myHop->atom2==atom2) || (myHop->atom1==atom2 && myHop->atom2==atom1)){
            return myHop->hop;
        }
	}
	 return 0;
}