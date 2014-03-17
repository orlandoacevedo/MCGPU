/*
	New version of SimBox

	Author: Nathan Coleman
	Last Changed: February 21, 2014
	
	-> February 26, 27, by Albert Wallace
*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "SerialBox.h"

using namespace std;

double randomFloat(const double start, const double end){return 0.0;}


SerialBox::SerialBox(IOUtilities configScan) : Box()
{
    configScan.pullInDataToConstructSimBox();
	// atoms = NULL;
	// molecules = NULL;
 //    stringstream ss;
	// environment = new Environment();
	// memcpy(environment, configScan.currentEnvironment, sizeof(Environment));

	// configScan.scanInOpls();
	// configScan.scanInZmatrix();

	// int moleculeIndex = 0;
	// int atomCount = 0;

	// vector<Molecule> molecVec = configScan.buildMolecule(atomCount);
	// int molecMod = environment->numOfMolecules % molecVec.size();

	// if (molecMod != 0)
 //   {
 //       environment->numOfMolecules += molecVec.size() - molecMod;
 //       std::cout << "Number of molecules not divisible by specified z-matrix. Changing number of molecules to: " << environment->numOfMolecules << std::endl;
 //    }

 //    molecules = (Molecule *)malloc(sizeof(Molecule) * environment->numOfMolecules);
    
 //    int molecDiv = environment->numOfMolecules / molecVec.size();
 //    molecTypenum=molecVec.size();
    
 //    int count[5];//sum up number of atoms,bonds,angles,dihedrals,hops

	// tables = new Table[molecVec.size()];
 //    memset(count,0,sizeof(count));
	// int currentAtomCount = 0;
	
 //    for(int j = 0; j < molecVec.size(); j++)
 //    {
 //  		Molecule molec1 = molecVec[j];   
 //         //Copy data from vector to molecule
 //  		count[0]+=molec1.numOfAtoms;
 //  		count[1]+=molec1.numOfBonds;
 //  		count[2]+=molec1.numOfAngles;
 //  		count[3]+=molec1.numOfDihedrals;
 //  		count[4]+=molec1.numOfHops;
  		
 //  		std::cout << "before table building. Number of atom "<< molec1.numOfAtoms<<std::endl;
  		
 //  		Hop *myHop = molec1.hops;
 //  		int **table;
 //  		table = new int*[molec1.numOfAtoms];
 //  		for(int k = 0; k< molec1.numOfAtoms;k++)
 //  			table[k] = new int[molec1.numOfAtoms];
 //  		//int table[molec1.numOfAtoms][molec1.numOfAtoms];
 //  		for(int test = 0; test< molec1.numOfAtoms;test++)
 //        {
 //  			for(int test1 = 0; test1 < molec1.numOfAtoms; test1++)
 //            {
 //  				table[test][test1] = 0;
 //  			}
	// 	}
		
	// 	for(int k2 = 0; k2<molec1.numOfHops;k2++)
 //        {
	// 		int atom1 = myHop->atom1;
	// 		std::cout << "atom1: "<< atom1-currentAtomCount <<std::endl;
	// 		int atom2 = myHop->atom2;
	// 		std::cout << "atom2: "<< atom2-currentAtomCount<<std::endl;
	// 		std::cout << "hop: " << myHop->hop <<std::endl;
	// 		table[atom1-currentAtomCount][atom2-currentAtomCount] = myHop->hop;
	// 		table[atom2-currentAtomCount][atom1-currentAtomCount] = myHop->hop;
	// 		myHop++;
	// 	}
	  
	//    std::cout << "after table building"<<std::endl;
	//    tables[j] = Table(table); //createTable is in metroUtil
	//    currentAtomCount += molec1.numOfAtoms;
	//    std::cout << "after table creation. Current atom cout: "<< currentAtomCount<<std::endl;
 //    }
     
 //    atoms     =(Atom *)malloc(sizeof(Atom)*molecDiv*count[0]);
 //    bonds     =(Bond *)malloc(sizeof(Bond)*molecDiv*count[1]);
 //    angles    =(Angle *)malloc(sizeof(Angle)*molecDiv*count[2]);
 //    dihedrals =(Dihedral *)malloc(sizeof(Dihedral)*molecDiv*count[3]);
 //    hops      =(Hop *)malloc(sizeof(Hop)*molecDiv*count[4]);  
 //    memset(atoms,0,sizeof(Atom)*molecDiv*count[0]);
 //    memset(bonds,0,sizeof(Bond)*molecDiv*count[1]);
 //    memset(angles,0,sizeof(Angle)*molecDiv*count[2]);
 //    memset(dihedrals,0,sizeof(Dihedral)*molecDiv*count[3]);
 //    memset(hops,0,sizeof(Hop)*molecDiv*count[4]);

 //    //arrange first part of molecules
 //    memset(count,0,sizeof(count));
 // 	for(int j = 0; j < molecVec.size(); j++)
 //    {
 // 	      //Copy data from vector to molecule
 //        Molecule molec1 = molecVec[j];   

 //        molecules[j].atoms = (Atom *)(atoms+count[0]);
 //        molecules[j].bonds = (Bond *)(bonds+count[1]);
 //        molecules[j].angles = (Angle *)(angles+count[2]);
 //        molecules[j].dihedrals = (Dihedral *)(dihedrals+count[3]);
 //        molecules[j].hops = (Hop *)(hops+count[4]);

 //        molecules[j].moleculeIdentificationNumber = molec1.moleculeIdentificationNumber;
 //        molecules[j].numOfAtoms = molec1.numOfAtoms;
 //        molecules[j].numOfBonds = molec1.numOfBonds;
 //        molecules[j].numOfDihedrals = molec1.numOfDihedrals;
 //        molecules[j].numOfAngles = molec1.numOfAngles;
 //        molecules[j].numOfHops = molec1.numOfHops;

 //        count[0]+=molec1.numOfAtoms;
 //        count[1]+=molec1.numOfBonds;
 //        count[2]+=molec1.numOfAngles;
 //        count[3]+=molec1.numOfDihedrals;
 //        count[4]+=molec1.numOfHops;

 //        //get the atoms from the vector molecule
 //        for(int k = 0; k < molec1.numOfAtoms; k++)
 //        {
 //            molecules[j].atoms[k] = molec1.atoms[k];
 //        }               
           
 //        //assign bonds
 //        for(int k = 0; k < molec1.numOfBonds; k++)
 //        {
 //            molecules[j].bonds[k] = molec1.bonds[k];
 //        }

 //        //assign angles
 //        for(int k = 0; k < molec1.numOfAngles; k++)
 //        {
 //            molecules[j].angles[k] = molec1.angles[k];
 //        }

 //        //assign dihedrals
 //        for(int k = 0; k < molec1.numOfDihedrals; k++)
 //        {
 //            molecules[j].dihedrals[k] = molec1.dihedrals[k];
 //        }

 //        //assign hops zx add
 //        for(int k = 0; k < molec1.numOfHops; k++)
 //        {
 //            molecules[j].hops[k] = molec1.hops[k];
 //        }
 //    }
   
 //    for(int m = 1; m < molecDiv; m++)
 //    {
 //        int offset=m*molecTypenum;
 //    	memcpy(&molecules[offset],molecules,sizeof(Molecule)*molecTypenum);
 //    	for(int n=0;n<molecTypenum;n++)
 //        {
 //    		molecules[offset+n].moleculeIdentificationNumber=offset+n;
 //            molecules[offset+n].atoms = molecules[n].atoms+count[0]*m;
 //            molecules[offset+n].bonds =  molecules[n].bonds+count[1]*m;
 //            molecules[offset+n].angles =  molecules[n].angles+count[2]*m;
 //            molecules[offset+n].dihedrals =  molecules[n].dihedrals+count[3]*m;
 //            molecules[offset+n].hops =  molecules[n].hops+count[4]*m;
 //        }
        
 //        memcpy(&atoms[offset*count[0]],atoms,sizeof(Atom)*count[0]);
 //        memcpy(&bonds[offset*count[1]],bonds,sizeof(Bond)*count[1]);
 //        memcpy(&angles[offset*count[2]],angles,sizeof(Angle)*count[2]);
 //        memcpy(&dihedrals[offset*count[3]],dihedrals,sizeof(Dihedral)*count[3]);
 //        memcpy(&hops[offset*count[4]],hops,sizeof(Hop)*count[4]);
        
 //        for(int k=0;k<count[0];k++)
 //        {
 //            atoms[offset*count[0]+k].atomIdentificationNumber=offset*count[0]+k;
 //        }
        
 //        for(int k=0;k<count[1];k++)
 //        {
 //            bonds[offset*count[1]+k].atom1+=m*count[0];
 //            bonds[offset*count[1]+k].atom2+=m*count[0];
 //        }
        
 //        for(int k=0;k<count[2];k++)
 //        {
 //            angles[offset*count[2]+k].atom1+=m*count[0];
 //            angles[offset*count[2]+k].atom2+=m*count[0];
 //        }
        
 //        for(int k=0;k<count[3];k++)
 //        {
 //            dihedrals[offset*count[3]+k].atom1+=m*count[0];
 //            dihedrals[offset*count[3]+k].atom2+=m*count[0];
 //        }
        
 //        for(int k=0;k<count[4];k++)
 //        {
 //            hops[offset*count[4]+k].atom1+=m*count[0];
 //            hops[offset*count[4]+k].atom2+=m*count[0];
 //        }
 //    }
     
 //    environment->numOfAtoms = count[0]*molecDiv;
	// ss << "Molecules Created into an Array" << std::endl;
 //    writeToLog(ss.str(),0);
     
 //    if (!configScan.statePath.empty())
 //    {
	// 	ss << "Reading State File \nPath: " << configScan.statePath << std::endl;
 //        std::cout<<ss.str()<<std::endl; writeToLog(ss.str(),0);
 //     	configScan.ReadStateFile(configScan.statePath.c_str(), environment, molecules); //Does this really need the statepath?
 //    }
 //    else
 //    {
 //        ss << "Assigning Molecule Positions..." << std::endl;
 //        writeToLog(ss.str(),0);
 //        //generatefccBox(molecules,environment);//generate fcc lattice box
 //        ss << "Finished Assigning Molecule Positions" << std::endl;
 //        writeToLog(ss.str(),0);
 //    }
}

SerialBox::~SerialBox()
{
	FREE(atoms);
	FREE(environment);
	FREE(molecules);
	FREE(energies);
}