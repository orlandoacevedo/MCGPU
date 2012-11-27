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
#include <assert.h>
#include "SimBox.h"

using namespace std;

#define FREE(ptr) if(ptr!=NULL) { free(ptr);ptr=NULL;}


double randomFloat(const double start, const double end)
{
    return (end-start) * (double(rand()) / RAND_MAX) + start;
}

SimBox::SimBox(Config_Scan configScan)
{
	molecules=NULL;
	enviro=NULL;
	stringstream ss;
	memset(&changedmole,0,sizeof(changedmole));
		
  ss << "Running simulation based on Z-Matrix File"<<endl;
  cout<<ss.str()<<endl; 
  writeToLog(ss);

  //get environment from the config file
  enviro=(Environment *)malloc(sizeof(Environment));
  memcpy(enviro,(configScan.getEnviro()),sizeof(Environment));
  
	ss << "Reading Configuation File \nPath: " << configScan.getConfigPath() << endl;
  cout<<ss.str()<<endl; writeToLog(ss);

  //set up Opls scan 
  ss << "Reading OPLS File \nPath: " << configScan.getOplsusaparPath() << endl;
  cout<<ss.str()<<endl; writeToLog(ss);
	
  string oplsPath = configScan.getOplsusaparPath();
  Opls_Scan oplsScan (oplsPath);
  oplsScan.scanInOpls(oplsPath);
	ss << "OplsScan and OPLS ref table Created " << endl;
  cout<<ss.str()<<endl; writeToLog(ss);

  //set up zMatrixScan
  ss << "Reading Z-Matrix File \nPath: " << configScan.getZmatrixPath() << endl;
  cout<<ss.str()<<endl; writeToLog(ss);
  Zmatrix_Scan zMatrixScan (configScan.getZmatrixPath(), &oplsScan);
  if (zMatrixScan.scanInZmatrix() == -1){
  	ss << "Error, Could not open: " << configScan.getZmatrixPath() << endl;
  	cerr << ss.str()<< endl;
  	writeToLog(ss);
    exit(1);
   }
   
   ss << "Opened Z-Matrix File \nBuilding "<< enviro->numOfMolecules << " Molecules..." << endl;
   cout<<ss.str()<<endl; writeToLog(ss);

   //Convert molecule vectors into an array
   int moleculeIndex = 0;
   int atomCount = 0;

   vector<Molecule> molecVec = zMatrixScan.buildMolecule(atomCount);
   int molecMod = enviro->numOfMolecules % molecVec.size();
   if (molecMod != 0){
       enviro->numOfMolecules += molecVec.size() - molecMod;
       cout << "Number of molecules not divisible by specified z-matrix. Changing number of molecules to: " << enviro->numOfMolecules << endl;
    }
    molecules = (Molecule *)malloc(sizeof(Molecule) * enviro->numOfMolecules);
    
    int molecDiv = enviro->numOfMolecules / molecVec.size();
    int molecTypenum=molecVec.size();
    
    int count[5];//sum up number of atoms,bonds,angles,dihedrals,hops
    memset(count,0,sizeof(count));
    for(int j = 0; j < molecVec.size(); j++){
       Molecule molec1 = molecVec[j];   
       //Copy data from vector to molecule
       count[0]+=molec1.numOfAtoms;
       count[1]+=molec1.numOfBonds;
       count[2]+=molec1.numOfAngles;
       count[3]+=molec1.numOfDihedrals;
       count[4]+=molec1.numOfHops;
     }
     
 	   atompool     =(Atom *)malloc(sizeof(Atom)*molecDiv*count[0]);
 	   bondpool     =(Bond *)malloc(sizeof(Bond)*molecDiv*count[1]);
 	   anglepool    =(Angle *)malloc(sizeof(Angle)*molecDiv*count[2]);
 	   dihedralpool =(Dihedral *)malloc(sizeof(Dihedral)*molecDiv*count[3]);
 	   hoppool      =(Hop *)malloc(sizeof(Hop)*molecDiv*count[4]);  
 	   memset(atompool,0,sizeof(Atom)*molecDiv*count[0]);
 	   memset(bondpool,0,sizeof(Bond)*molecDiv*count[1]);
 	   memset(anglepool,0,sizeof(Angle)*molecDiv*count[2]);
 	   memset(dihedralpool,0,sizeof(Dihedral)*molecDiv*count[3]);
 	   memset(hoppool,0,sizeof(Hop)*molecDiv*count[4]);
 	   
 	   //arrange first part of molecoles
 	   memset(count,0,sizeof(count));
 	   for(int j = 0; j < molecVec.size(); j++){
 	      //Copy data from vector to molecule
        Molecule molec1 = molecVec[j];   

         molecules[j].atoms = (Atom *)(atompool+count[0]);
         molecules[j].bonds = (Bond *)(bondpool+count[1]);
         molecules[j].angles = (Angle *)(anglepool+count[2]);
         molecules[j].dihedrals = (Dihedral *)(dihedralpool+count[3]);
         molecules[j].hops = (Hop *)(hoppool+count[4]);

         molecules[j].id = molec1.id;
         molecules[j].numOfAtoms = molec1.numOfAtoms;
         molecules[j].numOfBonds = molec1.numOfBonds;
         molecules[j].numOfDihedrals = molec1.numOfDihedrals;
         molecules[j].numOfAngles = molec1.numOfAngles;
         molecules[j].numOfHops = molec1.numOfHops;

         count[0]+=molec1.numOfAtoms;
         count[1]+=molec1.numOfBonds;
         count[2]+=molec1.numOfAngles;
         count[3]+=molec1.numOfDihedrals;
         count[4]+=molec1.numOfHops;

         //get the atoms from the vector molecule
         for(int k = 0; k < molec1.numOfAtoms; k++){
             molecules[j].atoms[k] = molec1.atoms[k];
         }               
               
         //assign bonds
         for(int k = 0; k < molec1.numOfBonds; k++){
             molecules[j].bonds[k] = molec1.bonds[k];
         }

         //assign angles
         for(int k = 0; k < molec1.numOfAngles; k++){
             molecules[j].angles[k] = molec1.angles[k];
         }

         //assign dihedrals
         for(int k = 0; k < molec1.numOfDihedrals; k++){
             molecules[j].dihedrals[k] = molec1.dihedrals[k];
         }
         
         //assign hops zx add
         for(int k = 0; k < molec1.numOfHops; k++){
             molecules[j].hops[k] = molec1.hops[k];
         }
      }
   
 	   for(int m = 1; m < molecDiv; m++){
 	      int offset=m*molecTypenum;
 	   		memcpy(&molecules[offset],molecules,sizeof(Molecule)*molecTypenum);
 	   		for(int n=0;n<molecTypenum;n++) {
 	   		    molecules[offset+n].id=offset+n;
            molecules[offset+n].atoms = molecules[n].atoms+count[0]*m;
            molecules[offset+n].bonds =  molecules[n].bonds+count[1]*m;
            molecules[offset+n].angles =  molecules[n].angles+count[2]*m;
            molecules[offset+n].dihedrals =  molecules[n].dihedrals+count[3]*m;
            molecules[offset+n].hops =  molecules[n].hops+count[4]*m;
        }
        
        memcpy(&atompool[offset*count[0]],atompool,sizeof(Atom)*count[0]);
        memcpy(&bondpool[offset*count[1]],bondpool,sizeof(Bond)*count[1]);
        memcpy(&anglepool[offset*count[2]],anglepool,sizeof(Angle)*count[2]);
        memcpy(&dihedralpool[offset*count[3]],dihedralpool,sizeof(Dihedral)*count[3]);
        memcpy(&hoppool[offset*count[4]],hoppool,sizeof(Hop)*count[4]);
        
        for(int k=0;k<count[0];k++)
          atompool[offset*count[0]+k].id=offset*count[0]+k;
        
     
        for(int k=0;k<count[1];k++){
          bondpool[offset*count[1]+k].atom1+=m*count[0];
          bondpool[offset*count[1]+k].atom2+=m*count[0];
        }
        
        for(int k=0;k<count[2];k++) {
          anglepool[offset*count[2]+k].atom1+=m*count[0];
          anglepool[offset*count[2]+k].atom2+=m*count[0];
        }
        
        for(int k=0;k<count[3];k++) {
          dihedralpool[offset*count[3]+k].atom1+=m*count[0];
          dihedralpool[offset*count[3]+k].atom2+=m*count[0];
        }
        
        for(int k=0;k<count[4];k++) {
          hoppool[offset*count[4]+k].atom1+=m*count[0];
          hoppool[offset*count[4]+k].atom2+=m*count[0];
        }
     }
     
     enviro->numOfAtoms = count[0]*molecDiv;
		 ss << "Molecules Created into an Array" << endl;
     writeToLog(ss);
     
     if (!configScan.getStatePath().empty()){
			  ss << "Reading State File \nPath: " << configScan.getStatePath() << endl;
        cout<<ss.str()<<endl; writeToLog(ss);
     	  ReadStateFile(configScan.getStatePath().c_str());
     } else {
        ss << "Assigning Molecule Positions..." << endl;
        writeToLog(ss);
        generatePoints(molecules,enviro);//generate points random
        ss << "Finished Assigning Molecule Positions" << endl;
        writeToLog(ss);
     }
}

SimBox::~SimBox()
{
  FREE(molecules);
  FREE(enviro);
  FREE(atompool);
  FREE(bondpool);
  FREE(anglepool);
  FREE(dihedralpool);
  FREE(hoppool);
 
  //free memory of changedmole
  FREE(changedmole.atoms);
  FREE(changedmole.bonds);
  FREE(changedmole.angles);
  FREE(changedmole.dihedrals);
  FREE(changedmole.hops);
}

int SimBox::ReadStateFile(char const* StateFile)
{
    ifstream inFile;
    Environment tmpenv;
    stringstream ss;
    char buf[250];
    
    cout<<"read state file "<<StateFile<<endl;
    //save current Enviroment to tmpenv at first
    memcpy(&tmpenv,enviro,sizeof(Environment));
    
    inFile.open(StateFile);
    
    //read and check the environment
    if (inFile.is_open()) {
      inFile>>tmpenv.x>>tmpenv.y>>tmpenv.z>>tmpenv.maxTranslation>>tmpenv.numOfAtoms>>tmpenv.temperature>>tmpenv.cutoff;
    }
    
    if (memcmp(&tmpenv,enviro,sizeof(Environment))!=0)
    {
       ss<<"Wrong state files,does not match other configfiles"<<endl;
       ss<<"x "<<tmpenv.x<<" "<<enviro->x<<endl;
       ss<<"y "<<tmpenv.y<<" "<<enviro->y<<endl;
       ss<<"z "<<tmpenv.z<<" "<<enviro->z<<endl;
       ss<<"numOfAtoms "<<tmpenv.numOfAtoms<<" "<<enviro->numOfAtoms<<endl;
       ss<<"temperature "<<tmpenv.temperature<<" "<<enviro->temperature<<endl;
       ss<<"cutoff "<<tmpenv.cutoff<<" "<<enviro->cutoff<<endl;
       ss<<ss.str()<<endl; writeToLog(ss);      
    } 
    inFile.getline(buf,sizeof(buf)); //ignore blank line
    int molecno,atomno;

    molecno=0;
    atomno=0;
    
    int no;
    Atom currentAtom;
   	Bond  currentBond;
 	  Angle currentAngle;
 	  Dihedral currentDi;
 	  Hop      currentHop;
 	  Molecule *ptr=molecules;;

    while(inFile.good()&&molecno<enviro->numOfMolecules){
      inFile>>no	;
      assert(ptr->id==no);
      inFile.getline(buf,sizeof(buf)); //bypass atom flag
      inFile.getline(buf,sizeof(buf));
      assert(strcmp(buf,"= Atoms")==0);
      
      for(int i=0;i<ptr->numOfAtoms;i++)
      {
      	inFile>>currentAtom.id >> currentAtom.x >> currentAtom.y >> currentAtom.z>> currentAtom.sigma >> currentAtom.epsilon >> currentAtom.charge;
      	assert(currentAtom.id==ptr->atoms[i].id);
      	//printf("id:%d,x:%f,y:%f\n",currentAtom.id,currentAtom.x,currentAtom.y);
      	memcpy(&ptr->atoms[i],&currentAtom,sizeof(Atom));
      }

      inFile.getline(buf,sizeof(buf)); //ignore bonds flag
      inFile.getline(buf,sizeof(buf));
      assert(strcmp(buf,"= Bonds")==0);
      for(int i=0;i<ptr->numOfBonds;i++)
      {
      	inFile>>currentBond.atom1 >>currentBond.atom2 >> currentBond.distance >> currentBond.variable;
      	assert(currentBond.atom1==ptr->bonds[i].atom1);
				assert(currentBond.atom2==ptr->bonds[i].atom2);      	
      	memcpy(&ptr->bonds[i],&currentBond,sizeof(Bond));
      }

      inFile.getline(buf,sizeof(buf)); //ignore Dihedrals flag
      inFile.getline(buf,sizeof(buf));
      assert(strcmp(buf,"= Dihedrals")==0);
      for(int i=0;i<ptr->numOfDihedrals;i++)
      {
      	inFile>>currentDi.atom1>>currentDi.atom2>>currentDi.value>>currentDi.variable;
      	assert(currentDi.atom1==ptr->dihedrals[i].atom1);
				assert(currentDi.atom2==ptr->dihedrals[i].atom2);      	
      	memcpy(&ptr->dihedrals[i],&currentDi,sizeof(Dihedral));
      }

      inFile.getline(buf,sizeof(buf)); //ignore hops flag
      inFile.getline(buf,sizeof(buf));
      assert(strcmp(buf,"=Hops")==0);
      for(int i=0;i<ptr->numOfHops;i++)
      {
      	inFile>>currentHop.atom1>>currentHop.atom2 >>currentHop.hop;
      	assert(currentHop.atom1==ptr->hops[i].atom1);
				assert(currentHop.atom2==ptr->hops[i].atom2);      	
      	memcpy(&ptr->hops[i],&currentHop,sizeof(Hop));
      }

      inFile.getline(buf,sizeof(buf)); //ignore angles flag
      inFile.getline(buf,sizeof(buf));
      assert(strcmp(buf,"= Angles")==0);
      for(int i=0;i<ptr->numOfAngles;i++)
      {
      	inFile>>currentAngle.atom1 >> currentAngle.atom2 >>currentAngle.value >>currentAngle.variable;
      	assert(currentAngle.atom1==ptr->angles[i].atom1);
				assert(currentAngle.atom2==ptr->angles[i].atom2);      	
      	memcpy(&ptr->angles[i],&currentAngle,sizeof(Angle));
      }       
      
      inFile.getline(buf,sizeof(buf)); //bypass == flag
      inFile.getline(buf,sizeof(buf));
      assert(strcmp(buf,"==")==0);   
      
      ptr++;                    
      molecno++;
    }
  inFile.close();
  WriteStateFile("Confirm.state");

	return 0;
}

int SimBox::WriteStateFile(char const* StateFile)
{
    ofstream outFile;
    int numOfMolecules=enviro->numOfMolecules;
    
    outFile.open(StateFile);
    
    //print the environment
    outFile << enviro->x << " " << enviro->y << " " << enviro->z << " " << enviro->maxTranslation<<" " << enviro->numOfAtoms
        << " " << enviro->temperature << " " << enviro->cutoff <<endl;
    outFile << endl; // blank line
    
    for(int i = 0; i < numOfMolecules; i++){
        Molecule currentMol = molecules[i];
        outFile << currentMol.id << endl;
        outFile << "= Atoms" << endl; // delimiter
    
        //write atoms
        for(int j = 0; j < currentMol.numOfAtoms; j++){
            Atom currentAtom = currentMol.atoms[j];
            outFile << currentAtom.id << " "
                << currentAtom.x << " " << currentAtom.y << " " << currentAtom.z
                << " " << currentAtom.sigma << " " << currentAtom.epsilon  << " "
                << currentAtom.charge << endl;
        }
        outFile << "= Bonds" << endl; // delimiter
        
        //write bonds
        for(int j = 0; j < currentMol.numOfBonds; j++){
            Bond currentBond = currentMol.bonds[j];
            outFile << currentBond.atom1 << " " << currentBond.atom2 << " "
                << currentBond.distance << " ";
            if(currentBond.variable)
                outFile << "1" << endl;
            else
                outFile << "0" << endl;

        }
        outFile << "= Dihedrals" << endl; // delimiter
        for(int j = 0; j < currentMol.numOfDihedrals; j++){
            Dihedral currentDi = currentMol.dihedrals[j];
            outFile << currentDi.atom1 << " " << currentDi.atom2 << " "
                << currentDi.value << " ";
            if(currentDi.variable)
                outFile << "1" << endl;
            else
                outFile << "0" << endl;
        }

        outFile << "=Hops" << endl;

        for(int j = 0; j < currentMol.numOfHops; j++){
            Hop currentHop = currentMol.hops[j];

            outFile << currentHop.atom1 << " " << currentHop.atom2 << " "
                << currentHop.hop << endl;
        }
        
        
        outFile << "= Angles" << endl; // delimiter

        //print angless
        for(int j = 0; j < currentMol.numOfAngles; j++){
            Angle currentAngle = currentMol.angles[j];

            outFile << currentAngle.atom1 << " " << currentAngle.atom2 << " "
                << currentAngle.value << " ";
            if(currentAngle.variable)
                outFile << "1" << endl;
            else
                outFile << "0" << endl;
        }


        //write a == line
        outFile << "==" << endl;
    }
    outFile.close();
	return 0;
}

int SimBox::writePDB(char const* pdbFile)
{
	printf("%s\n",pdbFile);
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

    //zx mod for global seed used srand((unsigned int) time(NULL));
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
    Molecule molec = molecules[currentIndex];
	int molecId = molec.atoms[0].id;
    while(atomId < molecId && currentIndex > 0){
        currentIndex -= 1;
		molec = molecules[currentIndex];
		molecId = molec.atoms[0].id;
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
    Hop *myHop = molecule->hops;
    for(int x=0; x< molecule->numOfHops; x++){        
        //compare atoms to each hop struct in molecule
		if((myHop->atom1==atom1 && myHop->atom2==atom2) || (myHop->atom1==atom2 && myHop->atom2==atom1)){
            return myHop->hop;
        }
    myHop++;
	}
	 return 0;
}

int SimBox::ChangeMolecule()
{
    double maxTranslation = enviro->maxTranslation;
    double maxRotation = enviro->maxRotation;

    //Pick a molecule to move
    int moleculeIndex = randomFloat(0, enviro->numOfMolecules);
        
    saveChangedMole(moleculeIndex);
        
   //Pick an atom in the molecule about which to rotate
   int atomIndex = randomFloat(0, molecules[moleculeIndex].numOfAtoms);
   Atom vertex = molecules[moleculeIndex].atoms[atomIndex];

   const double deltaX = randomFloat(-maxTranslation, maxTranslation);
   const double deltaY = randomFloat(-maxTranslation, maxTranslation);
   const double deltaZ = randomFloat(-maxTranslation, maxTranslation);

   const double degreesX = randomFloat(-maxRotation, maxRotation);
   const double degreesY = randomFloat(-maxRotation, maxRotation);
   const double degreesZ = randomFloat(-maxRotation, maxRotation); 

   moveMolecule(molecules[moleculeIndex], vertex, deltaX, deltaY, deltaZ,
        degreesX, degreesY, degreesZ);

   keepMoleculeInBox(&molecules[moleculeIndex], enviro);

   return moleculeIndex;
}

int SimBox::Rollback(int moleno)
{
	 return copyMolecule(&molecules[moleno],&changedmole);
}

int SimBox::saveChangedMole(int moleno)
{
  Molecule *mole_src=&molecules[moleno];
  
  //free memory of changedmole before allocate memory
  FREE(changedmole.atoms);
  FREE(changedmole.bonds);
  FREE(changedmole.angles);
  FREE(changedmole.dihedrals);
  FREE(changedmole.hops);

  memcpy(&changedmole,mole_src,sizeof(changedmole));
  
  changedmole.atoms = (Atom *)malloc(sizeof(Atom) * mole_src->numOfAtoms);
  changedmole.bonds = (Bond *)malloc(sizeof(Bond) * mole_src->numOfBonds);
  changedmole.angles = (Angle *)malloc(sizeof(Angle) * mole_src->numOfAngles);
  changedmole.dihedrals = (Dihedral *)malloc(sizeof(Dihedral) * mole_src->numOfDihedrals);
	changedmole.hops = (Hop *)malloc(sizeof(Hop) * mole_src->numOfHops);
	
	copyMolecule(&changedmole,mole_src);
	
	return 0;
}

int SimBox::copyMolecule(Molecule *mole_dst, Molecule *mole_src)
{
    mole_dst->numOfAtoms = mole_src->numOfAtoms;
    mole_dst->numOfBonds = mole_src->numOfBonds;
    mole_dst->numOfAngles = mole_src->numOfAngles;
    mole_dst->numOfDihedrals = mole_src->numOfDihedrals;
    mole_dst->numOfHops =  mole_src->numOfHops;
    mole_dst->id = mole_src->id;
    
  for(int i = 0; i < mole_src->numOfAtoms; i++){
        mole_dst->atoms[i] = mole_src->atoms[i];
  }

  for(int i = 0; i < mole_src->numOfBonds; i++){
        mole_dst->bonds[i] = mole_src->bonds[i];
  }
  
  for(int i = 0; i < mole_src->numOfAngles; i++){
        mole_dst->angles[i] = mole_src->angles[i];
  }
  
  for(int i = 0; i < mole_src->numOfDihedrals; i++){
        mole_dst->dihedrals[i] = mole_src->dihedrals[i];
  }
	
	for(int i = 0; i < mole_src->numOfHops; i++){
        mole_dst->hops[i] = mole_src->hops[i];
  }
  
  return 0;  	
}
