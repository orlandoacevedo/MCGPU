/*Intended goal: support read, parse, and extract operations on configuration files to properly initialize 
*  a simulation environment.
* [1] This file/class will, ideally, replace Config_Scan.cpp & will augment MetroUtil.cpp. (Started on 19 Feb. Beta completion on 28 Feb.)
* [2] This file/class will, ideally, replace Opls_Scan and Zmatrix_Scan. (Started, 26 Feb. -Albert)
*Created 19 February 2014. Albert Wallace
*/
/*
--Changes made on:
		->Sun, 23 Feb 2014 (Albert)
		->Wed, 26 Feb (Albert)
		->Thu, 27 Feb (Albert, then Tavis)
		->Fri, 28 Feb (Albert)
*/
/*Based on work from earlier sessions by Alexander Luchs, Riley Spahn, Seth Wooten, and Orlando Acevedo*/

//Dear Self: Tavis says we might be better off importing an instance of a box and filling it with the necessary information.
// He has good ideas and we should listen to him.
// Also, do destructors. Please.

//_________________________________________________________________________________________________________________
//  INCLUDE statements
//_________________________________________________________________________________________________________________
#include <assert.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <errno.h>
#include <exception>
#include <stdexcept>
#include <vector>

#include "StructLibrary.h"
#include "IOUtilities.h"

//_________________________________________________________________________________________________________________
//  DEFINE statements
//_________________________________________________________________________________________________________________
// These are required for the associated SWITCH statement during the writing of the log file
#define DEFAULT 0
#define START 1
#define END 2
#define OPLS 3
#define Z_MATRIX 4
#define GEOM 5

//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//  Initial construction & reading of configuration path information.
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________

/*
*Constructor for the entirety of the IOUtilities class. Will call various methods toward the end to facilitate an easy
*	simulation environment setup.
*
*@params: configPath: the path to the main configuration file, which itself points to other specialized configuration files
* 	and contains other bits of information to set up a proper Environment for simulation.
*/
IOUtilities::IOUtilities(std::string configPath){


	//UtilitiesInfo filePathsEtc; //all the variables used for this class are stuck in this struct for easy, yet unsafe, access
					//there were apparently getters and setters to access them all, so if necessary, we can have one getter later for the entire struct
	
	//note to people/myself: //enviro = currentEnvironment
	
    //memset(&filePathsEtc,0,sizeof(UtilitiesInfo)); //filePathsEtc is a struct of type UtilitiesInfo, and this is apparently the best way to instantiate the struct
    			//except that may not be required? but it's left in for legacy reasons
    filePathsEtc = new UtilitiesInfo();
    filePathsEtc->configPath = configPath;
    filePathsEtc->numOfSteps=0;
    readInConfigAlreadyDone = false;
    readInConfig(); //do the rest of the construction
    readInConfigAlreadyDone = true; //setting it this way will prevent unnecessarily running the entirety of readInConfig
    			//though this setup means we cannot change the config file during the run, but this is by design.
}

/*
*Primary method called to read from the main configuration file. Parses all known information in the predefined format
*	and stores it into an Environment structure, or otherwise in the UtilitiesInfo structure which holds various file paths, etc.
*
*@params: [none; uses variables within the class to pass information]
*@return: [none; uses variables within the class to pass information]
*/
void IOUtilities::readInConfig()
{
	if (readInConfigAlreadyDone)
	{
		//it's already been done during construction; do nothing
		}
	else
	{
		std::ifstream configscanner(filePathsEtc->configPath.c_str());
		if (! configscanner.is_open())
		{
			throwScanError("Configuration file failed to open.");
			return;
		}
		else
		{
			std::string line;
			int currentLine = 1;
			while (configscanner.good())
			{
				std::getline(configscanner,line);
			
				//assigns attributes based on line number
				//current line = actual line in the configuration file
				switch(currentLine)
				{
					case 2:
						if(line.length() > 0)
						{
							filePathsEtc->currentEnvironment->x = atof(line.c_str());
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing environment x value.");
							return;
						}
						break;
					case 3:
						if(line.length() > 0)
						{
							filePathsEtc->currentEnvironment->y = atof(line.c_str());
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing environment y value.");
							return;
						}
						break;
					case 4:
						if(line.length() > 0)
						{
							filePathsEtc->currentEnvironment->z = atof(line.c_str());
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing environment z value.");
							return;
						}
						break;
					case 6:
						if(line.length() > 0)
						{
							filePathsEtc->currentEnvironment->temp = atof(line.c_str());
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing environment temperature value.");
							return;
						}
						break;
					case 8:
						if(line.length() > 0)
						{
							filePathsEtc->currentEnvironment->maxTranslation = atof(line.c_str());
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing environment max translation value.");
							return;
						}
						break;
					case 10:
						if(line.length() > 0)
						{
							filePathsEtc->numOfSteps = atoi(line.c_str());
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing number of steps value.");
							return;
						}
						break;
					case 12:
						if(line.length() > 0)
						{
							filePathsEtc->currentEnvironment->numOfMolecules = atoi(line.c_str());
							//printf("number is %d",enviro.numOfMolecules);
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing number of molecules value.");
							return;
						}
						break;
					case 14:
						if(line.length() > 0)
						{
							filePathsEtc->oplsuaparPath = line;
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing oplsuapar path value.");
							return;
						}
						break;
					case 16:
						if(line.length() > 0)
						{
							filePathsEtc->zmatrixPath = line;
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing z-matrix path value.");
							return;
						}
						break;
					case 18:
						if(line.length() > 0)
						{
							filePathsEtc->statePath = line;
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing state file output path value.");
							return;
						}
						break;
					case 20:
						if(line.length() > 0){
							filePathsEtc->stateOutputPath = line;
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing state file output path value.");
							return;
						}
						break;
					case 22:
						if(line.length() > 0){
							filePathsEtc->pdbOutputPath = line;
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing PDB output path value.");
							return;
						}
						break;
					case 24:
						if(line.length() > 0)
						{
							filePathsEtc->currentEnvironment->cutoff = atof(line.c_str());
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing environment cutoff value.");
							return;
						}
						break;
					case 26:
						if(line.length() > 0)
						{
							filePathsEtc->currentEnvironment->maxRotation = atof(line.c_str());
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing environment max rotation value.");
							return;
						}
						break;
					case 28:
						if(line.length() > 0)
						{
							filePathsEtc->currentEnvironment->randomseed=atoi(line.c_str());
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing random seed value.");
							return;
						}
						break;
					case 30:
						if(line.length() > 0)
						{
							// Convert to a zero-based index
							filePathsEtc->currentEnvironment->primaryAtomIndex=atoi(line.c_str()) - 1;
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing environment primary atom index value.");
							return;
						}
						break;
						//end of disabled configuration file code.
				}
			
				currentLine++;
			}
		}
	}
}

void IOUtilities::throwScanError(std::string message)
{

	std::cerr << std::endl << message << std::endl << "	Error Number: " << errno << std::endl << std::endl;

	return;
}


/**
	This method is used to read in from a state file.
	//copied over from...SimBox, so retooling is necessary (variable references point to SimBox locations)
	
	@param StateFile - takes the location of the state file to be read in
*/
/*
int IOUtilities::ReadStateFile(char const* StateFile, Environment * destinationEnvironment, Molecule * destinationMoleculeCollection)
{
    ifstream inFile;
    Environment tmpenv;
    stringstream ss;
    char buf[250];
    
    cout<<"read state file "<<StateFile<<endl;
    //save current Enviroment to tmpenv at first
    memcpy(&tmpenv,destinationEnvironment,sizeof(Environment));
    
    inFile.open(StateFile);
    
    //read and check the environment
    if (inFile.is_open())
    {
      inFile>>tmpenv.x>>tmpenv.y>>tmpenv.z>>tmpenv.maxTranslation>>tmpenv.numOfAtoms>>tmpenv.temp>>tmpenv.cutoff;
    }
    
    if (memcmp(&tmpenv,destinationEnvironment,sizeof(Environment))!=0)
    {
       ss<<"Wrong state files,does not match other configfiles"<<endl;
       ss<<"x "<<tmpenv.x<<" "<<destinationEnvironment->x<<endl;
       ss<<"y "<<tmpenv.y<<" "<<destinationEnvironment->y<<endl;
       ss<<"z "<<tmpenv.z<<" "<<destinationEnvironment->z<<endl;
       ss<<"numOfAtoms "<<tmpenv.numOfAtoms<<" "<<destinationEnvironment->numOfAtoms<<endl;
       ss<<"temperature "<<tmpenv.temp<<" "<<destinationEnvironment->temp<<endl;
       ss<<"cutoff "<<tmpenv.cutoff<<" "<<destinationEnvironment->cutoff<<endl;
       ss<<ss.str()<<endl; writeToLog(ss);      
    } 
    inFile.getline(buf,sizeof(buf)); //ignore blank line
    int molecno = 0;
    //int atomno = 0; //reported at compile-time as being unused, so commented out

    int no;
    Atom currentAtom;
   	Bond  currentBond;
 	Angle currentAngle;
    Dihedral currentDi;
 	Hop      currentHop;
 	Molecule *ptr=destinationMoleculeCollection;

    while(inFile.good()&&molecno<destinationEnvironment->numOfMolecules)
    {
        inFile>>no;
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
        // known BUG - if molecule has no hops (3 atoms or less) state file gives error crashing simulation
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
    WriteStateFile("Confirm.state", destinationEnvironment, ptr);

	return 0;
}

*/


/**
	Used to write to a state file.
	//copied over from...Simbox. Retooling needed. (Variable references point to SimBox locations)
	
	@param StateFile - writes to a state file at the location given
*/
/*
int IOUtilities::WriteStateFile(char const* StateFile, Environment * sourceEnvironment, Molecule * sourceMoleculeCollection)
{
    ofstream outFile;
    int numOfMolecules=sourceEnvironment->numOfMolecules;
    
    outFile.open(StateFile);
    
    //print the environment
    outFile << sourceEnvironment->x << " " << sourceEnvironment->y << " " << sourceEnvironment->z << " " << sourceEnvironment->maxTranslation<<" " << sourceEnvironment->numOfAtoms
        << " " << sourceEnvironment->temperature << " " << sourceEnvironment->cutoff <<endl;
    outFile << endl; // blank line
    
    for(int i = 0; i < numOfMolecules; i++)
    {
        Molecule currentMol = sourceMoleculeCollections[i];
        outFile << currentMol.id << endl;
        outFile << "= Atoms" << endl; // delimiter
    
        //write atoms
        for(int j = 0; j < currentMol.numOfAtoms; j++)
        {
            Atom currentAtom = currentMol.atoms[j];
            outFile << currentAtom.id << " "
                << currentAtom.x << " " << currentAtom.y << " " << currentAtom.z
                << " " << currentAtom.sigma << " " << currentAtom.epsilon  << " "
                << currentAtom.charge << endl;
        }
        outFile << "= Bonds" << endl; // delimiter
        
        //write bonds
        for(int j = 0; j < currentMol.numOfBonds; j++)
        {
            Bond currentBond = currentMol.bonds[j];
            outFile << currentBond.atom1 << " " << currentBond.atom2 << " "
                << currentBond.distance << " ";
            if(currentBond.variable)
                outFile << "1" << endl;
            else
                outFile << "0" << endl;

        }
        outFile << "= Dihedrals" << endl; // delimiter
        for(int j = 0; j < currentMol.numOfDihedrals; j++)
        {
            Dihedral currentDi = currentMol.dihedrals[j];
            outFile << currentDi.atom1 << " " << currentDi.atom2 << " "
                << currentDi.value << " ";
            if(currentDi.variable)
            {
                outFile << "1" << endl;
            }
            else
            {
                outFile << "0" << endl;
            }
        }

        outFile << "=Hops" << endl;

        for(int j = 0; j < currentMol.numOfHops; j++)
        {
            Hop currentHop = currentMol.hops[j];

            outFile << currentHop.atom1 << " " << currentHop.atom2 << " "
                << currentHop.hop << endl;
        }
        
        
        outFile << "= Angles" << endl; // delimiter

        //print angless
        for(int j = 0; j < currentMol.numOfAngles; j++)
        {
            Angle currentAngle = currentMol.angles[j];

            outFile << currentAngle.atom1 << " " << currentAngle.atom2 << " "
                << currentAngle.value << " ";
            if(currentAngle.variable)
            {
                outFile << "1" << endl;
            }
            else
            {
                outFile << "0" << endl;
            }
        }


        //write a == line
        outFile << "==" << endl;
    }
    outFile.close();
	return 0;
}
*/




/**
	writes to a PDB file for visualizing the box
	//this method copied over from...Simbox. Retooling needed. (Variable references point to SimBox locations)
	
	//other versions found in metroUtil were not used; it relied only on SimBox
	
	//to do this, the main simulation controller--linearSimulationExample in the old data structure--will have to manage the data in the
	// PDB, by directly modifying it in the instance of the IOUtilities class we will have to make,
	// rather than SimBox managing it.
	
	// the downside is that this implementation, when done for parallel,
	//  will likely require writing back to the host, rather than storing the info on the GPU until the very end...
	//  as the SimBox holds all the current information being accessed by this method.
	
	@param pdbFile - Location of the pdbFile to be written to
*/
/*
int IOUtilities::writePDB(char const* pdbFile)
{
    ofstream outputFile;
    outputFile.open(pdbFile);
    int numOfMolecules=enviro->numOfMolecules;
    outputFile << "REMARK Created by MCGPU" << endl;
    //int atomIndex = 0;
    for (int i = 0; i < numOfMolecules; i++)
    {
    	Molecule currentMol = molecules[i];    	
        for (int j = 0; j < currentMol.numOfAtoms; j++)
        {
        	Atom currentAtom = currentMol.atoms[j];
            outputFile.setf(ios_base::left,ios_base::adjustfield);
            outputFile.width(6);
            outputFile << "ATOM";
            outputFile.setf(ios_base::right,ios_base::adjustfield);
            outputFile.width(5);
            outputFile << currentAtom.id + 1;
            outputFile.width(3); // change from 5
            outputFile << currentAtom.name;
            outputFile.width(6); // change from 4
            outputFile << "UNK";
            outputFile.width(6);
            outputFile << i + 1;
            outputFile.setf(ios_base::fixed, ios_base::floatfield);
            outputFile.precision(3);
            outputFile.width(12);
            outputFile << currentAtom.x;
            outputFile.width(8);
            outputFile << currentAtom.y;
            outputFile.width(8);
            outputFile << currentAtom.z << endl;
        }
        outputFile << "TER" << endl;
    }
    outputFile << "END" << endl;
    outputFile.close();

	return 0;
}
*/

//potential overloading of the above function

//@param: sourceMoleculeCollection: Is an array of molecules, dynamically allocated & created elsewhere [such as in the SimBox]

int IOUtilities::writePDB(char const* pdbFile, Environment sourceEnvironment, Molecule * sourceMoleculeCollection)
{
	//molecules = (Molecule *)malloc(sizeof(Molecule) * enviro->numOfMolecules);
	
    std::ofstream outputFile;
    outputFile.open(pdbFile);
    int numOfMolecules=sourceEnvironment.numOfMolecules;
    outputFile << "REMARK Created by MCGPU" << std::endl;
    //int atomIndex = 0;
    for (int i = 0; i < numOfMolecules; i++)
    {
    	Molecule currentMol = sourceMoleculeCollection[i];    	
        for (int j = 0; j < currentMol.numAtoms; j++)
        {
        	Atom currentAtom = currentMol.atoms[j];
            outputFile.setf(std::ios_base::left,std::ios_base::adjustfield);
            outputFile.width(6);
            outputFile << "ATOM";
            outputFile.setf(std::ios_base::right,std::ios_base::adjustfield);
            outputFile.width(5);
            outputFile << currentAtom.id + 1;
            outputFile.width(3); // change from 5
            outputFile << currentAtom.name;
            outputFile.width(6); // change from 4
            outputFile << "UNK";
            outputFile.width(6);
            outputFile << i + 1;
            outputFile.setf(std::ios_base::fixed, std::ios_base::floatfield);
            outputFile.precision(3);
            outputFile.width(12);
            outputFile << currentAtom.x;
            outputFile.width(8);
            outputFile << currentAtom.y;
            outputFile.width(8);
            outputFile << currentAtom.z << std::endl;
        }
        outputFile << "TER" << std::endl;
    }
    outputFile << "END" << std::endl;
    outputFile.close();

	return 0;
}

//_________________________________________________________________________________________________________________
//  OPLS configuration scanning. [from Opls_Scan.cpp]
//_________________________________________________________________________________________________________________

/*Opls_Scan::Opls_Scan(string filename)
{
   fileName = filename;
}*/

void IOUtilities::deleteOpls_Scan()
{
    oplsTable.clear();
}

int IOUtilities::scanInOpls()
{
    int numOfLines=0;
    std::ifstream oplsScanner(filePathsEtc->oplsuaparPath.c_str()); //##
    if( !oplsScanner.is_open() ) //##
        return -1; //##
    else { //##
        std::string line;  //##
        while( oplsScanner.good() ) //##
        { //##
            numOfLines++; //##
            std::getline(oplsScanner,line); //##

            //check if it is a commented line,
            //or if it is a title line
            try{ //##
                if(line.at(0) != '#' && numOfLines > 1)
                    addLineToTable(line,numOfLines);
            }
            catch (std::out_of_range& e){}
        }
        oplsScanner.close();
        logErrors();
    }
    return 0;
}

void IOUtilities::addLineToTable(string line, int numOfLines) //##
{
    std::string hashNum;
    int secCol;
    double charge,sigma,epsilon;
    std::string name, extra;
    std::stringstream ss(line);

    //check to see what format it is opls, V value, or neither
    int format = checkFormat(line);
              
    if(format == 1)
    {      	
        ss >> hashNum >> secCol >> name >> charge >> sigma >> epsilon;
        char *atomtype = (char*)name.c_str(); 
         
        Atom temp = createAtom(0, -1, -1, -1, sigma, epsilon, charge, *atomtype);
        std::pair<map<string,Atom>::iterator,bool> ret;
        ret = oplsTable.insert( std::pair<string,Atom>(hashNum,temp) );

        if (ret.second==false)
        {
            errHashes.push_back(hashNum);
        }
    }
    else if(format == 2)
    {
        double v0,v1,v2,v3;
        ss >> hashNum >> v0 >> v1 >> v2 >> v3 ;
        Fourier vValues = {v0,v1,v2,v3};
        std::pair<map<string,Fourier>::iterator,bool> ret2;
        ret2 = fourierTable.insert( std::pair<string,Fourier>(hashNum,vValues) );

        if (ret2.second==false)
        {
            errHashesFourier.push_back(hashNum);	
        }	  
    }
    else
    {
	     errLinesOPLS.push_back(numOfLines);
    }
}

int IOUtilities::checkFormat(std::string line)
{   	 
    int hashNum, secCol;
    double charge,sigma,epsilon;
    std::string name, extra;
    std::stringstream iss(line);

    double v1,v2,v3,v4;
    std::stringstream issw(line);

    //see if format is the V values for the diherdral format
    if((issw >> hashNum >> v1 >> v2 >> v3 >> v4) )
    {
        return 2;
    }
    //else see if format is normal opls format
    else if((iss >> hashNum >> secCol >> name >> charge >> sigma >> epsilon ))
    {
        return 1;
    }
    //if neither return -1
    else
    {
        return -1;
    }
}

void IOUtilities::logErrors()
{
    std::stringstream output;
    // See if there were any errors
    if(errLinesOPLS.empty() || errHashes.empty()|| errHashesFourier.empty())
    {
	    //Errors in the format
		output<<"Errors found in the OPLS file: "<< filePathsEtc->oplsuaparPath.c_str() <<std::endl;
        if(!errLinesOPLS.empty())
        {
		      output << "Found Errors in the Format of the following Lines: " << std::endl;
				for(int a=0; a<errLinesOPLS.size(); a++)
                {
				    if(a%10==0 && a!=0) //ten per line
                    {
					     output << std::endl;
                    }
				    output << errLinesOPLS[a]<< " ";
				}
				output << std::endl << std::endl;
		  }
		  if(!errHashes.empty())
          {
		      output << "Error - The following OPLS values existed more than once: " << std::endl;
				for(int a=0; a<errHashes.size(); a++)
                {
				    if(a%10==0 && a!=0) //ten per line
                    {
					     output << std::endl;
                    }
				    output << errHashes[a]<< " ";
				}
				output << std::endl << std::endl;
		  }
		  if(!errHashesFourier.empty())
          {
		      output << "Error - The following Fourier Coefficent values existed more than once: " << std::endl;
				for(int a=0; a<errHashesFourier.size(); a++)
                {
				    if(a%10==0 && a!=0) //ten per line
                    {
					     output << std::endl;
                    }
				    output << errHashesFourier[a]<< " ";
				}
				output << std::endl << std::endl;
		  }
		  writeToLog(output,OPLS);
	}
}

Atom IOUtilities::getAtom(std::string hashNum)
{
    if(oplsTable.count(hashNum)>0 )
    {
        return oplsTable[hashNum];
	}
	else
    {
	    std::cerr << "Index does not exist: "<< hashNum << std::endl;
		return createAtom(0, -1, -1, -1, -1, -1, -1, NULL);
	}
}

double IOUtilities::getSigma(std::string hashNum)
{
    if(oplsTable.count(hashNum)>0 )
    {
        Atom temp = oplsTable[hashNum];
        return temp.sigma;
    }
    else
    {
        std::cerr << "Index does not exist: "<< hashNum << std::endl;
        return -1;
    }
}

double IOUtilities::getEpsilon(std::string hashNum)
{
    if(oplsTable.count(hashNum)>0 )
    {
        Atom temp = oplsTable[hashNum];
        return temp.epsilon;
    }
    else
    {
        std::cerr << "Index does not exist: "<< hashNum << std::endl;
        return -1;
    }
}

double IOUtilities::getCharge(std::string hashNum)
{
    if(oplsTable.count(hashNum)>0 )
    {
        Atom temp = oplsTable[hashNum];
        return temp.charge;
    }
    else
    {
        std::cerr << "Index does not exist: "<< hashNum << std::endl;
        return -1;
    }
}

Fourier IOUtilities::getFourier(std::string hashNum)
{
    if(fourierTable.count(hashNum)>0 )
    {
        Fourier temp = fourierTable[hashNum];
        return temp;
    }
    else
    {	    
        std::cerr << "Index does not exist: "<< hashNum << std::endl;
        Fourier temp ={-1,-1,-1,-1};
        return temp;
    }
}




/*
//_________________________________________________________________________________________________________________
//  Zmatrix opening & parsing/Zmatrix configuration. [from Zmatrix_Scan.cpp]
//_________________________________________________________________________________________________________________
Zmatrix_Scan::Zmatrix_Scan(string filename, Opls_Scan* oplsScannerRef)
{
    fileName = filename;
    oplsScanner = oplsScannerRef;
    startNewMolecule = false;
}

Zmatrix_Scan::~Zmatrix_Scan(){}

int Zmatrix_Scan::scanInZmatrix()
{
    stringstream output;
    int numOfLines=0;
    ifstream zmatrixScanner(fileName.c_str());

    if( !zmatrixScanner.is_open() )
    {
        return -1;
    }
    else
    {
        string line; 
        while( zmatrixScanner.good() )
        {
            numOfLines++;
            getline(zmatrixScanner,line);

            Molecule workingMolecule;

            //check if it is a commented line,
            //or if it is a title line
            try
            {
                if(line.at(0) != '#' && numOfLines > 1)
                {
                    parseLine(line,numOfLines);
                }
            }
            catch(std::out_of_range& e){}

            if (startNewMolecule)
            {
                Atom* atomArray;
                Bond* bondArray;
                Angle* angleArray;
                Dihedral* dihedralArray;
                
                atomArray = (Atom*) malloc(sizeof(Atom) * atomVector.size());
                bondArray = (Bond*) malloc(sizeof(Bond) * bondVector.size());
                angleArray = (Angle*) malloc(sizeof(Angle) * angleVector.size());
                dihedralArray = (Dihedral*) malloc(sizeof(Dihedral) * dihedralVector.size());

                for (int i = 0; i < atomVector.size(); i++)
                {
                    atomArray[i] = atomVector[i];
                }
                for (int i = 0; i < bondVector.size(); i++)
                {
                    bondArray[i] = bondVector[i];
                }
                for (int i = 0; i < angleVector.size(); i++)
                {
                    angleArray[i] = angleVector[i];
                }
                for (int i = 0; i < dihedralVector.size(); i++)
                {
                    dihedralArray[i] = dihedralVector[i];
                }

                moleculePattern.push_back(createMolecule(-1, atomArray, angleArray, bondArray, dihedralArray, 
                     atomVector.size(), angleVector.size(), bondVector.size(), dihedralVector.size()));

                atomVector.clear();
                bondVector.clear();
                angleVector.clear();
                dihedralVector.clear();

                startNewMolecule = false;
            } 
        }

        zmatrixScanner.close();
    	 
    }

    return 0;
}

void Zmatrix_Scan::parseLine(string line, int numOfLines)
{

    string atomID, atomType, oplsA, oplsB, bondWith, bondDistance, angleWith, angleMeasure, dihedralWith, dihedralMeasure;

    stringstream ss;

    //check if line contains correct format
    int format = checkFormat(line);

    if(format == 1)
    {
        //read in strings in columns and store the data in temporary variables
        ss << line;    	
        ss >> atomID >> atomType >> oplsA >> oplsB >> bondWith >> bondDistance >> angleWith >> angleMeasure >> dihedralWith >> dihedralMeasure;

        //setup structures for permanent encapsulation
        Atom lineAtom;
        Bond lineBond;
        Angle lineAngle;
        Dihedral lineDihedral;
		  
        if (oplsA.compare("-1") != 0)
        {
            lineAtom = oplsScanner->getAtom(oplsA);
            lineAtom.id = atoi(atomID.c_str());
            lineAtom.x = 0;
            lineAtom.y = 0;
            lineAtom.z = 0;
        }
        else//dummy atom
        {
        	char dummy = 'X';
            lineAtom = createAtom(atoi(atomID.c_str()), -1, -1, -1, -1, -1, -1, dummy);
        }
		  atomVector.push_back(lineAtom);

        if (bondWith.compare("0") != 0)
        {
            lineBond.atom1 = lineAtom.id;
            lineBond.atom2 = atoi(bondWith.c_str());
            lineBond.distance = atof(bondDistance.c_str());
            lineBond.variable = false;
            bondVector.push_back(lineBond);
        }

        if (angleWith.compare("0") != 0)
        {
            lineAngle.atom1 = lineAtom.id;
            lineAngle.atom2 = atoi(angleWith.c_str());
            lineAngle.value = atof(angleMeasure.c_str());
            lineAngle.variable = false;
            angleVector.push_back(lineAngle);
        }

        if (dihedralWith.compare("0") != 0)
        {
            lineDihedral.atom1 = lineAtom.id;
            lineDihedral.atom2 = atoi(dihedralWith.c_str());
            lineDihedral.value = atof(dihedralMeasure.c_str());
            lineDihedral.variable = false;
            dihedralVector.push_back(lineDihedral);
        }
    } //end if format == 1

    else if(format == 2)
    {
        startNewMolecule = true;
	}
    else if(format == 3)
    {
        startNewMolecule = true;
    }
    if (previousFormat >= 3 && format == -1)
    {
        handleZAdditions(line, previousFormat);
    }

    previousFormat = format;
}
    
int Zmatrix_Scan::checkFormat(string line)
{
    int format =-1; 
    stringstream iss(line);
    stringstream iss2(line);
    string atomType, someLine;
    int atomID, oplsA, oplsB, bondWith, angleWith,dihedralWith,extra;
    double bondDistance, angleMeasure, dihedralMeasure;	 

    // check if it is the normal 11 line format
    if( iss >> atomID >> atomType >> 
        oplsA >> oplsB >> 
        bondWith >> bondDistance >> 
        angleWith >> angleMeasure >> 
        dihedralWith >> dihedralMeasure >> extra)
    {
        format = 1;
    }
    else
    {
        someLine = line;
        if(someLine.find("TERZ")!=string::npos)
        {
            format = 2;
		}
        else if(someLine.find("Geometry Variations")!=string::npos)
        {
            format = 3;
        }
        else if(someLine.find("Variable Bonds")!=string::npos)
        {
            format = 4;
        }
        else if(someLine.find("Additional Bonds")!=string::npos)
        {
            format = 5;
        }
        else if(someLine.find("Harmonic Constraints")!=string::npos)
        {
            format = 6;
        }
        else if(someLine.find("Variable Bond Angles")!=string::npos)
        {
            format = 7;
        }
        else if(someLine.find("Additional Bond Angles")!=string::npos)
        {
            format = 8;
        }
        else if(someLine.find("Variable Dihedrals")!=string::npos)
        {
            format = 9;
        }
        else if(someLine.find("Additional Dihedrals")!=string::npos)
        {
            format = 10;
        }
        else if(someLine.find("Domain Definitions")!=string::npos)
        {
            format = 11;
        }
        else if(someLine.find("Final blank line")!=string::npos)
        {
            format = -2;  
        }
    }

    return format;
}

void Zmatrix_Scan::handleZAdditions(string line, int cmdFormat)
{
    vector<int> atomIds;
    int id;
    stringstream tss(line.substr(0,15) );

    if(line.find("AUTO")!=string::npos)
    {
	     //Do stuff for AUTO
    }
    else
    {
        while(tss >> id)
        {
            atomIds.push_back(id);
            if(tss.peek()=='-'||tss.peek()==','||tss.peek()==' ')
            {
                tss.ignore();
            }
        }

        int start, end=0;
        if( atomIds.size()== 1)
        {
            start = atomIds[0];
            end = atomIds[0]; 
        }
        else if(atomIds.size() == 2)
        {
            start = atomIds[0]; end = atomIds[1];
        }

        switch(cmdFormat)
        {
            case 3:
                // Geometry Variations follow 
            break;
            case 4:
                // Variable Bonds follow
                for(int i=0; i< moleculePattern[0].numOfBonds; i++)
                {
                    if(  moleculePattern[0].bonds[i].atom1 >= start &&  moleculePattern[0].bonds[i].atom1 <= end)
                    {
                        //cout << "Bond Atom1: "<<  moleculePattern[0].bonds[i].atom1 << " : " <<  moleculePattern[0].bonds[i].variable<<endl;//DEBUG
                        moleculePattern[0].bonds[i].variable = true;
                    }
                }
                break;
            case 5:
                //  Additional Bonds follow 
                break;
            case 6:
                // Harmonic Constraints follow 
                break;
            case 7:
                //  Variable Bond Angles follow
                for(int i=0; i<  moleculePattern[0].numOfAngles; i++)
                {
                    if(  moleculePattern[0].angles[i].atom1 >= start && moleculePattern[0].angles[i].atom1 <= end)
                    {
                    //cout << "Angle Atom1: "<<  moleculePattern[0].angles[i].atom1 << " : " << moleculePattern[0].angles[i].variable << endl;//DEBUG
                    moleculePattern[0].angles[i].variable = true;
                    }
                }
                break;
            case 8:
                // Additional Bond Angles follow
                break;
            case 9:
                // Variable Dihedrals follow
                for(int i=0; i< moleculePattern[0].numOfDihedrals; i++)
                {
                    if(  moleculePattern[0].dihedrals[i].atom1 >= start &&  moleculePattern[0].dihedrals[i].atom1 <= end )
                    {
                        //cout << "Dihedral Atom1: "<<  moleculePattern[0].dihedrals[i].atom1 << " : " <<   moleculePattern[0].dihedrals[i].variable << endl;//DEBUG
                        moleculePattern[0].dihedrals[i].variable = true;
                    }
                }
                break;
            case 10:
                //  Domain Definitions follow
                break;
            default:
                //Do nothing
                break;
        }
    }
}

vector<Hop> Zmatrix_Scan::calculateHops(Molecule molec)
{
    vector<Hop> newHops;
    int **graph;
    int size = molec.numOfAtoms;
	 int startId = molec.atoms[0].id;

    buildAdjacencyMatrix(graph,molec);

    for(int atom1=0; atom1<size; atom1++)
    {
        for(int atom2=atom1+1; atom2<size; atom2++)
        {
            int distance = findHopDistance(atom1,atom2,size,graph);
            if(distance >=3)
            {
				Hop tempHop = createHop(atom1+startId,atom2+startId,distance); //+startId because atoms may not start at 1
                newHops.push_back(tempHop);					
            }  		      
        }
    }

    return newHops; 
}

bool Zmatrix_Scan::contains(vector<int> &vect, int item)
{
    for(int i=0; i<vect.size(); i++)
    {
        if(vect[i]==item)
        {
            return true;
        }
    }

    return false;
}


int Zmatrix_Scan::findHopDistance(int atom1,int atom2,int size, int **graph)
{
    map<int,int> distance;
    queue<int> Queue;
    vector<int> checked;
    vector<int> bonds;


    Queue.push(atom1);
    checked.push_back(atom1);
    distance.insert( pair<int,int>(atom1,0) );	

    while(!Queue.empty())
    {
        int target = Queue.front();
        Queue.pop();
        if(target == atom2)
        {
            return distance[target];
        }

        bonds.clear();
        for(int col=0;col<size;col++)
        {
            if( graph[target][col]==1 )
            {
                bonds.push_back(col);
            }
        }

        for(int x=0; x<bonds.size();x++)
        {
            int currentBond = bonds[x];
            if(!contains(checked,currentBond) )
            {
                checked.push_back(currentBond);
                int newDistance = distance[target]+1;
                distance.insert(pair<int,int>(currentBond, newDistance));
                Queue.push(currentBond);
            }
        }
    }

    //Needs a return value
    return -1; //Temp fill
}

void Zmatrix_Scan::buildAdjacencyMatrix(int **&graph, Molecule molec)
{
    int size = molec.numOfAtoms;
	int startId = molec.atoms[0].id; //the first atom ID in the molecule
	int lastId = startId + molec.numOfAtoms -1; //the last atom ID in the molecule
    graph =  new int*[size]; //create colums
    for(int i=0; i<size; i++) //create rows
    {
        graph[i]=new int[size];	
    }

    //fill with zero
    for(int c=0; c<size; c++)
    {
        for(int r=0; r<size; r++)
        {
            graph[c][r]=0;
        }
    }

    //fill with adjacent array with bonds
    for(int x=0; x<molec.numOfBonds; x++)
    {
        Bond bond = molec.bonds[x];
		  //make sure the bond is intermolecular
		if( (bond.atom1 >= startId && bond.atom1 <= lastId) &&
		   (bond.atom2 >= startId && bond.atom2 <= lastId) )
        {
			graph[bond.atom1-startId][bond.atom2-startId]=1;
            graph[bond.atom2-startId][bond.atom1-startId]=1;
		}       
    }
}

vector<Molecule> Zmatrix_Scan::buildMolecule(int startingID)
{
	int numOfMolec = moleculePattern.size();
	Molecule newMolecules[numOfMolec];
	 
    //need a deep copy of molecule pattern incase it is modified.
    for (int i = 0; i < moleculePattern.size(); i++)
    {
        Atom *atomCopy = new Atom[ moleculePattern[i].numOfAtoms] ;
        for(int a=0; a <  moleculePattern[i].numOfAtoms ; a++)
        {
            atomCopy[a]=  moleculePattern[i].atoms[a];
        }

        Bond *bondCopy = new Bond[ moleculePattern[i].numOfBonds] ;
        for(int a=0; a <  moleculePattern[i].numOfBonds ; a++)
        {
            bondCopy[a]=  moleculePattern[i].bonds[a];
        }

        Angle *angleCopy = new Angle[ moleculePattern[i].numOfAngles] ;
        for(int a=0; a <  moleculePattern[i].numOfAngles ; a++)
        {
            angleCopy[a]=  moleculePattern[i].angles[a];
        }

        Dihedral *dihedCopy = new Dihedral[ moleculePattern[i].numOfDihedrals];
        for(int a=0; a <  moleculePattern[i].numOfDihedrals ; a++)
        {
            dihedCopy[a]=  moleculePattern[i].dihedrals[a];
        }

        //calculate and add array of Hops to the molecule
        vector<Hop> calculatedHops;
        calculatedHops = calculateHops(moleculePattern[i]);
        int numOfHops = calculatedHops.size();
        Hop *hopCopy = new Hop[numOfHops];
        for(int a=0; a < numOfHops; a++)
        {
            hopCopy[a] = calculatedHops[a];
        }


        Molecule molecCopy = createMolecule(-1,atomCopy, angleCopy, bondCopy, dihedCopy, hopCopy, 
                                    moleculePattern[i].numOfAtoms, 
                                    moleculePattern[i].numOfAngles,
                                    moleculePattern[i].numOfBonds,
                                    moleculePattern[i].numOfDihedrals,
                                    numOfHops);	
		  newMolecules[i] = molecCopy; 
             
    }
				
	//Assign/calculate the appropiate x,y,z positions to the molecules. 									
	//buildMoleculeInSpace(newMolecules, numOfMolec);
	buildMoleculeXYZ(newMolecules, numOfMolec);
	 

    for (int i = 0; i < numOfMolec; i++)
    {
        if(i == 0)
        {
            newMolecules[i].id = startingID;
        }
        else
        {
            newMolecules[i].id = newMolecules[i-1].id + newMolecules[i-1].numOfAtoms; 
        }
    }
	 
    for (int j = 0; j < numOfMolec; j++)
    {
        Molecule newMolecule = newMolecules[j];
        //map unique IDs to atoms within structs based on startingID
        for(int i = 0; i < newMolecules[j].numOfAtoms; i++)
        {
            int atomID = newMolecule.atoms[i].id - 1;
            //newMolecule.atoms[i].id = atomID + newMolecule.id;
				newMolecule.atoms[i].id = atomID + startingID;

        }
        for (int i = 0; i < newMolecule.numOfBonds; i++)
        {
            int atom1ID = newMolecule.bonds[i].atom1 - 1;
            int atom2ID = newMolecule.bonds[i].atom2 - 1;

            //newMolecule.bonds[i].atom1 = atom1ID + newMolecule.id;
            //newMolecule.bonds[i].atom2 = atom2ID + newMolecule.id;
				newMolecule.bonds[i].atom1 = atom1ID + startingID;
            newMolecule.bonds[i].atom2 = atom2ID + startingID;
        }
        for (int i = 0; i < newMolecule.numOfAngles; i++)
        {
            int atom1ID = newMolecule.angles[i].atom1 - 1;
            int atom2ID = newMolecule.angles[i].atom2 - 1;

            //newMolecule.angles[i].atom1 = atom1ID + newMolecule.id;
            //newMolecule.angles[i].atom2 = atom2ID + newMolecule.id;
				newMolecule.angles[i].atom1 = atom1ID + startingID;
            newMolecule.angles[i].atom2 = atom2ID + startingID;
        }
        for (int i = 0; i < newMolecule.numOfDihedrals; i++)
        {
            int atom1ID = newMolecule.dihedrals[i].atom1 - 1;
            int atom2ID = newMolecule.dihedrals[i].atom2 - 1;

            //newMolecule.dihedrals[i].atom1 = atom1ID + newMolecule.id;
            //newMolecule.dihedrals[i].atom2 = atom2ID + newMolecule.id;
				newMolecule.dihedrals[i].atom1 = atom1ID + startingID;
            newMolecule.dihedrals[i].atom2 = atom2ID + startingID;
        }
        for (int i = 0; i < newMolecule.numOfHops; i++)
        {
            int atom1ID = newMolecule.hops[i].atom1 - 1;
            int atom2ID = newMolecule.hops[i].atom2 - 1;

            //newMolecule.hops[i].atom1 = atom1ID + newMolecule.id;
            //newMolecule.hops[i].atom2 = atom2ID + newMolecule.id;
				newMolecule.hops[i].atom1 = atom1ID + startingID;
            newMolecule.hops[i].atom2 = atom2ID + startingID;
        }
    }

    return vector<Molecule>(newMolecules,newMolecules+sizeof(newMolecules)/sizeof(Molecule) );
}
*/


//_________________________________________________________________________________________________________________
//  Methods to handle logging of status & errors.
//_________________________________________________________________________________________________________________

/*
//this comes from metroUtil...
//this supports writing to the log with relative ease.
//These functions are namespace-less and class-less.
//said string contains almost all of the text you wish to be written to the file.
//opening and closing of the file will be done on the fly, and should be guaranteed once this method reaches its end.
//If stringstream is needed, you may call the overloaded version below, which will still call this version of the method
*@param: text: the text to be written to the output file
*@param: stamp: under which category we log this text: start of simulation [START], end of simulation [END],
*		error while handling the OPLS files [OPLS] Zmatrix files [Z_MATRIX] or geometric config files [GEO], or generic [DEFAULT].
*
*@returns: [none]
*/
void writeToLog(std::string text,int stamp){
    std::string filename = "OutputLog";
	std::ofstream logFile;
	logFile.open(filename.c_str(),std::ios::out|std::ios::app);
	 
	std::string hash ="";
	time_t current_time;
    struct tm * time_info;
    char timeString[9];  // space for "HH:MM:SS\0"
	 
	 switch(stamp){
	     case START:
		      //The start of a new simulation
		      logFile << "\n\n\n\n\n\n" << std::endl;
				logFile << "======================================================================"<< std::endl;
				logFile << "                       Starting Simulation: ";				
            time(&current_time);
            time_info = localtime(&current_time);
            strftime(timeString, sizeof(timeString), "%H:%M:%S", time_info);
				logFile << timeString << std::endl;
				logFile << "----------------------------------------------------------------------"<< std::endl;
				break;
			case END: 
			   //The end of a running simulation
				logFile << "----------------------------------------------------------------------"<< std::endl;
				logFile << "                       Ending Simulation: ";
            time(&current_time);
            time_info = localtime(&current_time);
            strftime(timeString, sizeof(timeString), "%H:%M:%S", time_info);
				logFile << timeString << std::endl;
				logFile << "======================================================================"<< std::endl;
				break;		
	     case OPLS:
		      //OPLS error
	         logFile << "--OPLS: ";
		      break;
	     case Z_MATRIX:
		      //Zmatrix error
	         logFile << "--Z_Matrix: ";
		      break;
		  case GEOM:
		      //GEOM error Geometric
				logFile << "--GEOM: ";
				break;
	     default:
	         logFile << "";
		      break;		
	 }
	 logFile << text << std::endl;
	 logFile.close();	 
}


//this method allows for writing to a given log file, with some measure of automation
//this overloaded version allows for a stringstream, instead of a normal string, to be input.
//said string contains almost all of the text you wish to be written to the file.
//opening and closing of the file will be done on the fly, and should be guaranteed once this method reaches its end.
void writeToLog(std::stringstream& ss, int stamp ){
    writeToLog(ss.str(),stamp);
	 ss.str(""); // clears the string steam...
	 ss.clear();
}