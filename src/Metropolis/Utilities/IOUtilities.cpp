/*Intended goal: support read, parse, and extract operations on configuration files to properly initialize 
*  a simulation environment.
* [1] This file/class will, ideally, replace Config_Scan.cpp & will augment MetroUtil.cpp. (Started, 19 Feb. Beta completion, 28 Feb.)
* [2] This file/class will, ideally, replace Opls_Scan and Zmatrix_Scan. (Started, 26 Feb. -Albert)
*Created 19 February 2014. Albert Wallace
*/
/*
--Changes made on:
		->Sun, 23 Feb 2014 (Albert)
		->Wed, 26 Feb (Albert)
		->Thu, 27 Feb (Albert, then Tavis)
		->Fri, 28 Feb (Albert)
		->Mon, 03 Mar; Wed, 05 Mar; Thur, 06 Mar; Fri, 07 Mar; Mon, 10 Mar; Wed, 12 Mar; Thur, 13 Mar (Albert)
*/
/*Based on work from earlier sessions by Alexander Luchs, Riley Spahn, Seth Wooten, and Orlando Acevedo*/

//Dear Self: Tavis says we might be better off importing an instance of a box and filling it with the necessary information.
// He has good ideas and we should listen to him.
// Also, do destructors. Please.

//Again, dear self: Also try to handle the Zmatrix and State file. If there is no state file but there is a Zmatrix, it's fine.
// if there is no Zmatrix but a state file, also fine. If both are there, that's fine. If neither is there, make a critical error out of it
// and halt the process.


//Search for generatefccbox. At this point, we give the box all the info it needs to create itself.



//_________________________________________________________________________________________________________________
//  INCLUDE statements
//_________________________________________________________________________________________________________________

#include "IOUtilities.h"
#include "StructLibrary.h"
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
IOUtilities::IOUtilities(std::string configPathIn){


	//UtilitiesInfo ; //all the variables used for this class are stuck in this struct for easy, yet unsafe, access
					//there were apparently getters and setters to access them all, so if necessary, we can have one getter later for the entire struct
	
	//note to people/myself: //enviro = currentEnvironment
	
	configPath = configPathIn; //the path to the primary configuration file, which holds all other potential file paths
	currentEnvironment = new Environment(); //The current working environment for the simulation
    unsigned int numOfSteps = 0; //The number of steps to run the simulation
	oplsuaparPath = ""; //The path to the opls files containing additional geometry data, to be used (eventually) during simulation
	zmatrixPath = ""; //The path to the Z-matrix files to be used during simulation
	statePath = ""; //The path to the state information file to be used in the simulation
	stateOutputPath = ""; //The path where we write the state output files after simulation
	pdbOutputPath = ""; //The path where we write the pdb output files after simulation
	cutoff = 0; //The nonbonded cutoff distance.
	
	criticalErrorEncountered = false; //set to true if, at any point, a critical error prevents the program from reliably setting up the environment
    
    pullInDataToConstructSimBox(); //do the rest of the construction with this driver method
}

/*
*Primary method called to read from the main configuration file. Parses all known information in the predefined format
*	and stores it into an Environment structure, or otherwise in the UtilitiesInfo structure which holds various file paths, etc.
*
*@params: [none; uses variables within the class to pass information]
*@return: [none; uses variables within the class to pass information]
*/
bool IOUtilities::readInConfig()
{
	bool isSafeToContinue = true;
	std::ifstream configscanner(configPath.c_str());
	if (! configscanner.is_open())
	{
		throwScanError("Configuration file failed to open.");
		return false;
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
						currentEnvironment->x = atof(line.c_str());
					}
					else
					{
						throwScanError("Configuration file not well formed. Missing environment x value.");
						return false;
					}
					break;
				case 3:
					if(line.length() > 0)
					{
						currentEnvironment->y = atof(line.c_str());
					}
					else
					{
						throwScanError("Configuration file not well formed. Missing environment y value.");
						return false;
					}
					break;
				case 4:
					if(line.length() > 0)
					{
						currentEnvironment->z = atof(line.c_str());
					}
					else
					{
						throwScanError("Configuration file not well formed. Missing environment z value.");
						return false;
					}
					break;
				case 6:
					if(line.length() > 0)
					{
						currentEnvironment->temp = atof(line.c_str());
					}
					else
					{
						throwScanError("Configuration file not well formed. Missing environment temperature value.");
						return false;
					}
					break;
				case 8:
					if(line.length() > 0)
					{
						currentEnvironment->maxTranslation = atof(line.c_str());
					}
					else
					{
						throwScanError("Configuration file not well formed. Missing environment max translation value.");
						return false;
					}
					break;
				case 10:
					if(line.length() > 0)
					{
						numOfSteps = atoi(line.c_str());
					}
					else
					{
						throwScanError("Configuration file not well formed. Missing number of steps value.");
						return false;
					}
					break;
				case 12:
					if(line.length() > 0)
					{
						currentEnvironment->numOfMolecules = atoi(line.c_str());
						//printf("number is %d",enviro.numOfMolecules);
					}
					else
					{
						throwScanError("Configuration file not well formed. Missing number of molecules value.");
						return false;
					}
					break;
				case 14:
					if(line.length() > 0)
					{
						oplsuaparPath = line;
					}
					else
					{
						throwScanError("Configuration file not well formed. Missing oplsuapar path value.");
						return false;
					}
					break;
				case 16:
					if(line.length() > 0)
					{
						zmatrixPath = line;
					}
					else
					{
						throwScanError("Configuration file not well formed. Missing z-matrix path value.");
						isSafeToContinue = false; //now that it is false, check if there is a state file
						exit(1);
					}
					break;
				case 18:
					if(line.length() > 0)
					{
						statePath = line;
					}
					else if (isSafeToContinue == false)
					{
						criticalErrorEncountered = true;
						throwScanError("Configuration file not well formed. Missing value pointing to prior state file path. Cannot safely continue with program execution.");
						return false;
					}
					else
					{	
						throwScanError("INFO: Value pointing to prior state file path not found in main config file. Environment setup defaulting to clean simulation.");
						isSafeToContinue = true;
					}
					break;
				case 20:
					if(line.length() > 0){
						stateOutputPath = line;
					}
					else
					{
						throwScanError("Configuration file not well formed. Missing state file output path value.");
						return false;
					}
					break;
				case 22:
					if(line.length() > 0){
						pdbOutputPath = line;
					}
					else
					{
						throwScanError("Configuration file not well formed. Missing PDB output path value.");
						return false;
					}
					break;
				case 24:
					if(line.length() > 0)
					{
						currentEnvironment->cutoff = atof(line.c_str());
					}
					else
					{
						throwScanError("Configuration file not well formed. Missing environment cutoff value.");
						return false;
					}
					break;
				case 26:
					if(line.length() > 0)
					{
						currentEnvironment->maxRotation = atof(line.c_str());
					}
					else
					{
						throwScanError("Configuration file not well formed. Missing environment max rotation value.");
						return false;
					}
					break;
				case 28:
					if(line.length() > 0)
					{
						currentEnvironment->randomseed=atoi(line.c_str());
					}
					else
					{
						throwScanError("Configuration file not well formed. Missing random seed value.");
						return false;
					}
					break;
				case 30:
					if(line.length() > 0)
					{
						// Convert to a zero-based index
						currentEnvironment->primaryAtomIndex=atoi(line.c_str()) - 1;
					}
					else
					{
						throwScanError("Configuration file not well formed. Missing environment primary atom index value.");
						return false;
					}
					break;
					//end of disabled configuration file code.
			}
		
			currentLine++;
		}
	}
	return true;
}


/*
	*Used to throw error messages to screen during the read of the primary config file, with a count of the errors encountered
	*@param: message: the message to be displayed on the screen
	*/
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
int IOUtilities::ReadStateFile(char const* StateFile, Environment * destinationEnvironment, Molecule * destinationMoleculeCollection)
{
    std::ifstream inFile;
    Environment tmpenv;
    std::stringstream ss;
    char buf[250];
    
    std::cout<<"read state file "<<StateFile<<std::endl;
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
       ss<<"Wrong state files,does not match other configfiles"<<std::endl;
       ss<<"x "<<tmpenv.x<<" "<<destinationEnvironment->x<<std::endl;
       ss<<"y "<<tmpenv.y<<" "<<destinationEnvironment->y<<std::endl;
       ss<<"z "<<tmpenv.z<<" "<<destinationEnvironment->z<<std::endl;
       ss<<"numOfAtoms "<<tmpenv.numOfAtoms<<" "<<destinationEnvironment->numOfAtoms<<std::endl;
       ss<<"temperature "<<tmpenv.temp<<" "<<destinationEnvironment->temp<<std::endl;
       ss<<"cutoff "<<tmpenv.cutoff<<" "<<destinationEnvironment->cutoff<<std::endl;
       ss<<ss.str()<<std::endl; 
       writeToLog(ss, DEFAULT);      
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
 	Molecule *moleculePtr=destinationMoleculeCollection;

    while(inFile.good()&&molecno<destinationEnvironment->numOfMolecules)
    {
        inFile>>no;
        assert(moleculePtr->moleculeIdentificationNumber==no);
        inFile.getline(buf,sizeof(buf)); //bypass atom flag
        inFile.getline(buf,sizeof(buf));
        assert(strcmp(buf,"= Atoms")==0);

        for(int i=0;i<moleculePtr->numOfAtoms;i++)
        {
        	inFile>>currentAtom.atomIdentificationNumber >> currentAtom.x >> currentAtom.y >> currentAtom.z>> currentAtom.sigma >> currentAtom.epsilon >> currentAtom.charge;
        	assert(currentAtom.atomIdentificationNumber==moleculePtr->atoms[i].atomIdentificationNumber);
        	//printf("id:%d,x:%f,y:%f\n",currentAtom.atomIdentificationNumber,currentAtom.x,currentAtom.y);
        	memcpy(&moleculePtr->atoms[i],&currentAtom,sizeof(Atom));
        }

        inFile.getline(buf,sizeof(buf)); //ignore bonds flag
        inFile.getline(buf,sizeof(buf));
        assert(strcmp(buf,"= Bonds")==0);
        for(int i=0;i<moleculePtr->numOfBonds;i++)
        {
        	inFile>>currentBond.atom1 >>currentBond.atom2 >> currentBond.lengthOfBond >> currentBond.canBeVaried;
        	assert(currentBond.atom1==moleculePtr->bonds[i].atom1);
        	assert(currentBond.atom2==moleculePtr->bonds[i].atom2);      	
        	memcpy(&moleculePtr->bonds[i],&currentBond,sizeof(Bond));
        }

        inFile.getline(buf,sizeof(buf)); //ignore Dihedrals flag
        inFile.getline(buf,sizeof(buf));
        assert(strcmp(buf,"= Dihedrals")==0);
        for(int i=0;i<moleculePtr->numOfDihedrals;i++)
        {
        	inFile>>currentDi.atom1>>currentDi.atom2>>currentDi.distanceBetweenAtoms>>currentDi.canBeVaried;
        	assert(currentDi.atom1==moleculePtr->dihedrals[i].atom1);
        	assert(currentDi.atom2==moleculePtr->dihedrals[i].atom2);      	
        	memcpy(&moleculePtr->dihedrals[i],&currentDi,sizeof(Dihedral));
        }

        inFile.getline(buf,sizeof(buf)); //ignore hops flag
        inFile.getline(buf,sizeof(buf));
        assert(strcmp(buf,"=Hops")==0);
        // known BUG - if molecule has no hops (3 atoms or less) state file gives error crashing simulation
        for(int i=0;i<moleculePtr->numOfHops;i++)
        {
        	inFile>>currentHop.atom1>>currentHop.atom2 >>currentHop.hop;
        	assert(currentHop.atom1==moleculePtr->hops[i].atom1);
        	assert(currentHop.atom2==moleculePtr->hops[i].atom2);      	
        	memcpy(&moleculePtr->hops[i],&currentHop,sizeof(Hop));
        }

        inFile.getline(buf,sizeof(buf)); //ignore angles flag
        inFile.getline(buf,sizeof(buf));
        assert(strcmp(buf,"= Angles")==0);
        for(int i=0;i<moleculePtr->numOfAngles;i++)
        {
        	inFile>>currentAngle.atom1 >> currentAngle.atom2 >>currentAngle.magnitudeOfAngle >>currentAngle.canBeVaried;
        	assert(currentAngle.atom1==moleculePtr->angles[i].atom1);
        	assert(currentAngle.atom2==moleculePtr->angles[i].atom2);      	
        	memcpy(&moleculePtr->angles[i],&currentAngle,sizeof(Angle));
        }       

        inFile.getline(buf,sizeof(buf)); //bypass == flag
        inFile.getline(buf,sizeof(buf));
        assert(strcmp(buf,"==")==0);   

        moleculePtr++;                    
        molecno++;
    }
    inFile.close();
    WriteStateFile("Confirm.state", destinationEnvironment, moleculePtr);

	return 0;
}




/**
	Used to write to a state file.
	//copied over from...Simbox. Retooling needed. (Variable references point to SimBox locations)
	
	@param StateFile - writes to a state file at the location given
*/

int IOUtilities::WriteStateFile(char const* StateFile, Environment * sourceEnvironment, Molecule * sourceMoleculeCollection)
{
    std::ofstream outFile;
    int numOfMolecules=sourceEnvironment->numOfMolecules;
    
    outFile.open(StateFile);
    
    //print the environment
    outFile << sourceEnvironment->x << " " << sourceEnvironment->y << " " << sourceEnvironment->z << " " << sourceEnvironment->maxTranslation<<" " << sourceEnvironment->numOfAtoms
        << " " << sourceEnvironment->temp << " " << sourceEnvironment->cutoff <<std::endl;
    outFile << std::endl; // blank line
    
    for(int i = 0; i < numOfMolecules; i++)
    {
        Molecule currentMol = sourceMoleculeCollection[i];
        outFile << currentMol.moleculeIdentificationNumber << std::endl;
        outFile << "= Atoms" << std::endl; // delimiter
    
        //write atoms
        for(int j = 0; j < currentMol.numOfAtoms; j++)
        {
            Atom currentAtom = currentMol.atoms[j];
            outFile << currentAtom.atomIdentificationNumber << " "
                << currentAtom.x << " " << currentAtom.y << " " << currentAtom.z
                << " " << currentAtom.sigma << " " << currentAtom.epsilon  << " "
                << currentAtom.charge << std::endl;
        }
        outFile << "= Bonds" << std::endl; // delimiter
        
        //write bonds
        for(int j = 0; j < currentMol.numOfBonds; j++)
        {
            Bond currentBond = currentMol.bonds[j];
            outFile << currentBond.atom1 << " " << currentBond.atom2 << " "
                << currentBond.lengthOfBond << " ";
            if(currentBond.canBeVaried)
                outFile << "1" << std::endl;
            else
                outFile << "0" << std::endl;

        }
        outFile << "= Dihedrals" << std::endl; // delimiter
        for(int j = 0; j < currentMol.numOfDihedrals; j++)
        {
            Dihedral currentDi = currentMol.dihedrals[j];
            outFile << currentDi.atom1 << " " << currentDi.atom2 << " "
                << currentDi.distanceBetweenAtoms << " ";
            if(currentDi.canBeVaried)
            {
                outFile << "1" << std::endl;
            }
            else
            {
                outFile << "0" << std::endl;
            }
        }

        outFile << "=Hops" << std::endl;

        for(int j = 0; j < currentMol.numOfHops; j++)
        {
            Hop currentHop = currentMol.hops[j];

            outFile << currentHop.atom1 << " " << currentHop.atom2 << " "
                << currentHop.hop << std::endl;
        }
        
        
        outFile << "= Angles" << std::endl; // delimiter

        //print angless
        for(int j = 0; j < currentMol.numOfAngles; j++)
        {
            Angle currentAngle = currentMol.angles[j];

            outFile << currentAngle.atom1 << " " << currentAngle.atom2 << " "
                << currentAngle.magnitudeOfAngle << " ";
            if(currentAngle.canBeVaried)
            {
                outFile << "1" << std::endl;
            }
            else
            {
                outFile << "0" << std::endl;
            }
        }


        //write a == line
        outFile << "==" << std::endl;
    }
    outFile.close();
	return 0;
}





/**
	*Writes to a PDB file for visualizing the box
	//this method copied over from...Simbox. Retooling needed. (Variable references point to SimBox locations)
	
	//other versions found in metroUtil were not used; it relied only on SimBox
	
	//to do this, the main simulation controller--linearSimulationExample in the old data structure--will have to manage the data in the
	// PDB, by directly modifying it in the instance of the IOUtilities class we will have to make,
	// rather than SimBox managing it.
	
	// the downside is that this implementation, when done for parallel,
	//  will likely require writing back to the host, rather than storing the info on the GPU until the very end...
	//  as the SimBox holds all the current information being accessed by this method.
	
	*@param pdbFile - Location of the pdbFile to be written to
	*@param sourceEnvironment: The environment as it currently exists in the SimBox to be visualized
	*@param: sourceMoleculeCollection: Is an array of molecules, dynamically allocated & created elsewhere [such as in the SimBox]
*/

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
        for (int j = 0; j < currentMol.numOfAtoms; j++)
        {
        	Atom currentAtom = currentMol.atoms[j];
            outputFile.setf(std::ios_base::left,std::ios_base::adjustfield);
            outputFile.width(6);
            outputFile << "ATOM";
            outputFile.setf(std::ios_base::right,std::ios_base::adjustfield);
            outputFile.width(5);
            outputFile << currentAtom.atomIdentificationNumber + 1;
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
void IOUtilities::destructOpls_Scan()
{
    oplsTable.clear();
}

int IOUtilities::scanInOpls()
{
    int numOfLines=0;
    std::ifstream oplsScanner(oplsuaparPath.c_str()); //##
    if( !oplsScanner.is_open() ) //##
        return -1; //##
    else { //##
        string line;  //##
        while( oplsScanner.good() ) //##
        { //##
            numOfLines++; //##
            getline(oplsScanner,line); //##

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
    string hashNum;
    int secCol;
    double charge,sigma,epsilon;
    string name, extra;
    std::stringstream ss(line);

    //check to see what format it is opls, V value, or neither
    int format = OPLScheckFormat(line);
              
    if(format == 1)
    {      	
        ss >> hashNum >> secCol >> name >> charge >> sigma >> epsilon;
        char *atomtype = (char*)name.c_str(); 
         
        Atom temp = Atom(0, -1, -1, -1, sigma, epsilon, charge, *atomtype);
        std::pair<std::map<std::string,Atom>::iterator,bool> ret;
        ret = oplsTable.insert(std::pair<std::string,Atom>(hashNum,temp) );

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

int IOUtilities::OPLScheckFormat(string line)
{   	 
    int hashNum, secCol;
    double charge,sigma,epsilon;
    string name, extra;
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
		  output<<"Errors found in the OPLS file: "<< oplsuaparPath<<std::endl;
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
				output << std::endl<< std::endl;
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
				output << std::endl<< std::endl;
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
				output << std::endl<< std::endl;
		  }
		  writeToLog(output,OPLS);
	}
}

Atom IOUtilities::getAtom(string hashNum)
{
    if(oplsTable.count(hashNum)>0 )
    {
        return oplsTable[hashNum];
	}
	else
    {
	    std::cerr << "Index does not exist: "<< hashNum <<std::endl;
		return Atom(0, -1, -1, -1, -1, -1, -1, '0');
	}
}

double IOUtilities::getSigma(string hashNum)
{
    if(oplsTable.count(hashNum)>0 )
    {
        Atom temp = oplsTable[hashNum];
        return temp.sigma;
    }
    else
    {
        std::cerr << "Index does not exist: "<< hashNum <<std::endl;
        return -1;
    }
}

double IOUtilities::getEpsilon(string hashNum)
{
    if(oplsTable.count(hashNum)>0 )
    {
        Atom temp = oplsTable[hashNum];
        return temp.epsilon;
    }
    else
    {
        std::cerr << "Index does not exist: "<< hashNum <<std::endl;
        return -1;
    }
}

double IOUtilities::getCharge(string hashNum)
{
    if(oplsTable.count(hashNum)>0 )
    {
        Atom temp = oplsTable[hashNum];
        return temp.charge;
    }
    else
    {
        std::cerr << "Index does not exist: "<< hashNum <<std::endl;
        return -1;
    }
}

Fourier IOUtilities::getFourier(string hashNum)
{
    if(fourierTable.count(hashNum)>0 )
    {
        Fourier temp = fourierTable[hashNum];
        return temp;
    }
    else
    {	    
        std::cerr << "Index does not exist: "<< hashNum <<std::endl;
        Fourier temp ={-1,-1,-1,-1};
        return temp;
    }
}





//_________________________________________________________________________________________________________________
//  Zmatrix opening & parsing/Zmatrix configuration. [from Zmatrix_Scan.cpp]
//_________________________________________________________________________________________________________________


// Zmatrix_Scan::Zmatrix_Scan(string filename, Opls_Scan* oplsScannerRef)
// {
//     fileName = filename;
//     oplsScanner = oplsScannerRef;
//     startNewMolecule_ZM = false;
// }

void IOUtilities::destructZmatrix_Scan(){}

int IOUtilities::scanInZmatrix()
{
    std::stringstream output;
    int numOfLines=0;
    std::ifstream zmatrixScanner(zmatrixPath.c_str());

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

            if (startNewMolecule_ZM)
            {
                Atom* atomArray;
                Bond* bondArray;
                Angle* angleArray;
                Dihedral* dihedralArray;
                
                atomArray = (Atom*) malloc(sizeof(Atom) * atomVector_ZM.size());
                bondArray = (Bond*) malloc(sizeof(Bond) * bondVector_ZM.size());
                angleArray = (Angle*) malloc(sizeof(Angle) * angleVector_ZM.size());
                dihedralArray = (Dihedral*) malloc(sizeof(Dihedral) * dihedralVector_ZM.size());

                for (int i = 0; i < atomVector_ZM.size(); i++)
                {
                    atomArray[i] = atomVector_ZM[i];
                }
                for (int i = 0; i < bondVector_ZM.size(); i++)
                {
                    bondArray[i] = bondVector_ZM[i];
                }
                for (int i = 0; i < angleVector_ZM.size(); i++)
                {
                    angleArray[i] = angleVector_ZM[i];
                }
                for (int i = 0; i < dihedralVector_ZM.size(); i++)
                {
                    dihedralArray[i] = dihedralVector_ZM[i];
                }
				
				Hop * dummyHop = new Hop();
				
				Molecule dummyMolecule(-1, atomArray, angleArray, bondArray, dihedralArray, dummyHop,
                     atomVector_ZM.size(), angleVector_ZM.size(), bondVector_ZM.size(), dihedralVector_ZM.size(), 0);
				
                moleculePattern_ZM.push_back(dummyMolecule);

                atomVector_ZM.clear();
                bondVector_ZM.clear();
                angleVector_ZM.clear();
                dihedralVector_ZM.clear();

                startNewMolecule_ZM = false;
            } 
        }

        zmatrixScanner.close();
    	 
    }

    return 0;
}

void IOUtilities::parseLine(string line, int numOfLines)
{

    string atomID, atomType, oplsA, oplsB, bondWith, bondDistance, angleWith, angleMeasure, dihedralWith, dihedralMeasure;

    std::stringstream ss;

    //check if line contains correct format
    int format = ZMcheckFormat(line);

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
            lineAtom = getAtom(oplsA);
            lineAtom.atomIdentificationNumber = atoi(atomID.c_str());
            lineAtom.x = 0;
            lineAtom.y = 0;
            lineAtom.z = 0;
        }
        else//dummy atom
        {
        	char dummy = 'X';
            lineAtom = Atom(atoi(atomID.c_str()), -1, -1, -1, -1, -1, -1, dummy);
        }
		  atomVector_ZM.push_back(lineAtom);

        if (bondWith.compare("0") != 0)
        {
            lineBond.atom1 = lineAtom.atomIdentificationNumber;
            lineBond.atom2 = atoi(bondWith.c_str());
            lineBond.lengthOfBond = atof(bondDistance.c_str());
            lineBond.canBeVaried = false;
            bondVector_ZM.push_back(lineBond);
        }

        if (angleWith.compare("0") != 0)
        {
            lineAngle.atom1 = lineAtom.atomIdentificationNumber;
            lineAngle.atom2 = atoi(angleWith.c_str());
            lineAngle.magnitudeOfAngle = atof(angleMeasure.c_str());
            lineAngle.canBeVaried = false;
            angleVector_ZM.push_back(lineAngle);
        }

        if (dihedralWith.compare("0") != 0)
        {
            lineDihedral.atom1 = lineAtom.atomIdentificationNumber;
            lineDihedral.atom2 = atoi(dihedralWith.c_str());
            lineDihedral.distanceBetweenAtoms = atof(dihedralMeasure.c_str());
            lineDihedral.canBeVaried = false;
            dihedralVector_ZM.push_back(lineDihedral);
        }
    } //end if format == 1

    else if(format == 2)
    {
        startNewMolecule_ZM = true;
	}
    else if(format == 3)
    {
        startNewMolecule_ZM = true;
    }
    if (previousFormat_ZM >= 3 && format == -1)
    {
        handleZAdditions(line, previousFormat_ZM);
    }

    previousFormat_ZM = format;
}
    
int IOUtilities::ZMcheckFormat(string line)
{
    int format =-1; 
    std::stringstream iss(line);
    std::stringstream iss2(line);
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

void IOUtilities::handleZAdditions(string line, int cmdFormat)
{
    vector<int> atomIds;
    int id;
    std::stringstream tss(line.substr(0,15) );

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
                for(int i=0; i< moleculePattern_ZM[0].numOfBonds; i++)
                {
                    if(  moleculePattern_ZM[0].bonds[i].atom1 >= start &&  moleculePattern_ZM[0].bonds[i].atom1 <= end)
                    {
                        //std::cout << "Bond Atom1: "<<  moleculePattern_ZM[0].bonds[i].atom1 << " : " <<  moleculePattern_ZM[0].bonds[i].variable<<std::endl;//DEBUG
                        moleculePattern_ZM[0].bonds[i].canBeVaried = true;
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
                for(int i=0; i<  moleculePattern_ZM[0].numOfAngles; i++)
                {
                    if(  moleculePattern_ZM[0].angles[i].atom1 >= start && moleculePattern_ZM[0].angles[i].atom1 <= end)
                    {
                    //std::cout << "Angle Atom1: "<<  moleculePattern_ZM[0].angles[i].atom1 << " : " << moleculePattern_ZM[0].angles[i].variable << std::endl;//DEBUG
                    moleculePattern_ZM[0].angles[i].canBeVaried = true;
                    }
                }
                break;
            case 8:
                // Additional Bond Angles follow
                break;
            case 9:
                // Variable Dihedrals follow
                for(int i=0; i< moleculePattern_ZM[0].numOfDihedrals; i++)
                {
                    if(  moleculePattern_ZM[0].dihedrals[i].atom1 >= start &&  moleculePattern_ZM[0].dihedrals[i].atom1 <= end )
                    {
                        //std::cout << "Dihedral Atom1: "<<  moleculePattern_ZM[0].dihedrals[i].atom1 << " : " <<   moleculePattern_ZM[0].dihedrals[i].variable << std::endl;//DEBUG
                        moleculePattern_ZM[0].dihedrals[i].canBeVaried = true;
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

vector<Hop> IOUtilities::calculateHops(Molecule molec)
{
    vector<Hop> newHops;
    int **graph;
    int size = molec.numOfAtoms;
	 int startId = molec.atoms[0].atomIdentificationNumber;

    buildAdjacencyMatrix(graph,molec);

    for(int atom1=0; atom1<size; atom1++)
    {
        for(int atom2=atom1+1; atom2<size; atom2++)
        {
            int distance = findHopDistance(atom1,atom2,size,graph);
            if(distance >=3)
            {
				Hop tempHop(atom1+startId,atom2+startId,distance); //+startId because atoms may not start at 1
                newHops.push_back(tempHop);					
            }  		      
        }
    }

    return newHops; 
}

bool IOUtilities::contains(vector<int> &vect, int item)
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


int IOUtilities::findHopDistance(int atom1,int atom2,int size, int **graph)
{
    map<int,int> distance;
    std::queue<int> Queue;
    vector<int> checked;
    vector<int> bonds;


    Queue.push(atom1);
    checked.push_back(atom1);
    distance.insert( std::pair<int,int>(atom1,0) );	

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
                distance.insert(std::pair<int,int>(currentBond, newDistance));
                Queue.push(currentBond);
            }
        }
    }

    //Needs a return value
    return -1; //Temp fill
}

void IOUtilities::buildAdjacencyMatrix(int **&graph, Molecule molec)
{
    int size = molec.numOfAtoms;
	int startId = molec.atoms[0].atomIdentificationNumber; //the first atom ID in the molecule
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


vector<Molecule> IOUtilities::buildMolecule(int startingID)
{
	int numOfMolec = moleculePattern_ZM.size();
	Molecule newMolecules[numOfMolec];
	 
    //need a deep copy of molecule pattern incase it is modified.
    for (int i = 0; i < moleculePattern_ZM.size(); i++)
    {
        Atom *atomCopy = new Atom[ moleculePattern_ZM[i].numOfAtoms] ;
        for(int a=0; a <  moleculePattern_ZM[i].numOfAtoms ; a++)
        {
            atomCopy[a]=  moleculePattern_ZM[i].atoms[a];
        }

        Bond *bondCopy = new Bond[ moleculePattern_ZM[i].numOfBonds] ;
        for(int a=0; a <  moleculePattern_ZM[i].numOfBonds ; a++)
        {
            bondCopy[a]=  moleculePattern_ZM[i].bonds[a];
        }

        Angle *angleCopy = new Angle[ moleculePattern_ZM[i].numOfAngles] ;
        for(int a=0; a <  moleculePattern_ZM[i].numOfAngles ; a++)
        {
            angleCopy[a]=  moleculePattern_ZM[i].angles[a];
        }

        Dihedral *dihedCopy = new Dihedral[ moleculePattern_ZM[i].numOfDihedrals];
        for(int a=0; a <  moleculePattern_ZM[i].numOfDihedrals ; a++)
        {
            dihedCopy[a]=  moleculePattern_ZM[i].dihedrals[a];
        }

        //calculate and add array of Hops to the molecule
        vector<Hop> calculatedHops;
        calculatedHops = calculateHops(moleculePattern_ZM[i]);
        int numOfHops = calculatedHops.size();
        Hop *hopCopy = new Hop[numOfHops];
        for(int a=0; a < numOfHops; a++)
        {
            hopCopy[a] = calculatedHops[a];
        }


        
        	/*Molecule molecCopy(-1,atomCopy, angleCopy, bondCopy, dihedCopy, hopCopy, 
                                    moleculePattern_ZM[i].numOfAtoms, 
                                    moleculePattern_ZM[i].numOfAngles,
                                    moleculePattern_ZM[i].numOfBonds,
                                    moleculePattern_ZM[i].numOfDihedrals,
                                    numOfHops);	
		  newMolecules[i] = molecCopy; */
		 
		  newMolecules[i] = Molecule(-1,atomCopy, angleCopy, bondCopy, dihedCopy, hopCopy, 
                                    moleculePattern_ZM[i].numOfAtoms, 
                                    moleculePattern_ZM[i].numOfAngles,
                                    moleculePattern_ZM[i].numOfBonds,
                                    moleculePattern_ZM[i].numOfDihedrals,
                                    numOfHops);	
                                   
             
    }
				
	//Assign/calculate the appropiate x,y,z positions to the molecules. 									
	//buildMoleculeInSpace(newMolecules, numOfMolec);
	buildMoleculeXYZ(newMolecules, numOfMolec);
	 

    for (int i = 0; i < numOfMolec; i++)
    {
        if(i == 0)
        {
            newMolecules[i].moleculeIdentificationNumber = startingID;
        }
        else
        {
            newMolecules[i].moleculeIdentificationNumber = newMolecules[i-1].moleculeIdentificationNumber + newMolecules[i-1].numOfAtoms; 
        }
    }
	 
    for (int j = 0; j < numOfMolec; j++)
    {
        Molecule newMolecule = newMolecules[j];
        //map unique IDs to atoms within structs based on startingID
        for(int i = 0; i < newMolecules[j].numOfAtoms; i++)
        {
            int atomID = newMolecule.atoms[i].atomIdentificationNumber - 1;
            //newMolecule.atoms[i].id = atomID + newMolecule.id;
				newMolecule.atoms[i].atomIdentificationNumber = atomID + startingID;

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
    std::string logFilename = "OutputLog";
	std::ofstream logFile;
	logFile.open(logFilename.c_str(),std::ios::out|std::ios::app);
	 
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

/*________________________________________________________________________________________________________________________
__________________________________________________________________________________________________________________________
Environment/Box setup driver; ripped from old SimBox
__________________________________________________________________________________________________________________________
__________________________________________________________________________________________________________________________
*/
/**
	Creates the Box for use in the simulation.
	Builds the environment based on either the configuration file or the state file
	
	@param configScan - instance of configScan based on the utility
*/
void IOUtilities::pullInDataToConstructSimBox()
{
	molecules=NULL;
	enviro=NULL;
	std::stringstream ss;
	memset(&changedmole,0,sizeof(changedmole));
		
  ss << "Running simulation based on Z-Matrix File"<<std::endl;
  std::cout<<ss.str()<<std::endl; 
  writeToLog(ss, DEFAULT);

  //get environment from the config file
  
  //enviro=(Environment *)malloc(sizeof(Environment)); //TBRe
  //memcpy(enviro,(configScan.getEnviro()),sizeof(Environment)); //TBRe
  // //Albert note: enviro = currentEnvironment, FYI; TBR
  readInConfig();
	ss << "Reading Configuation File \nPath: " << configPath << std::endl;
  std::cout<<ss.str()<<std::endl; 
  writeToLog(ss, DEFAULT);

  //set up Opls scan 
  ss << "Reading OPLS File \nPath: " << oplsuaparPath << std::endl;
  std::cout<<ss.str()<<std::endl; 
  writeToLog(ss, DEFAULT);

  scanInOpls();
	ss << "OplsScan and OPLS ref table Created " << std::endl;
  std::cout<<ss.str()<<std::endl; 
  writeToLog(ss, DEFAULT);

  //set up zMatrixScan
  ss << "Reading Z-Matrix File \nPath: " << zmatrixPath << std::endl;
  std::cout<<ss.str()<<std::endl; 
  writeToLog(ss, DEFAULT);
  //Zmatrix_Scan zMatrixScan (configScan.getZmatrixPath(), &oplsScan); //TBRe
  if (scanInZmatrix() == -1)
  {
  	ss << "Error, Could not open: " << zmatrixPath << std::endl;
  	std::cerr << ss.str()<< std::endl;
  	writeToLog(ss, DEFAULT);
    exit(1);
   }
   
   ss << "Opened Z-Matrix File \nBuilding "<< currentEnvironment->numOfMolecules << " Molecules..." << std::endl;
   std::cout<<ss.str()<<std::endl; 
   writeToLog(ss, DEFAULT);

   //Convert molecule vectors into an array
   int moleculeIndex = 0;
   int atomCount = 0;

   vector<Molecule> molecVec = buildMolecule(atomCount);
   //###FLOATING POINT ERROR NEXT LINE!!!###
   int molecMod = currentEnvironment->numOfMolecules % molecVec.size();
   if (molecMod != 0)
   {
       currentEnvironment->numOfMolecules += molecVec.size() - molecMod;
       std::cout << "Number of molecules not divisible by specified z-matrix. Changing number of molecules to: " << currentEnvironment->numOfMolecules << std::endl;
    }
    molecules = (Molecule *)malloc(sizeof(Molecule) * currentEnvironment->numOfMolecules);
    
    int molecDiv = currentEnvironment->numOfMolecules / molecVec.size();
    //int molecTypenum;
    int molecTypenum=molecVec.size();
    
    int count[5];//sum up number of atoms,bonds,angles,dihedrals,hops
	//int * 
	//molecTables = new int
	//molecTables[molecVec.size()];
	Table * tables = new Table[molecVec.size()];
    memset(count,0,sizeof(count));
	int currentAtomCount = 0;
	
    for(int j = 0; j < molecVec.size(); j++)
    {
  		Molecule molec1 = molecVec[j];   
         //Copy data from vector to molecule
  		count[0]+=molec1.numOfAtoms;
  		count[1]+=molec1.numOfBonds;
  		count[2]+=molec1.numOfAngles;
  		count[3]+=molec1.numOfDihedrals;
  		count[4]+=molec1.numOfHops;
  		
  		std::cout << "before table building. Number of atom "<< molec1.numOfAtoms << std::endl;
  		
  		Hop *myHop = molec1.hops;
  		int **table;
  		table = new int*[molec1.numOfAtoms];
  		for(int k = 0; k< molec1.numOfAtoms;k++)
  			table[k] = new int[molec1.numOfAtoms];
  		//int table[molec1.numOfAtoms][molec1.numOfAtoms];
  		for(int test = 0; test< molec1.numOfAtoms;test++)
        {
  			for(int test1 = 0; test1 < molec1.numOfAtoms; test1++)
            {
  				table[test][test1] = 0;
  			}
		}
		
		for(int k2 = 0; k2<molec1.numOfHops;k2++)
        {
			int atom1 = myHop->atom1;
			std::cout << "atom1: "<< atom1-currentAtomCount << std::endl;
			int atom2 = myHop->atom2;
			std::cout << "atom2: "<< atom2-currentAtomCount<< std::endl;
			std::cout << "hop: " << myHop->hop << std::endl;
			table[atom1-currentAtomCount][atom2-currentAtomCount] = myHop->hop;
			table[atom2-currentAtomCount][atom1-currentAtomCount] = myHop->hop;
			myHop++;
		}
	  
	   std::cout << "after table building"<< std::endl;
	   //memset(table,0,sizeof(table));
	   //int table[molec1.numOfAtoms][molec1.numOfAtoms];
	   //std::cout << " this is " << j << std::endl;
	   tables[j] = Table(table); //createTable is in metroUtil
	   currentAtomCount += molec1.numOfAtoms;
	   std::cout << "after table creation. Current atom cout: "<< currentAtomCount << std::endl;
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

    //arrange first part of molecules
    memset(count,0,sizeof(count));
 	for(int j = 0; j < molecVec.size(); j++)
    {
 	      //Copy data from vector to molecule
        Molecule molec1 = molecVec[j];   

        molecules[j].atoms = (Atom *)(atompool+count[0]);
        molecules[j].bonds = (Bond *)(bondpool+count[1]);
        molecules[j].angles = (Angle *)(anglepool+count[2]);
        molecules[j].dihedrals = (Dihedral *)(dihedralpool+count[3]);
        molecules[j].hops = (Hop *)(hoppool+count[4]);

        molecules[j].moleculeIdentificationNumber = molec1.moleculeIdentificationNumber;
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
        for(int k = 0; k < molec1.numOfAtoms; k++)
        {
            molecules[j].atoms[k] = molec1.atoms[k];
        }               
           
        //assign bonds
        for(int k = 0; k < molec1.numOfBonds; k++)
        {
            molecules[j].bonds[k] = molec1.bonds[k];
        }

        //assign angles
        for(int k = 0; k < molec1.numOfAngles; k++)
        {
            molecules[j].angles[k] = molec1.angles[k];
        }

        //assign dihedrals
        for(int k = 0; k < molec1.numOfDihedrals; k++)
        {
            molecules[j].dihedrals[k] = molec1.dihedrals[k];
        }

        //assign hops zx add
        for(int k = 0; k < molec1.numOfHops; k++)
        {
            molecules[j].hops[k] = molec1.hops[k];
        }
    }
   
    for(int m = 1; m < molecDiv; m++)
    {
        int offset=m*molecTypenum;
    	memcpy(&molecules[offset],molecules,sizeof(Molecule)*molecTypenum);
    	for(int n=0;n<molecTypenum;n++)
        {
    		molecules[offset+n].moleculeIdentificationNumber=offset+n;
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
        {
            atompool[offset*count[0]+k].atomIdentificationNumber=offset*count[0]+k;
        }
        
        for(int k=0;k<count[1];k++)
        {
            bondpool[offset*count[1]+k].atom1+=m*count[0];
            bondpool[offset*count[1]+k].atom2+=m*count[0];
        }
        
        for(int k=0;k<count[2];k++)
        {
            anglepool[offset*count[2]+k].atom1+=m*count[0];
            anglepool[offset*count[2]+k].atom2+=m*count[0];
        }
        
        for(int k=0;k<count[3];k++)
        {
            dihedralpool[offset*count[3]+k].atom1+=m*count[0];
            dihedralpool[offset*count[3]+k].atom2+=m*count[0];
        }
        
        for(int k=0;k<count[4];k++)
        {
            hoppool[offset*count[4]+k].atom1+=m*count[0];
            hoppool[offset*count[4]+k].atom2+=m*count[0];
        }
    }
    std::cout << "ERROR AFTER THIS" << std::endl;
    //currentEnvironment->numOfAtoms = count[0]*molecDiv;
	//ss << "Molecules Created into an Array" << std::endl;
    //writeToLog(ss, DEFAULT);
     
    if (!statePath.empty())
    {
		//ss << "Reading State File \nPath: " << statePath << std::endl;
        //std::cout << ss.str() << std::endl; 
        //writeToLog(ss, DEFAULT);
     	//ReadStateFile(configScan.getStatePath().c_str());
     	//ReadStateFile(statePath.c_str(), currentEnvironment, molecules);
    }
    else
    {
        //ss << "Assigning Molecule Positions..." << std::endl;
        //writeToLog(ss, DEFAULT);
        //generatefccBox(molecules,enviro);//generate fcc lattice box
        //ss << "Finished Assigning Molecule Positions" << std::endl;
        //writeToLog(ss, DEFAULT);
    }
}