/*Intended goal: support read, parse, and extract operations on configuration files to properly initialize 
*  a simulation environment.
* //This file/class will, ideally, replace Config_Scan.cpp & will augment MetroUtil.cpp
*Created 19 February 2014. N. Coleman, A. Wallace
*/
//Changes on:
//	Sun, 23 Feb 2014. 1530PM to 1558PM, 1611 to 1655PM, 1757PM to 2031PM


/*!\file
  \brief Class used to read and write configuration files.
  \author Alexander Luchs, Riley Spahn, Seth Wooten, and Orlando Acevedo
 
 */
 
 
#include <assert.h>
#include "IOUtilities.cuh"
#include "errno.h"

IOUtilities::IOUtilities(string configPath){


	//UtilitiesInfo filePathsEtc; //all the variables used for this class are stuck in this struct for easy, yet unsafe, access
					//there were apparently getters and setters to access them all, so if necessary, we can have one getter later for the entire struct
	
	//note to people/myself: //enviro = currentEnvironment
	
    filePathsEtc.configPath = configPath;
    //memset(&enviro,0,sizeof(Environment)); //the Environment is a struct, and this is apparently the best way to instantiate the struct
    			//except that may not be required?
    filePathsEtc.numOfSteps=0;
    
    ifstream configscanner(configPath.c_str());
    if (! configscanner.is_open())
    {
        throwScanError("Configuration file failed to open.");
        return;
    }
    else
    {
        string line;
        int currentLine = 1;
        while (configscanner.good())
        {
            getline(configscanner,line);
			
            //assigns attributes based on line number
            switch(currentLine)
            {
                case 2:
					if(line.length() > 0)
					{
						filePathsEtc.currentEnvironment.x = atof(line.c_str());
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
						filePathsEtc.currentEnvironment.y = atof(line.c_str());
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
						filePathsEtc.currentEnvironment.z = atof(line.c_str());
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
						filePathsEtc.currentEnvironment.temperature = atof(line.c_str());
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
						filePathsEtc.currentEnvironment.maxTranslation = atof(line.c_str());
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
						filePathsEtc.numOfSteps = atoi(line.c_str());
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
						filePathsEtc.currentEnvironment.numOfMolecules = atoi(line.c_str());
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
						filePathsEtc.oplsuaparPath = line;
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
						filePathsEtc.zmatrixPath = line;
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
                        filePathsEtc.statePath = line;
                    }
                    else
					{
						throwScanError("Configuration file not well formed. Missing state file output path value.");
						return;
					}
                    break;
                case 20:
                    if(line.length() > 0){
                        filePathsEtc.stateOutputPath = line;
                    }
					else
					{
						throwScanError("Configuration file not well formed. Missing state file output path value.");
						return;
					}
                    break;
                case 22:
                    if(line.length() > 0){
                        filePathsEtc.pdbOutputPath = line;
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
						filePathsEtc.currentEnvironment.cutoff = atof(line.c_str());
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
						filePathsEtc.currentEnvironment.maxRotation = atof(line.c_str());
                    }
					else
					{
						throwScanError("Configuration file not well formed. Missing environment max rotation value.");
						return;
					}
                    break;
                case 28:
                    if(line.length() > 0){
						filePathsEtc.currentEnvironment.randomseed=atoi(line.c_str());
                    }
                	break;
                case 30:
                    if(line.length() > 0){
						// Convert to a zero-based index
						filePathsEtc.currentEnvironment.primaryAtomIndex=atoi(line.c_str()) - 1;
                    }
					else
					{
						throwScanError("Configuration file not well formed. Missing environment primary atom index value.");
						return;
					}
                	break;
            }
			
			currentLine++;
        }
    }
    
    readInConfigAlreadyDone = true;
}

/*
struct UtilitiesInfo
{
	Environment currentEnvironment; //The current working environment for the simulation
    string configPath; //The path to the main configuration file read in for the simulation
    unsigned int numOfSteps; //The number of steps to run the simulation
    string oplsuaparPath; //The path to the opls files containing additional geometry data, to be used (eventually) during simulation
    string zmatrixPath; //The path to the Z-matrix files to be used during simulation
    string statePath; //The path to the state information file to be used in the simulation
    string stateOutputPath; //The path where we write the state output files after simulation
    string pdbOutputPath; //The path where we write the pdb output files after simulation
    unsigned int cutoff; //The nonbonded cutoff distance.
};  */
  
void IOUtilities::readInConfig()
{
	if (readInConfigAlreadyDone)
	{
		//it's already been done during construction; do nothing
		}
	else
	{
		ifstream configscanner(filePathsEtc.configPath.c_str());
		if (! configscanner.is_open())
		{
			throwScanError("Configuration file failed to open.");
			return;
		}
		else
		{
			string line;
			int currentLine = 1;
			while (configscanner.good())
			{
				getline(configscanner,line);
			
				//assigns attributes based on line number
				switch(currentLine)
				{
					case 2:
						if(line.length() > 0)
						{
							filePathsEtc.currentEnvironment.x = atof(line.c_str());
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
							filePathsEtc.currentEnvironment.y = atof(line.c_str());
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
							filePathsEtc.currentEnvironment.z = atof(line.c_str());
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
							filePathsEtc.currentEnvironment.temperature = atof(line.c_str());
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
							filePathsEtc.currentEnvironment.maxTranslation = atof(line.c_str());
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
							filePathsEtc.numOfSteps = atoi(line.c_str());
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
							filePathsEtc.currentEnvironment.numOfMolecules = atoi(line.c_str());
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
							filePathsEtc.oplsuaparPath = line;
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
							filePathsEtc.zmatrixPath = line;
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
							filePathsEtc.statePath = line;
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing state file output path value.");
							return;
						}
						break;
					case 20:
						if(line.length() > 0){
							filePathsEtc.stateOutputPath = line;
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing state file output path value.");
							return;
						}
						break;
					case 22:
						if(line.length() > 0){
							filePathsEtc.pdbOutputPath = line;
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
							filePathsEtc.currentEnvironment.cutoff = atof(line.c_str());
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
							filePathsEtc.currentEnvironment.maxRotation = atof(line.c_str());
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing environment max rotation value.");
							return;
						}
						break;
					case 28:
						if(line.length() > 0){
							filePathsEtc.currentEnvironment.randomseed=atoi(line.c_str());
						}
						break;
					case 30:
						if(line.length() > 0){
							// Convert to a zero-based index
							filePathsEtc.currentEnvironment.primaryAtomIndex=atoi(line.c_str()) - 1;
						}
						else
						{
							throwScanError("Configuration file not well formed. Missing environment primary atom index value.");
							return;
						}
						break;
				}
			
				currentLine++;
			}
		}
	}
}

void IOUtilities::throwScanError(string message)
{

	cerr << endl << message << endl << "	Error Number: " << errno << endl <<endl;

	return;
}


/**
	This method is used to read in from a state file.
	//copied over from...SimBox, so retooling is necessary (variable references point to SimBox locations)
	
	@param StateFile - takes the location of the state file to be read in
*/

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
      inFile>>tmpenv.x>>tmpenv.y>>tmpenv.z>>tmpenv.maxTranslation>>tmpenv.numOfAtoms>>tmpenv.temperature>>tmpenv.cutoff;
    }
    
    if (memcmp(&tmpenv,destinationEnvironment,sizeof(Environment))!=0)
    {
       ss<<"Wrong state files,does not match other configfiles"<<endl;
       ss<<"x "<<tmpenv.x<<" "<<destinationEnvironment->x<<endl;
       ss<<"y "<<tmpenv.y<<" "<<destinationEnvironment->y<<endl;
       ss<<"z "<<tmpenv.z<<" "<<destinationEnvironment->z<<endl;
       ss<<"numOfAtoms "<<tmpenv.numOfAtoms<<" "<<destinationEnvironment->numOfAtoms<<endl;
       ss<<"temperature "<<tmpenv.temperature<<" "<<destinationEnvironment->temperature<<endl;
       ss<<"cutoff "<<tmpenv.cutoff<<" "<<destinationEnvironment->cutoff<<endl;
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
	
    ofstream outputFile;
    outputFile.open(pdbFile);
    int numOfMolecules=sourceEnvironment.numOfMolecules;
    outputFile << "REMARK Created by MCGPU" << endl;
    //int atomIndex = 0;
    for (int i = 0; i < numOfMolecules; i++)
    {
    	Molecule currentMol = sourceMoleculeCollection[i];    	
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


//this comes from metroUtil...
//this supports writing to the log with relative ease.
//said string contains almost all of the text you wish to be written to the file.
//opening and closing of the file will be done on the fly, and should be guaranteed once this method reaches its end.
//If stringstream is needed, you may call the overloaded version below, which will still call this version of the method
void writeToLog(string text,int stamp){
    string filename = "OutputLog";
	 ofstream logFile;
	 logFile.open(filename.c_str(),ios::out|ios::app);
	 
	 string hash ="";
	 time_t current_time;
    struct tm * time_info;
    char timeString[9];  // space for "HH:MM:SS\0"
	 
	 switch(stamp){
	     case START:
		      //The start of a new simulation
		      logFile << "\n\n\n\n\n\n" << endl;
				logFile << "======================================================================"<<endl;
				logFile << "                       Starting Simulation: ";				
            time(&current_time);
            time_info = localtime(&current_time);
            strftime(timeString, sizeof(timeString), "%H:%M:%S", time_info);
				logFile << timeString << endl;
				logFile << "----------------------------------------------------------------------"<<endl;
				break;
			case END: 
			   //The end of a running simulation
				logFile << "----------------------------------------------------------------------"<<endl;
				logFile << "                       Ending Simulation: ";
            time(&current_time);
            time_info = localtime(&current_time);
            strftime(timeString, sizeof(timeString), "%H:%M:%S", time_info);
				logFile << timeString << endl;
				logFile << "======================================================================"<<endl;
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
	 logFile << text << endl;
	 logFile.close();	 
}


//this method allows for writing to a given log file, with some measure of automation
//this overloaded version allows for a stringstream, instead of a normal string, to be input.
//said string contains almost all of the text you wish to be written to the file.
//opening and closing of the file will be done on the fly, and should be guaranteed once this method reaches its end.
void writeToLog(stringstream& ss, int stamp ){
    writeToLog(ss.str(),stamp);
	 ss.str(""); // clears the string steam...
	 ss.clear();
}
