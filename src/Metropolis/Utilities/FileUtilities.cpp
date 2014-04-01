#include "FileUtilities.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <exception>
#include <stdexcept>
#include <sstream>
#include "StructLibrary.h"
#include "Metropolis/Box.h"

using std::string;
using std::ifstream;


bool loadBoxData(string configpath, Box* box, long* steps)
{
	if (box == NULL)
	{
		std::cerr << "Error: loadBoxData(): Box is NULL" << std::endl;
		return false;
	}

	ConfigScanner config_scanner = ConfigScanner();
	if (!config_scanner.readInConfig(configpath))
	{
		std::cerr << "Error: loadBoxData(): Could not read config file" << std::endl;
		return false;
	}

    *steps = config_scanner.getSteps();
	box->environment = new Environment(config_scanner.getEnviro());

	OplsScanner opls_scanner = OplsScanner();
	if (!opls_scanner.readInOpls(config_scanner.getOplsusaparPath()))
	{
		std::cerr << "Error: loadBoxData(): Could not read OPLS file" << std::endl;
		return false;
	}

	ZmatrixScanner zmatrix_scanner = ZmatrixScanner();
	if (!zmatrix_scanner.readInZmatrix(config_scanner.getZmatrixPath(), &opls_scanner))
	{
		std::cerr << "Error: loadBoxData(): Could not read Z-Matrix file" << std::endl;
		return false;
	}

	// Fill the box with config stuff.
	if (!buildBoxData(&zmatrix_scanner, box))
	{
		std::cerr << "Error: loadBoxData(): Could not build box data" << std::endl;
		return false;
	}

	return true;
}

bool buildBoxData(ZmatrixScanner* zMatrixScan, Box* box)
{
	//Convert molecule vectors into an array
   int moleculeIndex = 0;
   int atomCount = 0;

   Environment* enviro = box->environment;

   vector<Molecule> molecVec = zMatrixScan->buildMolecule(atomCount);
   int molecMod = enviro->numOfMolecules % molecVec.size();
   if (molecMod != 0)
   {
       enviro->numOfMolecules += molecVec.size() - molecMod;
       cout << "Number of molecules not divisible by specified z-matrix. Changing number of molecules to: " << enviro->numOfMolecules << endl;
    }

    box->moleculeCount = enviro->numOfMolecules;
    box->molecules = (Molecule *)malloc(sizeof(Molecule) * box->moleculeCount);
    
    int molecDiv = enviro->numOfMolecules / molecVec.size();
    int molecTypenum=molecVec.size();
    
    int count[5];//sum up number of atoms,bonds,angles,dihedrals,hops
	//int * 
	//molecTables = new int
	//molecTables[molecVec.size()];
	//Table * tables;
	Table* tables = new Table[molecVec.size()];
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
  		
  		cout << "before table building. Number of atom "<< molec1.numOfAtoms<<endl;
  		
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
			cout << "atom1: "<< atom1-currentAtomCount <<endl;
			int atom2 = myHop->atom2;
			cout << "atom2: "<< atom2-currentAtomCount<<endl;
			cout << "hop: " << myHop->hop <<endl;
			table[atom1-currentAtomCount][atom2-currentAtomCount] = myHop->hop;
			table[atom2-currentAtomCount][atom1-currentAtomCount] = myHop->hop;
			myHop++;
		}
	  
	   cout << "after table building"<<endl;
	   //memset(table,0,sizeof(table));
	   //int table[molec1.numOfAtoms][molec1.numOfAtoms];
	   //cout << " this is " << j <<endl;
	   tables[j] = Table(table); //createTable is in metroUtil
	   currentAtomCount += molec1.numOfAtoms;
	   cout << "after table creation. Current atom cout: "<< currentAtomCount<<endl;
    }

    box->atomCount = molecDiv * count[0];
    box->bondCount = molecDiv * count[1];
    box->angleCount = molecDiv * count[2];
    box->dihedralCount = molecDiv * count[3];
    box->hopCount = molecDiv * count[4];

    std::cout << "Molecule Count: " << box->moleculeCount << std::endl;
    std::cout << "Atom Count: " << box->atomCount << std::endl;
    std::cout << "Bond Count: " << box->bondCount << std::endl;
    std::cout << "Angle Count: " << box->angleCount << std::endl;
    std::cout << "Dihedral Count: " << box->dihedralCount << std::endl;
    std::cout << "Hop Count: " << box->hopCount << std::endl;
     
    box->atoms 	   = (Atom *)malloc(sizeof(Atom)*box->atomCount);
    box->bonds     = (Bond *)malloc(sizeof(Bond)*box->bondCount);
    box->angles    = (Angle *)malloc(sizeof(Angle)*box->angleCount);
    box->dihedrals = (Dihedral *)malloc(sizeof(Dihedral)*box->dihedralCount);
    box->hops      = (Hop *)malloc(sizeof(Hop)*box->hopCount);

    memset(box->atoms,0,sizeof(Atom)*box->atomCount);
    memset(box->bonds,0,sizeof(Bond)*box->bondCount);
    memset(box->angles,0,sizeof(Angle)*box->angleCount);
    memset(box->dihedrals,0,sizeof(Dihedral)*box->dihedralCount);
    memset(box->hops,0,sizeof(Hop)*box->hopCount);

    //arrange first part of molecules
    memset(count,0,sizeof(count));
 	for(int j = 0; j < molecVec.size(); j++)
    {
 	      //Copy data from vector to molecule
        Molecule molec1 = molecVec[j];   

        box->molecules[j].atoms = (Atom *)(box->atoms+count[0]);
        box->molecules[j].bonds = (Bond *)(box->bonds+count[1]);
        box->molecules[j].angles = (Angle *)(box->angles+count[2]);
        box->molecules[j].dihedrals = (Dihedral *)(box->dihedrals+count[3]);
        box->molecules[j].hops = (Hop *)(box->hops+count[4]);

        box->molecules[j].id = molec1.id;
        box->molecules[j].numOfAtoms = molec1.numOfAtoms;
        box->molecules[j].numOfBonds = molec1.numOfBonds;
        box->molecules[j].numOfDihedrals = molec1.numOfDihedrals;
        box->molecules[j].numOfAngles = molec1.numOfAngles;
        box->molecules[j].numOfHops = molec1.numOfHops;

        count[0]+=molec1.numOfAtoms;
        count[1]+=molec1.numOfBonds;
        count[2]+=molec1.numOfAngles;
        count[3]+=molec1.numOfDihedrals;
        count[4]+=molec1.numOfHops;

        //get the atoms from the vector molecule
        for(int k = 0; k < molec1.numOfAtoms; k++)
        {
            box->molecules[j].atoms[k] = molec1.atoms[k];
        }               
           
        //assign bonds
        for(int k = 0; k < molec1.numOfBonds; k++)
        {
            box->molecules[j].bonds[k] = molec1.bonds[k];
        }

        //assign angles
        for(int k = 0; k < molec1.numOfAngles; k++)
        {
            box->molecules[j].angles[k] = molec1.angles[k];
        }

        //assign dihedrals
        for(int k = 0; k < molec1.numOfDihedrals; k++)
        {
            box->molecules[j].dihedrals[k] = molec1.dihedrals[k];
        }

        //assign hops zx add
        for(int k = 0; k < molec1.numOfHops; k++)
        {
            box->molecules[j].hops[k] = molec1.hops[k];
        }
    }
   
    for(int m = 1; m < molecDiv; m++)
    {
        int offset=m*molecTypenum;
    	memcpy(&(box->molecules[offset]),box->molecules,sizeof(Molecule)*molecTypenum);
    	for(int n=0;n<molecTypenum;n++)
        {
    		box->molecules[offset+n].id=offset+n;
            box->molecules[offset+n].atoms = box->molecules[n].atoms+count[0]*m;
            box->molecules[offset+n].bonds =  box->molecules[n].bonds+count[1]*m;
            box->molecules[offset+n].angles =  box->molecules[n].angles+count[2]*m;
            box->molecules[offset+n].dihedrals =  box->molecules[n].dihedrals+count[3]*m;
            box->molecules[offset+n].hops =  box->molecules[n].hops+count[4]*m;
        }
        
        memcpy(&(box->atoms[offset*count[0]]),box->atoms,sizeof(Atom)*count[0]);
        memcpy(&(box->bonds[offset*count[1]]),box->bonds,sizeof(Bond)*count[1]);
        memcpy(&(box->angles[offset*count[2]]),box->angles,sizeof(Angle)*count[2]);
        memcpy(&(box->dihedrals[offset*count[3]]),box->dihedrals,sizeof(Dihedral)*count[3]);
        memcpy(&(box->hops[offset*count[4]]),box->hops,sizeof(Hop)*count[4]);
        
        for(int k=0;k<count[0];k++)
        {
            box->atoms[offset*count[0]+k].id=offset*count[0]+k;
        }
        
        for(int k=0;k<count[1];k++)
        {
            box->bonds[offset*count[1]+k].atom1+=m*count[0];
            box->bonds[offset*count[1]+k].atom2+=m*count[0];
        }
        
        for(int k=0;k<count[2];k++)
        {
            box->angles[offset*count[2]+k].atom1+=m*count[0];
            box->angles[offset*count[2]+k].atom2+=m*count[0];
        }
        
        for(int k=0;k<count[3];k++)
        {
            box->dihedrals[offset*count[3]+k].atom1+=m*count[0];
            box->dihedrals[offset*count[3]+k].atom2+=m*count[0];
        }
        
        for(int k=0;k<count[4];k++)
        {
            box->hops[offset*count[4]+k].atom1+=m*count[0];
            box->hops[offset*count[4]+k].atom2+=m*count[0];
        }
    }
     
    enviro->numOfAtoms = count[0]*molecDiv;
	std::cout << "Molecules Created into an Array" << std::endl;
     
    // TODO: Right now we can only read from a zmatrix file and not a state file
    // In the future we need to make it possible.
    
    std::cout << "Assigning Molecule Positions..." << std::endl;
    if (!generatefccBox(box)) //generate fcc lattice box
    {
    	std::cerr << "Error: buildBoxData(): Could not generate FCC box" << std::endl;
    	return false;
    }
    std::cout << "Finished Assigning Molecule Positions" << std::endl;

    return true;
}

bool generatefccBox(Box* box)
{
	if (box->environment == NULL || box->molecules == NULL)
	{
		std::cerr << "Error: generatefccBox(): Box environment or molecules array NULL" << std::endl;
		return false;
	}
	
	double cells, dcells, cellL, halfcellL;
	Environment* enviro = box->environment;
	Molecule* molecules = box->molecules;
	
	//Determine the number of unit cells in each coordinate direction
	dcells = pow(0.25 * (double) enviro->numOfMolecules, 1.0/3.0);
	cells = (int)(dcells + 0.5);
		
	//Check if numOfMolecules is a non-fcc number of molecules
	//and increase the number of cells if necessary
	while((4 * cells * cells * cells) < enviro->numOfMolecules)
    {
		cells++;
    }
			
	//Determine length of unit cell
	cellL = enviro->x/ (double) cells;
	halfcellL = 0.5 * cellL;
	
	//Construct the unit cell
	for (int j = 0; j < molecules[0].numOfAtoms; j++)
    {
    	molecules[0].atoms[j].x += 0.0;
        molecules[0].atoms[j].y += 0.0;
        molecules[0].atoms[j].z += 0.0;
	}
	
	for (int j = 0; j < molecules[1].numOfAtoms; j++)
    {
    	molecules[1].atoms[j].x += halfcellL;
        molecules[1].atoms[j].y += halfcellL;
        molecules[1].atoms[j].z += 0.0;
    }
    
    for (int j = 0; j < molecules[2].numOfAtoms; j++)
    {	
        molecules[2].atoms[j].x += 0.0;
        molecules[2].atoms[j].y += halfcellL;
        molecules[2].atoms[j].z += halfcellL;
    }
    
    for (int j = 0; j < molecules[3].numOfAtoms; j++)
    {
        molecules[3].atoms[j].x += halfcellL;
        molecules[3].atoms[j].y += 0.0;
        molecules[3].atoms[j].z += halfcellL;
    }
    
	//Init all other molecules to initial coordinates
	//Build the lattice from the unit cell by repeatedly translating
	//the four vectors of the unit cell through a distance cellL in
	//the x, y, and z directions
	for(int i = 4; i < enviro->numOfMolecules; i++)
    {
		for (int j = 0; j < molecules[i].numOfAtoms; j++)
        {
			molecules[i].atoms[j].x += 0.0;
    		molecules[i].atoms[j].y += 0.0;
   	 		molecules[i].atoms[j].z += 0.0;
   	 	}		
	}
	
	int offset = 0;
	for(int z = 1; z <= cells; z++)
		for(int y = 1; y <= cells; y++)
			for(int x = 1; x <= cells; x++)
            {
				for(int a = 0; a < 4; a++)
                {
					int i = a + offset;
					if(i < enviro->numOfMolecules)
                    {								
						for (int j = 0; j < molecules[i].numOfAtoms; j++)
                        {
							molecules[i].atoms[j].x = molecules[a].atoms[j].x + cellL * (x-1);
							molecules[i].atoms[j].y = molecules[a].atoms[j].y + cellL * (y-1);
							molecules[i].atoms[j].z = molecules[a].atoms[j].z + cellL * (z-1);
						}
					}
				}
				offset += 4;
			}
	
	//Shift center of box to the origin
	for(int i = 0; i < enviro->numOfMolecules; i++)
    {
		for (int j = 0; j < molecules[i].numOfAtoms; j++)
        {
			molecules[i].atoms[j].x -= halfcellL;
			molecules[i].atoms[j].y -= halfcellL;
			molecules[i].atoms[j].z -= halfcellL;
		}
	}

	return true;
}

ConfigScanner::ConfigScanner()
{
	enviro = Environment();
	numOfSteps = 0;
}

bool ConfigScanner::readInConfig(string configpath)
{
	enviro = Environment();
	numOfSteps = 0;

	configPath = configpath;
	if (configpath.empty())
	{
		std::cerr << "Configuration File path is empty" << std::endl;
		return false;
	}

    ifstream configscanner(configpath.c_str());
    if (!configscanner.is_open())
    {
        std::cerr << "Unable to open configuration file (" << configpath << ")" << std::endl;
        return false;
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
						enviro.x = atof(line.c_str());
                    }
					else
					{
						std::cerr << "Error: Config File: Missing environment X value." << std::endl;
						return false;
					}
                    break;
                case 3:
					if(line.length() > 0)
					{
						enviro.y = atof(line.c_str());
                    }
					else
					{
						std::cerr << "Error: Config File: Missing environment Y value." << std::endl;
						return false;
					}
                    break;
                case 4:
					if(line.length() > 0)
					{
						enviro.z = atof(line.c_str());
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
						enviro.temp = atof(line.c_str());
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
						enviro.maxTranslation = atof(line.c_str());
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
						enviro.numOfMolecules = atoi(line.c_str());
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
						return false;
					}
                    break;
                case 18:
                    if(line.length() > 0)
					{
                        statePath = line;
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
						enviro.cutoff = atof(line.c_str());
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
						enviro.maxRotation = atof(line.c_str());
                    }
					else
					{
						throwScanError("Configuration file not well formed. Missing environment max rotation value.");
						return false;
					}
                    break;
                case 28:
                    if(line.length() > 0){
						enviro.randomseed=atoi(line.c_str());
                    }
                	break;
                case 30:
                    if(line.length() > 0){
						// Convert to a zero-based index
						enviro.primaryAtomIndex=atoi(line.c_str()) - 1;
                    }
					else
					{
						throwScanError("Configuration file not well formed. Missing environment primary atom index value.");
						return false;
					}
                	break;
            }
			
			currentLine++;
        }
    }

    configscanner.close();

    return true;
}

void ConfigScanner::throwScanError(string message)
{

	std::cerr << message << std::endl;
}

Environment* ConfigScanner::getEnviro()
{
    return &enviro;
}

string ConfigScanner::getConfigPath()
{
    return configPath;
}

long ConfigScanner::getSteps()
{
    return numOfSteps;
}
    
string ConfigScanner::getOplsusaparPath()
{
    return oplsuaparPath;
}

string ConfigScanner::getZmatrixPath()
{
    return zmatrixPath;
}

string ConfigScanner::getStatePath()
{
    return statePath;
}

string ConfigScanner::getStateOutputPath()
{
    return stateOutputPath;
}

string ConfigScanner::getPdbOutputPath()
{
    return pdbOutputPath;
}

// ============================================================================
// ============================= OPLS SCANNER =================================
// ============================================================================


OplsScanner::OplsScanner()
{
}

OplsScanner::~OplsScanner()
{
    oplsTable.clear();
    fourierTable.clear();
}

bool OplsScanner::readInOpls(string filename)
{
	fileName = filename;
	if (filename.empty())
	{
		std::cerr << "Error: readInOpls(): empty filename given" << std::endl;
		return false;
	}

    int numOfLines=0;
    ifstream oplsScanner(filename.c_str());
    if (!oplsScanner.is_open())
    {
    	std::cerr << "Error: readInOpls(): could not open file (" << filename << ")" << std::endl;
        return false;
    }
    else 
    {
        string line; 
        while (oplsScanner.good())
        {
            numOfLines++;
            getline(oplsScanner,line);

            //check if it is a commented line,
            //or if it is a title line
            try
            {
                if(line.at(0) != '#' && numOfLines > 1)
                    addLineToTable(line,numOfLines);
            }
            catch (std::out_of_range& e)
            {
            	// Eat the exception and continue...why aren't we failing?
            }
        }
        oplsScanner.close();
        logErrors();
    }

    return true;
}

void OplsScanner::addLineToTable(string line, int numOfLines)
{
    string hashNum;
    int secCol;
    double charge,sigma,epsilon;
    string name, extra;
    stringstream ss(line);

    //check to see what format it is opls, V value, or neither
    int format = checkFormat(line);
              
    if(format == 1)
    {      	
        ss >> hashNum >> secCol >> name >> charge >> sigma >> epsilon;
        char *atomtype = (char*)name.c_str(); 
         
        Atom temp = createAtom(0, -1, -1, -1, sigma, epsilon, charge, *atomtype);
        pair<map<string,Atom>::iterator,bool> ret;
        ret = oplsTable.insert( pair<string,Atom>(hashNum,temp) );

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
        pair<map<string,Fourier>::iterator,bool> ret2;
        ret2 = fourierTable.insert( pair<string,Fourier>(hashNum,vValues) );

        if (ret2.second==false)
        {
            errHashesFourier.push_back(hashNum);	
        }	  
    }
    else
    {
	     errLines.push_back(numOfLines);
    }
}

int OplsScanner::checkFormat(string line)
{   	 
    int hashNum, secCol;
    double charge,sigma,epsilon;
    string name, extra;
    stringstream iss(line);

    double v1,v2,v3,v4;
    stringstream issw(line);

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

void OplsScanner::logErrors()
{
    stringstream output;
    // See if there were any errors
    if(errLines.empty() || errHashes.empty()|| errHashesFourier.empty())
    {
	     //Errors in the format
		  output<<"Errors found in the OPLS file: "<< fileName<<endl;
        if(!errLines.empty())
        {
		      output << "Found Errors in the Format of the following Lines: " << endl;
				for(int a=0; a<errLines.size(); a++)
                {
				    if(a%10==0 && a!=0) //ten per line
                    {
					     output << endl;
                    }
				    output << errLines[a]<< " ";
				}
				output << endl<< endl;
		  }
		  if(!errHashes.empty())
          {
		      output << "Error - The following OPLS values existed more than once: " << endl;
				for(int a=0; a<errHashes.size(); a++)
                {
				    if(a%10==0 && a!=0) //ten per line
                    {
					     output << endl;
                    }
				    output << errHashes[a]<< " ";
				}
				output << endl<< endl;
		  }
		  if(!errHashesFourier.empty())
          {
		      output << "Error - The following Fourier Coefficent values existed more than once: " << endl;
				for(int a=0; a<errHashesFourier.size(); a++)
                {
				    if(a%10==0 && a!=0) //ten per line
                    {
					     output << endl;
                    }
				    output << errHashesFourier[a]<< " ";
				}
				output << endl<< endl;
		  }
		  writeToLog(output,OPLS);
	}
}

Atom OplsScanner::getAtom(string hashNum)
{
    if(oplsTable.count(hashNum)>0 )
    {
        return oplsTable[hashNum];
	}
	else
    {
	    cerr << "Index does not exist: "<< hashNum <<endl;
		return createAtom(0, -1, -1, -1, -1, -1, -1, NULL);
	}
}

double OplsScanner::getSigma(string hashNum)
{
    if(oplsTable.count(hashNum)>0 )
    {
        Atom temp = oplsTable[hashNum];
        return temp.sigma;
    }
    else
    {
        cerr << "Index does not exist: "<< hashNum <<endl;
        return -1;
    }
}

double OplsScanner::getEpsilon(string hashNum)
{
    if(oplsTable.count(hashNum)>0 )
    {
        Atom temp = oplsTable[hashNum];
        return temp.epsilon;
    }
    else
    {
        cerr << "Index does not exist: "<< hashNum <<endl;
        return -1;
    }
}

double OplsScanner::getCharge(string hashNum)
{
    if(oplsTable.count(hashNum)>0 )
    {
        Atom temp = oplsTable[hashNum];
        return temp.charge;
    }
    else
    {
        cerr << "Index does not exist: "<< hashNum <<endl;
        return -1;
    }
}

Fourier OplsScanner::getFourier(string hashNum)
{
    if(fourierTable.count(hashNum)>0 )
    {
        Fourier temp = fourierTable[hashNum];
        return temp;
    }
    else
    {	    
        cerr << "Index does not exist: "<< hashNum <<endl;
        Fourier temp ={-1,-1,-1,-1};
        return temp;
    }
}

// ============================================================================
// ======================== Z-Matrix Scanner ==================================
// ============================================================================

ZmatrixScanner::ZmatrixScanner()
{
    oplsScanner = NULL;
    startNewMolecule = false;
}

ZmatrixScanner::~ZmatrixScanner()
{
}

bool ZmatrixScanner::readInZmatrix(string filename, OplsScanner* scanner)
{
	fileName = filename;
	oplsScanner = scanner;
	startNewMolecule = false;

	if (filename.empty())
	{
		std::cerr << "Error: readInZmatrix(): Given filename is NULL" << std::endl;
		return false;
	}

    stringstream output;
    int numOfLines=0;

    ifstream zmatrixScanner(fileName.c_str());
    if( !zmatrixScanner.is_open() )
    {
    	std::cerr << "Error: Unable to open Z-Matrix file (" << fileName << ")" << std::endl;
        return false;
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
    }

    zmatrixScanner.close();    	 

    return true;
}

void ZmatrixScanner::parseLine(string line, int numOfLines)
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
    
int ZmatrixScanner::checkFormat(string line)
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

void ZmatrixScanner::handleZAdditions(string line, int cmdFormat)
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

vector<Hop> ZmatrixScanner::calculateHops(Molecule molec)
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
				Hop tempHop = Hop(atom1+startId,atom2+startId,distance); //+startId because atoms may not start at 1
                newHops.push_back(tempHop);			
            }  		      
        }
    }

    return newHops; 
}

bool ZmatrixScanner::contains(vector<int> &vect, int item)
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


int ZmatrixScanner::findHopDistance(int atom1,int atom2,int size, int **graph)
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

void ZmatrixScanner::buildAdjacencyMatrix(int **&graph, Molecule molec)
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

vector<Molecule> ZmatrixScanner::buildMolecule(int startingID)
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


        Molecule molecCopy = Molecule(-1,atomCopy, angleCopy, bondCopy, dihedCopy, hopCopy, 
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

// ============================================================================
// ======================= Logging Functions ==================================
// ============================================================================

void writeToLog(string text,int stamp)
{
    string filename = "OutputLog";
	ofstream logFile;
	logFile.open(filename.c_str(),ios::out|ios::app);
	 
	string hash ="";
	time_t current_time;
    struct tm * time_info;
    char timeString[9];  // space for "HH:MM:SS\0"
	 
	switch(stamp)
    {
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

void writeToLog(stringstream& ss, int stamp)
{
    writeToLog(ss.str(),stamp);
	ss.str(""); // clears the string steam...
	ss.clear();
}