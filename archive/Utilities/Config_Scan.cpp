
/*!\file
  \brief Class used to read and write configuration files.
  \author Alexander Luchs, Riley Spahn, Seth Wooten, and Orlando Acevedo
 
 */
#include "Config_Scan.h"
#include "errno.h"

Config_Scan::Config_Scan(string configPath){
    configpath = configPath;
    memset(&enviro,0,sizeof(Environment));
    numOfSteps=0;
}

void Config_Scan::readInConfig()
{
    ifstream configscanner(configpath.c_str());
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
						enviro.x = atof(line.c_str());
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
						enviro.y = atof(line.c_str());
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
						enviro.z = atof(line.c_str());
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
						enviro.temperature = atof(line.c_str());
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
						enviro.maxTranslation = atof(line.c_str());
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
						numOfSteps = atoi(line.c_str());
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
						enviro.numOfMolecules = atoi(line.c_str());
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
						oplsuaparPath = line;
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
						zmatrixPath = line;
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
						return;
					}
                    break;
                case 22:
                    if(line.length() > 0){
                        pdbOutputPath = line;
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
						enviro.cutoff = atof(line.c_str());
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
						enviro.maxRotation = atof(line.c_str());
                    }
					else
					{
						throwScanError("Configuration file not well formed. Missing environment max rotation value.");
						return;
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
						return;
					}
                	break;
            }
			
			currentLine++;
        }
    }
}

void Config_Scan::throwScanError(string message)
{

	cerr << endl << message << endl << "	Error Number: " << errno << endl <<endl;

	return;
}

Environment *Config_Scan::getEnviro()
{
    return &enviro;
}

string Config_Scan::getConfigPath()
{
    return configpath;
}

long Config_Scan::getSteps()
{
    return numOfSteps;
}
    
string Config_Scan::getOplsusaparPath()
{
    return oplsuaparPath;
}

string Config_Scan::getZmatrixPath()
{
    return zmatrixPath;
}

string Config_Scan::getStatePath()
{
    return statePath;
}

string Config_Scan::getStateOutputPath()
{
    return stateOutputPath;
}

string Config_Scan::getPdbOutputPath()
{
    return pdbOutputPath;
}
