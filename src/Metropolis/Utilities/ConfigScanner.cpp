#include "FileUtilities.h"

#include <exception>
#include <stdexcept>
#include <time.h>
//#include <sstream>
#include "Parsing.h"
#include "StructLibrary.h"
#include "Metropolis/Box.h"
#include "Metropolis/SimulationArgs.h"


ConfigScanner::ConfigScanner() {
	enviro = Environment();
	numOfSteps = 0;
	useStatefileSetup = false;
	useZMatrixSetup = false;
}


bool ConfigScanner::readInConfig(string configpath) {
	isSafeToContinue = true;
	enviro = Environment();
	numOfSteps = 0;

	configPath = configpath;
	if (configPath.empty()) {
		std::cerr << "Configuration File path is empty" << std::endl;
		return false;
	}

  ifstream configscanner(configPath.c_str());
  if (!configscanner.is_open()) {
    std::cerr << "Unable to open configuration file (" << configPath << ")"
              << std::endl;
    return false;
  } else {
    string line;
    int currentLine = 1;

    while (configscanner.good()) {
      getline(configscanner,line);
			trim(line);
      //assigns attributes based on line number
      switch(currentLine) {
        case 2:     // read in x size of environment
					if (line.length() > 0) {
						enviro.x = atof(line.c_str());
          } else {
						throwScanError("Error: Config File: Missing environment X value.");
						return false;
					}
          break;
        case 3:     // read in y size of environment
					if (line.length() > 0) {
						enviro.y = atof(line.c_str());
          } else {
						throwScanError("Error: Config File: Missing environment Y value.");
						return false;
					}
          break;
        case 4:     // read in z size of environment
					if (line.length() > 0) {
						enviro.z = atof(line.c_str());
          } else {
						throwScanError("Configuration file not well formed. Missing environment z value.");
						return false;
					}
          break;
        case 6:     // read in the temperature of the environment
					if (line.length() > 0) {
						enviro.temp = atof(line.c_str());
          } else {
						throwScanError("Configuration file not well formed. Missing environment temperature value.");
						return false;
					}
          break;
        case 8:     // read in the maximum translation an atom can have.
					if (line.length() > 0) {
						enviro.maxTranslation = atof(line.c_str());
          } else {
						throwScanError("Configuration file not well formed. Missing environment max translation value.");
						return false;
					}
          break;
        case 10:    // read in the number of steps to perform.
					if (line.length() > 0) {
						numOfSteps = atoi(line.c_str());
          } else {
						throwScanError("Configuration file not well formed. Missing number of steps value.");
						return false;
					}
          break;
        case 12:    // read in the number of molecules.
					if (line.length() > 0) {
						enviro.numOfMolecules = atoi(line.c_str());
          } else {
						throwScanError("Configuration file not well formed. Missing number of molecules value.");
						return false;
					}
          break;
        case 14:    // read in the oplsua.par path.
					if (line.length() > 0) {
						oplsuaparPath = line;
          }	else {
						throwScanError("Configuration file not well formed. Missing oplsuapar path value.");
						return false;
					}
          break;
        case 16:    // read in the z matrix path.
					if (line.length() > 0) {
						zmatrixPath = line;
						useZMatrixSetup = true;
						useStatefileSetup = true;
					} else {
						throwScanError("INFO: Configuration file not well formed. Missing z-matrix path value. Attempting to rely on prior state information...");
						isSafeToContinue = false; //now that it is false, check if there is a state file
						useZMatrixSetup = false;
						useStatefileSetup = true;
					}
					break;
        case 18:  // read in the state file path.
          if (line.length() > 0) {
            statePath = line;
            useZMatrixSetup = false; //we will go ahead and use the statefile, since the path was provided
            useStatefileSetup = true;
          } else if (isSafeToContinue == false) { // If no z matrix was provided, and no state file, there is a problem.
						throwScanError("Configuration file not well formed. Missing value pointing to prior state file path. Cannot safely continue with program execution.");
						useStatefileSetup = false; //we can't even hope to find the state file, so don't use it
						useZMatrixSetup = false;
						return false; //preferable to simply exiting, as we want to give the program a chance to do...whatever?
					} else { // If we don't have a state file, but do have a z-matrix, use the z-matrix.
						throwScanError("INFO: Value pointing to prior state file path not found in main config file. Environment setup defaulting to clean simulation.");
						useZMatrixSetup = true; //revert to using the Zmatrix file for all setup
						useStatefileSetup = false;
						isSafeToContinue = true;
					}
					break;
        case 20:    // read in the path to output state files to.
          if (line.length() > 0) {
            stateOutputPath = line;
          } else {
						throwScanError("Configuration file not well formed. Missing state file output path value.");
						return false;
					}
          break;
        case 22:    // read in the path to output pdb files to.
          if (line.length() > 0) {
            pdbOutputPath = line;
          } else {
						throwScanError("Configuration file not well formed. Missing PDB output path value.");
						return false;
					}
          break;
        case 24:    // read in the energy calculation cutoff value.
					if (line.length() > 0) {
						enviro.cutoff = atof(line.c_str());
          } else {
						throwScanError("Configuration file not well formed. Missing environment cutoff value.");
						return false;
					}
          break;
        case 26:  // read in the max rotation value.
					if (line.length() > 0) {
				    enviro.maxRotation = atof(line.c_str());
          } else {
						throwScanError("Configuration file not well formed. Missing environment max rotation value.");
						return false;
					}
          break;
        case 28:    // read in the random seed.
          if (line.length() > 0) {
						enviro.randomseed=atoi(line.c_str());
          } else {
						enviro.randomseed = (unsigned int) time(NULL);
						//throwScanError("Configuration file not well formed. Missing random seed value.");
						//return false;
					}
          break;
        case 30:  // Get the primary atom index definitions.
          if (line.length() > 0) {
						// Convert to a zero-based index
						parsePrimaryIndexDefinitions(line);
          } else {
						throwScanError("Configuration file not well formed. Missing environment primary atom index value.");
						return false;
					}
          break;
		case 32: // Simulation Name
		  if (line.length() > 0) {
		  	simName = line;
		  }
      }

			currentLine++;
    }
  }

  configscanner.close();

  return true;
}

void ConfigScanner::parsePrimaryIndexDefinitions(string definitions) {
  *enviro.primaryAtomIndexConfigLine = definitions;
  char *indexVector;
  char *charLine = (char *) malloc(sizeof(char) * (definitions.size()+1));
  strcpy(charLine, definitions.c_str());
  indexVector = strtok(charLine, ",");

  for (int i = 0; indexVector ; i++) {
    std::string sIndexVector = indexVector;
	  enviro.primaryAtomIndexDefinitions++;
    char* c;
    if ((c=strchr(indexVector, '[')) != NULL) {

      indexVector = (indexVector + 1);
	    while ((c=strchr(indexVector, ']')) == NULL) {
		    int currentPrimaryIndex = atoi(indexVector);
		    (*(enviro.primaryAtomIndexArray)).push_back(new std::vector<int>);
		    (*(*(enviro.primaryAtomIndexArray))[i]).push_back(currentPrimaryIndex - 1);
		    indexVector = strtok(NULL, ",");
	    }

	    *c = '\0';
	    int currentPrimaryIndex = atoi(indexVector);
	    (*(enviro.primaryAtomIndexArray)).push_back(new std::vector<int>);
	    (*(*(enviro.primaryAtomIndexArray))[i]).push_back(currentPrimaryIndex - 1);
	    indexVector =strtok(NULL, ",");
	    continue;
	  }

    int currentPrimaryIndex = atoi(indexVector);
	  (*(enviro.primaryAtomIndexArray)).push_back(new std::vector<int>);
	  (*(*(enviro.primaryAtomIndexArray))[i]).push_back(currentPrimaryIndex - 1);
	  indexVector = strtok(NULL, ",");
  }

  free (charLine);
}

void ConfigScanner::throwScanError(string message) {
	std::cerr << message << std::endl;
}

Environment* ConfigScanner::getEnviro() {
  return &enviro;
}

string ConfigScanner::getConfigPath() {
  return configPath;
}

long ConfigScanner::getSteps() {
  return numOfSteps;
}

string ConfigScanner::getOplsusaparPath() {
  return oplsuaparPath;
}

string ConfigScanner::getZmatrixPath() {
  return zmatrixPath;
}

string ConfigScanner::getStatePath() {
  return statePath;
}

string ConfigScanner::getStateOutputPath() {
  return stateOutputPath;
}

string ConfigScanner::getPdbOutputPath() {
  return pdbOutputPath;
}

string ConfigScanner::getSimulationName() {
  return simName;
}
