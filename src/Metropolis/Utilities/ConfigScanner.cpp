#include <exception>
#include <stdexcept>
#include <time.h>

#include "FileUtilities.h"
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
  }
  string line;

  while (configscanner.good()) {
    getline(configscanner, line);
    trim(line);

    // Ignore blank lines and comments
    if (line.empty() || line[0] == '#' || line[0] == ';') {
      continue;
    }

    // Sections are supported, but don't have any additional functionality
    if (line[0] == '[' && line.find(']') == line.length()) {
      continue;
    }

    int splitPos = line.find('=');
    if (splitPos == -1) {
      throwScanError("Error: Configuration File not well formed.\n"
        "Offending line: \"" + line + "\"");
      return false;
    }
    string key = line.substr(0, splitPos);
    string value = line.substr(splitPos+1, line.length());

    // Assign attributes based on key
    if (key == "x") {
      if (value.length() > 0) {
        enviro.x = atof(value.c_str());
      } else {
        throwScanError("Error: Config File: Missing environment X value.");
        return false;
      }
    } else if (key == "y") {
      if (value.length() > 0) {
        enviro.y = atof(value.c_str());
      } else {
        throwScanError("Error: Config File: Missing environment Y value.");
        return false;
      }
    } else if (key == "z") {
      if (value.length() > 0) {
        enviro.z = atof(value.c_str());
      } else {
        throwScanError("Configuration file not well formed. Missing environment z value.");
        return false;
      }
    } else if (key == "temp") {
      if (value.length() > 0) {
        enviro.temp = atof(value.c_str());
      } else {
        throwScanError("Configuration file not well formed. Missing environment temperature value.");
        return false;
      }
    } else if (key == "max-translation") {
      if (value.length() > 0) {
        enviro.maxTranslation = atof(value.c_str());
      } else {
        throwScanError("Configuration file not well formed. Missing environment max translation value.");
        return false;
        }
    } else if (key == "steps") {
        if (value.length() > 0) {
          numOfSteps = atoi(value.c_str());
        } else {
          throwScanError("Configuration file not well formed. Missing number of steps value.");
          return false;
        }
    } else if (key == "molecules") {
      if (value.length() > 0) {
        enviro.numOfMolecules = atoi(value.c_str());
      } else {
        throwScanError("Configuration file not well formed. Missing number of molecules value.");
        return false;
      }
    } else if (key == "opla.par") {
      if (value.length() > 0) {
        oplsuaparPath = value;
      } else {
        throwScanError("Configuration file not well formed. Missing oplsuapar path value.");
        return false;
      }
    } else if (key == "z-matrix") {
      if (value.length() > 0) {
        zmatrixPath = value;
        useZMatrixSetup = true;
        useStatefileSetup = true;
      } else {
        throwScanError("INFO: Configuration file not well formed. Missing z-matrix path value. Attempting to rely on prior state information...");
        isSafeToContinue = false; //now that it is false, check if there is a state file
        useZMatrixSetup = false;
        useStatefileSetup = true;
      }
    } else if (key == "state-input") {
      if (value.length() > 0) {
        statePath = value;
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
    } else if (key == "state-output") {
      if (value.length() > 0) {
        stateOutputPath = value;
      } else {
        throwScanError("Configuration file not well formed. Missing state file output path value.");
        return false;
      }
    } else if (key == "pdb-output") {
      if (value.length() > 0) {
        pdbOutputPath = value;
      } else {
        throwScanError("Configuration file not well formed. Missing PDB output path value.");
        return false;
      }
    } else if (key == "cutoff") {
      if (value.length() > 0) {
        enviro.cutoff = atof(value.c_str());
      } else {
        throwScanError("Configuration file not well formed. Missing environment cutoff value.");
        return false;
      }
    } else if (key == "max-rotation") {
      if (value.length() > 0) {
        enviro.maxRotation = atof(value.c_str());
      } else {
        throwScanError("Configuration file not well formed. Missing environment max rotation value.");
        return false;
      }
    } else if (key == "random-seed") {
      if (value.length() > 0) {
        enviro.randomseed=atoi(value.c_str());
      } else {
        enviro.randomseed = (unsigned int) time(NULL);
      }
    } else if (key == "primary-atom") {
      if (value.length() > 0) {
        // Convert to a zero-based index
        parsePrimaryIndexDefinitions(value);
      } else {
        throwScanError("Configuration file not well formed. Missing "
                       "environment primary atom index value.");
        return false;
      }
    } else if (key == "sim-name") {
      if (value.length() > 0) {
        simName = value;
      }
    } else if (key == "strategy") {
      if (value.length() > 0) {
        strategy = value;
      }
    } else if (key == "max-bond-delta") {
      enviro.maxBondDelta = atof(value.c_str());
    } else if (key == "max-angle-delta") {
      enviro.maxAngleDelta = atof(value.c_str());
    } else {
      throwScanError("Unexpected key encountered: " + key);
      return false;
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

string ConfigScanner::getStrategy() {
  return strategy;
}
