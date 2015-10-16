#include "FileUtilities.h"

#include <exception>
#include <stdexcept>
//#include <sstream>
#include "Parsing.h"
#include "StructLibrary.h"
#include "Metropolis/Box.h"
#include "Metropolis/SimulationArgs.h"

StateScanner::StateScanner(string filename) {
	universal_filename = filename;
}

StateScanner::~StateScanner() {
}

Environment* StateScanner::readInEnvironment() {
	std::string filename = universal_filename;

  ifstream inFile;
  inFile.open(filename.c_str());
  string line;
  Environment* environment;

  if (inFile.is_open()) {
    getline(inFile, line);
    environment = getEnvironmentFromLine(line);
  }

  return environment;
}

long StateScanner::readInStepNumber() {
  std::string filename = universal_filename;

  long startingStep = 0;
  ifstream inFile(filename.c_str());
  string line;

  if (inFile.is_open()) {
    getline(inFile, line); // Environment
    getline(inFile, line); // Step number

    if (!fromString<long>(line, startingStep))
      return -1;

    return startingStep;
  }

  return -1;
}

vector<Molecule> StateScanner::readInMolecules() {
	std::string filename = universal_filename;

  vector<Molecule> molecules;
  ifstream inFile(filename.c_str());
  string line;

  if (inFile.is_open()) {
    vector<Bond> bonds;
    vector<Angle> angles;
    vector<Atom> atoms;
    vector<Dihedral> dihedrals;
    vector<Hop> hops;

    //the third line starts the first molecule
    getline(inFile, line); // envrionment
    getline(inFile, line); // step number
    getline(inFile, line); //blank

    Molecule currentMol;
    int section = 0; // 0 = id, 1 = type, 2 = atom, 3 = bond, 4 = dihedral, 5 = hop, 6 = angle 7 = name
    int molNum = 0;

    while (inFile.good()) {
      getline(inFile, line);
      string hold = line.substr(0, 2);

      switch (section) {
        case 0: // id
          if (hold.compare("= ") == 0) {
            section++;
          } else {
            currentMol.id = atoi(line.c_str());
          }
          break;
        case 1: // type
          if (hold.compare("= ") == 0) {
            section++;
          } else {
            currentMol.type = atoi(line.c_str());
          }
          break;
        case 2: // atom
          if (hold.compare("= ") == 0) {
            section++;
          } else {
            atoms.push_back(getAtomFromLine(line));
          }
          break;
        case 3: // bond
          if (hold.compare("= ") == 0) {
            section++;
          } else {
					  bonds.push_back(getBondFromLine(line));
          }
          break;
        case 4: // dihedral
          if (hold.compare("= ") == 0) {
            section++;
          } else {
            dihedrals.push_back(getDihedralFromLine(line));
          }
          break;
        case 5: // hop
          if (hold.compare("= ") == 0) {
            section++;
          } else {
            hops.push_back(getHopFromLine(line));
          }
          break;
        case 6: // angle
          if (hold.compare("==") == 0) {
            section = 0;
            molNum++;

            //convert all vectors to arrays
            Bond *bondArray = (Bond *) malloc(sizeof(Bond) * bonds.size());
            Angle *angleArray = (Angle *) malloc(sizeof(Angle) * angles.size());
            Atom *atomArray = (Atom *) malloc(sizeof(Atom) * atoms.size());
            Dihedral *dihedralArray = (Dihedral *) malloc(sizeof(Dihedral) * dihedrals.size());
            Hop *hopArray = (Hop *) malloc(sizeof(Hop) * hops.size());

            for (int i = 0; i < bonds.size(); i++) {
              bondArray[i] = bonds[i];
            }

						for (int i = 0; i < angles.size(); i++) {
              angleArray[i] = angles[i];
            }

            for (int i = 0; i < atoms.size(); i++) {
              atomArray[i] = atoms[i];
            }

            for (int i = 0; i < dihedrals.size(); i++) {
              dihedralArray[i] = dihedrals[i];
            }

            for (int i = 0; i < hops.size(); i++) {
              hopArray[i] = hops[i];
            }

            //assign arrays to molecule
            currentMol.atoms = atomArray;
            currentMol.numOfAtoms = atoms.size();

            currentMol.bonds = bondArray;
            currentMol.numOfBonds = bonds.size();

            currentMol.angles = angleArray;
            currentMol.numOfAngles = angles.size();

            currentMol.dihedrals = dihedralArray;
            currentMol.numOfDihedrals = dihedrals.size();

            currentMol.hops = hopArray;
            currentMol.numOfHops = hops.size();

            //add molecule to array of returned molecules
            molecules.push_back(currentMol);

            Molecule newMolec;
            currentMol = newMolec;

            dihedrals.clear();
            atoms.clear();
            bonds.clear();
            angles.clear();
            hops.clear();
          } else {
            angles.push_back(getAngleFromLine(line));
          }
          break;
      }
    }

  } else {
    return molecules;
  }

  return molecules;
}


Angle StateScanner::getAngleFromLine(string line) {
  Angle angle = Angle();

  char *tokens;
  char *charLine = (char *) malloc(sizeof(char) * (line.size()+1));
  strcpy(charLine, line.c_str());

  tokens = strtok(charLine, " ");
  int tokenNumber = 0;

  while (tokens != NULL) {
    switch (tokenNumber) {
      case 0:
        angle.atom1 = atoi(tokens);
        break;
      case 1:
        angle.atom2 = atoi(tokens);
        break;
      case 2:
        angle.value = atof(tokens);
        break;
      case 3:
        if (atoi(tokens) == 1)
          angle.variable = true;
        else
          angle.variable = false;
    }

    tokens = strtok(NULL, " ");
    tokenNumber++;
  }

  free (charLine);

  return angle;
}

Atom StateScanner::getAtomFromLine(string line) {
	unsigned long inval = (unsigned long) -1;
  Atom atom = createAtom(inval, -1,-1,-1,-1,-1);

  char *tokens;
  char *charLine = (char *) malloc(sizeof(char) * (line.size()+1));
  strcpy(charLine, line.c_str());

  tokens = strtok(charLine, " ");
  int tokenNumber = 0;

  while (tokens != NULL) {
    switch (tokenNumber) {
      case 0:
        atom.id = atoi(tokens);
        break;
      case 1:
        atom.x = atof(tokens);
        break;
      case 2:
        atom.y = atof(tokens);
        break;
      case 3:
        atom.z = atof(tokens);
        break;
      case 4:
        atom.sigma = atof(tokens);
        break;
      case 5:
        atom.epsilon = atof(tokens);
        break;
      case 6:
        atom.charge = atof(tokens);
        break;
      case 7:
	      *atom.name = std::string(tokens);
 	      break;
    }

    tokens = strtok(NULL, " ");
    tokenNumber++;
  }

  free (charLine);

  return atom;
}

Bond StateScanner::getBondFromLine(string line) {
  Bond bond = Bond();
  char *tokens;

  char *charLine = (char *) malloc(sizeof(char) * (line.size()+1));
  strcpy(charLine, line.c_str());

  tokens = strtok(charLine, " ");
  int tokenNumber = 0;

  while (tokens != NULL) {
    switch (tokenNumber) {
      case 0: // atom 1
        bond.atom1 = atoi(tokens);
        break;
      case 1: // atom 2
        bond.atom2 = atoi(tokens);
        break;
      case 2: // value
        bond.distance = atof(tokens);
        break;
      case 3: // variable
        if (atoi(tokens) == 1)
          bond.variable = true;
        else
          bond.variable = false;
        break;
    }

    tokens = strtok(NULL, " ");
    tokenNumber++;
  }

  free (charLine);

  return bond;
}

Dihedral StateScanner::getDihedralFromLine(string line) {
  Dihedral dihedral = Dihedral();

  char *tokens;
  char *charLine = (char *) malloc(sizeof(char) * (line.size()+1));
  strcpy(charLine, line.c_str());

  tokens = strtok(charLine, " ");
  int tokenNumber = 0;

  while (tokens != NULL) {
    switch (tokenNumber) {
      case 0:
        dihedral.atom1 = atoi(tokens);
        break;
      case 1:
        dihedral.atom2 = atoi(tokens);
        break;
      case 2:
        dihedral.value = atof(tokens);
        break;
      case 3:
        if (atoi(tokens) == 1)
          dihedral.variable = true;
        else
          dihedral.variable = false;
        break;
    }

    tokens = strtok(NULL, " ");
    tokenNumber++;
  }

  free (charLine);

  return dihedral;
}

Hop StateScanner::getHopFromLine(string line) {
  Hop hop = Hop();

  char *tokens;
  char *charLine = (char *) malloc(sizeof(char) * (line.size()+1));
  strcpy(charLine, line.c_str());

  tokens = strtok(charLine, " ");
  int tokenNumber = 0;

  while (tokens != NULL) {
    switch (tokenNumber) {
      case 0:
        hop.atom1 = atoi(tokens);
        break;
      case 1:
        hop.atom2 = atoi(tokens);
        break;
      case 2:
        hop.hop = atoi(tokens);
        break;
    }

    tokens = strtok(NULL, " ");
    tokenNumber++;
  }

  free (charLine);

  return hop;
}

Environment* StateScanner::getEnvironmentFromLine(string line) {
  Environment* environment = new Environment();

  char *tokens;
  char *charLine = (char *) malloc(sizeof(char) * (line.size()+1));
  strcpy(charLine, line.c_str());
  tokens = strtok(charLine, " ");

  int tokenNumber = 0;

  while (tokens != NULL) {
    switch (tokenNumber) {
      case 0:
        environment->x = atof(tokens);
        break;
      case 1:
        environment->y = atof(tokens);
        break;
      case 2:
        environment->z = atof(tokens);
        break;
      case 3:
        environment->numOfMolecules = atoi(tokens);
        break;
      case 4:
        environment->numOfAtoms = atoi(tokens);
        break;
      case 5:
        environment->temp = atof(tokens);
        break;
      case 6:
        environment->cutoff = atof(tokens);
        break;
      case 7:
        environment->maxTranslation = atof(tokens);
        break;
      case 8:
        environment->maxRotation = atof(tokens);
        break;
      case 9:
        parsePrimaryIndexDefinitions(environment, std::string(tokens));
		    environment->randomseed =  atoi(tokens + (strlen(tokens) + 1));
        break;
      case 10:
        environment->randomseed = atoi(tokens);
        break;
    }

    tokens = strtok(NULL, " ");
    tokenNumber++;
  }

  free (charLine);

  return environment;
}

void StateScanner::parsePrimaryIndexDefinitions(Environment* enviro, string definitions) {
  *enviro->primaryAtomIndexConfigLine = definitions;

  char *indexVector;
  char *charLine = (char *) malloc(sizeof(char) * (definitions.size()+1));
  strcpy(charLine, definitions.c_str());

  indexVector = strtok(charLine, ",");

  for (int i = 0; indexVector; i++) {
  	std::string sIndexVector = indexVector;
  	enviro->primaryAtomIndexDefinitions++;
  	char* c;

  	if ((c=strchr(indexVector, '[')) != NULL) {
    	indexVector = (indexVector + 1);

      while ((c=strchr(indexVector, ']')) == NULL) {
        int currentPrimaryIndex = atoi(indexVector);

        (*(enviro->primaryAtomIndexArray)).push_back(new std::vector<int>);
        (*(*(enviro->primaryAtomIndexArray))[i]).push_back(currentPrimaryIndex - 1);
        indexVector = strtok(NULL, ",");
    	}

      *c = '\0';
    	int currentPrimaryIndex = atoi(indexVector);

      (*(enviro->primaryAtomIndexArray)).push_back(new std::vector<int>);
    	(*(*(enviro->primaryAtomIndexArray))[i]).push_back(currentPrimaryIndex - 1);
    	indexVector =strtok(NULL, ",");

      continue;
    }

    int currentPrimaryIndex = atoi(indexVector);

    (*(enviro->primaryAtomIndexArray)).push_back(new std::vector<int>);
    (*(*(enviro->primaryAtomIndexArray))[i]).push_back(currentPrimaryIndex - 1);
    indexVector = strtok(NULL, ",");
  }

  free (charLine);
}



void StateScanner::outputState(Environment *environment, Molecule *molecules, int numOfMolecules, int step, string filename) {
  ofstream outFile;
  outFile.open(filename.c_str());
  outFile << environment->x << " " << environment->y << " "
          << environment->z << " " << environment->numOfMolecules << " "
          << environment->numOfAtoms << " " << environment->temp << " "
          << environment->cutoff << " " << environment->maxTranslation << " "
          << environment->maxRotation << " "
          << *environment->primaryAtomIndexConfigLine << " "
          << environment->randomseed << " "
          << std::endl << step << std::endl << std::endl;

  for (int i = 0; i < numOfMolecules; i++) {
    Molecule currentMol = molecules[i];
    outFile << currentMol.id << std::endl;

    outFile << "= Type" << std::endl;
    outFile << currentMol.type << std::endl;

    // Write atoms
    outFile << "= Atoms" << std::endl;

    for (int j = 0; j < currentMol.numOfAtoms; j++) {
      Atom currentAtom = currentMol.atoms[j];

      outFile << currentAtom.id << " "
              << currentAtom.x << " " << currentAtom.y
              << " " << currentAtom.z << " "
              << currentAtom.sigma << " "
              << currentAtom.epsilon << " "
              << currentAtom.charge << " "
		          << *currentAtom.name << " " << std::endl;
    }

    // Write bonds
    outFile << "= Bonds" << std::endl;

    for (int j = 0; j < currentMol.numOfBonds; j++) {
      Bond currentBond = currentMol.bonds[j];
      outFile << currentBond.atom1 << " "
              << currentBond.atom2 << " "
              << currentBond.distance << " ";

      if (currentBond.variable) {
        outFile << "1" << std::endl;
      } else {
        outFile << "0" << std:: endl;
      }
    }

    // Write dihedrals
    outFile << "= Dihedrals" << std::endl;

    for (int j = 0; j < currentMol.numOfDihedrals; j++) {
      Dihedral currentDi = currentMol.dihedrals[j];
      outFile << currentDi.atom1 << " "
              << currentDi.atom2 << " "
              << currentDi.value << " ";

      if (currentDi.variable) {
        outFile << "1" << std::endl;
      } else {
        outFile << "0" << std::endl;
      }
    }

    // Write hops
    outFile << "= Hops" << std::endl;

    for (int j = 0; j < currentMol.numOfHops; j++) {
      Hop currentHop = currentMol.hops[j];

      outFile << currentHop.atom1 << " "
              << currentHop.atom2 << " "
              << currentHop.hop << std::endl;
    }

    // Write angles
    outFile << "= Angles" << std::endl;

    for (int j = 0; j < currentMol.numOfAngles; j++) {
      Angle currentAngle = currentMol.angles[j];

      outFile << currentAngle.atom1 << " "
              << currentAngle.atom2 << " "
              << currentAngle.value << " ";

      if(currentAngle.variable) {
        outFile << "1" << std::endl;
      } else {
        outFile << "0" << std::endl;
      }
    }

    outFile << "==" << std::endl;
  }

  outFile.close();
}
