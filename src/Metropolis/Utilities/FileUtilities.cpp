#include "FileUtilities.h"

#include "Parsing.h"
#include "StructLibrary.h"
#include "Metropolis/Box.h"
#include "Metropolis/SimulationArgs.h"
#include "MathLibrary.h"

#include <exception>
#include <stdexcept>
#include <time.h>


using std::string;
using std::ifstream;

#define DEFAULT_STEP_COUNT 100

bool loadBoxData(SimulationArgs& simArgs, Box* box, long* startStep,
                 long* steps, SBScanner* sbScanner) {
  if (box == NULL) { // If box is null, print an error and return.
    std::cerr << "Error: loadBoxData(): Box is NULL" << std::endl;
    return false;
  }

  Environment* enviro;
  vector<Molecule> moleculeVector;

  // Build from config/z-matrix.
  if (simArgs.fileType == InputFile::Configuration) {
    // Read in config information from config file.
    ConfigScanner config_scanner = ConfigScanner();
    if (!config_scanner.readInConfig(simArgs.filePath)) {
      std::cerr << "Error: loadBoxData(): Could not read config file"
                << std::endl;
      return false;
    }
    simArgs.pdbOutputPath = config_scanner.getPdbOutputPath();
    simArgs.stateOutputPath = config_scanner.getStateOutputPath();
    if (!config_scanner.getSimulationName().empty() &&
        simArgs.simulationName.empty()) {
      simArgs.simulationName = config_scanner.getSimulationName();
    }

    if (!config_scanner.getStrategy().empty() &&
        simArgs.strategy == Strategy::Default) {
      simArgs.strategy = Strategy::fromString(config_scanner.getStrategy());
      if (simArgs.strategy == Strategy::Unknown) {
        std::cerr << "Error: Unknown energy calculation strategy specified "
                     "in config file"
                  << std::endl;
        return false;
      }
    }

    // Getting bond and angle data from oplsaa.sb file.
    std::string sb_path = config_scanner.getOplsusaparPath();
    std::size_t slash_index= sb_path.find_last_of('/');
    sb_path = sb_path.substr(0, slash_index);

    if(!sbScanner->readInSB(sb_path + "/oplsaa.sb")) {
      std::cerr << "Error: loadBoxData(): Could not read OPLS SB file"
                << std::endl;
      return false;
    }

    // Getting atom data from the oplsaa.par / oplsua.par files.
    OplsScanner opls_scanner = OplsScanner();
    if (!opls_scanner.readInOpls(config_scanner.getOplsusaparPath())) {
      std::cerr << "Error: loadBoxData(): Could not read OPLS file"
                << std::endl;
      return false;
    }

    // Reading in molecule information from the z-matrix.
    ZmatrixScanner zmatrix_scanner = ZmatrixScanner();
    if (!zmatrix_scanner.readInZmatrix(config_scanner.getZmatrixPath(),
                                       &opls_scanner)) {
      std::cerr << "Error: loadBoxData(): Could not read Z-Matrix file"
                << std::endl;
      return false;
    }

    // Initialize moleculeVector, environment and steps.
    moleculeVector = zmatrix_scanner.buildMolecule(0);
    enviro = config_scanner.getEnviro();
    *steps = config_scanner.getSteps();
    *startStep = 0;

    // Verify that the # of molecules is consitent in z-matrix and config
    if (moleculeVector.size() != enviro->primaryAtomIndexDefinitions) {
      std::cerr << "Error: loadBoxData(): The number of molecules read "
                << "from the Z Matrix file (" << moleculeVector.size()
                << ") does not equal the number of primary index "
                << "definitions in the config file ("
                << enviro->primaryAtomIndexDefinitions << ")"
                << std::endl;
      return false;
    }

    // Instantiate the box's environment.
    box->environment = new Environment(enviro);

    // Build the box's data.
    if (!buildBoxData(enviro, moleculeVector, box, *sbScanner)) {
      std::cerr << "Error: loadBoxData(): Could not build box data" << std::endl;
      return false;
    }

    return true;
  }

  // Build from state file.
  if (simArgs.fileType == InputFile::State) {
    // Instantiate a state scanner and read in the environment.
    StateScanner state_scanner = StateScanner(simArgs.filePath);
    enviro = state_scanner.readInEnvironment();

    // Validate the environment
    if (enviro == NULL) {
      std::cerr << "Error: Unable to read environment from State File"
                << std::endl;
      return false;
    }

    // Get the molecule vector from the state scanner.
    moleculeVector = state_scanner.readInMolecules();

    // Validate the molecule vector.
    if (moleculeVector.size() == 0) {
      std::cerr << "Error: Unable to read molecule data from State file"
                << std::endl;
      return false;
    }

    // Instantiate steps / start step.
    *steps = DEFAULT_STEP_COUNT;
    *startStep = state_scanner.readInStepNumber();

    // Validate start step.
    if (*startStep < 0) {
      std::cerr << "Error: State file does not have a starting step number"
                << std::endl;
      return false;
    }

    // If environment is validated, set the box's environment.
    box->environment = new Environment(enviro);

    // Populate the box with data.
    if (!fillBoxData(enviro, moleculeVector, box, *sbScanner)) {
      std::cerr << "Error: loadBoxData(): Could not build box data"
                << std::endl;
      return false;
    }

    return true;
  }

  // If we are neither building from config/z-matrix or state file, print error.
  std::cout << "Error: Could not recognize input file type" << std::endl;
  return false;
}


bool fillBoxData(Environment* enviro, vector<Molecule>& molecVec, Box* box,
                 SBScanner& sbScanner) {
  // If the vector of molecules has no contents, print an error and return
  if (!enviro || !box || molecVec.size() < 1) {
    std::cerr << "Error: fillBoxData(): Could not fill molecule data."
              << std::endl;
    return false;
  }

  // Holds the running number of atoms, bonds, angles, dihedrals, and hops.
  int count[5];
  memset(count,0,sizeof(count));  // Instantiate count.

  // Instantiate the box's relevant variables.
  box->moleculeCount = 0;
  box->atomCount = 0;
  box->bondCount = 0;
  box->angleCount = 0;
  box->dihedralCount = 0;
  box->hopCount = 0;

  // For every molecule in the vector, increment the box's count for each
  // sub-value accordingly.
  for (int i = 0; i < molecVec.size(); ++i) {
    Molecule mol = molecVec[i];
    box->moleculeCount += 1;
    box->atomCount += mol.numOfAtoms;
    box->bondCount += mol.numOfBonds;
    box->angleCount += mol.numOfAngles;
    box->dihedralCount += mol.numOfDihedrals;
    box->hopCount += mol.numOfHops;
  }

  // Allocate space for the box's various structures.
  box->molecules = (Molecule *) malloc(sizeof(Molecule) * box->moleculeCount);
  box->atoms     = (Atom *) malloc(sizeof(Atom) * box->atomCount);
  box->bonds     = (Bond *) malloc(sizeof(Bond) *  box->bondCount);
  box->angles    = (Angle *) malloc(sizeof(Angle) * box->angleCount);
  box->dihedrals = (Dihedral *) malloc(sizeof(Dihedral) * box->dihedralCount);
  box->hops      = (Hop *) malloc(sizeof(Hop) * box->hopCount);

  // Clear the space allocated for the box's structures.
  // Why do we not clear box->molecules?
  memset(box->atoms, 0, sizeof(Atom) * box->atomCount);
  memset(box->bonds, 0, sizeof(Bond) * box->bondCount);
  memset(box->angles, 0, sizeof(Angle) * box->angleCount);
  memset(box->dihedrals, 0, sizeof(Dihedral) * box->dihedralCount);
  memset(box->hops, 0, sizeof(Hop) * box->hopCount);

  std::vector<Bond> bondVector; // bondVector will hold a list of bonds.
  // These will be used to get the common atom in every angle.

  //Copy data from vector to molecule
  for (int j = 0; j < molecVec.size(); j++) {
    Molecule molec1 = molecVec[j];

    // Point to memory allocated to atoms, bonds, etc in the box with the
    // atoms, bonds, etc. arrays in each molecule.
    box->molecules[j].atoms = (Atom *) (box->atoms + count[0]);
    box->molecules[j].bonds = (Bond *) (box->bonds + count[1]);
    box->molecules[j].angles = (Angle *) (box->angles + count[2]);
    box->molecules[j].dihedrals = (Dihedral *) (box->dihedrals + count[3]);
    box->molecules[j].hops = (Hop *) (box->hops + count[4]);

    // Copy over constants from each element in the molecule vector.
    box->molecules[j].id = molec1.id;
    box->molecules[j].type = molec1.type;
    box->molecules[j].numOfAtoms = molec1.numOfAtoms;
    box->molecules[j].numOfBonds = molec1.numOfBonds;
    box->molecules[j].numOfDihedrals = molec1.numOfDihedrals;
    box->molecules[j].numOfAngles = molec1.numOfAngles;
    box->molecules[j].numOfHops = molec1.numOfHops;

    // Increment the offset to point data to.
    count[0] += molec1.numOfAtoms;
    count[1] += molec1.numOfBonds;
    count[2] += molec1.numOfAngles;
    count[3] += molec1.numOfDihedrals;
    count[4] += molec1.numOfHops;

    // Copy over atoms from the molecule vector.
    for (int k = 0; k < molec1.numOfAtoms; k++) {
      box->molecules[j].atoms[k] = molec1.atoms[k];
    }

    // Copy over bonds from the molecule vector.
    for (int k = 0; k < molec1.numOfBonds; k++) {
      box->molecules[j].bonds[k] = molec1.bonds[k];

      // Add the bond to bondVector
      bondVector.push_back(molec1.bonds[k]);
    }

    // Copy over angles from the molecule vector.
    for (int k = 0; k < molec1.numOfAngles; k++) {
      box->molecules[j].angles[k] = molec1.angles[k];
    }

    // Copy over dihedrals from the molecule vector.
    for(int k = 0; k < molec1.numOfDihedrals; k++) {
      box->molecules[j].dihedrals[k] = molec1.dihedrals[k];
    }

    // Copy over hops from the molecule vector.
    for(int k = 0; k < molec1.numOfHops; k++) {
      box->molecules[j].hops[k] = molec1.hops[k];
    }

  } // Repeat for every molecule in the vector.

  enviro->numOfAtoms = box->atomCount;
  return true;
}

bool buildBoxData(Environment* enviro, vector<Molecule>& molecVec, Box* box,
                  SBScanner& sbScanner) {

  // Convert molecule vectors into an array
  //int moleculeIndex = 0;
  //int atomCount = 0;

  // If the vector of molecules has no contents, print an error and return.
  if (!enviro || !box || molecVec.size() < 1) {
    std::cerr << "Error: buildBoxData(): Could not load molecule data."
              << std::endl;
    return false;
  }

  // Make nunOfMolecules divisible by the size of the molecule vector.
  int molecMod = enviro->numOfMolecules % molecVec.size();
  if (molecMod != 0)  {
    enviro->numOfMolecules += molecVec.size() - molecMod;
  }

  // Allocate memory for the molecues in the box.
  box->moleculeCount = enviro->numOfMolecules;
  box->molecules = (Molecule *) malloc(sizeof(Molecule) * box->moleculeCount);

  // If there aren't any molecules, return an error.
  if(box->moleculeCount < 1 || molecVec.size() < 1) {
    std::cerr << "The simulation environment should have at least 1 "
              << "molecule. (The current count? " << box->moleculeCount
              <<")." << std::endl;
    return false;
  }


  int molecDiv = enviro->numOfMolecules / molecVec.size();
  int molecTypeNum = molecVec.size();

  // Holds the running number of atoms, bonds, angles, dihedrals, and hops.
  int count[5];
  memset(count,0,sizeof(count)); // Instantiate count.

  // Allocate memory for tables.
  Table* tables = new Table[molecVec.size()];
  int currentAtomCount = 0;

  //Copy data from vector to molecule
  for (int j = 0; j < molecVec.size(); j++) {
    Molecule molec1 = molecVec[j];

    // Count the number of atoms, bonds, etc.
    count[0] += molec1.numOfAtoms;
    count[1] += molec1.numOfBonds;
    count[2] += molec1.numOfAngles;
    count[3] += molec1.numOfDihedrals;
    count[4] += molec1.numOfHops;

    // Begin copying hops.
    Hop *myHop = molec1.hops;
    int **table;

    // Allocate memory to tables.
    table = new int*[molec1.numOfAtoms];
    for(int k = 0; k < molec1.numOfAtoms;k++) {
      table[k] = new int[molec1.numOfAtoms];
    }

    // Instantiate table to 0.
    for(int test = 0; test< molec1.numOfAtoms;test++) {
      for(int test1 = 0; test1 < molec1.numOfAtoms; test1++) {
          table[test][test1] = 0;
      }
    }

    // Fill table with hops.
    for(int k2 = 0; k2<molec1.numOfHops;k2++) {
      int atom1 = myHop->atom1;
      int atom2 = myHop->atom2;
      table[atom1-currentAtomCount][atom2-currentAtomCount] = myHop->hop;
      table[atom2-currentAtomCount][atom1-currentAtomCount] = myHop->hop;
      myHop++;
    }

    // Put the table in tables.
    tables[j] = Table(table);
    // Increase the number atoms in the molecule.
    currentAtomCount += molec1.numOfAtoms;
  }

  // Fill box variables.
  box->atomCount = molecDiv * count[0];
  box->bondCount = molecDiv * count[1];
  box->angleCount = molecDiv * count[2];
  box->dihedralCount = molecDiv * count[3];
  box->hopCount = molecDiv * count[4];

  // Allocate memory to box.
  box->atoms     = (Atom *) malloc(sizeof(Atom) * box->atomCount);
  box->bonds     = (Bond *) malloc(sizeof(Bond) * box->bondCount);
  box->angles    = (Angle *) malloc(sizeof(Angle) * box->angleCount);
  box->dihedrals = (Dihedral *) malloc(sizeof(Dihedral) * box->dihedralCount);
  box->hops      = (Hop *) malloc(sizeof(Hop) * box->hopCount);

  // Clear memory in the box.
  memset(box->atoms, 0, sizeof(Atom) * box->atomCount);
  memset(box->bonds, 0, sizeof(Bond) * box->bondCount);
  memset(box->angles, 0, sizeof(Angle) * box->angleCount);
  memset(box->dihedrals, 0, sizeof(Dihedral) * box->dihedralCount);
  memset(box->hops, 0, sizeof(Hop) * box->hopCount);

  //Clear count.
  memset(count, 0, sizeof(count));

  std::vector<Bond> bondVector; // bondVector will hold a list of bonds.
  // These will be used to get the common atom in every angle.

  //Copy data from vector to molecule.
  for(int j = 0; j < molecVec.size(); j++) {

    Molecule molec1 = molecVec[j]; // Get the molecule from the molecule vector.

    // Point to memory allocated to atoms, bonds, etc in the box with the
    // atoms, bonds, etc. arrays in each molecule.
    box->molecules[j].atoms = (Atom *) (box->atoms + count[0]);
    box->molecules[j].bonds = (Bond *) (box->bonds + count[1]);
    box->molecules[j].angles = (Angle *) (box->angles + count[2]);
    box->molecules[j].dihedrals = (Dihedral *) (box->dihedrals + count[3]);
    box->molecules[j].hops = (Hop *) (box->hops + count[4]);

    // Copy over molecule values from the molecule in the vector.
    box->molecules[j].id = molec1.id;
    box->molecules[j].type = molec1.type;
    box->molecules[j].numOfAtoms = molec1.numOfAtoms;
    box->molecules[j].numOfBonds = molec1.numOfBonds;
    box->molecules[j].numOfDihedrals = molec1.numOfDihedrals;
    box->molecules[j].numOfAngles = molec1.numOfAngles;
    box->molecules[j].numOfHops = molec1.numOfHops;

    // Increment the count of atoms, bonds, etc.
    // This will set the offset for the next pointer assignment.
    count[0] += molec1.numOfAtoms;
    count[1] += molec1.numOfBonds;
    count[2] += molec1.numOfAngles;
    count[3] += molec1.numOfDihedrals;
    count[4] += molec1.numOfHops;

    // Copy atoms from the molecule vector.
    for (int k = 0; k < molec1.numOfAtoms; k++) {
      box->molecules[j].atoms[k] = molec1.atoms[k];
    }

    // Copy bonds from the molecule vector.
    for (int k = 0; k < molec1.numOfBonds; k++) {
      box->molecules[j].bonds[k] = molec1.bonds[k];

      // Add the bond to bondVector
      bondVector.push_back(molec1.bonds[k]);
    }

    // Copy angles from the molecule vector.
    for (int k = 0; k < molec1.numOfAngles; k++) {
      box->molecules[j].angles[k] = molec1.angles[k];
      int a1Idx = molec1.angles[k].atom1;
      int a2Idx = molec1.angles[k].atom2;
      if (box->molecules[j].angles[k].commonAtom != 0) {
        box->molecules[j].angles[k].commonAtom = getCommonAtom(bondVector,
                                                               a1Idx, a2Idx);
      }
    }


    // Copy dihedrals from the molecule vector.
    for (int k = 0; k < molec1.numOfDihedrals; k++) {
      box->molecules[j].dihedrals[k] = molec1.dihedrals[k];
    }

    // Copy hops from the molecule vector.
    for (int k = 0; k < molec1.numOfHops; k++) {
      box->molecules[j].hops[k] = molec1.hops[k];
    }

  } // Repeat for every molecule in the molecule vector.

  // Repeat for each copy of every molecule in the vector.
  for (int m = 1; m < molecDiv; m++) {
    // Copy the molecule into the corresponding spot in the box.
    int offset = m * molecTypeNum;
    memcpy(&(box->molecules[offset]), box->molecules,
           sizeof(Molecule) * molecTypeNum);

    // Copy atoms, bonds, etc. for each molecule.
    for (int n = 0; n < molecTypeNum; n++) {
      box->molecules[offset+n].id=offset+n;
      box->molecules[offset+n].atoms = box->molecules[n].atoms + count[0] * m;
      box->molecules[offset+n].bonds = box->molecules[n].bonds + count[1] * m;
      box->molecules[offset+n].angles = box->molecules[n].angles 
        + count[2] * m;
      box->molecules[offset+n].dihedrals = box->molecules[n].dihedrals
        + count[3] * m;
      box->molecules[offset+n].hops =  box->molecules[n].hops + count[4] * m;
    }

    memcpy(&(box->atoms[m * count[0]]), box->atoms, sizeof(Atom) * count[0]);
    memcpy(&(box->bonds[m * count[1]]), box->bonds, sizeof(Bond) * count[1]);
    memcpy(&(box->angles[m * count[2]]),
           box->angles, sizeof(Angle) * count[2]);
    memcpy(&(box->dihedrals[m * count[3]]), box->dihedrals,
           sizeof(Dihedral) * count[3]);
    memcpy(&(box->hops[m * count[4]]), box->hops, sizeof(Hop) * count[4]);

    // Copy atoms for each molecule.
    for (int k = 0; k < count[0]; k++)  {
      box->atoms[m * count[0] + k].id = m * count[0] + k;
    }

    // Copy bonds for each molecule.
    for (int k = 0; k < count[1]; k++) {
      box->bonds[m * count[1] + k].atom1 += m * count[0];
      box->bonds[m * count[1] + k].atom2 += m * count[0];
    }

    // Copy angles for each molecule.
    for(int k = 0; k < count[2]; k++) {
      box->angles[m * count[2] + k].atom1 += m * count[0];
      box->angles[m * count[2] + k].atom2 += m * count[0];
      box->angles[m * count[2] + k].commonAtom += m * count[0];
    }

    // Copy dihedrals for each molecule.
    for(int k = 0; k < count[3]; k++) {
      box->dihedrals[m * count[3] + k].atom1 += m * count[0];
      box->dihedrals[m * count[3] + k].atom2 += m * count[0];
    }

    // Copy hops for each molecule.
    for(int k = 0; k < count[4]; k++) {
      box->hops[m * count[4] + k].atom1 += m * count[0];
      box->hops[m * count[4] + k].atom2 += m * count[0];
    }

  } // Repeat for every molecule.

  enviro->numOfAtoms = count[0] * molecDiv;

  //generate fcc lattice box
  if (!generatefccBox(box)) {
    std::cerr << "Error: buildBoxData(): Could not generate FCC box"
              << std::endl;
    return false;
  }

  // Now that we have XYZ coordinates, compute the values of any implicit
  // angles.
  for (int angle = 0; angle < box->angleCount; angle++) {
    if (box->angles[angle].value == 0) {
      // Calculate the correct angle value base on the XYZ coordinates
      Angle a = box->angles[angle];
      Atom *atoms = box->atoms;
      box->angles[angle].value = getAngle(atoms[a.atom1], atoms[a.commonAtom],
                                          atoms[a.atom2]);
    }
  }


  return true;
}

bool generatefccBox(Box* box) {
  // Return an error if the box is missing an environment or molecules.
  if (box->environment == NULL || box->molecules == NULL) {
    std::cerr << "Error: generatefccBox(): Box environment or molecules "
                 "array NULL" << std::endl;
    return false;
  }

  Real cells, dcells, cellL, halfcellL;
  Environment* enviro = box->environment;
  Molecule* molecules = box->molecules;

  // Determine the number of unit cells in each coordinate direction
  dcells = pow(0.25 * (Real) enviro->numOfMolecules, 1.0/3.0);
  cells = (int)(dcells + 0.5);

  // Check if numOfMolecules is a non-fcc number of molecules and increase
  // the number of cells if necessary
  while((4 * cells * cells * cells) < enviro->numOfMolecules) {
    cells++;
  }

  //Determine length of unit cell
  cellL = enviro->x / (Real) cells;
  halfcellL = 0.5 * cellL;

  //Construct the unit cell
  for (int j = 0; j < molecules[0].numOfAtoms; j++) {
    molecules[0].atoms[j].x += 0.0;
    molecules[0].atoms[j].y += 0.0;
    molecules[0].atoms[j].z += 0.0;
  }

  for (int j = 0; j < molecules[1].numOfAtoms; j++) {
    molecules[1].atoms[j].x += halfcellL;
    molecules[1].atoms[j].y += halfcellL;
    molecules[1].atoms[j].z += 0.0;
  }

  for (int j = 0; j < molecules[2].numOfAtoms; j++) {
    molecules[2].atoms[j].x += 0.0;
    molecules[2].atoms[j].y += halfcellL;
    molecules[2].atoms[j].z += halfcellL;
  }

  for (int j = 0; j < molecules[3].numOfAtoms; j++) {
    molecules[3].atoms[j].x += halfcellL;
    molecules[3].atoms[j].y += 0.0;
    molecules[3].atoms[j].z += halfcellL;
  }

  // Init all other molecules to initial coordinates
  // Build the lattice from the unit cell by repeatedly translating
  // the four vectors of the unit cell through a distance cellL in
  // the x, y, and z directions
  for(int i = 4; i < enviro->numOfMolecules; i++) {
    for (int j = 0; j < molecules[i].numOfAtoms; j++) {
      molecules[i].atoms[j].x += 0.0;
      molecules[i].atoms[j].y += 0.0;
      molecules[i].atoms[j].z += 0.0;
    }
  }

  int offset = 0;
  for (int z = 1; z <= cells; z++) {
    for (int y = 1; y <= cells; y++) {
      for (int x = 1; x <= cells; x++) {
        for (int a = 0; a < 4; a++) {
          int i = a + offset;
          if (i < enviro->numOfMolecules) {
            for (int j = 0; j < molecules[i].numOfAtoms; j++) {
              molecules[i].atoms[j].x = molecules[a].atoms[j].x + cellL * (x-1);
              molecules[i].atoms[j].y = molecules[a].atoms[j].y + cellL * (y-1);
              molecules[i].atoms[j].z = molecules[a].atoms[j].z + cellL * (z-1);
            }
          }
        }
        offset += 4;
      }
    }
  }

  //Shift center of box to the origin
  for (int i = 0; i < enviro->numOfMolecules; i++) {
    for (int j = 0; j < molecules[i].numOfAtoms; j++) {
      molecules[i].atoms[j].x -= halfcellL;
      molecules[i].atoms[j].y -= halfcellL;
      molecules[i].atoms[j].z -= halfcellL;
    }
  }

  return true;
}

// ============================================================================
// ======================= Logging Functions ==================================
// ============================================================================

void writeToLog(string text,int stamp) {
  string filename = "OutputLog";
  ofstream logFile;
  logFile.open(filename.c_str(),ios::out|ios::app);

  string hash ="";
  time_t current_time;
  struct tm * time_info;
  char timeString[9];  // space for "HH:MM:SS\0"

  switch(stamp) {
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

void writeToLog(stringstream& ss, int stamp) {
  writeToLog(ss.str(),stamp);
  ss.str(""); // clears the string steam...
  ss.clear();
}

std::string& rtrim(std::string& s) {
  s.erase(s.find_last_not_of(TRIMMED_CHARS) + 1);
  return s;
}

std::string& ltrim(std::string& s) {
  std::size_t len = s.find_first_not_of(TRIMMED_CHARS);
  if (len != std::string::npos)
    s.erase(0, len);
  return s;
}

std::string& trim(std::string& s) {
  return ltrim(rtrim(s));
}
