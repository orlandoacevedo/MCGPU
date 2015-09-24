#include "FileUtilities.h"

//#include <iostream>
//#include <fstream>
//#include <stdlib.h>
#include <exception>
#include <stdexcept>
//#include <sstream>
#include "Parsing.h"
#include "StructLibrary.h"
#include "Metropolis/Box.h"
#include "Metropolis/SimulationArgs.h"

using std::string;
using std::ifstream;

#define DEFAULT_STEP_COUNT 100


bool loadBoxData(string inputPath, InputFileType inputType, Box* box, long* startStep, long* steps) {

	if (box == NULL) {			// If box is null, print an error and return.
		std::cerr << "Error: loadBoxData(): Box is NULL" << std::endl;
		return false;
	}

  Environment* enviro;
  vector<Molecule> moleculeVector;

  if (inputType == InputFile::Configuration) {	// Build from config/z-matrix.

		// Reading in config information from config file.
		ConfigScanner config_scanner = ConfigScanner();
    if (!config_scanner.readInConfig(inputPath)) {
  		std::cerr << "Error: loadBoxData(): Could not read config file" << std::endl;
      return false;
    }

		// Getting bond and angle data from oplsaa.sb file.
		SBScanner sb_scanner = SBScanner();
		std::string sb_path = config_scanner.getOplsusaparPath();
		std::size_t slash_index= sb_path.find_last_of('/');
		sb_path = sb_path.substr(0, slash_index);

		if(!sb_scanner.readInSB(sb_path + "/oplsaa.sb")) {
			std::cerr << "Error: loadBoxData(): Could not read OPLS SB file" << std::endl;
			return false;
		}

		// Getting atom data from the oplsaa.par / oplsua.par files.
    OplsScanner opls_scanner = OplsScanner();
    if (!opls_scanner.readInOpls(config_scanner.getOplsusaparPath())) {
    	std::cerr << "Error: loadBoxData(): Could not read OPLS file" << std::endl;
      return false;
    }

		// Reading in molecule information from the z-matrix.
    ZmatrixScanner zmatrix_scanner = ZmatrixScanner();
    if (!zmatrix_scanner.readInZmatrix(config_scanner.getZmatrixPath(), &opls_scanner)) {
    	std::cerr << "Error: loadBoxData(): Could not read Z-Matrix file" << std::endl;
      return false;
    }

		// Initialize moleculeVector, environment and steps.
    moleculeVector = zmatrix_scanner.buildMolecule(0);
    enviro = config_scanner.getEnviro();
    *steps = config_scanner.getSteps();
    *startStep = 0;

		// Verify that the # of molecules is consitent in z-matrix and config files.
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
    if (!buildBoxData(enviro, moleculeVector, box, sb_scanner)) {
    	std::cerr << "Error: loadBoxData(): Could not build box data" << std::endl;
      return false;
    }

  	return true;

  } else if (inputType == InputFile::State) {		// Build from state file.

		// Instantiate a state scanner and read in the environment.
		StateScanner state_scanner = StateScanner(inputPath);
    enviro = state_scanner.readInEnvironment();

		// Validate the environment
    if (enviro == NULL) {
    	std::cerr << "Error: Unable to read environment from State File" << std::endl;
      return false;
    }

		// Get the molecule vector from the state scanner.
		moleculeVector = state_scanner.readInMolecules();

		// Validate the molecule vector.
    if (moleculeVector.size() == 0) {
    	std::cerr << "Error: Unable to read molecule data from State file" << std::endl;
      return false;
    }

		// Instantiate steps / start step.
    *steps = DEFAULT_STEP_COUNT;
    *startStep = state_scanner.readInStepNumber();

		// Validate start step.
    if (*startStep < 0) {
  		std::cerr << "Error: State file does not have a starting step number" << std::endl;
      return false;
  	}

		// If environment is validated, set the box's environment.
    box->environment = new Environment(enviro);

		// Populate the box with data.
    if (!fillBoxData(enviro, moleculeVector, box)) {
    	std::cerr << "Error: loadBoxData(): Could not build box data" << std::endl;
      return false;
    }

    return true;
  }

	// If we are neither building from config/z-matrix or state file, print error.
  std::cout << "Error: Could not recognize input file type" << std::endl;
	return false;
}


bool fillBoxData(Environment* enviro, vector<Molecule>& molecVec, Box* box) {

	//if the vector of molecules has no contents, print an error and return.
	if (!enviro || !box || molecVec.size() < 1) {
  	std::cerr << "Error: fillBoxData(): Could not fill molecule data." << std::endl;
    return false;
  }

	// Holds the running number of atoms, bonds, angles, dihedrals, and hops.
  int count[5];
  memset(count,0,sizeof(count));	// Instantiate count.

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

bool buildBoxData(Environment* enviro, vector<Molecule>& molecVec, Box* box, SBScanner& sb_scanner) {

	// Convert molecule vectors into an array
  int moleculeIndex = 0;
  int atomCount = 0;

	// If the vector of molecules has no contents, print an error and return.
  if (!enviro || !box || molecVec.size() < 1) {
    std::cerr << "Error: buildBoxData(): Could not load molecule data." << std::endl;
    return false;
  }

	// Make nunOfMolecules divisible by the size of the molecule vector.
  int molecMod = enviro->numOfMolecules % molecVec.size();
  if (molecMod != 0)	{
  	enviro->numOfMolecules += molecVec.size() - molecMod;
  }

	// Allocate memory for the molecues in the box.
  box->moleculeCount = enviro->numOfMolecules;
  box->molecules = (Molecule *) malloc(sizeof(Molecule) * box->moleculeCount);

	// If there aren't any molecules, return an error.
  if(box->moleculeCount < 1 || molecVec.size() < 1) {
  	std::cerr << "The simulation environment should have at least 1 " <<
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
  box->atoms 	   = (Atom *) malloc(sizeof(Atom) * box->atomCount);
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

			int a1Idx = molec1.bonds[k].atom1;	// Get the indexes of both atoms
			int a2Idx = molec1.bonds[k].atom2;	// in the bond.

			std::string a1Name = *(molec1.atoms[a1Idx].name); // Get the names of both
			std::string a2Name = *(molec1.atoms[a2Idx].name); // atoms in the bond.

			// Assign the forceConstant and bondDistance values for the bond
			// based on name lookup from the oplsaa.sb file using sb_scanner.
			molec1.bonds[k].forceConstant = sb_scanner.getKBond(a1Name, a2Name);
			molec1.bonds[k].eqBondDist = sb_scanner.getEqBondDist(a1Name, a2Name);

			// Copy the bond from the molecule in the vector to the box.
			box->molecules[j].bonds[k] = molec1.bonds[k];
    }

    // Copy angles from the molecule vector.
    for (int k = 0; k < molec1.numOfAngles; k++) {
  		box->molecules[j].angles[k] = molec1.angles[k];
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
    memcpy(&(box->molecules[offset]), box->molecules, sizeof(Molecule) * molecTypeNum);

		// Copy atoms, bonds, etc. for each molecule.
	  for (int n = 0; n < molecTypeNum; n++) {
    	box->molecules[offset+n].id=offset+n;
      box->molecules[offset+n].atoms = box->molecules[n].atoms + count[0] * m;
      box->molecules[offset+n].bonds =  box->molecules[n].bonds + count[1] * m;
      box->molecules[offset+n].angles =  box->molecules[n].angles + count[2] * m;
      box->molecules[offset+n].dihedrals =  box->molecules[n].dihedrals + count[3] * m;
      box->molecules[offset+n].hops =  box->molecules[n].hops + count[4] * m;
    }

		memcpy(&(box->atoms[m * count[0]]), box->atoms, sizeof(Atom) * count[0]);
	  memcpy(&(box->bonds[m * count[1]]), box->bonds, sizeof(Bond) * count[1]);
	  memcpy(&(box->angles[m * count[2]]), box->angles, sizeof(Angle) * count[2]);
	  memcpy(&(box->dihedrals[m * count[3]]), box->dihedrals, sizeof(Dihedral) * count[3]);
	  memcpy(&(box->hops[m * count[4]]), box->hops, sizeof(Hop) * count[4]);

		// Copy atoms for each molecule.
    for (int k = 0; k < count[0]; k++)	{
    	box->atoms[m * count[0] + k].id = offset * count[0] + k;
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
    }

		// Copy dihedrals for each molecule.
    for(int k = 0; k < count[3]; k++) {
    	box->dihedrals[m * count[3] + k].atom1+= m * count[0];
      box->dihedrals[m * count[3] + k].atom2+= m * count[0];
  	}

		// Copy hops for each molecule.
    for(int k = 0; k < count[4]; k++) {
    	box->hops[m * count[4] + k].atom1 += m * count[0];
      box->hops[m * count[4] + k].atom2 += m * count[0];
    }

  }	// Repeat for every molecule.

  enviro->numOfAtoms = count[0] * molecDiv;

	//generate fcc lattice box
  if (!generatefccBox(box)) {
  	std::cerr << "Error: buildBoxData(): Could not generate FCC box" << std::endl;
    return false;
  }

  return true;
}

bool generatefccBox(Box* box) {
	// Return an error if the box is missing an environment or molecules.
	if (box->environment == NULL || box->molecules == NULL) {
		std::cerr << "Error: generatefccBox(): Box environment or molecules array NULL" << std::endl;
		return false;
	}

	Real cells, dcells, cellL, halfcellL;
	Environment* enviro = box->environment;
	Molecule* molecules = box->molecules;

	//Determine the number of unit cells in each coordinate direction
	dcells = pow(0.25 * (Real) enviro->numOfMolecules, 1.0/3.0);
	cells = (int)(dcells + 0.5);

	//Check if numOfMolecules is a non-fcc number of molecules
	//and increase the number of cells if necessary
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

	//Init all other molecules to initial coordinates
	//Build the lattice from the unit cell by repeatedly translating
	//the four vectors of the unit cell through a distance cellL in
	//the x, y, and z directions
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
// ========================= State Scanner ====================================
// ============================================================================

StateScanner::StateScanner(string filename)
{
	universal_filename = filename;
}

StateScanner::~StateScanner()
{

}

Environment* StateScanner::readInEnvironment()
{
	std::string filename = universal_filename;

    ifstream inFile;
    inFile.open(filename.c_str());
    string line;
    Environment* environment;

    if(inFile.is_open())
    {
        getline(inFile, line);
        environment = getEnvironmentFromLine(line);
    }

    return environment;
}

long StateScanner::readInStepNumber()
{
    std::string filename = universal_filename;

    long startingStep = 0;
    ifstream inFile(filename.c_str());
    string line;

    if (inFile.is_open())
    {
        getline(inFile, line); // Environment
        getline(inFile, line); // Step number

        if (!fromString<long>(line, startingStep))
            return -1;
        return startingStep;
    }

    return -1;
}

vector<Molecule> StateScanner::readInMolecules()
{
	std::string filename = universal_filename;

    vector<Molecule> molecules;
    ifstream inFile(filename.c_str());
    string line;

    if(inFile.is_open())
    {
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
        while(inFile.good())
        {
            getline(inFile, line);
            string hold = line.substr(0, 2);
            switch(section)
            {
                case 0: // id
                    if(hold.compare("= ") == 0)
                    {
                        section++;
                    }
                    else
                    {
                        currentMol.id = atoi(line.c_str());
                    }
                    break;
                case 1: // type
                    if(hold.compare("= ") == 0)
                    {
                        section++;
                    }
                    else
                    {
                        currentMol.type = atoi(line.c_str());
                    }
                    break;
                case 2: // atom
                    if(hold.compare("= ") == 0)
                    {
                        section++;
                    }
                    else
                    {
                       atoms.push_back(getAtomFromLine(line));
                    }
                    break;
                case 3: // bond
                    if(hold.compare("= ") == 0)
                    {
                        section++;
                    }
                    else
                    {
					   bonds.push_back(getBondFromLine(line));
                    }
                    break;
                case 4: // dihedral
                    if(hold.compare("= ") == 0)
                    {
                        section++;
                    }
                    else
                    {
                       dihedrals.push_back(getDihedralFromLine(line));
                    }
                    break;
                case 5: // hop
                    if(hold.compare("= ") == 0)
                    {
                        section++;
                    }
                    else
                    {
                        hops.push_back(getHopFromLine(line));
                    }
                    break;
                case 6: // angle
                    if(hold.compare("==") == 0)
                    {
                        section = 0;
                        molNum++;

                        //convert all vectors to arrays
                        Bond *bondArray = (Bond *) malloc(sizeof(Bond) * bonds.size());
                        Angle *angleArray = (Angle *) malloc(sizeof(Angle) * angles.size());
                        Atom *atomArray = (Atom *) malloc(sizeof(Atom) * atoms.size());
                        Dihedral *dihedralArray = (Dihedral *) malloc(sizeof(Dihedral) * dihedrals.size());
                        Hop *hopArray = (Hop *) malloc(sizeof(Hop) * hops.size());

                        for(int i = 0; i < bonds.size(); i++)
                        {
                            bondArray[i] = bonds[i];
                        }
						for(int i = 0; i < angles.size(); i++)
                        {
                            angleArray[i] = angles[i];
                        }
						for(int i = 0; i < atoms.size(); i++)
                        {
                            atomArray[i] = atoms[i];
                        }
                        for(int i = 0; i < dihedrals.size(); i++)
                        {
                            dihedralArray[i] = dihedrals[i];
                        }
                        for(int i = 0; i < hops.size(); i++)
                        {
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
                    }
                    else
                    {
                       angles.push_back(getAngleFromLine(line));
                    }
                    break;
            }
        }
    }
    else
    {
        return molecules;
    }

   return molecules;
}

Angle StateScanner::getAngleFromLine(string line)
{
    Angle angle = Angle();
    char *tokens;
    char *charLine = (char *) malloc(sizeof(char) * (line.size()+1));
    strcpy(charLine, line.c_str());
    tokens = strtok(charLine, " ");
    int tokenNumber = 0;

    while(tokens != NULL)
    {
       switch(tokenNumber)
       {
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
               if(atoi(tokens) == 1)
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

Atom StateScanner::getAtomFromLine(string line)
{
    Atom atom = createAtom(-1,-1,-1,-1,-1,-1);
    char *tokens;
    char *charLine = (char *) malloc(sizeof(char) * (line.size()+1));
    strcpy(charLine, line.c_str());
    tokens = strtok(charLine, " ");
    int tokenNumber = 0;

    while(tokens != NULL)
    {
        switch(tokenNumber)
        {
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

Bond StateScanner::getBondFromLine(string line)
{
    Bond bond = Bond();
    char *tokens;
    char *charLine = (char *) malloc(sizeof(char) * (line.size()+1));
    strcpy(charLine, line.c_str());
    tokens = strtok(charLine, " ");
    int tokenNumber = 0;

    while(tokens != NULL)
    {
        switch(tokenNumber)
        {
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
                if(atoi(tokens) == 1)
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

Dihedral StateScanner::getDihedralFromLine(string line)
{
    Dihedral dihedral = Dihedral();

    char *tokens;
    char *charLine = (char *) malloc(sizeof(char) * (line.size()+1));
    strcpy(charLine, line.c_str());
    tokens = strtok(charLine, " ");
    int tokenNumber = 0;
    while(tokens != NULL)
    {
        switch(tokenNumber)
        {
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
                if(atoi(tokens) == 1)
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

Environment* StateScanner::getEnvironmentFromLine(string line)
{
    Environment* environment = new Environment();
    char *tokens;
    char *charLine = (char *) malloc(sizeof(char) * (line.size()+1));
    strcpy(charLine, line.c_str());
    tokens = strtok(charLine, " ");
    int numOfAtoms, tokenNumber = 0;
    Real x,y,z,cutoff;

    while(tokens != NULL)
    {
        switch(tokenNumber)
        {
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

void StateScanner::parsePrimaryIndexDefinitions(Environment* enviro, string definitions)
{
    *enviro->primaryAtomIndexConfigLine = definitions;
    char *indexVector;
    char *charLine = (char *) malloc(sizeof(char) * (definitions.size()+1));
    strcpy(charLine, definitions.c_str());
    indexVector = strtok(charLine, ",");

    for (int i = 0; indexVector ; i++)
    {
    	std::string sIndexVector = indexVector;
    	enviro->primaryAtomIndexDefinitions++;
    	char* c;
    	if ((c=strchr(indexVector, '['))!=NULL)
    	{
    	    indexVector = (indexVector + 1);
    	    while ((c=strchr(indexVector, ']'))==NULL)
    	    {
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

Hop StateScanner::getHopFromLine(string line)
{
    Hop hop = Hop();
    char *tokens;
    char *charLine = (char *) malloc(sizeof(char) * (line.size()+1));
    strcpy(charLine, line.c_str());
    tokens = strtok(charLine, " ");
    int tokenNumber = 0;

    while(tokens != NULL)
    {
        switch(tokenNumber)
        {
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

void StateScanner::outputState(Environment *environment, Molecule *molecules, int numOfMolecules, int step, string filename)
{
    ofstream outFile;
    outFile.open(filename.c_str());
    outFile << environment->x << " " << environment->y << " "
        << environment->z << " " << environment->numOfMolecules << " "
        << environment->numOfAtoms << " " << environment->temp << " "
        << environment->cutoff << " " << environment->maxTranslation << " "
        << environment->maxRotation << " " << *environment->primaryAtomIndexConfigLine << " "
        << environment->randomseed << " "
        << std::endl << step << std::endl << std::endl;

    for(int i=0; i<numOfMolecules; i++)
    {
        Molecule currentMol = molecules[i];
        outFile << currentMol.id << std::endl;

        outFile << "= Type" << std::endl;
        outFile << currentMol.type << std::endl;

        // Write atoms
        outFile << "= Atoms" << std::endl;
        for(int j=0; j<currentMol.numOfAtoms; j++)
        {
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
        for(int j=0; j<currentMol.numOfBonds; j++)
        {
            Bond currentBond = currentMol.bonds[j];
            outFile << currentBond.atom1 << " "
                << currentBond.atom2 << " "
                << currentBond.distance << " ";

            if(currentBond.variable)
            {
                outFile << "1" << std::endl;
            }
            else
            {
                outFile << "0" << std:: endl;
            }
        }

        // Write dihedrals
        outFile << "= Dihedrals" << std::endl;
        for(int j=0; j<currentMol.numOfDihedrals; j++)
        {
            Dihedral currentDi = currentMol.dihedrals[j];
            outFile << currentDi.atom1 << " "
                << currentDi.atom2 << " "
                << currentDi.value << " ";

            if(currentDi.variable)
            {
                outFile << "1" << std::endl;
            }
            else
            {
                outFile << "0" << std::endl;
            }
        }

        // Write hops
        outFile << "= Hops" << std::endl;
        for(int j=0; j<currentMol.numOfHops; j++)
        {
            Hop currentHop = currentMol.hops[j];

            outFile << currentHop.atom1 << " "
                << currentHop.atom2 << " "
                << currentHop.hop << std::endl;
        }

        // Write angles
        outFile << "= Angles" << std::endl;
        for(int j=0; j<currentMol.numOfAngles; j++)
        {
            Angle currentAngle = currentMol.angles[j];

            outFile << currentAngle.atom1 << " "
                << currentAngle.atom2 << " "
                << currentAngle.value << " ";

            if(currentAngle.variable)
            {
                outFile << "1" << std::endl;
            }
            else
            {
                outFile << "0" << std::endl;
            }
        }

        outFile << "==" << std::endl;
    }

    outFile.close();
}

//=============================================================================
//======================== SB Scanner - Reads OPLSAA.sb =======================
//=============================================================================
void SBScanner::processBond(string line) {
	string atom1, atom2;
	if(line.substr(1,1) == " ") {
		atom1 = line.substr(0,1);
	} else {
		atom1 = line.substr(0,2);
	}
	if(line.substr(4,1) == " ") {
		atom2 = line.substr(3,1);
	} else {
		atom2 = line.substr(3,2);
	}
	double forceK, bondDist;
	string rest_of_line = line.substr(6, line.size()-6);
	stringstream ss(rest_of_line);
	ss >> forceK >> bondDist;
	bondDataMap[atom1][atom2].kBond = forceK;
	bondDataMap[atom2][atom1].kBond = forceK;
	bondDataMap[atom1][atom2].eqBondDist = bondDist;
	bondDataMap[atom2][atom1].eqBondDist = bondDist;
}

void SBScanner::processAngle(string line) {
	string end1 = line.substr(0, 2);
	string midAtom = line.substr(3, 2);
	string end2 = line.substr(6, 2);
	string rest_of_line = line.substr(8, line.size()-8);
	double angleK, angle;
	stringstream ss(rest_of_line);
	ss >> angleK >> angle;
	angleDataMap[end1][midAtom][end2].kAngle = angleK;
	angleDataMap[end2][midAtom][end1].kAngle = angleK;
	angleDataMap[end1][midAtom][end2].eqAngle = angle;
	angleDataMap[end2][midAtom][end1].eqAngle = angle;
}

void SBScanner::processLine(string line) {
	if(line.at(5) == '-')
		processAngle(line);
	else
		processBond(line);
}

SBScanner::SBScanner() {

}

SBScanner::~SBScanner() {
	bondDataMap.clear();
	angleDataMap.clear();
}

bool SBScanner::readInSB(string filename) {
	fileName = filename;

	if (filename.empty()) {
		std::cerr << "Error: readInOpls(): empty filename given" << std::endl;
		return false;
	}

	int numOfLines=0;
    ifstream oplsScanner(filename.c_str());
    if (!oplsScanner.is_open()) {
    	std::cerr << "Error: readInOpls(): could not open file (" << filename << ")" << std::endl;
        return false;
    } else {
        string line;
        while (oplsScanner.good()) {
            numOfLines++;
            getline(oplsScanner,line);

            //check if it is a commented line,
            //or if it is a title lines
            try {
                if(line.at(0) != '*' && line.at(0) != '\"' && line.at(0) != '#')
				{
					processLine(line);
                    //cout <<(line) << endl;
				}
            } catch (std::out_of_range& e) {
            	// Eat the exception and continue...why aren't we failing?
            }
        }
        oplsScanner.close();
    }
	return true;
}

map<string, BondData> SBScanner::getBondMap(string atom1) {
	if(bondDataMap.count(atom1) == 0) {
		map<string, BondData> blankMap;
		return blankMap;
	} else {
		return bondDataMap[atom1];
	}
}

map<string, map<string, AngleData> > SBScanner::getAngleMap(string endpoint1) {
	if(angleDataMap.count(endpoint1) == 0) {
		map<string, map<string, AngleData> > blankMap;
		return blankMap;
	} else {
		return angleDataMap[endpoint1];
	}
}


Real SBScanner::getKBond(string atom1, string atom2) {
	if(bondDataMap.count(atom1) == 0)
		return -1;
	if(bondDataMap[atom1].count(atom2) == 0)
		return -1;
	return bondDataMap[atom1][atom2].kBond;
}

Real SBScanner::getEqBondDist(string atom1, string atom2) {
	if(bondDataMap.count(atom1) == 0)
		return -1;
	if(bondDataMap[atom1].count(atom2) == 0)
		return -1;
	return bondDataMap[atom1][atom2].eqBondDist;
}

Real SBScanner::getKAngle(string endpoint1, string middleAtom, string endpoint2) {
	if(angleDataMap.count(endpoint1) == 0)
		return -1;
	if(angleDataMap[endpoint1].count(middleAtom) == 0)
		return -1;
	if(angleDataMap[endpoint1][middleAtom].count(endpoint2) == 0)
		return -1;
	return angleDataMap[endpoint1][middleAtom][endpoint2].kAngle;
}

Real SBScanner::getEqAngle(string endpoint1, string middleAtom, string endpoint2) {
	if(angleDataMap.count(endpoint1) == 0)
		return -1;
	if(angleDataMap[endpoint1].count(middleAtom) == 0)
		return -1;
	if(angleDataMap[endpoint1][middleAtom].count(endpoint2) == 0)
		return -1;
	return angleDataMap[endpoint1][middleAtom][endpoint2].eqAngle;
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

static inline std::string& rtrim(std::string& s) {
	s.erase(s.find_last_not_of(TRIMMED_CHARS) + 1);
  return s;
}

static inline std::string& ltrim(std::string& s) {
  std::size_t len = s.find_first_not_of(TRIMMED_CHARS);
  if (len != std::string::npos)
    s.erase(0, len);
  return s;
}

static inline std::string& trim(std::string& s) {
  return ltrim(rtrim(s));
}
