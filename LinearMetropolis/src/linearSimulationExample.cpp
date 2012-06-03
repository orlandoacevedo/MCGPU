#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <sstream>
#include "../../Utilities/src/Opls_Scan.h"
#include "../../Utilities/src/Config_Scan.h"
#include "../../Utilities/src/metroUtil.h"
#include "../../Utilities/src/Zmatrix_Scan.h"
#include "../../Utilities/src/State_Scan.h"
#include "../src/metroLinearUtil.h"
   
using namespace std;

// boltzman constant
const double kBoltz = 1.987206504191549E-003;
const double maxRotation = 15.0; // degrees

stringstream ss;

double randomFloat(const double start, const double end)
{
    return (end-start) * (double(rand()) / RAND_MAX) + start;
}

void runLinear(Molecule *molecules, Environment *enviro, int numberOfSteps, string stateFile, string pdbFile){
    int accepted = 0; // number of accepted moves
    int rejected = 0; // number of rejected moves
    double maxTranslation = enviro->maxTranslation;
    double temperature = enviro->temperature;
    double kT = kBoltz * temperature;

    ss << "Assigning Molecule Positions..." << endl;
    writeToLog(ss);
	generatefccBox(molecules, enviro);
   
    ss << "Finished Assigning Molecule Positions" << endl;
    writeToLog(ss);
		
	double currentEnergy = 0.0;
	double oldEnergy = 0.0;
	double newEnergy = 0.0;
	
	// Compute long-range correction to LJ energy,
    // only need to do once for NVT simulation
    double energyLRC = Energy_LRC(molecules, enviro);
    ss << "Long-range cutoff energy: "<< energyLRC << endl;
    writeToLog(ss);

    printState(enviro, molecules, enviro->numOfMolecules, "initialState");

    ss << "Starting first step in Simulation" <<endl;
    writeToLog(ss);
	 
    for(int move = 0; move < numberOfSteps; move++){
        if(move < 1){
        	oldEnergy = calcEnergyWrapper_NLC(molecules, enviro);
        	// Include long-range correction to LJ energy below
        	oldEnergy += energyLRC;
        }
        else{
        	oldEnergy = currentEnergy;
        }

        //Pick a molecule to move
        int moleculeIndex = randomFloat(0, enviro->numOfMolecules);
        Molecule toMove = molecules[moleculeIndex];
        Molecule oldToMove;
        copyMolecule(&oldToMove, &toMove);

        //Pick an atom in the molecule about which to rotate
        int atomIndex = randomFloat(0, molecules[moleculeIndex].numOfAtoms);
        Atom vertex = molecules[moleculeIndex].atoms[atomIndex];

        const double deltaX = randomFloat(-maxTranslation, maxTranslation);
        const double deltaY = randomFloat(-maxTranslation, maxTranslation);
        const double deltaZ = randomFloat(-maxTranslation, maxTranslation);

        const double degreesX = randomFloat(-maxRotation, maxRotation);
        const double degreesY = randomFloat(-maxRotation, maxRotation);
        const double degreesZ = randomFloat(-maxRotation, maxRotation); 

        toMove = moveMolecule(toMove, vertex, deltaX, deltaY, deltaZ,
        degreesX, degreesY, degreesZ);

        keepMoleculeInBox(&toMove, enviro);

        molecules[moleculeIndex] = toMove;

        newEnergy = calcEnergyWrapper_NLC(molecules, enviro);
        // Include long-range correction to LJ energy below
        newEnergy += energyLRC;

        bool accept = false;

        if(newEnergy < oldEnergy){
            accept = true;
        }
        else{
            double x = exp(-(newEnergy - oldEnergy) / kT);

            if(x >= randomFloat(0.0, 1.0)){
                accept = true;
            }
            else{
                accept = false;
            }
        }

        if(accept){
            accepted++;
            currentEnergy = newEnergy;
        }
        else{
            rejected++;
            currentEnergy = oldEnergy;
            //restore previous configuration
            molecules[moleculeIndex] = oldToMove;
        }


        //Print the state every 100 moves.
        if(move % 100 == 0){
            char fileName[50];
            sprintf(fileName, "%dState.state", move);
            string fileNameStr(fileName);
            printState(enviro, molecules, enviro->numOfMolecules, fileNameStr);

			ss << "Step Number: "<< move <<  endl;
			ss << "Current Energy: " << currentEnergy << endl;
			cout << ss.str();

            ss << "Accepted: "<< accepted << endl <<"Rejected: "<< rejected << endl;
            writeToLog(ss);
        }
    }
    char fileName[50];
    sprintf(fileName, "FinalState.state");
    string fileNameStr(fileName);
    printState(enviro, molecules, enviro->numOfMolecules, fileNameStr);
    writePDB(molecules, *enviro, pdbFile);
    
    ss << "Steps Complete"<<endl;        
    ss << "Final Energy: " << currentEnergy << endl;
    ss << "Accepted Moves: " << accepted << endl;
    ss << "Rejected Moves: " << rejected << endl;
    ss << "Acceptance Rate: " << (int) ((float) accepted/ (float) numberOfSteps*100) << "%" << endl;
	 cout << ss.str();
	 writeToLog(ss);
}

/**
  ./bin/linearSim flag path/to/config/file
    flags:
        -z run from z matrix file spcecified in the configuration file
        -s run from state input file specified in the configuration file
*/
int main(int argc, char ** argv){
    writeToLog("",START);
    clock_t startTime, endTime;
    startTime = clock();

    if(argc != 3){
        printf("Error.  Expected usage: run_gpu bin/linearSim {-z | -s} {configPath}\n");
        exit(1);
    }

    // z matrix or state flag
    string flag = argv[1];

    //path to the configuration file
    string configPath = argv[2];

    //demo configuration path: "bin/demoConfiguration.txt";
    //Configuration file scanner
    Config_Scan configScan(configPath);
    configScan.readInConfig();

    //Environment for the simulation
    Environment enviro;
    unsigned long simulationSteps = configScan.getSteps();
    Molecule *molecules;

    //Simulation will run based on the zMatrix and configuration Files
    if(flag.compare("-z") == 0){
        ss << "Running simulation based on Z-Matrix File"<<endl;
        cout<<ss.str()<<endl; writeToLog(ss);

        //get environment from the config file
        enviro = configScan.getEnviro();
		  ss << "Reading Configuation File \nPath: " << configScan.getConfigPath() << endl;
        cout<<ss.str()<<endl; writeToLog(ss);

        //set up Opls scan 
        ss << "Reading OPLS File \nPath: " << configScan.getOplsusaparPath() << endl;
        cout<<ss.str()<<endl; writeToLog(ss);
			
        string oplsPath = configScan.getOplsusaparPath();
        Opls_Scan oplsScan (oplsPath);
        oplsScan.scanInOpls(oplsPath);
		  ss << "OplsScan and OPLS ref table Created " << endl;
         cout<<ss.str()<<endl; writeToLog(ss);

        //set up zMatrixScan
        ss << "Reading Z-Matrix File \nPath: " << configScan.getZmatrixPath() << endl;
        cout<<ss.str()<<endl; writeToLog(ss);
        Zmatrix_Scan zMatrixScan (configScan.getZmatrixPath(), &oplsScan);
        if (zMatrixScan.scanInZmatrix() == -1){
            ss << "Error, Could not open: " << configScan.getZmatrixPath() << endl;
            cerr << ss.str()<< endl;
            writeToLog(ss);
            exit(1);
        }
        ss << "Opened Z-Matrix File \nBuilding "<< enviro.numOfMolecules << " Molecules..." << endl;
        cout<<ss.str()<<endl; writeToLog(ss);

        //Convert molecule vectors into an array
        int moleculeIndex = 0;
        int atomCount = 0;

        vector<Molecule> molecVec = zMatrixScan.buildMolecule(atomCount);
        int molecMod = enviro.numOfMolecules % molecVec.size();
        if (molecMod != 0){
            enviro.numOfMolecules += molecVec.size() - molecMod;
            cout << "Number of molecules not divisible by specified z-matrix. Changing number of molecules to: " << enviro.numOfMolecules << endl;
        }
        molecules = (Molecule *)malloc(sizeof(Molecule) * enviro.numOfMolecules);
            
        while(moleculeIndex < enviro.numOfMolecules){
            molecVec = zMatrixScan.buildMolecule(atomCount);
            //cycle through the number of molecules from the zMatrix
            for(int j = 0; j < molecVec.size(); j++){
                //Copy data from vector to molecule
                Molecule molec1 = molecVec[j];

                molecules[moleculeIndex].atoms = (Atom *)malloc(sizeof(Atom) * molec1.numOfAtoms);
                molecules[moleculeIndex].bonds = (Bond *)malloc(sizeof(Bond) * molec1.numOfBonds);
                molecules[moleculeIndex].angles = (Angle *)malloc(sizeof(Angle) * molec1.numOfAngles);
                molecules[moleculeIndex].dihedrals = (Dihedral *)malloc(sizeof(Dihedral) * molec1.numOfDihedrals);
                molecules[moleculeIndex].hops = (Hop *)malloc(sizeof(Hop) * molec1.numOfHops);

                molecules[moleculeIndex].id = molec1.id;
                molecules[moleculeIndex].numOfAtoms = molec1.numOfAtoms;
                molecules[moleculeIndex].numOfBonds = molec1.numOfBonds;
                molecules[moleculeIndex].numOfDihedrals = molec1.numOfDihedrals;
                molecules[moleculeIndex].numOfAngles = molec1.numOfAngles;
                molecules[moleculeIndex].numOfHops = molec1.numOfHops;

                //get the atoms from the vector molecule
                for(int k = 0; k < molec1.numOfAtoms; k++){
                    molecules[moleculeIndex].atoms[k] = molec1.atoms[k];
                }               
               
                //assign bonds
                for(int k = 0; k < molec1.numOfBonds; k++){
                    molecules[moleculeIndex].bonds[k] = molec1.bonds[k];
                }

                //assign angles
                for(int k = 0; k < molec1.numOfAngles; k++){
                    molecules[moleculeIndex].angles[k] = molec1.angles[k];
                }

                //assign dihedrals
                for(int k = 0; k < molec1.numOfDihedrals; k++){
                    molecules[moleculeIndex].dihedrals[k] = molec1.dihedrals[k];
                }

                atomCount += molecules[moleculeIndex].numOfAtoms;
               
                moleculeIndex++;
            }
        }
        enviro.numOfAtoms = atomCount;
		  ss << "Molecules Created into an Array" << endl;
        writeToLog(ss);
    }       
    //Simulation will run based on the state file
    else if(flag.compare("-s") == 0){
        ss << "Running simulation based on State File"<<endl;
        cout<<ss.str()<<endl; writeToLog(ss);  
      	
        //get path for the state file
        string statePath = configScan.getStatePath();
        ss << "Reading State File \nPath: " << statePath << endl;
        cout<<ss.str()<<endl; writeToLog(ss);  
      	
        //get environment from the state file
        enviro = readInEnvironment(statePath);
        ss << "Scanning in Enviorment " << endl;
        cout<<ss.str()<<endl; writeToLog(ss);
        
        //get vector of molecules from the state file
        vector<Molecule> molecVec = readInMolecules(statePath);
        enviro.numOfMolecules = molecVec.size();
        
        //convert vector of molecules to array
        molecules = (Molecule *)malloc(sizeof(Molecule) * molecVec.size());
        for(int i = 0; i < molecVec.size(); i++){
            molecules[i] = molecVec[i];
        }
    }
    else{
        ss << "Error, Unknown flag"<<endl;
        cout << ss.str();
		  writeToLog(ss);
        exit(1);
    }

    ss << "\nBeginning simulation with: " << endl;
    ss << "\tmolecules "<< enviro.numOfMolecules << endl;
    ss << "\tatoms: "<< enviro.numOfAtoms << endl;
    ss << "\tsteps: "<< simulationSteps << endl;
    cout << ss.str() <<endl;
    writeToLog(ss);
    runLinear(molecules, &enviro, simulationSteps, configScan.getStateOutputPath(),
    configScan.getPdbOutputPath());
         
    for (int i = 0; i < enviro.numOfMolecules; i++){
        free(molecules[i].atoms);
        free(molecules[i].bonds);
        free(molecules[i].angles);
        free(molecules[i].dihedrals);
        free(molecules[i].hops);
    } 

    free(molecules);
    endTime = clock();
    double diffTime = difftime(endTime, startTime) / CLOCKS_PER_SEC;
    ss << "\nSimulation Complete \nRun Time: " << diffTime << endl;
    cout << ss.str() <<endl;
    writeToLog(ss);      
    writeToLog("",END);
}
