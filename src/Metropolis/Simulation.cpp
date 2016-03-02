/*
	Driver for the simulation. Takes in a SimulationArgs object and creates the
	the necessary Box type, state file output path, etc.

	Author: Nathan Coleman
	Created: February 21, 2014

	-> February 26, by Albert Wallace
	-> March 28, by Joshua Mosby
	-> April 21, by Nathan Coleman
*/

#include <string>
#include <iostream>
#include <fstream>
#include <time.h>
#include <stdio.h>
#include <omp.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <sstream>

#include "Simulation.h"
#include "SimulationArgs.h"
#include "Box.h"
#include "Metropolis/Utilities/MathLibrary.h"
#include "Metropolis/Utilities/Parsing.h"
#include "SerialSim/SerialBox.h"
#include "SerialSim/SerialCalcs.h"
#include "SerialSim/NeighborList.h"
//#include "ParallelSim/ParallelCalcs.h"
#include "Utilities/FileUtilities.h"
#include "SimBox.h"
#include "SimBoxBuilder.h"
#include "GPUCopy.h"
#include "SimCalcs.h"

#define RESULTS_FILE_DEFAULT "run"
#define RESULTS_FILE_EXT ".results"

Simulation::Simulation(SimulationArgs simArgs) {
	args = simArgs;
	stepStart = 0;

	if (false) {
		/*We need to set this to 1 in parallel mode because it is irrelevant BUT is used in the
		  runtime calculation. See summary equation for explanation.*/
		threadsToSpawn = 1;
	//	box = ParallelCalcs::createBox(args.filePath, args.fileType, &stepStart, &simSteps);
	} else {
		int processorCount = omp_get_num_procs();
		/*We seem to get reasonable performance if we use half as many threads as there are 'logical' processors.
		  We could do performance tuning to get more information on what the ideal number might be.*/
		threadsToSpawn = max(processorCount / 2, 1);

		if (simArgs.threadCount > 0)
			threadsToSpawn = min(omp_get_max_threads(), simArgs.threadCount);

		std::cout << processorCount << " processors detected by OpenMP; using " << threadsToSpawn << " threads." << endl;
		omp_set_num_threads(threadsToSpawn);
		omp_set_dynamic(0); //forces OpenMP to use the exact number of threads specified above (no less)
		box = SerialCalcs::createBox(args.filePath, args.fileType, &stepStart, &simSteps);
	}

	if (box == NULL) {
		std::cerr << "Error: Unable to initialize simulation Box" << std::endl;
		exit(EXIT_FAILURE);
	} else {
		std::cout << "Using seed: " << box->environment->randomseed << std::endl;
		seed(box->environment->randomseed);
	}

	if (args.stepCount > 0)
		simSteps = args.stepCount;
}

Simulation::~Simulation() {
	if (box != NULL) {
		delete box;
		box = NULL;
	}
}

void Simulation::run() {
	std::cout << "Simulation Name: " << args.simulationName << std::endl;

	for (int molIdx = 0; molIdx < box->environment->numOfMolecules; molIdx++) {
		box->keepMoleculeInBox(molIdx);
	}

	if (args.useNeighborList) {
		box->createNeighborList();
	}

	Real oldEnergy_sb = 0;
	Real oldEnergy = 0, currentEnergy = 0;
	Real newEnergyCont = 0, oldEnergyCont = 0;
	Real lj_energy = 0, charge_energy = 0;
	Real new_lj = 0, old_lj = 0;
	Real new_charge = 0, old_charge = 0;

	Real energy_LRC = SerialCalcs::calcEnergy_LRC(box);
	//Real intraMolEnergy = SerialCalcs::calcIntraMolecularEnergy(box, lj_energy, charge_energy);

	Real kT = kBoltz * box->getEnvironment()->temp;
	int accepted = 0;
	int rejected = 0;

	string directory = getcwd(NULL, 0);

	std::string mc ("MCGPU");
	std::size_t found = directory.rfind(mc);

	if (found != std::string::npos)
		directory = directory.substr(0,found+6);

	std::string MCGPU = directory;

	MCGPU.append("/bin/pdbFiles");

	mkdir(MCGPU.data(), 0777);

	clock_t startTime, endTime;
	startTime = clock();

  clock_t function_time_start, function_time_end;
	function_time_start = clock();

	if(args.verboseOutput) {
		log = Logger(VERBOSE);
	} else {
		log = Logger(NON_VERBOSE);
	}

	// Build SimBox below
	SimBoxBuilder builder = SimBoxBuilder(args.useNeighborList, new SBScanner());
	SimBox* sb = builder.build(box);
	GPUCopy::copyIn(sb);
	GPUCopy::setParallel(true);
	SimCalcs::setSB(sb);
	//Calculate original starting energy for the entire system
	if (oldEnergy == 0) {
		if (false) {
			if (args.useNeighborList) {
				log.verbose("Using neighbor-list for energy calculation");
				//oldEnergy = ParallelCalcs::calcSystemEnergy_NLC(box);
			} else {
				log.verbose("Using original parallel energy calculation");
			//	oldEnergy = ParallelCalcs::calcSystemEnergy(box);
			}
		} else {
			if (args.useNeighborList) {
				log.verbose("Using neighbor-list for energy calculation");
				sb->useNLC = true;
				//oldEnergy = SerialCalcs::calcSystemEnergy_NLC(box, lj_energy, charge_energy);
			} else {
				log.verbose("Using original system energy calculation");
				//oldEnergy = SerialCalcs::calcSystemEnergy(box, lj_energy, charge_energy);
			}
		}
		oldEnergy_sb = sb->calcSystemEnergy(lj_energy, charge_energy);
		oldEnergy_sb += energy_LRC;
	}
	function_time_end = clock();
	GPUCopy::copyOut(sb);
	double duration = difftime(function_time_end, function_time_start) / (CLOCKS_PER_SEC * threadsToSpawn);

	//for testing/debug purposes
	std::stringstream durationConv;
	durationConv << "Duration of system energy calculation function: " << duration << " seconds.";
	log.verbose(durationConv.str());
	log.verbose("Threads to spawn: " + threadsToSpawn);

	std::stringstream simStepsConv;
	simStepsConv << "\nRunning " << (simSteps) << " steps\n";
	log.verbose(simStepsConv.str());

	std::string baseStateFile = "";
	if (!args.simulationName.empty()) {
		baseStateFile.append(args.simulationName);
	} else {
		baseStateFile.append("untitled");
	}

	//Loop for each individual step
	for (int move = stepStart; move < (stepStart + simSteps); move++) {
		std::vector<int> neighbors;
		new_lj = 0, old_lj = 0;
		new_charge = 0, old_charge = 0;

		//provide printouts at each pre-determined interval (not at each step)
		if (args.statusInterval > 0 && (move - stepStart) % args.statusInterval == 0) {
			stringstream moveConv;
			moveConv << "Step " << (move) << ":\n--Current Energy: " << oldEnergy_sb << "\n";
			log.verbose(moveConv.str());
		}

		if (args.stateInterval > 0 && move > stepStart && (move - stepStart) % args.stateInterval == 0) {
			log.verbose("");
			saveState(baseStateFile, move, sb);
			log.verbose("");
		}

		//Randomly select index of a molecule for changing
		//int changeIdx = box->chooseMolecule();
		int changeIdx = sb->chooseMolecule();

		//Calculate the current/original/old energy contribution for the current molecule
		if (false) {
			if (args.useNeighborList) {
			//	oldEnergyCont = ParallelCalcs::calcMolecularEnergyContribution_NLC(box, changeIdx, neighbors);
			} else {
			//	oldEnergyCont = ParallelCalcs::calcMolecularEnergyContribution(box, changeIdx);
			}
		} else {
			if (args.useNeighborList) {
				oldEnergyCont = sb->calcMolecularEnergyContribution(old_lj, old_charge, changeIdx, 0);
			} else {
			  oldEnergyCont = sb->calcMolecularEnergyContribution(old_lj, old_charge, changeIdx, 0);
			}
		}

		//Actually translate the molecule at the preselected index
   	sb->changeMolecule(changeIdx);

		//Calculate the new energy after translation
		if (false) {
			if (args.useNeighborList) {
			//	newEnergyCont = ParallelCalcs::calcMolecularEnergyContribution_NLC(box, changeIdx, neighbors);
			} else {
			//	newEnergyCont = ParallelCalcs::calcMolecularEnergyContribution(box, changeIdx);
			}
		} else {
			if (args.useNeighborList) {
				newEnergyCont = sb->calcMolecularEnergyContribution(new_lj, new_charge, changeIdx, 0);
			} else {
				newEnergyCont = sb->calcMolecularEnergyContribution(new_lj, new_charge, changeIdx, 0);
			}
		}

		// Compare new energy and old energy to decide if we should accept or not
		bool accept = false;

		if (newEnergyCont < oldEnergyCont) {
			// Always accept decrease in energy
			accept = true;
		} else {
			// Otherwise use statistics+random number to determine weather to accept increase in energy
			Real x = exp(-(newEnergyCont - oldEnergyCont) / kT);

			if (x >= randomReal(0.0, 1.0)) {
				accept = true;
			} else {
				accept = false;
			}
		}

		if (accept) {
			accepted++;
			oldEnergy_sb += newEnergyCont - oldEnergyCont;
			lj_energy += new_lj - old_lj;
			charge_energy += new_charge - old_charge;
		} else {
			rejected++;
			sb->rollback(changeIdx);
		}

	}
	endTime = clock();
	writePDB(box->getEnvironment(), box->getMolecules(), sb);

	/*This number will understate 'true' time the more threads we have, since not all parts of the program are threaded.
	  However, it is a good enough estimation without adding unnecessary complexity.*/
	double diffTime = difftime(endTime, startTime) / (CLOCKS_PER_SEC * threadsToSpawn);

	currentEnergy = oldEnergy_sb;
	stringstream startConv;
	startConv << "Step " << (stepStart + simSteps) << ":\r\n--Current Energy: " << currentEnergy;
	log.verbose(startConv.str());

	if (args.stateInterval >= 0)
		saveState(baseStateFile, (stepStart + simSteps), sb);

	fprintf(stdout, "\nFinished running %ld steps\n", simSteps);

  fprintf(stdout, "LJ-Energy Subtotal: %.3f\n", lj_energy);
	fprintf(stdout, "Charge Energy Subtotal: %.3f\n", charge_energy);

	fprintf(stdout, "Energy Long-range Correction: %.3f\n", energy_LRC);
	//fprintf(stdout, "Intramolecular Energy: %.3f\n", intraMolEnergy);

	fprintf(stdout, "Final Energy: %.3f\n", currentEnergy);
	fprintf(stdout, "Run Time: %.3f seconds\n", diffTime);
	fprintf(stdout, "Accepted Moves: %d\n", accepted);
	fprintf(stdout, "Rejected Moves: %d\n", rejected);
	fprintf(stdout, "Acceptance Ratio: %.2f%%\n", 100.0 * accepted / (accepted + rejected));

	std::string resultsName;

	if (args.simulationName.empty())
		resultsName = RESULTS_FILE_DEFAULT;
	else
		resultsName = args.simulationName;

	resultsName.append(RESULTS_FILE_EXT);

	// Save the simulation results.
	std::ofstream resultsFile;
	resultsFile.open(resultsName.c_str());

	resultsFile << "######### MCGPU Results File #############" << std::endl;
	resultsFile << "[Information]" << std::endl;
	resultsFile << "Timestamp = " << currentDateTime() << std::endl;

	if (!args.simulationName.empty())
		resultsFile << "Simulation-Name = " << args.simulationName << std::endl;

	if (false) {
		resultsFile << "Simulation-Mode = GPU" << std::endl;
	} else {
		resultsFile << "Simulation-Mode = CPU" << std::endl;
		resultsFile << "Threads-Used = " << threadsToSpawn << std::endl;
	}

	resultsFile << "Starting-Step = " << stepStart << std::endl;
	resultsFile << "Steps = " << simSteps << std::endl;
	resultsFile << "Molecule-Count = " << box->environment->numOfMolecules << std::endl << std::endl;
	resultsFile << "[Results]" << std::endl;

	resultsFile << "LJ-Energy Subtotal: " << lj_energy << std::endl;
	resultsFile << "Charge Energy Subtotal: " << charge_energy << std::endl;

	resultsFile << "Final-Energy = " << currentEnergy << std::endl;
	resultsFile << "Run-Time = " << diffTime << " seconds" << std::endl;
	resultsFile << "Accepted-Moves = " << accepted << std::endl;
	resultsFile << "Rejected-Moves = " << rejected << std::endl;
	resultsFile << "Acceptance-Rate = " << 100.0f * accepted / (float) (accepted + rejected) << "%" << std::endl;

	resultsFile.close();
}

void Simulation::saveState(const std::string& baseFileName, int simStep, const SimBox* sb) {
	StateScanner statescan = StateScanner("");
	std::string stateOutputPath = baseFileName;
	std::string stepCount;

	if (!toString<int>(simStep, stepCount))
		return;

	stateOutputPath.append("_");
	stateOutputPath.append(stepCount); //add the step number to the name of the output file
	stateOutputPath.append(".state");

	log.verbose("Saving state file " + stateOutputPath );

	statescan.outputState(box->getEnvironment(), box->getMolecules(), box->getMoleculeCount(), simStep, stateOutputPath, sb->atomCoordinates);
}

int Simulation::writePDB(Environment sourceEnvironment, Molecule * sourceMoleculeCollection, SimBox* sb) {
  std::string pdbName;

  if (args.simulationName.empty())
    pdbName = RESULTS_FILE_DEFAULT;
  else
    pdbName = args.simulationName;

  pdbName.append(".pdb");
  std::ofstream pdbFile;
  pdbFile.open(pdbName.c_str());

	int numOfMolecules = sourceEnvironment.numOfMolecules;
	pdbFile << "REMARK Created by MCGPU" << std::endl;
	int atomIdx = 0;

	Real** atomCoords = sb->atomCoordinates;

	for (int i = 0; i < numOfMolecules; i++) {
		Molecule currentMol = sourceMoleculeCollection[i];
    for (int j = 0; j < currentMol.numOfAtoms; j++) {
			Atom currentAtom = currentMol.atoms[j];
			pdbFile.setf(std::ios_base::left,std::ios_base::adjustfield);
			pdbFile.width(6);
			pdbFile << "ATOM";
			pdbFile.setf(std::ios_base::right,std::ios_base::adjustfield);
			pdbFile.width(5);
			pdbFile << currentAtom.id + 1;
			pdbFile.width(3); // change from 5
			pdbFile << *currentAtom.name;
			pdbFile.width(6); // change from 4
			pdbFile << "UNK";
			pdbFile.width(6);
			pdbFile << i + 1;
			pdbFile.setf(std::ios_base::fixed, std::ios_base::floatfield);
			pdbFile.precision(3);
			pdbFile.width(12);
			pdbFile << atomCoords[0][atomIdx];
			pdbFile.width(8);
			pdbFile << atomCoords[1][atomIdx];
			pdbFile.width(8);
			pdbFile << atomCoords[2][atomIdx] << std::endl;
			atomIdx++;
    }
  	pdbFile << "TER" << std::endl;
  }
  pdbFile << "END" << std::endl;
  pdbFile.close();

	return 0;
}

const std::string Simulation::currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);

    return buf;
}
