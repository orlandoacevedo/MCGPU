/*
	New Simulation to replace linearSim and parallelSim

	Author: Nathan Coleman
	Last Changed: February 21, 2014
	
	-> February 26, by Albert Wallace
	-> March 28, by Joshua Mosby
*/

#include <string>
#include <iostream>
#include <time.h>

#include "Simulation.h"
#include "SimulationArgs.h"
#include "Box.h"
#include "Metropolis/Utilities/MathLibrary.h"
#include "Metropolis/Utilities/Parsing.h"
#include "SerialSim/SerialBox.h"
#include "SerialSim/SerialCalcs.h"
#include "ParallelSim/ParallelCalcs.h"
#include "Utilities/FileUtilities.h"


//Constructor & Destructor
Simulation::Simulation(SimulationArgs simArgs)
{
	args = simArgs;

	if (simArgs.simulationMode == SimulationMode::Parallel)
		box = ParallelCalcs::createBox(args.filePath, args.fileType, &simSteps);
	else
		box = SerialCalcs::createBox(args.filePath, args.fileType, &simSteps);

	if (box == NULL)
	{
		std::cerr << "Error: Unable to initialize simulation Box" << std::endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		std::cout << "Using seed: " << box->environment->randomseed << std::endl;
		seed(box->environment->randomseed);
	}

	if (args.stepCount > 0)
		simSteps = args.stepCount;
}

Simulation::~Simulation()
{
	if(box != NULL)
	{
		delete box;
		box = NULL;
	}
}

//Utility
void Simulation::run()
{
	std::cout << "Simulation Name: " << args.simulationName << std::endl;
	//declare variables common to both parallel and serial
	Molecule *molecules = box->getMolecules();
	Environment *enviro = box->getEnvironment();
	
	Real oldEnergy = 0, currentEnergy = 0;
	Real newEnergyCont, oldEnergyCont;
	Real  kT = kBoltz * enviro->temp;
	int accepted = 0;
	int rejected = 0;

	clock_t startTime, endTime;
    startTime = clock();
	
	//calculate old energy
	if (oldEnergy == 0)
	{
		if (args.simulationMode == SimulationMode::Parallel)
		{
			oldEnergy = ParallelCalcs::calcSystemEnergy(box);
		}
		else
		{
			oldEnergy = SerialCalcs::calcSystemEnergy(molecules, enviro);
		}
	}
	
	std::cout << std::endl << "Running " << simSteps << " steps" << std::endl << std::endl;
	
	//determine where we want the state file to go
	std::string baseStateFile = "";	
	if (!args.simulationName.empty())
	{
		baseStateFile.append(args.simulationName);
	}
	else
	{
		baseStateFile.append("untitled");
	}
	
	for(int move = 0; move < simSteps; move++)
	{
		if (args.statusInterval > 0 && move % args.statusInterval == 0)
		{
			std::cout << "Step " << move << ":\n--Current Energy: " << oldEnergy << std::endl;	
		}
		
		if (args.stateInterval > 0 && move > 0 && move % args.stateInterval == 0)
		{
			saveState(baseStateFile, move);
		}
		
		int changeIdx = box->chooseMolecule();
		
		if (args.simulationMode == SimulationMode::Parallel)
		{
			oldEnergyCont = ParallelCalcs::calcMolecularEnergyContribution(box, changeIdx);
		}
		else
		{
			oldEnergyCont = SerialCalcs::calcMolecularEnergyContribution(molecules, enviro, changeIdx);
		}
			
		box->changeMolecule(changeIdx);
		
		if (args.simulationMode == SimulationMode::Parallel)
		{
			newEnergyCont = ParallelCalcs::calcMolecularEnergyContribution(box, changeIdx);
		}
		else
		{
			newEnergyCont = SerialCalcs::calcMolecularEnergyContribution(molecules, enviro, changeIdx);
		}
		
		bool accept = false;
		
		if(newEnergyCont < oldEnergyCont)
		{
			accept = true;
		}
		else
		{
			Real x = exp(-(newEnergyCont - oldEnergyCont) / kT);
			
			if(x >= randomReal(0.0, 1.0))
			{
				accept = true;
			}
			else
			{
				accept = false;
			}
		}
		
		if(accept)
		{
			accepted++;
			oldEnergy += newEnergyCont - oldEnergyCont;
		}
		else
		{
			rejected++;
			//restore previous configuration
			box->rollback(changeIdx);
		}
	}

	endTime = clock();
    double diffTime = difftime(endTime, startTime) / CLOCKS_PER_SEC;

	std::cout << "Step " << simSteps << ":\r\n--Current Energy: " << oldEnergy << std::endl;
	currentEnergy = oldEnergy;

	// Save the final state of the simulation
	if (args.stateInterval >= 0)
	{
		saveState(baseStateFile, simSteps);
	}
	
	std::cout << std::endl << "Finished running " << simSteps << " steps" << std::endl;
	std::cout << "Final Energy: " << currentEnergy << std::endl;
	std::cout << "Run Time: " << diffTime << " seconds" << std::endl;
	std::cout << "Accepted Moves: " << accepted << std::endl;
	std::cout << "Rejected Moves: " << rejected << std::endl;
	std::cout << "Acceptance Ratio: " << 100.0 * accepted / (accepted + rejected) << '\%' << std::endl;
}

void Simulation::saveState(const std::string& baseFileName, int simStep)
{
	StateScanner statescan = StateScanner("");
	std::string stateOutputPath = baseFileName;
	std::string stepCount;

	if (!toString<int>(simStep, stepCount))
		return;

	stateOutputPath.append("_");
	stateOutputPath.append(stepCount); //add the step number to the name of the output file
	stateOutputPath.append(".state");

	std::cout << "Saving state file " << stateOutputPath << std::endl;

	statescan.outputState(box->getEnvironment(), box->getMolecules(), box->getMoleculeCount(), stateOutputPath);
}