/*
	New Simulation to replace linearSim and parallelSim

	Author: Nathan Coleman
	Last Changed: February 21, 2014
	
	-> February 26, by Albert Wallace
   -> March 28, by Joshua Mosby
*/
#include <stdio.h>

#include "Simulation.h"
#include "SimulationArgs.h"
#include "Box.h"
#include "Metropolis/Utilities/MathLibrary.h"
#include "SerialSim/SerialBox.h"
#include "SerialSim/SerialCalcs.h"
#include "ParallelSim/ParallelCalcs.h"
#include "Utilities/IOUtilities.h"


//Constructor & Destructor
Simulation::Simulation(SimulationArgs simArgs)
{
   IOUtilities ioUtil = IOUtilities(args.configPath);

   if (!ioUtil.readInConfig() )
   {
      fprintf(stderr, "Terminating Simulation...\n\n");
      exit(1);
   }
   args = simArgs;
   box = new Box(ioUtil);
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
   //declare variables common to both parallel and serial
   Molecule *molecules = box->getMolecules();
   Environment *enviro = box->getEnvironment();
   double oldEnergy = 0, currentEnergy = 0;
   double newEnergyCont, oldEnergyCont;
   double  kT = kBoltz * enviro->temp;
   int accepted = 0;
   int rejected = 0;
   
   //calculate old energy
   if (oldEnergy == 0)
   {
      if (args.simulationMode == SimulationMode::Parallel) {
         oldEnergy = ParallelCalcs::calcSystemEnergy();
      }
      else {
         oldEnergy = SerialCalcs::calcSystemEnergy(molecules, enviro);
      }
   }
   
   for(int move = 0; move < simSteps; move++)
   {
      int changeIdx = box->chooseMolecule();
      
      if (args.simulationMode == SimulationMode::Parallel) {
         oldEnergyCont = ParallelCalcs::calcMolecularEnergyContribution(changeIdx);
      }
      else {
         oldEnergyCont = SerialCalcs::calcMolecularEnergyContribution(molecules, enviro, changeIdx);
      }
         
      box->changeMolecule(changeIdx);
      
      if (args.simulationMode == SimulationMode::Parallel) {
         newEnergyCont = ParallelCalcs::calcMolecularEnergyContribution(changeIdx);
      }
      else {
         newEnergyCont = SerialCalcs::calcMolecularEnergyContribution(molecules, enviro, changeIdx);
      }
      
      bool accept = false;
      
      if(newEnergyCont < oldEnergyCont)
      {
         accept = true;
      }
      else
      {
         double x = exp(-(newEnergyCont - oldEnergyCont) / kT);
         
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
   currentEnergy = oldEnergy;
}


