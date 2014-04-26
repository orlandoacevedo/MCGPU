//this is simply used to instantiate an instance of an IOUtilities class (which is the class that reads in config info)
//this will help pinpoint errors to be purely during the import process, as it only does one thing, and does not modify the contents. 

#include <stdlib.h>
#include <stdio.h>

#include "IOUtilities.cpp"
#include "../Box.cpp"
#include "MathLibrary.h"
#include "../SerialSim/SerialBox.cpp"
#include "../SerialSim/SerialCalcs.cpp"

const double kBoltz = 0.00198717;

void run()
{
	IOUtilities testUtils = IOUtilities("../../../stuff/demoConfigurationAlbert.txt");
	
	Box * box;
	
	box = new Box(testUtils);
	
	std::cout << "IOUtil CZM0001"<< std::endl;

	int simSteps = testUtils.numOfSteps;
	
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
         oldEnergy = SerialCalcs::calcSystemEnergy(molecules, enviro);
   }
   
   for(int move = 0; move < simSteps; move++)
   {
      int changeIdx = box->chooseMolecule();
      
         oldEnergyCont = SerialCalcs::calcMolecularEnergyContribution(molecules, enviro, changeIdx);
         
      box->changeMolecule(changeIdx);
      
         newEnergyCont = SerialCalcs::calcMolecularEnergyContribution(molecules, enviro, changeIdx);
      
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


int main(int argc, char** argv)
{
	run();
	return 0;
	
		
	
}