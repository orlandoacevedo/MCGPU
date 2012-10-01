/*!\file
  \Class for simulation Box, including Enviroments and points to molocoles.
  \author David(Xiao Zhang).
 
  This file contains implement of SimBox that are used to handle enviroments and common function
  for box.
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <time.h>
#include "SimBox.h"


SimBox::SimBox(char const* ConfigFile)
{
	molecules=NULL;
	enviro=NULL;
}

Molecule * SimBox::GetMolecule()
{
	return molecules;
}

Environment * SimBox::GetEnviro()
{
	return enviro;
}

int SimBox::ReadStateFile(char const* StateFile)
{
	printf("%s\n",StateFile);
	return 0;
}

int SimBox::WriteStateFile(char const* StateFile)
{
	printf("%s\n",StateFile);
	return 0;
}
 
