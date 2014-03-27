//this is simply used to instantiate an instance of an IOUtilities class (which is the class that reads in config info)
//this will help pinpoint errors to be purely during the import process, as it only does one thing, and does not modify the contents. 

#include <stdlib.h>
#include <stdio.h>

#include "IOUtilities.cpp"



int main(int argc, char** argv)
{
	IOUtilities testUtils = IOUtilities("../../../stuff/demoConfigurationAlbert.txt");
	
	return 0;
	
}