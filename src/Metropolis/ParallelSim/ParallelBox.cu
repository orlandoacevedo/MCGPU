/*
	New version of GPUSimBox
	Serves as a wrapper for SimBox

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

#include "ParallelBox.cuh"

using namespace std;

//Constructor & Destructor
ParallelBox::ParallelBox(Config_Scan configScan){}
ParallelBox::~ParallelBox(){}

//Utility
int copyBoxToHost(){return 0;}
int copyBoxToDevice(){return 0;}