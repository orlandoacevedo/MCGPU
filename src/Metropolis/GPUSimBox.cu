/*
	New version of GPUSimBox
	Serves as a wrapper for SimBox

	Author: Nathan Coleman
	Last Changed: February 21, 2014
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "GPUSimBox.cuh"

using namespace std;

//Constructor & Destructor
GPUSimBox::GPUSimBox(Config_Scan configScan) : SimBox(configScan){}
GPUSimBox::~GPUSimBox(){}

//Utility
int copyBoxToHost(SimBox *hostBox){return 0;}
int copyBoxToDevice(SimBox *hostBox){return 0;}

//Cuda necessities
void cudaAssert(const cudaError err, const char *file, const int line){}