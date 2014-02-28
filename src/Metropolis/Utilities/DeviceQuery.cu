//This program checks if there is a CUDA capable graphics card 
//and selects the best one
#include <stdio.h>
#include <stdlib.h>

//This function checks the device (devProp) against the specifications
//It returns true if the device meets specifications, false otherwise
bool matchSpecs(cudaDeviceProp devProp, int specMajor, int specMinor) {
	//if device major is greater, return true
	if (devProp.major > specMajor) {
		return true;
	}
	//if device major is equal, look at minor
	else if (devProp.major == specMajor) {
		//if minor is less, return false
		if (devProp.minor < specMinor) {
			return false;
		}
		//if minor is greater or equal, return true
		else{
			return true;
		}
	}
	//if device major is less, return false
	else {
		return false;
	}
}

//This function checks for the Device and chooses the best one
bool chooseBestDevice(int specMajor, int specMinor) {
	//declare variables
	int devCount;
	cudaDeviceProp devProp;
	bool match;
	//get the number of CUDA devices
   cudaGetDeviceCount(&devCount);
	printf("There are %d CUDA device(s)\n", devCount);

	//take appropriate action based on number of devices
	if (devCount == 0) {
		printf("No CUDA capable cards found\n");
		return false;		
	}
	else if (devCount == 1) {
		printf("One CUDA capable card found\n");
		cudaGetDeviceProperties(&devProp, 0);
		//make sure card matches minimum specifications
		printf("%s has capability %d.%d\n", devProp.name, devProp.major, devProp.minor);
		printf("Minimum Capability: %d.%d\n", specMajor, specMinor);
		match = matchSpecs(devProp, specMajor, specMinor);		
		if (match) {
			printf("%s matches specifications\n", devProp.name);
			return true;
		}
		else {
			printf("%s does not match specifications\n", devProp.name);
			return false;
		}
	}	
	else {
		cudaDeviceProp devPropArray[devCount];//Optimization issue due to unused memory
      int devIndexArray[devCount]; // List of device indexes
		int newDevCount = 0; 
      //Loop through available devices 
		for (int i = 0; i < devCount; i++) {
			cudaGetDeviceProperties(&devProp, i);
			if (matchSpecs(devProp, specMajor, specMinor)) {
            devPropArray[newDevCount] = devProp;
            devIndexArray[newDevCount] = i;
            newDevCount++;
         }
		}
      if (newDevCount == 0) {
         printf("GPU cards detected but none meet minumum specifications.\n");
         return false;
      }
      
		devID = 0;
      for (i = 1; i < newDevCount; i++) {
         if (matchSpecs(devPropArray[i], devPropArray[devID].major, devPropArray[devID].minor)) {
            devID = i;
         }
      }
      
      int bestDeviceIndex = devIndexArray[devID];
      cudaSetDevice(bestDeviceIndex);
		return true;	
	}
}
/*
//Arguments are Major and Minor
int main(int argc, char** argv) {
	if (argc != 3)
	{
		printf("need major and minor version\n");
		exit(0);
	}
	int major = strtol(argv[1], NULL, 10);
	int minor = strtol(argv[2], NULL, 10);
	chooseBestDevice(major, minor);
}*/
