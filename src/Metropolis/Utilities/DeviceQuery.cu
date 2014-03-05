/// @file DeviceQuery.cu
///
/// Contains functions to query the current system device capabilities
/// and set the current CUDA device that meets the minimum user
/// specifications.
///
/// @author Joshua Mosby
/// @date Created 2/28/2014
/// @date Updated 3/2/2014

//This program checks if there is a CUDA capable graphics card 
//and selects the best one
#include <stdio.h>
#include <stdlib.h>

#include "DeviceQuery.h"
#include "DeviceQuery.cuh"

bool openDeviceContext(struct DeviceContext* context, int major, int minor, int deviceIndex)
{
	cudaDeviceProp properties;
	int device;

	if (!context)
		return false;

	if (deviceIndex >= 0)
	{
		device = deviceIndex;
		if (!findDevice(major, minor, device, &properties))
			return false;
	}
	else
	{
		if (!findBestDevice(major, minor, &device, &properties))
			return false;
	}

	cudaSetDevice(device);

	fprintf(stdout, "Currently using CUDA capable device %d\n", device);

	context->index = device;
	context->major = properties.major;
	context->minor = properties.minor;

	return true;
}

bool closeDeviceContext(struct DeviceContext* context)
{
	if (!context)
		return false;

	if (!context->isOpen())
	{
		fprintf(stderr, "Error: Could not clean up opened device context.\n");
		return false;
	}
	
	cudaDeviceSynchronize();
	cudaDeviceReset();
	*context = DeviceContext();

	return true;
}

void printDeviceInformation()
{
	int devCount = 0;
	cudaGetDeviceCount(&devCount);

	if (devCount <= 0)
	{
		fprintf(stdout, "There are no available devices detected.\n");
		return;
	}

	fprintf(stdout, "There are %d device(s) detected.\n\n", devCount);

	for (int i = 0; i < devCount; ++i)
	{
		cudaDeviceProp info;
		cudaGetDeviceProperties(&info, i);

		fprintf(stdout, "Device %d: %s\n", i, info.name);
		fprintf(stdout, "Clock Rate: %d\n", info.clockRate);
		fprintf(stdout, "CUDA Version: %d.%d\n", info.major, info.minor);
		fprintf(stdout, "Processor Count: %d\n", info.multiProcessorCount);
		fprintf(stdout, "Global Memory: %d Bytes\n", info.totalGlobalMem);
		fprintf(stdout, "Constant Memory: %d Bytes\n", info.totalConstMem);
		fprintf(stdout, "Warp Size: %d\n", info.warpSize);
		fprintf(stdout, "\n");
	}
}

bool findDevice(int major, int minor, int device, cudaDeviceProp* properties)
{
	int devCount = 0;
	cudaGetDeviceCount(&devCount);

	if (devCount <= 0)
	{
		fprintf(stderr, "Error: There are no available devices detected.\n");
		return false;
	}

	fprintf(stdout, "There are %d device(s) detected.\n", devCount);
	fprintf(stdout, "Attempting to initialize device %d...\n", device);

	if (device < 0 || device >= devCount)
	{
		fprintf(stderr, "Error: The specified device index is not a valid device.\n");
		return false;
	}

	cudaGetDeviceProperties(properties, device);

	if (properties->computeMode == cudaComputeModeProhibited)
	{
		fprintf(stderr, "Error: Device %d is running in <Compute Mode Prohibited> and cannot be used.\n", device);
        return false;
	}

	if (properties->major < 1)
	{
		fprintf(stderr, "Error: Device %d does not support CUDA operations.\n", device);
		return false;
	}

	if (!matchSpecs(*properties, major, minor))
	{
		fprintf(stderr, "Error: Device %d does not meet the minimum specs.\n", device);
		return false;
	}

	return true;
}

bool findBestDevice(int major, int minor, int* device, cudaDeviceProp* properties)
{
	int devCount = 0;
	cudaGetDeviceCount(&devCount);

	if (devCount <= 0)
	{
		fprintf(stderr, "Error: There are no available devices detected.\n");
		return false;
	}

	fprintf(stdout, "There are %d device(s) detected.\n", devCount);
	fprintf(stdout, "Attempting to find the best device...\n", device);

	cudaDeviceProp devPropArray[devCount];
	int devIndexArray[devCount];
	int matchedDevCount = 0; 
	
	for (int i = 0; i < devCount; i++)
	{
		cudaDeviceProp* currentProp = &(devPropArray[matchedDevCount]);
		int* currentIndex = &(devIndexArray[matchedDevCount]);

		cudaGetDeviceProperties(currentProp, i);
		*currentIndex = i;

		if (currentProp->computeMode != cudaComputeModeProhibited && currentProp->major >= 1)
		{
			if (matchSpecs(*currentProp, major, minor))
			{
				matchedDevCount++;
			}
		}
	}

	if (matchedDevCount == 0)
	{
		printf("Error: No CUDA capable devices are available that meet the minimum specs detected.\n");
		return false;
	}

	int devID = 0;
	for (int i = 1; i < matchedDevCount; i++)
	{
		if (matchSpecs(devPropArray[i], devPropArray[devID].major, devPropArray[devID].minor))
		{
			devID = i;
		}
	}

	*properties = devPropArray[devID];
	*device = devIndexArray[devID];

	return true;
}

//This function checks the device (devProp) against the specifications
//It returns true if the device meets specifications, false otherwise
bool matchSpecs(cudaDeviceProp devProp, int specMajor, int specMinor)
{
	//if device major is greater, return true
	if (devProp.major > specMajor)
	{
		return true;
	}
	//if device major is equal, look at minor
	else if (devProp.major == specMajor)
	{
		//if minor is less, return false
		if (devProp.minor < specMinor)
		{
			return false;
		}
		//if minor is greater or equal, return true
		else
		{
			return true;
		}
	}
	//if device major is less, return false
	else
	{
		return false;
	}

	/* return (devProp.major > specMajor || (devProp.major == specMajor && devProp.minor < specMinor)) */
}