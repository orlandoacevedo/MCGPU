/// @file DeviceQuery.cuh
///
/// Contains device querying declarations and constants.
///
/// @author Joshua Mosby
/// @date Created 2/28/2014
/// @date Updated 3/2/2014

#ifndef DEVICE_QUERY_H
#define DEVICE_QUERY_H

#define DEVICE_ANY -1
#define MIN_MAJOR_VER 2
#define MIN_MINOR_VER 1

struct DeviceContext
{
	int index;
	int major;
	int minor;

	DeviceContext()
	{ 
		index = -1;
		major = 0;
		minor = 0;
	}

	bool isOpen() 
	{
		return index >= 0;
	}
};

bool openDeviceContext(DeviceContext* context, int major, int minor, int deviceIndex);

bool closeDeviceContext(DeviceContext* context);

void printDeviceInformation();

#endif