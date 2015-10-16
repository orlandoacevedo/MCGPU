/// @file DeviceQuery.cpp
///
/// This serves as a temporary file to replace DeviceQuery.cu
/// while openACC refactoring occurs.
///
/// @author Mitch Price
/// @date Created 10/10/2015
/// @date Updated 10/10/2015

#include "DeviceQuery.h"

bool openDeviceContext(struct DeviceContext* context, int major, int minor, int deviceIndex)
{
	return false;
}

bool closeDeviceContext(struct DeviceContext* context)
{
	return false;
}

void printDeviceInformation()
{
  std::cout << "CUDA is currently unsupported" << std::endl;
}
