/// @file DeviceQuery.cuh
///
/// Contains device querying declarations and constants.
///
/// @author Joshua Mosby
/// @date Created 2/28/2014
/// @date Updated 2/28/2014

#ifndef DEVICE_QUERY_H
#define DEVICE_QUERY_H

bool matchSpecs(cudaDeviceProp devProp, int specMajor, int specMinor);
bool chooseBestDevice(int specMajor, int specMinor);

#endif