/// @file DeviceQuery.cuh
///
/// Contains device querying declarations and constants.
///
/// @author Joshua Mosby
/// @date Created 2/28/2014
/// @date Updated 3/2/2014

#ifndef DEVICE_QUERY_CUH
#define DEVICE_QUERY_CUH

bool findDevice(int major, int minor, int device, cudaDeviceProp* properties);

bool findBestDevice(int major, int minor, int* device, cudaDeviceProp* properties);

/// Determines if a CUDA device matches the given minimum specifictions.
///
/// @param devProp The set of properties that describe the device capabilities.
/// @param specMajor The minimum allowed version major
/// @param specMinor The minimum allowed version minor
/// @returns True if the device met the minimum requirements; false if the
///     device is invalid or does not meet the minimum requirements.
bool matchSpecs(cudaDeviceProp devProp, int specMajor, int specMinor);

#endif