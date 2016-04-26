/// @file CudaUtilities.cuh
///
/// Collection of functions and macros that augment the base CUDA
/// libraries.
///
/// @author Tavis Maclellan

#ifndef CUDAUTILITIES_CUH
#define CUDAUTILITIES_CUH

/// Handles error codes returned from CUDA function calls.
#define cudaCheck(call) { cudaAssert((call), #call, __FILE__, __LINE__); }

/// Releases memory allocated on the GPU device.
#define cudaFREE(ptr) if(ptr!=NULL) { cudaFree(ptr);ptr=NULL;}

/// Asserts that a CUDA error code is not fatal. If the error code is not
/// a success, the function will output the error message along with the 
/// file name and line number of the error. By default the function will
/// terminate the program if there is an error; override this behavior
/// by setting the abort parameter to false.
inline void cudaAssert(const cudaError_t code, char const *const func, const char* file, const int line, bool abort=true)
{
	if (code != cudaSuccess) 
	{
		fprintf(stderr,"CUDA Error: %s(%d): %s: %s\n", file, line, func, cudaGetErrorString(code));
		if (abort)
		{
			cudaDeviceReset();
			exit(code);
		}
	}
}

#endif