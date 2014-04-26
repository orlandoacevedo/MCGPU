This file contains instructions for how to profile the program
Date created: 3/27/14
Date updated: 3/28/14
Author: Josh Mosby

NOTE: The ./ command is not needed for the app_name

For CUDA functions
nvprof		- collects events/metrics for CUDA kernels
		- actually runs the program, and thus needs full arguments
	nvprof [options] app_name [app_options]

For C/C++ functions
gprof		- profiles serial functions
		- Whenever the program is run, a file called gmon.out is created with the metrics from the last run
		- Gprof does not run the program, it simply analyzes these metrics, thus no command line arguments are needed
	gprof [options] app_name gmon.out


Other tools
CUDA-MEMCHECK 	- detects memory access errors in CUDA applications
	      	- works with all SM architectures
	      	- does not require any special compilation settings
		- should support dynamic parallelism
	cuda-memcheck [options] app_name [app_options]

CUDA-RACECHECK 	- helps identify memory access race conditions in CUDA applications that use shared memory
	       	- only looks at on-chip shared access memory (defined with the __shared__ flag)
		- does not require any special compilation settigns
		- only works on SM architectures 2.0 and above
		- should support dynamic parallelism
	cuda-memcheck --tool racecheck [memcheck-options] app_name [app_options]
