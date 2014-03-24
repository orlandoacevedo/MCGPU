This file contains instructions for how to profile the program
Date: 3/24/14
Author: Josh Mosby

NOTE: at the time of writing, there are still two code bases. This may change slightly
when we have a single code base

For CUDA functions, use nvprof 
nvprof actually runs the program, and thus requires all command line arguments
	nvprof [nvprof options] executable_name [executable arguments]
The ./ command is not needed for the executable_name

For C/C++ functions, use gprof
Whenever the program is run, a file called gmon.out is created with the metrics
from the last run. Gprof does not run the program, it simply analyzes these metrics.
Thus no command line arguments are needed
	gprof [gprof options] executable_name gmon.out
The ./ command is not needed for the executable_name
