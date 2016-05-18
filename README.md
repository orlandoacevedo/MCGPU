MCGPU (Monte Carlo on Graphics Processing Units)
===============================================================

##Requirements
 * [PGI Accelerator C/C++ Compiler with OpenACC](https://www.pgroup.com/resources/accel.htm) *or*
   [OpenACC Toolkit](https://developer.nvidia.com/openacc-toolkit) (free for academic use)
    * *Note*: If you are using the Alabama Supercomputer Center's Dense Memory Cluster (DMC), type ```module load pgi``` to load the PGI compilers.
 * For GPU Execution: NVIDIA CUDA-capable graphics card
    * Tested on NVIDIA Kepler K20m and K40
    * By default, the PGI compilers target Fermi-generation GPUs (Compute Capability 2.0) and higher
 * Linux operating system

##Build
```
git clone git://github.com/orlandoacevedo/MCGPU.git
cd MCGPU/
make
```

*Note*: To build in debug mode (so the executable can be debugged using cuda-gdb), use BUILD=debug:
```
make BUILD=debug
```

##Run
###To Run a Simulation on a local machine:
```
cd /path/to/MCGPU/
cd bin/
./metrosim ./[configuration file] [options]

Examples:
./metrosim ./indole4000.config --threads 6 -s --name indole4000 -n 1000 -i 100 -k
runs job using 6 OpenMP CPU threads, for 1000 steps, printing out every 100 intervals

./metrosim ./indole4000.config -p --name indole4000 -n 1000 -i 100 -k
runs job on best available GPU, for 1000 steps, printing out every 100 intervals
```
Where `[configuration file]` is a .config file containing configuration information, and `[options]` are command-line options. An example demo.config can be found in the resources folder. See below for specific .config file documentation and all command-line options available.

###To Run a Simulation on the Alabama Supercomputer Center's DMC:
```
cd /path/to/MCGPU/
run_gpu demo_script.txt
```
Choose a batch job queue:
```
Queue                 CPU    Mem # CPUs
-------------- ---------- ------ ------
small-serial     40:00:00    4gb      1 
medium-serial    90:00:00   16gb      1 
large-serial    240:00:00  120gb      1 
class             2:00:00   64gb   1-64 
daytime           4:00:00   16gb    1-4 
express          01:00:00  500mb      1
```

```
Enter Queue Name (default <cr>: small-serial) <must be a serial queue>
Enter Time Limit (default <cr>: 40:00:00 HH:MM:SS) <enter time limit>
Enter memory limit (default <cr>: 500mb) <enter required memory>
Enter GPU architecture [t10/fermi/kepler/any] (default <cr>: any) <kepler>
```

Your standard out for your job will be written to 
```
<jobname>.o<job number>
```

###To Run a Simulation on the Alabama Supercomputer Center in debug mode:
```
gpu_interactive

What architecture GPU do you want [any,t10,fermi,kepler]: <kepler>
Do you want to use X-windows [y/n]: <n>

cd /path/to/MCGPU/
cd bin/
./metrosim ./[configuration file]

```

For more information, see the Alabama Supercomputer Center manual.

##Visualizing PDB Files

At the end of the simulation, MCGPU produces a PDB file, which can be loaded
into a visualization program to view the box.  Recommended PDB viewers
include:

* [Jmol](http://jmol.sourceforge.net/)

* [Chimera](https://www.cgl.ucsf.edu/chimera/)

* [RasMol](http://www.openrasmol.org/)

##Running Automated Tests
```
cd /path/to/MCGPU/
make          # Build the metrosim binary first (required by metrotest)
make tests    # Then build the metrotest binary
cd bin/
./metrotest   # Will take several minutes
```

##Profiling
### For CPU profiling:
For CPU profiling, build metrosim with profiling enabled, run it (which will
produce a file called gmon.out), and then use gprof to view the resulting
profile data.
```
make BUILD=profile
bin/metrosim -s resources/exampleFiles/indole4000.config   # or another config file
gprof bin/metrosim
```

### For GPU profiling:
For GPU profiling, build metrosim in release mode, and run it using nvprof.
```
make
nvprof --print-gpu-summary bin/metrosim -s resources/exampleFiles/indole4000.config   # or another config file
```
Alternatively, ```nvprof --print-gpu-trace``` will print information about every
kernel launch (so only run this with a very few time steps).

The NVIDIA Visual Profiler, nvvp, is highly recommended and provides much more
detailed information that the nvprof commands above.

##Running With Multiple Solvents
MCGPU currently supports the simulaton of two solvents within one z-matrix file where separate solvents are separated by TERZ.

When using multiple solvents, the primary index array (Configuration File line 30), must contan at least one primary index array for each molecule, with the arrays enclosed in brackets and comma separated. For example [2],[1,3] represents the primary index structure for two molecules where the first molecule (defined above TERZ) has the primary index of '2' and the second molecule (defined below TERZ) has the primary indexes of '1' and '3'.

##Available Command-line Options
 * `--serial (-s)`: Runs simulation on CPU (default)
 * `--parallel (-p)`: Runs simulation on GPU (requries CUDA)
 * `--name <title>`: Specifies the name of the simulation that will be run.
 * `--steps <count> (-n)`: Specifies how many simulation steps to execute in the Monte Carlo Metropolis algorithm. Ignores steps to run in config file, if present (line 10).
 * `--verbose (-k)`: Enables real time energy printouts
 * `--neighbor <interval> (-l)`: Specifies to use the neighborlist structure for molecular organization. interval is optional and refers to how many steps between updating the neighborlist (default is 100).
 * `--status-interval <interval> (-i)`: Specifies the number of simulation steps between status updates.
 * `--state-interval <interval> (-I)`: Specifies the number of simulation steps between state file snapshots of the current simulation run.

To view documentation for all command-line flags available, use the --help flag:
```
./metrosim --help
```

##Configuration File
Configuration files are used to configure a simulation. Command-line options override values given in this file.
Line numbers are important

```
[1]     #line 1 is a comment and will be ignored
[2]     <x dimension>
[3]     <y dimension>
[4]     <z dimension>
[5]     #line 5 is a comment and will be ignored
[6]     <temperature in kelvin>
[7]     #line seven is a comment and will be ignored
[8]     <max translation for a molecule in angstroms>
[9]     #line 9 is a comment and is ignored
[10]    <number of steps for which to run the simulation> 
[11]    #line 11 is a comment and is ignored
[12]    <number of molecules in the simulation>
[13]    #line 13 is a comment and is ignored
[14]    <path to the opls-aa.par file>
[15]    #line 15 is a comment and is ignored
[16]    <path to the z-matrix file to use>
[17]    #line 17 is a comment and is ignored.
[18]    <path to the state input file; overrides z-matrix setting if present>
[19]    #line 19 is a comment and is ignored.
[20]    <path to the state output file>
[21]    #line 21 is a comment and is ignored.
[22]    <path to the pdb output file>
[23]    #line 23 is a comment and is ignored.
[24]    <nonbonded cutoff distance in angstroms>
[25]    #line 25 is a comment and is ignored.
[26]    <max rotation for a molecule in degrees>
[27]    #line 27 is a comment and is ignored.
[28]    <random number seed input as integer value>
[29]    #line 29 is a comment and is ignored.
[30]    <primary atom index array to be used during cutoff as integer indexes of z-matrix atom in molecule, comma separated, starting from zero>
```

**Contributing Authors**: Guillermo Aguirre, Scott Aldige, James Bass, Jared Brown, Matt Campbell, William Champion, Nathan Coleman, Yitong Dai, Seth Denney, Matthew Hardwick, Andrew Lewis, Alexander Luchs, Jennifer Lynch, Tavis Maclellan, Joshua Mosby, Jeffrey Overbey, Mitchell Price, Robert Sanek, Jonathan Sligh, Riley Spahn, Kalan Stowe, Ashley Tolbert, Albert Wallace, Jay Whaley, Seth Wooten, James Young, Francis Zayek, Xiao (David) Zhang, and Orlando Acevedo*

**Software License**:
MCGPU. Computational Chemistry: Highly Parallel Monte Carlo Simulations on CPUs and GPUs.
Copyright (C) 2016  Orlando Acevedo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. <http://www.gnu.org/licenses/>
