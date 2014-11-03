Computational Chemistry: Monte Carlo Simulations on CPU and GPU 
===============================================================

##Requirements
###Required Hardware:
 * Nvidia graphics card with compute Compute Capability 2.0 (or greater)
    * Tested on Nvidia Tesla M2070

###Required Software:
 * Linux Operating System
    * Tested on Ubuntu 14.04 LTS, SUSE Linux Enterprise Server 11 SP2
 * [Nvidia Developer Toolkit](http://developer.nvidia.com/cuda-downloads)
    * Tested with CUDA 4.2, 6.5

*Note*: If you are using the Alabama Supercomputer Center (ASC), configure based on the instructions received when you set up your account.

##Build
```
git clone git://github.com/orlandoacevedo/MCGPU.git
cd MCGPU/
make
```

*Note*: To build in debug mode, use BUILD=debug:
```
make BUILD=debug
```

##Run
###To Run a Simulation on a local machine:
```
cd /path/to/MCGPU/
cd bin/
./metrosim ./[configuration file] [options]
```
Where `[configuration file]` is a .config file containing configuration information, and `[options]` are command-line options. An example demo.config can be found in the resources folder. See below for specific .config file documentation and all command-line options available.

###To Run a Simulation on the Alabama Supercomputer Center:
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
Enter GPU architecture [t10/fermi/any] (default <cr>: any) <fermi>
```

You standard out for your job will be written to 
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

##Available Command-line Options
 * `--serial`: Run simulation serially (on CPU; default)
 * `--parallel`: Run simulation in parallel (on GPU; requries CUDA)
 * `--list-devices`: List available CUDA-capable devices (requires CUDA)
 * `--device <index>`: Specify what device to use when running a simulation. Index refers to one given in --list-devices. (requires CUDA)
 * `--name <title>`: Specifies the name of the simulation that will be run.
 * `--steps <count>`: Specifies how many simulation steps to execute in the Monte Carlo Metropolis algorithm. Ignores steps to run in config file, if present (line 10).
 
To view documentation for all command-line flags available, use the --help flag:
```
./metrosim --help
```

##Configuration File
Configuration files are used to configure a simulation. 
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
[11]    #line 11 is a commend and is ignored
[12]    <number of molecules in the simulation>
[13]    #line 13 is a comment and is ignored
[14]    <path to the opls-aa.par file>
[15]    #line 15 is a comment and is ignored
[16]    <path to the z-matrix file to use>
[17]    #line 17 is a comment and is ignored.
[18]    <path to the state input file; overrides z-matrix setting is present>
[19]    #line 19 is a comment and is ignored.
[20]    <path to the stateoutput file>
[21]    #line 21 is a comment and is ignored.
[22]    <path to the pdb output file>
[23]    #line 23 is a comment and is ignored.
[24]    <nonbonded cutoff distance in angstroms>
[25]    #line 25 is a comment and is ignored.
[26]    <max rotation for a molecule in degrees>
[27]    #line 27 is a comment and is ignored.
[28]    <random number seed input as integer value>
[29]    #line 29 is a comment and is ignored.
[30]    <primary atom index to be used during cutoff as integer index of z-matrix atom in molecule>
```

**Contributing Authors**: Scott Aldige, Matt Campbell, William Champion, Matthew Hardwick, Andrew Lewis, Alexander Luchs, Robert Sanek, Riley Spahn, Kalan Stowe, Ashley Tolbert, Seth Wooten, Xiao (David) Zhang, and Orlando Acevedo
