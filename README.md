Computational Chemistry: Highly Parallel Monte Carlo Simulations on CPUs and GPUs 
MCGPU (Monte Carlo on Graphics Processing Units)
===============================================================

##Requirements
###Required Hardware:
 * For Parallel Execution - Nvidia graphics card with compute Compute Capability 2.0 (or greater)
    * Tested on Nvidia Fermi M2070, Kepler K20m


###Required Software:
 * Linux or OSX Operating System
    * Tested on Ubuntu 14.04 LTS, SUSE Linux Enterprise Server 11 SP2, OSX Yosemite 
 * [Nvidia Developer Toolkit](http://developer.nvidia.com/cuda-downloads)
    * Tested with CUDA 5.5, 6.5
 * [OpenMp](http://www.openmp.org)

*Note*: If you are using the Alabama Supercomputer Center (ASC), configure based on the instructions received when you set up your account.

##Build
```
git clone git://github.com/orlandoacevedo/MCGPU.git
cd MCGPU/
make
```

*Note*: To build on a local machine, use LOCAL_INSTALL=1
```
make LOCAL_INSTALL=1
```

*Note*: To build on an OSX machine, use MAC=1
```
make MAC=1
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
Enter GPU architecture [t10/fermi/kepler/any] (default <cr>: any) <kepler>
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


##Running With Multiple Solvents
MCGPU currently supports the simulaton of two solvents within one z-matrix file where separate solvents are separated by TERZ.

When using multiple solvents, the primary index array (Configuration File line 30), must contan at least one primary index array for each molecule, with the arrays enclosed in brackets and comma separated. For example [2],[0,3] represents the primary index structure for two molecules where the first molecule (defined above TERZ) has the primary index of '2' and the second molecule (defined below TERZ) has the primary indexes of '0' and '3'.


##Available Command-line Options
 * `--serial`: Runs simulation on CPU (default)
 * `--parallel`: Runs simulation on GPU (requries CUDA)
 * `--list-devices`: Lists available CUDA-capable devices (requires CUDA)
 * `--device <index>`: Specifies what device to use when running a simulation. Index refers to one given in --list-devices. (requires CUDA)
 * `--threads <count>`: Specifies number of threads to use when running on a multiprocessor CPU (--serial only)
 * `--name <title>`: Specifies the name of the simulation that will be run.
 * `--steps <count>`: Specifies how many simulation steps to execute in the Monte Carlo Metropolis algorithm. Ignores steps to run in config file, if present (line 10).
 * `--verbose`: Enables real time energy printouts
 * `--neighbor <interval>`: Specifies to use the neighborlist structure for molecular organization. interval is optional and refers to how many steps between updating the neighborlist (default is 100).

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
[30]    <primary atom index array to be used during cutoff as integer indexes of z-matrix atom in molecule, comma separated, starting from zero>
```

**Contributing Authors**: Scott Aldige, Jared Brown, Matt Campbell, William Champion, Nathan Coleman, Yitong Dai, Seth Denney, Matthew Hardwick, Andrew Lewis, Alexander Luchs, Jennifer Lynch, Tavis Maclellan, Joshua Mosby, Robert Sanek, Riley Spahn, Kalan Stowe, Ashley Tolbert, Albert Wallace, Seth Wooten, James Young, Francis Zayek, Xiao (David) Zhang, and Orlando Acevedo*
