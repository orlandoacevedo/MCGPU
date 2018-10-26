MCGPU (Monte Carlo on Graphics Processing Units)
===============================================================

## Requirements
 * [PGI Accelerator C/C++ Compiler with OpenACC](https://www.pgroup.com/resources/accel.htm) *or*
   [OpenACC Toolkit](https://developer.nvidia.com/openacc-toolkit) (free for academic use)
    * *Note*: If you are using the Alabama Supercomputer Center's Dense Memory Cluster (DMC), type ```module load pgi``` to load the PGI compilers.
 * For GPU Execution: NVIDIA CUDA-capable graphics card
    * Tested on NVIDIA Kepler K20m and K40
    * By default, the PGI compilers target Fermi-generation GPUs (Compute Capability 2.0) and higher
 * Linux operating system

## Build
```
git clone git://github.com/orlandoacevedo/MCGPU.git
cd MCGPU/
make
```

*Note*: To build in debug mode (so the executable can be debugged using cuda-gdb), use BUILD=debug:
```
make BUILD=debug
```

*Note*: MCGPU can also be compiled with GCC, though it doesn't support GPU
offloading and is therefore substantially slower. To use GCC, compile with

```
make CC=g++
```

If you compile with GCC, you **cannot** run in parallel mode.

## Run

To Run a Simulation on a Local Machine:
```
cd /path/to/MCGPU/
bin/metrosim [configuration file] [options]
```
where `[configuration file]` is a .config file containing configuration
information and `[options]` are command-line options. An example demo.config
can be found in the resources folder. See below for specific .config file
documentation and all command-line options available.

*Example 1:*
```
bin/metrosim resources/exampleFiles/indole4000.config -k
```
runs a simulation on the GPU if possible and on the CPU otherwise.  The ```-k``` option enables
verbose output: status information will be printed every 1000 time steps.

*Example 2:*
```
bin/metrosim resources/exampleFiles/indole4000.config -p --name indole4000 -n 5000 -i 1000 -k
```
runs a simulation on the GPU (```-p```), for 5000 steps, printing status information every 1000 intervals.

*Example 3:*
```
bin/metrosim resources/exampleFiles/indole4000.config -s --name indole4000 -n 1000 -i 100 -k
```
runs a simulation on the CPU (```-s```), for 1000 steps, printing status information every 100 steps.

## To Run a Simulation on the Alabama Supercomputer Center's DMC:
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

## To Run a Simulation on the Alabama Supercomputer Center in debug mode:
```
gpu_interactive

What architecture GPU do you want [any,t10,fermi,kepler]: <kepler>
Do you want to use X-windows [y/n]: <n>

cd /path/to/MCGPU/
cd bin/
./metrosim ./[configuration file]

```

For more information, see the Alabama Supercomputer Center manual.

## Visualizing PDB Files

At the end of the simulation, MCGPU produces a PDB file, which can be loaded
into a visualization program to view the box.  Recommended PDB viewers
include:

* [Jmol](http://jmol.sourceforge.net/)

* [Chimera](https://www.cgl.ucsf.edu/chimera/)

* [RasMol](http://www.openrasmol.org/)

## Running Automated Tests
```
cd /path/to/MCGPU/
make          # Build the metrosim binary first (required by metrotest)
make tests    # Then build the metrotest binary
cd bin/
./metrotest   # Will take several minutes
```

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

## Running With Multiple Solvents
MCGPU currently supports the simulation of two solvents within one z-matrix file where separate solvents are separated by TERZ.

When using multiple solvents, the primary index array (Configuration File line 30), must contain at least one primary index array for each molecule, with the arrays enclosed in brackets and comma separated. For example [2],[1,3] represents the primary index structure for two molecules where the first molecule (defined above TERZ) has the primary index of '2' and the second molecule (defined below TERZ) has the primary indexes of '1' and '3'.

## Available Command-line Options
 * `--serial (-s)`: Runs simulation on CPU
 * `--parallel (-p)`: Runs simulation on GPU
 * `--name <title>`: Specifies the name of the simulation that will be run.
 * `--steps <count> (-n)`: Specifies how many simulation steps to execute in the Monte Carlo Metropolis algorithm. Ignores steps to run in config file, if present (line 10).
 * `--verbose (-k)`: Enables real time energy printouts
 * `--status-interval <interval> (-i)`: Specifies the number of simulation steps between status updates.
 * `--state-interval <interval> (-I)`: Specifies the number of simulation steps between state file snapshots of the current simulation run.
 * `--strategy <strategy-name> (-S)`: Specifies the energy calculation strategy to utilize. Current options include `brute-force` and `proximity-matrix`

To view documentation for all command-line flags available, use the --help flag:
```
./metrosim --help
```

## Configuration File
Configuration files are used to configure a simulation. They are formatted
according to the [INI](https://en.wikipedia.org/wiki/INI_file) format.
Command-line options override values given in this file. An example
configuration file is shown below.

```
# Name for the simulation
sim-name=MyTestSimulation

# The dimensions of the periodic simulation box (in angstroms)
x=55
y=55
z=55

# Temperature (in Kelvin)
temp=298.15

# Maximum translation for a molecule during the simulation
max-translation=0.15

# Number of steps to run in the simulation
steps=1000

# Number of molecules
molecules=5120

# Path to opla.par file
opla.par=/absolute/path/to/oplsaa.par

# Path to z-matrix file
z-matrix=/absolute/path/to/matrix.z

# Path to input state *directory*
state-input=/absolute/path/to/input/dir

# Path to state output *directory*
state-output=/absolute/path/to/output/dir

# Path to pdb output *directory*
pdb-output=/absolute/path/to/output/dir

# Cutoff distance (in angstroms)
cutoff=25

# Maximum rotation of any particle
max-rotation=15

# Seed for random generator
random-seed=12345

# Primary atom index (integer indexes of z-matrix atom in molecule, comma
# separated, starting from zero)
primary-atom=1

# Strategy for energy calculations
strategy=brute-force
```

The order of the attributes is not significant. Comments can begin with `#` or
`;`. Empty lines will be ignored.

**Contributing Authors**: Guillermo Aguirre, Scott Aldige, James Bass, Jared Brown, Matt Campbell, William Champion, Nathan Coleman, Yitong Dai, Seth Denney, Matthew Hardwick, Andrew Lewis, Alexander Luchs, Jennifer Lynch, Tavis Maclellan, Brandon Morris, Joshua Mosby, Jeffrey Overbey, Mitchell Price, Robert Sanek, Jonathan Sligh, Riley Spahn, Kalan Stowe, Ashley Tolbert, Albert Wallace, Jay Whaley, Seth Wooten, James Young, Francis Zayek, Xiao (David) Zhang, and Orlando Acevedo*

**Funding**: Gratitude is expressed to the National Science Foundation (NSF CHE-1562205) for funding the project.

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
