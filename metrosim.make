# Relative search paths for Include Files
IncPaths := src test .

# Automatically generate the compiler flags for included search paths.
Include := $(addprefix -I,$(IncPaths))

# Compiler specific flags for the C++ compiler when generating .o files
# and when generating .d files for dependency information
CxxFlags := -c -g -pg

# Compiler specific flags for the CUDA compiler when generating .o files
# and when generating .d files for dependency information
CuFlags := -c -g -G -pg -arch=sm_35 -rdc=true

# Linker specific flags when the compiler is linking the executable
LFlags := -g -pg

# Linker specific flags for the CUDA compiler
LCuFlags := -g -G -pg -arch=sm_35 -rdc=true

Single := -DUSE_SINGLE_PRECISION

AppName := metrosim

all : metrosim_double

metrosim_double : build
	g++ src/Applications/Application.cpp $(Include) $(CxxFlags) -o obj/src/Applications/Application.o
	g++ src/Metropolis/Box.cpp $(Include) $(CxxFlags) -o obj/src/Metropolis/Box.o
	g++ src/Applications/CommandParsing.cpp $(Include) $(CxxFlags) -o obj/src/Applications/CommandParsing.o
	g++ src/Metropolis/Simulation.cpp $(Include) $(CxxFlags) -o obj/src/Metropolis/Simulation.o
	g++ src/Metropolis/Box.cpp $(Include) $(CxxFlags) -o obj/src/Metropolis/Box.o
	g++ src/Metropolis/SerialSim/SerialBox.cpp $(Include) $(CxxFlags) -o obj/src/Metropolis/SerialSim/SerialBox.o
	g++ src/Metropolis/SerialSim/SerialCalcs.cpp $(Include) $(CxxFlags) -o obj/src/Metropolis/SerialSim/SerialCalcs.o
	g++ src/Metropolis/Utilities/StructLibrary.cpp $(Include) $(CxxFlags) -o obj/src/Metropolis/Utilities/StructLibrary.o
	g++ src/Metropolis/Utilities/MathLibrary.cpp $(Include) $(CxxFlags) -o obj/src/Metropolis/Utilities/MathLibrary.o
	g++ src/Metropolis/Utilities/FileUtilities.cpp $(Include) $(CxxFlags) -o obj/src/Metropolis/Utilities/FileUtilities.o
	nvcc src/Metropolis/Utilities/DeviceQuery.cu $(Include) $(CuFlags) -o obj/src/Metropolis/Utilities/DeviceQuery.o
	nvcc src/Metropolis/ParallelSim/ParallelBox.cu $(Include) $(CuFlags) -o obj/src/Metropolis/ParallelSim/ParallelBox.o
	nvcc src/Metropolis/ParallelSim/ParallelCalcs.cu $(Include) $(CuFlags) -o obj/src/Metropolis/ParallelSim/ParallelCalcs.o
	nvcc -o bin/$(AppName) $(LCuFlags) obj/src/Applications/Application.o obj/src/Applications/CommandParsing.o \
								  obj/src/Metropolis/Simulation.o obj/src/Metropolis/SerialSim/SerialBox.o \
								  obj/src/Metropolis/Utilities/FileUtilities.o obj/src/Metropolis/Utilities/StructLibrary.o \
								  obj/src/Metropolis/Utilities/MathLibrary.o obj/src/Metropolis/Utilities/DeviceQuery.o \
								  obj/src/Metropolis/Box.o obj/src/Metropolis/SerialSim/SerialCalcs.o \
								  obj/src/Metropolis/ParallelSim/ParallelBox.o obj/src/Metropolis/ParallelSim/ParallelCalcs.o

metrosim_single : build
	g++ src/Applications/Application.cpp $(Single) $(Include) $(CxxFlags) -o obj/src/Applications/Application.o
	g++ src/Metropolis/Box.cpp $(Single) $(Include) $(CxxFlags) -o obj/src/Metropolis/Box.o
	g++ src/Applications/CommandParsing.cpp $(Single) $(Include) $(CxxFlags) -o obj/src/Applications/CommandParsing.o
	g++ src/Metropolis/Simulation.cpp $(Single) $(Include) $(CxxFlags) -o obj/src/Metropolis/Simulation.o
	g++ src/Metropolis/Box.cpp $(Single) $(Include) $(CxxFlags) -o obj/src/Metropolis/Box.o
	g++ src/Metropolis/SerialSim/SerialBox.cpp $(Single) $(Include) $(CxxFlags) -o obj/src/Metropolis/SerialSim/SerialBox.o
	g++ src/Metropolis/SerialSim/SerialCalcs.cpp $(Single) $(Include) $(CxxFlags) -o obj/src/Metropolis/SerialSim/SerialCalcs.o
	g++ src/Metropolis/Utilities/IOUtilities.cpp $(Single) $(Include) $(CxxFlags) -o obj/src/Metropolis/Utilities/IOUtilities.o
	g++ src/Metropolis/Utilities/StructLibrary.cpp $(Single) $(Include) $(CxxFlags) -o obj/src/Metropolis/Utilities/StructLibrary.o
	g++ src/Metropolis/Utilities/MathLibrary.cpp $(Single) $(Include) $(CxxFlags) -o obj/src/Metropolis/Utilities/MathLibrary.o
	g++ src/Metropolis/Utilities/FileUtilities.cpp $(Include) $(CxxFlags) -o obj/src/Metropolis/Utilities/FileUtilities.o
	nvcc src/Metropolis/Utilities/DeviceQuery.cu $(Single) $(Include) $(CuFlags) -o obj/src/Metropolis/Utilities/DeviceQuery.o
	nvcc src/Metropolis/ParallelSim/ParallelBox.cu $(Single) $(Include) $(CuFlags) -o obj/src/Metropolis/ParallelSim/ParallelBox.o
	nvcc src/Metropolis/ParallelSim/ParallelCalcs.cu $(Single) $(Include) $(CuFlags) -o obj/src/Metropolis/ParallelSim/ParallelCalcs.o
	nvcc -o bin/$(AppName) $(Single) $(LCuFlags) obj/src/Applications/Application.o obj/src/Applications/CommandParsing.o \
								  obj/src/Metropolis/Simulation.o obj/src/Metropolis/SerialSim/SerialBox.o \
								  obj/src/Metropolis/Utilities/FileUtilities.o obj/src/Metropolis/Utilities/StructLibrary.o \
								  obj/src/Metropolis/Utilities/MathLibrary.o obj/src/Metropolis/Utilities/DeviceQuery.o
								  obj/src/Metropolis/Box.o obj/src/Metropolis/SerialSim/SerialCalcs.o \
								  obj/src/Metropolis/ParallelSim/ParallelBox.o obj/src/Metropolis/ParallelSim/ParallelCalcs.o

build :
	@mkdir -p bin obj obj/src obj/src/Applications obj/src/Metropolis obj/src/Metropolis/SerialSim obj/src/Metropolis/ParallelSim obj/src/Metropolis/Utilities