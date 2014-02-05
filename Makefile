#######################################
# MCGPU Makefile
# Version 1.0
# Authors: Tavis Maclellan
# 
# This Makefile depends on GNU make.
########################
# Makefile Target List #
########################
#
# make [all] : Creates all output files 
# make clean : Cleans out the objects folder
# make remove : Cleans out the objects and bin folder
# make compile : Compiles all the code, but no programs are linked
# make recompile : Deletes the object folder and forces the computer
#      to recompile the source code
#
#############################
# Project Structure Details #
#############################
#
# Defines the main folder names in the project root directory:
#
#   SourceDir : Contains the source code and header files
#   TestDir : Contains the testing code and input files
#   ObjDir : This folder is created by the Makefile and will be
#            populated with compiled object files and binaries
#   BinDir : This folder is created by the Makefile and will
#            hold the generated program executable files
SourceDir := src
TestDir := test
ObjDir := obj
BinDir := bin

# Defines the modules that exist in the source directory that should be
# included by the compiler. All files within these directories that are
# valid file types for the Makefile to handle will automatically be
# compiled into object files in the object directory. Make sure to add
# any new modules to this list, or else they will not be compiled.
Modules := LinearMetropolis ParallelMetropolis Utilities

##############################
# Compiler Specific Settings #
##############################

# Defines the compilers used to compile and link the source files.
# CC will be used to compile C++ files, and NVCC will be used to
# compile CUDA files.
CC := g++
NVCC := nvcc

# Defines the types of files that the Makefile knows how to compile
# and link. Specify the filetype by using a modulus (percent sign),
# followed by a dot and the file extension (e.g. %.java, %.txt).
FileTypes := %.cpp %.cu

# Relative search paths for Include Files
IncPaths := $(SourceDir) $(TestDir) .

# Compiler specific flags for the C++ compiler when generating .o files
# and when generating .d files for dependency information
CxxFlags := -c -g -pg
CxxDepFlags := -MM

# Compiler specific flags for the CUDA compiler when generating .o files
# and when generating .d files for dependency information
CuFlags := -c -g -arch sm_20
CuDepFlags := -M

# Linker specific flags when the compiler is linking the executable
LFlags := -g

###################
# Program Outputs #
###################

# The names of the programs to generate executable output.
LinearSim := linearSim
ParallelSim := parallelSim

# The modules in the source directory that each program is linked with. If 
# the module name is given, then the program will be linked with all of the
# source files in that folder. If the program only needs a specific object
# file from the module, then specify the specific source file needed with
# a .o extension added to the end of the source file name.
LinearSimModules := LinearMetropolis Utilities
ParallelSimModules := 	ParallelMetropolis Utilities \
			LinearMetropolis/SimBox.o

#############################
# Automated Testing Details #
#############################

# The name of the unit test executable program.
UnitTestProg := unittests

# The relative path to the testing module containing the unit test source.
UnitTestDir := $(TestDir)/unittests

# The relative path to the Google Test module that contains the source
# code and libraries for the Google Test framework.
GTestDir := $(TestDir)/gtest-1.7.0

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTestHeaders = $(GTestDir)/include/gtest/*.h \
               $(GTestDir)/include/gtest/internal/*.h

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
GTestFlags := -isystem $(GTestDir)/include $(Include)
GTestFlags += -g -pthread #-Wall -Wextra

# Builds gtest.a and gtest_main.a.
# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTestDir)/src/*.cc $(GTestDir)/src/*.h $(GTestHeaders)

######################
# Internal Variables #
######################

# Derives the paths to each of the source modules.
SourceModules := $(addprefix $(SourceDir)/,$(Modules))

# Creates a list of folders inside the object output directory that need
# to be created for the compiled files.
ObjFolders := $(addprefix $(ObjDir)/,$(SourceModules))

# Searches through the specified Modules list for all of the valid
# files that it can find and compile. Once all of the files are 
# found, they are appended with an .o and prefixed with the object
# directory path. This allows the compiled object files to be routed
# to the proper output directory.
Sources := $(filter $(FileTypes),$(wildcard $(addsuffix /*,$(SourceModules))))
Objects := $(patsubst %,$(ObjDir)/%.o,$(basename $(Sources)))

## The unit testing objects are all gathered seperately because they are 
## included all at once from the testing directory and are compiled into the
## output program alongside the source objects.
#UnitTestingSources := $(filter $(FileTypes),$(wildcard $(UnitTestDir)/*))
#UnitTestingObjects := $(patsubst %,$(ObjDir)/%.o,\
#		      $(basename $(UnitTestingSources)))

# Contains the list of directories to be added into the include search
# path used to located included header files.
Include := $(addprefix -I,$(IncPaths))

# This is the directory and matching criteria that the dependency files use
# to figure out if files have been updated or not.
DepDir := $(ObjDir)
DF = $(DepDir)/$*

#################
# Template Code #
#################

# Name: find-objects
# Description: Search through all of the object files that are to be compiled,
#	and finds all object files that exist inside the given folder.
# Parameters:
#	folder : This variable should be set inside the foreach loop used
#	         to call this function.
find-objects = $$(filter $$(folder)/%.o,$$(Objects))


# Name: program-template
# Description: Using the given program name and modules list, return
# 	a formatted target that links the program with its objects and
#	saves it to the bin directory.
# Remarks: This function is finding all of the object files that need
#	to be linked with the program, and it uses the modules list specified
#	at the top of the Makefile to do so. The final list is sorted in
#	order to remove duplicate values only (we don't actually need to
#	sort anything).
# Parameters (Called values):
# 	(1) : The name of the executable program. This value should match
#	      the value in the Programs variable.
#	(2) : The list of module names that the program needs to link with
#	(3) : The compiler to use to link the program.
define program-template
$(1)_LIST := $(addprefix $(ObjDir)/$(SourceDir)/,$(2))
$(1)_MODULE_LIST := $$(filter-out %.o,$$($(1)_LIST))
$(1)_OBJS := $$(filter %.o,$$($(1)_LIST))
$(1)_OBJS += $$(foreach folder,$$($(1)_MODULE_LIST),$(find-objects))

$(1): $$(sort $$($(1)_OBJS))
	$(3) $(LFlags) $$^ -o $(BinDir)/$$@
	
endef

##############################
# Makefile Rules and Targets #
##############################

# Specifies that these make targets are not actual files and therefore will
# not break if a similar named file exists in the directory.
.PHONY : all build compile recompile clean remove

# The house-keeping targets that provide the basic functionality of the 
# Makefile.

all : build $(LinearSim) $(ParallelSim)

compile : build $(Objects)

recompile : remove build $(Objects)

build :
	@mkdir -p $(ObjDir) $(ObjFolders) $(BinDir)

clean : 
	rm -rf $(ObjDir)

remove :
	rm -rf $(ObjDir) $(BinDir)

# The targets and rules used to build the program output files that are
# placed in the bin directory. The target and rules are being generated
# by a function that takes the program name, the program modules list, and
# the linker as arguments.

$(eval $(call program-template,$(LinearSim),$(LinearSimModules), $(CC)))

$(eval $(call program-template,$(ParallelSim),$(ParallelSimModules), $(NVCC)))


# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
$(ObjDir)/gtest-all.o : $(GTEST_SRCS_)
	$(CC) $(GTestFlags) -I$(GTestDir) -c \
            $(GTestDir)/src/gtest-all.cc -o $@

$(ObjDir)/gtest_main.o : $(GTEST_SRCS_)
	$(CC) $(GTestFlags) -I$(GTestDir) -c \
	  $(GTestDir)/src/gtest_main.cc -o $@

$(ObjDir)/gtest.a : $(ObjDir)/gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

$(ObjDir)/gtest_main.a : $(ObjDir)/gtest-all.o $(ObjDir)/gtest_main.o
	$(AR) $(ARFLAGS) $@ $^


# Here are the Rules that determine how to compile a CUDA and a C++ source 
# file into an object file. Also, this rule will generate the file's
# dependecies and format the file into a format that allows for easy and
# effecient dependency resolution. The CUDA code must be compiled twice (once
# for the object file and once for the dependencies), whereas the C++ code
# can accomplish both actions with one compile (by using the -MMD flag).

$(ObjDir)/%.o : %.cu
	$(NVCC) $(Include) -M $< -o $(ObjDir)/$*.d
	@cp $(DF).d $(DF).P
	@sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	    -e '/^$$/ d' -e 's/$$/ :/' < $(DF).d >> $(DF).P
	@rm -f $(DF).d
	$(NVCC) $(CuFlags) $(Include) $< -o $@

$(ObjDir)/%.o : %.cpp
	$(CC) $(CxxFlags) $(Include) -MMD $< -o $@
	@cp $(DF).d $(DF).P
	@sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	    -e '/^$$/ d' -e 's/$$/ :/' < $(DF).d >> $(DF).P
	@rm -f $(DF).d

######################
# Dependency Include #
######################

# This conditional statement will attempt to include all of the dependency
# files located in the object directory. If the files exist, then their
# dependency information is loaded, and each source file checks to see if
# it needs to be recompiled. The if statements are used to make sure that
# the dependency info isn't rebuilt when the object directory is being
# cleaned.

ifneq ($(MAKECMDGOALS),remove)
ifneq ($(MAKECMDGOALS),clean)
-include $(Objects:.o=.P)
endif
endif

#####################
# TODO: Future Work #
#####################

# While looking at several example Makefile resources when building Version
# 1.0 I came across an implementation that split the output programs into
# either Debug or Release mode. This is definitely a feature that we should
# probably implement in the future. The main difference is a different set
# of compiler flags are used when compiling and linking the source files and
# object files. These flags specify how much to optimize the code and whether
# or not to include debugging information.
#
# The Debug mode will compile faster and allow for debugger and profilers to
# monitor the executable program, but there is a performance penalty. The
# Release mode, on the other hand, is highly optimized to run very fast at
# the cost of compile time and debuggin/profiling support.
#
# This would not take too long to implement, although I believe that a
# template method would work best that would generate the Makefile rules
# for each mode (Debug and Release) for each executable output file. See
# the eval() and call() functions for more detail on dynamic generation
# of Makefile rules.
#
# - (Tavis Maclellan) : 2014-02-01
