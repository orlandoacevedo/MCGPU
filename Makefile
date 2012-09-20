NV=nvcc

SRC=src/
TST=test/
UTIL=Utilities/
PARA=parallelMetropolis/
LIN=LinearMetropolis/
BIN=bin/

FLAGS=-arch sm_20 -lcurand -g -c
FLAGSALT=-arch sm_20 -lcurand -g

TSTEXE=parallelTest
EXE=parallelExample
LINEXE=linearExample
UTILTESTEXE=utilTest

CommonUtilHeader=$(UTIL)$(SRC)geometricUtil.h $(UTIL)$(SRC)metroUtil.h
CommonUtilOBJ=$(UTIL)$(SRC)geometricUtil.o $(UTIL)$(SRC)metroUtil.o

ScannerHeader=$(UTIL)$(SRC)Zmatrix_Scan.h $(UTIL)$(SRC)Opls_Scan.h $(UTIL)$(SRC)State_Scan.h $(UTIL)$(SRC)Config_Scan.h
ScannerOBJ=$(UTIL)$(SRC)State_Scan.o $(UTIL)$(SRC)Config_Scan.o $(UTIL)$(SRC)Opls_Scan.o $(UTIL)$(SRC)Zmatrix_Scan.o

linearSimHeader= $(CommonUtilHeader) $(ScannerHeader) $(LIN)$(SRC)metroLinearUtil.h  $(PARA)$(TST)baseTests.h
linearSimOBJ=$(CommonUtilOBJ) $(ScannerOBJ) $(LIN)$(SRC)linearSimulationExample.o  $(LIN)$(SRC)metroLinearUtil.o  $(PARA)$(TST)baseTests.o

metroSimHeader= $(CommonUtilHeader) $(ScannerHeader) $(PARA)$(SRC)metroCudaUtil.cuh
metroSimOBJ=$(CommonUtilOBJ) $(ScannerOBJ) $(PARA)$(SRC)metropolisSimulationExample.o  $(PARA)$(SRC)metroCudaUtil.o  $(PARA)$(TST)baseTests.o

TestsHeader= $(CommonUtilHeader) $(PARA)$(TST)parallelTest.cuh $(PARA)$(TST)copyTests.cuh	$(PARA)$(TST)baseTests.h	\
	$(PARA)$(SRC)metroCudaUtil.cuh 
TestsOBJ=$(CommonUtilOBJ)	$(PARA)$(TST)parallelTest.o	$(PARA)$(TST)copyTests.o $(PARA)$(TST)baseTests.o	\
	 $(PARA)$(SRC)metroCudaUtil.o 
	 
UtiltestHeader=$(CommonUtilHeader) $(ScannerHeader) $(UTIL)$(TST)stateTest.h $(UTIL)$(TST)configurationTest.h $(UTIL)$(TST)zMatrixTest.h
UtiltestOBJ=$(CommonUtilOBJ) $(ScannerOBJ) $(UTIL)$(TST)stateTest.o $(UTIL)$(TST)configurationTest.o $(UTIL)$(TST)zMatrixTest.o $(UTIL)$(TST)utilityTests.o

all: dir tests utilTests metroSim linearSim

metroSim:$(metroSimHeader) $(metroSimOBJ)
	$(NV) $(FLAGSALT) $(metroSimOBJ) -o $(BIN)$(EXE) 

linearSim: $(linearSimHeader) $(linearSimOBJ)
	$(NV) $(FLAGSALT) $(linearSimOBJ) -o $(BIN)$(LINEXE) 
	
tests: $(TestsHeader) $(TestsOBJ)
	$(NV) $(FLAGSALT) $(TestsOBJ) -o $(BIN)$(TSTEXE) 
	
utilTests: $(UtiltestOBJ)
	$(NV) $(UtiltestOBJ) -o $(BIN)$(UTILTESTEXE)

.cpp.o:
	$(NV) $(FLAGS) $< -o $@
	
%.o : %.cu
	$(NV) $(FLAGS) $< -o $@
	
dir:
	mkdir -p $(BIN)

clean:
	rm -f $(UTIL)$(SRC)*.o
	rm -f $(PARA)$(SRC)*.o
	rm -f $(LIN)$(SRC)*.o
	rm -f $(UTIL)$(TST)*.o
	rm -f $(PARA)$(TST)*.o
	rm -f $(BIN)$(UTILTESTEXE)
	rm -f $(BIN)$(EXE)
	rm -f $(BIN)$(TSTEXE)
	rm -f $(BIN)$(LINEXE)
