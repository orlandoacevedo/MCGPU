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
COMMONUTILO=$(BIN)geometricUtil.o $(BIN)metroUtil.o
SCANNERO=$(BIN)State_Scan.o $(BIN)Config_Scan.o $(BIN)Opls_Scan.o $(BIN)Zmatrix_Scan.o

all: dir tests utilTests metroSim linearSim

metroSim: cudaUtil metroUtil stateScan configScan zMatrix OPLSScan baseTest $(PARA)$(SRC)metropolisSimulationExample.cu
	$(NV) $(FLAGS) $(PARA)$(SRC)metropolisSimulationExample.cu -o $(BIN)metropolisSimulationExample.o 
	$(NV) $(FLAGSALT) $(COMMONUTILO) $(SCANNERO) $(BIN)metroCudaUtil.o $(BIN)metropolisSimulationExample.o $(BIN)baseTests.o -o $(BIN)$(EXE) 

linearSim: linearUtil metroUtil stateScan configScan zMatrix OPLSScan baseTest $(LIN)$(SRC)linearSimulationExample.cpp
	$(NV) $(FLAGS) $(LIN)$(SRC)linearSimulationExample.cpp -o $(BIN)linearSimulationExample.o 
	$(NV) $(FLAGSALT) $(COMMONUTILO) $(SCANNERO) $(BIN)metroLinearUtil.o $(BIN)linearSimulationExample.o $(BIN)baseTests.o -o $(BIN)$(LINEXE) 

tests: cudaUtil geoUtil metroUtil copyTest baseTest $(PARA)$(TST)parallelTest.cu $(PARA)$(TST)parallelTest.cuh
	$(NV) $(FLAGSALT) $(COMMONUTILO) $(BIN)copyTests.o $(BIN)baseTests.o $(BIN)metroCudaUtil.o $(PARA)$(TST)parallelTest.cu -o $(BIN)$(TSTEXE)

copyTest: $(PARA)$(TST)copyTests.cuh $(PARA)$(TST)copyTests.cu
	$(NV) $(FLAGS) $(PARA)$(TST)copyTests.cu -o $(BIN)copyTests.o

baseTest: $(PARA)$(TST)baseTests.h $(PARA)$(TST)baseTests.cpp
	$(NV) $(FLAGS) $(PARA)$(TST)baseTests.cpp -o $(BIN)baseTests.o

cudaUtil: metroUtil baseTest $(PARA)$(SRC)metroCudaUtil.cuh $(PARA)$(SRC)metroCudaUtil.cu
	$(NV) $(FLAGS) $(PARA)$(SRC)metroCudaUtil.cu -o $(BIN)metroCudaUtil.o

linearUtil: metroUtil baseTest $(LIN)$(SRC)metroLinearUtil.h $(LIN)$(SRC)metroLinearUtil.cpp
	$(NV) $(FLAGS) $(LIN)$(SRC)metroLinearUtil.cpp -o $(BIN)metroLinearUtil.o

metroUtil: OPLSScan zMatrix geoUtil $(UTIL)$(SRC)metroUtil.h $(UTIL)$(SRC)metroUtil.cpp
	$(NV) $(FLAGS) $(UTIL)$(SRC)metroUtil.cpp -o $(BIN)metroUtil.o

utilTests: stateScan stateTest configTest metroUtil zMatrix OPLSScan geoTest zMatrixTest
	$(NV) $(COMMONUTILO) $(SCANNERO) $(BIN)configurationTest.o $(BIN)stateTest.o $(BIN)geometricTest.o $(BIN)zMatrixTest.o Utilities/test/utilityTests.cpp -o $(BIN)$(UTILTESTEXE)

zMatrix: geoUtil $(UTIL)$(SRC)Zmatrix_Scan.cpp $(UTIL)$(SRC)Zmatrix_Scan.h
	$(NV) $(FLAGS) $(UTIL)$(SRC)Zmatrix_Scan.cpp -o $(BIN)Zmatrix_Scan.o

OPLSScan: $(UTIL)$(SRC)Opls_Scan.cpp $(UTIL)$(SRC)Opls_Scan.h
	$(NV) $(FLAGS) $(UTIL)$(SRC)Opls_Scan.cpp -o $(BIN)Opls_Scan.o

geoUtil: $(UTIL)$(SRC)geometricUtil.cpp $(UTIL)$(SRC)geometricUtil.h
	$(NV) $(FLAGS) $(UTIL)$(SRC)geometricUtil.cpp -o $(BIN)geometricUtil.o

stateTest: metroUtil $(UTIL)$(TST)stateTest.cpp $(UTIL)$(TST)stateTest.h 
	$(NV) $(FLAGS) $(UTIL)$(TST)stateTest.cpp -o $(BIN)stateTest.o

configScan: $(UTIL)$(SRC)Config_Scan.cpp $(UTIL)$(SRC)Config_Scan.h
	$(NV) $(FLAGS) $(UTIL)$(SRC)Config_Scan.cpp -o $(BIN)Config_Scan.o

stateScan: $(UTIL)$(SRC)State_Scan.cpp $(UTIL)$(SRC)State_Scan.h
	$(NV) $(FLAGS) $(UTIL)$(SRC)State_Scan.cpp -o $(BIN)State_Scan.o

configTest: configScan $(UTIL)$(TST)configurationTest.h $(UTIL)$(TST)configurationTest.cpp
	$(NV) $(FLAGS) $(UTIL)$(TST)configurationTest.cpp -o $(BIN)configurationTest.o

geoTest: geoUtil $(UTIL)$(TST)geometricTest.cpp $(UTIL)$(TST)geometricTest.h
	$(NV) $(FLAGS) $(UTIL)$(TST)geometricTest.cpp -o $(BIN)geometricTest.o

zMatrixTest: metroUtil zMatrix OPLSScan $(UTIL)$(TST)zMatrixTest.cpp $(UTIL)$(TST)zMatrixTest.h 
	$(NV) $(FLAGS) $(UTIL)$(TST)zMatrixTest.cpp -o $(BIN)zMatrixTest.o


dir:
	mkdir -p $(BIN)

clean:
	rm -f $(BIN)*.o
	rm -f $(BIN)$(UTILTESTEXE)
	rm -f $(BIN)$(EXE)
	rm -f $(BIN)$(TSTEXE)
	rm -f $(BIN)$(LINEXE)
