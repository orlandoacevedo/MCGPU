CC=gcc
NVCC=nvcc

UTILSRC=Utilities/src/
UTILTST=Utilities/test/
LINSRC=LinearMetropolis/src/
PARASRC=parallelMetropolis/src/
PARATST=parallelMetropolis/test/
CLASSSRC=src/
TESTSRC=test/
BIN=bin/

CFLAGS=-g -pg -c
FLAGSALT=-lstdc++ -pg -g

CUFLAGS=-arch sm_20 -g -c

TESTCLASS=linearSim
PARASIM=parallelSim

TSTEXE=parallelTest
EXE=parallelExample
LINEXE=linearExample
UTILTESTEXE=utilTests

CommonUtilOBJ=$(UTILSRC)geometricUtil.o $(UTILSRC)metroUtil.o

ScannerOBJ=$(UTILSRC)State_Scan.o $(UTILSRC)Config_Scan.o $(UTILSRC)Opls_Scan.o $(UTILSRC)Zmatrix_Scan.o

TESTCLASSOBJ=$(CommonUtilOBJ)	$(ScannerOBJ) $(LINSRC)SimBox.o	$(LINSRC)linearSim.o $(LINSRC)linearSimulationExample.o 

ParallelOBJS=$(CommonUtilOBJ)	$(ScannerOBJ) $(LINSRC)SimBox.o	$(PARASRC)parallelSim.o $(PARASRC)parallelSimulationExample.o $(PARASRC)GPUSimBox.o $(PARASRC)parallelUtil.o
	 
all: $(TESTCLASS) ${PARASIM}

$(TESTCLASS):$(TESTCLASSOBJ)
	$(CC) $(FLAGSALT) $(TESTCLASSOBJ) -o $(BIN)/$(TESTCLASS) 
${PARASIM}:${ParallelOBJS}
	$(NVCC) $(FLAGSALT) $(ParallelOBJS) -o $(BIN)/$(PARASIM)

-include $(OBJS:.o=.d) 
	
%.d:%.cpp
	@echo "Generate depend file for " $<;                      \
	dirname $< >$@.p.$$$$;echo "\\" >>$@.p.$$$$;tr -d "\n" <$@.p.$$$$ >$@.$$$$;     \
	$(NV) -M $(CPPFLAGS) $< >> $@.$$$$;                      \
  sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;     \
  rm -f $@.$$$$ $@.p.$$$$
  
%.d:%.cu
	@echo "Generate depend file for " $<;                      \
	dirname $< >$@.p.$$$$;echo "\\" >>$@.p.$$$$;tr -d "\n" <$@.p.$$$$ >$@.$$$$;     \
	$(CC) -M $(CPPFLAGS) $< >> $@.$$$$;                      \
  sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;     \
  rm -f $@.$$$$ $@.p.$$$$

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
%.o:%.cu
	$(NVCC) $(CUFLAGS) $< -o $@
	
dir:
	mkdir -p $(BIN)

clean:
	rm -f $(PARASRC)*.o $(LINSRC)*.o $(TESTSRC)*.o 
	rm -f $(BIN)/$(TESTCLASS) $(BIN)/$(PARASIM)
