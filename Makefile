NV=nvcc

UTILSRC=Utilities/src/
UTILTST=Utilities/test/
LINSRC=LinearMetropolis/src/
PARASRC=parallelMetropolis/src/
PARATST=parallelMetropolis/test/
BIN=bin/

FLAGS=-arch sm_20 -lcurand -g -c
FLAGSALT=-arch sm_20 -lcurand -g

TSTEXE=parallelTest
EXE=parallelExample
LINEXE=linearExample
UTILTESTEXE=utilTest

CommonUtilOBJ=$(UTILSRC)geometricUtil.o $(UTILSRC)metroUtil.o

ScannerOBJ=$(UTILSRC)State_Scan.o $(UTILSRC)Config_Scan.o $(UTILSRC)Opls_Scan.o $(UTILSRC)Zmatrix_Scan.o

linearSimOBJ=$(CommonUtilOBJ) $(ScannerOBJ) $(LINSRC)linearSimulationExample.o  $(LINSRC)metroLinearUtil.o  $(PARATST)baseTests.o

metroSimOBJ=$(CommonUtilOBJ) $(ScannerOBJ) $(PARASRC)metropolisSimulationExample.o  $(PARASRC)metroCudaUtil.o  $(PARATST)baseTests.o

TestsOBJ=$(CommonUtilOBJ)	$(PARATST)parallelTest.o	$(PARATST)copyTests.o $(PARATST)baseTests.o	\
	 $(PARASRC)metroCudaUtil.o 
	 
UtiltestOBJ=$(CommonUtilOBJ) $(ScannerOBJ) $(UTILTST)stateTest.o $(UTILTST)configurationTest.o $(UTILTST)zMatrixTest.o \
	$(UTILTST)utilityTests.o $(UTILTST)geometricTest.o 

all: dir $(BIN)$(TSTEXE)  $(BIN)$(UTILTESTEXE) $(BIN)$(EXE) $(BIN)$(LINEXE)

$(BIN)$(EXE):$(metroSimOBJ)
	$(NV) $(FLAGSALT) $(metroSimOBJ) -o $(BIN)$(EXE) 

$(BIN)$(LINEXE): $(linearSimOBJ)
	$(NV) $(FLAGSALT) $(linearSimOBJ) -o $(BIN)$(LINEXE) 
	
$(BIN)$(TSTEXE) : $(TestsOBJ)
	$(NV) $(FLAGSALT) $(TestsOBJ) -o $(BIN)$(TSTEXE) 
	
$(BIN)$(UTILTESTEXE): $(UtiltestOBJ)
	$(NV) $(UtiltestOBJ) -o $(BIN)$(UTILTESTEXE)

OBJS=$(CommonUtilOBJ) $(ScannerOBJ) $(metroSimOBJ) $(linearSimOBJ) $(TestsOBJ) $(UtiltestOBJ)
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
	$(NV) -M $(CPPFLAGS) $< >> $@.$$$$;                      \
  sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;     \
  rm -f $@.$$$$ $@.p.$$$$

.cpp.o:
	$(NV) $(FLAGS) $< -o $@
	
%.o:%.cu
	$(NV) $(FLAGS) $< -o $@
	
dir:
	mkdir -p $(BIN)

clean:
	rm -f $(UTILSRC)*.o $(PARASRC)*.o $(LINSRC)*.o
	rm -f $(UTILTST)*.o $(PARATST)*.o

	rm -f $(UTILSRC)*.d $(PARASRC)*.d $(LINSRC)*.d
	rm -f $(UTILTST)*.d $(PARATST)*.d

	rm -f $(BIN)$(UTILTESTEXE)
	rm -f $(BIN)$(EXE)
	rm -f $(BIN)$(TSTEXE)
	rm -f $(BIN)$(LINEXE)
