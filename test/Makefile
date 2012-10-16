CC=gcc

UTILSRC=../Utilities/src/
UTILTST=../Utilities/test/
LINSRC=../LinearMetropolis/src/
PARASRC=../parallelMetropolis/src/
PARATST=../parallelMetropolis/test/
CLASSSRC=src/
TESTSRC=test/
BIN=bin/

FLAGS=-g -pg -c
FLAGSALT=-lstdc++ -pg -g

TESTCLASS=testclass
TSTEXE=parallelTest
EXE=parallelExample
LINEXE=linearExample
UTILTESTEXE=utilTests

CommonUtilOBJ=$(UTILSRC)geometricUtil.o $(UTILSRC)metroUtil.o

ScannerOBJ=$(UTILSRC)State_Scan.o $(UTILSRC)Config_Scan.o $(UTILSRC)Opls_Scan.o $(UTILSRC)Zmatrix_Scan.o

TESTCLASSOBJ=$(CommonUtilOBJ)	$(ScannerOBJ) $(CLASSSRC)SimBox.o	$(CLASSSRC)linearSim.o $(CLASSSRC)linearSimulationExample.o 
	 
all: $(TESTCLASS)

$(TESTCLASS):$(TESTCLASSOBJ)
	$(CC) $(FLAGSALT) $(TESTCLASSOBJ) -o $(TESTCLASS) 

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
	$(CC) $(FLAGS) $< -o $@
	
%.o:%.cu
	$(CC) $(FLAGS) $< -o $@
	
dir:
	mkdir -p $(BIN)

clean:
	rm -f $(CLASSSRC)*.o $(TESTSRC)*.o 
	rm -f $(TESTCLASS)
