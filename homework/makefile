OS = linux
#OS = windows

GCCOPT = gcc -O2 -funroll-loops
GCCDEB = gcc -ggdb -Wall -pedantic

#CCC = $(GCCOPT)
CCC = $(GCCDEB)

CFLAGS=

ifeq ($(OS),linux)
  LIBS = -lm -lpthread
else ifeq ($(OS),windows)
  LIBS = -lpthreadGC2
else
  LIBS = 
endif


CCCLNFLAGS = $(LIBS)

CCCFLAGS = 

RPOWER = rpower
RPOWER_OBJ = bin/mainrpower.o bin/utilities.o bin/mwf.o bin/power.o

MYOPT = myopt
MYOPT_OBJ = bin/myopt.o bin/myobasic.o bin/myoalgo.o bin/utilities.o 

COVPCA = covpca
COVPCA_OBJ = bin/pca.o bin/utilities.o bin/mwf.o bin/power.o

all: bin/$(RPOWER) bin/$(MYOPT) bin/$(COVPCA)

bin/%.o: src/%.c
	@echo compiling $*.c with $(CCC) $(CCCFLAGS)
	@$(CCC) $(CCCFLAGS) -c $< -o $@

bin/$(RPOWER): $(RPOWER_OBJ)
	$(CCC) $(CCCFLAGS) -o bin/$(RPOWER) $(RPOWER_OBJ) $(CCCLNFLAGS)

bin/$(MYOPT): $(MYOPT_OBJ)
	$(CCC) $(CCCFLAGS) -o bin/$(MYOPT) $(MYOPT_OBJ) $(CCCLNFLAGS)

bin/$(COVPCA): $(COVPCA_OBJ)
	$(CCC) $(CCCFLAGS) -o bin/$(COVPCA) $(COVPCA_OBJ) $(CCCLNFLAGS)

clean:
	rm bin/*

