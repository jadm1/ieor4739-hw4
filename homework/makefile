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

BDIR = bin
ODIR = obj



RPOWER = $(BDIR)/rpower
RPOWER_OBJ = $(ODIR)/mainrpower.o $(ODIR)/utilities.o $(ODIR)/mwf.o $(ODIR)/power.o

COVPCA = $(BDIR)/covpca
COVPCA_OBJ = $(ODIR)/pca.o $(ODIR)/utilities.o $(ODIR)/mwf.o $(ODIR)/power.o

MYOPT = $(BDIR)/myopt
MYOPT_OBJ = $(ODIR)/myopt.o $(ODIR)/myoalgo.o $(ODIR)/utilities.o 

PFSIMUL = $(BDIR)/pfsimul
PFSIMUL_OBJ = $(ODIR)/utilities.o $(ODIR)/portfoliosimul.o $(ODIR)/mwf.o



all: $(RPOWER) $(COVPCA) $(MYOPT) $(PFSIMUL)

$(ODIR)/%.o: src/%.c
	@echo compiling $*.c with $(CCC) $(CCCFLAGS)
	@$(CCC) $(CCCFLAGS) -c $< -o $@

$(RPOWER): $(RPOWER_OBJ)
	$(CCC) $(CCCFLAGS) -o $(RPOWER) $(RPOWER_OBJ) $(CCCLNFLAGS)

$(COVPCA): $(COVPCA_OBJ)
	$(CCC) $(CCCFLAGS) -o $(COVPCA) $(COVPCA_OBJ) $(CCCLNFLAGS)

$(MYOPT): $(MYOPT_OBJ)
	$(CCC) $(CCCFLAGS) -o $(MYOPT) $(MYOPT_OBJ) $(CCCLNFLAGS)

$(PFSIMUL): $(PFSIMUL_OBJ)
	$(CCC) $(CCCFLAGS) -o $(PFSIMUL) $(PFSIMUL_OBJ) $(CCCLNFLAGS)

clean:
	rm $(RPOWER) $(COVPCA) $(MYOPT) $(PFSIMUL) -f
	rm $(ODIR)/* -f

