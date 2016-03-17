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

PROG = rpower
PROG_OBJ = bin/mainrpower.o bin/utilities.o bin/mwf.o bin/power.o


all: bin/$(PROG)

bin/%.o: src/%.c
	@echo compiling $*.c with $(CCC) $(CCCFLAGS)
	@$(CCC) $(CCCFLAGS) -c $< -o $@

bin/$(PROG): $(PROG_OBJ)
	$(CCC) $(CCCFLAGS) -o bin/$(PROG) $(PROG_OBJ) $(CCCLNFLAGS)

clean:
	rm bin/*
	
