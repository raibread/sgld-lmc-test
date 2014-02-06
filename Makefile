#   Makefile Langevin SGD bias tests

#   The author gives permission for anyone to use this publicly posted 
#   code for any purpose.  The code was written for teaching, not research 
#   or commercial use.  It has not been tested thoroughly and probably has
#   serious bugs.  Results may be inaccurate, incorrect, or just wrong.

#   Author: Raiden Hasegawa
#   E-mail: raiden [dot] hasegawa [at] gmail [dot] com

SUFFIXES :=
%.c:
%.cpp:
%.o:
%.h:
%.py:

#        Compilers, linkers, compile options, link options 

CC      = gcc             #  C
CPP     = g++             #  C++ 

CCFLAGS = -O2

OTHERFLAGS = -std=c++11   

LINKER  =  g++       # linker, use C++ linker to link C and C++ object files.

#          Lists of files

C_SOURCES   = rng.c 

CPP_SOURCES = main.cpp toyModel.cpp

#P_SOURCES   = histogram.py 

C_OBJECTS   = $(patsubst %.c,   %.o, $(C_SOURCES) )
CPP_OBJECTS = $(patsubst %.cpp, %.o, $(CPP_SOURCES) )


#           Stuff that every single program depends on

FOR_ALL   = Makefile header.h

ALL_SOURCES = $(CPP_SOURCES) $(C_SOURCES) $(P_SOURCES) Makefile header.h README.md


#        compiling instructions, all manual, no implicit rules

main.o: main.cpp $(FOR_ALL)
	$(CPP) -c $(CCFLAGS) $(OTHERFLAGS) main.cpp

toyModel.o: toyModel.cpp $(FOR_ALL)
	$(CPP) -c $(CCFLAGS) $(OTHERFLAGS) toyModel.cpp

rng.o:   rng.c $(FOR_ALL)
	$(CC) -c $(CCFLAGS) rng.c

#      Build and run executables

#      "main program" sgldTest

sgldTestExecutable: main.o toyModel.o rng.o
	$(LINKER)  -o sgldTestExecutable main.o toyModel.o rng.o

genResults: sgldTestExecutable 
	./sgldTestExecutable

sgldTest: genResults 
	python plotting.py

tarball: $( All_SOURCES)  
	tar -cvf sgldTestSuite.tar $(ALL_SOURCES) 



