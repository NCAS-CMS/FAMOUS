# ----------------------------- COPYRIGHT ---------------------------- #
#           (c) BRITISH CROWN COPYRIGHT 2000, THE MET.OFFICE           #
#                                                                      #
# Use, duplication or disclosure of this code is subject to the        #
# restrictions as set forth in the contract. If no contract has been   #
# raised with this copy of the code, use, duplication or disclosure    #
# of it is strictly prohibited. Permission to do so must be obtained   #
# in writing from the Head of Numerical Modelling at The Met. Office.  #
# ----------------------------- COPYRIGHT ---------------------------- #
# Unset ENV to make sure that .kshrc (or anything else) is not run
ENV =

# Number of things to do in parallel
NPROC = 1

# Fortran compiler names
FORT = mpif90

# Fortran compiler options
FORTOPTS = -i8 -r8 -fp-model source -O2

# C Compiler name
CC = gcc

# C Compiler options
CCOPTS = 

# Object files
OBJS=\
NORMAL/aencor1a.o NORMAL/cemcor1a.o NORMAL/diaeng1a.o NORMAL/emdiag1a.o\
NORMAL/engsum1a.o NORMAL/fldiag1a.o NORMAL/massum1a.o

all: $(OBJS)


NORMAL/aencor1a.o: aencor1a.f
	$(FORT) $(FORTOPTS)  -c aencor1a.f
	@if test -f aencor1a.T;then echo "	mv -f aencor1a.T NORMAL/.";	mv -f aencor1a.T NORMAL/. ; fi
	mv -f aencor1a.o NORMAL/aencor1a.o

NORMAL/cemcor1a.o: cemcor1a.f
	$(FORT) $(FORTOPTS)  -c cemcor1a.f
	@if test -f cemcor1a.T;then echo "	mv -f cemcor1a.T NORMAL/.";	mv -f cemcor1a.T NORMAL/. ; fi
	mv -f cemcor1a.o NORMAL/cemcor1a.o

NORMAL/diaeng1a.o: diaeng1a.f
	$(FORT) $(FORTOPTS)  -c diaeng1a.f
	@if test -f diaeng1a.T;then echo "	mv -f diaeng1a.T NORMAL/.";	mv -f diaeng1a.T NORMAL/. ; fi
	mv -f diaeng1a.o NORMAL/diaeng1a.o

NORMAL/emdiag1a.o: emdiag1a.f
	$(FORT) $(FORTOPTS)  -c emdiag1a.f
	@if test -f emdiag1a.T;then echo "	mv -f emdiag1a.T NORMAL/.";	mv -f emdiag1a.T NORMAL/. ; fi
	mv -f emdiag1a.o NORMAL/emdiag1a.o

NORMAL/engsum1a.o: engsum1a.f
	$(FORT) $(FORTOPTS)  -c engsum1a.f
	@if test -f engsum1a.T;then echo "	mv -f engsum1a.T NORMAL/.";	mv -f engsum1a.T NORMAL/. ; fi
	mv -f engsum1a.o NORMAL/engsum1a.o

NORMAL/fldiag1a.o: fldiag1a.f
	$(FORT) $(FORTOPTS)  -c fldiag1a.f
	@if test -f fldiag1a.T;then echo "	mv -f fldiag1a.T NORMAL/.";	mv -f fldiag1a.T NORMAL/. ; fi
	mv -f fldiag1a.o NORMAL/fldiag1a.o

NORMAL/massum1a.o: massum1a.f
	$(FORT) $(FORTOPTS)  -c massum1a.f
	@if test -f massum1a.T;then echo "	mv -f massum1a.T NORMAL/.";	mv -f massum1a.T NORMAL/. ; fi
	mv -f massum1a.o NORMAL/massum1a.o

