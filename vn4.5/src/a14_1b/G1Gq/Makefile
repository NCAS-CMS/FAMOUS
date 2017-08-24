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
NORMAL/aencor1a.o NORMAL/cemcor1a.o NORMAL/emcalc1b.o NORMAL/emdiag1b.o\
NORMAL/fldiag1a.o

all: $(OBJS)


NORMAL/aencor1a.o: aencor1a.f
	$(FORT) $(FORTOPTS)  -c aencor1a.f
	@if test -f aencor1a.T;then echo "	mv -f aencor1a.T NORMAL/.";	mv -f aencor1a.T NORMAL/. ; fi
	mv -f aencor1a.o NORMAL/aencor1a.o

NORMAL/cemcor1a.o: cemcor1a.f
	$(FORT) $(FORTOPTS)  -c cemcor1a.f
	@if test -f cemcor1a.T;then echo "	mv -f cemcor1a.T NORMAL/.";	mv -f cemcor1a.T NORMAL/. ; fi
	mv -f cemcor1a.o NORMAL/cemcor1a.o

NORMAL/emcalc1b.o: emcalc1b.f
	$(FORT) $(FORTOPTS)  -c emcalc1b.f
	@if test -f emcalc1b.T;then echo "	mv -f emcalc1b.T NORMAL/.";	mv -f emcalc1b.T NORMAL/. ; fi
	mv -f emcalc1b.o NORMAL/emcalc1b.o

NORMAL/emdiag1b.o: emdiag1b.f
	$(FORT) $(FORTOPTS)  -c emdiag1b.f
	@if test -f emdiag1b.T;then echo "	mv -f emdiag1b.T NORMAL/.";	mv -f emdiag1b.T NORMAL/. ; fi
	mv -f emdiag1b.o NORMAL/emdiag1b.o

NORMAL/fldiag1a.o: fldiag1a.f
	$(FORT) $(FORTOPTS)  -c fldiag1a.f
	@if test -f fldiag1a.T;then echo "	mv -f fldiag1a.T NORMAL/.";	mv -f fldiag1a.T NORMAL/. ; fi
	mv -f fldiag1a.o NORMAL/fldiag1a.o

