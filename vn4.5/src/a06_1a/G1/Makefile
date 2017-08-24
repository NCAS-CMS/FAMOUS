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
NORMAL/gwave1a.o NORMAL/gwrich1a.o NORMAL/gwsurf1a.o NORMAL/gwvict1a.o

all: $(OBJS)


NORMAL/gwave1a.o: gwave1a.f
	$(FORT) $(FORTOPTS)  -c gwave1a.f
	@if test -f gwave1a.T;then echo "	mv -f gwave1a.T NORMAL/.";	mv -f gwave1a.T NORMAL/. ; fi
	mv -f gwave1a.o NORMAL/gwave1a.o

NORMAL/gwrich1a.o: gwrich1a.f
	$(FORT) $(FORTOPTS)  -c gwrich1a.f
	@if test -f gwrich1a.T;then echo "	mv -f gwrich1a.T NORMAL/.";	mv -f gwrich1a.T NORMAL/. ; fi
	mv -f gwrich1a.o NORMAL/gwrich1a.o

NORMAL/gwsurf1a.o: gwsurf1a.f
	$(FORT) $(FORTOPTS)  -c gwsurf1a.f
	@if test -f gwsurf1a.T;then echo "	mv -f gwsurf1a.T NORMAL/.";	mv -f gwsurf1a.T NORMAL/. ; fi
	mv -f gwsurf1a.o NORMAL/gwsurf1a.o

NORMAL/gwvict1a.o: gwvict1a.f
	$(FORT) $(FORTOPTS)  -c gwvict1a.f
	@if test -f gwvict1a.T;then echo "	mv -f gwvict1a.T NORMAL/.";	mv -f gwvict1a.T NORMAL/. ; fi
	mv -f gwvict1a.o NORMAL/gwvict1a.o

