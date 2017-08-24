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
NORMAL/fromla1a.o NORMAL/p21bit1a.o NORMAL/toland1a.o

all: $(OBJS)


NORMAL/fromla1a.o: fromla1a.f
	$(FORT) $(FORTOPTS)  -c fromla1a.f
	@if test -f fromla1a.T;then echo "	mv -f fromla1a.T NORMAL/.";	mv -f fromla1a.T NORMAL/. ; fi
	mv -f fromla1a.o NORMAL/fromla1a.o

NORMAL/p21bit1a.o: p21bit1a.f
	$(FORT) $(FORTOPTS)  -c p21bit1a.f
	@if test -f p21bit1a.T;then echo "	mv -f p21bit1a.T NORMAL/.";	mv -f p21bit1a.T NORMAL/. ; fi
	mv -f p21bit1a.o NORMAL/p21bit1a.o

NORMAL/toland1a.o: toland1a.f
	$(FORT) $(FORTOPTS)  -c toland1a.f
	@if test -f toland1a.T;then echo "	mv -f toland1a.T NORMAL/.";	mv -f toland1a.T NORMAL/. ; fi
	mv -f toland1a.o NORMAL/toland1a.o

