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
NORMAL/propags.o NORMAL/propdot.o

all: $(OBJS)


NORMAL/propags.o: propags.f
	$(FORT) $(FORTOPTS)  -c propags.f
	@if test -f propags.T;then echo "	mv -f propags.T NORMAL/.";	mv -f propags.T NORMAL/. ; fi
	mv -f propags.o NORMAL/propags.o

NORMAL/propdot.o: propdot.f
	$(FORT) $(FORTOPTS)  -c propdot.f
	@if test -f propdot.T;then echo "	mv -f propdot.T NORMAL/.";	mv -f propdot.T NORMAL/. ; fi
	mv -f propdot.o NORMAL/propdot.o

