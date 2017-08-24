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
NORMAL/nlweigt.o NORMAL/snonlin.o

all: $(OBJS)


NORMAL/nlweigt.o: nlweigt.f
	$(FORT) $(FORTOPTS)  -c nlweigt.f
	@if test -f nlweigt.T;then echo "	mv -f nlweigt.T NORMAL/.";	mv -f nlweigt.T NORMAL/. ; fi
	mv -f nlweigt.o NORMAL/nlweigt.o

NORMAL/snonlin.o: snonlin.f
	$(FORT) $(FORTOPTS)  -c snonlin.f
	@if test -f snonlin.T;then echo "	mv -f snonlin.T NORMAL/.";	mv -f snonlin.T NORMAL/. ; fi
	mv -f snonlin.o NORMAL/snonlin.o

