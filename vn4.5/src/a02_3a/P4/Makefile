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
NORMAL/lwrad3a.o

all: $(OBJS)


NORMAL/lwrad3a.o: lwrad3a.f
	$(FORT) $(FORTOPTS)  -c lwrad3a.f
	@if test -f lwrad3a.T;then echo "	mv -f lwrad3a.T NORMAL/.";	mv -f lwrad3a.T NORMAL/. ; fi
	mv -f lwrad3a.o NORMAL/lwrad3a.o

