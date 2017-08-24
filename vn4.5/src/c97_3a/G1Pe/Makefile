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
NORMAL/timefn2a.o NORMAL/timer3a.o

all: $(OBJS)


NORMAL/timefn2a.o: timefn2a.f
	$(FORT) $(FORTOPTS)  -c timefn2a.f
	@if test -f timefn2a.T;then echo "	mv -f timefn2a.T NORMAL/.";	mv -f timefn2a.T NORMAL/. ; fi
	mv -f timefn2a.o NORMAL/timefn2a.o

NORMAL/timer3a.o: timer3a.f
	$(FORT) $(FORTOPTS)  -c timer3a.f
	@if test -f timer3a.T;then echo "	mv -f timer3a.T NORMAL/.";	mv -f timer3a.T NORMAL/. ; fi
	mv -f timer3a.o NORMAL/timer3a.o

