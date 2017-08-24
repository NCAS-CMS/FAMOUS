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
NORMAL/drintfa1.o NORMAL/fldstat1.o NORMAL/thl2th1.o

all: $(OBJS)


NORMAL/drintfa1.o: drintfa1.f
	$(FORT) $(FORTOPTS)  -c drintfa1.f
	@if test -f drintfa1.T;then echo "	mv -f drintfa1.T NORMAL/.";	mv -f drintfa1.T NORMAL/. ; fi
	mv -f drintfa1.o NORMAL/drintfa1.o

NORMAL/fldstat1.o: fldstat1.f
	$(FORT) $(FORTOPTS)  -c fldstat1.f
	@if test -f fldstat1.T;then echo "	mv -f fldstat1.T NORMAL/.";	mv -f fldstat1.T NORMAL/. ; fi
	mv -f fldstat1.o NORMAL/fldstat1.o

NORMAL/thl2th1.o: thl2th1.f
	$(FORT) $(FORTOPTS)  -c thl2th1.f
	@if test -f thl2th1.T;then echo "	mv -f thl2th1.T NORMAL/.";	mv -f thl2th1.T NORMAL/. ; fi
	mv -f thl2th1.o NORMAL/thl2th1.o

