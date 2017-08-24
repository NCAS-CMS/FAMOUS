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
NORMAL/gwave3a.o NORMAL/gwjump3a.o NORMAL/gwlee3a.o NORMAL/gwsatn3a.o\
NORMAL/gwscor3a.o NORMAL/gwsurf3a.o NORMAL/gwvert3a.o NORMAL/gwvict3a.o

all: $(OBJS)


NORMAL/gwave3a.o: gwave3a.f
	$(FORT) $(FORTOPTS)  -c gwave3a.f
	@if test -f gwave3a.T;then echo "	mv -f gwave3a.T NORMAL/.";	mv -f gwave3a.T NORMAL/. ; fi
	mv -f gwave3a.o NORMAL/gwave3a.o

NORMAL/gwjump3a.o: gwjump3a.f
	$(FORT) $(FORTOPTS)  -c gwjump3a.f
	@if test -f gwjump3a.T;then echo "	mv -f gwjump3a.T NORMAL/.";	mv -f gwjump3a.T NORMAL/. ; fi
	mv -f gwjump3a.o NORMAL/gwjump3a.o

NORMAL/gwlee3a.o: gwlee3a.f
	$(FORT) $(FORTOPTS)  -c gwlee3a.f
	@if test -f gwlee3a.T;then echo "	mv -f gwlee3a.T NORMAL/.";	mv -f gwlee3a.T NORMAL/. ; fi
	mv -f gwlee3a.o NORMAL/gwlee3a.o

NORMAL/gwsatn3a.o: gwsatn3a.f
	$(FORT) $(FORTOPTS)  -c gwsatn3a.f
	@if test -f gwsatn3a.T;then echo "	mv -f gwsatn3a.T NORMAL/.";	mv -f gwsatn3a.T NORMAL/. ; fi
	mv -f gwsatn3a.o NORMAL/gwsatn3a.o

NORMAL/gwscor3a.o: gwscor3a.f
	$(FORT) $(FORTOPTS)  -c gwscor3a.f
	@if test -f gwscor3a.T;then echo "	mv -f gwscor3a.T NORMAL/.";	mv -f gwscor3a.T NORMAL/. ; fi
	mv -f gwscor3a.o NORMAL/gwscor3a.o

NORMAL/gwsurf3a.o: gwsurf3a.f
	$(FORT) $(FORTOPTS)  -c gwsurf3a.f
	@if test -f gwsurf3a.T;then echo "	mv -f gwsurf3a.T NORMAL/.";	mv -f gwsurf3a.T NORMAL/. ; fi
	mv -f gwsurf3a.o NORMAL/gwsurf3a.o

NORMAL/gwvert3a.o: gwvert3a.f
	$(FORT) $(FORTOPTS)  -c gwvert3a.f
	@if test -f gwvert3a.T;then echo "	mv -f gwvert3a.T NORMAL/.";	mv -f gwvert3a.T NORMAL/. ; fi
	mv -f gwvert3a.o NORMAL/gwvert3a.o

NORMAL/gwvict3a.o: gwvict3a.f
	$(FORT) $(FORTOPTS)  -c gwvict3a.f
	@if test -f gwvict3a.T;then echo "	mv -f gwvict3a.T NORMAL/.";	mv -f gwvict3a.T NORMAL/. ; fi
	mv -f gwvict3a.o NORMAL/gwvict3a.o

