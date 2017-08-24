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
NORMAL/frunof1a.o NORMAL/hyd_intct1a.o NORMAL/hydrol1a.o NORMAL/sfsnow1a.o\
NORMAL/sieve1a.o NORMAL/soilhy1a.o NORMAL/surfhy1a.o

all: $(OBJS)


NORMAL/frunof1a.o: frunof1a.f
	$(FORT) $(FORTOPTS)  -c frunof1a.f
	@if test -f frunof1a.T;then echo "	mv -f frunof1a.T NORMAL/.";	mv -f frunof1a.T NORMAL/. ; fi
	mv -f frunof1a.o NORMAL/frunof1a.o

NORMAL/hyd_intct1a.o: hyd_intct1a.f
	$(FORT) $(FORTOPTS)  -c hyd_intct1a.f
	@if test -f hyd_intct1a.T;then echo "	mv -f hyd_intct1a.T NORMAL/.";	mv -f hyd_intct1a.T NORMAL/. ; fi
	mv -f hyd_intct1a.o NORMAL/hyd_intct1a.o

NORMAL/hydrol1a.o: hydrol1a.f
	$(FORT) $(FORTOPTS)  -c hydrol1a.f
	@if test -f hydrol1a.T;then echo "	mv -f hydrol1a.T NORMAL/.";	mv -f hydrol1a.T NORMAL/. ; fi
	mv -f hydrol1a.o NORMAL/hydrol1a.o

NORMAL/sfsnow1a.o: sfsnow1a.f
	$(FORT) $(FORTOPTS)  -c sfsnow1a.f
	@if test -f sfsnow1a.T;then echo "	mv -f sfsnow1a.T NORMAL/.";	mv -f sfsnow1a.T NORMAL/. ; fi
	mv -f sfsnow1a.o NORMAL/sfsnow1a.o

NORMAL/sieve1a.o: sieve1a.f
	$(FORT) $(FORTOPTS)  -c sieve1a.f
	@if test -f sieve1a.T;then echo "	mv -f sieve1a.T NORMAL/.";	mv -f sieve1a.T NORMAL/. ; fi
	mv -f sieve1a.o NORMAL/sieve1a.o

NORMAL/soilhy1a.o: soilhy1a.f
	$(FORT) $(FORTOPTS)  -c soilhy1a.f
	@if test -f soilhy1a.T;then echo "	mv -f soilhy1a.T NORMAL/.";	mv -f soilhy1a.T NORMAL/. ; fi
	mv -f soilhy1a.o NORMAL/soilhy1a.o

NORMAL/surfhy1a.o: surfhy1a.f
	$(FORT) $(FORTOPTS)  -c surfhy1a.f
	@if test -f surfhy1a.T;then echo "	mv -f surfhy1a.T NORMAL/.";	mv -f surfhy1a.T NORMAL/. ; fi
	mv -f surfhy1a.o NORMAL/surfhy1a.o

