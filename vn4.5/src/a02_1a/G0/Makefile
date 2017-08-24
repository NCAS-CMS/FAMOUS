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
NORMAL/lwcld1a.o NORMAL/lwdcsf1a.o NORMAL/lwmast1a.o NORMAL/lwplan1a.o\
NORMAL/lwptsc1a.o NORMAL/lwrad1a.o NORMAL/lwtran1a.o

all: $(OBJS)


NORMAL/lwcld1a.o: lwcld1a.f
	$(FORT) $(FORTOPTS)  -c lwcld1a.f
	@if test -f lwcld1a.T;then echo "	mv -f lwcld1a.T NORMAL/.";	mv -f lwcld1a.T NORMAL/. ; fi
	mv -f lwcld1a.o NORMAL/lwcld1a.o

NORMAL/lwdcsf1a.o: lwdcsf1a.f
	$(FORT) $(FORTOPTS)  -c lwdcsf1a.f
	@if test -f lwdcsf1a.T;then echo "	mv -f lwdcsf1a.T NORMAL/.";	mv -f lwdcsf1a.T NORMAL/. ; fi
	mv -f lwdcsf1a.o NORMAL/lwdcsf1a.o

NORMAL/lwmast1a.o: lwmast1a.f
	$(FORT) $(FORTOPTS)  -c lwmast1a.f
	@if test -f lwmast1a.T;then echo "	mv -f lwmast1a.T NORMAL/.";	mv -f lwmast1a.T NORMAL/. ; fi
	mv -f lwmast1a.o NORMAL/lwmast1a.o

NORMAL/lwplan1a.o: lwplan1a.f
	$(FORT) $(FORTOPTS)  -c lwplan1a.f
	@if test -f lwplan1a.T;then echo "	mv -f lwplan1a.T NORMAL/.";	mv -f lwplan1a.T NORMAL/. ; fi
	mv -f lwplan1a.o NORMAL/lwplan1a.o

NORMAL/lwptsc1a.o: lwptsc1a.f
	$(FORT) $(FORTOPTS)  -c lwptsc1a.f
	@if test -f lwptsc1a.T;then echo "	mv -f lwptsc1a.T NORMAL/.";	mv -f lwptsc1a.T NORMAL/. ; fi
	mv -f lwptsc1a.o NORMAL/lwptsc1a.o

NORMAL/lwrad1a.o: lwrad1a.f
	$(FORT) $(FORTOPTS)  -c lwrad1a.f
	@if test -f lwrad1a.T;then echo "	mv -f lwrad1a.T NORMAL/.";	mv -f lwrad1a.T NORMAL/. ; fi
	mv -f lwrad1a.o NORMAL/lwrad1a.o

NORMAL/lwtran1a.o: lwtran1a.f
	$(FORT) $(FORTOPTS)  -c lwtran1a.f
	@if test -f lwtran1a.T;then echo "	mv -f lwtran1a.T NORMAL/.";	mv -f lwtran1a.T NORMAL/. ; fi
	mv -f lwtran1a.o NORMAL/lwtran1a.o

