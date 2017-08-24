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
NORMAL/lwcld1a.o NORMAL/lwdcsf1a.o NORMAL/lwmast1a.o NORMAL/lwplan1b.o\
NORMAL/lwptsc1b.o NORMAL/lwrad1a.o NORMAL/lwtran1b.o

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

NORMAL/lwplan1b.o: lwplan1b.f
	$(FORT) $(FORTOPTS)  -c lwplan1b.f
	@if test -f lwplan1b.T;then echo "	mv -f lwplan1b.T NORMAL/.";	mv -f lwplan1b.T NORMAL/. ; fi
	mv -f lwplan1b.o NORMAL/lwplan1b.o

NORMAL/lwptsc1b.o: lwptsc1b.f
	$(FORT) $(FORTOPTS)  -c lwptsc1b.f
	@if test -f lwptsc1b.T;then echo "	mv -f lwptsc1b.T NORMAL/.";	mv -f lwptsc1b.T NORMAL/. ; fi
	mv -f lwptsc1b.o NORMAL/lwptsc1b.o

NORMAL/lwrad1a.o: lwrad1a.f
	$(FORT) $(FORTOPTS)  -c lwrad1a.f
	@if test -f lwrad1a.T;then echo "	mv -f lwrad1a.T NORMAL/.";	mv -f lwrad1a.T NORMAL/. ; fi
	mv -f lwrad1a.o NORMAL/lwrad1a.o

NORMAL/lwtran1b.o: lwtran1b.f
	$(FORT) $(FORTOPTS)  -c lwtran1b.f
	@if test -f lwtran1b.T;then echo "	mv -f lwtran1b.T NORMAL/.";	mv -f lwtran1b.T NORMAL/. ; fi
	mv -f lwtran1b.o NORMAL/lwtran1b.o

