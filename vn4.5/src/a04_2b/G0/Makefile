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
NORMAL/glspp_2b.o NORMAL/lspeva2a.o NORMAL/lspfcw1a.o NORMAL/lspfor2b.o\
NORMAL/lspfrm1a.o NORMAL/lsppn2b.o NORMAL/lspsc1a.o

all: $(OBJS)


NORMAL/glspp_2b.o: glspp_2b.f
	$(FORT) $(FORTOPTS)  -c glspp_2b.f
	@if test -f glspp_2b.T;then echo "	mv -f glspp_2b.T NORMAL/.";	mv -f glspp_2b.T NORMAL/. ; fi
	mv -f glspp_2b.o NORMAL/glspp_2b.o

NORMAL/lspeva2a.o: lspeva2a.f
	$(FORT) $(FORTOPTS)  -c lspeva2a.f
	@if test -f lspeva2a.T;then echo "	mv -f lspeva2a.T NORMAL/.";	mv -f lspeva2a.T NORMAL/. ; fi
	mv -f lspeva2a.o NORMAL/lspeva2a.o

NORMAL/lspfcw1a.o: lspfcw1a.f
	$(FORT) $(FORTOPTS)  -c lspfcw1a.f
	@if test -f lspfcw1a.T;then echo "	mv -f lspfcw1a.T NORMAL/.";	mv -f lspfcw1a.T NORMAL/. ; fi
	mv -f lspfcw1a.o NORMAL/lspfcw1a.o

NORMAL/lspfor2b.o: lspfor2b.f
	$(FORT) $(FORTOPTS)  -c lspfor2b.f
	@if test -f lspfor2b.T;then echo "	mv -f lspfor2b.T NORMAL/.";	mv -f lspfor2b.T NORMAL/. ; fi
	mv -f lspfor2b.o NORMAL/lspfor2b.o

NORMAL/lspfrm1a.o: lspfrm1a.f
	$(FORT) $(FORTOPTS)  -c lspfrm1a.f
	@if test -f lspfrm1a.T;then echo "	mv -f lspfrm1a.T NORMAL/.";	mv -f lspfrm1a.T NORMAL/. ; fi
	mv -f lspfrm1a.o NORMAL/lspfrm1a.o

NORMAL/lsppn2b.o: lsppn2b.f
	$(FORT) $(FORTOPTS)  -c lsppn2b.f
	@if test -f lsppn2b.T;then echo "	mv -f lsppn2b.T NORMAL/.";	mv -f lsppn2b.T NORMAL/. ; fi
	mv -f lsppn2b.o NORMAL/lsppn2b.o

NORMAL/lspsc1a.o: lspsc1a.f
	$(FORT) $(FORTOPTS)  -c lspsc1a.f
	@if test -f lspsc1a.T;then echo "	mv -f lspsc1a.T NORMAL/.";	mv -f lspsc1a.T NORMAL/. ; fi
	mv -f lspsc1a.o NORMAL/lspsc1a.o

