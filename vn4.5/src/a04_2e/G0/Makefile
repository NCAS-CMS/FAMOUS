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
NORMAL/glspp_2d.o NORMAL/lspeva2e.o NORMAL/lspfcw2d.o NORMAL/lspfor2d.o\
NORMAL/lspfrm1a.o NORMAL/lsppn2d.o NORMAL/lspsc1a.o NORMAL/rainou1a.o\
NORMAL/slspsc1a.o

all: $(OBJS)


NORMAL/glspp_2d.o: glspp_2d.f
	$(FORT) $(FORTOPTS)  -c glspp_2d.f
	@if test -f glspp_2d.T;then echo "	mv -f glspp_2d.T NORMAL/.";	mv -f glspp_2d.T NORMAL/. ; fi
	mv -f glspp_2d.o NORMAL/glspp_2d.o

NORMAL/lspeva2e.o: lspeva2e.f
	$(FORT) $(FORTOPTS)  -c lspeva2e.f
	@if test -f lspeva2e.T;then echo "	mv -f lspeva2e.T NORMAL/.";	mv -f lspeva2e.T NORMAL/. ; fi
	mv -f lspeva2e.o NORMAL/lspeva2e.o

NORMAL/lspfcw2d.o: lspfcw2d.f
	$(FORT) $(FORTOPTS)  -c lspfcw2d.f
	@if test -f lspfcw2d.T;then echo "	mv -f lspfcw2d.T NORMAL/.";	mv -f lspfcw2d.T NORMAL/. ; fi
	mv -f lspfcw2d.o NORMAL/lspfcw2d.o

NORMAL/lspfor2d.o: lspfor2d.f
	$(FORT) $(FORTOPTS)  -c lspfor2d.f
	@if test -f lspfor2d.T;then echo "	mv -f lspfor2d.T NORMAL/.";	mv -f lspfor2d.T NORMAL/. ; fi
	mv -f lspfor2d.o NORMAL/lspfor2d.o

NORMAL/lspfrm1a.o: lspfrm1a.f
	$(FORT) $(FORTOPTS)  -c lspfrm1a.f
	@if test -f lspfrm1a.T;then echo "	mv -f lspfrm1a.T NORMAL/.";	mv -f lspfrm1a.T NORMAL/. ; fi
	mv -f lspfrm1a.o NORMAL/lspfrm1a.o

NORMAL/lsppn2d.o: lsppn2d.f
	$(FORT) $(FORTOPTS)  -c lsppn2d.f
	@if test -f lsppn2d.T;then echo "	mv -f lsppn2d.T NORMAL/.";	mv -f lsppn2d.T NORMAL/. ; fi
	mv -f lsppn2d.o NORMAL/lsppn2d.o

NORMAL/lspsc1a.o: lspsc1a.f
	$(FORT) $(FORTOPTS)  -c lspsc1a.f
	@if test -f lspsc1a.T;then echo "	mv -f lspsc1a.T NORMAL/.";	mv -f lspsc1a.T NORMAL/. ; fi
	mv -f lspsc1a.o NORMAL/lspsc1a.o

NORMAL/rainou1a.o: rainou1a.f
	$(FORT) $(FORTOPTS)  -c rainou1a.f
	@if test -f rainou1a.T;then echo "	mv -f rainou1a.T NORMAL/.";	mv -f rainou1a.T NORMAL/. ; fi
	mv -f rainou1a.o NORMAL/rainou1a.o

NORMAL/slspsc1a.o: slspsc1a.f
	$(FORT) $(FORTOPTS)  -c slspsc1a.f
	@if test -f slspsc1a.T;then echo "	mv -f slspsc1a.T NORMAL/.";	mv -f slspsc1a.T NORMAL/. ; fi
	mv -f slspsc1a.o NORMAL/slspsc1a.o

