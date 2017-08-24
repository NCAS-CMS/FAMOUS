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
NORMAL/empty.o NORMAL/inicmc.o NORMAL/inioasis.o NORMAL/inizptr.o\
NORMAL/locrd.o NORMAL/locwt.o NORMAL/oasisdiag.o NORMAL/oasisdiagi.o\
NORMAL/oasistep.o NORMAL/pdefmr.o NORMAL/pdefmw.o NORMAL/pinim.o

all: $(OBJS)


NORMAL/empty.o: empty.f
	$(FORT) $(FORTOPTS)  -c empty.f
	@if test -f empty.T;then echo "	mv -f empty.T NORMAL/.";	mv -f empty.T NORMAL/. ; fi
	mv -f empty.o NORMAL/empty.o

NORMAL/inicmc.o: inicmc.f
	$(FORT) $(FORTOPTS)  -c inicmc.f
	@if test -f inicmc.T;then echo "	mv -f inicmc.T NORMAL/.";	mv -f inicmc.T NORMAL/. ; fi
	mv -f inicmc.o NORMAL/inicmc.o

NORMAL/inioasis.o: inioasis.f
	$(FORT) $(FORTOPTS)  -c inioasis.f
	@if test -f inioasis.T;then echo "	mv -f inioasis.T NORMAL/.";	mv -f inioasis.T NORMAL/. ; fi
	mv -f inioasis.o NORMAL/inioasis.o

NORMAL/inizptr.o: inizptr.f
	$(FORT) $(FORTOPTS)  -c inizptr.f
	@if test -f inizptr.T;then echo "	mv -f inizptr.T NORMAL/.";	mv -f inizptr.T NORMAL/. ; fi
	mv -f inizptr.o NORMAL/inizptr.o

NORMAL/locrd.o: locrd.f
	$(FORT) $(FORTOPTS)  -c locrd.f
	@if test -f locrd.T;then echo "	mv -f locrd.T NORMAL/.";	mv -f locrd.T NORMAL/. ; fi
	mv -f locrd.o NORMAL/locrd.o

NORMAL/locwt.o: locwt.f
	$(FORT) $(FORTOPTS)  -c locwt.f
	@if test -f locwt.T;then echo "	mv -f locwt.T NORMAL/.";	mv -f locwt.T NORMAL/. ; fi
	mv -f locwt.o NORMAL/locwt.o

NORMAL/oasisdiag.o: oasisdiag.f
	$(FORT) $(FORTOPTS)  -c oasisdiag.f
	@if test -f oasisdiag.T;then echo "	mv -f oasisdiag.T NORMAL/.";	mv -f oasisdiag.T NORMAL/. ; fi
	mv -f oasisdiag.o NORMAL/oasisdiag.o

NORMAL/oasisdiagi.o: oasisdiagi.f
	$(FORT) $(FORTOPTS)  -c oasisdiagi.f
	@if test -f oasisdiagi.T;then echo "	mv -f oasisdiagi.T NORMAL/.";	mv -f oasisdiagi.T NORMAL/. ; fi
	mv -f oasisdiagi.o NORMAL/oasisdiagi.o

NORMAL/oasistep.o: oasistep.f
	$(FORT) $(FORTOPTS)  -c oasistep.f
	@if test -f oasistep.T;then echo "	mv -f oasistep.T NORMAL/.";	mv -f oasistep.T NORMAL/. ; fi
	mv -f oasistep.o NORMAL/oasistep.o

NORMAL/pdefmr.o: pdefmr.f
	$(FORT) $(FORTOPTS)  -c pdefmr.f
	@if test -f pdefmr.T;then echo "	mv -f pdefmr.T NORMAL/.";	mv -f pdefmr.T NORMAL/. ; fi
	mv -f pdefmr.o NORMAL/pdefmr.o

NORMAL/pdefmw.o: pdefmw.f
	$(FORT) $(FORTOPTS)  -c pdefmw.f
	@if test -f pdefmw.T;then echo "	mv -f pdefmw.T NORMAL/.";	mv -f pdefmw.T NORMAL/. ; fi
	mv -f pdefmw.o NORMAL/pdefmw.o

NORMAL/pinim.o: pinim.f
	$(FORT) $(FORTOPTS)  -c pinim.f
	@if test -f pinim.T;then echo "	mv -f pinim.T NORMAL/.";	mv -f pinim.T NORMAL/. ; fi
	mv -f pinim.o NORMAL/pinim.o

