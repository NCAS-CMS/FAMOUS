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
NORMAL/cllvd3a.o NORMAL/ftsa1a.o NORMAL/solang1a.o NORMAL/solpos1a.o\
NORMAL/swrad3a.o

all: $(OBJS)


NORMAL/cllvd3a.o: cllvd3a.f
	$(FORT) $(FORTOPTS)  -c cllvd3a.f
	@if test -f cllvd3a.T;then echo "	mv -f cllvd3a.T NORMAL/.";	mv -f cllvd3a.T NORMAL/. ; fi
	mv -f cllvd3a.o NORMAL/cllvd3a.o

NORMAL/ftsa1a.o: ftsa1a.f
	$(FORT) $(FORTOPTS)  -c ftsa1a.f
	@if test -f ftsa1a.T;then echo "	mv -f ftsa1a.T NORMAL/.";	mv -f ftsa1a.T NORMAL/. ; fi
	mv -f ftsa1a.o NORMAL/ftsa1a.o

NORMAL/solang1a.o: solang1a.f
	$(FORT) $(FORTOPTS)  -c solang1a.f
	@if test -f solang1a.T;then echo "	mv -f solang1a.T NORMAL/.";	mv -f solang1a.T NORMAL/. ; fi
	mv -f solang1a.o NORMAL/solang1a.o

NORMAL/solpos1a.o: solpos1a.f
	$(FORT) $(FORTOPTS)  -c solpos1a.f
	@if test -f solpos1a.T;then echo "	mv -f solpos1a.T NORMAL/.";	mv -f solpos1a.T NORMAL/. ; fi
	mv -f solpos1a.o NORMAL/solpos1a.o

NORMAL/swrad3a.o: swrad3a.f
	$(FORT) $(FORTOPTS)  -c swrad3a.f
	@if test -f swrad3a.T;then echo "	mv -f swrad3a.T NORMAL/.";	mv -f swrad3a.T NORMAL/. ; fi
	mv -f swrad3a.o NORMAL/swrad3a.o

