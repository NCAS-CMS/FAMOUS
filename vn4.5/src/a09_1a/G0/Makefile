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
NORMAL/cld2rh1a.o NORMAL/cldbas1a.o NORMAL/cldcvr1a.o NORMAL/gcld_1a.o\
NORMAL/lscld1a.o

all: $(OBJS)


NORMAL/cld2rh1a.o: cld2rh1a.f
	$(FORT) $(FORTOPTS)  -c cld2rh1a.f
	@if test -f cld2rh1a.T;then echo "	mv -f cld2rh1a.T NORMAL/.";	mv -f cld2rh1a.T NORMAL/. ; fi
	mv -f cld2rh1a.o NORMAL/cld2rh1a.o

NORMAL/cldbas1a.o: cldbas1a.f
	$(FORT) $(FORTOPTS)  -c cldbas1a.f
	@if test -f cldbas1a.T;then echo "	mv -f cldbas1a.T NORMAL/.";	mv -f cldbas1a.T NORMAL/. ; fi
	mv -f cldbas1a.o NORMAL/cldbas1a.o

NORMAL/cldcvr1a.o: cldcvr1a.f
	$(FORT) $(FORTOPTS)  -c cldcvr1a.f
	@if test -f cldcvr1a.T;then echo "	mv -f cldcvr1a.T NORMAL/.";	mv -f cldcvr1a.T NORMAL/. ; fi
	mv -f cldcvr1a.o NORMAL/cldcvr1a.o

NORMAL/gcld_1a.o: gcld_1a.f
	$(FORT) $(FORTOPTS)  -c gcld_1a.f
	@if test -f gcld_1a.T;then echo "	mv -f gcld_1a.T NORMAL/.";	mv -f gcld_1a.T NORMAL/. ; fi
	mv -f gcld_1a.o NORMAL/gcld_1a.o

NORMAL/lscld1a.o: lscld1a.f
	$(FORT) $(FORTOPTS)  -c lscld1a.f
	@if test -f lscld1a.T;then echo "	mv -f lscld1a.T NORMAL/.";	mv -f lscld1a.T NORMAL/. ; fi
	mv -f lscld1a.o NORMAL/lscld1a.o

