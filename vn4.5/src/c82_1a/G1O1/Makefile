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
NORMAL/inanca1a.o NORMAL/inanco1a.o NORMAL/rpanca1a.o NORMAL/rpanco1a.o

all: $(OBJS)


NORMAL/inanca1a.o: inanca1a.f
	$(FORT) $(FORTOPTS)  -c inanca1a.f
	@if test -f inanca1a.T;then echo "	mv -f inanca1a.T NORMAL/.";	mv -f inanca1a.T NORMAL/. ; fi
	mv -f inanca1a.o NORMAL/inanca1a.o

NORMAL/inanco1a.o: inanco1a.f
	$(FORT) $(FORTOPTS)  -c inanco1a.f
	@if test -f inanco1a.T;then echo "	mv -f inanco1a.T NORMAL/.";	mv -f inanco1a.T NORMAL/. ; fi
	mv -f inanco1a.o NORMAL/inanco1a.o

NORMAL/rpanca1a.o: rpanca1a.f
	$(FORT) $(FORTOPTS)  -c rpanca1a.f
	@if test -f rpanca1a.T;then echo "	mv -f rpanca1a.T NORMAL/.";	mv -f rpanca1a.T NORMAL/. ; fi
	mv -f rpanca1a.o NORMAL/rpanca1a.o

NORMAL/rpanco1a.o: rpanco1a.f
	$(FORT) $(FORTOPTS)  -c rpanco1a.f
	@if test -f rpanco1a.T;then echo "	mv -f rpanco1a.T NORMAL/.";	mv -f rpanco1a.T NORMAL/. ; fi
	mv -f rpanco1a.o NORMAL/rpanco1a.o

