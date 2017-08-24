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
NORMAL/settra1a.o NORMAL/tracad1a.o NORMAL/travad1a.o NORMAL/trbdry1a.o\
NORMAL/trsrce1a.o

all: $(OBJS)


NORMAL/settra1a.o: settra1a.f
	$(FORT) $(FORTOPTS)  -c settra1a.f
	@if test -f settra1a.T;then echo "	mv -f settra1a.T NORMAL/.";	mv -f settra1a.T NORMAL/. ; fi
	mv -f settra1a.o NORMAL/settra1a.o

NORMAL/tracad1a.o: tracad1a.f
	$(FORT) $(FORTOPTS)  -c tracad1a.f
	@if test -f tracad1a.T;then echo "	mv -f tracad1a.T NORMAL/.";	mv -f tracad1a.T NORMAL/. ; fi
	mv -f tracad1a.o NORMAL/tracad1a.o

NORMAL/travad1a.o: travad1a.f
	$(FORT) $(FORTOPTS)  -c travad1a.f
	@if test -f travad1a.T;then echo "	mv -f travad1a.T NORMAL/.";	mv -f travad1a.T NORMAL/. ; fi
	mv -f travad1a.o NORMAL/travad1a.o

NORMAL/trbdry1a.o: trbdry1a.f
	$(FORT) $(FORTOPTS)  -c trbdry1a.f
	@if test -f trbdry1a.T;then echo "	mv -f trbdry1a.T NORMAL/.";	mv -f trbdry1a.T NORMAL/. ; fi
	mv -f trbdry1a.o NORMAL/trbdry1a.o

NORMAL/trsrce1a.o: trsrce1a.f
	$(FORT) $(FORTOPTS)  -c trsrce1a.f
	@if test -f trsrce1a.T;then echo "	mv -f trsrce1a.T NORMAL/.";	mv -f trsrce1a.T NORMAL/. ; fi
	mv -f trsrce1a.o NORMAL/trsrce1a.o

