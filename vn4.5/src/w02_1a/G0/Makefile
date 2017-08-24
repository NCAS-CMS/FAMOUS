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
NORMAL/airsea.o NORMAL/sinput.o NORMAL/stress.o NORMAL/stresso.o\
NORMAL/tauhf.o

all: $(OBJS)


NORMAL/airsea.o: airsea.f
	$(FORT) $(FORTOPTS)  -c airsea.f
	@if test -f airsea.T;then echo "	mv -f airsea.T NORMAL/.";	mv -f airsea.T NORMAL/. ; fi
	mv -f airsea.o NORMAL/airsea.o

NORMAL/sinput.o: sinput.f
	$(FORT) $(FORTOPTS)  -c sinput.f
	@if test -f sinput.T;then echo "	mv -f sinput.T NORMAL/.";	mv -f sinput.T NORMAL/. ; fi
	mv -f sinput.o NORMAL/sinput.o

NORMAL/stress.o: stress.f
	$(FORT) $(FORTOPTS)  -c stress.f
	@if test -f stress.T;then echo "	mv -f stress.T NORMAL/.";	mv -f stress.T NORMAL/. ; fi
	mv -f stress.o NORMAL/stress.o

NORMAL/stresso.o: stresso.f
	$(FORT) $(FORTOPTS)  -c stresso.f
	@if test -f stresso.T;then echo "	mv -f stresso.T NORMAL/.";	mv -f stresso.T NORMAL/. ; fi
	mv -f stresso.o NORMAL/stresso.o

NORMAL/tauhf.o: tauhf.f
	$(FORT) $(FORTOPTS)  -c tauhf.f
	@if test -f tauhf.T;then echo "	mv -f tauhf.T NORMAL/.";	mv -f tauhf.T NORMAL/. ; fi
	mv -f tauhf.o NORMAL/tauhf.o

