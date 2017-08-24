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
NORMAL/initacc1.o NORMAL/initmin1.o NORMAL/initveg1.o NORMAL/phenol1a.o\
NORMAL/pparm1a.o NORMAL/sparm1a.o NORMAL/tpts1a.o NORMAL/veg1a.o\
NORMAL/veg_ic1a.o

all: $(OBJS)


NORMAL/initacc1.o: initacc1.f
	$(FORT) $(FORTOPTS)  -c initacc1.f
	@if test -f initacc1.T;then echo "	mv -f initacc1.T NORMAL/.";	mv -f initacc1.T NORMAL/. ; fi
	mv -f initacc1.o NORMAL/initacc1.o

NORMAL/initmin1.o: initmin1.f
	$(FORT) $(FORTOPTS)  -c initmin1.f
	@if test -f initmin1.T;then echo "	mv -f initmin1.T NORMAL/.";	mv -f initmin1.T NORMAL/. ; fi
	mv -f initmin1.o NORMAL/initmin1.o

NORMAL/initveg1.o: initveg1.f
	$(FORT) $(FORTOPTS)  -c initveg1.f
	@if test -f initveg1.T;then echo "	mv -f initveg1.T NORMAL/.";	mv -f initveg1.T NORMAL/. ; fi
	mv -f initveg1.o NORMAL/initveg1.o

NORMAL/phenol1a.o: phenol1a.f
	$(FORT) $(FORTOPTS)  -c phenol1a.f
	@if test -f phenol1a.T;then echo "	mv -f phenol1a.T NORMAL/.";	mv -f phenol1a.T NORMAL/. ; fi
	mv -f phenol1a.o NORMAL/phenol1a.o

NORMAL/pparm1a.o: pparm1a.f
	$(FORT) $(FORTOPTS)  -c pparm1a.f
	@if test -f pparm1a.T;then echo "	mv -f pparm1a.T NORMAL/.";	mv -f pparm1a.T NORMAL/. ; fi
	mv -f pparm1a.o NORMAL/pparm1a.o

NORMAL/sparm1a.o: sparm1a.f
	$(FORT) $(FORTOPTS)  -c sparm1a.f
	@if test -f sparm1a.T;then echo "	mv -f sparm1a.T NORMAL/.";	mv -f sparm1a.T NORMAL/. ; fi
	mv -f sparm1a.o NORMAL/sparm1a.o

NORMAL/tpts1a.o: tpts1a.f
	$(FORT) $(FORTOPTS)  -c tpts1a.f
	@if test -f tpts1a.T;then echo "	mv -f tpts1a.T NORMAL/.";	mv -f tpts1a.T NORMAL/. ; fi
	mv -f tpts1a.o NORMAL/tpts1a.o

NORMAL/veg1a.o: veg1a.f
	$(FORT) $(FORTOPTS)  -c veg1a.f
	@if test -f veg1a.T;then echo "	mv -f veg1a.T NORMAL/.";	mv -f veg1a.T NORMAL/. ; fi
	mv -f veg1a.o NORMAL/veg1a.o

NORMAL/veg_ic1a.o: veg_ic1a.f
	$(FORT) $(FORTOPTS)  -c veg_ic1a.f
	@if test -f veg_ic1a.T;then echo "	mv -f veg_ic1a.T NORMAL/.";	mv -f veg_ic1a.T NORMAL/. ; fi
	mv -f veg_ic1a.o NORMAL/veg_ic1a.o

