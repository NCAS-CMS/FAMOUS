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
NORMAL/calcp21a.o NORMAL/calcpv1a.o NORMAL/cat1a.o NORMAL/dthdp1a.o\
NORMAL/dyndia1a.o NORMAL/omegdi1a.o NORMAL/pvpin1a.o NORMAL/pvthin1a.o\
NORMAL/testdi1a.o NORMAL/thetpv1a.o NORMAL/thpvin1a.o NORMAL/vorti11a.o\
NORMAL/vorti21a.o NORMAL/vorti31a.o NORMAL/vorti41a.o NORMAL/vorti51a.o\
NORMAL/windmx1a.o

all: $(OBJS)


NORMAL/calcp21a.o: calcp21a.f
	$(FORT) $(FORTOPTS)  -c calcp21a.f
	@if test -f calcp21a.T;then echo "	mv -f calcp21a.T NORMAL/.";	mv -f calcp21a.T NORMAL/. ; fi
	mv -f calcp21a.o NORMAL/calcp21a.o

NORMAL/calcpv1a.o: calcpv1a.f
	$(FORT) $(FORTOPTS)  -c calcpv1a.f
	@if test -f calcpv1a.T;then echo "	mv -f calcpv1a.T NORMAL/.";	mv -f calcpv1a.T NORMAL/. ; fi
	mv -f calcpv1a.o NORMAL/calcpv1a.o

NORMAL/cat1a.o: cat1a.f
	$(FORT) $(FORTOPTS)  -c cat1a.f
	@if test -f cat1a.T;then echo "	mv -f cat1a.T NORMAL/.";	mv -f cat1a.T NORMAL/. ; fi
	mv -f cat1a.o NORMAL/cat1a.o

NORMAL/dthdp1a.o: dthdp1a.f
	$(FORT) $(FORTOPTS)  -c dthdp1a.f
	@if test -f dthdp1a.T;then echo "	mv -f dthdp1a.T NORMAL/.";	mv -f dthdp1a.T NORMAL/. ; fi
	mv -f dthdp1a.o NORMAL/dthdp1a.o

NORMAL/dyndia1a.o: dyndia1a.f
	$(FORT) $(FORTOPTS)  -c dyndia1a.f
	@if test -f dyndia1a.T;then echo "	mv -f dyndia1a.T NORMAL/.";	mv -f dyndia1a.T NORMAL/. ; fi
	mv -f dyndia1a.o NORMAL/dyndia1a.o

NORMAL/omegdi1a.o: omegdi1a.f
	$(FORT) $(FORTOPTS)  -c omegdi1a.f
	@if test -f omegdi1a.T;then echo "	mv -f omegdi1a.T NORMAL/.";	mv -f omegdi1a.T NORMAL/. ; fi
	mv -f omegdi1a.o NORMAL/omegdi1a.o

NORMAL/pvpin1a.o: pvpin1a.f
	$(FORT) $(FORTOPTS)  -c pvpin1a.f
	@if test -f pvpin1a.T;then echo "	mv -f pvpin1a.T NORMAL/.";	mv -f pvpin1a.T NORMAL/. ; fi
	mv -f pvpin1a.o NORMAL/pvpin1a.o

NORMAL/pvthin1a.o: pvthin1a.f
	$(FORT) $(FORTOPTS)  -c pvthin1a.f
	@if test -f pvthin1a.T;then echo "	mv -f pvthin1a.T NORMAL/.";	mv -f pvthin1a.T NORMAL/. ; fi
	mv -f pvthin1a.o NORMAL/pvthin1a.o

NORMAL/testdi1a.o: testdi1a.f
	$(FORT) $(FORTOPTS)  -c testdi1a.f
	@if test -f testdi1a.T;then echo "	mv -f testdi1a.T NORMAL/.";	mv -f testdi1a.T NORMAL/. ; fi
	mv -f testdi1a.o NORMAL/testdi1a.o

NORMAL/thetpv1a.o: thetpv1a.f
	$(FORT) $(FORTOPTS)  -c thetpv1a.f
	@if test -f thetpv1a.T;then echo "	mv -f thetpv1a.T NORMAL/.";	mv -f thetpv1a.T NORMAL/. ; fi
	mv -f thetpv1a.o NORMAL/thetpv1a.o

NORMAL/thpvin1a.o: thpvin1a.f
	$(FORT) $(FORTOPTS)  -c thpvin1a.f
	@if test -f thpvin1a.T;then echo "	mv -f thpvin1a.T NORMAL/.";	mv -f thpvin1a.T NORMAL/. ; fi
	mv -f thpvin1a.o NORMAL/thpvin1a.o

NORMAL/vorti11a.o: vorti11a.f
	$(FORT) $(FORTOPTS)  -c vorti11a.f
	@if test -f vorti11a.T;then echo "	mv -f vorti11a.T NORMAL/.";	mv -f vorti11a.T NORMAL/. ; fi
	mv -f vorti11a.o NORMAL/vorti11a.o

NORMAL/vorti21a.o: vorti21a.f
	$(FORT) $(FORTOPTS)  -c vorti21a.f
	@if test -f vorti21a.T;then echo "	mv -f vorti21a.T NORMAL/.";	mv -f vorti21a.T NORMAL/. ; fi
	mv -f vorti21a.o NORMAL/vorti21a.o

NORMAL/vorti31a.o: vorti31a.f
	$(FORT) $(FORTOPTS)  -c vorti31a.f
	@if test -f vorti31a.T;then echo "	mv -f vorti31a.T NORMAL/.";	mv -f vorti31a.T NORMAL/. ; fi
	mv -f vorti31a.o NORMAL/vorti31a.o

NORMAL/vorti41a.o: vorti41a.f
	$(FORT) $(FORTOPTS)  -c vorti41a.f
	@if test -f vorti41a.T;then echo "	mv -f vorti41a.T NORMAL/.";	mv -f vorti41a.T NORMAL/. ; fi
	mv -f vorti41a.o NORMAL/vorti41a.o

NORMAL/vorti51a.o: vorti51a.f
	$(FORT) $(FORTOPTS)  -c vorti51a.f
	@if test -f vorti51a.T;then echo "	mv -f vorti51a.T NORMAL/.";	mv -f vorti51a.T NORMAL/. ; fi
	mv -f vorti51a.o NORMAL/vorti51a.o

NORMAL/windmx1a.o: windmx1a.f
	$(FORT) $(FORTOPTS)  -c windmx1a.f
	@if test -f windmx1a.T;then echo "	mv -f windmx1a.T NORMAL/.";	mv -f windmx1a.T NORMAL/. ; fi
	mv -f windmx1a.o NORMAL/windmx1a.o

