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
NORMAL/adjctl1a.o NORMAL/angmom1a.o NORMAL/dia10a1a.o NORMAL/dia10b1a.o\
NORMAL/dia10q1a.o NORMAL/filtfl1a.o NORMAL/pthadj1a.o NORMAL/setfil1b.o\
NORMAL/thetlq1a.o NORMAL/uvadj1a.o NORMAL/vertve1b.o

all: $(OBJS)


NORMAL/adjctl1a.o: adjctl1a.f
	$(FORT) $(FORTOPTS)  -c adjctl1a.f
	@if test -f adjctl1a.T;then echo "	mv -f adjctl1a.T NORMAL/.";	mv -f adjctl1a.T NORMAL/. ; fi
	mv -f adjctl1a.o NORMAL/adjctl1a.o

NORMAL/angmom1a.o: angmom1a.f
	$(FORT) $(FORTOPTS)  -c angmom1a.f
	@if test -f angmom1a.T;then echo "	mv -f angmom1a.T NORMAL/.";	mv -f angmom1a.T NORMAL/. ; fi
	mv -f angmom1a.o NORMAL/angmom1a.o

NORMAL/dia10a1a.o: dia10a1a.f
	$(FORT) $(FORTOPTS)  -c dia10a1a.f
	@if test -f dia10a1a.T;then echo "	mv -f dia10a1a.T NORMAL/.";	mv -f dia10a1a.T NORMAL/. ; fi
	mv -f dia10a1a.o NORMAL/dia10a1a.o

NORMAL/dia10b1a.o: dia10b1a.f
	$(FORT) $(FORTOPTS)  -c dia10b1a.f
	@if test -f dia10b1a.T;then echo "	mv -f dia10b1a.T NORMAL/.";	mv -f dia10b1a.T NORMAL/. ; fi
	mv -f dia10b1a.o NORMAL/dia10b1a.o

NORMAL/dia10q1a.o: dia10q1a.f
	$(FORT) $(FORTOPTS)  -c dia10q1a.f
	@if test -f dia10q1a.T;then echo "	mv -f dia10q1a.T NORMAL/.";	mv -f dia10q1a.T NORMAL/. ; fi
	mv -f dia10q1a.o NORMAL/dia10q1a.o

NORMAL/filtfl1a.o: filtfl1a.f
	$(FORT) $(FORTOPTS)  -c filtfl1a.f
	@if test -f filtfl1a.T;then echo "	mv -f filtfl1a.T NORMAL/.";	mv -f filtfl1a.T NORMAL/. ; fi
	mv -f filtfl1a.o NORMAL/filtfl1a.o

NORMAL/pthadj1a.o: pthadj1a.f
	$(FORT) $(FORTOPTS)  -c pthadj1a.f
	@if test -f pthadj1a.T;then echo "	mv -f pthadj1a.T NORMAL/.";	mv -f pthadj1a.T NORMAL/. ; fi
	mv -f pthadj1a.o NORMAL/pthadj1a.o

NORMAL/setfil1b.o: setfil1b.f
	$(FORT) $(FORTOPTS)  -c setfil1b.f
	@if test -f setfil1b.T;then echo "	mv -f setfil1b.T NORMAL/.";	mv -f setfil1b.T NORMAL/. ; fi
	mv -f setfil1b.o NORMAL/setfil1b.o

NORMAL/thetlq1a.o: thetlq1a.f
	$(FORT) $(FORTOPTS)  -c thetlq1a.f
	@if test -f thetlq1a.T;then echo "	mv -f thetlq1a.T NORMAL/.";	mv -f thetlq1a.T NORMAL/. ; fi
	mv -f thetlq1a.o NORMAL/thetlq1a.o

NORMAL/uvadj1a.o: uvadj1a.f
	$(FORT) $(FORTOPTS)  -c uvadj1a.f
	@if test -f uvadj1a.T;then echo "	mv -f uvadj1a.T NORMAL/.";	mv -f uvadj1a.T NORMAL/. ; fi
	mv -f uvadj1a.o NORMAL/uvadj1a.o

NORMAL/vertve1b.o: vertve1b.f
	$(FORT) $(FORTOPTS)  -c vertve1b.f
	@if test -f vertve1b.T;then echo "	mv -f vertve1b.T NORMAL/.";	mv -f vertve1b.T NORMAL/. ; fi
	mv -f vertve1b.o NORMAL/vertve1b.o

