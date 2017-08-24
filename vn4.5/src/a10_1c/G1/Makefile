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
NORMAL/adjctl1c.o NORMAL/angmom1a.o NORMAL/dia10a1a.o NORMAL/dia10b1a.o\
NORMAL/dia10q1a.o NORMAL/filtfl1a.o NORMAL/pthadj1a.o NORMAL/setfil1a.o\
NORMAL/thetlq1c.o NORMAL/uvadj1c.o NORMAL/vertve1c.o

all: $(OBJS)


NORMAL/adjctl1c.o: adjctl1c.f
	$(FORT) $(FORTOPTS)  -c adjctl1c.f
	@if test -f adjctl1c.T;then echo "	mv -f adjctl1c.T NORMAL/.";	mv -f adjctl1c.T NORMAL/. ; fi
	mv -f adjctl1c.o NORMAL/adjctl1c.o

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

NORMAL/setfil1a.o: setfil1a.f
	$(FORT) $(FORTOPTS)  -c setfil1a.f
	@if test -f setfil1a.T;then echo "	mv -f setfil1a.T NORMAL/.";	mv -f setfil1a.T NORMAL/. ; fi
	mv -f setfil1a.o NORMAL/setfil1a.o

NORMAL/thetlq1c.o: thetlq1c.f
	$(FORT) $(FORTOPTS)  -c thetlq1c.f
	@if test -f thetlq1c.T;then echo "	mv -f thetlq1c.T NORMAL/.";	mv -f thetlq1c.T NORMAL/. ; fi
	mv -f thetlq1c.o NORMAL/thetlq1c.o

NORMAL/uvadj1c.o: uvadj1c.f
	$(FORT) $(FORTOPTS)  -c uvadj1c.f
	@if test -f uvadj1c.T;then echo "	mv -f uvadj1c.T NORMAL/.";	mv -f uvadj1c.T NORMAL/. ; fi
	mv -f uvadj1c.o NORMAL/uvadj1c.o

NORMAL/vertve1c.o: vertve1c.f
	$(FORT) $(FORTOPTS)  -c vertve1c.f
	@if test -f vertve1c.T;then echo "	mv -f vertve1c.T NORMAL/.";	mv -f vertve1c.T NORMAL/. ; fi
	mv -f vertve1c.o NORMAL/vertve1c.o

