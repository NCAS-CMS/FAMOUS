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
NORMAL/calcrs1a.o NORMAL/calcts1a.o NORMAL/chgdim1a.o NORMAL/copd3d1a.o\
NORMAL/copdia1a.o NORMAL/dewpnt1a.o NORMAL/extdia1a.o NORMAL/polar1a.o\
NORMAL/polauv1a.o NORMAL/ptocu1a.o NORMAL/ptocv1a.o NORMAL/ptouv1a.o\
NORMAL/qsat2b.o NORMAL/twbulb1a.o NORMAL/uvtop1a.o NORMAL/uvtopf1a.o\
NORMAL/visbty1a.o

all: $(OBJS)


NORMAL/calcrs1a.o: calcrs1a.f
	$(FORT) $(FORTOPTS)  -c calcrs1a.f
	@if test -f calcrs1a.T;then echo "	mv -f calcrs1a.T NORMAL/.";	mv -f calcrs1a.T NORMAL/. ; fi
	mv -f calcrs1a.o NORMAL/calcrs1a.o

NORMAL/calcts1a.o: calcts1a.f
	$(FORT) $(FORTOPTS)  -c calcts1a.f
	@if test -f calcts1a.T;then echo "	mv -f calcts1a.T NORMAL/.";	mv -f calcts1a.T NORMAL/. ; fi
	mv -f calcts1a.o NORMAL/calcts1a.o

NORMAL/chgdim1a.o: chgdim1a.f
	$(FORT) $(FORTOPTS)  -c chgdim1a.f
	@if test -f chgdim1a.T;then echo "	mv -f chgdim1a.T NORMAL/.";	mv -f chgdim1a.T NORMAL/. ; fi
	mv -f chgdim1a.o NORMAL/chgdim1a.o

NORMAL/copd3d1a.o: copd3d1a.f
	$(FORT) $(FORTOPTS)  -c copd3d1a.f
	@if test -f copd3d1a.T;then echo "	mv -f copd3d1a.T NORMAL/.";	mv -f copd3d1a.T NORMAL/. ; fi
	mv -f copd3d1a.o NORMAL/copd3d1a.o

NORMAL/copdia1a.o: copdia1a.f
	$(FORT) $(FORTOPTS)  -c copdia1a.f
	@if test -f copdia1a.T;then echo "	mv -f copdia1a.T NORMAL/.";	mv -f copdia1a.T NORMAL/. ; fi
	mv -f copdia1a.o NORMAL/copdia1a.o

NORMAL/dewpnt1a.o: dewpnt1a.f
	$(FORT) $(FORTOPTS)  -c dewpnt1a.f
	@if test -f dewpnt1a.T;then echo "	mv -f dewpnt1a.T NORMAL/.";	mv -f dewpnt1a.T NORMAL/. ; fi
	mv -f dewpnt1a.o NORMAL/dewpnt1a.o

NORMAL/extdia1a.o: extdia1a.f
	$(FORT) $(FORTOPTS)  -c extdia1a.f
	@if test -f extdia1a.T;then echo "	mv -f extdia1a.T NORMAL/.";	mv -f extdia1a.T NORMAL/. ; fi
	mv -f extdia1a.o NORMAL/extdia1a.o

NORMAL/polar1a.o: polar1a.f
	$(FORT) $(FORTOPTS)  -c polar1a.f
	@if test -f polar1a.T;then echo "	mv -f polar1a.T NORMAL/.";	mv -f polar1a.T NORMAL/. ; fi
	mv -f polar1a.o NORMAL/polar1a.o

NORMAL/polauv1a.o: polauv1a.f
	$(FORT) $(FORTOPTS)  -c polauv1a.f
	@if test -f polauv1a.T;then echo "	mv -f polauv1a.T NORMAL/.";	mv -f polauv1a.T NORMAL/. ; fi
	mv -f polauv1a.o NORMAL/polauv1a.o

NORMAL/ptocu1a.o: ptocu1a.f
	$(FORT) $(FORTOPTS)  -c ptocu1a.f
	@if test -f ptocu1a.T;then echo "	mv -f ptocu1a.T NORMAL/.";	mv -f ptocu1a.T NORMAL/. ; fi
	mv -f ptocu1a.o NORMAL/ptocu1a.o

NORMAL/ptocv1a.o: ptocv1a.f
	$(FORT) $(FORTOPTS)  -c ptocv1a.f
	@if test -f ptocv1a.T;then echo "	mv -f ptocv1a.T NORMAL/.";	mv -f ptocv1a.T NORMAL/. ; fi
	mv -f ptocv1a.o NORMAL/ptocv1a.o

NORMAL/ptouv1a.o: ptouv1a.f
	$(FORT) $(FORTOPTS)  -c ptouv1a.f
	@if test -f ptouv1a.T;then echo "	mv -f ptouv1a.T NORMAL/.";	mv -f ptouv1a.T NORMAL/. ; fi
	mv -f ptouv1a.o NORMAL/ptouv1a.o

NORMAL/qsat2b.o: qsat2b.f
	$(FORT) $(FORTOPTS)  -c qsat2b.f
	@if test -f qsat2b.T;then echo "	mv -f qsat2b.T NORMAL/.";	mv -f qsat2b.T NORMAL/. ; fi
	mv -f qsat2b.o NORMAL/qsat2b.o

NORMAL/twbulb1a.o: twbulb1a.f
	$(FORT) $(FORTOPTS)  -c twbulb1a.f
	@if test -f twbulb1a.T;then echo "	mv -f twbulb1a.T NORMAL/.";	mv -f twbulb1a.T NORMAL/. ; fi
	mv -f twbulb1a.o NORMAL/twbulb1a.o

NORMAL/uvtop1a.o: uvtop1a.f
	$(FORT) $(FORTOPTS)  -c uvtop1a.f
	@if test -f uvtop1a.T;then echo "	mv -f uvtop1a.T NORMAL/.";	mv -f uvtop1a.T NORMAL/. ; fi
	mv -f uvtop1a.o NORMAL/uvtop1a.o

NORMAL/uvtopf1a.o: uvtopf1a.f
	$(FORT) $(FORTOPTS)  -c uvtopf1a.f
	@if test -f uvtopf1a.T;then echo "	mv -f uvtopf1a.T NORMAL/.";	mv -f uvtopf1a.T NORMAL/. ; fi
	mv -f uvtopf1a.o NORMAL/uvtopf1a.o

NORMAL/visbty1a.o: visbty1a.f
	$(FORT) $(FORTOPTS)  -c visbty1a.f
	@if test -f visbty1a.T;then echo "	mv -f visbty1a.T NORMAL/.";	mv -f visbty1a.T NORMAL/. ; fi
	mv -f visbty1a.o NORMAL/visbty1a.o

