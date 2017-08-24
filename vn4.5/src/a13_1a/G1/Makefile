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
NORMAL/cofthq1a.o NORMAL/cofuv1a.o NORMAL/difctl1a.o NORMAL/filtuv1a.o\
NORMAL/qtpos1a.o NORMAL/thqdif1a.o NORMAL/uvdif1a.o

all: $(OBJS)


NORMAL/cofthq1a.o: cofthq1a.f
	$(FORT) $(FORTOPTS)  -c cofthq1a.f
	@if test -f cofthq1a.T;then echo "	mv -f cofthq1a.T NORMAL/.";	mv -f cofthq1a.T NORMAL/. ; fi
	mv -f cofthq1a.o NORMAL/cofthq1a.o

NORMAL/cofuv1a.o: cofuv1a.f
	$(FORT) $(FORTOPTS)  -c cofuv1a.f
	@if test -f cofuv1a.T;then echo "	mv -f cofuv1a.T NORMAL/.";	mv -f cofuv1a.T NORMAL/. ; fi
	mv -f cofuv1a.o NORMAL/cofuv1a.o

NORMAL/difctl1a.o: difctl1a.f
	$(FORT) $(FORTOPTS)  -c difctl1a.f
	@if test -f difctl1a.T;then echo "	mv -f difctl1a.T NORMAL/.";	mv -f difctl1a.T NORMAL/. ; fi
	mv -f difctl1a.o NORMAL/difctl1a.o

NORMAL/filtuv1a.o: filtuv1a.f
	$(FORT) $(FORTOPTS)  -c filtuv1a.f
	@if test -f filtuv1a.T;then echo "	mv -f filtuv1a.T NORMAL/.";	mv -f filtuv1a.T NORMAL/. ; fi
	mv -f filtuv1a.o NORMAL/filtuv1a.o

NORMAL/qtpos1a.o: qtpos1a.f
	$(FORT) $(FORTOPTS)  -c qtpos1a.f
	@if test -f qtpos1a.T;then echo "	mv -f qtpos1a.T NORMAL/.";	mv -f qtpos1a.T NORMAL/. ; fi
	mv -f qtpos1a.o NORMAL/qtpos1a.o

NORMAL/thqdif1a.o: thqdif1a.f
	$(FORT) $(FORTOPTS)  -c thqdif1a.f
	@if test -f thqdif1a.T;then echo "	mv -f thqdif1a.T NORMAL/.";	mv -f thqdif1a.T NORMAL/. ; fi
	mv -f thqdif1a.o NORMAL/thqdif1a.o

NORMAL/uvdif1a.o: uvdif1a.f
	$(FORT) $(FORTOPTS)  -c uvdif1a.f
	@if test -f uvdif1a.T;then echo "	mv -f uvdif1a.T NORMAL/.";	mv -f uvdif1a.T NORMAL/. ; fi
	mv -f uvdif1a.o NORMAL/uvdif1a.o

