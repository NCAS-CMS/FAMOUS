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
NORMAL/cofthq1c.o NORMAL/cofuv1c.o NORMAL/difctl1c.o NORMAL/filtuv1c.o\
NORMAL/qtpos1a.o NORMAL/thqdif1c.o NORMAL/uvdif1c.o

all: $(OBJS)


NORMAL/cofthq1c.o: cofthq1c.f
	$(FORT) $(FORTOPTS)  -c cofthq1c.f
	@if test -f cofthq1c.T;then echo "	mv -f cofthq1c.T NORMAL/.";	mv -f cofthq1c.T NORMAL/. ; fi
	mv -f cofthq1c.o NORMAL/cofthq1c.o

NORMAL/cofuv1c.o: cofuv1c.f
	$(FORT) $(FORTOPTS)  -c cofuv1c.f
	@if test -f cofuv1c.T;then echo "	mv -f cofuv1c.T NORMAL/.";	mv -f cofuv1c.T NORMAL/. ; fi
	mv -f cofuv1c.o NORMAL/cofuv1c.o

NORMAL/difctl1c.o: difctl1c.f
	$(FORT) $(FORTOPTS)  -c difctl1c.f
	@if test -f difctl1c.T;then echo "	mv -f difctl1c.T NORMAL/.";	mv -f difctl1c.T NORMAL/. ; fi
	mv -f difctl1c.o NORMAL/difctl1c.o

NORMAL/filtuv1c.o: filtuv1c.f
	$(FORT) $(FORTOPTS)  -c filtuv1c.f
	@if test -f filtuv1c.T;then echo "	mv -f filtuv1c.T NORMAL/.";	mv -f filtuv1c.T NORMAL/. ; fi
	mv -f filtuv1c.o NORMAL/filtuv1c.o

NORMAL/qtpos1a.o: qtpos1a.f
	$(FORT) $(FORTOPTS)  -c qtpos1a.f
	@if test -f qtpos1a.T;then echo "	mv -f qtpos1a.T NORMAL/.";	mv -f qtpos1a.T NORMAL/. ; fi
	mv -f qtpos1a.o NORMAL/qtpos1a.o

NORMAL/thqdif1c.o: thqdif1c.f
	$(FORT) $(FORTOPTS)  -c thqdif1c.f
	@if test -f thqdif1c.T;then echo "	mv -f thqdif1c.T NORMAL/.";	mv -f thqdif1c.T NORMAL/. ; fi
	mv -f thqdif1c.o NORMAL/thqdif1c.o

NORMAL/uvdif1c.o: uvdif1c.f
	$(FORT) $(FORTOPTS)  -c uvdif1c.f
	@if test -f uvdif1c.T;then echo "	mv -f uvdif1c.T NORMAL/.";	mv -f uvdif1c.T NORMAL/. ; fi
	mv -f uvdif1c.o NORMAL/uvdif1c.o

