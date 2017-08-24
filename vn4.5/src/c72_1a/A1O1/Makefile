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
NORMAL/chk1box1.o NORMAL/copya2o1.o NORMAL/copyo2a1.o NORMAL/cyclcbc1.o\
NORMAL/doarav1.o NORMAL/postint1.o NORMAL/prearav1.o NORMAL/rowswap1.o

all: $(OBJS)


NORMAL/chk1box1.o: chk1box1.f
	$(FORT) $(FORTOPTS)  -c chk1box1.f
	@if test -f chk1box1.T;then echo "	mv -f chk1box1.T NORMAL/.";	mv -f chk1box1.T NORMAL/. ; fi
	mv -f chk1box1.o NORMAL/chk1box1.o

NORMAL/copya2o1.o: copya2o1.f
	$(FORT) $(FORTOPTS)  -c copya2o1.f
	@if test -f copya2o1.T;then echo "	mv -f copya2o1.T NORMAL/.";	mv -f copya2o1.T NORMAL/. ; fi
	mv -f copya2o1.o NORMAL/copya2o1.o

NORMAL/copyo2a1.o: copyo2a1.f
	$(FORT) $(FORTOPTS)  -c copyo2a1.f
	@if test -f copyo2a1.T;then echo "	mv -f copyo2a1.T NORMAL/.";	mv -f copyo2a1.T NORMAL/. ; fi
	mv -f copyo2a1.o NORMAL/copyo2a1.o

NORMAL/cyclcbc1.o: cyclcbc1.f
	$(FORT) $(FORTOPTS)  -c cyclcbc1.f
	@if test -f cyclcbc1.T;then echo "	mv -f cyclcbc1.T NORMAL/.";	mv -f cyclcbc1.T NORMAL/. ; fi
	mv -f cyclcbc1.o NORMAL/cyclcbc1.o

NORMAL/doarav1.o: doarav1.f
	$(FORT) $(FORTOPTS)  -c doarav1.f
	@if test -f doarav1.T;then echo "	mv -f doarav1.T NORMAL/.";	mv -f doarav1.T NORMAL/. ; fi
	mv -f doarav1.o NORMAL/doarav1.o

NORMAL/postint1.o: postint1.f
	$(FORT) $(FORTOPTS)  -c postint1.f
	@if test -f postint1.T;then echo "	mv -f postint1.T NORMAL/.";	mv -f postint1.T NORMAL/. ; fi
	mv -f postint1.o NORMAL/postint1.o

NORMAL/prearav1.o: prearav1.f
	$(FORT) $(FORTOPTS)  -c prearav1.f
	@if test -f prearav1.T;then echo "	mv -f prearav1.T NORMAL/.";	mv -f prearav1.T NORMAL/. ; fi
	mv -f prearav1.o NORMAL/prearav1.o

NORMAL/rowswap1.o: rowswap1.f
	$(FORT) $(FORTOPTS)  -c rowswap1.f
	@if test -f rowswap1.T;then echo "	mv -f rowswap1.T NORMAL/.";	mv -f rowswap1.T NORMAL/. ; fi
	mv -f rowswap1.o NORMAL/rowswap1.o

