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
NORMAL/colm1a.o NORMAL/d2lamp1a.o NORMAL/del2p1a.o NORMAL/del2u1a.o\
NORMAL/evalsp1a.o NORMAL/filter1a.o NORMAL/glbm1a.o NORMAL/hgrad1a.o\
NORMAL/matinv1a.o NORMAL/mltgrd1a.o NORMAL/setf1a.o NORMAL/spline1a.o\
NORMAL/zonm1a.o

all: $(OBJS)


NORMAL/colm1a.o: colm1a.f
	$(FORT) $(FORTOPTS)  -c colm1a.f
	@if test -f colm1a.T;then echo "	mv -f colm1a.T NORMAL/.";	mv -f colm1a.T NORMAL/. ; fi
	mv -f colm1a.o NORMAL/colm1a.o

NORMAL/d2lamp1a.o: d2lamp1a.f
	$(FORT) $(FORTOPTS)  -c d2lamp1a.f
	@if test -f d2lamp1a.T;then echo "	mv -f d2lamp1a.T NORMAL/.";	mv -f d2lamp1a.T NORMAL/. ; fi
	mv -f d2lamp1a.o NORMAL/d2lamp1a.o

NORMAL/del2p1a.o: del2p1a.f
	$(FORT) $(FORTOPTS)  -c del2p1a.f
	@if test -f del2p1a.T;then echo "	mv -f del2p1a.T NORMAL/.";	mv -f del2p1a.T NORMAL/. ; fi
	mv -f del2p1a.o NORMAL/del2p1a.o

NORMAL/del2u1a.o: del2u1a.f
	$(FORT) $(FORTOPTS)  -c del2u1a.f
	@if test -f del2u1a.T;then echo "	mv -f del2u1a.T NORMAL/.";	mv -f del2u1a.T NORMAL/. ; fi
	mv -f del2u1a.o NORMAL/del2u1a.o

NORMAL/evalsp1a.o: evalsp1a.f
	$(FORT) $(FORTOPTS)  -c evalsp1a.f
	@if test -f evalsp1a.T;then echo "	mv -f evalsp1a.T NORMAL/.";	mv -f evalsp1a.T NORMAL/. ; fi
	mv -f evalsp1a.o NORMAL/evalsp1a.o

NORMAL/filter1a.o: filter1a.f
	$(FORT) $(FORTOPTS)  -c filter1a.f
	@if test -f filter1a.T;then echo "	mv -f filter1a.T NORMAL/.";	mv -f filter1a.T NORMAL/. ; fi
	mv -f filter1a.o NORMAL/filter1a.o

NORMAL/glbm1a.o: glbm1a.f
	$(FORT) $(FORTOPTS)  -c glbm1a.f
	@if test -f glbm1a.T;then echo "	mv -f glbm1a.T NORMAL/.";	mv -f glbm1a.T NORMAL/. ; fi
	mv -f glbm1a.o NORMAL/glbm1a.o

NORMAL/hgrad1a.o: hgrad1a.f
	$(FORT) $(FORTOPTS)  -c hgrad1a.f
	@if test -f hgrad1a.T;then echo "	mv -f hgrad1a.T NORMAL/.";	mv -f hgrad1a.T NORMAL/. ; fi
	mv -f hgrad1a.o NORMAL/hgrad1a.o

NORMAL/matinv1a.o: matinv1a.f
	$(FORT) $(FORTOPTS)  -c matinv1a.f
	@if test -f matinv1a.T;then echo "	mv -f matinv1a.T NORMAL/.";	mv -f matinv1a.T NORMAL/. ; fi
	mv -f matinv1a.o NORMAL/matinv1a.o

NORMAL/mltgrd1a.o: mltgrd1a.f
	$(FORT) $(FORTOPTS)  -c mltgrd1a.f
	@if test -f mltgrd1a.T;then echo "	mv -f mltgrd1a.T NORMAL/.";	mv -f mltgrd1a.T NORMAL/. ; fi
	mv -f mltgrd1a.o NORMAL/mltgrd1a.o

NORMAL/setf1a.o: setf1a.f
	$(FORT) $(FORTOPTS)  -c setf1a.f
	@if test -f setf1a.T;then echo "	mv -f setf1a.T NORMAL/.";	mv -f setf1a.T NORMAL/. ; fi
	mv -f setf1a.o NORMAL/setf1a.o

NORMAL/spline1a.o: spline1a.f
	$(FORT) $(FORTOPTS)  -c spline1a.f
	@if test -f spline1a.T;then echo "	mv -f spline1a.T NORMAL/.";	mv -f spline1a.T NORMAL/. ; fi
	mv -f spline1a.o NORMAL/spline1a.o

NORMAL/zonm1a.o: zonm1a.f
	$(FORT) $(FORTOPTS)  -c zonm1a.f
	@if test -f zonm1a.T;then echo "	mv -f zonm1a.T NORMAL/.";	mv -f zonm1a.T NORMAL/. ; fi
	mv -f zonm1a.o NORMAL/zonm1a.o

