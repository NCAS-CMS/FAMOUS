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
NORMAL/contra1a.o NORMAL/diatha1a.o NORMAL/duct1a.o NORMAL/freeze1a.o\
NORMAL/icaoht1a.o NORMAL/phydia1a.o NORMAL/qctrop1a.o NORMAL/snowpr1a.o\
NORMAL/thetaw1a.o

all: $(OBJS)


NORMAL/contra1a.o: contra1a.f
	$(FORT) $(FORTOPTS)  -c contra1a.f
	@if test -f contra1a.T;then echo "	mv -f contra1a.T NORMAL/.";	mv -f contra1a.T NORMAL/. ; fi
	mv -f contra1a.o NORMAL/contra1a.o

NORMAL/diatha1a.o: diatha1a.f
	$(FORT) $(FORTOPTS)  -c diatha1a.f
	@if test -f diatha1a.T;then echo "	mv -f diatha1a.T NORMAL/.";	mv -f diatha1a.T NORMAL/. ; fi
	mv -f diatha1a.o NORMAL/diatha1a.o

NORMAL/duct1a.o: duct1a.f
	$(FORT) $(FORTOPTS)  -c duct1a.f
	@if test -f duct1a.T;then echo "	mv -f duct1a.T NORMAL/.";	mv -f duct1a.T NORMAL/. ; fi
	mv -f duct1a.o NORMAL/duct1a.o

NORMAL/freeze1a.o: freeze1a.f
	$(FORT) $(FORTOPTS)  -c freeze1a.f
	@if test -f freeze1a.T;then echo "	mv -f freeze1a.T NORMAL/.";	mv -f freeze1a.T NORMAL/. ; fi
	mv -f freeze1a.o NORMAL/freeze1a.o

NORMAL/icaoht1a.o: icaoht1a.f
	$(FORT) $(FORTOPTS)  -c icaoht1a.f
	@if test -f icaoht1a.T;then echo "	mv -f icaoht1a.T NORMAL/.";	mv -f icaoht1a.T NORMAL/. ; fi
	mv -f icaoht1a.o NORMAL/icaoht1a.o

NORMAL/phydia1a.o: phydia1a.f
	$(FORT) $(FORTOPTS)  -c phydia1a.f
	@if test -f phydia1a.T;then echo "	mv -f phydia1a.T NORMAL/.";	mv -f phydia1a.T NORMAL/. ; fi
	mv -f phydia1a.o NORMAL/phydia1a.o

NORMAL/qctrop1a.o: qctrop1a.f
	$(FORT) $(FORTOPTS)  -c qctrop1a.f
	@if test -f qctrop1a.T;then echo "	mv -f qctrop1a.T NORMAL/.";	mv -f qctrop1a.T NORMAL/. ; fi
	mv -f qctrop1a.o NORMAL/qctrop1a.o

NORMAL/snowpr1a.o: snowpr1a.f
	$(FORT) $(FORTOPTS)  -c snowpr1a.f
	@if test -f snowpr1a.T;then echo "	mv -f snowpr1a.T NORMAL/.";	mv -f snowpr1a.T NORMAL/. ; fi
	mv -f snowpr1a.o NORMAL/snowpr1a.o

NORMAL/thetaw1a.o: thetaw1a.f
	$(FORT) $(FORTOPTS)  -c thetaw1a.f
	@if test -f thetaw1a.T;then echo "	mv -f thetaw1a.T NORMAL/.";	mv -f thetaw1a.T NORMAL/. ; fi
	mv -f thetaw1a.o NORMAL/thetaw1a.o

