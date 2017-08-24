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
NORMAL/grvset1a.o NORMAL/soot1a.o NORMAL/sulph1a.o NORMAL/vgrav1a.o

all: $(OBJS)


NORMAL/grvset1a.o: grvset1a.f
	$(FORT) $(FORTOPTS)  -c grvset1a.f
	@if test -f grvset1a.T;then echo "	mv -f grvset1a.T NORMAL/.";	mv -f grvset1a.T NORMAL/. ; fi
	mv -f grvset1a.o NORMAL/grvset1a.o

NORMAL/soot1a.o: soot1a.f
	$(FORT) $(FORTOPTS)  -c soot1a.f
	@if test -f soot1a.T;then echo "	mv -f soot1a.T NORMAL/.";	mv -f soot1a.T NORMAL/. ; fi
	mv -f soot1a.o NORMAL/soot1a.o

NORMAL/sulph1a.o: sulph1a.f
	$(FORT) $(FORTOPTS)  -c sulph1a.f
	@if test -f sulph1a.T;then echo "	mv -f sulph1a.T NORMAL/.";	mv -f sulph1a.T NORMAL/. ; fi
	mv -f sulph1a.o NORMAL/sulph1a.o

NORMAL/vgrav1a.o: vgrav1a.f
	$(FORT) $(FORTOPTS)  -c vgrav1a.f
	@if test -f vgrav1a.T;then echo "	mv -f vgrav1a.T NORMAL/.";	mv -f vgrav1a.T NORMAL/. ; fi
	mv -f vgrav1a.o NORMAL/vgrav1a.o

