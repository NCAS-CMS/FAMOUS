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
NORMAL/advctl1e.o NORMAL/advpgd1e.o NORMAL/advugd1e.o NORMAL/divdmp1a.o\
NORMAL/masuvw1a.o NORMAL/masuwt1a.o NORMAL/thadv1e.o NORMAL/uvadv1e.o\
NORMAL/vcorio1a.o

all: $(OBJS)


NORMAL/advctl1e.o: advctl1e.f
	$(FORT) $(FORTOPTS)  -c advctl1e.f
	@if test -f advctl1e.T;then echo "	mv -f advctl1e.T NORMAL/.";	mv -f advctl1e.T NORMAL/. ; fi
	mv -f advctl1e.o NORMAL/advctl1e.o

NORMAL/advpgd1e.o: advpgd1e.f
	$(FORT) $(FORTOPTS)  -c advpgd1e.f
	@if test -f advpgd1e.T;then echo "	mv -f advpgd1e.T NORMAL/.";	mv -f advpgd1e.T NORMAL/. ; fi
	mv -f advpgd1e.o NORMAL/advpgd1e.o

NORMAL/advugd1e.o: advugd1e.f
	$(FORT) $(FORTOPTS)  -c advugd1e.f
	@if test -f advugd1e.T;then echo "	mv -f advugd1e.T NORMAL/.";	mv -f advugd1e.T NORMAL/. ; fi
	mv -f advugd1e.o NORMAL/advugd1e.o

NORMAL/divdmp1a.o: divdmp1a.f
	$(FORT) $(FORTOPTS)  -c divdmp1a.f
	@if test -f divdmp1a.T;then echo "	mv -f divdmp1a.T NORMAL/.";	mv -f divdmp1a.T NORMAL/. ; fi
	mv -f divdmp1a.o NORMAL/divdmp1a.o

NORMAL/masuvw1a.o: masuvw1a.f
	$(FORT) $(FORTOPTS)  -c masuvw1a.f
	@if test -f masuvw1a.T;then echo "	mv -f masuvw1a.T NORMAL/.";	mv -f masuvw1a.T NORMAL/. ; fi
	mv -f masuvw1a.o NORMAL/masuvw1a.o

NORMAL/masuwt1a.o: masuwt1a.f
	$(FORT) $(FORTOPTS)  -c masuwt1a.f
	@if test -f masuwt1a.T;then echo "	mv -f masuwt1a.T NORMAL/.";	mv -f masuwt1a.T NORMAL/. ; fi
	mv -f masuwt1a.o NORMAL/masuwt1a.o

NORMAL/thadv1e.o: thadv1e.f
	$(FORT) $(FORTOPTS)  -c thadv1e.f
	@if test -f thadv1e.T;then echo "	mv -f thadv1e.T NORMAL/.";	mv -f thadv1e.T NORMAL/. ; fi
	mv -f thadv1e.o NORMAL/thadv1e.o

NORMAL/uvadv1e.o: uvadv1e.f
	$(FORT) $(FORTOPTS)  -c uvadv1e.f
	@if test -f uvadv1e.T;then echo "	mv -f uvadv1e.T NORMAL/.";	mv -f uvadv1e.T NORMAL/. ; fi
	mv -f uvadv1e.o NORMAL/uvadv1e.o

NORMAL/vcorio1a.o: vcorio1a.f
	$(FORT) $(FORTOPTS)  -c vcorio1a.f
	@if test -f vcorio1a.T;then echo "	mv -f vcorio1a.T NORMAL/.";	mv -f vcorio1a.T NORMAL/. ; fi
	mv -f vcorio1a.o NORMAL/vcorio1a.o

