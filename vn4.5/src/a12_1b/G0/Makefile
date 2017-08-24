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
NORMAL/advctl1a.o NORMAL/advpgd1a.o NORMAL/advugd1a.o NORMAL/divdmp1a.o\
NORMAL/masuvw1a.o NORMAL/masuwt1a.o NORMAL/qtadv1a.o NORMAL/thadv1a.o\
NORMAL/uvadv1b.o NORMAL/vcorio1a.o

all: $(OBJS)


NORMAL/advctl1a.o: advctl1a.f
	$(FORT) $(FORTOPTS)  -c advctl1a.f
	@if test -f advctl1a.T;then echo "	mv -f advctl1a.T NORMAL/.";	mv -f advctl1a.T NORMAL/. ; fi
	mv -f advctl1a.o NORMAL/advctl1a.o

NORMAL/advpgd1a.o: advpgd1a.f
	$(FORT) $(FORTOPTS)  -c advpgd1a.f
	@if test -f advpgd1a.T;then echo "	mv -f advpgd1a.T NORMAL/.";	mv -f advpgd1a.T NORMAL/. ; fi
	mv -f advpgd1a.o NORMAL/advpgd1a.o

NORMAL/advugd1a.o: advugd1a.f
	$(FORT) $(FORTOPTS)  -c advugd1a.f
	@if test -f advugd1a.T;then echo "	mv -f advugd1a.T NORMAL/.";	mv -f advugd1a.T NORMAL/. ; fi
	mv -f advugd1a.o NORMAL/advugd1a.o

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

NORMAL/qtadv1a.o: qtadv1a.f
	$(FORT) $(FORTOPTS)  -c qtadv1a.f
	@if test -f qtadv1a.T;then echo "	mv -f qtadv1a.T NORMAL/.";	mv -f qtadv1a.T NORMAL/. ; fi
	mv -f qtadv1a.o NORMAL/qtadv1a.o

NORMAL/thadv1a.o: thadv1a.f
	$(FORT) $(FORTOPTS)  -c thadv1a.f
	@if test -f thadv1a.T;then echo "	mv -f thadv1a.T NORMAL/.";	mv -f thadv1a.T NORMAL/. ; fi
	mv -f thadv1a.o NORMAL/thadv1a.o

NORMAL/uvadv1b.o: uvadv1b.f
	$(FORT) $(FORTOPTS)  -c uvadv1b.f
	@if test -f uvadv1b.T;then echo "	mv -f uvadv1b.T NORMAL/.";	mv -f uvadv1b.T NORMAL/. ; fi
	mv -f uvadv1b.o NORMAL/uvadv1b.o

NORMAL/vcorio1a.o: vcorio1a.f
	$(FORT) $(FORTOPTS)  -c vcorio1a.f
	@if test -f vcorio1a.T;then echo "	mv -f vcorio1a.T NORMAL/.";	mv -f vcorio1a.T NORMAL/. ; fi
	mv -f vcorio1a.o NORMAL/vcorio1a.o

