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
NORMAL/ftsa2b.o NORMAL/solang1a.o NORMAL/solpos1a.o NORMAL/swclop2b.o\
NORMAL/swdkdi1a.o NORMAL/swdtca1a.o NORMAL/swmast2b.o NORMAL/swmsal1a.o\
NORMAL/swptsc1a.o NORMAL/swrad2b.o NORMAL/swtran1a.o

all: $(OBJS)


NORMAL/ftsa2b.o: ftsa2b.f
	$(FORT) $(FORTOPTS)  -c ftsa2b.f
	@if test -f ftsa2b.T;then echo "	mv -f ftsa2b.T NORMAL/.";	mv -f ftsa2b.T NORMAL/. ; fi
	mv -f ftsa2b.o NORMAL/ftsa2b.o

NORMAL/solang1a.o: solang1a.f
	$(FORT) $(FORTOPTS)  -c solang1a.f
	@if test -f solang1a.T;then echo "	mv -f solang1a.T NORMAL/.";	mv -f solang1a.T NORMAL/. ; fi
	mv -f solang1a.o NORMAL/solang1a.o

NORMAL/solpos1a.o: solpos1a.f
	$(FORT) $(FORTOPTS)  -c solpos1a.f
	@if test -f solpos1a.T;then echo "	mv -f solpos1a.T NORMAL/.";	mv -f solpos1a.T NORMAL/. ; fi
	mv -f solpos1a.o NORMAL/solpos1a.o

NORMAL/swclop2b.o: swclop2b.f
	$(FORT) $(FORTOPTS)  -c swclop2b.f
	@if test -f swclop2b.T;then echo "	mv -f swclop2b.T NORMAL/.";	mv -f swclop2b.T NORMAL/. ; fi
	mv -f swclop2b.o NORMAL/swclop2b.o

NORMAL/swdkdi1a.o: swdkdi1a.f
	$(FORT) $(FORTOPTS)  -c swdkdi1a.f
	@if test -f swdkdi1a.T;then echo "	mv -f swdkdi1a.T NORMAL/.";	mv -f swdkdi1a.T NORMAL/. ; fi
	mv -f swdkdi1a.o NORMAL/swdkdi1a.o

NORMAL/swdtca1a.o: swdtca1a.f
	$(FORT) $(FORTOPTS)  -c swdtca1a.f
	@if test -f swdtca1a.T;then echo "	mv -f swdtca1a.T NORMAL/.";	mv -f swdtca1a.T NORMAL/. ; fi
	mv -f swdtca1a.o NORMAL/swdtca1a.o

NORMAL/swmast2b.o: swmast2b.f
	$(FORT) $(FORTOPTS)  -c swmast2b.f
	@if test -f swmast2b.T;then echo "	mv -f swmast2b.T NORMAL/.";	mv -f swmast2b.T NORMAL/. ; fi
	mv -f swmast2b.o NORMAL/swmast2b.o

NORMAL/swmsal1a.o: swmsal1a.f
	$(FORT) $(FORTOPTS)  -c swmsal1a.f
	@if test -f swmsal1a.T;then echo "	mv -f swmsal1a.T NORMAL/.";	mv -f swmsal1a.T NORMAL/. ; fi
	mv -f swmsal1a.o NORMAL/swmsal1a.o

NORMAL/swptsc1a.o: swptsc1a.f
	$(FORT) $(FORTOPTS)  -c swptsc1a.f
	@if test -f swptsc1a.T;then echo "	mv -f swptsc1a.T NORMAL/.";	mv -f swptsc1a.T NORMAL/. ; fi
	mv -f swptsc1a.o NORMAL/swptsc1a.o

NORMAL/swrad2b.o: swrad2b.f
	$(FORT) $(FORTOPTS)  -c swrad2b.f
	@if test -f swrad2b.T;then echo "	mv -f swrad2b.T NORMAL/.";	mv -f swrad2b.T NORMAL/. ; fi
	mv -f swrad2b.o NORMAL/swrad2b.o

NORMAL/swtran1a.o: swtran1a.f
	$(FORT) $(FORTOPTS)  -c swtran1a.f
	@if test -f swtran1a.T;then echo "	mv -f swtran1a.T NORMAL/.";	mv -f swtran1a.T NORMAL/. ; fi
	mv -f swtran1a.o NORMAL/swtran1a.o

