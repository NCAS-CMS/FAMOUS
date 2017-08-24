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
NORMAL/aloghf2a.o NORMAL/exphf2a.o NORMAL/fourie3a.o NORMAL/ieeepk1a.o\
NORMAL/isamax2a.o NORMAL/ismax2a.o NORMAL/ismin2a.o NORMAL/isrchf2a.o\
NORMAL/p_decode.o NORMAL/p_orders.o NORMAL/sqrthf2a.o NORMAL/ssum2a.o

all: $(OBJS)


NORMAL/aloghf2a.o: aloghf2a.f
	$(FORT) $(FORTOPTS)  -c aloghf2a.f
	@if test -f aloghf2a.T;then echo "	mv -f aloghf2a.T NORMAL/.";	mv -f aloghf2a.T NORMAL/. ; fi
	mv -f aloghf2a.o NORMAL/aloghf2a.o

NORMAL/exphf2a.o: exphf2a.f
	$(FORT) $(FORTOPTS)  -c exphf2a.f
	@if test -f exphf2a.T;then echo "	mv -f exphf2a.T NORMAL/.";	mv -f exphf2a.T NORMAL/. ; fi
	mv -f exphf2a.o NORMAL/exphf2a.o

NORMAL/fourie3a.o: fourie3a.f
	$(FORT) $(FORTOPTS)  -c fourie3a.f
	@if test -f fourie3a.T;then echo "	mv -f fourie3a.T NORMAL/.";	mv -f fourie3a.T NORMAL/. ; fi
	mv -f fourie3a.o NORMAL/fourie3a.o

NORMAL/ieeepk1a.o: ieeepk1a.f
	$(FORT) $(FORTOPTS)  -c ieeepk1a.f
	@if test -f ieeepk1a.T;then echo "	mv -f ieeepk1a.T NORMAL/.";	mv -f ieeepk1a.T NORMAL/. ; fi
	mv -f ieeepk1a.o NORMAL/ieeepk1a.o

NORMAL/isamax2a.o: isamax2a.f
	$(FORT) $(FORTOPTS)  -c isamax2a.f
	@if test -f isamax2a.T;then echo "	mv -f isamax2a.T NORMAL/.";	mv -f isamax2a.T NORMAL/. ; fi
	mv -f isamax2a.o NORMAL/isamax2a.o

NORMAL/ismax2a.o: ismax2a.f
	$(FORT) $(FORTOPTS)  -c ismax2a.f
	@if test -f ismax2a.T;then echo "	mv -f ismax2a.T NORMAL/.";	mv -f ismax2a.T NORMAL/. ; fi
	mv -f ismax2a.o NORMAL/ismax2a.o

NORMAL/ismin2a.o: ismin2a.f
	$(FORT) $(FORTOPTS)  -c ismin2a.f
	@if test -f ismin2a.T;then echo "	mv -f ismin2a.T NORMAL/.";	mv -f ismin2a.T NORMAL/. ; fi
	mv -f ismin2a.o NORMAL/ismin2a.o

NORMAL/isrchf2a.o: isrchf2a.f
	$(FORT) $(FORTOPTS)  -c isrchf2a.f
	@if test -f isrchf2a.T;then echo "	mv -f isrchf2a.T NORMAL/.";	mv -f isrchf2a.T NORMAL/. ; fi
	mv -f isrchf2a.o NORMAL/isrchf2a.o

NORMAL/p_decode.o: p_decode.f
	$(FORT) $(FORTOPTS)  -c p_decode.f
	@if test -f p_decode.T;then echo "	mv -f p_decode.T NORMAL/.";	mv -f p_decode.T NORMAL/. ; fi
	mv -f p_decode.o NORMAL/p_decode.o

NORMAL/p_orders.o: p_orders.f
	$(FORT) $(FORTOPTS)  -c p_orders.f
	@if test -f p_orders.T;then echo "	mv -f p_orders.T NORMAL/.";	mv -f p_orders.T NORMAL/. ; fi
	mv -f p_orders.o NORMAL/p_orders.o

NORMAL/sqrthf2a.o: sqrthf2a.f
	$(FORT) $(FORTOPTS)  -c sqrthf2a.f
	@if test -f sqrthf2a.T;then echo "	mv -f sqrthf2a.T NORMAL/.";	mv -f sqrthf2a.T NORMAL/. ; fi
	mv -f sqrthf2a.o NORMAL/sqrthf2a.o

NORMAL/ssum2a.o: ssum2a.f
	$(FORT) $(FORTOPTS)  -c ssum2a.f
	@if test -f ssum2a.T;then echo "	mv -f ssum2a.T NORMAL/.";	mv -f ssum2a.T NORMAL/. ; fi
	mv -f ssum2a.o NORMAL/ssum2a.o

