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
NORMAL/darcy5a.o NORMAL/frunof7a.o NORMAL/htcond5a.o NORMAL/hyd_ic7a.o\
NORMAL/hydcon5a.o NORMAL/hydrol7a.o NORMAL/icehtc5a.o NORMAL/infilt5a.o\
NORMAL/sfsnow7a.o NORMAL/sieve7a.o NORMAL/soilht7a.o NORMAL/soilhy5a.o\
NORMAL/soilmc7a.o NORMAL/surfhy7a.o

all: $(OBJS)


NORMAL/darcy5a.o: darcy5a.f
	$(FORT) $(FORTOPTS)  -c darcy5a.f
	@if test -f darcy5a.T;then echo "	mv -f darcy5a.T NORMAL/.";	mv -f darcy5a.T NORMAL/. ; fi
	mv -f darcy5a.o NORMAL/darcy5a.o

NORMAL/frunof7a.o: frunof7a.f
	$(FORT) $(FORTOPTS)  -c frunof7a.f
	@if test -f frunof7a.T;then echo "	mv -f frunof7a.T NORMAL/.";	mv -f frunof7a.T NORMAL/. ; fi
	mv -f frunof7a.o NORMAL/frunof7a.o

NORMAL/htcond5a.o: htcond5a.f
	$(FORT) $(FORTOPTS)  -c htcond5a.f
	@if test -f htcond5a.T;then echo "	mv -f htcond5a.T NORMAL/.";	mv -f htcond5a.T NORMAL/. ; fi
	mv -f htcond5a.o NORMAL/htcond5a.o

NORMAL/hyd_ic7a.o: hyd_ic7a.f
	$(FORT) $(FORTOPTS)  -c hyd_ic7a.f
	@if test -f hyd_ic7a.T;then echo "	mv -f hyd_ic7a.T NORMAL/.";	mv -f hyd_ic7a.T NORMAL/. ; fi
	mv -f hyd_ic7a.o NORMAL/hyd_ic7a.o

NORMAL/hydcon5a.o: hydcon5a.f
	$(FORT) $(FORTOPTS)  -c hydcon5a.f
	@if test -f hydcon5a.T;then echo "	mv -f hydcon5a.T NORMAL/.";	mv -f hydcon5a.T NORMAL/. ; fi
	mv -f hydcon5a.o NORMAL/hydcon5a.o

NORMAL/hydrol7a.o: hydrol7a.f
	$(FORT) $(FORTOPTS)  -c hydrol7a.f
	@if test -f hydrol7a.T;then echo "	mv -f hydrol7a.T NORMAL/.";	mv -f hydrol7a.T NORMAL/. ; fi
	mv -f hydrol7a.o NORMAL/hydrol7a.o

NORMAL/icehtc5a.o: icehtc5a.f
	$(FORT) $(FORTOPTS)  -c icehtc5a.f
	@if test -f icehtc5a.T;then echo "	mv -f icehtc5a.T NORMAL/.";	mv -f icehtc5a.T NORMAL/. ; fi
	mv -f icehtc5a.o NORMAL/icehtc5a.o

NORMAL/infilt5a.o: infilt5a.f
	$(FORT) $(FORTOPTS)  -c infilt5a.f
	@if test -f infilt5a.T;then echo "	mv -f infilt5a.T NORMAL/.";	mv -f infilt5a.T NORMAL/. ; fi
	mv -f infilt5a.o NORMAL/infilt5a.o

NORMAL/sfsnow7a.o: sfsnow7a.f
	$(FORT) $(FORTOPTS)  -c sfsnow7a.f
	@if test -f sfsnow7a.T;then echo "	mv -f sfsnow7a.T NORMAL/.";	mv -f sfsnow7a.T NORMAL/. ; fi
	mv -f sfsnow7a.o NORMAL/sfsnow7a.o

NORMAL/sieve7a.o: sieve7a.f
	$(FORT) $(FORTOPTS)  -c sieve7a.f
	@if test -f sieve7a.T;then echo "	mv -f sieve7a.T NORMAL/.";	mv -f sieve7a.T NORMAL/. ; fi
	mv -f sieve7a.o NORMAL/sieve7a.o

NORMAL/soilht7a.o: soilht7a.f
	$(FORT) $(FORTOPTS)  -c soilht7a.f
	@if test -f soilht7a.T;then echo "	mv -f soilht7a.T NORMAL/.";	mv -f soilht7a.T NORMAL/. ; fi
	mv -f soilht7a.o NORMAL/soilht7a.o

NORMAL/soilhy5a.o: soilhy5a.f
	$(FORT) $(FORTOPTS)  -c soilhy5a.f
	@if test -f soilhy5a.T;then echo "	mv -f soilhy5a.T NORMAL/.";	mv -f soilhy5a.T NORMAL/. ; fi
	mv -f soilhy5a.o NORMAL/soilhy5a.o

NORMAL/soilmc7a.o: soilmc7a.f
	$(FORT) $(FORTOPTS)  -c soilmc7a.f
	@if test -f soilmc7a.T;then echo "	mv -f soilmc7a.T NORMAL/.";	mv -f soilmc7a.T NORMAL/. ; fi
	mv -f soilmc7a.o NORMAL/soilmc7a.o

NORMAL/surfhy7a.o: surfhy7a.f
	$(FORT) $(FORTOPTS)  -c surfhy7a.f
	@if test -f surfhy7a.T;then echo "	mv -f surfhy7a.T NORMAL/.";	mv -f surfhy7a.T NORMAL/. ; fi
	mv -f surfhy7a.o NORMAL/surfhy7a.o

