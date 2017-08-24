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
NORMAL/darcy5a.o NORMAL/frunof1a.o NORMAL/heatcp5a.o NORMAL/htcond5a.o\
NORMAL/hyd_ic5a.o NORMAL/hydcon5a.o NORMAL/hydrol5a.o NORMAL/icehtc5a.o\
NORMAL/infilt5a.o NORMAL/sfsnow2a.o NORMAL/sieve1a.o NORMAL/smroot5a.o\
NORMAL/soilht5a.o NORMAL/soilhy5a.o NORMAL/surfhy1a.o

all: $(OBJS)


NORMAL/darcy5a.o: darcy5a.f
	$(FORT) $(FORTOPTS)  -c darcy5a.f
	@if test -f darcy5a.T;then echo "	mv -f darcy5a.T NORMAL/.";	mv -f darcy5a.T NORMAL/. ; fi
	mv -f darcy5a.o NORMAL/darcy5a.o

NORMAL/frunof1a.o: frunof1a.f
	$(FORT) $(FORTOPTS)  -c frunof1a.f
	@if test -f frunof1a.T;then echo "	mv -f frunof1a.T NORMAL/.";	mv -f frunof1a.T NORMAL/. ; fi
	mv -f frunof1a.o NORMAL/frunof1a.o

NORMAL/heatcp5a.o: heatcp5a.f
	$(FORT) $(FORTOPTS)  -c heatcp5a.f
	@if test -f heatcp5a.T;then echo "	mv -f heatcp5a.T NORMAL/.";	mv -f heatcp5a.T NORMAL/. ; fi
	mv -f heatcp5a.o NORMAL/heatcp5a.o

NORMAL/htcond5a.o: htcond5a.f
	$(FORT) $(FORTOPTS)  -c htcond5a.f
	@if test -f htcond5a.T;then echo "	mv -f htcond5a.T NORMAL/.";	mv -f htcond5a.T NORMAL/. ; fi
	mv -f htcond5a.o NORMAL/htcond5a.o

NORMAL/hyd_ic5a.o: hyd_ic5a.f
	$(FORT) $(FORTOPTS)  -c hyd_ic5a.f
	@if test -f hyd_ic5a.T;then echo "	mv -f hyd_ic5a.T NORMAL/.";	mv -f hyd_ic5a.T NORMAL/. ; fi
	mv -f hyd_ic5a.o NORMAL/hyd_ic5a.o

NORMAL/hydcon5a.o: hydcon5a.f
	$(FORT) $(FORTOPTS)  -c hydcon5a.f
	@if test -f hydcon5a.T;then echo "	mv -f hydcon5a.T NORMAL/.";	mv -f hydcon5a.T NORMAL/. ; fi
	mv -f hydcon5a.o NORMAL/hydcon5a.o

NORMAL/hydrol5a.o: hydrol5a.f
	$(FORT) $(FORTOPTS)  -c hydrol5a.f
	@if test -f hydrol5a.T;then echo "	mv -f hydrol5a.T NORMAL/.";	mv -f hydrol5a.T NORMAL/. ; fi
	mv -f hydrol5a.o NORMAL/hydrol5a.o

NORMAL/icehtc5a.o: icehtc5a.f
	$(FORT) $(FORTOPTS)  -c icehtc5a.f
	@if test -f icehtc5a.T;then echo "	mv -f icehtc5a.T NORMAL/.";	mv -f icehtc5a.T NORMAL/. ; fi
	mv -f icehtc5a.o NORMAL/icehtc5a.o

NORMAL/infilt5a.o: infilt5a.f
	$(FORT) $(FORTOPTS)  -c infilt5a.f
	@if test -f infilt5a.T;then echo "	mv -f infilt5a.T NORMAL/.";	mv -f infilt5a.T NORMAL/. ; fi
	mv -f infilt5a.o NORMAL/infilt5a.o

NORMAL/sfsnow2a.o: sfsnow2a.f
	$(FORT) $(FORTOPTS)  -c sfsnow2a.f
	@if test -f sfsnow2a.T;then echo "	mv -f sfsnow2a.T NORMAL/.";	mv -f sfsnow2a.T NORMAL/. ; fi
	mv -f sfsnow2a.o NORMAL/sfsnow2a.o

NORMAL/sieve1a.o: sieve1a.f
	$(FORT) $(FORTOPTS)  -c sieve1a.f
	@if test -f sieve1a.T;then echo "	mv -f sieve1a.T NORMAL/.";	mv -f sieve1a.T NORMAL/. ; fi
	mv -f sieve1a.o NORMAL/sieve1a.o

NORMAL/smroot5a.o: smroot5a.f
	$(FORT) $(FORTOPTS)  -c smroot5a.f
	@if test -f smroot5a.T;then echo "	mv -f smroot5a.T NORMAL/.";	mv -f smroot5a.T NORMAL/. ; fi
	mv -f smroot5a.o NORMAL/smroot5a.o

NORMAL/soilht5a.o: soilht5a.f
	$(FORT) $(FORTOPTS)  -c soilht5a.f
	@if test -f soilht5a.T;then echo "	mv -f soilht5a.T NORMAL/.";	mv -f soilht5a.T NORMAL/. ; fi
	mv -f soilht5a.o NORMAL/soilht5a.o

NORMAL/soilhy5a.o: soilhy5a.f
	$(FORT) $(FORTOPTS)  -c soilhy5a.f
	@if test -f soilhy5a.T;then echo "	mv -f soilhy5a.T NORMAL/.";	mv -f soilhy5a.T NORMAL/. ; fi
	mv -f soilhy5a.o NORMAL/soilhy5a.o

NORMAL/surfhy1a.o: surfhy1a.f
	$(FORT) $(FORTOPTS)  -c surfhy1a.f
	@if test -f surfhy1a.T;then echo "	mv -f surfhy1a.T NORMAL/.";	mv -f surfhy1a.T NORMAL/. ; fi
	mv -f surfhy1a.o NORMAL/surfhy1a.o

