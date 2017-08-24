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
NORMAL/findpic.o NORMAL/fldout.o NORMAL/regroup.o NORMAL/trhou4.o\
NORMAL/vagdirt.o NORMAL/vintpar.o NORMAL/vtott.o NORMAL/waveh.o\
NORMAL/wavetr.o NORMAL/wtrain.o NORMAL/wtrain1.o NORMAL/wtrain2.o\
NORMAL/wtreorg.o

all: $(OBJS)


NORMAL/findpic.o: findpic.f
	$(FORT) $(FORTOPTS)  -c findpic.f
	@if test -f findpic.T;then echo "	mv -f findpic.T NORMAL/.";	mv -f findpic.T NORMAL/. ; fi
	mv -f findpic.o NORMAL/findpic.o

NORMAL/fldout.o: fldout.f
	$(FORT) $(FORTOPTS)  -c fldout.f
	@if test -f fldout.T;then echo "	mv -f fldout.T NORMAL/.";	mv -f fldout.T NORMAL/. ; fi
	mv -f fldout.o NORMAL/fldout.o

NORMAL/regroup.o: regroup.f
	$(FORT) $(FORTOPTS)  -c regroup.f
	@if test -f regroup.T;then echo "	mv -f regroup.T NORMAL/.";	mv -f regroup.T NORMAL/. ; fi
	mv -f regroup.o NORMAL/regroup.o

NORMAL/trhou4.o: trhou4.f
	$(FORT) $(FORTOPTS)  -c trhou4.f
	@if test -f trhou4.T;then echo "	mv -f trhou4.T NORMAL/.";	mv -f trhou4.T NORMAL/. ; fi
	mv -f trhou4.o NORMAL/trhou4.o

NORMAL/vagdirt.o: vagdirt.f
	$(FORT) $(FORTOPTS)  -c vagdirt.f
	@if test -f vagdirt.T;then echo "	mv -f vagdirt.T NORMAL/.";	mv -f vagdirt.T NORMAL/. ; fi
	mv -f vagdirt.o NORMAL/vagdirt.o

NORMAL/vintpar.o: vintpar.f
	$(FORT) $(FORTOPTS)  -c vintpar.f
	@if test -f vintpar.T;then echo "	mv -f vintpar.T NORMAL/.";	mv -f vintpar.T NORMAL/. ; fi
	mv -f vintpar.o NORMAL/vintpar.o

NORMAL/vtott.o: vtott.f
	$(FORT) $(FORTOPTS)  -c vtott.f
	@if test -f vtott.T;then echo "	mv -f vtott.T NORMAL/.";	mv -f vtott.T NORMAL/. ; fi
	mv -f vtott.o NORMAL/vtott.o

NORMAL/waveh.o: waveh.f
	$(FORT) $(FORTOPTS)  -c waveh.f
	@if test -f waveh.T;then echo "	mv -f waveh.T NORMAL/.";	mv -f waveh.T NORMAL/. ; fi
	mv -f waveh.o NORMAL/waveh.o

NORMAL/wavetr.o: wavetr.f
	$(FORT) $(FORTOPTS)  -c wavetr.f
	@if test -f wavetr.T;then echo "	mv -f wavetr.T NORMAL/.";	mv -f wavetr.T NORMAL/. ; fi
	mv -f wavetr.o NORMAL/wavetr.o

NORMAL/wtrain.o: wtrain.f
	$(FORT) $(FORTOPTS)  -c wtrain.f
	@if test -f wtrain.T;then echo "	mv -f wtrain.T NORMAL/.";	mv -f wtrain.T NORMAL/. ; fi
	mv -f wtrain.o NORMAL/wtrain.o

NORMAL/wtrain1.o: wtrain1.f
	$(FORT) $(FORTOPTS)  -c wtrain1.f
	@if test -f wtrain1.T;then echo "	mv -f wtrain1.T NORMAL/.";	mv -f wtrain1.T NORMAL/. ; fi
	mv -f wtrain1.o NORMAL/wtrain1.o

NORMAL/wtrain2.o: wtrain2.f
	$(FORT) $(FORTOPTS)  -c wtrain2.f
	@if test -f wtrain2.T;then echo "	mv -f wtrain2.T NORMAL/.";	mv -f wtrain2.T NORMAL/. ; fi
	mv -f wtrain2.o NORMAL/wtrain2.o

NORMAL/wtreorg.o: wtreorg.f
	$(FORT) $(FORTOPTS)  -c wtreorg.f
	@if test -f wtreorg.T;then echo "	mv -f wtreorg.T NORMAL/.";	mv -f wtreorg.T NORMAL/. ; fi
	mv -f wtreorg.o NORMAL/wtreorg.o

