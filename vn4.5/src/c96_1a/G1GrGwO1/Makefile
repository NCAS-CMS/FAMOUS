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
NORMAL/chdcmp1a.o NORMAL/cpyfld1a.o NORMAL/decmpo1a.o NORMAL/decomp1a.o\
NORMAL/flhalo1a.o NORMAL/gatfld1a.o NORMAL/gatzon1a.o NORMAL/gengat1a.o\
NORMAL/gensct1a.o NORMAL/gltolo1.o NORMAL/gtalbc1a.o NORMAL/gtolbc1a.o\
NORMAL/ieeepk1a.o NORMAL/mppmfl1a.o NORMAL/pbfin1a.o NORMAL/pbfout1a.o\
NORMAL/pclose1a.o NORMAL/popen1a.o NORMAL/pstpos1a.o NORMAL/rdlsm1a.o\
NORMAL/rdmult1a.o NORMAL/scalbc1a.o NORMAL/scolbc1a.o NORMAL/sctfld1a.o\
NORMAL/sctzon1a.o NORMAL/setneb1a.o NORMAL/stsfld1a.o NORMAL/stside1a.o\
NORMAL/swpbnd1a.o NORMAL/wtmult1a.o

all: $(OBJS)


NORMAL/chdcmp1a.o: chdcmp1a.f
	$(FORT) $(FORTOPTS)  -c chdcmp1a.f
	@if test -f chdcmp1a.T;then echo "	mv -f chdcmp1a.T NORMAL/.";	mv -f chdcmp1a.T NORMAL/. ; fi
	mv -f chdcmp1a.o NORMAL/chdcmp1a.o

NORMAL/cpyfld1a.o: cpyfld1a.f
	$(FORT) $(FORTOPTS)  -c cpyfld1a.f
	@if test -f cpyfld1a.T;then echo "	mv -f cpyfld1a.T NORMAL/.";	mv -f cpyfld1a.T NORMAL/. ; fi
	mv -f cpyfld1a.o NORMAL/cpyfld1a.o

NORMAL/decmpo1a.o: decmpo1a.f
	$(FORT) $(FORTOPTS)  -c decmpo1a.f
	@if test -f decmpo1a.T;then echo "	mv -f decmpo1a.T NORMAL/.";	mv -f decmpo1a.T NORMAL/. ; fi
	mv -f decmpo1a.o NORMAL/decmpo1a.o

NORMAL/decomp1a.o: decomp1a.f
	$(FORT) $(FORTOPTS)  -c decomp1a.f
	@if test -f decomp1a.T;then echo "	mv -f decomp1a.T NORMAL/.";	mv -f decomp1a.T NORMAL/. ; fi
	mv -f decomp1a.o NORMAL/decomp1a.o

NORMAL/flhalo1a.o: flhalo1a.f
	$(FORT) $(FORTOPTS)  -c flhalo1a.f
	@if test -f flhalo1a.T;then echo "	mv -f flhalo1a.T NORMAL/.";	mv -f flhalo1a.T NORMAL/. ; fi
	mv -f flhalo1a.o NORMAL/flhalo1a.o

NORMAL/gatfld1a.o: gatfld1a.f
	$(FORT) $(FORTOPTS)  -c gatfld1a.f
	@if test -f gatfld1a.T;then echo "	mv -f gatfld1a.T NORMAL/.";	mv -f gatfld1a.T NORMAL/. ; fi
	mv -f gatfld1a.o NORMAL/gatfld1a.o

NORMAL/gatzon1a.o: gatzon1a.f
	$(FORT) $(FORTOPTS)  -c gatzon1a.f
	@if test -f gatzon1a.T;then echo "	mv -f gatzon1a.T NORMAL/.";	mv -f gatzon1a.T NORMAL/. ; fi
	mv -f gatzon1a.o NORMAL/gatzon1a.o

NORMAL/gengat1a.o: gengat1a.f
	$(FORT) $(FORTOPTS)  -c gengat1a.f
	@if test -f gengat1a.T;then echo "	mv -f gengat1a.T NORMAL/.";	mv -f gengat1a.T NORMAL/. ; fi
	mv -f gengat1a.o NORMAL/gengat1a.o

NORMAL/gensct1a.o: gensct1a.f
	$(FORT) $(FORTOPTS)  -c gensct1a.f
	@if test -f gensct1a.T;then echo "	mv -f gensct1a.T NORMAL/.";	mv -f gensct1a.T NORMAL/. ; fi
	mv -f gensct1a.o NORMAL/gensct1a.o

NORMAL/gltolo1.o: gltolo1.f
	$(FORT) $(FORTOPTS)  -c gltolo1.f
	@if test -f gltolo1.T;then echo "	mv -f gltolo1.T NORMAL/.";	mv -f gltolo1.T NORMAL/. ; fi
	mv -f gltolo1.o NORMAL/gltolo1.o

NORMAL/gtalbc1a.o: gtalbc1a.f
	$(FORT) $(FORTOPTS)  -c gtalbc1a.f
	@if test -f gtalbc1a.T;then echo "	mv -f gtalbc1a.T NORMAL/.";	mv -f gtalbc1a.T NORMAL/. ; fi
	mv -f gtalbc1a.o NORMAL/gtalbc1a.o

NORMAL/gtolbc1a.o: gtolbc1a.f
	$(FORT) $(FORTOPTS)  -c gtolbc1a.f
	@if test -f gtolbc1a.T;then echo "	mv -f gtolbc1a.T NORMAL/.";	mv -f gtolbc1a.T NORMAL/. ; fi
	mv -f gtolbc1a.o NORMAL/gtolbc1a.o

NORMAL/ieeepk1a.o: ieeepk1a.f
	$(FORT) $(FORTOPTS)  -c ieeepk1a.f
	@if test -f ieeepk1a.T;then echo "	mv -f ieeepk1a.T NORMAL/.";	mv -f ieeepk1a.T NORMAL/. ; fi
	mv -f ieeepk1a.o NORMAL/ieeepk1a.o

NORMAL/mppmfl1a.o: mppmfl1a.f
	$(FORT) $(FORTOPTS)  -c mppmfl1a.f
	@if test -f mppmfl1a.T;then echo "	mv -f mppmfl1a.T NORMAL/.";	mv -f mppmfl1a.T NORMAL/. ; fi
	mv -f mppmfl1a.o NORMAL/mppmfl1a.o

NORMAL/pbfin1a.o: pbfin1a.f
	$(FORT) $(FORTOPTS)  -c pbfin1a.f
	@if test -f pbfin1a.T;then echo "	mv -f pbfin1a.T NORMAL/.";	mv -f pbfin1a.T NORMAL/. ; fi
	mv -f pbfin1a.o NORMAL/pbfin1a.o

NORMAL/pbfout1a.o: pbfout1a.f
	$(FORT) $(FORTOPTS)  -c pbfout1a.f
	@if test -f pbfout1a.T;then echo "	mv -f pbfout1a.T NORMAL/.";	mv -f pbfout1a.T NORMAL/. ; fi
	mv -f pbfout1a.o NORMAL/pbfout1a.o

NORMAL/pclose1a.o: pclose1a.f
	$(FORT) $(FORTOPTS)  -c pclose1a.f
	@if test -f pclose1a.T;then echo "	mv -f pclose1a.T NORMAL/.";	mv -f pclose1a.T NORMAL/. ; fi
	mv -f pclose1a.o NORMAL/pclose1a.o

NORMAL/popen1a.o: popen1a.f
	$(FORT) $(FORTOPTS)  -c popen1a.f
	@if test -f popen1a.T;then echo "	mv -f popen1a.T NORMAL/.";	mv -f popen1a.T NORMAL/. ; fi
	mv -f popen1a.o NORMAL/popen1a.o

NORMAL/pstpos1a.o: pstpos1a.f
	$(FORT) $(FORTOPTS)  -c pstpos1a.f
	@if test -f pstpos1a.T;then echo "	mv -f pstpos1a.T NORMAL/.";	mv -f pstpos1a.T NORMAL/. ; fi
	mv -f pstpos1a.o NORMAL/pstpos1a.o

NORMAL/rdlsm1a.o: rdlsm1a.f
	$(FORT) $(FORTOPTS)  -c rdlsm1a.f
	@if test -f rdlsm1a.T;then echo "	mv -f rdlsm1a.T NORMAL/.";	mv -f rdlsm1a.T NORMAL/. ; fi
	mv -f rdlsm1a.o NORMAL/rdlsm1a.o

NORMAL/rdmult1a.o: rdmult1a.f
	$(FORT) $(FORTOPTS)  -c rdmult1a.f
	@if test -f rdmult1a.T;then echo "	mv -f rdmult1a.T NORMAL/.";	mv -f rdmult1a.T NORMAL/. ; fi
	mv -f rdmult1a.o NORMAL/rdmult1a.o

NORMAL/scalbc1a.o: scalbc1a.f
	$(FORT) $(FORTOPTS)  -c scalbc1a.f
	@if test -f scalbc1a.T;then echo "	mv -f scalbc1a.T NORMAL/.";	mv -f scalbc1a.T NORMAL/. ; fi
	mv -f scalbc1a.o NORMAL/scalbc1a.o

NORMAL/scolbc1a.o: scolbc1a.f
	$(FORT) $(FORTOPTS)  -c scolbc1a.f
	@if test -f scolbc1a.T;then echo "	mv -f scolbc1a.T NORMAL/.";	mv -f scolbc1a.T NORMAL/. ; fi
	mv -f scolbc1a.o NORMAL/scolbc1a.o

NORMAL/sctfld1a.o: sctfld1a.f
	$(FORT) $(FORTOPTS)  -c sctfld1a.f
	@if test -f sctfld1a.T;then echo "	mv -f sctfld1a.T NORMAL/.";	mv -f sctfld1a.T NORMAL/. ; fi
	mv -f sctfld1a.o NORMAL/sctfld1a.o

NORMAL/sctzon1a.o: sctzon1a.f
	$(FORT) $(FORTOPTS)  -c sctzon1a.f
	@if test -f sctzon1a.T;then echo "	mv -f sctzon1a.T NORMAL/.";	mv -f sctzon1a.T NORMAL/. ; fi
	mv -f sctzon1a.o NORMAL/sctzon1a.o

NORMAL/setneb1a.o: setneb1a.f
	$(FORT) $(FORTOPTS)  -c setneb1a.f
	@if test -f setneb1a.T;then echo "	mv -f setneb1a.T NORMAL/.";	mv -f setneb1a.T NORMAL/. ; fi
	mv -f setneb1a.o NORMAL/setneb1a.o

NORMAL/stsfld1a.o: stsfld1a.f
	$(FORT) $(FORTOPTS)  -c stsfld1a.f
	@if test -f stsfld1a.T;then echo "	mv -f stsfld1a.T NORMAL/.";	mv -f stsfld1a.T NORMAL/. ; fi
	mv -f stsfld1a.o NORMAL/stsfld1a.o

NORMAL/stside1a.o: stside1a.f
	$(FORT) $(FORTOPTS)  -c stside1a.f
	@if test -f stside1a.T;then echo "	mv -f stside1a.T NORMAL/.";	mv -f stside1a.T NORMAL/. ; fi
	mv -f stside1a.o NORMAL/stside1a.o

NORMAL/swpbnd1a.o: swpbnd1a.f
	$(FORT) $(FORTOPTS)  -c swpbnd1a.f
	@if test -f swpbnd1a.T;then echo "	mv -f swpbnd1a.T NORMAL/.";	mv -f swpbnd1a.T NORMAL/. ; fi
	mv -f swpbnd1a.o NORMAL/swpbnd1a.o

NORMAL/wtmult1a.o: wtmult1a.f
	$(FORT) $(FORTOPTS)  -c wtmult1a.f
	@if test -f wtmult1a.T;then echo "	mv -f wtmult1a.T NORMAL/.";	mv -f wtmult1a.T NORMAL/. ; fi
	mv -f wtmult1a.o NORMAL/wtmult1a.o

