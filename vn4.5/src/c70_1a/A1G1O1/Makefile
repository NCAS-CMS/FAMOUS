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
NORMAL/abcalc1.o NORMAL/addrchk1.o NORMAL/blkdata.o NORMAL/ctoi1.o\
NORMAL/day2cha1.o NORMAL/delhist1.o NORMAL/diagdes1.o NORMAL/ereport1.o\
NORMAL/errrep1.o NORMAL/errtrce1.o NORMAL/exitpro1.o NORMAL/exppx1.o\
NORMAL/f_type1.o NORMAL/intfnfl1.o NORMAL/intfout1.o NORMAL/levsrt1.o\
NORMAL/mergfld1.o NORMAL/readstm1.o NORMAL/readwrit.o NORMAL/setperl1.o\
NORMAL/umsetup1.o NORMAL/umsubin1.o

all: $(OBJS)


NORMAL/abcalc1.o: abcalc1.f
	$(FORT) $(FORTOPTS)  -c abcalc1.f
	@if test -f abcalc1.T;then echo "	mv -f abcalc1.T NORMAL/.";	mv -f abcalc1.T NORMAL/. ; fi
	mv -f abcalc1.o NORMAL/abcalc1.o

NORMAL/addrchk1.o: addrchk1.f
	$(FORT) $(FORTOPTS)  -c addrchk1.f
	@if test -f addrchk1.T;then echo "	mv -f addrchk1.T NORMAL/.";	mv -f addrchk1.T NORMAL/. ; fi
	mv -f addrchk1.o NORMAL/addrchk1.o

NORMAL/blkdata.o: blkdata.f
	$(FORT) $(FORTOPTS)  -c blkdata.f
	@if test -f blkdata.T;then echo "	mv -f blkdata.T NORMAL/.";	mv -f blkdata.T NORMAL/. ; fi
	mv -f blkdata.o NORMAL/blkdata.o

NORMAL/ctoi1.o: ctoi1.f
	$(FORT) $(FORTOPTS)  -c ctoi1.f
	@if test -f ctoi1.T;then echo "	mv -f ctoi1.T NORMAL/.";	mv -f ctoi1.T NORMAL/. ; fi
	mv -f ctoi1.o NORMAL/ctoi1.o

NORMAL/day2cha1.o: day2cha1.f
	$(FORT) $(FORTOPTS)  -c day2cha1.f
	@if test -f day2cha1.T;then echo "	mv -f day2cha1.T NORMAL/.";	mv -f day2cha1.T NORMAL/. ; fi
	mv -f day2cha1.o NORMAL/day2cha1.o

NORMAL/delhist1.o: delhist1.f
	$(FORT) $(FORTOPTS)  -c delhist1.f
	@if test -f delhist1.T;then echo "	mv -f delhist1.T NORMAL/.";	mv -f delhist1.T NORMAL/. ; fi
	mv -f delhist1.o NORMAL/delhist1.o

NORMAL/diagdes1.o: diagdes1.f
	$(FORT) $(FORTOPTS)  -c diagdes1.f
	@if test -f diagdes1.T;then echo "	mv -f diagdes1.T NORMAL/.";	mv -f diagdes1.T NORMAL/. ; fi
	mv -f diagdes1.o NORMAL/diagdes1.o

NORMAL/ereport1.o: ereport1.f
	$(FORT) $(FORTOPTS)  -c ereport1.f
	@if test -f ereport1.T;then echo "	mv -f ereport1.T NORMAL/.";	mv -f ereport1.T NORMAL/. ; fi
	mv -f ereport1.o NORMAL/ereport1.o

NORMAL/errrep1.o: errrep1.f
	$(FORT) $(FORTOPTS)  -c errrep1.f
	@if test -f errrep1.T;then echo "	mv -f errrep1.T NORMAL/.";	mv -f errrep1.T NORMAL/. ; fi
	mv -f errrep1.o NORMAL/errrep1.o

NORMAL/errtrce1.o: errtrce1.f
	$(FORT) $(FORTOPTS)  -c errtrce1.f
	@if test -f errtrce1.T;then echo "	mv -f errtrce1.T NORMAL/.";	mv -f errtrce1.T NORMAL/. ; fi
	mv -f errtrce1.o NORMAL/errtrce1.o

NORMAL/exitpro1.o: exitpro1.f
	$(FORT) $(FORTOPTS)  -c exitpro1.f
	@if test -f exitpro1.T;then echo "	mv -f exitpro1.T NORMAL/.";	mv -f exitpro1.T NORMAL/. ; fi
	mv -f exitpro1.o NORMAL/exitpro1.o

NORMAL/exppx1.o: exppx1.f
	$(FORT) $(FORTOPTS)  -c exppx1.f
	@if test -f exppx1.T;then echo "	mv -f exppx1.T NORMAL/.";	mv -f exppx1.T NORMAL/. ; fi
	mv -f exppx1.o NORMAL/exppx1.o

NORMAL/f_type1.o: f_type1.f
	$(FORT) $(FORTOPTS)  -c f_type1.f
	@if test -f f_type1.T;then echo "	mv -f f_type1.T NORMAL/.";	mv -f f_type1.T NORMAL/. ; fi
	mv -f f_type1.o NORMAL/f_type1.o

NORMAL/intfnfl1.o: intfnfl1.f
	$(FORT) $(FORTOPTS)  -c intfnfl1.f
	@if test -f intfnfl1.T;then echo "	mv -f intfnfl1.T NORMAL/.";	mv -f intfnfl1.T NORMAL/. ; fi
	mv -f intfnfl1.o NORMAL/intfnfl1.o

NORMAL/intfout1.o: intfout1.f
	$(FORT) $(FORTOPTS)  -c intfout1.f
	@if test -f intfout1.T;then echo "	mv -f intfout1.T NORMAL/.";	mv -f intfout1.T NORMAL/. ; fi
	mv -f intfout1.o NORMAL/intfout1.o

NORMAL/levsrt1.o: levsrt1.f
	$(FORT) $(FORTOPTS)  -c levsrt1.f
	@if test -f levsrt1.T;then echo "	mv -f levsrt1.T NORMAL/.";	mv -f levsrt1.T NORMAL/. ; fi
	mv -f levsrt1.o NORMAL/levsrt1.o

NORMAL/mergfld1.o: mergfld1.f
	$(FORT) $(FORTOPTS)  -c mergfld1.f
	@if test -f mergfld1.T;then echo "	mv -f mergfld1.T NORMAL/.";	mv -f mergfld1.T NORMAL/. ; fi
	mv -f mergfld1.o NORMAL/mergfld1.o

NORMAL/readstm1.o: readstm1.f
	$(FORT) $(FORTOPTS)  -c readstm1.f
	@if test -f readstm1.T;then echo "	mv -f readstm1.T NORMAL/.";	mv -f readstm1.T NORMAL/. ; fi
	mv -f readstm1.o NORMAL/readstm1.o

NORMAL/readwrit.o: readwrit.f
	$(FORT) $(FORTOPTS)  -c readwrit.f
	@if test -f readwrit.T;then echo "	mv -f readwrit.T NORMAL/.";	mv -f readwrit.T NORMAL/. ; fi
	mv -f readwrit.o NORMAL/readwrit.o

NORMAL/setperl1.o: setperl1.f
	$(FORT) $(FORTOPTS)  -c setperl1.f
	@if test -f setperl1.T;then echo "	mv -f setperl1.T NORMAL/.";	mv -f setperl1.T NORMAL/. ; fi
	mv -f setperl1.o NORMAL/setperl1.o

NORMAL/umsetup1.o: umsetup1.f
	$(FORT) $(FORTOPTS)  -c umsetup1.f
	@if test -f umsetup1.T;then echo "	mv -f umsetup1.T NORMAL/.";	mv -f umsetup1.T NORMAL/. ; fi
	mv -f umsetup1.o NORMAL/umsetup1.o

NORMAL/umsubin1.o: umsubin1.f
	$(FORT) $(FORTOPTS)  -c umsubin1.f
	@if test -f umsubin1.T;then echo "	mv -f umsubin1.T NORMAL/.";	mv -f umsubin1.T NORMAL/. ; fi
	mv -f umsubin1.o NORMAL/umsubin1.o

