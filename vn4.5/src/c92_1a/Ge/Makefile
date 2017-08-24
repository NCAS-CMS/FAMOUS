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
NORMAL/boxsum1.o NORMAL/coasaj1a.o NORMAL/eqtoll1a.o NORMAL/fractm1a.o\
NORMAL/hint_aw1.o NORMAL/hint_bl1.o NORMAL/hintco1a.o NORMAL/hintctl1.o\
NORMAL/horcub1a.o NORMAL/horint1a.o NORMAL/horlin1a.o NORMAL/hormon1a.o\
NORMAL/horqui1a.o NORMAL/intfca1a.o NORMAL/lltoeq1a.o NORMAL/nearpt1a.o\
NORMAL/pmsl1a.o NORMAL/sprals1a.o NORMAL/tint1a.o NORMAL/tintc1a.o\
NORMAL/trop1a.o NORMAL/vercub1a.o NORMAL/verint1a.o NORMAL/verlin1a.o\
NORMAL/verqui1a.o NORMAL/vint1a.o NORMAL/vintt1a.o NORMAL/vintth1a.o\
NORMAL/vinttp1a.o NORMAL/vintz1a.o NORMAL/vintzh1a.o NORMAL/wcoeff1a.o\
NORMAL/weqtol1a.o NORMAL/wlltoe1a.o

all: $(OBJS)


NORMAL/boxsum1.o: boxsum1.f
	$(FORT) $(FORTOPTS)  -c boxsum1.f
	@if test -f boxsum1.T;then echo "	mv -f boxsum1.T NORMAL/.";	mv -f boxsum1.T NORMAL/. ; fi
	mv -f boxsum1.o NORMAL/boxsum1.o

NORMAL/coasaj1a.o: coasaj1a.f
	$(FORT) $(FORTOPTS)  -c coasaj1a.f
	@if test -f coasaj1a.T;then echo "	mv -f coasaj1a.T NORMAL/.";	mv -f coasaj1a.T NORMAL/. ; fi
	mv -f coasaj1a.o NORMAL/coasaj1a.o

NORMAL/eqtoll1a.o: eqtoll1a.f
	$(FORT) $(FORTOPTS)  -c eqtoll1a.f
	@if test -f eqtoll1a.T;then echo "	mv -f eqtoll1a.T NORMAL/.";	mv -f eqtoll1a.T NORMAL/. ; fi
	mv -f eqtoll1a.o NORMAL/eqtoll1a.o

NORMAL/fractm1a.o: fractm1a.f
	$(FORT) $(FORTOPTS)  -c fractm1a.f
	@if test -f fractm1a.T;then echo "	mv -f fractm1a.T NORMAL/.";	mv -f fractm1a.T NORMAL/. ; fi
	mv -f fractm1a.o NORMAL/fractm1a.o

NORMAL/hint_aw1.o: hint_aw1.f
	$(FORT) $(FORTOPTS)  -c hint_aw1.f
	@if test -f hint_aw1.T;then echo "	mv -f hint_aw1.T NORMAL/.";	mv -f hint_aw1.T NORMAL/. ; fi
	mv -f hint_aw1.o NORMAL/hint_aw1.o

NORMAL/hint_bl1.o: hint_bl1.f
	$(FORT) $(FORTOPTS)  -c hint_bl1.f
	@if test -f hint_bl1.T;then echo "	mv -f hint_bl1.T NORMAL/.";	mv -f hint_bl1.T NORMAL/. ; fi
	mv -f hint_bl1.o NORMAL/hint_bl1.o

NORMAL/hintco1a.o: hintco1a.f
	$(FORT) $(FORTOPTS)  -c hintco1a.f
	@if test -f hintco1a.T;then echo "	mv -f hintco1a.T NORMAL/.";	mv -f hintco1a.T NORMAL/. ; fi
	mv -f hintco1a.o NORMAL/hintco1a.o

NORMAL/hintctl1.o: hintctl1.f
	$(FORT) $(FORTOPTS)  -c hintctl1.f
	@if test -f hintctl1.T;then echo "	mv -f hintctl1.T NORMAL/.";	mv -f hintctl1.T NORMAL/. ; fi
	mv -f hintctl1.o NORMAL/hintctl1.o

NORMAL/horcub1a.o: horcub1a.f
	$(FORT) $(FORTOPTS)  -c horcub1a.f
	@if test -f horcub1a.T;then echo "	mv -f horcub1a.T NORMAL/.";	mv -f horcub1a.T NORMAL/. ; fi
	mv -f horcub1a.o NORMAL/horcub1a.o

NORMAL/horint1a.o: horint1a.f
	$(FORT) $(FORTOPTS)  -c horint1a.f
	@if test -f horint1a.T;then echo "	mv -f horint1a.T NORMAL/.";	mv -f horint1a.T NORMAL/. ; fi
	mv -f horint1a.o NORMAL/horint1a.o

NORMAL/horlin1a.o: horlin1a.f
	$(FORT) $(FORTOPTS)  -c horlin1a.f
	@if test -f horlin1a.T;then echo "	mv -f horlin1a.T NORMAL/.";	mv -f horlin1a.T NORMAL/. ; fi
	mv -f horlin1a.o NORMAL/horlin1a.o

NORMAL/hormon1a.o: hormon1a.f
	$(FORT) $(FORTOPTS)  -c hormon1a.f
	@if test -f hormon1a.T;then echo "	mv -f hormon1a.T NORMAL/.";	mv -f hormon1a.T NORMAL/. ; fi
	mv -f hormon1a.o NORMAL/hormon1a.o

NORMAL/horqui1a.o: horqui1a.f
	$(FORT) $(FORTOPTS)  -c horqui1a.f
	@if test -f horqui1a.T;then echo "	mv -f horqui1a.T NORMAL/.";	mv -f horqui1a.T NORMAL/. ; fi
	mv -f horqui1a.o NORMAL/horqui1a.o

NORMAL/intfca1a.o: intfca1a.f
	$(FORT) $(FORTOPTS)  -c intfca1a.f
	@if test -f intfca1a.T;then echo "	mv -f intfca1a.T NORMAL/.";	mv -f intfca1a.T NORMAL/. ; fi
	mv -f intfca1a.o NORMAL/intfca1a.o

NORMAL/lltoeq1a.o: lltoeq1a.f
	$(FORT) $(FORTOPTS)  -c lltoeq1a.f
	@if test -f lltoeq1a.T;then echo "	mv -f lltoeq1a.T NORMAL/.";	mv -f lltoeq1a.T NORMAL/. ; fi
	mv -f lltoeq1a.o NORMAL/lltoeq1a.o

NORMAL/nearpt1a.o: nearpt1a.f
	$(FORT) $(FORTOPTS)  -c nearpt1a.f
	@if test -f nearpt1a.T;then echo "	mv -f nearpt1a.T NORMAL/.";	mv -f nearpt1a.T NORMAL/. ; fi
	mv -f nearpt1a.o NORMAL/nearpt1a.o

NORMAL/pmsl1a.o: pmsl1a.f
	$(FORT) $(FORTOPTS)  -c pmsl1a.f
	@if test -f pmsl1a.T;then echo "	mv -f pmsl1a.T NORMAL/.";	mv -f pmsl1a.T NORMAL/. ; fi
	mv -f pmsl1a.o NORMAL/pmsl1a.o

NORMAL/sprals1a.o: sprals1a.f
	$(FORT) $(FORTOPTS)  -c sprals1a.f
	@if test -f sprals1a.T;then echo "	mv -f sprals1a.T NORMAL/.";	mv -f sprals1a.T NORMAL/. ; fi
	mv -f sprals1a.o NORMAL/sprals1a.o

NORMAL/tint1a.o: tint1a.f
	$(FORT) $(FORTOPTS)  -c tint1a.f
	@if test -f tint1a.T;then echo "	mv -f tint1a.T NORMAL/.";	mv -f tint1a.T NORMAL/. ; fi
	mv -f tint1a.o NORMAL/tint1a.o

NORMAL/tintc1a.o: tintc1a.f
	$(FORT) $(FORTOPTS)  -c tintc1a.f
	@if test -f tintc1a.T;then echo "	mv -f tintc1a.T NORMAL/.";	mv -f tintc1a.T NORMAL/. ; fi
	mv -f tintc1a.o NORMAL/tintc1a.o

NORMAL/trop1a.o: trop1a.f
	$(FORT) $(FORTOPTS)  -c trop1a.f
	@if test -f trop1a.T;then echo "	mv -f trop1a.T NORMAL/.";	mv -f trop1a.T NORMAL/. ; fi
	mv -f trop1a.o NORMAL/trop1a.o

NORMAL/vercub1a.o: vercub1a.f
	$(FORT) $(FORTOPTS)  -c vercub1a.f
	@if test -f vercub1a.T;then echo "	mv -f vercub1a.T NORMAL/.";	mv -f vercub1a.T NORMAL/. ; fi
	mv -f vercub1a.o NORMAL/vercub1a.o

NORMAL/verint1a.o: verint1a.f
	$(FORT) $(FORTOPTS)  -c verint1a.f
	@if test -f verint1a.T;then echo "	mv -f verint1a.T NORMAL/.";	mv -f verint1a.T NORMAL/. ; fi
	mv -f verint1a.o NORMAL/verint1a.o

NORMAL/verlin1a.o: verlin1a.f
	$(FORT) $(FORTOPTS)  -c verlin1a.f
	@if test -f verlin1a.T;then echo "	mv -f verlin1a.T NORMAL/.";	mv -f verlin1a.T NORMAL/. ; fi
	mv -f verlin1a.o NORMAL/verlin1a.o

NORMAL/verqui1a.o: verqui1a.f
	$(FORT) $(FORTOPTS)  -c verqui1a.f
	@if test -f verqui1a.T;then echo "	mv -f verqui1a.T NORMAL/.";	mv -f verqui1a.T NORMAL/. ; fi
	mv -f verqui1a.o NORMAL/verqui1a.o

NORMAL/vint1a.o: vint1a.f
	$(FORT) $(FORTOPTS)  -c vint1a.f
	@if test -f vint1a.T;then echo "	mv -f vint1a.T NORMAL/.";	mv -f vint1a.T NORMAL/. ; fi
	mv -f vint1a.o NORMAL/vint1a.o

NORMAL/vintt1a.o: vintt1a.f
	$(FORT) $(FORTOPTS)  -c vintt1a.f
	@if test -f vintt1a.T;then echo "	mv -f vintt1a.T NORMAL/.";	mv -f vintt1a.T NORMAL/. ; fi
	mv -f vintt1a.o NORMAL/vintt1a.o

NORMAL/vintth1a.o: vintth1a.f
	$(FORT) $(FORTOPTS)  -c vintth1a.f
	@if test -f vintth1a.T;then echo "	mv -f vintth1a.T NORMAL/.";	mv -f vintth1a.T NORMAL/. ; fi
	mv -f vintth1a.o NORMAL/vintth1a.o

NORMAL/vinttp1a.o: vinttp1a.f
	$(FORT) $(FORTOPTS)  -c vinttp1a.f
	@if test -f vinttp1a.T;then echo "	mv -f vinttp1a.T NORMAL/.";	mv -f vinttp1a.T NORMAL/. ; fi
	mv -f vinttp1a.o NORMAL/vinttp1a.o

NORMAL/vintz1a.o: vintz1a.f
	$(FORT) $(FORTOPTS)  -c vintz1a.f
	@if test -f vintz1a.T;then echo "	mv -f vintz1a.T NORMAL/.";	mv -f vintz1a.T NORMAL/. ; fi
	mv -f vintz1a.o NORMAL/vintz1a.o

NORMAL/vintzh1a.o: vintzh1a.f
	$(FORT) $(FORTOPTS)  -c vintzh1a.f
	@if test -f vintzh1a.T;then echo "	mv -f vintzh1a.T NORMAL/.";	mv -f vintzh1a.T NORMAL/. ; fi
	mv -f vintzh1a.o NORMAL/vintzh1a.o

NORMAL/wcoeff1a.o: wcoeff1a.f
	$(FORT) $(FORTOPTS)  -c wcoeff1a.f
	@if test -f wcoeff1a.T;then echo "	mv -f wcoeff1a.T NORMAL/.";	mv -f wcoeff1a.T NORMAL/. ; fi
	mv -f wcoeff1a.o NORMAL/wcoeff1a.o

NORMAL/weqtol1a.o: weqtol1a.f
	$(FORT) $(FORTOPTS)  -c weqtol1a.f
	@if test -f weqtol1a.T;then echo "	mv -f weqtol1a.T NORMAL/.";	mv -f weqtol1a.T NORMAL/. ; fi
	mv -f weqtol1a.o NORMAL/weqtol1a.o

NORMAL/wlltoe1a.o: wlltoe1a.f
	$(FORT) $(FORTOPTS)  -c wlltoe1a.f
	@if test -f wlltoe1a.T;then echo "	mv -f wlltoe1a.T NORMAL/.";	mv -f wlltoe1a.T NORMAL/. ; fi
	mv -f wlltoe1a.o NORMAL/wlltoe1a.o

