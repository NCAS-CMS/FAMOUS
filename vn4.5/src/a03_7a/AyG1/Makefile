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
NORMAL/bdylyr7a.o NORMAL/bl_ic7a.o NORMAL/bl_lsp1a.o NORMAL/bouytq5b.o\
NORMAL/btqint5b.o NORMAL/excoef7a.o NORMAL/exfxtq5b.o NORMAL/exfxuv5b.o\
NORMAL/fcdch7a.o NORMAL/htcond5a.o NORMAL/imcltq7a.o NORMAL/imcluv5b.o\
NORMAL/impmix2c.o NORMAL/kmkh7a.o NORMAL/leaf7a.o NORMAL/llit7a.o\
NORMAL/microb7a.o NORMAL/physio7a.o NORMAL/radmoses.o NORMAL/raero7a.o\
NORMAL/rootfr7a.o NORMAL/screen7a.o NORMAL/sfevap7a.o NORMAL/sfexch7a.o\
NORMAL/sfflux7a.o NORMAL/sflint7a.o NORMAL/sfmelt7a.o NORMAL/sforog7a.o\
NORMAL/sfrest7a.o NORMAL/sfrib7a.o NORMAL/sfstom7a.o NORMAL/siceht5b.o\
NORMAL/smcext7a.o NORMAL/stdev7a.o NORMAL/trmix2c.o NORMAL/trsrce1a.o\
NORMAL/vshrz7a.o NORMAL/z1a.o

all: $(OBJS)


NORMAL/bdylyr7a.o: bdylyr7a.f
	$(FORT) $(FORTOPTS)  -c bdylyr7a.f
	@if test -f bdylyr7a.T;then echo "	mv -f bdylyr7a.T NORMAL/.";	mv -f bdylyr7a.T NORMAL/. ; fi
	mv -f bdylyr7a.o NORMAL/bdylyr7a.o

NORMAL/bl_ic7a.o: bl_ic7a.f
	$(FORT) $(FORTOPTS)  -c bl_ic7a.f
	@if test -f bl_ic7a.T;then echo "	mv -f bl_ic7a.T NORMAL/.";	mv -f bl_ic7a.T NORMAL/. ; fi
	mv -f bl_ic7a.o NORMAL/bl_ic7a.o

NORMAL/bl_lsp1a.o: bl_lsp1a.f
	$(FORT) $(FORTOPTS)  -c bl_lsp1a.f
	@if test -f bl_lsp1a.T;then echo "	mv -f bl_lsp1a.T NORMAL/.";	mv -f bl_lsp1a.T NORMAL/. ; fi
	mv -f bl_lsp1a.o NORMAL/bl_lsp1a.o

NORMAL/bouytq5b.o: bouytq5b.f
	$(FORT) $(FORTOPTS)  -c bouytq5b.f
	@if test -f bouytq5b.T;then echo "	mv -f bouytq5b.T NORMAL/.";	mv -f bouytq5b.T NORMAL/. ; fi
	mv -f bouytq5b.o NORMAL/bouytq5b.o

NORMAL/btqint5b.o: btqint5b.f
	$(FORT) $(FORTOPTS)  -c btqint5b.f
	@if test -f btqint5b.T;then echo "	mv -f btqint5b.T NORMAL/.";	mv -f btqint5b.T NORMAL/. ; fi
	mv -f btqint5b.o NORMAL/btqint5b.o

NORMAL/excoef7a.o: excoef7a.f
	$(FORT) $(FORTOPTS)  -c excoef7a.f
	@if test -f excoef7a.T;then echo "	mv -f excoef7a.T NORMAL/.";	mv -f excoef7a.T NORMAL/. ; fi
	mv -f excoef7a.o NORMAL/excoef7a.o

NORMAL/exfxtq5b.o: exfxtq5b.f
	$(FORT) $(FORTOPTS)  -c exfxtq5b.f
	@if test -f exfxtq5b.T;then echo "	mv -f exfxtq5b.T NORMAL/.";	mv -f exfxtq5b.T NORMAL/. ; fi
	mv -f exfxtq5b.o NORMAL/exfxtq5b.o

NORMAL/exfxuv5b.o: exfxuv5b.f
	$(FORT) $(FORTOPTS)  -c exfxuv5b.f
	@if test -f exfxuv5b.T;then echo "	mv -f exfxuv5b.T NORMAL/.";	mv -f exfxuv5b.T NORMAL/. ; fi
	mv -f exfxuv5b.o NORMAL/exfxuv5b.o

NORMAL/fcdch7a.o: fcdch7a.f
	$(FORT) $(FORTOPTS)  -c fcdch7a.f
	@if test -f fcdch7a.T;then echo "	mv -f fcdch7a.T NORMAL/.";	mv -f fcdch7a.T NORMAL/. ; fi
	mv -f fcdch7a.o NORMAL/fcdch7a.o

NORMAL/htcond5a.o: htcond5a.f
	$(FORT) $(FORTOPTS)  -c htcond5a.f
	@if test -f htcond5a.T;then echo "	mv -f htcond5a.T NORMAL/.";	mv -f htcond5a.T NORMAL/. ; fi
	mv -f htcond5a.o NORMAL/htcond5a.o

NORMAL/imcltq7a.o: imcltq7a.f
	$(FORT) $(FORTOPTS)  -c imcltq7a.f
	@if test -f imcltq7a.T;then echo "	mv -f imcltq7a.T NORMAL/.";	mv -f imcltq7a.T NORMAL/. ; fi
	mv -f imcltq7a.o NORMAL/imcltq7a.o

NORMAL/imcluv5b.o: imcluv5b.f
	$(FORT) $(FORTOPTS)  -c imcluv5b.f
	@if test -f imcluv5b.T;then echo "	mv -f imcluv5b.T NORMAL/.";	mv -f imcluv5b.T NORMAL/. ; fi
	mv -f imcluv5b.o NORMAL/imcluv5b.o

NORMAL/impmix2c.o: impmix2c.f
	$(FORT) $(FORTOPTS)  -c impmix2c.f
	@if test -f impmix2c.T;then echo "	mv -f impmix2c.T NORMAL/.";	mv -f impmix2c.T NORMAL/. ; fi
	mv -f impmix2c.o NORMAL/impmix2c.o

NORMAL/kmkh7a.o: kmkh7a.f
	$(FORT) $(FORTOPTS)  -c kmkh7a.f
	@if test -f kmkh7a.T;then echo "	mv -f kmkh7a.T NORMAL/.";	mv -f kmkh7a.T NORMAL/. ; fi
	mv -f kmkh7a.o NORMAL/kmkh7a.o

NORMAL/leaf7a.o: leaf7a.f
	$(FORT) $(FORTOPTS)  -c leaf7a.f
	@if test -f leaf7a.T;then echo "	mv -f leaf7a.T NORMAL/.";	mv -f leaf7a.T NORMAL/. ; fi
	mv -f leaf7a.o NORMAL/leaf7a.o

NORMAL/llit7a.o: llit7a.f
	$(FORT) $(FORTOPTS)  -c llit7a.f
	@if test -f llit7a.T;then echo "	mv -f llit7a.T NORMAL/.";	mv -f llit7a.T NORMAL/. ; fi
	mv -f llit7a.o NORMAL/llit7a.o

NORMAL/microb7a.o: microb7a.f
	$(FORT) $(FORTOPTS)  -c microb7a.f
	@if test -f microb7a.T;then echo "	mv -f microb7a.T NORMAL/.";	mv -f microb7a.T NORMAL/. ; fi
	mv -f microb7a.o NORMAL/microb7a.o

NORMAL/physio7a.o: physio7a.f
	$(FORT) $(FORTOPTS)  -c physio7a.f
	@if test -f physio7a.T;then echo "	mv -f physio7a.T NORMAL/.";	mv -f physio7a.T NORMAL/. ; fi
	mv -f physio7a.o NORMAL/physio7a.o

NORMAL/radmoses.o: radmoses.f
	$(FORT) $(FORTOPTS)  -c radmoses.f
	@if test -f radmoses.T;then echo "	mv -f radmoses.T NORMAL/.";	mv -f radmoses.T NORMAL/. ; fi
	mv -f radmoses.o NORMAL/radmoses.o

NORMAL/raero7a.o: raero7a.f
	$(FORT) $(FORTOPTS)  -c raero7a.f
	@if test -f raero7a.T;then echo "	mv -f raero7a.T NORMAL/.";	mv -f raero7a.T NORMAL/. ; fi
	mv -f raero7a.o NORMAL/raero7a.o

NORMAL/rootfr7a.o: rootfr7a.f
	$(FORT) $(FORTOPTS)  -c rootfr7a.f
	@if test -f rootfr7a.T;then echo "	mv -f rootfr7a.T NORMAL/.";	mv -f rootfr7a.T NORMAL/. ; fi
	mv -f rootfr7a.o NORMAL/rootfr7a.o

NORMAL/screen7a.o: screen7a.f
	$(FORT) $(FORTOPTS)  -c screen7a.f
	@if test -f screen7a.T;then echo "	mv -f screen7a.T NORMAL/.";	mv -f screen7a.T NORMAL/. ; fi
	mv -f screen7a.o NORMAL/screen7a.o

NORMAL/sfevap7a.o: sfevap7a.f
	$(FORT) $(FORTOPTS)  -c sfevap7a.f
	@if test -f sfevap7a.T;then echo "	mv -f sfevap7a.T NORMAL/.";	mv -f sfevap7a.T NORMAL/. ; fi
	mv -f sfevap7a.o NORMAL/sfevap7a.o

NORMAL/sfexch7a.o: sfexch7a.f
	$(FORT) $(FORTOPTS)  -c sfexch7a.f
	@if test -f sfexch7a.T;then echo "	mv -f sfexch7a.T NORMAL/.";	mv -f sfexch7a.T NORMAL/. ; fi
	mv -f sfexch7a.o NORMAL/sfexch7a.o

NORMAL/sfflux7a.o: sfflux7a.f
	$(FORT) $(FORTOPTS)  -c sfflux7a.f
	@if test -f sfflux7a.T;then echo "	mv -f sfflux7a.T NORMAL/.";	mv -f sfflux7a.T NORMAL/. ; fi
	mv -f sfflux7a.o NORMAL/sfflux7a.o

NORMAL/sflint7a.o: sflint7a.f
	$(FORT) $(FORTOPTS)  -c sflint7a.f
	@if test -f sflint7a.T;then echo "	mv -f sflint7a.T NORMAL/.";	mv -f sflint7a.T NORMAL/. ; fi
	mv -f sflint7a.o NORMAL/sflint7a.o

NORMAL/sfmelt7a.o: sfmelt7a.f
	$(FORT) $(FORTOPTS)  -c sfmelt7a.f
	@if test -f sfmelt7a.T;then echo "	mv -f sfmelt7a.T NORMAL/.";	mv -f sfmelt7a.T NORMAL/. ; fi
	mv -f sfmelt7a.o NORMAL/sfmelt7a.o

NORMAL/sforog7a.o: sforog7a.f
	$(FORT) $(FORTOPTS)  -c sforog7a.f
	@if test -f sforog7a.T;then echo "	mv -f sforog7a.T NORMAL/.";	mv -f sforog7a.T NORMAL/. ; fi
	mv -f sforog7a.o NORMAL/sforog7a.o

NORMAL/sfrest7a.o: sfrest7a.f
	$(FORT) $(FORTOPTS)  -c sfrest7a.f
	@if test -f sfrest7a.T;then echo "	mv -f sfrest7a.T NORMAL/.";	mv -f sfrest7a.T NORMAL/. ; fi
	mv -f sfrest7a.o NORMAL/sfrest7a.o

NORMAL/sfrib7a.o: sfrib7a.f
	$(FORT) $(FORTOPTS)  -c sfrib7a.f
	@if test -f sfrib7a.T;then echo "	mv -f sfrib7a.T NORMAL/.";	mv -f sfrib7a.T NORMAL/. ; fi
	mv -f sfrib7a.o NORMAL/sfrib7a.o

NORMAL/sfstom7a.o: sfstom7a.f
	$(FORT) $(FORTOPTS)  -c sfstom7a.f
	@if test -f sfstom7a.T;then echo "	mv -f sfstom7a.T NORMAL/.";	mv -f sfstom7a.T NORMAL/. ; fi
	mv -f sfstom7a.o NORMAL/sfstom7a.o

NORMAL/siceht5b.o: siceht5b.f
	$(FORT) $(FORTOPTS)  -c siceht5b.f
	@if test -f siceht5b.T;then echo "	mv -f siceht5b.T NORMAL/.";	mv -f siceht5b.T NORMAL/. ; fi
	mv -f siceht5b.o NORMAL/siceht5b.o

NORMAL/smcext7a.o: smcext7a.f
	$(FORT) $(FORTOPTS)  -c smcext7a.f
	@if test -f smcext7a.T;then echo "	mv -f smcext7a.T NORMAL/.";	mv -f smcext7a.T NORMAL/. ; fi
	mv -f smcext7a.o NORMAL/smcext7a.o

NORMAL/stdev7a.o: stdev7a.f
	$(FORT) $(FORTOPTS)  -c stdev7a.f
	@if test -f stdev7a.T;then echo "	mv -f stdev7a.T NORMAL/.";	mv -f stdev7a.T NORMAL/. ; fi
	mv -f stdev7a.o NORMAL/stdev7a.o

NORMAL/trmix2c.o: trmix2c.f
	$(FORT) $(FORTOPTS)  -c trmix2c.f
	@if test -f trmix2c.T;then echo "	mv -f trmix2c.T NORMAL/.";	mv -f trmix2c.T NORMAL/. ; fi
	mv -f trmix2c.o NORMAL/trmix2c.o

NORMAL/trsrce1a.o: trsrce1a.f
	$(FORT) $(FORTOPTS)  -c trsrce1a.f
	@if test -f trsrce1a.T;then echo "	mv -f trsrce1a.T NORMAL/.";	mv -f trsrce1a.T NORMAL/. ; fi
	mv -f trsrce1a.o NORMAL/trsrce1a.o

NORMAL/vshrz7a.o: vshrz7a.f
	$(FORT) $(FORTOPTS)  -c vshrz7a.f
	@if test -f vshrz7a.T;then echo "	mv -f vshrz7a.T NORMAL/.";	mv -f vshrz7a.T NORMAL/. ; fi
	mv -f vshrz7a.o NORMAL/vshrz7a.o

NORMAL/z1a.o: z1a.f
	$(FORT) $(FORTOPTS)  -c z1a.f
	@if test -f z1a.T;then echo "	mv -f z1a.T NORMAL/.";	mv -f z1a.T NORMAL/. ; fi
	mv -f z1a.o NORMAL/z1a.o

