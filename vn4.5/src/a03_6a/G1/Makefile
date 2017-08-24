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
NORMAL/bdylyr6a.o NORMAL/bl_ic6a.o NORMAL/bl_lsp1a.o NORMAL/bouytq6a.o\
NORMAL/btqint6a.o NORMAL/excfnl6a.o NORMAL/excoef6a.o NORMAL/exfxtq6a.o\
NORMAL/exfxuv5b.o NORMAL/fcdch6a.o NORMAL/htcond5a.o NORMAL/imcltq5b.o\
NORMAL/imcluv5b.o NORMAL/impmix2c.o NORMAL/kmkh6a.o NORMAL/kmkhz6a.o\
NORMAL/leaf5b.o NORMAL/phimh6a.o NORMAL/sfevap5b.o NORMAL/sfexch6a.o\
NORMAL/sfflux5b.o NORMAL/sflbes5b.o NORMAL/sflint6a.o NORMAL/sfmelt5b.o\
NORMAL/sfrest6a.o NORMAL/sfrib6a.o NORMAL/sfrugh6a.o NORMAL/sfstom5b.o\
NORMAL/siceht5b.o NORMAL/trmix2c.o NORMAL/trsrce1a.o NORMAL/z1a.o

all: $(OBJS)


NORMAL/bdylyr6a.o: bdylyr6a.f
	$(FORT) $(FORTOPTS)  -c bdylyr6a.f
	@if test -f bdylyr6a.T;then echo "	mv -f bdylyr6a.T NORMAL/.";	mv -f bdylyr6a.T NORMAL/. ; fi
	mv -f bdylyr6a.o NORMAL/bdylyr6a.o

NORMAL/bl_ic6a.o: bl_ic6a.f
	$(FORT) $(FORTOPTS)  -c bl_ic6a.f
	@if test -f bl_ic6a.T;then echo "	mv -f bl_ic6a.T NORMAL/.";	mv -f bl_ic6a.T NORMAL/. ; fi
	mv -f bl_ic6a.o NORMAL/bl_ic6a.o

NORMAL/bl_lsp1a.o: bl_lsp1a.f
	$(FORT) $(FORTOPTS)  -c bl_lsp1a.f
	@if test -f bl_lsp1a.T;then echo "	mv -f bl_lsp1a.T NORMAL/.";	mv -f bl_lsp1a.T NORMAL/. ; fi
	mv -f bl_lsp1a.o NORMAL/bl_lsp1a.o

NORMAL/bouytq6a.o: bouytq6a.f
	$(FORT) $(FORTOPTS)  -c bouytq6a.f
	@if test -f bouytq6a.T;then echo "	mv -f bouytq6a.T NORMAL/.";	mv -f bouytq6a.T NORMAL/. ; fi
	mv -f bouytq6a.o NORMAL/bouytq6a.o

NORMAL/btqint6a.o: btqint6a.f
	$(FORT) $(FORTOPTS)  -c btqint6a.f
	@if test -f btqint6a.T;then echo "	mv -f btqint6a.T NORMAL/.";	mv -f btqint6a.T NORMAL/. ; fi
	mv -f btqint6a.o NORMAL/btqint6a.o

NORMAL/excfnl6a.o: excfnl6a.f
	$(FORT) $(FORTOPTS)  -c excfnl6a.f
	@if test -f excfnl6a.T;then echo "	mv -f excfnl6a.T NORMAL/.";	mv -f excfnl6a.T NORMAL/. ; fi
	mv -f excfnl6a.o NORMAL/excfnl6a.o

NORMAL/excoef6a.o: excoef6a.f
	$(FORT) $(FORTOPTS)  -c excoef6a.f
	@if test -f excoef6a.T;then echo "	mv -f excoef6a.T NORMAL/.";	mv -f excoef6a.T NORMAL/. ; fi
	mv -f excoef6a.o NORMAL/excoef6a.o

NORMAL/exfxtq6a.o: exfxtq6a.f
	$(FORT) $(FORTOPTS)  -c exfxtq6a.f
	@if test -f exfxtq6a.T;then echo "	mv -f exfxtq6a.T NORMAL/.";	mv -f exfxtq6a.T NORMAL/. ; fi
	mv -f exfxtq6a.o NORMAL/exfxtq6a.o

NORMAL/exfxuv5b.o: exfxuv5b.f
	$(FORT) $(FORTOPTS)  -c exfxuv5b.f
	@if test -f exfxuv5b.T;then echo "	mv -f exfxuv5b.T NORMAL/.";	mv -f exfxuv5b.T NORMAL/. ; fi
	mv -f exfxuv5b.o NORMAL/exfxuv5b.o

NORMAL/fcdch6a.o: fcdch6a.f
	$(FORT) $(FORTOPTS)  -c fcdch6a.f
	@if test -f fcdch6a.T;then echo "	mv -f fcdch6a.T NORMAL/.";	mv -f fcdch6a.T NORMAL/. ; fi
	mv -f fcdch6a.o NORMAL/fcdch6a.o

NORMAL/htcond5a.o: htcond5a.f
	$(FORT) $(FORTOPTS)  -c htcond5a.f
	@if test -f htcond5a.T;then echo "	mv -f htcond5a.T NORMAL/.";	mv -f htcond5a.T NORMAL/. ; fi
	mv -f htcond5a.o NORMAL/htcond5a.o

NORMAL/imcltq5b.o: imcltq5b.f
	$(FORT) $(FORTOPTS)  -c imcltq5b.f
	@if test -f imcltq5b.T;then echo "	mv -f imcltq5b.T NORMAL/.";	mv -f imcltq5b.T NORMAL/. ; fi
	mv -f imcltq5b.o NORMAL/imcltq5b.o

NORMAL/imcluv5b.o: imcluv5b.f
	$(FORT) $(FORTOPTS)  -c imcluv5b.f
	@if test -f imcluv5b.T;then echo "	mv -f imcluv5b.T NORMAL/.";	mv -f imcluv5b.T NORMAL/. ; fi
	mv -f imcluv5b.o NORMAL/imcluv5b.o

NORMAL/impmix2c.o: impmix2c.f
	$(FORT) $(FORTOPTS)  -c impmix2c.f
	@if test -f impmix2c.T;then echo "	mv -f impmix2c.T NORMAL/.";	mv -f impmix2c.T NORMAL/. ; fi
	mv -f impmix2c.o NORMAL/impmix2c.o

NORMAL/kmkh6a.o: kmkh6a.f
	$(FORT) $(FORTOPTS)  -c kmkh6a.f
	@if test -f kmkh6a.T;then echo "	mv -f kmkh6a.T NORMAL/.";	mv -f kmkh6a.T NORMAL/. ; fi
	mv -f kmkh6a.o NORMAL/kmkh6a.o

NORMAL/kmkhz6a.o: kmkhz6a.f
	$(FORT) $(FORTOPTS)  -c kmkhz6a.f
	@if test -f kmkhz6a.T;then echo "	mv -f kmkhz6a.T NORMAL/.";	mv -f kmkhz6a.T NORMAL/. ; fi
	mv -f kmkhz6a.o NORMAL/kmkhz6a.o

NORMAL/leaf5b.o: leaf5b.f
	$(FORT) $(FORTOPTS)  -c leaf5b.f
	@if test -f leaf5b.T;then echo "	mv -f leaf5b.T NORMAL/.";	mv -f leaf5b.T NORMAL/. ; fi
	mv -f leaf5b.o NORMAL/leaf5b.o

NORMAL/phimh6a.o: phimh6a.f
	$(FORT) $(FORTOPTS)  -c phimh6a.f
	@if test -f phimh6a.T;then echo "	mv -f phimh6a.T NORMAL/.";	mv -f phimh6a.T NORMAL/. ; fi
	mv -f phimh6a.o NORMAL/phimh6a.o

NORMAL/sfevap5b.o: sfevap5b.f
	$(FORT) $(FORTOPTS)  -c sfevap5b.f
	@if test -f sfevap5b.T;then echo "	mv -f sfevap5b.T NORMAL/.";	mv -f sfevap5b.T NORMAL/. ; fi
	mv -f sfevap5b.o NORMAL/sfevap5b.o

NORMAL/sfexch6a.o: sfexch6a.f
	$(FORT) $(FORTOPTS)  -c sfexch6a.f
	@if test -f sfexch6a.T;then echo "	mv -f sfexch6a.T NORMAL/.";	mv -f sfexch6a.T NORMAL/. ; fi
	mv -f sfexch6a.o NORMAL/sfexch6a.o

NORMAL/sfflux5b.o: sfflux5b.f
	$(FORT) $(FORTOPTS)  -c sfflux5b.f
	@if test -f sfflux5b.T;then echo "	mv -f sfflux5b.T NORMAL/.";	mv -f sfflux5b.T NORMAL/. ; fi
	mv -f sfflux5b.o NORMAL/sfflux5b.o

NORMAL/sflbes5b.o: sflbes5b.f
	$(FORT) $(FORTOPTS)  -c sflbes5b.f
	@if test -f sflbes5b.T;then echo "	mv -f sflbes5b.T NORMAL/.";	mv -f sflbes5b.T NORMAL/. ; fi
	mv -f sflbes5b.o NORMAL/sflbes5b.o

NORMAL/sflint6a.o: sflint6a.f
	$(FORT) $(FORTOPTS)  -c sflint6a.f
	@if test -f sflint6a.T;then echo "	mv -f sflint6a.T NORMAL/.";	mv -f sflint6a.T NORMAL/. ; fi
	mv -f sflint6a.o NORMAL/sflint6a.o

NORMAL/sfmelt5b.o: sfmelt5b.f
	$(FORT) $(FORTOPTS)  -c sfmelt5b.f
	@if test -f sfmelt5b.T;then echo "	mv -f sfmelt5b.T NORMAL/.";	mv -f sfmelt5b.T NORMAL/. ; fi
	mv -f sfmelt5b.o NORMAL/sfmelt5b.o

NORMAL/sfrest6a.o: sfrest6a.f
	$(FORT) $(FORTOPTS)  -c sfrest6a.f
	@if test -f sfrest6a.T;then echo "	mv -f sfrest6a.T NORMAL/.";	mv -f sfrest6a.T NORMAL/. ; fi
	mv -f sfrest6a.o NORMAL/sfrest6a.o

NORMAL/sfrib6a.o: sfrib6a.f
	$(FORT) $(FORTOPTS)  -c sfrib6a.f
	@if test -f sfrib6a.T;then echo "	mv -f sfrib6a.T NORMAL/.";	mv -f sfrib6a.T NORMAL/. ; fi
	mv -f sfrib6a.o NORMAL/sfrib6a.o

NORMAL/sfrugh6a.o: sfrugh6a.f
	$(FORT) $(FORTOPTS)  -c sfrugh6a.f
	@if test -f sfrugh6a.T;then echo "	mv -f sfrugh6a.T NORMAL/.";	mv -f sfrugh6a.T NORMAL/. ; fi
	mv -f sfrugh6a.o NORMAL/sfrugh6a.o

NORMAL/sfstom5b.o: sfstom5b.f
	$(FORT) $(FORTOPTS)  -c sfstom5b.f
	@if test -f sfstom5b.T;then echo "	mv -f sfstom5b.T NORMAL/.";	mv -f sfstom5b.T NORMAL/. ; fi
	mv -f sfstom5b.o NORMAL/sfstom5b.o

NORMAL/siceht5b.o: siceht5b.f
	$(FORT) $(FORTOPTS)  -c siceht5b.f
	@if test -f siceht5b.T;then echo "	mv -f siceht5b.T NORMAL/.";	mv -f siceht5b.T NORMAL/. ; fi
	mv -f siceht5b.o NORMAL/siceht5b.o

NORMAL/trmix2c.o: trmix2c.f
	$(FORT) $(FORTOPTS)  -c trmix2c.f
	@if test -f trmix2c.T;then echo "	mv -f trmix2c.T NORMAL/.";	mv -f trmix2c.T NORMAL/. ; fi
	mv -f trmix2c.o NORMAL/trmix2c.o

NORMAL/trsrce1a.o: trsrce1a.f
	$(FORT) $(FORTOPTS)  -c trsrce1a.f
	@if test -f trsrce1a.T;then echo "	mv -f trsrce1a.T NORMAL/.";	mv -f trsrce1a.T NORMAL/. ; fi
	mv -f trsrce1a.o NORMAL/trsrce1a.o

NORMAL/z1a.o: z1a.f
	$(FORT) $(FORTOPTS)  -c z1a.f
	@if test -f z1a.T;then echo "	mv -f z1a.T NORMAL/.";	mv -f z1a.T NORMAL/. ; fi
	mv -f z1a.o NORMAL/z1a.o

