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
NORMAL/bdylyr5a.o NORMAL/bl_ic5a.o NORMAL/bl_lsp1a.o NORMAL/excoef3a.o\
NORMAL/fcdch3a.o NORMAL/htcond5a.o NORMAL/implca4a.o NORMAL/impmix2c.o\
NORMAL/kmkh3a.o NORMAL/leaf5a.o NORMAL/sfevap5a.o NORMAL/sfexch5a.o\
NORMAL/sfflux5a.o NORMAL/sflint3a.o NORMAL/sfmelt5a.o NORMAL/sfrest5a.o\
NORMAL/sfrib5a.o NORMAL/sfrugh5a.o NORMAL/sfstom5a.o NORMAL/siceht4a.o\
NORMAL/trmix2c.o NORMAL/z1a.o

all: $(OBJS)


NORMAL/bdylyr5a.o: bdylyr5a.f
	$(FORT) $(FORTOPTS)  -c bdylyr5a.f
	@if test -f bdylyr5a.T;then echo "	mv -f bdylyr5a.T NORMAL/.";	mv -f bdylyr5a.T NORMAL/. ; fi
	mv -f bdylyr5a.o NORMAL/bdylyr5a.o

NORMAL/bl_ic5a.o: bl_ic5a.f
	$(FORT) $(FORTOPTS)  -c bl_ic5a.f
	@if test -f bl_ic5a.T;then echo "	mv -f bl_ic5a.T NORMAL/.";	mv -f bl_ic5a.T NORMAL/. ; fi
	mv -f bl_ic5a.o NORMAL/bl_ic5a.o

NORMAL/bl_lsp1a.o: bl_lsp1a.f
	$(FORT) $(FORTOPTS)  -c bl_lsp1a.f
	@if test -f bl_lsp1a.T;then echo "	mv -f bl_lsp1a.T NORMAL/.";	mv -f bl_lsp1a.T NORMAL/. ; fi
	mv -f bl_lsp1a.o NORMAL/bl_lsp1a.o

NORMAL/excoef3a.o: excoef3a.f
	$(FORT) $(FORTOPTS)  -c excoef3a.f
	@if test -f excoef3a.T;then echo "	mv -f excoef3a.T NORMAL/.";	mv -f excoef3a.T NORMAL/. ; fi
	mv -f excoef3a.o NORMAL/excoef3a.o

NORMAL/fcdch3a.o: fcdch3a.f
	$(FORT) $(FORTOPTS)  -c fcdch3a.f
	@if test -f fcdch3a.T;then echo "	mv -f fcdch3a.T NORMAL/.";	mv -f fcdch3a.T NORMAL/. ; fi
	mv -f fcdch3a.o NORMAL/fcdch3a.o

NORMAL/htcond5a.o: htcond5a.f
	$(FORT) $(FORTOPTS)  -c htcond5a.f
	@if test -f htcond5a.T;then echo "	mv -f htcond5a.T NORMAL/.";	mv -f htcond5a.T NORMAL/. ; fi
	mv -f htcond5a.o NORMAL/htcond5a.o

NORMAL/implca4a.o: implca4a.f
	$(FORT) $(FORTOPTS)  -c implca4a.f
	@if test -f implca4a.T;then echo "	mv -f implca4a.T NORMAL/.";	mv -f implca4a.T NORMAL/. ; fi
	mv -f implca4a.o NORMAL/implca4a.o

NORMAL/impmix2c.o: impmix2c.f
	$(FORT) $(FORTOPTS)  -c impmix2c.f
	@if test -f impmix2c.T;then echo "	mv -f impmix2c.T NORMAL/.";	mv -f impmix2c.T NORMAL/. ; fi
	mv -f impmix2c.o NORMAL/impmix2c.o

NORMAL/kmkh3a.o: kmkh3a.f
	$(FORT) $(FORTOPTS)  -c kmkh3a.f
	@if test -f kmkh3a.T;then echo "	mv -f kmkh3a.T NORMAL/.";	mv -f kmkh3a.T NORMAL/. ; fi
	mv -f kmkh3a.o NORMAL/kmkh3a.o

NORMAL/leaf5a.o: leaf5a.f
	$(FORT) $(FORTOPTS)  -c leaf5a.f
	@if test -f leaf5a.T;then echo "	mv -f leaf5a.T NORMAL/.";	mv -f leaf5a.T NORMAL/. ; fi
	mv -f leaf5a.o NORMAL/leaf5a.o

NORMAL/sfevap5a.o: sfevap5a.f
	$(FORT) $(FORTOPTS)  -c sfevap5a.f
	@if test -f sfevap5a.T;then echo "	mv -f sfevap5a.T NORMAL/.";	mv -f sfevap5a.T NORMAL/. ; fi
	mv -f sfevap5a.o NORMAL/sfevap5a.o

NORMAL/sfexch5a.o: sfexch5a.f
	$(FORT) $(FORTOPTS)  -c sfexch5a.f
	@if test -f sfexch5a.T;then echo "	mv -f sfexch5a.T NORMAL/.";	mv -f sfexch5a.T NORMAL/. ; fi
	mv -f sfexch5a.o NORMAL/sfexch5a.o

NORMAL/sfflux5a.o: sfflux5a.f
	$(FORT) $(FORTOPTS)  -c sfflux5a.f
	@if test -f sfflux5a.T;then echo "	mv -f sfflux5a.T NORMAL/.";	mv -f sfflux5a.T NORMAL/. ; fi
	mv -f sfflux5a.o NORMAL/sfflux5a.o

NORMAL/sflint3a.o: sflint3a.f
	$(FORT) $(FORTOPTS)  -c sflint3a.f
	@if test -f sflint3a.T;then echo "	mv -f sflint3a.T NORMAL/.";	mv -f sflint3a.T NORMAL/. ; fi
	mv -f sflint3a.o NORMAL/sflint3a.o

NORMAL/sfmelt5a.o: sfmelt5a.f
	$(FORT) $(FORTOPTS)  -c sfmelt5a.f
	@if test -f sfmelt5a.T;then echo "	mv -f sfmelt5a.T NORMAL/.";	mv -f sfmelt5a.T NORMAL/. ; fi
	mv -f sfmelt5a.o NORMAL/sfmelt5a.o

NORMAL/sfrest5a.o: sfrest5a.f
	$(FORT) $(FORTOPTS)  -c sfrest5a.f
	@if test -f sfrest5a.T;then echo "	mv -f sfrest5a.T NORMAL/.";	mv -f sfrest5a.T NORMAL/. ; fi
	mv -f sfrest5a.o NORMAL/sfrest5a.o

NORMAL/sfrib5a.o: sfrib5a.f
	$(FORT) $(FORTOPTS)  -c sfrib5a.f
	@if test -f sfrib5a.T;then echo "	mv -f sfrib5a.T NORMAL/.";	mv -f sfrib5a.T NORMAL/. ; fi
	mv -f sfrib5a.o NORMAL/sfrib5a.o

NORMAL/sfrugh5a.o: sfrugh5a.f
	$(FORT) $(FORTOPTS)  -c sfrugh5a.f
	@if test -f sfrugh5a.T;then echo "	mv -f sfrugh5a.T NORMAL/.";	mv -f sfrugh5a.T NORMAL/. ; fi
	mv -f sfrugh5a.o NORMAL/sfrugh5a.o

NORMAL/sfstom5a.o: sfstom5a.f
	$(FORT) $(FORTOPTS)  -c sfstom5a.f
	@if test -f sfstom5a.T;then echo "	mv -f sfstom5a.T NORMAL/.";	mv -f sfstom5a.T NORMAL/. ; fi
	mv -f sfstom5a.o NORMAL/sfstom5a.o

NORMAL/siceht4a.o: siceht4a.f
	$(FORT) $(FORTOPTS)  -c siceht4a.f
	@if test -f siceht4a.T;then echo "	mv -f siceht4a.T NORMAL/.";	mv -f siceht4a.T NORMAL/. ; fi
	mv -f siceht4a.o NORMAL/siceht4a.o

NORMAL/trmix2c.o: trmix2c.f
	$(FORT) $(FORTOPTS)  -c trmix2c.f
	@if test -f trmix2c.T;then echo "	mv -f trmix2c.T NORMAL/.";	mv -f trmix2c.T NORMAL/. ; fi
	mv -f trmix2c.o NORMAL/trmix2c.o

NORMAL/z1a.o: z1a.f
	$(FORT) $(FORTOPTS)  -c z1a.f
	@if test -f z1a.T;then echo "	mv -f z1a.T NORMAL/.";	mv -f z1a.T NORMAL/. ; fi
	mv -f z1a.o NORMAL/z1a.o

