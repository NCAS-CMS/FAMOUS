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
NORMAL/bdylyr3a.o NORMAL/bl_ic3a.o NORMAL/bl_lsp1a.o NORMAL/excoef3a.o\
NORMAL/fcdch3a.o NORMAL/implca2c.o NORMAL/impmix2c.o NORMAL/kmkh3a.o\
NORMAL/sfevap2b.o NORMAL/sfexch3a.o NORMAL/sflint3a.o NORMAL/siceht1a.o\
NORMAL/soilht1a.o NORMAL/trmix2c.o NORMAL/z1a.o

all: $(OBJS)


NORMAL/bdylyr3a.o: bdylyr3a.f
	$(FORT) $(FORTOPTS)  -c bdylyr3a.f
	@if test -f bdylyr3a.T;then echo "	mv -f bdylyr3a.T NORMAL/.";	mv -f bdylyr3a.T NORMAL/. ; fi
	mv -f bdylyr3a.o NORMAL/bdylyr3a.o

NORMAL/bl_ic3a.o: bl_ic3a.f
	$(FORT) $(FORTOPTS)  -c bl_ic3a.f
	@if test -f bl_ic3a.T;then echo "	mv -f bl_ic3a.T NORMAL/.";	mv -f bl_ic3a.T NORMAL/. ; fi
	mv -f bl_ic3a.o NORMAL/bl_ic3a.o

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

NORMAL/implca2c.o: implca2c.f
	$(FORT) $(FORTOPTS)  -c implca2c.f
	@if test -f implca2c.T;then echo "	mv -f implca2c.T NORMAL/.";	mv -f implca2c.T NORMAL/. ; fi
	mv -f implca2c.o NORMAL/implca2c.o

NORMAL/impmix2c.o: impmix2c.f
	$(FORT) $(FORTOPTS)  -c impmix2c.f
	@if test -f impmix2c.T;then echo "	mv -f impmix2c.T NORMAL/.";	mv -f impmix2c.T NORMAL/. ; fi
	mv -f impmix2c.o NORMAL/impmix2c.o

NORMAL/kmkh3a.o: kmkh3a.f
	$(FORT) $(FORTOPTS)  -c kmkh3a.f
	@if test -f kmkh3a.T;then echo "	mv -f kmkh3a.T NORMAL/.";	mv -f kmkh3a.T NORMAL/. ; fi
	mv -f kmkh3a.o NORMAL/kmkh3a.o

NORMAL/sfevap2b.o: sfevap2b.f
	$(FORT) $(FORTOPTS)  -c sfevap2b.f
	@if test -f sfevap2b.T;then echo "	mv -f sfevap2b.T NORMAL/.";	mv -f sfevap2b.T NORMAL/. ; fi
	mv -f sfevap2b.o NORMAL/sfevap2b.o

NORMAL/sfexch3a.o: sfexch3a.f
	$(FORT) $(FORTOPTS)  -c sfexch3a.f
	@if test -f sfexch3a.T;then echo "	mv -f sfexch3a.T NORMAL/.";	mv -f sfexch3a.T NORMAL/. ; fi
	mv -f sfexch3a.o NORMAL/sfexch3a.o

NORMAL/sflint3a.o: sflint3a.f
	$(FORT) $(FORTOPTS)  -c sflint3a.f
	@if test -f sflint3a.T;then echo "	mv -f sflint3a.T NORMAL/.";	mv -f sflint3a.T NORMAL/. ; fi
	mv -f sflint3a.o NORMAL/sflint3a.o

NORMAL/siceht1a.o: siceht1a.f
	$(FORT) $(FORTOPTS)  -c siceht1a.f
	@if test -f siceht1a.T;then echo "	mv -f siceht1a.T NORMAL/.";	mv -f siceht1a.T NORMAL/. ; fi
	mv -f siceht1a.o NORMAL/siceht1a.o

NORMAL/soilht1a.o: soilht1a.f
	$(FORT) $(FORTOPTS)  -c soilht1a.f
	@if test -f soilht1a.T;then echo "	mv -f soilht1a.T NORMAL/.";	mv -f soilht1a.T NORMAL/. ; fi
	mv -f soilht1a.o NORMAL/soilht1a.o

NORMAL/trmix2c.o: trmix2c.f
	$(FORT) $(FORTOPTS)  -c trmix2c.f
	@if test -f trmix2c.T;then echo "	mv -f trmix2c.T NORMAL/.";	mv -f trmix2c.T NORMAL/. ; fi
	mv -f trmix2c.o NORMAL/trmix2c.o

NORMAL/z1a.o: z1a.f
	$(FORT) $(FORTOPTS)  -c z1a.f
	@if test -f z1a.T;then echo "	mv -f z1a.T NORMAL/.";	mv -f z1a.T NORMAL/. ; fi
	mv -f z1a.o NORMAL/z1a.o

