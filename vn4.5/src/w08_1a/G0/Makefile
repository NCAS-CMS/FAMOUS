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
NORMAL/blokfld.o NORMAL/femean.o NORMAL/gradi.o NORMAL/implsch.o\
NORMAL/mblock.o NORMAL/mgrid.o NORMAL/mtabs.o NORMAL/mubuf.o\
NORMAL/prewind1.o NORMAL/semean.o NORMAL/setupwv.o NORMAL/unblok.o\
NORMAL/wav_to_d11.o NORMAL/wavfrst1.o NORMAL/wv_check.o NORMAL/wv_shape.o

all: $(OBJS)


NORMAL/blokfld.o: blokfld.f
	$(FORT) $(FORTOPTS)  -c blokfld.f
	@if test -f blokfld.T;then echo "	mv -f blokfld.T NORMAL/.";	mv -f blokfld.T NORMAL/. ; fi
	mv -f blokfld.o NORMAL/blokfld.o

NORMAL/femean.o: femean.f
	$(FORT) $(FORTOPTS)  -c femean.f
	@if test -f femean.T;then echo "	mv -f femean.T NORMAL/.";	mv -f femean.T NORMAL/. ; fi
	mv -f femean.o NORMAL/femean.o

NORMAL/gradi.o: gradi.f
	$(FORT) $(FORTOPTS)  -c gradi.f
	@if test -f gradi.T;then echo "	mv -f gradi.T NORMAL/.";	mv -f gradi.T NORMAL/. ; fi
	mv -f gradi.o NORMAL/gradi.o

NORMAL/implsch.o: implsch.f
	$(FORT) $(FORTOPTS)  -c implsch.f
	@if test -f implsch.T;then echo "	mv -f implsch.T NORMAL/.";	mv -f implsch.T NORMAL/. ; fi
	mv -f implsch.o NORMAL/implsch.o

NORMAL/mblock.o: mblock.f
	$(FORT) $(FORTOPTS)  -c mblock.f
	@if test -f mblock.T;then echo "	mv -f mblock.T NORMAL/.";	mv -f mblock.T NORMAL/. ; fi
	mv -f mblock.o NORMAL/mblock.o

NORMAL/mgrid.o: mgrid.f
	$(FORT) $(FORTOPTS)  -c mgrid.f
	@if test -f mgrid.T;then echo "	mv -f mgrid.T NORMAL/.";	mv -f mgrid.T NORMAL/. ; fi
	mv -f mgrid.o NORMAL/mgrid.o

NORMAL/mtabs.o: mtabs.f
	$(FORT) $(FORTOPTS)  -c mtabs.f
	@if test -f mtabs.T;then echo "	mv -f mtabs.T NORMAL/.";	mv -f mtabs.T NORMAL/. ; fi
	mv -f mtabs.o NORMAL/mtabs.o

NORMAL/mubuf.o: mubuf.f
	$(FORT) $(FORTOPTS)  -c mubuf.f
	@if test -f mubuf.T;then echo "	mv -f mubuf.T NORMAL/.";	mv -f mubuf.T NORMAL/. ; fi
	mv -f mubuf.o NORMAL/mubuf.o

NORMAL/prewind1.o: prewind1.f
	$(FORT) $(FORTOPTS)  -c prewind1.f
	@if test -f prewind1.T;then echo "	mv -f prewind1.T NORMAL/.";	mv -f prewind1.T NORMAL/. ; fi
	mv -f prewind1.o NORMAL/prewind1.o

NORMAL/semean.o: semean.f
	$(FORT) $(FORTOPTS)  -c semean.f
	@if test -f semean.T;then echo "	mv -f semean.T NORMAL/.";	mv -f semean.T NORMAL/. ; fi
	mv -f semean.o NORMAL/semean.o

NORMAL/setupwv.o: setupwv.f
	$(FORT) $(FORTOPTS)  -c setupwv.f
	@if test -f setupwv.T;then echo "	mv -f setupwv.T NORMAL/.";	mv -f setupwv.T NORMAL/. ; fi
	mv -f setupwv.o NORMAL/setupwv.o

NORMAL/unblok.o: unblok.f
	$(FORT) $(FORTOPTS)  -c unblok.f
	@if test -f unblok.T;then echo "	mv -f unblok.T NORMAL/.";	mv -f unblok.T NORMAL/. ; fi
	mv -f unblok.o NORMAL/unblok.o

NORMAL/wav_to_d11.o: wav_to_d11.f
	$(FORT) $(FORTOPTS)  -c wav_to_d11.f
	@if test -f wav_to_d11.T;then echo "	mv -f wav_to_d11.T NORMAL/.";	mv -f wav_to_d11.T NORMAL/. ; fi
	mv -f wav_to_d11.o NORMAL/wav_to_d11.o

NORMAL/wavfrst1.o: wavfrst1.f
	$(FORT) $(FORTOPTS)  -c wavfrst1.f
	@if test -f wavfrst1.T;then echo "	mv -f wavfrst1.T NORMAL/.";	mv -f wavfrst1.T NORMAL/. ; fi
	mv -f wavfrst1.o NORMAL/wavfrst1.o

NORMAL/wv_check.o: wv_check.f
	$(FORT) $(FORTOPTS)  -c wv_check.f
	@if test -f wv_check.T;then echo "	mv -f wv_check.T NORMAL/.";	mv -f wv_check.T NORMAL/. ; fi
	mv -f wv_check.o NORMAL/wv_check.o

NORMAL/wv_shape.o: wv_shape.f
	$(FORT) $(FORTOPTS)  -c wv_shape.f
	@if test -f wv_shape.T;then echo "	mv -f wv_shape.T NORMAL/.";	mv -f wv_shape.T NORMAL/. ; fi
	mv -f wv_shape.o NORMAL/wv_shape.o

