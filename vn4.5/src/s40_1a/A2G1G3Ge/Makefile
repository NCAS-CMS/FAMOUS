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
NORMAL/cu_to_h.o NORMAL/cv_to_h.o NORMAL/h_to_cu.o NORMAL/h_to_cv.o\
NORMAL/slbctl1a.o NORMAL/slbdif1a.o NORMAL/slbhc1a.o NORMAL/slbhca1a.o\
NORMAL/slbiadv.o NORMAL/slbicav.o NORMAL/slbice1a.o NORMAL/slbidrif.o\
NORMAL/slbidyn.o NORMAL/slbifree.o NORMAL/slbtadv.o NORMAL/slbtuv.o\
NORMAL/trna2s1a.o NORMAL/trns2a1a.o NORMAL/umslab1a.o NORMAL/uv_to_cu.o\
NORMAL/uv_to_cv.o NORMAL/uv_to_h.o

all: $(OBJS)


NORMAL/cu_to_h.o: cu_to_h.f
	$(FORT) $(FORTOPTS)  -c cu_to_h.f
	@if test -f cu_to_h.T;then echo "	mv -f cu_to_h.T NORMAL/.";	mv -f cu_to_h.T NORMAL/. ; fi
	mv -f cu_to_h.o NORMAL/cu_to_h.o

NORMAL/cv_to_h.o: cv_to_h.f
	$(FORT) $(FORTOPTS)  -c cv_to_h.f
	@if test -f cv_to_h.T;then echo "	mv -f cv_to_h.T NORMAL/.";	mv -f cv_to_h.T NORMAL/. ; fi
	mv -f cv_to_h.o NORMAL/cv_to_h.o

NORMAL/h_to_cu.o: h_to_cu.f
	$(FORT) $(FORTOPTS)  -c h_to_cu.f
	@if test -f h_to_cu.T;then echo "	mv -f h_to_cu.T NORMAL/.";	mv -f h_to_cu.T NORMAL/. ; fi
	mv -f h_to_cu.o NORMAL/h_to_cu.o

NORMAL/h_to_cv.o: h_to_cv.f
	$(FORT) $(FORTOPTS)  -c h_to_cv.f
	@if test -f h_to_cv.T;then echo "	mv -f h_to_cv.T NORMAL/.";	mv -f h_to_cv.T NORMAL/. ; fi
	mv -f h_to_cv.o NORMAL/h_to_cv.o

NORMAL/slbctl1a.o: slbctl1a.f
	$(FORT) $(FORTOPTS)  -c slbctl1a.f
	@if test -f slbctl1a.T;then echo "	mv -f slbctl1a.T NORMAL/.";	mv -f slbctl1a.T NORMAL/. ; fi
	mv -f slbctl1a.o NORMAL/slbctl1a.o

NORMAL/slbdif1a.o: slbdif1a.f
	$(FORT) $(FORTOPTS)  -c slbdif1a.f
	@if test -f slbdif1a.T;then echo "	mv -f slbdif1a.T NORMAL/.";	mv -f slbdif1a.T NORMAL/. ; fi
	mv -f slbdif1a.o NORMAL/slbdif1a.o

NORMAL/slbhc1a.o: slbhc1a.f
	$(FORT) $(FORTOPTS)  -c slbhc1a.f
	@if test -f slbhc1a.T;then echo "	mv -f slbhc1a.T NORMAL/.";	mv -f slbhc1a.T NORMAL/. ; fi
	mv -f slbhc1a.o NORMAL/slbhc1a.o

NORMAL/slbhca1a.o: slbhca1a.f
	$(FORT) $(FORTOPTS)  -c slbhca1a.f
	@if test -f slbhca1a.T;then echo "	mv -f slbhca1a.T NORMAL/.";	mv -f slbhca1a.T NORMAL/. ; fi
	mv -f slbhca1a.o NORMAL/slbhca1a.o

NORMAL/slbiadv.o: slbiadv.f
	$(FORT) $(FORTOPTS)  -c slbiadv.f
	@if test -f slbiadv.T;then echo "	mv -f slbiadv.T NORMAL/.";	mv -f slbiadv.T NORMAL/. ; fi
	mv -f slbiadv.o NORMAL/slbiadv.o

NORMAL/slbicav.o: slbicav.f
	$(FORT) $(FORTOPTS)  -c slbicav.f
	@if test -f slbicav.T;then echo "	mv -f slbicav.T NORMAL/.";	mv -f slbicav.T NORMAL/. ; fi
	mv -f slbicav.o NORMAL/slbicav.o

NORMAL/slbice1a.o: slbice1a.f
	$(FORT) $(FORTOPTS)  -c slbice1a.f
	@if test -f slbice1a.T;then echo "	mv -f slbice1a.T NORMAL/.";	mv -f slbice1a.T NORMAL/. ; fi
	mv -f slbice1a.o NORMAL/slbice1a.o

NORMAL/slbidrif.o: slbidrif.f
	$(FORT) $(FORTOPTS)  -c slbidrif.f
	@if test -f slbidrif.T;then echo "	mv -f slbidrif.T NORMAL/.";	mv -f slbidrif.T NORMAL/. ; fi
	mv -f slbidrif.o NORMAL/slbidrif.o

NORMAL/slbidyn.o: slbidyn.f
	$(FORT) $(FORTOPTS)  -c slbidyn.f
	@if test -f slbidyn.T;then echo "	mv -f slbidyn.T NORMAL/.";	mv -f slbidyn.T NORMAL/. ; fi
	mv -f slbidyn.o NORMAL/slbidyn.o

NORMAL/slbifree.o: slbifree.f
	$(FORT) $(FORTOPTS)  -c slbifree.f
	@if test -f slbifree.T;then echo "	mv -f slbifree.T NORMAL/.";	mv -f slbifree.T NORMAL/. ; fi
	mv -f slbifree.o NORMAL/slbifree.o

NORMAL/slbtadv.o: slbtadv.f
	$(FORT) $(FORTOPTS)  -c slbtadv.f
	@if test -f slbtadv.T;then echo "	mv -f slbtadv.T NORMAL/.";	mv -f slbtadv.T NORMAL/. ; fi
	mv -f slbtadv.o NORMAL/slbtadv.o

NORMAL/slbtuv.o: slbtuv.f
	$(FORT) $(FORTOPTS)  -c slbtuv.f
	@if test -f slbtuv.T;then echo "	mv -f slbtuv.T NORMAL/.";	mv -f slbtuv.T NORMAL/. ; fi
	mv -f slbtuv.o NORMAL/slbtuv.o

NORMAL/trna2s1a.o: trna2s1a.f
	$(FORT) $(FORTOPTS)  -c trna2s1a.f
	@if test -f trna2s1a.T;then echo "	mv -f trna2s1a.T NORMAL/.";	mv -f trna2s1a.T NORMAL/. ; fi
	mv -f trna2s1a.o NORMAL/trna2s1a.o

NORMAL/trns2a1a.o: trns2a1a.f
	$(FORT) $(FORTOPTS)  -c trns2a1a.f
	@if test -f trns2a1a.T;then echo "	mv -f trns2a1a.T NORMAL/.";	mv -f trns2a1a.T NORMAL/. ; fi
	mv -f trns2a1a.o NORMAL/trns2a1a.o

NORMAL/umslab1a.o: umslab1a.f
	$(FORT) $(FORTOPTS)  -c umslab1a.f
	@if test -f umslab1a.T;then echo "	mv -f umslab1a.T NORMAL/.";	mv -f umslab1a.T NORMAL/. ; fi
	mv -f umslab1a.o NORMAL/umslab1a.o

NORMAL/uv_to_cu.o: uv_to_cu.f
	$(FORT) $(FORTOPTS)  -c uv_to_cu.f
	@if test -f uv_to_cu.T;then echo "	mv -f uv_to_cu.T NORMAL/.";	mv -f uv_to_cu.T NORMAL/. ; fi
	mv -f uv_to_cu.o NORMAL/uv_to_cu.o

NORMAL/uv_to_cv.o: uv_to_cv.f
	$(FORT) $(FORTOPTS)  -c uv_to_cv.f
	@if test -f uv_to_cv.T;then echo "	mv -f uv_to_cv.T NORMAL/.";	mv -f uv_to_cv.T NORMAL/. ; fi
	mv -f uv_to_cv.o NORMAL/uv_to_cv.o

NORMAL/uv_to_h.o: uv_to_h.f
	$(FORT) $(FORTOPTS)  -c uv_to_h.f
	@if test -f uv_to_h.T;then echo "	mv -f uv_to_h.T NORMAL/.";	mv -f uv_to_h.T NORMAL/. ; fi
	mv -f uv_to_h.o NORMAL/uv_to_h.o

