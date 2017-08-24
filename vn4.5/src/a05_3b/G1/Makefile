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
NORMAL/cal3dcca.o NORMAL/chgphs3b.o NORMAL/cloudw1a.o NORMAL/conrad1a.o\
NORMAL/consc1a.o NORMAL/conv23a.o NORMAL/convec3a.o NORMAL/coreng1a.o\
NORMAL/crsfrz2a.o NORMAL/ddcall3a.o NORMAL/ddenv3a.o NORMAL/ddinit3a.o\
NORMAL/ddraug3a.o NORMAL/detrai1a.o NORMAL/detrat1a.o NORMAL/devap3a.o\
NORMAL/downd3b.o NORMAL/dqsdth1a.o NORMAL/enviro3a.o NORMAL/evp3a.o\
NORMAL/flagw1a.o NORMAL/flxini2a.o NORMAL/gconv_3b.o NORMAL/latent1a.o\
NORMAL/laycn1a.o NORMAL/layerd2a.o NORMAL/lifpar3a.o NORMAL/parcel1a.o\
NORMAL/pevbcb3a.o NORMAL/satcal2a.o NORMAL/sconsc1a.o NORMAL/termco1a.o\
NORMAL/termdd2a.o NORMAL/thetar1a.o NORMAL/thpdet1a.o

all: $(OBJS)


NORMAL/cal3dcca.o: cal3dcca.f
	$(FORT) $(FORTOPTS)  -c cal3dcca.f
	@if test -f cal3dcca.T;then echo "	mv -f cal3dcca.T NORMAL/.";	mv -f cal3dcca.T NORMAL/. ; fi
	mv -f cal3dcca.o NORMAL/cal3dcca.o

NORMAL/chgphs3b.o: chgphs3b.f
	$(FORT) $(FORTOPTS)  -c chgphs3b.f
	@if test -f chgphs3b.T;then echo "	mv -f chgphs3b.T NORMAL/.";	mv -f chgphs3b.T NORMAL/. ; fi
	mv -f chgphs3b.o NORMAL/chgphs3b.o

NORMAL/cloudw1a.o: cloudw1a.f
	$(FORT) $(FORTOPTS)  -c cloudw1a.f
	@if test -f cloudw1a.T;then echo "	mv -f cloudw1a.T NORMAL/.";	mv -f cloudw1a.T NORMAL/. ; fi
	mv -f cloudw1a.o NORMAL/cloudw1a.o

NORMAL/conrad1a.o: conrad1a.f
	$(FORT) $(FORTOPTS)  -c conrad1a.f
	@if test -f conrad1a.T;then echo "	mv -f conrad1a.T NORMAL/.";	mv -f conrad1a.T NORMAL/. ; fi
	mv -f conrad1a.o NORMAL/conrad1a.o

NORMAL/consc1a.o: consc1a.f
	$(FORT) $(FORTOPTS)  -c consc1a.f
	@if test -f consc1a.T;then echo "	mv -f consc1a.T NORMAL/.";	mv -f consc1a.T NORMAL/. ; fi
	mv -f consc1a.o NORMAL/consc1a.o

NORMAL/conv23a.o: conv23a.f
	$(FORT) $(FORTOPTS)  -c conv23a.f
	@if test -f conv23a.T;then echo "	mv -f conv23a.T NORMAL/.";	mv -f conv23a.T NORMAL/. ; fi
	mv -f conv23a.o NORMAL/conv23a.o

NORMAL/convec3a.o: convec3a.f
	$(FORT) $(FORTOPTS)  -c convec3a.f
	@if test -f convec3a.T;then echo "	mv -f convec3a.T NORMAL/.";	mv -f convec3a.T NORMAL/. ; fi
	mv -f convec3a.o NORMAL/convec3a.o

NORMAL/coreng1a.o: coreng1a.f
	$(FORT) $(FORTOPTS)  -c coreng1a.f
	@if test -f coreng1a.T;then echo "	mv -f coreng1a.T NORMAL/.";	mv -f coreng1a.T NORMAL/. ; fi
	mv -f coreng1a.o NORMAL/coreng1a.o

NORMAL/crsfrz2a.o: crsfrz2a.f
	$(FORT) $(FORTOPTS)  -c crsfrz2a.f
	@if test -f crsfrz2a.T;then echo "	mv -f crsfrz2a.T NORMAL/.";	mv -f crsfrz2a.T NORMAL/. ; fi
	mv -f crsfrz2a.o NORMAL/crsfrz2a.o

NORMAL/ddcall3a.o: ddcall3a.f
	$(FORT) $(FORTOPTS)  -c ddcall3a.f
	@if test -f ddcall3a.T;then echo "	mv -f ddcall3a.T NORMAL/.";	mv -f ddcall3a.T NORMAL/. ; fi
	mv -f ddcall3a.o NORMAL/ddcall3a.o

NORMAL/ddenv3a.o: ddenv3a.f
	$(FORT) $(FORTOPTS)  -c ddenv3a.f
	@if test -f ddenv3a.T;then echo "	mv -f ddenv3a.T NORMAL/.";	mv -f ddenv3a.T NORMAL/. ; fi
	mv -f ddenv3a.o NORMAL/ddenv3a.o

NORMAL/ddinit3a.o: ddinit3a.f
	$(FORT) $(FORTOPTS)  -c ddinit3a.f
	@if test -f ddinit3a.T;then echo "	mv -f ddinit3a.T NORMAL/.";	mv -f ddinit3a.T NORMAL/. ; fi
	mv -f ddinit3a.o NORMAL/ddinit3a.o

NORMAL/ddraug3a.o: ddraug3a.f
	$(FORT) $(FORTOPTS)  -c ddraug3a.f
	@if test -f ddraug3a.T;then echo "	mv -f ddraug3a.T NORMAL/.";	mv -f ddraug3a.T NORMAL/. ; fi
	mv -f ddraug3a.o NORMAL/ddraug3a.o

NORMAL/detrai1a.o: detrai1a.f
	$(FORT) $(FORTOPTS)  -c detrai1a.f
	@if test -f detrai1a.T;then echo "	mv -f detrai1a.T NORMAL/.";	mv -f detrai1a.T NORMAL/. ; fi
	mv -f detrai1a.o NORMAL/detrai1a.o

NORMAL/detrat1a.o: detrat1a.f
	$(FORT) $(FORTOPTS)  -c detrat1a.f
	@if test -f detrat1a.T;then echo "	mv -f detrat1a.T NORMAL/.";	mv -f detrat1a.T NORMAL/. ; fi
	mv -f detrat1a.o NORMAL/detrat1a.o

NORMAL/devap3a.o: devap3a.f
	$(FORT) $(FORTOPTS)  -c devap3a.f
	@if test -f devap3a.T;then echo "	mv -f devap3a.T NORMAL/.";	mv -f devap3a.T NORMAL/. ; fi
	mv -f devap3a.o NORMAL/devap3a.o

NORMAL/downd3b.o: downd3b.f
	$(FORT) $(FORTOPTS)  -c downd3b.f
	@if test -f downd3b.T;then echo "	mv -f downd3b.T NORMAL/.";	mv -f downd3b.T NORMAL/. ; fi
	mv -f downd3b.o NORMAL/downd3b.o

NORMAL/dqsdth1a.o: dqsdth1a.f
	$(FORT) $(FORTOPTS)  -c dqsdth1a.f
	@if test -f dqsdth1a.T;then echo "	mv -f dqsdth1a.T NORMAL/.";	mv -f dqsdth1a.T NORMAL/. ; fi
	mv -f dqsdth1a.o NORMAL/dqsdth1a.o

NORMAL/enviro3a.o: enviro3a.f
	$(FORT) $(FORTOPTS)  -c enviro3a.f
	@if test -f enviro3a.T;then echo "	mv -f enviro3a.T NORMAL/.";	mv -f enviro3a.T NORMAL/. ; fi
	mv -f enviro3a.o NORMAL/enviro3a.o

NORMAL/evp3a.o: evp3a.f
	$(FORT) $(FORTOPTS)  -c evp3a.f
	@if test -f evp3a.T;then echo "	mv -f evp3a.T NORMAL/.";	mv -f evp3a.T NORMAL/. ; fi
	mv -f evp3a.o NORMAL/evp3a.o

NORMAL/flagw1a.o: flagw1a.f
	$(FORT) $(FORTOPTS)  -c flagw1a.f
	@if test -f flagw1a.T;then echo "	mv -f flagw1a.T NORMAL/.";	mv -f flagw1a.T NORMAL/. ; fi
	mv -f flagw1a.o NORMAL/flagw1a.o

NORMAL/flxini2a.o: flxini2a.f
	$(FORT) $(FORTOPTS)  -c flxini2a.f
	@if test -f flxini2a.T;then echo "	mv -f flxini2a.T NORMAL/.";	mv -f flxini2a.T NORMAL/. ; fi
	mv -f flxini2a.o NORMAL/flxini2a.o

NORMAL/gconv_3b.o: gconv_3b.f
	$(FORT) $(FORTOPTS)  -c gconv_3b.f
	@if test -f gconv_3b.T;then echo "	mv -f gconv_3b.T NORMAL/.";	mv -f gconv_3b.T NORMAL/. ; fi
	mv -f gconv_3b.o NORMAL/gconv_3b.o

NORMAL/latent1a.o: latent1a.f
	$(FORT) $(FORTOPTS)  -c latent1a.f
	@if test -f latent1a.T;then echo "	mv -f latent1a.T NORMAL/.";	mv -f latent1a.T NORMAL/. ; fi
	mv -f latent1a.o NORMAL/latent1a.o

NORMAL/laycn1a.o: laycn1a.f
	$(FORT) $(FORTOPTS)  -c laycn1a.f
	@if test -f laycn1a.T;then echo "	mv -f laycn1a.T NORMAL/.";	mv -f laycn1a.T NORMAL/. ; fi
	mv -f laycn1a.o NORMAL/laycn1a.o

NORMAL/layerd2a.o: layerd2a.f
	$(FORT) $(FORTOPTS)  -c layerd2a.f
	@if test -f layerd2a.T;then echo "	mv -f layerd2a.T NORMAL/.";	mv -f layerd2a.T NORMAL/. ; fi
	mv -f layerd2a.o NORMAL/layerd2a.o

NORMAL/lifpar3a.o: lifpar3a.f
	$(FORT) $(FORTOPTS)  -c lifpar3a.f
	@if test -f lifpar3a.T;then echo "	mv -f lifpar3a.T NORMAL/.";	mv -f lifpar3a.T NORMAL/. ; fi
	mv -f lifpar3a.o NORMAL/lifpar3a.o

NORMAL/parcel1a.o: parcel1a.f
	$(FORT) $(FORTOPTS)  -c parcel1a.f
	@if test -f parcel1a.T;then echo "	mv -f parcel1a.T NORMAL/.";	mv -f parcel1a.T NORMAL/. ; fi
	mv -f parcel1a.o NORMAL/parcel1a.o

NORMAL/pevbcb3a.o: pevbcb3a.f
	$(FORT) $(FORTOPTS)  -c pevbcb3a.f
	@if test -f pevbcb3a.T;then echo "	mv -f pevbcb3a.T NORMAL/.";	mv -f pevbcb3a.T NORMAL/. ; fi
	mv -f pevbcb3a.o NORMAL/pevbcb3a.o

NORMAL/satcal2a.o: satcal2a.f
	$(FORT) $(FORTOPTS)  -c satcal2a.f
	@if test -f satcal2a.T;then echo "	mv -f satcal2a.T NORMAL/.";	mv -f satcal2a.T NORMAL/. ; fi
	mv -f satcal2a.o NORMAL/satcal2a.o

NORMAL/sconsc1a.o: sconsc1a.f
	$(FORT) $(FORTOPTS)  -c sconsc1a.f
	@if test -f sconsc1a.T;then echo "	mv -f sconsc1a.T NORMAL/.";	mv -f sconsc1a.T NORMAL/. ; fi
	mv -f sconsc1a.o NORMAL/sconsc1a.o

NORMAL/termco1a.o: termco1a.f
	$(FORT) $(FORTOPTS)  -c termco1a.f
	@if test -f termco1a.T;then echo "	mv -f termco1a.T NORMAL/.";	mv -f termco1a.T NORMAL/. ; fi
	mv -f termco1a.o NORMAL/termco1a.o

NORMAL/termdd2a.o: termdd2a.f
	$(FORT) $(FORTOPTS)  -c termdd2a.f
	@if test -f termdd2a.T;then echo "	mv -f termdd2a.T NORMAL/.";	mv -f termdd2a.T NORMAL/. ; fi
	mv -f termdd2a.o NORMAL/termdd2a.o

NORMAL/thetar1a.o: thetar1a.f
	$(FORT) $(FORTOPTS)  -c thetar1a.f
	@if test -f thetar1a.T;then echo "	mv -f thetar1a.T NORMAL/.";	mv -f thetar1a.T NORMAL/. ; fi
	mv -f thetar1a.o NORMAL/thetar1a.o

NORMAL/thpdet1a.o: thpdet1a.f
	$(FORT) $(FORTOPTS)  -c thpdet1a.f
	@if test -f thpdet1a.T;then echo "	mv -f thpdet1a.T NORMAL/.";	mv -f thpdet1a.T NORMAL/. ; fi
	mv -f thpdet1a.o NORMAL/thpdet1a.o

