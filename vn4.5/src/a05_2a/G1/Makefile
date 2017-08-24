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
NORMAL/chgphs2a.o NORMAL/cloudw1a.o NORMAL/conrad1a.o NORMAL/consc1a.o\
NORMAL/conv21a.o NORMAL/convec2a.o NORMAL/coreng1a.o NORMAL/crsfrz2a.o\
NORMAL/ddcall2a.o NORMAL/ddenv2a.o NORMAL/ddinit2a.o NORMAL/ddraug2a.o\
NORMAL/detrai1a.o NORMAL/detrat1a.o NORMAL/devap2a.o NORMAL/downd2a.o\
NORMAL/dqsdth1a.o NORMAL/enviro1a.o NORMAL/evp2a.o NORMAL/flagw1a.o\
NORMAL/flxini2a.o NORMAL/gconv_2a.o NORMAL/latent1a.o NORMAL/laycn1a.o\
NORMAL/layerd2a.o NORMAL/lifpar1a.o NORMAL/parcel1a.o NORMAL/pevbcb2a.o\
NORMAL/satcal2a.o NORMAL/termco1a.o NORMAL/termdd2a.o NORMAL/thetar1a.o\
NORMAL/thpdet1a.o

all: $(OBJS)


NORMAL/chgphs2a.o: chgphs2a.f
	$(FORT) $(FORTOPTS)  -c chgphs2a.f
	@if test -f chgphs2a.T;then echo "	mv -f chgphs2a.T NORMAL/.";	mv -f chgphs2a.T NORMAL/. ; fi
	mv -f chgphs2a.o NORMAL/chgphs2a.o

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

NORMAL/conv21a.o: conv21a.f
	$(FORT) $(FORTOPTS)  -c conv21a.f
	@if test -f conv21a.T;then echo "	mv -f conv21a.T NORMAL/.";	mv -f conv21a.T NORMAL/. ; fi
	mv -f conv21a.o NORMAL/conv21a.o

NORMAL/convec2a.o: convec2a.f
	$(FORT) $(FORTOPTS)  -c convec2a.f
	@if test -f convec2a.T;then echo "	mv -f convec2a.T NORMAL/.";	mv -f convec2a.T NORMAL/. ; fi
	mv -f convec2a.o NORMAL/convec2a.o

NORMAL/coreng1a.o: coreng1a.f
	$(FORT) $(FORTOPTS)  -c coreng1a.f
	@if test -f coreng1a.T;then echo "	mv -f coreng1a.T NORMAL/.";	mv -f coreng1a.T NORMAL/. ; fi
	mv -f coreng1a.o NORMAL/coreng1a.o

NORMAL/crsfrz2a.o: crsfrz2a.f
	$(FORT) $(FORTOPTS)  -c crsfrz2a.f
	@if test -f crsfrz2a.T;then echo "	mv -f crsfrz2a.T NORMAL/.";	mv -f crsfrz2a.T NORMAL/. ; fi
	mv -f crsfrz2a.o NORMAL/crsfrz2a.o

NORMAL/ddcall2a.o: ddcall2a.f
	$(FORT) $(FORTOPTS)  -c ddcall2a.f
	@if test -f ddcall2a.T;then echo "	mv -f ddcall2a.T NORMAL/.";	mv -f ddcall2a.T NORMAL/. ; fi
	mv -f ddcall2a.o NORMAL/ddcall2a.o

NORMAL/ddenv2a.o: ddenv2a.f
	$(FORT) $(FORTOPTS)  -c ddenv2a.f
	@if test -f ddenv2a.T;then echo "	mv -f ddenv2a.T NORMAL/.";	mv -f ddenv2a.T NORMAL/. ; fi
	mv -f ddenv2a.o NORMAL/ddenv2a.o

NORMAL/ddinit2a.o: ddinit2a.f
	$(FORT) $(FORTOPTS)  -c ddinit2a.f
	@if test -f ddinit2a.T;then echo "	mv -f ddinit2a.T NORMAL/.";	mv -f ddinit2a.T NORMAL/. ; fi
	mv -f ddinit2a.o NORMAL/ddinit2a.o

NORMAL/ddraug2a.o: ddraug2a.f
	$(FORT) $(FORTOPTS)  -c ddraug2a.f
	@if test -f ddraug2a.T;then echo "	mv -f ddraug2a.T NORMAL/.";	mv -f ddraug2a.T NORMAL/. ; fi
	mv -f ddraug2a.o NORMAL/ddraug2a.o

NORMAL/detrai1a.o: detrai1a.f
	$(FORT) $(FORTOPTS)  -c detrai1a.f
	@if test -f detrai1a.T;then echo "	mv -f detrai1a.T NORMAL/.";	mv -f detrai1a.T NORMAL/. ; fi
	mv -f detrai1a.o NORMAL/detrai1a.o

NORMAL/detrat1a.o: detrat1a.f
	$(FORT) $(FORTOPTS)  -c detrat1a.f
	@if test -f detrat1a.T;then echo "	mv -f detrat1a.T NORMAL/.";	mv -f detrat1a.T NORMAL/. ; fi
	mv -f detrat1a.o NORMAL/detrat1a.o

NORMAL/devap2a.o: devap2a.f
	$(FORT) $(FORTOPTS)  -c devap2a.f
	@if test -f devap2a.T;then echo "	mv -f devap2a.T NORMAL/.";	mv -f devap2a.T NORMAL/. ; fi
	mv -f devap2a.o NORMAL/devap2a.o

NORMAL/downd2a.o: downd2a.f
	$(FORT) $(FORTOPTS)  -c downd2a.f
	@if test -f downd2a.T;then echo "	mv -f downd2a.T NORMAL/.";	mv -f downd2a.T NORMAL/. ; fi
	mv -f downd2a.o NORMAL/downd2a.o

NORMAL/dqsdth1a.o: dqsdth1a.f
	$(FORT) $(FORTOPTS)  -c dqsdth1a.f
	@if test -f dqsdth1a.T;then echo "	mv -f dqsdth1a.T NORMAL/.";	mv -f dqsdth1a.T NORMAL/. ; fi
	mv -f dqsdth1a.o NORMAL/dqsdth1a.o

NORMAL/enviro1a.o: enviro1a.f
	$(FORT) $(FORTOPTS)  -c enviro1a.f
	@if test -f enviro1a.T;then echo "	mv -f enviro1a.T NORMAL/.";	mv -f enviro1a.T NORMAL/. ; fi
	mv -f enviro1a.o NORMAL/enviro1a.o

NORMAL/evp2a.o: evp2a.f
	$(FORT) $(FORTOPTS)  -c evp2a.f
	@if test -f evp2a.T;then echo "	mv -f evp2a.T NORMAL/.";	mv -f evp2a.T NORMAL/. ; fi
	mv -f evp2a.o NORMAL/evp2a.o

NORMAL/flagw1a.o: flagw1a.f
	$(FORT) $(FORTOPTS)  -c flagw1a.f
	@if test -f flagw1a.T;then echo "	mv -f flagw1a.T NORMAL/.";	mv -f flagw1a.T NORMAL/. ; fi
	mv -f flagw1a.o NORMAL/flagw1a.o

NORMAL/flxini2a.o: flxini2a.f
	$(FORT) $(FORTOPTS)  -c flxini2a.f
	@if test -f flxini2a.T;then echo "	mv -f flxini2a.T NORMAL/.";	mv -f flxini2a.T NORMAL/. ; fi
	mv -f flxini2a.o NORMAL/flxini2a.o

NORMAL/gconv_2a.o: gconv_2a.f
	$(FORT) $(FORTOPTS)  -c gconv_2a.f
	@if test -f gconv_2a.T;then echo "	mv -f gconv_2a.T NORMAL/.";	mv -f gconv_2a.T NORMAL/. ; fi
	mv -f gconv_2a.o NORMAL/gconv_2a.o

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

NORMAL/lifpar1a.o: lifpar1a.f
	$(FORT) $(FORTOPTS)  -c lifpar1a.f
	@if test -f lifpar1a.T;then echo "	mv -f lifpar1a.T NORMAL/.";	mv -f lifpar1a.T NORMAL/. ; fi
	mv -f lifpar1a.o NORMAL/lifpar1a.o

NORMAL/parcel1a.o: parcel1a.f
	$(FORT) $(FORTOPTS)  -c parcel1a.f
	@if test -f parcel1a.T;then echo "	mv -f parcel1a.T NORMAL/.";	mv -f parcel1a.T NORMAL/. ; fi
	mv -f parcel1a.o NORMAL/parcel1a.o

NORMAL/pevbcb2a.o: pevbcb2a.f
	$(FORT) $(FORTOPTS)  -c pevbcb2a.f
	@if test -f pevbcb2a.T;then echo "	mv -f pevbcb2a.T NORMAL/.";	mv -f pevbcb2a.T NORMAL/. ; fi
	mv -f pevbcb2a.o NORMAL/pevbcb2a.o

NORMAL/satcal2a.o: satcal2a.f
	$(FORT) $(FORTOPTS)  -c satcal2a.f
	@if test -f satcal2a.T;then echo "	mv -f satcal2a.T NORMAL/.";	mv -f satcal2a.T NORMAL/. ; fi
	mv -f satcal2a.o NORMAL/satcal2a.o

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

