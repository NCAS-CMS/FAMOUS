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
NORMAL/cal3dcca.o NORMAL/chgphs3c.o NORMAL/cloudw1a.o NORMAL/conrad1a.o\
NORMAL/consc1a.o NORMAL/conv23a.o NORMAL/convec3c.o NORMAL/coreng1a.o\
NORMAL/crsfrz3c.o NORMAL/ddcall3c.o NORMAL/ddenv3a.o NORMAL/ddinit3a.o\
NORMAL/ddraug3c.o NORMAL/detrai1a.o NORMAL/detrat3c.o NORMAL/devap3a.o\
NORMAL/downd3c.o NORMAL/dqsdth3c.o NORMAL/enviro3a.o NORMAL/evp3a.o\
NORMAL/flagw3c.o NORMAL/flxini2a.o NORMAL/gconv_3a.o NORMAL/latent3c.o\
NORMAL/laycn3c.o NORMAL/layerd3c.o NORMAL/lifpar3a.o NORMAL/parcel3c.o\
NORMAL/pevbcb3a.o NORMAL/satcal3c.o NORMAL/sconsc1a.o NORMAL/termco1a.o\
NORMAL/termdd2a.o NORMAL/thetar1a.o NORMAL/thpdet1a.o

all: $(OBJS)


NORMAL/cal3dcca.o: cal3dcca.f
	$(FORT) $(FORTOPTS)  -c cal3dcca.f
	@if test -f cal3dcca.T;then echo "	mv -f cal3dcca.T NORMAL/.";	mv -f cal3dcca.T NORMAL/. ; fi
	mv -f cal3dcca.o NORMAL/cal3dcca.o

NORMAL/chgphs3c.o: chgphs3c.f
	$(FORT) $(FORTOPTS)  -c chgphs3c.f
	@if test -f chgphs3c.T;then echo "	mv -f chgphs3c.T NORMAL/.";	mv -f chgphs3c.T NORMAL/. ; fi
	mv -f chgphs3c.o NORMAL/chgphs3c.o

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

NORMAL/convec3c.o: convec3c.f
	$(FORT) $(FORTOPTS)  -c convec3c.f
	@if test -f convec3c.T;then echo "	mv -f convec3c.T NORMAL/.";	mv -f convec3c.T NORMAL/. ; fi
	mv -f convec3c.o NORMAL/convec3c.o

NORMAL/coreng1a.o: coreng1a.f
	$(FORT) $(FORTOPTS)  -c coreng1a.f
	@if test -f coreng1a.T;then echo "	mv -f coreng1a.T NORMAL/.";	mv -f coreng1a.T NORMAL/. ; fi
	mv -f coreng1a.o NORMAL/coreng1a.o

NORMAL/crsfrz3c.o: crsfrz3c.f
	$(FORT) $(FORTOPTS)  -c crsfrz3c.f
	@if test -f crsfrz3c.T;then echo "	mv -f crsfrz3c.T NORMAL/.";	mv -f crsfrz3c.T NORMAL/. ; fi
	mv -f crsfrz3c.o NORMAL/crsfrz3c.o

NORMAL/ddcall3c.o: ddcall3c.f
	$(FORT) $(FORTOPTS)  -c ddcall3c.f
	@if test -f ddcall3c.T;then echo "	mv -f ddcall3c.T NORMAL/.";	mv -f ddcall3c.T NORMAL/. ; fi
	mv -f ddcall3c.o NORMAL/ddcall3c.o

NORMAL/ddenv3a.o: ddenv3a.f
	$(FORT) $(FORTOPTS)  -c ddenv3a.f
	@if test -f ddenv3a.T;then echo "	mv -f ddenv3a.T NORMAL/.";	mv -f ddenv3a.T NORMAL/. ; fi
	mv -f ddenv3a.o NORMAL/ddenv3a.o

NORMAL/ddinit3a.o: ddinit3a.f
	$(FORT) $(FORTOPTS)  -c ddinit3a.f
	@if test -f ddinit3a.T;then echo "	mv -f ddinit3a.T NORMAL/.";	mv -f ddinit3a.T NORMAL/. ; fi
	mv -f ddinit3a.o NORMAL/ddinit3a.o

NORMAL/ddraug3c.o: ddraug3c.f
	$(FORT) $(FORTOPTS)  -c ddraug3c.f
	@if test -f ddraug3c.T;then echo "	mv -f ddraug3c.T NORMAL/.";	mv -f ddraug3c.T NORMAL/. ; fi
	mv -f ddraug3c.o NORMAL/ddraug3c.o

NORMAL/detrai1a.o: detrai1a.f
	$(FORT) $(FORTOPTS)  -c detrai1a.f
	@if test -f detrai1a.T;then echo "	mv -f detrai1a.T NORMAL/.";	mv -f detrai1a.T NORMAL/. ; fi
	mv -f detrai1a.o NORMAL/detrai1a.o

NORMAL/detrat3c.o: detrat3c.f
	$(FORT) $(FORTOPTS)  -c detrat3c.f
	@if test -f detrat3c.T;then echo "	mv -f detrat3c.T NORMAL/.";	mv -f detrat3c.T NORMAL/. ; fi
	mv -f detrat3c.o NORMAL/detrat3c.o

NORMAL/devap3a.o: devap3a.f
	$(FORT) $(FORTOPTS)  -c devap3a.f
	@if test -f devap3a.T;then echo "	mv -f devap3a.T NORMAL/.";	mv -f devap3a.T NORMAL/. ; fi
	mv -f devap3a.o NORMAL/devap3a.o

NORMAL/downd3c.o: downd3c.f
	$(FORT) $(FORTOPTS)  -c downd3c.f
	@if test -f downd3c.T;then echo "	mv -f downd3c.T NORMAL/.";	mv -f downd3c.T NORMAL/. ; fi
	mv -f downd3c.o NORMAL/downd3c.o

NORMAL/dqsdth3c.o: dqsdth3c.f
	$(FORT) $(FORTOPTS)  -c dqsdth3c.f
	@if test -f dqsdth3c.T;then echo "	mv -f dqsdth3c.T NORMAL/.";	mv -f dqsdth3c.T NORMAL/. ; fi
	mv -f dqsdth3c.o NORMAL/dqsdth3c.o

NORMAL/enviro3a.o: enviro3a.f
	$(FORT) $(FORTOPTS)  -c enviro3a.f
	@if test -f enviro3a.T;then echo "	mv -f enviro3a.T NORMAL/.";	mv -f enviro3a.T NORMAL/. ; fi
	mv -f enviro3a.o NORMAL/enviro3a.o

NORMAL/evp3a.o: evp3a.f
	$(FORT) $(FORTOPTS)  -c evp3a.f
	@if test -f evp3a.T;then echo "	mv -f evp3a.T NORMAL/.";	mv -f evp3a.T NORMAL/. ; fi
	mv -f evp3a.o NORMAL/evp3a.o

NORMAL/flagw3c.o: flagw3c.f
	$(FORT) $(FORTOPTS)  -c flagw3c.f
	@if test -f flagw3c.T;then echo "	mv -f flagw3c.T NORMAL/.";	mv -f flagw3c.T NORMAL/. ; fi
	mv -f flagw3c.o NORMAL/flagw3c.o

NORMAL/flxini2a.o: flxini2a.f
	$(FORT) $(FORTOPTS)  -c flxini2a.f
	@if test -f flxini2a.T;then echo "	mv -f flxini2a.T NORMAL/.";	mv -f flxini2a.T NORMAL/. ; fi
	mv -f flxini2a.o NORMAL/flxini2a.o

NORMAL/gconv_3a.o: gconv_3a.f
	$(FORT) $(FORTOPTS)  -c gconv_3a.f
	@if test -f gconv_3a.T;then echo "	mv -f gconv_3a.T NORMAL/.";	mv -f gconv_3a.T NORMAL/. ; fi
	mv -f gconv_3a.o NORMAL/gconv_3a.o

NORMAL/latent3c.o: latent3c.f
	$(FORT) $(FORTOPTS)  -c latent3c.f
	@if test -f latent3c.T;then echo "	mv -f latent3c.T NORMAL/.";	mv -f latent3c.T NORMAL/. ; fi
	mv -f latent3c.o NORMAL/latent3c.o

NORMAL/laycn3c.o: laycn3c.f
	$(FORT) $(FORTOPTS)  -c laycn3c.f
	@if test -f laycn3c.T;then echo "	mv -f laycn3c.T NORMAL/.";	mv -f laycn3c.T NORMAL/. ; fi
	mv -f laycn3c.o NORMAL/laycn3c.o

NORMAL/layerd3c.o: layerd3c.f
	$(FORT) $(FORTOPTS)  -c layerd3c.f
	@if test -f layerd3c.T;then echo "	mv -f layerd3c.T NORMAL/.";	mv -f layerd3c.T NORMAL/. ; fi
	mv -f layerd3c.o NORMAL/layerd3c.o

NORMAL/lifpar3a.o: lifpar3a.f
	$(FORT) $(FORTOPTS)  -c lifpar3a.f
	@if test -f lifpar3a.T;then echo "	mv -f lifpar3a.T NORMAL/.";	mv -f lifpar3a.T NORMAL/. ; fi
	mv -f lifpar3a.o NORMAL/lifpar3a.o

NORMAL/parcel3c.o: parcel3c.f
	$(FORT) $(FORTOPTS)  -c parcel3c.f
	@if test -f parcel3c.T;then echo "	mv -f parcel3c.T NORMAL/.";	mv -f parcel3c.T NORMAL/. ; fi
	mv -f parcel3c.o NORMAL/parcel3c.o

NORMAL/pevbcb3a.o: pevbcb3a.f
	$(FORT) $(FORTOPTS)  -c pevbcb3a.f
	@if test -f pevbcb3a.T;then echo "	mv -f pevbcb3a.T NORMAL/.";	mv -f pevbcb3a.T NORMAL/. ; fi
	mv -f pevbcb3a.o NORMAL/pevbcb3a.o

NORMAL/satcal3c.o: satcal3c.f
	$(FORT) $(FORTOPTS)  -c satcal3c.f
	@if test -f satcal3c.T;then echo "	mv -f satcal3c.T NORMAL/.";	mv -f satcal3c.T NORMAL/. ; fi
	mv -f satcal3c.o NORMAL/satcal3c.o

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

