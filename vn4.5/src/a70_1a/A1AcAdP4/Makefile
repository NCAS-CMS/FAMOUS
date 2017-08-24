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
NORMAL/agfx3a.o NORMAL/aggcld3a.o NORMAL/agtfx3a.o NORMAL/asfx3a.o\
NORMAL/bnsol3a.o NORMAL/caden3a.o NORMAL/clcol3a.o NORMAL/clrsp3a.o\
NORMAL/colsf3a.o NORMAL/dfpln3a.o NORMAL/diag3a.o NORMAL/e3a013a.o\
NORMAL/fill3a.o NORMAL/fxca3a.o NORMAL/gascal1a.o NORMAL/gasop3a.o\
NORMAL/glbct3a.o NORMAL/greyk3a.o NORMAL/gssan3a.o NORMAL/intfx3a.o\
NORMAL/irs3a.o NORMAL/lcldn3a.o NORMAL/mngfx3a.o NORMAL/mnirr3a.o\
NORMAL/monfx3a.o NORMAL/mxaps3a.o NORMAL/mxcof3a.o NORMAL/mxcol3a.o\
NORMAL/mxmte3a.o NORMAL/mxssr3a.o NORMAL/opicl3a.o NORMAL/opwcl3a.o\
NORMAL/ovltrp3a.o NORMAL/ovmxm3a.o NORMAL/ovmxr3a.o NORMAL/rank3a.o\
NORMAL/rscas3a.o NORMAL/rscnt3a.o NORMAL/rsctw3a.o NORMAL/sb1g3a.o\
NORMAL/sbcf3a.o NORMAL/sbf3a.o NORMAL/sbke3a.o NORMAL/sbrv3a.o\
NORMAL/sbwg3a.o NORMAL/sclab3a.o NORMAL/sclge3a.o NORMAL/sclpt3a.o\
NORMAL/slhgdr3a.o NORMAL/slmxdr3a.o NORMAL/slvtpa3a.o NORMAL/slvtrp3a.o\
NORMAL/smap3a.o NORMAL/smtpn3a.o NORMAL/snclp3a.o NORMAL/snscf3a.o\
NORMAL/solcf3a.o NORMAL/solsr3a.o NORMAL/spin3a.o NORMAL/splmx3a.o\
NORMAL/ssct3a.o NORMAL/sscta3a.o NORMAL/stmtf3a.o NORMAL/stsct3a.o\
NORMAL/surfp3a.o NORMAL/tcf3a.o NORMAL/tcfb3a.o NORMAL/tcfc3a.o\
NORMAL/tropin1a.o NORMAL/trpclm3a.o NORMAL/trpssr3a.o NORMAL/trsfc3a.o\
NORMAL/tstr3a.o NORMAL/twcfrg3a.o

all: $(OBJS)


NORMAL/agfx3a.o: agfx3a.f
	$(FORT) $(FORTOPTS)  -c agfx3a.f
	@if test -f agfx3a.T;then echo "	mv -f agfx3a.T NORMAL/.";	mv -f agfx3a.T NORMAL/. ; fi
	mv -f agfx3a.o NORMAL/agfx3a.o

NORMAL/aggcld3a.o: aggcld3a.f
	$(FORT) $(FORTOPTS)  -c aggcld3a.f
	@if test -f aggcld3a.T;then echo "	mv -f aggcld3a.T NORMAL/.";	mv -f aggcld3a.T NORMAL/. ; fi
	mv -f aggcld3a.o NORMAL/aggcld3a.o

NORMAL/agtfx3a.o: agtfx3a.f
	$(FORT) $(FORTOPTS)  -c agtfx3a.f
	@if test -f agtfx3a.T;then echo "	mv -f agtfx3a.T NORMAL/.";	mv -f agtfx3a.T NORMAL/. ; fi
	mv -f agtfx3a.o NORMAL/agtfx3a.o

NORMAL/asfx3a.o: asfx3a.f
	$(FORT) $(FORTOPTS)  -c asfx3a.f
	@if test -f asfx3a.T;then echo "	mv -f asfx3a.T NORMAL/.";	mv -f asfx3a.T NORMAL/. ; fi
	mv -f asfx3a.o NORMAL/asfx3a.o

NORMAL/bnsol3a.o: bnsol3a.f
	$(FORT) $(FORTOPTS)  -c bnsol3a.f
	@if test -f bnsol3a.T;then echo "	mv -f bnsol3a.T NORMAL/.";	mv -f bnsol3a.T NORMAL/. ; fi
	mv -f bnsol3a.o NORMAL/bnsol3a.o

NORMAL/caden3a.o: caden3a.f
	$(FORT) $(FORTOPTS)  -c caden3a.f
	@if test -f caden3a.T;then echo "	mv -f caden3a.T NORMAL/.";	mv -f caden3a.T NORMAL/. ; fi
	mv -f caden3a.o NORMAL/caden3a.o

NORMAL/clcol3a.o: clcol3a.f
	$(FORT) $(FORTOPTS)  -c clcol3a.f
	@if test -f clcol3a.T;then echo "	mv -f clcol3a.T NORMAL/.";	mv -f clcol3a.T NORMAL/. ; fi
	mv -f clcol3a.o NORMAL/clcol3a.o

NORMAL/clrsp3a.o: clrsp3a.f
	$(FORT) $(FORTOPTS)  -c clrsp3a.f
	@if test -f clrsp3a.T;then echo "	mv -f clrsp3a.T NORMAL/.";	mv -f clrsp3a.T NORMAL/. ; fi
	mv -f clrsp3a.o NORMAL/clrsp3a.o

NORMAL/colsf3a.o: colsf3a.f
	$(FORT) $(FORTOPTS)  -c colsf3a.f
	@if test -f colsf3a.T;then echo "	mv -f colsf3a.T NORMAL/.";	mv -f colsf3a.T NORMAL/. ; fi
	mv -f colsf3a.o NORMAL/colsf3a.o

NORMAL/dfpln3a.o: dfpln3a.f
	$(FORT) $(FORTOPTS)  -c dfpln3a.f
	@if test -f dfpln3a.T;then echo "	mv -f dfpln3a.T NORMAL/.";	mv -f dfpln3a.T NORMAL/. ; fi
	mv -f dfpln3a.o NORMAL/dfpln3a.o

NORMAL/diag3a.o: diag3a.f
	$(FORT) $(FORTOPTS)  -c diag3a.f
	@if test -f diag3a.T;then echo "	mv -f diag3a.T NORMAL/.";	mv -f diag3a.T NORMAL/. ; fi
	mv -f diag3a.o NORMAL/diag3a.o

NORMAL/e3a013a.o: e3a013a.f
	$(FORT) $(FORTOPTS)  -c e3a013a.f
	@if test -f e3a013a.T;then echo "	mv -f e3a013a.T NORMAL/.";	mv -f e3a013a.T NORMAL/. ; fi
	mv -f e3a013a.o NORMAL/e3a013a.o

NORMAL/fill3a.o: fill3a.f
	$(FORT) $(FORTOPTS)  -c fill3a.f
	@if test -f fill3a.T;then echo "	mv -f fill3a.T NORMAL/.";	mv -f fill3a.T NORMAL/. ; fi
	mv -f fill3a.o NORMAL/fill3a.o

NORMAL/fxca3a.o: fxca3a.f
	$(FORT) $(FORTOPTS)  -c fxca3a.f
	@if test -f fxca3a.T;then echo "	mv -f fxca3a.T NORMAL/.";	mv -f fxca3a.T NORMAL/. ; fi
	mv -f fxca3a.o NORMAL/fxca3a.o

NORMAL/gascal1a.o: gascal1a.f
	$(FORT) $(FORTOPTS)  -c gascal1a.f
	@if test -f gascal1a.T;then echo "	mv -f gascal1a.T NORMAL/.";	mv -f gascal1a.T NORMAL/. ; fi
	mv -f gascal1a.o NORMAL/gascal1a.o

NORMAL/gasop3a.o: gasop3a.f
	$(FORT) $(FORTOPTS)  -c gasop3a.f
	@if test -f gasop3a.T;then echo "	mv -f gasop3a.T NORMAL/.";	mv -f gasop3a.T NORMAL/. ; fi
	mv -f gasop3a.o NORMAL/gasop3a.o

NORMAL/glbct3a.o: glbct3a.f
	$(FORT) $(FORTOPTS)  -c glbct3a.f
	@if test -f glbct3a.T;then echo "	mv -f glbct3a.T NORMAL/.";	mv -f glbct3a.T NORMAL/. ; fi
	mv -f glbct3a.o NORMAL/glbct3a.o

NORMAL/greyk3a.o: greyk3a.f
	$(FORT) $(FORTOPTS)  -c greyk3a.f
	@if test -f greyk3a.T;then echo "	mv -f greyk3a.T NORMAL/.";	mv -f greyk3a.T NORMAL/. ; fi
	mv -f greyk3a.o NORMAL/greyk3a.o

NORMAL/gssan3a.o: gssan3a.f
	$(FORT) $(FORTOPTS)  -c gssan3a.f
	@if test -f gssan3a.T;then echo "	mv -f gssan3a.T NORMAL/.";	mv -f gssan3a.T NORMAL/. ; fi
	mv -f gssan3a.o NORMAL/gssan3a.o

NORMAL/intfx3a.o: intfx3a.f
	$(FORT) $(FORTOPTS)  -c intfx3a.f
	@if test -f intfx3a.T;then echo "	mv -f intfx3a.T NORMAL/.";	mv -f intfx3a.T NORMAL/. ; fi
	mv -f intfx3a.o NORMAL/intfx3a.o

NORMAL/irs3a.o: irs3a.f
	$(FORT) $(FORTOPTS)  -c irs3a.f
	@if test -f irs3a.T;then echo "	mv -f irs3a.T NORMAL/.";	mv -f irs3a.T NORMAL/. ; fi
	mv -f irs3a.o NORMAL/irs3a.o

NORMAL/lcldn3a.o: lcldn3a.f
	$(FORT) $(FORTOPTS)  -c lcldn3a.f
	@if test -f lcldn3a.T;then echo "	mv -f lcldn3a.T NORMAL/.";	mv -f lcldn3a.T NORMAL/. ; fi
	mv -f lcldn3a.o NORMAL/lcldn3a.o

NORMAL/mngfx3a.o: mngfx3a.f
	$(FORT) $(FORTOPTS)  -c mngfx3a.f
	@if test -f mngfx3a.T;then echo "	mv -f mngfx3a.T NORMAL/.";	mv -f mngfx3a.T NORMAL/. ; fi
	mv -f mngfx3a.o NORMAL/mngfx3a.o

NORMAL/mnirr3a.o: mnirr3a.f
	$(FORT) $(FORTOPTS)  -c mnirr3a.f
	@if test -f mnirr3a.T;then echo "	mv -f mnirr3a.T NORMAL/.";	mv -f mnirr3a.T NORMAL/. ; fi
	mv -f mnirr3a.o NORMAL/mnirr3a.o

NORMAL/monfx3a.o: monfx3a.f
	$(FORT) $(FORTOPTS)  -c monfx3a.f
	@if test -f monfx3a.T;then echo "	mv -f monfx3a.T NORMAL/.";	mv -f monfx3a.T NORMAL/. ; fi
	mv -f monfx3a.o NORMAL/monfx3a.o

NORMAL/mxaps3a.o: mxaps3a.f
	$(FORT) $(FORTOPTS)  -c mxaps3a.f
	@if test -f mxaps3a.T;then echo "	mv -f mxaps3a.T NORMAL/.";	mv -f mxaps3a.T NORMAL/. ; fi
	mv -f mxaps3a.o NORMAL/mxaps3a.o

NORMAL/mxcof3a.o: mxcof3a.f
	$(FORT) $(FORTOPTS)  -c mxcof3a.f
	@if test -f mxcof3a.T;then echo "	mv -f mxcof3a.T NORMAL/.";	mv -f mxcof3a.T NORMAL/. ; fi
	mv -f mxcof3a.o NORMAL/mxcof3a.o

NORMAL/mxcol3a.o: mxcol3a.f
	$(FORT) $(FORTOPTS)  -c mxcol3a.f
	@if test -f mxcol3a.T;then echo "	mv -f mxcol3a.T NORMAL/.";	mv -f mxcol3a.T NORMAL/. ; fi
	mv -f mxcol3a.o NORMAL/mxcol3a.o

NORMAL/mxmte3a.o: mxmte3a.f
	$(FORT) $(FORTOPTS)  -c mxmte3a.f
	@if test -f mxmte3a.T;then echo "	mv -f mxmte3a.T NORMAL/.";	mv -f mxmte3a.T NORMAL/. ; fi
	mv -f mxmte3a.o NORMAL/mxmte3a.o

NORMAL/mxssr3a.o: mxssr3a.f
	$(FORT) $(FORTOPTS)  -c mxssr3a.f
	@if test -f mxssr3a.T;then echo "	mv -f mxssr3a.T NORMAL/.";	mv -f mxssr3a.T NORMAL/. ; fi
	mv -f mxssr3a.o NORMAL/mxssr3a.o

NORMAL/opicl3a.o: opicl3a.f
	$(FORT) $(FORTOPTS)  -c opicl3a.f
	@if test -f opicl3a.T;then echo "	mv -f opicl3a.T NORMAL/.";	mv -f opicl3a.T NORMAL/. ; fi
	mv -f opicl3a.o NORMAL/opicl3a.o

NORMAL/opwcl3a.o: opwcl3a.f
	$(FORT) $(FORTOPTS)  -c opwcl3a.f
	@if test -f opwcl3a.T;then echo "	mv -f opwcl3a.T NORMAL/.";	mv -f opwcl3a.T NORMAL/. ; fi
	mv -f opwcl3a.o NORMAL/opwcl3a.o

NORMAL/ovltrp3a.o: ovltrp3a.f
	$(FORT) $(FORTOPTS)  -c ovltrp3a.f
	@if test -f ovltrp3a.T;then echo "	mv -f ovltrp3a.T NORMAL/.";	mv -f ovltrp3a.T NORMAL/. ; fi
	mv -f ovltrp3a.o NORMAL/ovltrp3a.o

NORMAL/ovmxm3a.o: ovmxm3a.f
	$(FORT) $(FORTOPTS)  -c ovmxm3a.f
	@if test -f ovmxm3a.T;then echo "	mv -f ovmxm3a.T NORMAL/.";	mv -f ovmxm3a.T NORMAL/. ; fi
	mv -f ovmxm3a.o NORMAL/ovmxm3a.o

NORMAL/ovmxr3a.o: ovmxr3a.f
	$(FORT) $(FORTOPTS)  -c ovmxr3a.f
	@if test -f ovmxr3a.T;then echo "	mv -f ovmxr3a.T NORMAL/.";	mv -f ovmxr3a.T NORMAL/. ; fi
	mv -f ovmxr3a.o NORMAL/ovmxr3a.o

NORMAL/rank3a.o: rank3a.f
	$(FORT) $(FORTOPTS)  -c rank3a.f
	@if test -f rank3a.T;then echo "	mv -f rank3a.T NORMAL/.";	mv -f rank3a.T NORMAL/. ; fi
	mv -f rank3a.o NORMAL/rank3a.o

NORMAL/rscas3a.o: rscas3a.f
	$(FORT) $(FORTOPTS)  -c rscas3a.f
	@if test -f rscas3a.T;then echo "	mv -f rscas3a.T NORMAL/.";	mv -f rscas3a.T NORMAL/. ; fi
	mv -f rscas3a.o NORMAL/rscas3a.o

NORMAL/rscnt3a.o: rscnt3a.f
	$(FORT) $(FORTOPTS)  -c rscnt3a.f
	@if test -f rscnt3a.T;then echo "	mv -f rscnt3a.T NORMAL/.";	mv -f rscnt3a.T NORMAL/. ; fi
	mv -f rscnt3a.o NORMAL/rscnt3a.o

NORMAL/rsctw3a.o: rsctw3a.f
	$(FORT) $(FORTOPTS)  -c rsctw3a.f
	@if test -f rsctw3a.T;then echo "	mv -f rsctw3a.T NORMAL/.";	mv -f rsctw3a.T NORMAL/. ; fi
	mv -f rsctw3a.o NORMAL/rsctw3a.o

NORMAL/sb1g3a.o: sb1g3a.f
	$(FORT) $(FORTOPTS)  -c sb1g3a.f
	@if test -f sb1g3a.T;then echo "	mv -f sb1g3a.T NORMAL/.";	mv -f sb1g3a.T NORMAL/. ; fi
	mv -f sb1g3a.o NORMAL/sb1g3a.o

NORMAL/sbcf3a.o: sbcf3a.f
	$(FORT) $(FORTOPTS)  -c sbcf3a.f
	@if test -f sbcf3a.T;then echo "	mv -f sbcf3a.T NORMAL/.";	mv -f sbcf3a.T NORMAL/. ; fi
	mv -f sbcf3a.o NORMAL/sbcf3a.o

NORMAL/sbf3a.o: sbf3a.f
	$(FORT) $(FORTOPTS)  -c sbf3a.f
	@if test -f sbf3a.T;then echo "	mv -f sbf3a.T NORMAL/.";	mv -f sbf3a.T NORMAL/. ; fi
	mv -f sbf3a.o NORMAL/sbf3a.o

NORMAL/sbke3a.o: sbke3a.f
	$(FORT) $(FORTOPTS)  -c sbke3a.f
	@if test -f sbke3a.T;then echo "	mv -f sbke3a.T NORMAL/.";	mv -f sbke3a.T NORMAL/. ; fi
	mv -f sbke3a.o NORMAL/sbke3a.o

NORMAL/sbrv3a.o: sbrv3a.f
	$(FORT) $(FORTOPTS)  -c sbrv3a.f
	@if test -f sbrv3a.T;then echo "	mv -f sbrv3a.T NORMAL/.";	mv -f sbrv3a.T NORMAL/. ; fi
	mv -f sbrv3a.o NORMAL/sbrv3a.o

NORMAL/sbwg3a.o: sbwg3a.f
	$(FORT) $(FORTOPTS)  -c sbwg3a.f
	@if test -f sbwg3a.T;then echo "	mv -f sbwg3a.T NORMAL/.";	mv -f sbwg3a.T NORMAL/. ; fi
	mv -f sbwg3a.o NORMAL/sbwg3a.o

NORMAL/sclab3a.o: sclab3a.f
	$(FORT) $(FORTOPTS)  -c sclab3a.f
	@if test -f sclab3a.T;then echo "	mv -f sclab3a.T NORMAL/.";	mv -f sclab3a.T NORMAL/. ; fi
	mv -f sclab3a.o NORMAL/sclab3a.o

NORMAL/sclge3a.o: sclge3a.f
	$(FORT) $(FORTOPTS)  -c sclge3a.f
	@if test -f sclge3a.T;then echo "	mv -f sclge3a.T NORMAL/.";	mv -f sclge3a.T NORMAL/. ; fi
	mv -f sclge3a.o NORMAL/sclge3a.o

NORMAL/sclpt3a.o: sclpt3a.f
	$(FORT) $(FORTOPTS)  -c sclpt3a.f
	@if test -f sclpt3a.T;then echo "	mv -f sclpt3a.T NORMAL/.";	mv -f sclpt3a.T NORMAL/. ; fi
	mv -f sclpt3a.o NORMAL/sclpt3a.o

NORMAL/slhgdr3a.o: slhgdr3a.f
	$(FORT) $(FORTOPTS)  -c slhgdr3a.f
	@if test -f slhgdr3a.T;then echo "	mv -f slhgdr3a.T NORMAL/.";	mv -f slhgdr3a.T NORMAL/. ; fi
	mv -f slhgdr3a.o NORMAL/slhgdr3a.o

NORMAL/slmxdr3a.o: slmxdr3a.f
	$(FORT) $(FORTOPTS)  -c slmxdr3a.f
	@if test -f slmxdr3a.T;then echo "	mv -f slmxdr3a.T NORMAL/.";	mv -f slmxdr3a.T NORMAL/. ; fi
	mv -f slmxdr3a.o NORMAL/slmxdr3a.o

NORMAL/slvtpa3a.o: slvtpa3a.f
	$(FORT) $(FORTOPTS)  -c slvtpa3a.f
	@if test -f slvtpa3a.T;then echo "	mv -f slvtpa3a.T NORMAL/.";	mv -f slvtpa3a.T NORMAL/. ; fi
	mv -f slvtpa3a.o NORMAL/slvtpa3a.o

NORMAL/slvtrp3a.o: slvtrp3a.f
	$(FORT) $(FORTOPTS)  -c slvtrp3a.f
	@if test -f slvtrp3a.T;then echo "	mv -f slvtrp3a.T NORMAL/.";	mv -f slvtrp3a.T NORMAL/. ; fi
	mv -f slvtrp3a.o NORMAL/slvtrp3a.o

NORMAL/smap3a.o: smap3a.f
	$(FORT) $(FORTOPTS)  -c smap3a.f
	@if test -f smap3a.T;then echo "	mv -f smap3a.T NORMAL/.";	mv -f smap3a.T NORMAL/. ; fi
	mv -f smap3a.o NORMAL/smap3a.o

NORMAL/smtpn3a.o: smtpn3a.f
	$(FORT) $(FORTOPTS)  -c smtpn3a.f
	@if test -f smtpn3a.T;then echo "	mv -f smtpn3a.T NORMAL/.";	mv -f smtpn3a.T NORMAL/. ; fi
	mv -f smtpn3a.o NORMAL/smtpn3a.o

NORMAL/snclp3a.o: snclp3a.f
	$(FORT) $(FORTOPTS)  -c snclp3a.f
	@if test -f snclp3a.T;then echo "	mv -f snclp3a.T NORMAL/.";	mv -f snclp3a.T NORMAL/. ; fi
	mv -f snclp3a.o NORMAL/snclp3a.o

NORMAL/snscf3a.o: snscf3a.f
	$(FORT) $(FORTOPTS)  -c snscf3a.f
	@if test -f snscf3a.T;then echo "	mv -f snscf3a.T NORMAL/.";	mv -f snscf3a.T NORMAL/. ; fi
	mv -f snscf3a.o NORMAL/snscf3a.o

NORMAL/solcf3a.o: solcf3a.f
	$(FORT) $(FORTOPTS)  -c solcf3a.f
	@if test -f solcf3a.T;then echo "	mv -f solcf3a.T NORMAL/.";	mv -f solcf3a.T NORMAL/. ; fi
	mv -f solcf3a.o NORMAL/solcf3a.o

NORMAL/solsr3a.o: solsr3a.f
	$(FORT) $(FORTOPTS)  -c solsr3a.f
	@if test -f solsr3a.T;then echo "	mv -f solsr3a.T NORMAL/.";	mv -f solsr3a.T NORMAL/. ; fi
	mv -f solsr3a.o NORMAL/solsr3a.o

NORMAL/spin3a.o: spin3a.f
	$(FORT) $(FORTOPTS)  -c spin3a.f
	@if test -f spin3a.T;then echo "	mv -f spin3a.T NORMAL/.";	mv -f spin3a.T NORMAL/. ; fi
	mv -f spin3a.o NORMAL/spin3a.o

NORMAL/splmx3a.o: splmx3a.f
	$(FORT) $(FORTOPTS)  -c splmx3a.f
	@if test -f splmx3a.T;then echo "	mv -f splmx3a.T NORMAL/.";	mv -f splmx3a.T NORMAL/. ; fi
	mv -f splmx3a.o NORMAL/splmx3a.o

NORMAL/ssct3a.o: ssct3a.f
	$(FORT) $(FORTOPTS)  -c ssct3a.f
	@if test -f ssct3a.T;then echo "	mv -f ssct3a.T NORMAL/.";	mv -f ssct3a.T NORMAL/. ; fi
	mv -f ssct3a.o NORMAL/ssct3a.o

NORMAL/sscta3a.o: sscta3a.f
	$(FORT) $(FORTOPTS)  -c sscta3a.f
	@if test -f sscta3a.T;then echo "	mv -f sscta3a.T NORMAL/.";	mv -f sscta3a.T NORMAL/. ; fi
	mv -f sscta3a.o NORMAL/sscta3a.o

NORMAL/stmtf3a.o: stmtf3a.f
	$(FORT) $(FORTOPTS)  -c stmtf3a.f
	@if test -f stmtf3a.T;then echo "	mv -f stmtf3a.T NORMAL/.";	mv -f stmtf3a.T NORMAL/. ; fi
	mv -f stmtf3a.o NORMAL/stmtf3a.o

NORMAL/stsct3a.o: stsct3a.f
	$(FORT) $(FORTOPTS)  -c stsct3a.f
	@if test -f stsct3a.T;then echo "	mv -f stsct3a.T NORMAL/.";	mv -f stsct3a.T NORMAL/. ; fi
	mv -f stsct3a.o NORMAL/stsct3a.o

NORMAL/surfp3a.o: surfp3a.f
	$(FORT) $(FORTOPTS)  -c surfp3a.f
	@if test -f surfp3a.T;then echo "	mv -f surfp3a.T NORMAL/.";	mv -f surfp3a.T NORMAL/. ; fi
	mv -f surfp3a.o NORMAL/surfp3a.o

NORMAL/tcf3a.o: tcf3a.f
	$(FORT) $(FORTOPTS)  -c tcf3a.f
	@if test -f tcf3a.T;then echo "	mv -f tcf3a.T NORMAL/.";	mv -f tcf3a.T NORMAL/. ; fi
	mv -f tcf3a.o NORMAL/tcf3a.o

NORMAL/tcfb3a.o: tcfb3a.f
	$(FORT) $(FORTOPTS)  -c tcfb3a.f
	@if test -f tcfb3a.T;then echo "	mv -f tcfb3a.T NORMAL/.";	mv -f tcfb3a.T NORMAL/. ; fi
	mv -f tcfb3a.o NORMAL/tcfb3a.o

NORMAL/tcfc3a.o: tcfc3a.f
	$(FORT) $(FORTOPTS)  -c tcfc3a.f
	@if test -f tcfc3a.T;then echo "	mv -f tcfc3a.T NORMAL/.";	mv -f tcfc3a.T NORMAL/. ; fi
	mv -f tcfc3a.o NORMAL/tcfc3a.o

NORMAL/tropin1a.o: tropin1a.f
	$(FORT) $(FORTOPTS)  -c tropin1a.f
	@if test -f tropin1a.T;then echo "	mv -f tropin1a.T NORMAL/.";	mv -f tropin1a.T NORMAL/. ; fi
	mv -f tropin1a.o NORMAL/tropin1a.o

NORMAL/trpclm3a.o: trpclm3a.f
	$(FORT) $(FORTOPTS)  -c trpclm3a.f
	@if test -f trpclm3a.T;then echo "	mv -f trpclm3a.T NORMAL/.";	mv -f trpclm3a.T NORMAL/. ; fi
	mv -f trpclm3a.o NORMAL/trpclm3a.o

NORMAL/trpssr3a.o: trpssr3a.f
	$(FORT) $(FORTOPTS)  -c trpssr3a.f
	@if test -f trpssr3a.T;then echo "	mv -f trpssr3a.T NORMAL/.";	mv -f trpssr3a.T NORMAL/. ; fi
	mv -f trpssr3a.o NORMAL/trpssr3a.o

NORMAL/trsfc3a.o: trsfc3a.f
	$(FORT) $(FORTOPTS)  -c trsfc3a.f
	@if test -f trsfc3a.T;then echo "	mv -f trsfc3a.T NORMAL/.";	mv -f trsfc3a.T NORMAL/. ; fi
	mv -f trsfc3a.o NORMAL/trsfc3a.o

NORMAL/tstr3a.o: tstr3a.f
	$(FORT) $(FORTOPTS)  -c tstr3a.f
	@if test -f tstr3a.T;then echo "	mv -f tstr3a.T NORMAL/.";	mv -f tstr3a.T NORMAL/. ; fi
	mv -f tstr3a.o NORMAL/tstr3a.o

NORMAL/twcfrg3a.o: twcfrg3a.f
	$(FORT) $(FORTOPTS)  -c twcfrg3a.f
	@if test -f twcfrg3a.T;then echo "	mv -f twcfrg3a.T NORMAL/.";	mv -f twcfrg3a.T NORMAL/. ; fi
	mv -f twcfrg3a.o NORMAL/twcfrg3a.o

