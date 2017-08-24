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
NORMAL/coex1a.o NORMAL/exptenc1.o NORMAL/extts1a.o NORMAL/grbtim1a.o\
NORMAL/grbwrt1a.o NORMAL/initpp1a.o NORMAL/mosgrd1a.o NORMAL/muspac1a.o\
NORMAL/pp2griba.o NORMAL/ppfile1a.o NORMAL/pphead1a.o NORMAL/setlst1a.o\
NORMAL/setpls1a.o NORMAL/spacia1a.o NORMAL/stacum1a.o NORMAL/stcolm1a.o\
NORMAL/stcomg1a.o NORMAL/stextc1a.o NORMAL/stfldm1a.o NORMAL/stglom1a.o\
NORMAL/stlev1a.o NORMAL/stmax1a.o NORMAL/stmerm1a.o NORMAL/stmin1a.o\
NORMAL/stwork1a.o NORMAL/stzonm1a.o NORMAL/temprl1a.o

all: $(OBJS)


NORMAL/coex1a.o: coex1a.f
	$(FORT) $(FORTOPTS)  -c coex1a.f
	@if test -f coex1a.T;then echo "	mv -f coex1a.T NORMAL/.";	mv -f coex1a.T NORMAL/. ; fi
	mv -f coex1a.o NORMAL/coex1a.o

NORMAL/exptenc1.o: exptenc1.f
	$(FORT) $(FORTOPTS)  -c exptenc1.f
	@if test -f exptenc1.T;then echo "	mv -f exptenc1.T NORMAL/.";	mv -f exptenc1.T NORMAL/. ; fi
	mv -f exptenc1.o NORMAL/exptenc1.o

NORMAL/extts1a.o: extts1a.f
	$(FORT) $(FORTOPTS)  -c extts1a.f
	@if test -f extts1a.T;then echo "	mv -f extts1a.T NORMAL/.";	mv -f extts1a.T NORMAL/. ; fi
	mv -f extts1a.o NORMAL/extts1a.o

NORMAL/grbtim1a.o: grbtim1a.f
	$(FORT) $(FORTOPTS)  -c grbtim1a.f
	@if test -f grbtim1a.T;then echo "	mv -f grbtim1a.T NORMAL/.";	mv -f grbtim1a.T NORMAL/. ; fi
	mv -f grbtim1a.o NORMAL/grbtim1a.o

NORMAL/grbwrt1a.o: grbwrt1a.f
	$(FORT) $(FORTOPTS)  -c grbwrt1a.f
	@if test -f grbwrt1a.T;then echo "	mv -f grbwrt1a.T NORMAL/.";	mv -f grbwrt1a.T NORMAL/. ; fi
	mv -f grbwrt1a.o NORMAL/grbwrt1a.o

NORMAL/initpp1a.o: initpp1a.f
	$(FORT) $(FORTOPTS)  -c initpp1a.f
	@if test -f initpp1a.T;then echo "	mv -f initpp1a.T NORMAL/.";	mv -f initpp1a.T NORMAL/. ; fi
	mv -f initpp1a.o NORMAL/initpp1a.o

NORMAL/mosgrd1a.o: mosgrd1a.f
	$(FORT) $(FORTOPTS)  -c mosgrd1a.f
	@if test -f mosgrd1a.T;then echo "	mv -f mosgrd1a.T NORMAL/.";	mv -f mosgrd1a.T NORMAL/. ; fi
	mv -f mosgrd1a.o NORMAL/mosgrd1a.o

NORMAL/muspac1a.o: muspac1a.f
	$(FORT) $(FORTOPTS)  -c muspac1a.f
	@if test -f muspac1a.T;then echo "	mv -f muspac1a.T NORMAL/.";	mv -f muspac1a.T NORMAL/. ; fi
	mv -f muspac1a.o NORMAL/muspac1a.o

NORMAL/pp2griba.o: pp2griba.f
	$(FORT) $(FORTOPTS)  -c pp2griba.f
	@if test -f pp2griba.T;then echo "	mv -f pp2griba.T NORMAL/.";	mv -f pp2griba.T NORMAL/. ; fi
	mv -f pp2griba.o NORMAL/pp2griba.o

NORMAL/ppfile1a.o: ppfile1a.f
	$(FORT) $(FORTOPTS)  -c ppfile1a.f
	@if test -f ppfile1a.T;then echo "	mv -f ppfile1a.T NORMAL/.";	mv -f ppfile1a.T NORMAL/. ; fi
	mv -f ppfile1a.o NORMAL/ppfile1a.o

NORMAL/pphead1a.o: pphead1a.f
	$(FORT) $(FORTOPTS)  -c pphead1a.f
	@if test -f pphead1a.T;then echo "	mv -f pphead1a.T NORMAL/.";	mv -f pphead1a.T NORMAL/. ; fi
	mv -f pphead1a.o NORMAL/pphead1a.o

NORMAL/setlst1a.o: setlst1a.f
	$(FORT) $(FORTOPTS)  -c setlst1a.f
	@if test -f setlst1a.T;then echo "	mv -f setlst1a.T NORMAL/.";	mv -f setlst1a.T NORMAL/. ; fi
	mv -f setlst1a.o NORMAL/setlst1a.o

NORMAL/setpls1a.o: setpls1a.f
	$(FORT) $(FORTOPTS)  -c setpls1a.f
	@if test -f setpls1a.T;then echo "	mv -f setpls1a.T NORMAL/.";	mv -f setpls1a.T NORMAL/. ; fi
	mv -f setpls1a.o NORMAL/setpls1a.o

NORMAL/spacia1a.o: spacia1a.f
	$(FORT) $(FORTOPTS)  -c spacia1a.f
	@if test -f spacia1a.T;then echo "	mv -f spacia1a.T NORMAL/.";	mv -f spacia1a.T NORMAL/. ; fi
	mv -f spacia1a.o NORMAL/spacia1a.o

NORMAL/stacum1a.o: stacum1a.f
	$(FORT) $(FORTOPTS)  -c stacum1a.f
	@if test -f stacum1a.T;then echo "	mv -f stacum1a.T NORMAL/.";	mv -f stacum1a.T NORMAL/. ; fi
	mv -f stacum1a.o NORMAL/stacum1a.o

NORMAL/stcolm1a.o: stcolm1a.f
	$(FORT) $(FORTOPTS)  -c stcolm1a.f
	@if test -f stcolm1a.T;then echo "	mv -f stcolm1a.T NORMAL/.";	mv -f stcolm1a.T NORMAL/. ; fi
	mv -f stcolm1a.o NORMAL/stcolm1a.o

NORMAL/stcomg1a.o: stcomg1a.f
	$(FORT) $(FORTOPTS)  -c stcomg1a.f
	@if test -f stcomg1a.T;then echo "	mv -f stcomg1a.T NORMAL/.";	mv -f stcomg1a.T NORMAL/. ; fi
	mv -f stcomg1a.o NORMAL/stcomg1a.o

NORMAL/stextc1a.o: stextc1a.f
	$(FORT) $(FORTOPTS)  -c stextc1a.f
	@if test -f stextc1a.T;then echo "	mv -f stextc1a.T NORMAL/.";	mv -f stextc1a.T NORMAL/. ; fi
	mv -f stextc1a.o NORMAL/stextc1a.o

NORMAL/stfldm1a.o: stfldm1a.f
	$(FORT) $(FORTOPTS)  -c stfldm1a.f
	@if test -f stfldm1a.T;then echo "	mv -f stfldm1a.T NORMAL/.";	mv -f stfldm1a.T NORMAL/. ; fi
	mv -f stfldm1a.o NORMAL/stfldm1a.o

NORMAL/stglom1a.o: stglom1a.f
	$(FORT) $(FORTOPTS)  -c stglom1a.f
	@if test -f stglom1a.T;then echo "	mv -f stglom1a.T NORMAL/.";	mv -f stglom1a.T NORMAL/. ; fi
	mv -f stglom1a.o NORMAL/stglom1a.o

NORMAL/stlev1a.o: stlev1a.f
	$(FORT) $(FORTOPTS)  -c stlev1a.f
	@if test -f stlev1a.T;then echo "	mv -f stlev1a.T NORMAL/.";	mv -f stlev1a.T NORMAL/. ; fi
	mv -f stlev1a.o NORMAL/stlev1a.o

NORMAL/stmax1a.o: stmax1a.f
	$(FORT) $(FORTOPTS)  -c stmax1a.f
	@if test -f stmax1a.T;then echo "	mv -f stmax1a.T NORMAL/.";	mv -f stmax1a.T NORMAL/. ; fi
	mv -f stmax1a.o NORMAL/stmax1a.o

NORMAL/stmerm1a.o: stmerm1a.f
	$(FORT) $(FORTOPTS)  -c stmerm1a.f
	@if test -f stmerm1a.T;then echo "	mv -f stmerm1a.T NORMAL/.";	mv -f stmerm1a.T NORMAL/. ; fi
	mv -f stmerm1a.o NORMAL/stmerm1a.o

NORMAL/stmin1a.o: stmin1a.f
	$(FORT) $(FORTOPTS)  -c stmin1a.f
	@if test -f stmin1a.T;then echo "	mv -f stmin1a.T NORMAL/.";	mv -f stmin1a.T NORMAL/. ; fi
	mv -f stmin1a.o NORMAL/stmin1a.o

NORMAL/stwork1a.o: stwork1a.f
	$(FORT) $(FORTOPTS)  -c stwork1a.f
	@if test -f stwork1a.T;then echo "	mv -f stwork1a.T NORMAL/.";	mv -f stwork1a.T NORMAL/. ; fi
	mv -f stwork1a.o NORMAL/stwork1a.o

NORMAL/stzonm1a.o: stzonm1a.f
	$(FORT) $(FORTOPTS)  -c stzonm1a.f
	@if test -f stzonm1a.T;then echo "	mv -f stzonm1a.T NORMAL/.";	mv -f stzonm1a.T NORMAL/. ; fi
	mv -f stzonm1a.o NORMAL/stzonm1a.o

NORMAL/temprl1a.o: temprl1a.f
	$(FORT) $(FORTOPTS)  -c temprl1a.f
	@if test -f temprl1a.T;then echo "	mv -f temprl1a.T NORMAL/.";	mv -f temprl1a.T NORMAL/. ; fi
	mv -f temprl1a.o NORMAL/temprl1a.o

