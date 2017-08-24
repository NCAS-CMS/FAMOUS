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
NORMAL/abort1a.o NORMAL/chkloo1a.o NORMAL/exp32b1a.o NORMAL/g_ch_l1a.o\
NORMAL/getdim1a.o NORMAL/ieeepk1a.o NORMAL/initfh1a.o NORMAL/ioerr1a.o\
NORMAL/newpac1a.o NORMAL/poserr1a.o NORMAL/pr_f_c1a.o NORMAL/prfixh1a.o\
NORMAL/prifld1a.o NORMAL/prihda1a.o NORMAL/prlfld1a.o NORMAL/prlook1a.o\
NORMAL/prrfld1a.o NORMAL/prrhda1a.o NORMAL/readdm1a.o NORMAL/readfh1a.o\
NORMAL/readfl1a.o NORMAL/readhe1a.o NORMAL/sdfadr1a.o NORMAL/writdm1a.o\
NORMAL/writfl1a.o NORMAL/writhe1a.o

all: $(OBJS)


NORMAL/abort1a.o: abort1a.f
	$(FORT) $(FORTOPTS)  -c abort1a.f
	@if test -f abort1a.T;then echo "	mv -f abort1a.T NORMAL/.";	mv -f abort1a.T NORMAL/. ; fi
	mv -f abort1a.o NORMAL/abort1a.o

NORMAL/chkloo1a.o: chkloo1a.f
	$(FORT) $(FORTOPTS)  -c chkloo1a.f
	@if test -f chkloo1a.T;then echo "	mv -f chkloo1a.T NORMAL/.";	mv -f chkloo1a.T NORMAL/. ; fi
	mv -f chkloo1a.o NORMAL/chkloo1a.o

NORMAL/exp32b1a.o: exp32b1a.f
	$(FORT) $(FORTOPTS)  -c exp32b1a.f
	@if test -f exp32b1a.T;then echo "	mv -f exp32b1a.T NORMAL/.";	mv -f exp32b1a.T NORMAL/. ; fi
	mv -f exp32b1a.o NORMAL/exp32b1a.o

NORMAL/g_ch_l1a.o: g_ch_l1a.f
	$(FORT) $(FORTOPTS)  -c g_ch_l1a.f
	@if test -f g_ch_l1a.T;then echo "	mv -f g_ch_l1a.T NORMAL/.";	mv -f g_ch_l1a.T NORMAL/. ; fi
	mv -f g_ch_l1a.o NORMAL/g_ch_l1a.o

NORMAL/getdim1a.o: getdim1a.f
	$(FORT) $(FORTOPTS)  -c getdim1a.f
	@if test -f getdim1a.T;then echo "	mv -f getdim1a.T NORMAL/.";	mv -f getdim1a.T NORMAL/. ; fi
	mv -f getdim1a.o NORMAL/getdim1a.o

NORMAL/ieeepk1a.o: ieeepk1a.f
	$(FORT) $(FORTOPTS)  -c ieeepk1a.f
	@if test -f ieeepk1a.T;then echo "	mv -f ieeepk1a.T NORMAL/.";	mv -f ieeepk1a.T NORMAL/. ; fi
	mv -f ieeepk1a.o NORMAL/ieeepk1a.o

NORMAL/initfh1a.o: initfh1a.f
	$(FORT) $(FORTOPTS)  -c initfh1a.f
	@if test -f initfh1a.T;then echo "	mv -f initfh1a.T NORMAL/.";	mv -f initfh1a.T NORMAL/. ; fi
	mv -f initfh1a.o NORMAL/initfh1a.o

NORMAL/ioerr1a.o: ioerr1a.f
	$(FORT) $(FORTOPTS)  -c ioerr1a.f
	@if test -f ioerr1a.T;then echo "	mv -f ioerr1a.T NORMAL/.";	mv -f ioerr1a.T NORMAL/. ; fi
	mv -f ioerr1a.o NORMAL/ioerr1a.o

NORMAL/newpac1a.o: newpac1a.f
	$(FORT) $(FORTOPTS)  -c newpac1a.f
	@if test -f newpac1a.T;then echo "	mv -f newpac1a.T NORMAL/.";	mv -f newpac1a.T NORMAL/. ; fi
	mv -f newpac1a.o NORMAL/newpac1a.o

NORMAL/poserr1a.o: poserr1a.f
	$(FORT) $(FORTOPTS)  -c poserr1a.f
	@if test -f poserr1a.T;then echo "	mv -f poserr1a.T NORMAL/.";	mv -f poserr1a.T NORMAL/. ; fi
	mv -f poserr1a.o NORMAL/poserr1a.o

NORMAL/pr_f_c1a.o: pr_f_c1a.f
	$(FORT) $(FORTOPTS)  -c pr_f_c1a.f
	@if test -f pr_f_c1a.T;then echo "	mv -f pr_f_c1a.T NORMAL/.";	mv -f pr_f_c1a.T NORMAL/. ; fi
	mv -f pr_f_c1a.o NORMAL/pr_f_c1a.o

NORMAL/prfixh1a.o: prfixh1a.f
	$(FORT) $(FORTOPTS)  -c prfixh1a.f
	@if test -f prfixh1a.T;then echo "	mv -f prfixh1a.T NORMAL/.";	mv -f prfixh1a.T NORMAL/. ; fi
	mv -f prfixh1a.o NORMAL/prfixh1a.o

NORMAL/prifld1a.o: prifld1a.f
	$(FORT) $(FORTOPTS)  -c prifld1a.f
	@if test -f prifld1a.T;then echo "	mv -f prifld1a.T NORMAL/.";	mv -f prifld1a.T NORMAL/. ; fi
	mv -f prifld1a.o NORMAL/prifld1a.o

NORMAL/prihda1a.o: prihda1a.f
	$(FORT) $(FORTOPTS)  -c prihda1a.f
	@if test -f prihda1a.T;then echo "	mv -f prihda1a.T NORMAL/.";	mv -f prihda1a.T NORMAL/. ; fi
	mv -f prihda1a.o NORMAL/prihda1a.o

NORMAL/prlfld1a.o: prlfld1a.f
	$(FORT) $(FORTOPTS)  -c prlfld1a.f
	@if test -f prlfld1a.T;then echo "	mv -f prlfld1a.T NORMAL/.";	mv -f prlfld1a.T NORMAL/. ; fi
	mv -f prlfld1a.o NORMAL/prlfld1a.o

NORMAL/prlook1a.o: prlook1a.f
	$(FORT) $(FORTOPTS)  -c prlook1a.f
	@if test -f prlook1a.T;then echo "	mv -f prlook1a.T NORMAL/.";	mv -f prlook1a.T NORMAL/. ; fi
	mv -f prlook1a.o NORMAL/prlook1a.o

NORMAL/prrfld1a.o: prrfld1a.f
	$(FORT) $(FORTOPTS)  -c prrfld1a.f
	@if test -f prrfld1a.T;then echo "	mv -f prrfld1a.T NORMAL/.";	mv -f prrfld1a.T NORMAL/. ; fi
	mv -f prrfld1a.o NORMAL/prrfld1a.o

NORMAL/prrhda1a.o: prrhda1a.f
	$(FORT) $(FORTOPTS)  -c prrhda1a.f
	@if test -f prrhda1a.T;then echo "	mv -f prrhda1a.T NORMAL/.";	mv -f prrhda1a.T NORMAL/. ; fi
	mv -f prrhda1a.o NORMAL/prrhda1a.o

NORMAL/readdm1a.o: readdm1a.f
	$(FORT) $(FORTOPTS)  -c readdm1a.f
	@if test -f readdm1a.T;then echo "	mv -f readdm1a.T NORMAL/.";	mv -f readdm1a.T NORMAL/. ; fi
	mv -f readdm1a.o NORMAL/readdm1a.o

NORMAL/readfh1a.o: readfh1a.f
	$(FORT) $(FORTOPTS)  -c readfh1a.f
	@if test -f readfh1a.T;then echo "	mv -f readfh1a.T NORMAL/.";	mv -f readfh1a.T NORMAL/. ; fi
	mv -f readfh1a.o NORMAL/readfh1a.o

NORMAL/readfl1a.o: readfl1a.f
	$(FORT) $(FORTOPTS)  -c readfl1a.f
	@if test -f readfl1a.T;then echo "	mv -f readfl1a.T NORMAL/.";	mv -f readfl1a.T NORMAL/. ; fi
	mv -f readfl1a.o NORMAL/readfl1a.o

NORMAL/readhe1a.o: readhe1a.f
	$(FORT) $(FORTOPTS)  -c readhe1a.f
	@if test -f readhe1a.T;then echo "	mv -f readhe1a.T NORMAL/.";	mv -f readhe1a.T NORMAL/. ; fi
	mv -f readhe1a.o NORMAL/readhe1a.o

NORMAL/sdfadr1a.o: sdfadr1a.f
	$(FORT) $(FORTOPTS)  -c sdfadr1a.f
	@if test -f sdfadr1a.T;then echo "	mv -f sdfadr1a.T NORMAL/.";	mv -f sdfadr1a.T NORMAL/. ; fi
	mv -f sdfadr1a.o NORMAL/sdfadr1a.o

NORMAL/writdm1a.o: writdm1a.f
	$(FORT) $(FORTOPTS)  -c writdm1a.f
	@if test -f writdm1a.T;then echo "	mv -f writdm1a.T NORMAL/.";	mv -f writdm1a.T NORMAL/. ; fi
	mv -f writdm1a.o NORMAL/writdm1a.o

NORMAL/writfl1a.o: writfl1a.f
	$(FORT) $(FORTOPTS)  -c writfl1a.f
	@if test -f writfl1a.T;then echo "	mv -f writfl1a.T NORMAL/.";	mv -f writfl1a.T NORMAL/. ; fi
	mv -f writfl1a.o NORMAL/writfl1a.o

NORMAL/writhe1a.o: writhe1a.f
	$(FORT) $(FORTOPTS)  -c writhe1a.f
	@if test -f writhe1a.T;then echo "	mv -f writhe1a.T NORMAL/.";	mv -f writhe1a.T NORMAL/. ; fi
	mv -f writhe1a.o NORMAL/writhe1a.o

