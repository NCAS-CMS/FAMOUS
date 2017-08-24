# ------------------------ FILE: MAKEFILE ------------------------

POBJS=update.o proc_dir.o utils.o proc_pl.o proc_if.o proc_ibd.o strfns.o
MOBJS=mklib.o utils.o strfns.o
ROBJS=rsort.o utils.o strfns.o
INCLUDE=defns.h structs.h

# Uncomment the following to use alternative string compare functions
STRCMP=-DALT_STRCMP
#STRCMP=

# for the Cray C compiler
 CFLAGS = -T cray-ymp -O -Aa -DSAVE_LOG -DHP $(STRCMP)
 CC=cc
 LINK=cc
 LFLAGS= -T cray-ymp

all: nupdate mklib rsort

nupdate : $(POBJS) $(INCLUDE)
	$(LINK) $(LFLAGS) $(POBJS) -o nupdate

mklib : $(MOBJS) $(INCLUDE)
	$(LINK) $(LFLAGS) $(MOBJS) -o mklib

rsort : $(ROBJS) $(INCLUDE)
	$(LINK) $(LFLAGS) $(ROBJS) -o rsort

utils.o : utils.c $(INCLUDE)

update.o : update.c $(INCLUDE)

proc_pl.o : proc_pl.c $(INCLUDE)

proc_dir.o : proc_dir.c $(INCLUDE)

proc_if.o : proc_if.c $(INCLUDE)

proc_ibd.o : proc_ibd.c $(INCLUDE)

mklib.o : mklib.c $(INCLUDE)

strfns.o : strfns.c

clean :
	rm nupdate mklib rsort
	rm *.o
	rm core


