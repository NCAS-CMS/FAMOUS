FORT=$(FCOM_CMD)
CC=$(CCOM_CMD)

FORTOPTS=$(FCOM_OPTS)
CCOPTS=$(CCOM_OPTS)

C_CODE=codec
F_CODE=codef
EXE=fc_test.exe
MAKEFILE=Makefile
PACKAGE=fc_test.tar
PACKAGE_LIST=$(MAKEFILE) $(F_CODE).f $(C_CODE).c

all:$(EXE)

clean:
	rm -f $(F_CODE).o $(C_CODE).o $(EXE) $(PACKAGE)

package: $(PACKAGE)

$(PACKAGE) : $(PACKAGE_LIST)
	tar cf $(PACKAGE) $(PACKAGE_LIST)

$(C_CODE).o:$(C_CODE).c
	$(CC) $(CCOPTS) -c $(C_CODE).c

$(EXE):$(F_CODE).f $(C_CODE).o
	$(FORT) $(FORTOPTS) $(F_CODE).f $(C_CODE).o -o $(EXE)
	./$(EXE)
