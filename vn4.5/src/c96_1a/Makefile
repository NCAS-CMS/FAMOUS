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

BUILDTYPES =  A1G1G3GrGw A1G1G3GrGwO1 A1G1GrGw G1G3GrGwO1 G1GrGwO1

all:
	@for buildtype in $(BUILDTYPES) ;\
	do ( cd $$buildtype ; pwd ; make ) ; done

install: all
	@for buildtype in $(BUILDTYPES) ;\
	do ( cd $$buildtype ; pwd ; make install ) ; done
