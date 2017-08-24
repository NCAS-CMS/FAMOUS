#                         GC - GCOM Makefile

# *****************************COPYRIGHT*******************************
# (c) CROWN COPYRIGHT 2001, Met Office, All Rights Reserved. 
# Please refer to Copyright file in top level GCOM directory
#                 for further details
# *****************************COPYRIGHT*******************************

# "user_inc.mk" contains all the macro definitions and anything user-editable
# There should be no need to edit this Makefile itself.


include ${GC_TOPDIR}/make_inc/user_inc.mk
include ${GC_TOPDIR}/make_inc/sys_inc.mk
include ${GC_TOPDIR}/make_inc/files.mk
include ${GC_TOPDIR}/make_inc/rules.mk

COMPONENT=GC

all : 
	@${ECHO} "========================================================="
	@${ECHO} \*\*\* Compiling ${UPACKAGE}-${COMPONENT}
	@${ECHO}
	@${MAKE} gcom-gc GC_VERSION=${GC_VERSION} \
                         GC_TOPDIR=${GC_TOPDIR} \
                         GC_BUILD_DATE="${GC_BUILD_DATE}"
	@${ECHO}
	@${ECHO} "========================================================="
	@${ECHO}

gcom-gc : ${GC_LIB}

${GC_LIB} : ${GC_OBJ} 

include ${GC_TOPDIR}/make_inc/gc_depend.mk

clean:
	@${ECHO} \*\*\* Deleting files from ${UPACKAGE}-${COMPONENT}
	${DELETE} ${GC_OBJ} *.f90 *.mod
# DO NOT DELETE
