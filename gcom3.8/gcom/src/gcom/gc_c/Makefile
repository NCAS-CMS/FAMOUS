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

${GC_LIB} : gc__abort.o


clean:
	@${ECHO} \*\*\* Deleting files from ${UPACKAGE}-${COMPONENT}
	${DELETE} gc__abort.o
# DO NOT DELETE
