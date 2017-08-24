#                      Top Level GCOM Makefile

# *****************************COPYRIGHT*******************************
# (c) CROWN COPYRIGHT 2001, Met Office, All Rights Reserved. 
# Please refer to Copyright file in top level GCOM directory
#                 for further details
# *****************************COPYRIGHT*******************************

# "make_inc/user_inc.mk" contains all user-editable bits and pieces
# There should be no need to edit this main Makefile itself.

include make_inc/user_inc.mk
include make_inc/sys_inc.mk

all : gcom
library : gcom

gcom :
	@${MAKE} do_makelib \
        GC_VERSION="3.8" \
        GC_TOPDIR=`pwd` \
        GC_BUILD_DATE="`date`"

do_makelib:
	@${ECHO} \*\*\* "Creating ${GC_DESCRIP} version of ${UPACKAGE} (version ${GC_VERSION})"
	@${ECHO} \*\*\* Library name is ${GC_LIB}
	@${ECHO} 
	@${ECHO} >> ${GC_LOG}
	@${ECHO} \*\*\* make started at `date` >> ${GC_LOG}
	@cd src/gcom/mpl ; ${MAKE} GC_VERSION=${GC_VERSION} \
                         GC_TOPDIR=${GC_TOPDIR} \
                         GC_BUILD_DATE="${GC_BUILD_DATE}"
	@cd src/gcom/gc ; ${MAKE} GC_VERSION=${GC_VERSION} \
                         GC_TOPDIR=${GC_TOPDIR} \
                         GC_BUILD_DATE="${GC_BUILD_DATE}"
	@cd src/gcom/gcg ; ${MAKE} GC_VERSION=${GC_VERSION} \
                          GC_TOPDIR=${GC_TOPDIR} \
                          GC_BUILD_DATE="${GC_BUILD_DATE}"
	@cd src/gcom/gc_c ; ${MAKE}

clean :
	@${MAKE} do_clean ${GC_EXTRA_BITS_CLEAN} \
          GC_VERSION="3.8" \
          GC_TOPDIR=`pwd`

do_clean :
	@${ECHO} \*\*\* "Cleaning ${UPACKAGE} (version ${GC_VERSION})"
	@${ECHO}
	@cd  src/gcom/gc ; ${MAKE} clean GC_VERSION=${GC_VERSION} GC_TOPDIR=${GC_TOPDIR}
	@cd  src/gcom/gcg ; ${MAKE} clean GC_VERSION=${GC_VERSION} GC_TOPDIR=${GC_TOPDIR}
	@cd src/gcom/mpl ; ${MAKE} clean GC_VERSION=${GC_VERSION} GC_TOPDIR=${GC_TOPDIR}
	@cd  src/gcom/gc_c ; ${MAKE} clean GC_VERSION=${GC_VERSION} GC_TOPDIR=${GC_TOPDIR}
	@${ECHO} \*\*\*\* Library was deleted at `date` >> ${GC_LOG}
	@${ECHO}
	@${ECHO} \*\*\* Deleting library `${BASENAME} ${GC_LIB}`
	${DELETE} ${GC_LIB}

