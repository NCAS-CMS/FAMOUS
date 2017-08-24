! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gcg_prolog.h"    

FUNCTION GCG_ME (GID)
!     ******************************************************************
!     * Purpose:
!     *  
!     *  Return the rank of my node in this group, in the 
!     *  range 0...groupsize(GID)-1 or -1 of invalid GID / not a member.
!     *
!     * Input:
!     *  GID     - processor group ID 
!     *
!     * Output:
!     *
!     * NOTES:       
!     * 
!     ******************************************************************

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gcg_constants.h"

INTEGER (KIND=GC_INT_KIND) :: GCG_ME
INTEGER (KIND=GC_INT_KIND) :: GID

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: ISTAT, IGID, GRANK
#endif

#if defined(MPI_SRC)
IF (GID  ==  GCG__ALLGROUP) THEN
   IGID = GC__MY_MPI_COMM_WORLD
ELSE
   IGID = GID
ENDIF
CALL MPL_COMM_RANK(IGID, GRANK, ISTAT)
GCG_ME = GRANK
#endif

#if defined(SERIAL_SRC)
GCG_ME = 0
#endif

RETURN
END
