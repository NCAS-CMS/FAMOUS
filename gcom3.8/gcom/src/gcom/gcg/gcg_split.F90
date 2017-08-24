! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gcg_prolog.h"    

SUBROUTINE GCG_SPLIT(ME, NPROC, COLOR, ISTAT, GID)
!     ******************************************************************
!     * Purpose:
!     *
!     *  GC groups routine. Splits processors in disjoint groups based 
!     *  on a color attribute.  
!     *
!     * Input:
!     *  ME      - my node ID
!     *  NPROC   - number of nodes
!     *  COLOR   - color attribute for this node
!     *
!     * Output:
!     *  GID     - group identifier
!     *  ISTAT   - status of split 0 is OK
!     *
!     * NOTES:
!     *  GCG_SPLIT() is a synchronizing call.
!     *
!     ******************************************************************

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"

INTEGER (KIND=GC_INT_KIND) :: ME, NPROC, COLOR, GID, ISTAT
LOGICAL (KIND=GC_LOG_KIND), SAVE :: GINITED = .FALSE.

ISTAT = GC__OK
IF (.NOT. GINITED) CALL GCG__INIT_INTERNALS(GINITED)

#if defined(MPI_SRC)
CALL MPL_COMM_SPLIT(GC__MY_MPI_COMM_WORLD, COLOR,                 &
                    0_GC_INT_KIND, GID, ISTAT)
#endif

RETURN
END
