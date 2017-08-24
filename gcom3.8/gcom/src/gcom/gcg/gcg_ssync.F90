! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gcg_prolog.h"    

SUBROUTINE GCG_SSYNC (GID, ISTAT)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Group synchronisation.
!     *
!     * Input:
!     *  GID     - processor group ID
!     *
!     * Output:
!     *  ISTAT   - status of rsum. 0 is OK (MPI_SRC only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

IMPLICIT NONE

#include "gc_kinds.h"

INTEGER (KIND=GC_INT_KIND) :: GID, ISTAT

#if defined(GCOM_TIMER)
#include "gc_timer.h"
#define THIS_TIMER TIMET_BARRIER
#define THIS_LENGTH 0
#endif


INTEGER (KIND=GC_INT_KIND) :: L, ME, NPROC, IGID, ILOC, GRANK
#include "gc_functions.h"
INTEGER (KIND=GC_INT_KIND) :: I


ISTAT = 0

IF (GID  ==  0) THEN
   NPROC = GC_NPROC()
   CALL GC_SSYNC(NPROC, ISTAT)
   RETURN
ENDIF

RETURN
END
