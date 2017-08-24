! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gcg_prolog.h"    

SUBROUTINE GCG_RSUM (LEN, GID, ISTAT, RSUM)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Calculate the real sum across all processors of a group and
!     *  distribute the result to all members of the group.
!     *
!     * Input:
!     *  LEN     - number of elements in message
!     *  GID     - processor group ID
!     *  RSUM    - array with elements to be added up across the nodes
!     *
!     * Output:
!     *  RSUM    - array containing the sums across the nodes
!     *  ISTAT   - status of rsum. 0 is OK (MPI_SRC only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

USE MPL, ONLY :                                                   &
    MPL_REAL,                                                     &
    MPL_SUM

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD

USE GCOM_MOD, ONLY :                                              &
    GC_FORCE_BITREP,                                              &
    GC_ON
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"
#include "gcg_constants.h"

INTEGER (KIND=GC_INT_KIND) :: LEN, GID, ISTAT
INTEGER (KIND=GC_INT_KIND) :: OPT
REAL (KIND=GC_REAL_KIND)   :: RSUM(LEN)

REAL (KIND=GC_REAL_KIND)   :: REDUCE_DATA_WRK(LEN)

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: IGID
#endif

#if defined(GCOM_TIMER)
#include "gc_timer.h"
#define THIS_TIMER TIMET_COLL
#define THIS_LENGTH LEN*GC__RSIZE
#endif
INTEGER (KIND=GC_INT_KIND) :: I

#include "gc_start_timer.h"

ISTAT = GC__OK

#if defined(MPI_SRC)
CALL GC_GETOPT(GC_FORCE_BITREP, OPT, ISTAT)
IF (OPT == GC_ON) THEN
  CALL GCG_RSUMR(LEN, GID, ISTAT, RSUM)
ELSE
  IF (GID  ==  GCG__ALLGROUP) THEN
     IGID = GC__MY_MPI_COMM_WORLD
  ELSE
     IGID = GID
  ENDIF
  DO I = 1,LEN
     REDUCE_DATA_WRK(I) = RSUM(I)
  ENDDO
  CALL MPL_ALLREDUCE(REDUCE_DATA_WRK, RSUM, LEN, MPL_REAL,        &
       MPL_SUM, IGID, ISTAT)
END IF
#endif


#include "gc_end_timer.h"

RETURN
END
