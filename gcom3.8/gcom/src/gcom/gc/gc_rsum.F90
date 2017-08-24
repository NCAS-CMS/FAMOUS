! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"    

SUBROUTINE GC_RSUM (LEN, NPROC, ISTAT, SSUM)
!     ******************************************************************
!     * Purpose:
!     *  Calculate the real sum across all processors and distribute
!     *  the result to all the processors.
!     *
!     * Input:
!     *  LEN     - number of elements in message
!     *  NPROC   - number of processors
!     *  SSUM    - array with elements to be added up across the nodes
!     *
!     * Output:
!     *  SSUM    - array containing the sums across the nodes
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

INTEGER (KIND=GC_INT_KIND) :: LEN, NPROC, ISTAT
INTEGER (KIND=GC_INT_KIND) :: OPT
REAL (KIND=GC_REAL_KIND)   :: SSUM(LEN)
REAL (KIND=GC_REAL_KIND)   :: REDUCE_DATA_WRK(LEN)

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: I
#endif

#if defined(GCOM_TIMER)
#include "gc_timer.h"
#define THIS_TIMER TIMET_COLL
#define THIS_LENGTH LEN*GC__RSIZE
#endif

#include "gc_start_timer.h"

#if defined(MPI_SRC)
CALL GC_GETOPT(GC_FORCE_BITREP, OPT, ISTAT)
IF (OPT == GC_ON) THEN
  CALL GC_RSUMR(LEN, NPROC, ISTAT, SSUM)
ELSE
  DO I = 1,LEN
     REDUCE_DATA_WRK(I) = SSUM(I)
  END DO
  CALL MPL_ALLREDUCE(REDUCE_DATA_WRK, SSUM, LEN, MPL_REAL,        &
       MPL_SUM, GC__MY_MPI_COMM_WORLD, ISTAT)
END IF
#endif

#if defined(SERIAL_SRC)
ISTAT = GC__OK
#endif

#include "gc_end_timer.h"

RETURN
END


