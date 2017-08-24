! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"    

SUBROUTINE GC_RSUMR (LEN, NPROC, ISTAT, SSUM)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Calculate in a reproducible way the real sum across all 
!     *  processors and distribute the result to all the processors.
!     *
!     * Input:
!     *  LEN     - number of elements in message
!     *  NPROC   - number of processors
!     *  SSUM    - array with elements to be added up across the nodes
!     *
!     * Output:
!     *  SSUM    - array containing the sums across the nodes
!     *  ISTAT   - status of rsumr. 0 is OK (MPI_SRC only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

USE MPL, ONLY :                                                   &
    MPL_STATUS_SIZE,                                              &
    MPL_REAL

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_mtags.h"
#include "gc_constants.h"

INTEGER (KIND=GC_INT_KIND) :: LEN, NPROC, ISTAT
REAL (KIND=GC_REAL_KIND)   :: SSUM(LEN)
REAL (KIND=GC_REAL_KIND)   :: REDUCE_DATA_WRK(LEN,NPROC)

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: I, L, ME
#include "gc_functions.h"
INTEGER (KIND=GC_INT_KIND) :: STATUS(MPL_STATUS_SIZE)
INTEGER (KIND=GC_INT_KIND) :: REQUESTS(NPROC)
#endif

#if defined(GCOM_TIMER)
#include "gc_timer.h"
#define THIS_TIMER TIMET_COLL
#define THIS_LENGTH LEN*GC__RSIZE
#endif

#include "gc_start_timer.h"
#if defined(MPI_SRC)
ME = GC_ME()
IF (ME  /=  GC__IONODE) THEN
  CALL MPL_SEND(SSUM, LEN, MPL_REAL, GC__IONODE,                 &
                GCID__RSUM0, GC__MY_MPI_COMM_WORLD, ISTAT)
ELSE
  ! Issue all receives - non-blocking
  DO I = 1, NPROC-1
    CALL MPL_IRECV(REDUCE_DATA_WRK(:,I), LEN, MPL_REAL,          &
                   I, GCID__RSUM0, GC__MY_MPI_COMM_WORLD,        &
                   REQUESTS(I), ISTAT)
  END DO

  ! Wait for receives to complete in order and then add in contribution
  DO I = 1, NPROC-1
    CALL MPL_WAIT(REQUESTS(I), STATUS, ISTAT)

    DO L = 1,LEN
      SSUM(L) = SSUM(L) + REDUCE_DATA_WRK(L,I)
    END DO
  END DO
END IF
CALL MPL_BCAST(SSUM, LEN, MPL_REAL, GC__IONODE,                   &
               GC__MY_MPI_COMM_WORLD, ISTAT)
#endif

#if defined(SERIAL_SRC)
ISTAT = GC__OK
#endif

#include "gc_end_timer.h"

RETURN
END
