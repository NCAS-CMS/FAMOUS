! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"    

SUBROUTINE GC_RMIN (LEN, NPROC, ISTAT, SMIN)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Finds the real minimum across all processors and distribute
!     *  the result to all the processors.
!     *
!     * Input:
!     *  LEN     - number of elements in message
!     *  NPROC   - number of processors
!     *  SMIN    - array with elements of which the elementwise minimum
!     *            across the nodes is to be found
!     *
!     * Output:
!     *  SMIN    - array containing the minimums across the nodes
!     *  ISTAT   - status of rsum. 0 is OK (MPI_SRC only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

USE MPL, ONLY :                                                   &
    MPL_REAL,                                                     &
    MPL_MIN

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"

INTEGER (KIND=GC_INT_KIND) :: LEN, NPROC, ISTAT
REAL (KIND=GC_REAL_KIND)   :: SMIN(LEN)
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
DO I = 1,LEN
   REDUCE_DATA_WRK(I) = SMIN(I)
ENDDO
CALL MPL_ALLREDUCE(REDUCE_DATA_WRK, SMIN, LEN, MPL_REAL, MPL_MIN, &
     GC__MY_MPI_COMM_WORLD, ISTAT)
#endif

#if defined(SERIAL_SRC)
ISTAT = GC__OK
#endif

#include "gc_end_timer.h"

RETURN
END
