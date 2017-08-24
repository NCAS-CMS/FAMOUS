! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"    

SUBROUTINE GC_ISUM (LEN, NPROC, ISTAT, ISUM)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Calculate the integer sum across all processors and distribute
!     *  the result to all the processors.
!     *
!     * Input:
!     *  LEN     - number of elements in message
!     *  NPROC   - number of processors
!     *  ISUM    - array with elements to be added up across the nodes
!     *
!     * Output:
!     *  ISUM    - array containing the sums across the nodes
!     *  ISTAT   - status of isum. 0 is OK (MPI_SRC only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

USE MPL, ONLY :                                                   &
    MPL_INTEGER,                                                  &
    MPL_SUM

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"

INTEGER (KIND=GC_INT_KIND) :: LEN, NPROC, ISTAT, ISUM(LEN)
INTEGER (KIND=GC_INT_KIND) :: REDUCE_DATA_IWRK(LEN)

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: I
#endif

#if defined(GCOM_TIMER)
#include "gc_timer.h"
#define THIS_TIMER TIMET_COLL
#define THIS_LENGTH LEN*GC__ISIZE
#endif

#include "gc_start_timer.h"

#if defined(MPI_SRC)
DO I = 1,LEN
   REDUCE_DATA_IWRK(I) = ISUM(I)
ENDDO
CALL MPL_ALLREDUCE(REDUCE_DATA_IWRK, ISUM, LEN, MPL_INTEGER,      &
     MPL_SUM, GC__MY_MPI_COMM_WORLD, ISTAT)

#endif

#if defined(SERIAL_SRC)
ISTAT = GC__OK
#endif

#include "gc_end_timer.h"

RETURN
END
