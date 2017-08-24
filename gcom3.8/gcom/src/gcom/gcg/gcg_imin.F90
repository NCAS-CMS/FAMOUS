! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gcg_prolog.h"    

SUBROUTINE GCG_IMIN (LEN, GID, ISTAT, IMIN)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Calculate the integer minimum across all processors of a group
!     *  and distribute the result to all members of the group.
!     *
!     * Input:
!     *  LEN     - number of elements in message
!     *  GID     - processor group ID
!     *  IMIN    - array with elements to be added up across the nodes
!     *
!     * Output:
!     *  IMIN    - array containing the sums across the nodes
!     *  ISTAT   - status of rsum. 0 is OK (MPI_SRC only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

USE MPL, ONLY :                                                   &
    MPL_INTEGER,                                                  &
    MPL_MIN

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"
#include "gcg_constants.h"

INTEGER (KIND=GC_INT_KIND) :: LEN, GID, ISTAT, IMIN(LEN)

INTEGER (KIND=GC_INT_KIND) :: REDUCE_DATA_IWRK(LEN)

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: IGID
#endif

#if defined(GCOM_TIMER)
#include "gc_timer.h"
#define THIS_TIMER TIMET_COLL
#define THIS_LENGTH LEN*GC__ISIZE
#endif
INTEGER (KIND=GC_INT_KIND) :: I

#include "gc_start_timer.h"

ISTAT = GC__OK

#if defined(MPI_SRC)
IF (GID  ==  GCG__ALLGROUP) THEN
   IGID = GC__MY_MPI_COMM_WORLD
ELSE
   IGID = GID
ENDIF
DO I = 1,LEN
   REDUCE_DATA_IWRK(I) = IMIN(I)
ENDDO
CALL MPL_ALLREDUCE(REDUCE_DATA_IWRK, IMIN, LEN, MPL_INTEGER,      &
     MPL_MIN, IGID, ISTAT)
#endif

#include "gc_end_timer.h"

RETURN
END
