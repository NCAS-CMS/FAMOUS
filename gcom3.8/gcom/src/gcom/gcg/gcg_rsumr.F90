! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gcg_prolog.h"   

SUBROUTINE GCG_RSUMR (LEN, GID, ISTAT, RSUM)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Calculate in a reproducible way the real sum across all 
!     *  processors of a group and distribute the result to all members 
!     *  of the group.
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
    MPL_STATUS_SIZE,                                              &
    MPL_BYTE

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"
#include "gcg_constants.h"
#include "gcg_mtags.h"

INTEGER (KIND=GC_INT_KIND) :: LEN, GID, ISTAT
REAL (KIND=GC_REAL_KIND)   :: RSUM(LEN)


#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: L, NPROC, IGID, GRANK, GSIZE
INTEGER (KIND=GC_INT_KIND) :: STATUS(MPL_STATUS_SIZE)

REAL (KIND=GC_REAL_KIND), ALLOCATABLE   :: REDUCE_DATA_WRK(:,:)
INTEGER (KIND=GC_INT_KIND), ALLOCATABLE :: REQUESTS(:)
#include "gc_functions.h"
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
IF (GID  ==  GCG__ALLGROUP) THEN
  NPROC = GC_NPROC()
  CALL GC_RSUMR(LEN, NPROC, ISTAT, RSUM)
  RETURN
ELSE
  IGID = GID
  CALL MPL_COMM_RANK(IGID, GRANK, ISTAT)
  CALL MPL_COMM_SIZE(IGID, GSIZE, ISTAT)

  ALLOCATE( REQUESTS(GSIZE) )
  ALLOCATE( REDUCE_DATA_WRK(LEN, GSIZE) )

  IF (GRANK  /=  0) THEN
    CALL MPL_SEND(RSUM, GC__RSIZE*LEN, MPL_BYTE,               &
                  0_GC_INT_KIND, GCGID__VEC0, IGID, ISTAT)
  ELSE
  ! Issue all receives - non-blocking
    DO I = 1, GSIZE-1
      CALL MPL_IRECV(REDUCE_DATA_WRK(:,I), GC__RSIZE*LEN,      &
                    MPL_BYTE, I, GCGID__VEC0, IGID,            &
                    REQUESTS(I), ISTAT)
    END DO

  ! Wait for receives to complete in order and then add in contribution
    DO I = 1, GSIZE-1
      CALL MPL_WAIT(REQUESTS(I), STATUS, ISTAT)
      
      DO L = 1,LEN
        RSUM(L) = RSUM(L) + REDUCE_DATA_WRK(L,I)
      END DO
    END DO
  END IF
  CALL MPL_BCAST(RSUM, GC__RSIZE*LEN, MPL_BYTE,                  &
                 0_GC_INT_KIND, IGID, ISTAT)

  DEALLOCATE( REQUESTS )
  DEALLOCATE( REDUCE_DATA_WRK )
END IF
#endif

#include "gc_end_timer.h"
RETURN
END
