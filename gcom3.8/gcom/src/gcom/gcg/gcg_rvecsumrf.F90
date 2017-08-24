! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gcg_prolog.h"    

SUBROUTINE GCG_RVECSUMRF(GVL, LVL, LSL, LSO, NV, FIELD, GID,      &
                        ISTAT, SUMS)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Calculate in a reproducible way the real sum of a set of
!     *  vectors across all members of a group, and distribute the
!     *  results to all members of the group.
!     *
!     * Input:
!     *  GVL       Global vector length
!     *  LVL     - Local Vector Length
!     *  LSL     - Local Sum Length (the length of the subsection to be
!     *            summed for each vector)
!     *  LSO     - Local Sum Offset (element where the summation start)
!     *  NV      - Number of Vectors
!     *  FIELD   - local array containing the vectors to be summed
!     *  GID     - processor group ID
!     *
!     * Output:
!     *  SUMS    - array containing the sums across the nodes
!     *  ISTAT   - status of rsum. 0 is OK (MPI_SRC only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

USE MPL, ONLY :                                                   &
    MPL_STATUS_SIZE,                                              &
    MPL_BYTE,                                                     &
    MPL_REAL, MPL_MAX, MPL_SUM, MPL_INTEGER

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"
#include "gcg_constants.h"

INTEGER (KIND=GC_INT_KIND) :: GVL, LVL, LSL, LSO, NV, GID, ISTAT
REAL (KIND=GC_REAL_KIND)   :: FIELD(LVL,NV), SUMS(NV)

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: STATUS(MPL_STATUS_SIZE)
#endif

#if defined(GCOM_TIMER)
#include "gc_timer.h"
#define THIS_TIMER TIMET_COLL
#define THIS_LENGTH LSL*NV*GC__RSIZE
#endif

INTEGER (KIND=GC_INT_KIND) :: I, J, ME, GRANK, GSIZE, IGID, ILOC
INTEGER (KIND=GC_INT_KIND) :: LSLT
INTEGER (KIND=GC_INT_KIND) :: ZERO_COUNT
INTEGER (KIND=GC_INT_KIND) :: ISUM(NV)
INTEGER (KIND=GC_INT_KIND) :: JSUM(NV)
REAL (KIND=GC_REAL_KIND) :: TMAX(NV)
REAL (KIND=GC_REAL_KIND) :: TTMAX(NV)
REAL (KIND=GC_REAL_KIND) :: FAC(NV)

#include "gc_start_timer.h"

ISTAT = GC__OK

!---  Set all sums to zero
DO J = 1, NV
   SUMS(J) = 0.0
END DO

#if defined(MPI_SRC)
IF (GID  ==  GCG__ALLGROUP) THEN
   IGID = GC__MY_MPI_COMM_WORLD
ELSE
   IGID = GID
END IF
#endif

#if defined(SERIAL_SRC)
DO J = 1, NV
   DO I = LSO, LSL+LSO-1
      SUMS(J) = SUMS(J) + FIELD(I,J)
   END DO
END DO
#else
! Algorithm obtains a reproducible sum by doing it in fixed
! precision (integer) arithmetic. This loses accuracy so needs
! to be used with care.

! Fist find the maximum data size over all processors
TMAX(:) = 0
DO J=1,NV
  DO I=LSO,LSL+LSO-1
    TMAX(J)=MAX(TMAX(J),ABS(FIELD(I,J)))
  END DO
END DO

LSLT = GVL
CALL MPL_ALLREDUCE(TMAX, TTMAX, NV, MPL_REAL, MPL_MAX,            &
                   IGID, ISTAT)

! And find scaling factors
DO J=1,NV
  FAC(J)=(TTMAX(J)*LSLT)/2**62
END DO

! Sum and convert to integer in 1 dimension in one go
ZERO_COUNT=0
ISUM(:)=0
DO J=1,NV
  IF(FAC(J) /= 0.0) THEN
    DO I=LSO,LSL+LSO-1
      ISUM(J)=ISUM(J)+INT(FIELD(I,J)*(1.0/FAC(J)))
    END DO
  ELSE
    ZERO_COUNT=ZERO_COUNT+1
  END IF
END DO

! If field isn't entirely zero, do the summation across processors
! and convert back to real.
IF (ZERO_COUNT /= NV) THEN
  CALL MPL_ALLREDUCE(ISUM, JSUM, NV, MPL_INTEGER, MPL_SUM,        &
                     IGID, ISTAT)
  DO J=1,NV
    SUMS(J)=JSUM(J)*FAC(J)
  END DO

ELSE

  DO J=1,NV
    SUMS(J)=0.
  END DO

END IF
#endif

#include "gc_end_timer.h"

RETURN
END
