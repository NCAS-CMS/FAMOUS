! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gcg_prolog.h"    

SUBROUTINE GCG_RVECSUMR(LVL, LSL, LSO, NV, FIELD, GID,            &
                        ISTAT, SUMS)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Calculate in a reproducible way the real sum of a set of
!     *  vectors across all members of a group, and distribute the
!     *  results to all members of the group.
!     *
!     * Input:
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
    MPL_BYTE

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"
#include "gcg_constants.h"
#include "gcg_mtags.h"

INTEGER (KIND=GC_INT_KIND) :: LVL, LSL, LSO, NV, GID, ISTAT
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
INTEGER (KIND=GC_INT_KIND) :: G0, GPREV, GNEXT, GLAST

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

CALL MPL_COMM_RANK(IGID, GRANK, ISTAT)
CALL MPL_COMM_SIZE(IGID, GSIZE, ISTAT)
GLAST = GSIZE - 1
GPREV = GRANK - 1
GNEXT = GRANK + 1
#endif

#if defined(SERIAL_SRC)
DO J = 1, NV
  DO I = LSO, LSL+LSO-1
    SUMS(J) = SUMS(J) + FIELD(I,J)
  END DO
END DO
#else

!---  Perform the reproducible global sums for this GID. The first 
!     group member (GRANK = 0) starts. He sums his elements, pass on 
!     the partial results to the next member, and so on until the last 
!     member, which broadcast the sums after adding his contributions.
IF (GRANK  /=  0) THEN
  CALL MPL_RECV(SUMS, GC__RSIZE*NV, MPL_BYTE, GPREV,             &
                GCGID__VEC0, IGID, STATUS, ISTAT)
END IF

DO J = 1, NV
!cdir novector
   DO I = LSO, LSL+LSO-1
      SUMS(J) = SUMS(J) + FIELD(I,J)
   END DO
END DO

IF (GRANK  <   GSIZE-1) THEN
  CALL MPL_SEND(SUMS, GC__RSIZE*NV, MPL_BYTE, GNEXT,            &
                GCGID__VEC0, IGID, ISTAT)
END IF

IF (GSIZE  >   1) THEN
  CALL MPL_BCAST(SUMS, GC__RSIZE*NV, MPL_BYTE, GLAST,            &
                 IGID, ISTAT)
END IF
#endif

#include "gc_end_timer.h"

RETURN
END
