! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gcg_prolog.h"    

SUBROUTINE GCG_RVECSUMF(LVL, LSL, LSO, NV, FIELD, GID,            &
           ISTAT, SUMS)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Calculate fast, but not necessarily reproducible, the real sum 
!     *  of a set of vectors across all members of a group, and 
!     *  distribute the results to all members of the group.
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

USE GCOM_MOD, ONLY :                                              &
    GC_ON,                                                        &
    GC_FORCE_BITREP

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"
#include "gcg_constants.h"

INTEGER (KIND=GC_INT_KIND) :: LVL, LSL, LSO, NV, GID, ISTAT
INTEGER (KIND=GC_INT_KIND) :: OPT
REAL (KIND=GC_REAL_KIND)   :: FIELD(LVL,NV), SUMS(NV)

INTEGER (KIND=GC_INT_KIND) :: I, J

CALL GC_GETOPT(GC_FORCE_BITREP, OPT, ISTAT)
IF (OPT == GC_ON) THEN
  CALL GCG_RVECSUMR(LVL, LSL, LSO, NV, FIELD, GID,                &
                    ISTAT, SUMS)
ELSE
  ISTAT = GC__OK

!---  Set all sums to zero
  DO J = 1, NV
     SUMS(J) = 0.0
  ENDDO

!---  Not necessary to verify GID in this routine, as this is done
!     in GCG_RSUM


#if defined(SERIAL_SRC)
  DO J = 1, NV
     DO I = LSO, LSL+LSO-1
        SUMS(J) = SUMS(J) + FIELD(I,J)
     ENDDO
  ENDDO
#else

!---  Find my sums
  DO J = 1, NV
     DO I = LSO, LSL+LSO-1
        SUMS(J) = SUMS(J) + FIELD(I,J)
     ENDDO
  ENDDO

!---  Call RSUM to do the global sums
  CALL GCG_RSUM(NV, GID, ISTAT, SUMS)

#endif

END IF

RETURN
END
