! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

!=======================================================================
!  This is an INTERNAL routine to be used within the GCG interface ONLY.
!=======================================================================

#include "gcg_prolog.h"    

#if defined(MPI_SRC)
INTEGER FUNCTION GCG__MPI_RANK(RANK, GID)

!     MPI Global to local rank translation routine. Translates
!     RANK in communicator GC__MY_MPI_COMM_WORLD to that of communicator
!     GID if RANK is a member of GID.

USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD

IMPLICIT NONE

#include "gc_kinds.h"

INTEGER (KIND=GC_INT_KIND) :: RANK, GID
INTEGER (KIND=GC_INT_KIND) :: I, SIZE, IGID, ISTAT
LOGICAL (KIND=GC_LOG_KIND), SAVE :: GINITED = .FALSE.
LOGICAL (KIND=GC_LOG_KIND), SAVE :: INITMPI = .FALSE.
INTEGER (KIND=GC_INT_KIND), SAVE :: GIDW

INTEGER (KIND=GC_INT_KIND), ALLOCATABLE, SAVE :: RANKS(:,:)

IF (.NOT. GINITED) CALL GCG__INIT_INTERNALS(GINITED)

IF (.NOT. INITMPI) THEN
   INITMPI = .TRUE.
   CALL MPL_COMM_GROUP(GC__MY_MPI_COMM_WORLD, GIDW, ISTAT)

   CALL MPL_COMM_SIZE(GC__MY_MPI_COMM_WORLD, SIZE, ISTAT)
   ALLOCATE (RANKS(0:SIZE-1, 2))
   DO I = 0, SIZE-1
      RANKS(I,1) = I
   ENDDO
ENDIF

!---  Translate ranks from GID to GC__MY_MPI_COMM_WORLD GID
CALL MPL_COMM_SIZE(GID, SIZE, ISTAT)
CALL MPL_COMM_GROUP(GID, IGID, ISTAT)
CALL MPL_GROUP_TRANSLATE_RANKS(IGID, SIZE, RANKS(0,1), GIDW,      &
     RANKS(0,2), ISTAT)
CALL MPL_GROUP_FREE(IGID, ISTAT)

!---  Search for specified global rank in GID group
DO I = 0, SIZE-1
   IF (RANKS(I,2)  ==  RANK) THEN
      GCG__MPI_RANK = I
      RETURN
   ENDIF
ENDDO

!---  Global rank RANK not found in GID:
GCG__MPI_RANK = -1

RETURN
END

#else

INTEGER FUNCTION GCG__MPI_RANK(dummy1,dummy2)

IMPLICIT NONE

#include "gc_kinds.h"

INTEGER (KIND=GC_INT_KIND) :: dummy1,dummy2

#include "gc_functions.h"

CALL GC_ABORT(GC_ME(),GC_NPROC(),                                 &
              'GCG__MPI_RANK called for non-MPI')
     
GCG__MPI_RANK = -1

RETURN
END
#endif
