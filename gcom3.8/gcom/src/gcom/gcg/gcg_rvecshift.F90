! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gcg_prolog.h"    

SUBROUTINE GCG_RVECSHIFT (LVL, LSL, LSO, NV, SHFT, WRAP, FIELD,   &
     GID, ISTAT)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Shift (rotate) the elements in a set of vectors distributed
!     *  across all members of a group.
!     *
!     * Input:
!     *  LVL     - Local Vector Length
!     *  LSL     - Local Shift Length (the length of the subsection to 
!     *            be shifted for each vector)
!     *  LSO     - Local Shift Offset (element where the summation 
!     *            start)
!     *  NV      - Number of Vectors
!     *  SHFT    - Number of Shifts to be done
!     *  WRAP    - Logical indicating whether the vectors should be
!     *            wrapped around on shifts
!     *  FIELD   - local array containing the vectors to be shifted
!     *  GID     - processor group ID
!     *
!     * Output:
!     *  FIELD   - Local array containing the shifted data
!     *  ISTAT   - status of rsum. 0 is OK (MPI_SRC only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

USE MPL, ONLY :                                                   &
    MPL_STATUS_SIZE,                                              &
    MPL_BYTE

USE GC_GLOBALS_MOD, ONLY :                                        &
    MAX_ROTATE

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"
#include "gcg_constants.h"
#include "gcg_mtags.h"


INTEGER (KIND=GC_INT_KIND)  :: LVL, LSL, LSO, NV, SHFT, GID, ISTAT
REAL (KIND=GC_REAL_KIND)    :: FIELD(LVL,NV)
LOGICAL (KIND=GC_INT_KIND)  :: WRAP

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND)  :: STATUS(MPL_STATUS_SIZE)
#endif
INTEGER (KIND=GC_INT_KIND)  :: I, J, K, ME, GRANK, GSIZE, GRANK0, &
                               IGID, ILOC, GJ
INTEGER (KIND=GC_INT_KIND)  :: LDLEN
INTEGER (KIND=GC_INT_KIND)  :: LST(2)
INTEGER (KIND=GC_INT_KIND), ALLOCATABLE  :: GLST(:,:)
REAL    (KIND=GC_REAL_KIND) :: GARR(MAX_ROTATE)
REAL    (KIND=GC_REAL_KIND) :: HARR(MAX_ROTATE)

#if defined(GCOM_TIMER)
#include "gc_timer.h"
#define THIS_TIMER TIMET_COLL
#define THIS_LENGTH LVL*NV*GC__RSIZE
#endif

#include "gc_start_timer.h"

ISTAT = GC__OK


!---  Return if no shifts
IF (SHFT  ==  0) RETURN

#if defined(MPI_SRC)
IF (GID  ==  GCG__ALLGROUP) THEN
   IGID = GC__MY_MPI_COMM_WORLD
ELSE
   IGID = GID
ENDIF
CALL MPL_COMM_RANK(IGID, GRANK, ISTAT)
CALL MPL_COMM_SIZE(IGID, GSIZE, ISTAT)
ALLOCATE( GLST(2,0:gsize-1) )
GRANK0 = 0
#endif


!---  Check if one or more processors in the group is asked to shift 
!     more elements than it has
LDLEN = LVL - (LSO + LSL - 1)
CALL GCG_IMIN(1_GC_INT_KIND, GID, ISTAT, LDLEN)
IF (LDLEN < 0 .OR. ISTAT < 0) THEN
   ISTAT = -1
   RETURN
ENDIF

#if defined(SERIAL_SRC)
!---  Loop over all vectors 
DO K = 1,NV

!---  Rotate by copying into a vector, copying this to another vector
!     and copying back (!) Can obviously be simplyfied, but how often
!     is this operation done in serial mode?
   DO J = 1, LSL
      GARR(J) = FIELD(LSO+J-1,K)
   ENDDO
   
   IF (SHFT  >   0) THEN
      DO I = 1, LSL-SHFT
         HARR(I+SHFT) = GARR(I)
      ENDDO
      IF (WRAP) THEN
         DO I = 1, SHFT
            HARR(I) = GARR(LSL-SHFT+I)
         ENDDO
      ELSE
         DO I = 1, SHFT
            HARR(I) = GARR(I)
         ENDDO
      ENDIF
   ELSE
      DO I = 1-SHFT, LSL
         HARR(I+SHFT) = GARR(I)
      ENDDO
      IF (WRAP) THEN
         DO I = 0,-1-SHFT
            HARR(LSL-I) = GARR(-I-SHFT)
         ENDDO
      ELSE
         DO I = 0,-1-SHFT
            HARR(LSL-I) = GARR(LSL-I)
         ENDDO
      ENDIF
   ENDIF

   DO J = 1, LSL
      FIELD(LSO+J-1,K) = HARR(J)
   ENDDO

ENDDO
#else

!---  Send the shift length and offset of all processors to the first
!     processor in the group. 
LST(1) = LSL
LST(2) = LSO
IF (GRANK  /=  0) THEN
   CALL MPL_BSEND(LST, GC__ISIZE*2, MPL_BYTE, GRANK0,             &
        GCGID__ROT0, IGID, ISTAT)
ELSE
   DO I = 1, GSIZE-1
      CALL MPL_RECV(GLST(1,I), GC__ISIZE*2, MPL_BYTE,             &
                    I, GCGID__ROT0,                               &
                    IGID, STATUS, ISTAT)
   ENDDO

!---  Exit if total length of the vectors is greater than MAX_ROTATE

   GLST(1,0) = LSL
   GLST(2,0) = LSO
   GJ = 0
   DO I = 0,GSIZE-1
      GJ = GJ + GLST(1,I)
   ENDDO
   IF (GJ  >   MAX_ROTATE)                                        &
      CALL GCG__ERRLIM(1_GC_INT_KIND, 'RVECSHIFT',                &
                       'ROTATE', MAX_ROTATE, GJ)

ENDIF

!---  Loop over all vectors 
DO K = 1,NV

!---  Send to first processor in the group, which do the rotate
!     and distribute back again
   IF (GRANK  /=  0) THEN
      CALL MPL_BSEND(FIELD(LSO,K), GC__RSIZE*LSL,                 &
                     MPL_BYTE, GRANK0,                            &
                     GCGID__ROT1, IGID, ISTAT)
      CALL MPL_RECV(FIELD(LSO,K), GC__RSIZE*LSL,                  &
                    MPL_BYTE, GRANK0,                             &
                    GCGID__ROT2, IGID, STATUS, ISTAT)

   ELSE
      
      DO J = 1, LSL
         GARR(J) = FIELD(LSO+J-1,K)
      ENDDO
      GJ = LSL
   
      DO  I = 1, GSIZE-1
         CALL MPL_RECV(GARR(GJ+1), GC__RSIZE*GLST(1,I),           &
                       MPL_BYTE, I, GCGID__ROT1,                  &
                       IGID, STATUS, ISTAT)
         GJ = GJ + GLST(1,I)
      ENDDO

      IF (SHFT  >   0) THEN
         IF (WRAP) THEN
            DO I = 1, SHFT
               HARR(I) = GARR(GJ-SHFT+I)
            ENDDO
         ELSE
            DO I = 1, SHFT
               HARR(I) = GARR(I)
            ENDDO
         ENDIF
         DO I = 1, GJ-SHFT
            HARR(I+SHFT) = GARR(I)
         ENDDO
      ELSE
         DO I = 1-SHFT, GJ
            HARR(I+SHFT) = GARR(I)
         ENDDO
         IF (WRAP) THEN
            DO I = 1,-SHFT
               HARR(GJ+SHFT+I) = GARR(I)
            ENDDO
         ELSE
            DO I = 1,-SHFT
               HARR(GJ+SHFT+I) = GARR(GJ+SHFT+I)
            ENDDO
         ENDIF
      ENDIF

      DO J = 1, LSL
         FIELD(LSO+J-1,K) = HARR(J)
      ENDDO
      GJ = LSL

      DO I = 1, GSIZE-1
         CALL MPL_BSEND(HARR(GJ+1), GC__RSIZE*GLST(1,I),          &
                        MPL_BYTE, I, GCGID__ROT2, IGID, ISTAT)
         GJ = GJ + GLST(1,I)
      ENDDO

   ENDIF
ENDDO

DEALLOCATE( GLST )
#endif

#include "gc_end_timer.h"

RETURN
END
