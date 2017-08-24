! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

!=======================================================================
!  This is an INTERNAL routine to be used within the GC interface ONLY.
!=======================================================================

#include "gc_prolog.h"    

SUBROUTINE GC__GET_MPI_TYPE(LENGTH,BLOCKS,STRIDE,THIS_TYPE)

!     Creates a new MPI derived datatype for the structure described
!     by LENGTH,BLOCKS,STRIDE or returns an existing one

USE MPL, ONLY:                                                    &
    MPL_BYTE,                                                     &
    MPL_ADDRESS_KIND

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY:                                         &
    MPI_TYPES_ARRAY,                                              &
    MAX_MPI_TYPES,                                                &
    N_MPI_TYPES_DEFINED,                                          &
    MPI_TYPE_SEQ
#endif

IMPLICIT NONE

#include "gc_kinds.h"

INTEGER (KIND=GC_INT_KIND) :: LENGTH      ! IN : Number of words in each data block
INTEGER (KIND=GC_INT_KIND) :: BLOCKS      ! IN : Number of blocks, each of length LENGTH
INTEGER (KIND=GC_INT_KIND) :: STRIDE      ! IN : Stride between the start of each block
INTEGER (KIND=GC_INT_KIND) :: THIS_TYPE   ! OUT : MPI reference ID for the new derived type

INTEGER (KIND=GC_INT_KIND) :: ISTAT       ! MPI return code
INTEGER (KIND=GC_INT_KIND) :: TYPE_COUNT  ! Counter through existing types
INTEGER (KIND=GC_INT_KIND) :: TYPE_INDEX  ! Index into MPI_TYPES_ARRAY
INTEGER (KIND=GC_INT_KIND) :: BEST_SEQ_VAL
INTEGER (KIND=GC_INT_KIND) :: INDEX
INTEGER (KIND=GC_INT_KIND) :: I

TYPE_INDEX=-1
#if defined(MPI_SRC)
MPI_TYPE_SEQ=MPI_TYPE_SEQ+1
TYPE_COUNT=0
DO
  TYPE_COUNT=TYPE_COUNT+1
  IF (TYPE_COUNT > N_MPI_TYPES_DEFINED) EXIT
  
  IF ( ( MPI_TYPES_ARRAY(2,TYPE_COUNT) == LENGTH) .AND.           &
       ( MPI_TYPES_ARRAY(3,TYPE_COUNT) == BLOCKS) .AND.           &
       ( MPI_TYPES_ARRAY(4,TYPE_COUNT) == STRIDE) ) THEN
    ! Found an existing type which matches
    TYPE_INDEX=TYPE_COUNT
    EXIT
  ENDIF
ENDDO

IF (TYPE_INDEX  ==  -1) THEN
  ! Need to add a new type
  IF (N_MPI_TYPES_DEFINED >= MAX_MPI_TYPES) THEN
    ! The table is filled up, so we must search 
    ! for least recently used entry
    BEST_SEQ_VAL=MPI_TYPE_SEQ+1
    INDEX=1
    DO I=1,N_MPI_TYPES_DEFINED
      IF (MPI_TYPES_ARRAY(5,I)  <   BEST_SEQ_VAL) THEN
        BEST_SEQ_VAL=MPI_TYPES_ARRAY(5,I)
        INDEX=I
      ENDIF
    ENDDO
  ELSE
    N_MPI_TYPES_DEFINED=N_MPI_TYPES_DEFINED+1
    INDEX=N_MPI_TYPES_DEFINED
  ENDIF
  TYPE_INDEX=INDEX
  
  ! And put the information into the table
  MPI_TYPES_ARRAY(2,TYPE_INDEX)=LENGTH
  MPI_TYPES_ARRAY(3,TYPE_INDEX)=BLOCKS
  MPI_TYPES_ARRAY(4,TYPE_INDEX)=STRIDE
  
  CALL MPL_TYPE_VECTOR(BLOCKS,GC__RSIZE*LENGTH,                   &
                       GC__RSIZE*STRIDE,                          &
                       MPL_BYTE,MPI_TYPES_ARRAY(1,TYPE_INDEX),    &
                       ISTAT)
  CALL MPL_TYPE_COMMIT(MPI_TYPES_ARRAY(1,TYPE_INDEX),ISTAT)
  
ENDIF

! Finally update the sequence number for this entry - this 
! enables us to spot which entries are rarely used and can
! be used for new entries if the table is full

MPI_TYPES_ARRAY(5,TYPE_INDEX)=MPI_TYPE_SEQ
THIS_TYPE=MPI_TYPES_ARRAY(1,TYPE_INDEX)
#endif

RETURN
END
     
SUBROUTINE GC__FREE_MPI_TYPES()

!     Frees MPI internals of the derived data types
#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    N_MPI_TYPES_DEFINED,                                          &
    MPI_TYPES_ARRAY,                                              &
    MPI_TYPE_SEQ
#endif

IMPLICIT NONE

#include "gc_kinds.h"

INTEGER (KIND=GC_INT_KIND) :: I,ISTAT

#if defined(MPI_SRC)      
DO I=1,N_MPI_TYPES_DEFINED
  CALL MPL_TYPE_FREE(MPI_TYPES_ARRAY(1,I),ISTAT)
ENDDO
N_MPI_TYPES_DEFINED=0
MPI_TYPE_SEQ=0
#endif

RETURN
END
