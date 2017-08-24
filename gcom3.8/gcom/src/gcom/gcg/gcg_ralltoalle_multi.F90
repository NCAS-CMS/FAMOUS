! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gcg_prolog.h"    

SUBROUTINE GCG_RALLTOALLE_MULTI(                                  &
   SEND_ARRAY, SEND_MAP, N_ITEMS_TO_SEND, SARR_LEN,               &
   RECV_ARRAY, RECV_MAP, N_ITEMS_TO_RECV, RARR_LEN,               &
   GID, FLAG, ISTAT)
!     ******************************************************************
!     * Purpose:
!     *
!     *  An extended all-to-all permutation of real data between the 
!     *  processors in a group. One processor may send several items 
!     *  of data to another processor. Similarly, a processor may 
!     *  receive several items from another processor. This routine 
!     *  may also be used for 1-to-all, all-to-1 and some-to-some
!     *  permutations.
!     *
!     * Input:
!     *  SEND_ARRAY       - array containing all data to be sent
!     *  SEND_MAP         - a map containing the following information
!     *                     for each of the items to be sent:
!     *                        1 - destination processor
!     *                        2 - base address in SEND_ARRAY
!     *                        3 - number of elements in this item 
!     *                        4 - stride between the elements in 
!     *                            SEND_ARRAY
!     *                        5 - element length
!     *                        6 - base address in the receiving 
!     *                            processor's RECV_ARRAY
!     *                        7 - stride between the elements in the
!     *                            receiving processor's RECV_ARRAY
!     *  N_ITEMS_TO_SEND  - total number of items to be sent from this
!     *                     processor
!     *  SARR_LEN         - length of SEND_ARRAY
!     *  RECV_MAP         - a map containing the following information
!     *                     for each of the items to be received:
!     *                        1 - source processor
!     *                        2 - base address in RECV_ARRAY
!     *                        3 - number of elements in this item 
!     *                        4 - stride between the elements in 
!     *                            RECV_ARRAY
!     *                        5 - element length
!     *                        6 - base address in the sending 
!     *                            processor's SEND_ARRAY
!     *                        7 - stride between the elements in the
!     *                            sending processor's SEND_ARRAY
!     *  N_ITEMS_TO_RECV  - total number of items to be received at this
!     *                     processor
!     *  RARR_LEN         - length of RECV_ARRAY
!     *  GID              - processor group ID
!     *  FLAG             - Not currently used. Expected to be used
!     *                     to characterize the permutation
!     *
!     * Output:
!     *  RECV_ARRAY       - array containing the received data, in
!     *                     the structure defined by RECV_MAP. 
!     *  ISTAT            - Status variable. 0 is OK, refer to the
!     *                     header files for nonzero status codes
!     *
!     ******************************************************************

IMPLICIT NONE

#include "gc_kinds.h"

INTEGER (KIND=GC_INT_KIND) :: N_ITEMS_TO_SEND, SARR_LEN,          &
                              SEND_MAP(7,N_ITEMS_TO_SEND)
INTEGER (KIND=GC_INT_KIND) :: N_ITEMS_TO_RECV, RARR_LEN,          &
                              RECV_MAP(7,N_ITEMS_TO_RECV)
INTEGER (KIND=GC_INT_KIND) :: GID, FLAG, ISTAT
REAL (KIND=GC_REAL_KIND)   :: SEND_ARRAY(SARR_LEN),               &
                              RECV_ARRAY(RARR_LEN)

#if defined(SERIAL_SRC)
INTEGER (KIND=GC_INT_KIND) :: I,J,K
#endif
     
#if defined(SERIAL_SRC)
DO K=1,N_ITEMS_TO_SEND
  DO J=1,SEND_MAP(3,K)
    DO I=1,SEND_MAP(5,K)
      RECV_ARRAY(I-1+SEND_MAP(6,K)+                               &
                 (J-1)*SEND_MAP(7,K))=                            &
      SEND_ARRAY(I-1+SEND_MAP(2,K)+                               &
                 (J-1)*SEND_MAP(4,K))
    END DO
  END DO
END DO
#else
CALL GCG__RALLTOALLE_MULTI(                                       &
   SEND_ARRAY, SEND_MAP, N_ITEMS_TO_SEND, SARR_LEN,               &
   RECV_ARRAY, RECV_MAP, N_ITEMS_TO_RECV, RARR_LEN,               &
   GID, FLAG, ISTAT)
#endif

RETURN
END


SUBROUTINE GCG__RALLTOALLE_MULTI(                                 &
   SEND_ARRAY, SEND_MAP, N_ITEMS_TO_SEND, SARR_LEN,               &
   RECV_ARRAY, RECV_MAP, N_ITEMS_TO_RECV, RARR_LEN,               &
   GID, FLAG, ISTAT)

USE MPL, ONLY :                                                   &
    MPL_STATUS_SIZE,                                              &
    MPL_REAL

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"

INTEGER (KIND=GC_INT_KIND) :: N_ITEMS_TO_SEND, SARR_LEN,          &
                              SEND_MAP(7,N_ITEMS_TO_SEND)
INTEGER (KIND=GC_INT_KIND) :: N_ITEMS_TO_RECV, RARR_LEN,          &
                              RECV_MAP(7,N_ITEMS_TO_RECV)
INTEGER (KIND=GC_INT_KIND) :: GID, MAX_BUF_SIZE, FLAG, ISTAT
REAL (KIND=GC_REAL_KIND)   :: SEND_ARRAY(SARR_LEN),               &
                              RECV_ARRAY(RARR_LEN)

INTEGER (KIND=GC_INT_KIND) :: I, J, K, B, L, SBASE,               &
                              SSTRIDE, RBASE, RSTRIDE, TAG

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: ME
INTEGER (KIND=GC_INT_KIND) :: NPROC
INTEGER (KIND=GC_INT_KIND) :: N_HANDLES
INTEGER (KIND=GC_INT_KIND) :: IERROR
INTEGER (KIND=GC_INT_KIND) :: SOURCE
INTEGER (KIND=GC_INT_KIND) :: DEST
INTEGER (KIND=GC_INT_KIND) :: POS

REAL (KIND=GC_REAL_KIND), ALLOCATABLE   :: SEND_DATA(:)
REAL (KIND=GC_REAL_KIND), ALLOCATABLE   :: RECV_DATA(:)

#include "gc_functions.h"

INTEGER (KIND=GC_INT_KIND) :: THIS_TYPE  ! Derived MPI data type for a send/recv
INTEGER (KIND=GC_INT_KIND) :: LENGTH    ! Size of a single block of data
INTEGER (KIND=GC_INT_KIND) :: BLOCKS    ! Number of blocks of data

INTEGER (KIND=GC_INT_KIND), ALLOCATABLE :: RECV_LEN(:)
INTEGER (KIND=GC_INT_KIND), ALLOCATABLE :: SEND_LEN(:)
INTEGER (KIND=GC_INT_KIND), ALLOCATABLE :: HANDLES(:)
INTEGER (KIND=GC_INT_KIND), ALLOCATABLE :: POS1(:)
INTEGER (KIND=GC_INT_KIND), ALLOCATABLE :: STATUS(:,:)
 
LOGICAL (KIND=GC_LOG_KIND) :: SEND_COMPLETED  ! All the outstanding sends have completed
LOGICAL (KIND=GC_LOG_KIND) :: RECV_COMPLETED  ! All the outstanding receives have completed
#if defined(GCOM_TIMER)
INTEGER (KIND=GC_INT_KIND) :: TOTAL_SENDLEN   ! Total size of data sent
INTEGER (KIND=GC_INT_KIND) :: TOTAL_RECVLEN   ! Total size of data received
#endif

#endif
#if defined(GCOM_TIMER)
#include "gc_timer.h"
#endif     

ISTAT = GC__OK

#if defined(MPI_SRC)
ME = GC_ME()
NPROC = GC_NPROC()
ALLOCATE(RECV_LEN(0:NPROC-1))
ALLOCATE(SEND_LEN(0:NPROC-1))
ALLOCATE(HANDLES(2*(NPROC-1)))
ALLOCATE(POS1(0:NPROC-1))

! Loop over all the items I have to receive  
#define THIS_TIMER TIMET_RECV
#define THIS_COUNT 0
#define THIS_LENGTH 0
#if defined(GCOM_TIMER)
TOTAL_RECVLEN=0
#endif
#include "gc_start_timer.h"

RECV_LEN(:) = 0
DO J=1,N_ITEMS_TO_RECV
  SOURCE=RECV_MAP(1,J)
  IF (SOURCE /= ME) THEN
    BLOCKS=RECV_MAP(3,J)
    LENGTH=RECV_MAP(5,J)
    IF (BLOCKS * LENGTH > 0) THEN
      RECV_LEN(SOURCE)=RECV_LEN(SOURCE)+1+BLOCKS*LENGTH
    END IF
  END IF
END DO 
ALLOCATE (RECV_DATA(SUM(RECV_LEN(:))))

POS=1
N_HANDLES=0    
DO I=0, NPROC-1
  IF (RECV_LEN(I) > 0) THEN
    TAG=ME*NPROC+I
    N_HANDLES=N_HANDLES+1
    CALL MPL_IRECV(RECV_DATA(POS),RECV_LEN(I),                    &
      MPL_REAL,I,                                                 &
      TAG,GC__MY_MPI_COMM_WORLD,                                  &
      HANDLES(N_HANDLES),ISTAT)
    IF (ISTAT  /=  0) RETURN
#if defined(GCOM_TIMER)
    TOTAL_RECVLEN=TOTAL_RECVLEN+LENGTH*BLOCKS*GC__RSIZE
#endif
    POS=POS+RECV_LEN(I)
  END IF
END DO
#include "gc_end_timer.h"     
 
! Loop over all the items I have to send
#define THIS_TIMER TIMET_SEND
#define THIS_COUNT 0
#define THIS_LENGTH 0
#if defined(GCOM_TIMER)
TOTAL_SENDLEN=0
#endif
#include "gc_start_timer.h"
SEND_LEN(:) = 0
DO J=1,N_ITEMS_TO_SEND
  DEST=SEND_MAP(1,J)
  IF (DEST /= ME) THEN
    BLOCKS=SEND_MAP(3,J)
    LENGTH=SEND_MAP(5,J)
    IF (BLOCKS * LENGTH > 0) THEN
      SEND_LEN(DEST)=SEND_LEN(DEST)+1+BLOCKS*LENGTH
    END IF
  END IF
END DO 

ALLOCATE (SEND_DATA(SUM(SEND_LEN(:))))

POS=1
DO I=0,NPROC-1
  POS1(I)=POS
  DO J=1,N_ITEMS_TO_SEND
    DEST=SEND_MAP(1,J)
    IF (DEST == I) THEN
      SBASE=SEND_MAP(2,J)
      BLOCKS=SEND_MAP(3,J)
      SSTRIDE=SEND_MAP(4,J)
      LENGTH=SEND_MAP(5,J)
      RBASE=SEND_MAP(6,J)
      RSTRIDE=SEND_MAP(7,J)
      IF (BLOCKS * LENGTH > 0) THEN
        IF (DEST == ME) THEN ! Sending to myself
          DO B=1,BLOCKS
            DO K=1,LENGTH
              RECV_ARRAY(RBASE + (B-1)*RSTRIDE + K-1) =           &
                SEND_ARRAY(SBASE + (B-1)*SSTRIDE + K-1)
            END DO
          END DO
        ELSE
          ! Store RBASE in send data in case maps have different
          ! orders on sender and recipient, this makes it 
          ! unabmiguious
          SEND_DATA(POS)=REAL(RBASE)
          POS=POS+1
          DO B=1,BLOCKS
            DO K=1,LENGTH
              SEND_DATA(POS)=                                     &
                SEND_ARRAY(SBASE + (B-1)*SSTRIDE + K-1)
              POS=POS+1
            END DO
          END DO
        END IF
      END IF
    END IF
  END DO
END DO

DO I=0, NPROC-1
  IF (SEND_LEN(I) > 0) THEN  

    TAG=I*NPROC+ME
    N_HANDLES=N_HANDLES+1
    CALL MPL_ISEND(SEND_DATA(POS1(I)),SEND_LEN(I),                &
      MPL_REAL,I,                                                 &
      TAG,GC__MY_MPI_COMM_WORLD,                                  &
      HANDLES(N_HANDLES),ISTAT)
    IF (ISTAT  /=  0) RETURN
#if defined(GCOM_TIMER)
    TOTAL_SENDLEN=TOTAL_SENDLEN+SEND_LEN(I)*GC__RSIZE
#endif
  END IF
END DO
#include "gc_end_timer.h"      
! Now all the sends and receives have been initiated, we 
! can just sit back and wait for them to complete

ALLOCATE(STATUS(MPL_STATUS_SIZE,N_HANDLES))
CALL MPL_WAITALL(N_HANDLES,HANDLES,STATUS,ierror)
DEALLOCATE(STATUS)

POS=1
DO I=0,NPROC-1
  DO J=1,N_ITEMS_TO_RECV
    SOURCE=RECV_MAP(1,J)
    IF (SOURCE == I .AND. SOURCE /= ME) THEN
      BLOCKS=RECV_MAP(3,J)
      RSTRIDE=RECV_MAP(4,J)
      LENGTH=RECV_MAP(5,J)
      RBASE=INT(RECV_DATA(POS))
      IF (BLOCKS * LENGTH > 0) THEN
        POS=POS+1
        DO B=1,BLOCKS
          DO K=1,LENGTH
            RECV_ARRAY(RBASE+(B-1)*RSTRIDE+K-1) = RECV_DATA(POS)
            POS=POS+1
          END DO
        END DO
      END IF
    END IF
  END DO
END DO

DEALLOCATE(RECV_LEN)
DEALLOCATE(SEND_LEN)
DEALLOCATE(RECV_DATA)
DEALLOCATE(SEND_DATA)
DEALLOCATE(HANDLES)
DEALLOCATE(POS1)

#endif

RETURN
END
