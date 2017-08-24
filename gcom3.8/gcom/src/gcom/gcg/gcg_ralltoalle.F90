! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gcg_prolog.h"    

SUBROUTINE GCG_RALLTOALLE(                                        &
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

#if defined(MPI_SRC)
USE GCOM_MOD, ONLY:                                               &
    GC_ALLTOALL_VERSION,                                          &
    GC_ALLTOALL_ORIG
#endif 

IMPLICIT NONE

#include "gc_kinds.h"

INTEGER (KIND=GC_INT_KIND) :: N_ITEMS_TO_SEND, SARR_LEN,          &
                              SEND_MAP(7,N_ITEMS_TO_SEND)
INTEGER (KIND=GC_INT_KIND) :: N_ITEMS_TO_RECV, RARR_LEN,          &
                              RECV_MAP(7,N_ITEMS_TO_RECV)
INTEGER (KIND=GC_INT_KIND) :: GID, FLAG, ISTAT
INTEGER (KIND=GC_INT_KIND) :: OPT
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
    ENDDO
  ENDDO
ENDDO
#else
CALL GC_GETOPT(GC_ALLTOALL_VERSION, OPT, ISTAT)
IF (OPT == GC_ALLTOALL_ORIG) THEN
  CALL GCG__RALLTOALLE(                                           &
     SEND_ARRAY, SEND_MAP, N_ITEMS_TO_SEND, SARR_LEN,             &
     RECV_ARRAY, RECV_MAP, N_ITEMS_TO_RECV, RARR_LEN,             &
     GID, FLAG, ISTAT)
ELSE
  CALL GCG__RALLTOALLE_MULTI(                                     &
     SEND_ARRAY, SEND_MAP, N_ITEMS_TO_SEND, SARR_LEN,             &
     RECV_ARRAY, RECV_MAP, N_ITEMS_TO_RECV, RARR_LEN,             &
     GID, FLAG, ISTAT)
END IF
#endif

RETURN
END


SUBROUTINE GCG__RALLTOALLE(                                       &
   SEND_ARRAY, SEND_MAP, N_ITEMS_TO_SEND, SARR_LEN,               &
   RECV_ARRAY, RECV_MAP, N_ITEMS_TO_RECV, RARR_LEN,               &
   GID, FLAG, ISTAT)

USE MPL, ONLY :                                                   &
    MPL_STATUS_SIZE

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD,                                        &
    GC__MPI_MAXTAG
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

INTEGER (KIND=GC_INT_KIND) :: I, J, K, L, LENGTH, LBASE, LSTRIDE, &
                              RBASE, RSTRIDE, TAG

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: ME

#include "gc_functions.h"

INTEGER (KIND=GC_INT_KIND) :: THIS_TYPE  ! Derived MPI data type for a send/recv
INTEGER (KIND=GC_INT_KIND) :: LLENGTH    ! Size of a single block of data
INTEGER (KIND=GC_INT_KIND) :: LBLOCKS    ! Number of blocks of data
INTEGER (KIND=GC_INT_KIND) :: SEND_HANDLE(N_ITEMS_TO_SEND)  ! Handles for sent data
INTEGER (KIND=GC_INT_KIND) :: RECV_HANDLE(N_ITEMS_TO_RECV)  ! Handles for received data
INTEGER (KIND=GC_INT_KIND) :: SEND_STATUS(MPL_STATUS_SIZE,        &
                                          N_ITEMS_TO_SEND)  ! Status (sends)
INTEGER (KIND=GC_INT_KIND) :: RECV_STATUS(MPL_STATUS_SIZE,        &
                                          N_ITEMS_TO_RECV)  ! Status (recvs)
INTEGER (KIND=GC_INT_KIND) :: N_SEND_HANDLES  ! Number of messages sent
INTEGER (KIND=GC_INT_KIND) :: N_RECV_HANDLES  ! Number of messages received
 
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

! Loop over all the items I have to receive  
#define THIS_TIMER TIMET_RECV
#define THIS_COUNT 0
#define THIS_LENGTH 0
#if defined(GCOM_TIMER)
TOTAL_RECVLEN=0
#endif
#include "gc_start_timer.h"
N_RECV_HANDLES=0    
DO J=1,N_ITEMS_TO_RECV
  
  IF (RECV_MAP(1,J)  /=  ME) THEN  ! Only receive if the message
                                   ! is not from myself
    LBASE=RECV_MAP(2,J)
    LBLOCKS=RECV_MAP(3,J)
    LSTRIDE=RECV_MAP(4,J)
    LLENGTH=RECV_MAP(5,J)
    
    ! Get an MPI derived type which describes this data structure
    CALL GC__GET_MPI_TYPE(LLENGTH,LBLOCKS,LSTRIDE,THIS_TYPE)
    
    TAG=MOD(LBASE,INT(GC__MPI_MAXTAG,GC_INT_KIND))
    N_RECV_HANDLES=N_RECV_HANDLES+1
    CALL MPL_IRECV(RECV_ARRAY(LBASE),1_GC_INT_KIND,               &
                   THIS_TYPE,RECV_MAP(1,J),                       &
                   TAG,GC__MY_MPI_COMM_WORLD,                     &
                   RECV_HANDLE(N_RECV_HANDLES),ISTAT)
    IF (ISTAT  /=  0) RETURN
#if defined(GCOM_TIMER)
    TOTAL_RECVLEN=TOTAL_RECVLEN+LLENGTH*LBLOCKS*GC__RSIZE
#endif
    
  ENDIF
ENDDO ! DO J=1,N_ITEMS_TO_RECV
#include "gc_end_timer.h"     
 
! Loop over all the items I have to send
#define THIS_TIMER TIMET_SEND
#define THIS_COUNT 0
#define THIS_LENGTH 0
#if defined(GCOM_TIMER)
TOTAL_SENDLEN=0
#endif
#include "gc_start_timer.h"
N_SEND_HANDLES=0      
DO J=1,N_ITEMS_TO_SEND

  LBASE=SEND_MAP(2,J)
  LBLOCKS=SEND_MAP(3,J)
  LSTRIDE=SEND_MAP(4,J)
  LLENGTH=SEND_MAP(5,J)
  
  RBASE=SEND_MAP(6,J)
  RSTRIDE=SEND_MAP(7,J)
  
  IF (SEND_MAP(1,J)  ==  ME) THEN ! Sending to myself
    DO I=1,LBLOCKS
      DO K=1,LLENGTH
        RECV_ARRAY(RBASE + (I-1)*RSTRIDE + K-1) =                 &
        SEND_ARRAY(LBASE + (I-1)*LSTRIDE + K-1)
      ENDDO
    ENDDO
  ELSE ! Sending to another processor
    
    CALL GC__GET_MPI_TYPE(LLENGTH,LBLOCKS,LSTRIDE,THIS_TYPE)
     
    TAG=MOD(RBASE,INT(GC__MPI_MAXTAG,GC_INT_KIND))
    N_SEND_HANDLES=N_SEND_HANDLES+1
    CALL MPL_ISEND(SEND_ARRAY(LBASE),1_GC_INT_KIND,               &
                  THIS_TYPE,SEND_MAP(1,J),                        &
                   TAG,GC__MY_MPI_COMM_WORLD,                     &
                   SEND_HANDLE(N_SEND_HANDLES),ISTAT)
    IF (ISTAT  /=  0) RETURN
#if defined(GCOM_TIMER)
    TOTAL_SENDLEN=TOTAL_SENDLEN+LLENGTH*LBLOCKS*GC__RSIZE
#endif          
  ENDIF
ENDDO
#include "gc_end_timer.h"      
! Now all the sends and receives have been initiated, we 
! can just sit back and wait for them to complete

#include "gc_start_timer.h"      
SEND_COMPLETED=.FALSE.
RECV_COMPLETED=.FALSE.
DO  ! Spin loop waiting for completion
  IF ((SEND_COMPLETED) .AND. (RECV_COMPLETED)) EXIT
  
  IF (.NOT. SEND_COMPLETED) THEN
    CALL MPL_TESTALL(N_SEND_HANDLES,SEND_HANDLE,                  &
                     SEND_COMPLETED,SEND_STATUS,ISTAT)
#if defined(GCOM_TIMER)
    IF (SEND_COMPLETED) THEN
#define THIS_TIMER TIMET_SEND
#define THIS_LENGTH TOTAL_SENDLEN
#define THIS_COUNT N_SEND_HANDLES
#include "gc_end_timer.h"
    ENDIF
#endif
  ENDIF
  
  IF (.NOT. RECV_COMPLETED) THEN
    CALL MPL_TESTALL(N_RECV_HANDLES,RECV_HANDLE,                  &
                     RECV_COMPLETED,RECV_STATUS,ISTAT)
#if defined(GCOM_TIMER)
    IF (RECV_COMPLETED) THEN
#define THIS_TIMER TIMET_RECV
#define THIS_LENGTH TOTAL_RECVLEN
#define THIS_COUNT N_RECV_HANDLES
#include "gc_end_timer.h"
    ENDIF
#endif
  ENDIF
ENDDO    
  
#endif

RETURN
END
