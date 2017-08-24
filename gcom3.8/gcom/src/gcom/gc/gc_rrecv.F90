! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"    

SUBROUTINE GC_RRECV (MSG, LEN, SEND, ISTAT, RARR, SARR)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Receive a real array from processor SEND.
!     *
!     * Input:
!     *  MSG     - message tag
!     *  LEN     - number of elements in message
!     *  SEND    - sender of the message (SEND = GC_ANY means any
!     *            processor)
!     *  SARR    - name of the array on the sending processor
!     *            (Obsolete)
!     *
!     * Output:
!     *  RARR    - array to be received
!     *  ISTAT   - status of send 0 is OK (MPI_SRC only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *  The use of ISTAT as an input argument is obsoleted. Use
!     *  GC_SETOPT().     
!     *
!     ******************************************************************

USE MPL, ONLY :                                                   &
    MPL_STATUS_SIZE,                                              &
    MPL_ANY_SOURCE,                                               &
    MPL_BYTE

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"

INTEGER (KIND=GC_INT_KIND) :: MSG, LEN, SEND, ISTAT
REAL (KIND=GC_REAL_KIND)   :: RARR(LEN), SARR(LEN)

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: STATUS(MPL_STATUS_SIZE)
INTEGER (KIND=GC_INT_KIND) :: ME,I
#include "gc_functions.h"
#endif
#if defined(GCOM_TIMER)
#include "gc_timer.h"
#define THIS_TIMER TIMET_RECV
#define THIS_LENGTH LEN*GC__RSIZE
#endif

#include "gc_start_timer.h"

#if defined(MPI_SRC)
ME=GC_ME()
IF (SEND  ==  ME) THEN  ! Receiving from myself
  DO I=1,LEN
    RARR(I)=SARR(I)
  ENDDO
ELSE ! Receiving from another processor
  IF (SEND  ==  GC__ANY) THEN
     CALL MPL_RECV(RARR, GC__RSIZE*LEN, MPL_BYTE,                 &
          MPL_ANY_SOURCE,                                         &
          MSG, GC__MY_MPI_COMM_WORLD, STATUS, ISTAT)
  ELSE
     CALL MPL_RECV(RARR, GC__RSIZE*LEN, MPL_BYTE,                 &
          SEND, MSG,                                              &
          GC__MY_MPI_COMM_WORLD, STATUS, ISTAT)
  ENDIF
ENDIF ! IF (SEND  ==  ME)
#endif

#if defined(SERIAL_SRC)
ISTAT = GC__OK
#endif

#include "gc_end_timer.h"

RETURN
END
