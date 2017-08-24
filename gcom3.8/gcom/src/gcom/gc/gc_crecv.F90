! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"    

SUBROUTINE GC_CRECV (MSG, LEN, SEND, ISTAT, RARR, SARR)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Receive a character array from processor SEND.
!     *
!     * Input:
!     *  MSG     - message tag
!     *  LEN     - number of BYTES in message
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
!     *
!     ******************************************************************

USE MPL, ONLY :                                                   &
    MPL_STATUS_SIZE,                                              &
    MPL_CHARACTER,                                                &
    MPL_ANY_SOURCE

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"

INTEGER (KIND=GC_INT_KIND) :: MSG, LEN, SEND, ISTAT
CHARACTER*(*)              :: RARR, SARR

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: STATUS(MPL_STATUS_SIZE)
INTEGER (KIND=GC_INT_KIND) :: ME,I
#include "gc_functions.h"
#endif

#if defined(GCOM_TIMER)
#include "gc_timer.h"
#define THIS_TIMER TIMET_RECV
#define THIS_LENGTH LEN
#endif

#include "gc_start_timer.h"

#if defined(MPI_SRC)
ME=GC_ME()
IF (SEND  ==  ME) THEN  ! Receiving from myself
  DO I=1,LEN
    RARR(I:I)=SARR(I:I)
  ENDDO
ELSE ! Receiving from another processor
  IF (SEND  ==  GC__ANY) THEN
     CALL MPL_RECV(RARR, LEN, MPL_CHARACTER, MPL_ANY_SOURCE, MSG, &
          GC__MY_MPI_COMM_WORLD, STATUS, ISTAT)
  ELSE
     CALL MPL_RECV(RARR, LEN, MPL_CHARACTER, SEND, MSG,           &
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
