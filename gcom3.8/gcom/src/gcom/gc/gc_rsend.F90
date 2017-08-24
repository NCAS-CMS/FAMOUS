! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"    

SUBROUTINE GC_RSEND (MSG, LEN, RECI, ISTAT, RARR, SARR)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Send a real array from this processor to processor RECI.
!     *
!     * Input:
!     *  MSG     - message tag
!     *  LEN     - number of elements in message
!     *  RECI    - receiver of the message
!     *  RARR    - name of the array on recieving processor
!     *            (Obsolete)
!     *  SARR    - array to be sent
!     *
!     * Output:
!     *  ISTAT   - status of send 0 is OK (MPI_SRC only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:  
!     *  The use of ISTAT as an input argument is obsoleted. Use
!     *  GC_SETOPT().     
!     *
!     ******************************************************************

USE MPL, ONLY :                                                   &
    MPL_BYTE

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"

INTEGER (KIND=GC_INT_KIND) :: MSG, LEN, RECI, ISTAT
REAL (KIND=GC_REAL_KIND)   :: RARR(LEN), SARR(LEN)

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: ME
#include "gc_functions.h"
#endif

#if defined(SERIAL_SRC)
INTEGER (KIND=GC_INT_KIND) :: I
#endif
#if defined(GCOM_TIMER)
#include "gc_timer.h"
#define THIS_TIMER TIMET_SEND
#define THIS_LENGTH LEN*GC__RSIZE
#endif

#include "gc_start_timer.h"

#if defined(MPI_SRC)
ME=GC_ME()
IF (RECI  /=  ME) THEN    ! If I'm not sending to myself
   CALL MPL_BSEND(SARR, GC__RSIZE*LEN, MPL_BYTE, RECI, MSG,       &
        GC__MY_MPI_COMM_WORLD, ISTAT)
ENDIF
#endif

#if defined(SERIAL_SRC)
DO I=1,LEN
  RARR(I)=SARR(I)
ENDDO

ISTAT = GC__OK
#endif

#include "gc_end_timer.h"

RETURN
END
