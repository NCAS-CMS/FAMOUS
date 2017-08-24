! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"    

SUBROUTINE GC_CSEND (MSG, LEN, RECI, ISTAT, RARR, SARR)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Send a character array from this processor to processor RECI.
!     *
!     * Input:
!     *  MSG     - message tag
!     *  LEN     - number of characters in message
!     *  RECI    - receiver of the message
!     *  RARR    - name of the array on recieving processor
!     *            (Obsolete)
!     *  SARR    - array to be sent
!     *
!     * Output:
!     *  ISTAT   - status of send. 0 is OK (MPI_SRC only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     *
!     ******************************************************************

USE MPL, ONLY :                                                   &
    MPL_CHARACTER

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"

INTEGER (KIND=GC_INT_KIND) :: MSG, LEN, RECI, ISTAT
CHARACTER*(*)              :: RARR, SARR

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: ME
#include "gc_functions.h"
#endif

#if defined(GCOM_TIMER)
#include "gc_timer.h"
#define THIS_TIMER TIMET_SEND
#define THIS_LENGTH LEN
#endif

#include "gc_start_timer.h"


#if defined(MPI_SRC)
ME=GC_ME()
IF (RECI  /=  ME) THEN  ! If I'm not sending to myself
   CALL MPL_BSEND(SARR, LEN, MPL_CHARACTER, RECI, MSG,            &
        GC__MY_MPI_COMM_WORLD, ISTAT)
ENDIF
#endif

#if defined(SERIAL_SRC)
RARR=SARR
ISTAT = GC__OK
#endif

#include "gc_end_timer.h"

RETURN
END
