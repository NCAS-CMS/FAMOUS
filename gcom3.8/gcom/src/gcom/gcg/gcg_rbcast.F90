! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gcg_prolog.h"    

SUBROUTINE GCG_RBCAST (MSG, LEN, SEND, GID, ISTAT, SARR)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Broadcast a real array to every processor of a group.
!     *
!     * Input:
!     *  MSG     - message tag
!     *  LEN     - number of elements in message
!     *  SEND    - sender of the message (global rank)
!     *  GID     - processor group ID
!     *  SARR    - array to be sent
!     *
!     * Output:
!     *  SARR    - array to be received (on nodes != SEND)
!     *  ISTAT   - status of bcast. 0 is OK (MPI_SRC only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
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
#include "gcg_constants.h"

INTEGER (KIND=GC_INT_KIND) :: MSG, LEN, SEND, GID, ISTAT
REAL (KIND=GC_REAL_KIND)   :: SARR(LEN)

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: IGID, RANK, GCG__MPI_RANK
#endif

#if defined(GCOM_TIMER)
#include "gc_timer.h"
#define THIS_LENGTH LEN*GC__RSIZE
#include "gc_functions.h"
#endif

#include "gc_start_timer.h"

ISTAT = GC__OK

#if defined(MPI_SRC)
IF (GID  ==  GCG__ALLGROUP) THEN
   IGID = GC__MY_MPI_COMM_WORLD
ELSE
   IGID = GID
ENDIF
ISTAT = GCG__MPI_RANK(SEND, IGID)
IF (ISTAT  ==  -1) RETURN
RANK = ISTAT
CALL MPL_BCAST(SARR, GC__RSIZE*LEN, MPL_BYTE, RANK,               &
               IGID, ISTAT)
#endif

#if defined(GCOM_TIMER)
IF (GC_ME()  ==  SEND) THEN
#define THIS_TIMER TIMET_BCAST_SEND
#include "gc_end_timer.h"
ELSE ! At the receiving end
#define THIS_TIMER TIMET_BCAST_RECV
#include "gc_end_timer.h"
ENDIF
#endif

RETURN
END
