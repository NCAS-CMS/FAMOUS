! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"    

SUBROUTINE GC_GSYNC (NPROC, ISTAT)                  
!     ******************************************************************
!     * Purpose:
!     *  
!     *  Synchronize the processors. Mainly used in front of
!     *  (asynchronous) communication and in connection with timing.
!     *
!     * Input:
!     *  NPROC   - number of nodes 
!     *
!     * Output:
!     *  ISTAT   - status of send 0 is OK (MPI_SRC only),
!     *            refer to the header files for nonzero status codes
!     *
!     * NOTES:       
!     * 
!     *  No node can continue execution before everybody have reached
!     *  this point.
!     ******************************************************************

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"

INTEGER (KIND=GC_INT_KIND) :: NPROC, ISTAT
#if defined(GCOM_TIMER)
#include "gc_timer.h"
#define THIS_TIMER TIMET_BARRIER
#define THIS_LENGTH 0
#endif

#include "gc_start_timer.h"

#if defined(MPI_SRC)
CALL MPL_BARRIER(GC__MY_MPI_COMM_WORLD, ISTAT)
#endif

#if defined(SERIAL_SRC)
ISTAT = GC__OK
#endif

#include "gc_end_timer.h"

RETURN
END
