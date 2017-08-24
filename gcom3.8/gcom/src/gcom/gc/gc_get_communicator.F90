! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"    

SUBROUTINE GC_GET_COMMUNICATOR(VALUE, ISTAT)
!     ******************************************************************
!     * Purpose:
!     *
!     *  For MPI_SRC, this routine returns the 
!     *  the active communicator used by GCOM. For other 
!     *  versions, returns GC_OK.
!     *
!     *
!     * Input:
!     *  none.
!     *
!     * Output:
!     *  VALUE   - present value of communicator (or GC_OK)
!     *  ISTAT   - status of option read, 0 is OK.
!     *
!     * NOTES:       
!     *
!     ******************************************************************

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"

INTEGER (KIND=GC_INT_KIND), INTENT(OUT) :: VALUE
INTEGER (KIND=GC_INT_KIND), INTENT(OUT) :: ISTAT

#if defined(MPI_SRC)
VALUE = GC__MY_MPI_COMM_WORLD
#else
VALUE=GC__OK
#endif

ISTAT=GC__OK

RETURN
END


