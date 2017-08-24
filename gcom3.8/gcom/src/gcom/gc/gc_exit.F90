! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"


SUBROUTINE GC_EXIT()
!     ******************************************************************
!     * Purpose:
!     *
!     *  Controlled cleanup of the parallel system.
!     *
!     * Input:
!     *
!     * Output:
!     *
!     * NOTES:
!     *  
!     ******************************************************************
IMPLICIT NONE

#include "gc_kinds.h"

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: INFO
CALL GC__FREE_MPI_TYPES()
CALL MPL_FINALIZE(INFO)
#endif

RETURN
END
