! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"    

SUBROUTINE GC_SSYNC (NPROC, ISTAT)                  
!     ******************************************************************
!     * Purpose:
!     *  
!     *  Obsolete - Shared memory synchronization of processors.
!     *
!     * Input:
!     *  NPROC   - number of nodes 
!     *
!     * Output:
!     *  ISTAT   - Return status, GC_OK if success.
!     *
!     * NOTES:       
!     * 
!     ******************************************************************

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"

INTEGER (KIND=GC_INT_KIND) :: NPROC, ISTAT

ISTAT = GC__OK

RETURN
END
