! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

FUNCTION GC_ME ()
!     ******************************************************************
!     * Purpose:
!     *  
!     *  Return my node id, in the range 0...nprocs-1.
!     *
!     * Input:
!     *
!     * Output:
!     *
!     * NOTES:       
!     * 
!     ******************************************************************

USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__ME

IMPLICIT NONE

#include "gc_kinds.h"
INTEGER (KIND=GC_INT_KIND) :: GC_ME

GC_ME = GC__ME

RETURN
END
