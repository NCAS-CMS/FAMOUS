! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

FUNCTION GC_NPROC ()
!     ******************************************************************
!     * Purpose:
!     *  
!     *  Return the total number of processors.
!     *
!     * Input:
!     *
!     * Output:
!     *
!     * NOTES:       
!     * 
!     ******************************************************************

USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__NPROC

IMPLICIT NONE

#include "gc_kinds.h"
INTEGER (KIND=GC_INT_KIND) :: GC_NPROC

GC_NPROC = GC__NPROC

RETURN
END
