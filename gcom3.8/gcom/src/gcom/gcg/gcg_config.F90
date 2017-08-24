! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************
#include "gcg_prolog.h"    

SUBROUTINE GCG_CONFIG (MXGRP, MXROT)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Return information about the GCG configuration.
!     *
!     * Output:
!     *  MXPROC    - maximum number of groups per processor (deprecated)
!     *  MXROT     - maximum number of elements in a rotate (shift)
!     *              operation 
!     *
!     * NOTES:
!     *  Refer to GC_CONFIG() for core GC configuration settings.
!     *    
!     ******************************************************************

USE GC_GLOBALS_MOD, ONLY :                                        &
    MAX_ROTATE

IMPLICIT NONE

#include "gc_kinds.h"

INTEGER (KIND=GC_INT_KIND) :: MXGRP, MXROT


MXGRP = 0                ! Deprecated
MXROT = MAX_ROTATE
END
