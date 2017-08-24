! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"    

SUBROUTINE GC_GETOPT(VAR, VALUE, ISTAT)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Get a runtime GC configuration option
!     *
!     * Input:
!     *  VAR     - configuration variable to read
!     *
!     * Output:
!     *  VALUE   - present value of VAR
!     *  ISTAT   - status of option read, 0 is OK.
!     *
!     * NOTES:       
!     *
!     ******************************************************************

USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MAX_OPTS,                                                 &
    GC__OPTIONS

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"
#include "gc_functions.h"

INTEGER (KIND=GC_INT_KIND) :: VAR, VALUE, ISTAT

IF (VAR > GC__MAX_OPTS .OR. VAR < 1) THEN
  ! VAR is out of range
  CALL GC_ABORT(GC_ME(), GC_NPROC(),                              &
                'Cannot get option - out of range')
END IF

VALUE = GC__OPTIONS(VAR)
ISTAT = GC__OK

RETURN
END
