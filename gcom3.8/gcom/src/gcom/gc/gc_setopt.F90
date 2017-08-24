! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"    

SUBROUTINE GC_SETOPT(VAR, VALUE, ISTAT)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Set a runtime GC configuration option
!     *
!     * Input:
!     *  VAR     - configuration variable to change
!     *  VALUE   - new value of VAR
!     *
!     * Output:
!     *  ISTAT   - status of option setting, 0 is OK.
!     *
!     * NOTES:
!     *  This is a synchronizing call for valid VAR/VALUE pairs, or a
!     *  NO-OP otherwise.
!     *
!     ******************************************************************

USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MAX_OPTS,                                                 &
    GC__OPTIONS

USE GCOM_MOD, ONLY :                                              &
    GC_FORCE_BITREP,                                              &
    GC_ON,                                                        &
    GC_OFF

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"
#include "gc_functions.h"

INTEGER (KIND=GC_INT_KIND) :: VAR, VALUE, ISTAT

! Ensure Option is one we recognise
IF (VAR > GC__MAX_OPTS .OR. VAR < 1) THEN
  CALL GC_ABORT(GC_ME(), GC_NPROC(),                              &
                'Cannot set option - unrecognised')
END IF

! Ensure Option values are recognised
IF (VAR == GC_FORCE_BITREP .AND.                                  &
   (VALUE /= GC_ON .AND. VALUE /= GC_OFF)) THEN
  CALL GC_ABORT(GC_ME(), GC_NPROC(),                              &
                'Cannot set GC_FORCE_BITREP - value unrecognised')
END IF
     

GC__OPTIONS(VAR) = VALUE

ISTAT = GC__OK

RETURN
END
