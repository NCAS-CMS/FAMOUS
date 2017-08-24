! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

!=======================================================================
!  This is an INTERNAL routine to be used within the GC interface ONLY.
!=======================================================================

#include "gc_prolog.h"   

SUBROUTINE GC__ERRLIM(IABRT, SUB, LIM, MVAL, AVAL)

!     Support function to exit or abort (if IABRT > 0) if an internal
!     limit is exceeded.

!     Currently, only abort is supported (RS 970408)
!
IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_functions.h"

INTEGER (KIND=GC_INT_KIND) :: IABRT, MVAL, AVAL
CHARACTER*(*) SUB, LIM

WRITE(*,*) 'GC_', SUB, '(): internal limit MAX_', LIM,         &
           ' exceeded on processor ', GC_ME()
WRITE(*,*) 'Maximum value is ', MVAL, '. Actual value is ',    &
           AVAL, '. Exiting.'

CALL GC_ABORT(GC_ME(), GC_NPROC(), '*** STATIC LIMIT EXCEEDED ***')

RETURN
END
