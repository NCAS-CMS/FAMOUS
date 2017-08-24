! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

!=======================================================================
!  This is an INTERNAL routine to be used within the GC interface ONLY.
!=======================================================================

#include "gcg_prolog.h"    

SUBROUTINE GCG__ERRLIM(IABRT, SUB, LIM, MVAL, AVAL)

!     Support function to exit or abort (if IABRT > 0) if an internal
!     limit is exceeded.

!     Currently, only abort is supported (RS 970408)

IMPLICIT NONE

#include "gc_kinds.h"

INTEGER (KIND=GC_INT_KIND) :: IABRT, MVAL, AVAL
CHARACTER*(*)              :: SUB, LIM
INTEGER (KIND=GC_INT_KIND) :: GCG__NPROC, GCG__ME
#include "gc_functions.h"
EXTERNAL GC_ME, GC_NPROC

GCG__ME = GC_ME()
GCG__NPROC = GC_NPROC()
WRITE(*,*) 'GCG_', SUB, '(): internal limit MAX_', LIM,        &
           ' exceeded on processor ', GCG__ME
WRITE(*,*) 'Maximum value is ', MVAL, '. Actual value is ',    &
           AVAL, '. Exiting.'

CALL GC_ABORT(GCG__ME, GCG__NPROC, '*** DEFINED LIMIT EXCEEDED ***')

RETURN
END
