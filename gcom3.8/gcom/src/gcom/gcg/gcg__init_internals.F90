! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

!=======================================================================
!  This is an INTERNAL routine to be used within the GCG interface ONLY.
!=======================================================================

#include "gcg_prolog.h"  
  
SUBROUTINE GCG__INIT_INTERNALS(INIT)

IMPLICIT NONE

#include "gc_kinds.h"

LOGICAL (KIND=GC_LOG_KIND) :: INIT
LOGICAL (KIND=GC_LOG_KIND), SAVE :: GINITED = .FALSE.

!     Set GINITED to make sure the subroutine body is processed once only.
INIT = .TRUE.
IF (GINITED) RETURN

GINITED = .TRUE.

RETURN
END
