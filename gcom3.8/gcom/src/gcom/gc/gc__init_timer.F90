! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

!=======================================================================
!  This is an INTERNAL routine to be used within the GC interface ONLY.
!=======================================================================

#include "gc_prolog.h"    

SUBROUTINE GC__INIT_TIMER()

! Support subroutine to initialise the GCOM timer internals

IMPLICIT NONE

#include "gc_kinds.h"

#if defined(GCOM_TIMER)
INTEGER (KIND=GC_INT_KIND) :: I
#include "gc_timer.h"

DO I=1,TIMER_NTYPES
  GCOM_TIMER_TIME(I)=0.0
  GCOM_TIMER_COUNT(I)=0
  GCOM_TIMER_DATA(I)=0
ENDDO
#endif

RETURN
END
