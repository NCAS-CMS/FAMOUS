! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

SUBROUTINE GC__STAMP()

IMPLICIT NONE

#include "gc_kinds.h"

#include "gc_precision.h"

WRITE(6,*)      
WRITE(6,*) '====================================================='
WRITE(6,*) 'GCOM Version ',                                       &
 GC_VERSION
WRITE(6,*)                                                        &
 GC_DESCRIP
WRITE(6,*) 'Using precision : ',                                  &
 GC_INT_TYPE , ' and ', GC_REAL_TYPE
WRITE(6,*) 'Built at ',                                           &
 GC_BUILD_DATE
#if defined(GCOM_TIMER)
WRITE(6,*) 'Timer/Performance Instrumentation is active'
#endif
WRITE(6,*) '====================================================='
WRITE(6,*)

RETURN
END
