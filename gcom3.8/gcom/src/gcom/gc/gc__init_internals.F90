! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

!=======================================================================
!  This is an INTERNAL routine to be used within the GC interface ONLY.
!=======================================================================

#include "gc_prolog.h"    

SUBROUTINE GC__INIT_INTERNALS(INIT)

!     Support function to initialize GC internals. Block data can
!     not be used with CRI MPP f90.

USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__NPROC,                                                    &
    GC__ME,                                                       &
    GC__INITED

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    N_MPI_TYPES_DEFINED,                                          &
    MPI_TYPE_SEQ
#endif

IMPLICIT NONE

#include "gc_kinds.h"

LOGICAL (KIND=GC_LOG_KIND) :: INIT

#if defined(MPI_SRC)
N_MPI_TYPES_DEFINED=0
MPI_TYPE_SEQ=0
#endif
GC__NPROC = -1
GC__ME = -1
GC__INITED = .FALSE.
INIT = .TRUE.

RETURN
END
