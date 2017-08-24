! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"    

SUBROUTINE GC_SET_COMMUNICATOR(VALUE, ME, NPROC, ISTAT)
!     ******************************************************************
!     * Purpose:
!     *
!     *  For MPI_SRC, this routine sets the active communicator used 
!     *  by GCOM and returns the number of processors and my process
!     *  id within this communicator.
!     *  
!     *
!     *
!     * Input:
!     *  VALUE   - value of communicator to use 
!     *
!     * Output:
!     *  ME      - processor ID in this communicator
!     *  NPROC   - number of processors in this communicator
!     *  ISTAT   - status of call, 0 is OK.
!     *
!     * NOTES:       
!     *
!     ******************************************************************

USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__ME,                                                       &
    GC__NPROC

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"
#include "gc_constants.h"
#include "gc_functions.h"

INTEGER (KIND=GC_INT_KIND), INTENT(IN ) :: VALUE
INTEGER (KIND=GC_INT_KIND), INTENT(OUT) :: ME
INTEGER (KIND=GC_INT_KIND), INTENT(OUT) :: NPROC
INTEGER (KIND=GC_INT_KIND), INTENT(OUT) :: ISTAT

INTEGER (KIND=GC_INT_KIND)              :: INFO 


#if defined(MPI_SRC)
GC__MY_MPI_COMM_WORLD = VALUE
CALL MPL_COMM_RANK(GC__MY_MPI_COMM_WORLD, ME, INFO)
CALL MPL_COMM_SIZE(GC__MY_MPI_COMM_WORLD, NPROC, INFO)
#elif defined(SERIAL_SRC)
ME = 0
NPROC = 1
#endif

GC__NPROC = NPROC
GC__ME    = ME
ISTAT     = GC__OK

RETURN
END
