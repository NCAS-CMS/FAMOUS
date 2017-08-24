! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"
#include "gc_assorted.h"

SUBROUTINE GC_ABORT (ME, NPROC, MESG)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Aborts program execution on all processors and clean up
!     *  the parallel system.
!     *
!     * Input:
!     *  ME      - my node ID, 0..NPROC-1
!     *  NPROC   - number of nodes
!     *  MESG    - abort message to be printed on stdout
!     *
!     * Output:
!     *
!     * NOTES:       
!     *
!     ******************************************************************

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD
#endif

IMPLICIT NONE

#include "gc_kinds.h"

INTEGER (KIND=GC_INT_KIND) :: ME, NPROC, I, INFO
CHARACTER*(*) MESG

WRITE(GC__FORTERRUNIT,*) 'gc_abort (Processor ',ME,'): ', MESG
CALL GC__FLUSH(6)

#if defined(IBM)
CALL xl__trbk()
#endif

#if defined(MPI_SRC)
CALL MPL_ABORT(GC__MY_MPI_COMM_WORLD , MPIABORT_ERRNO , INFO)
#endif

!---  Should only get through to here if no other abort facility
!---  has been found. Use the C stdlib abort via gc__abort.
CALL GC__ABORT()

RETURN
END

