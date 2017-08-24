! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"

SUBROUTINE GC_INIT (PATH, ME, NPROC)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Initialize all (machine dependent) variables used in the
!     *  communication.
!     *
!     * Input:
!     *  NPROC   - OBSOLETE 
!     *  PATH    - OBSOLETE
!     *                   
!     *                  
!     *
!     * Output:
!     *  NPROC   - number of nodes 
!     *  ME      - my node ID, 0..NPROC-1
!     *
!     * NOTES:       
!     *
!     * Implementation:
!     *  Split into GC_INIT_INTRO(COMM) 
!     *         and GC_INIT_FINAL(ME,NPROC,COMM)
!     *  COMM (INTEGER) is output from INTRO and input into FINAL
!     *  and is the MPI global communicator (undefined for other
!     *  communications systems) - this allows the user to call
!     *  INTRO & FINAL seperately rather than the normal GC_INIT
!     *  and intercept the global communicator and replace with their
!     *  own MPI communicator if required.
!     *            
!     *
!     ******************************************************************

IMPLICIT NONE

#include "gc_kinds.h"

CHARACTER*(*) PATH
INTEGER (KIND=GC_INT_KIND) :: ME, NPROC

INTEGER (KIND=GC_INT_KIND) :: COMM

CALL GC_INIT_INTRO(COMM)
CALL GC_INIT_FINAL(ME,NPROC,COMM)

RETURN
END

SUBROUTINE GC_INIT_INTRO (COMM)

USE MPL, ONLY :                                                   &
    MPL_COMM_WORLD

IMPLICIT NONE

#include "gc_kinds.h"

INTEGER (KIND=GC_INT_KIND) :: COMM ! OUT : Global communicator from MPI (otherwise 0)

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND) :: INFO
LOGICAL (KIND=GC_LOG_KIND) :: FLAG
#endif

COMM=0

#if defined(MPI_SRC)
CALL MPL_INITIALIZED(FLAG,INFO)
! Only initialise MPI if is isn't already done
IF (.NOT. FLAG) THEN
  CALL MPL_INIT(INFO)
END IF

COMM=MPL_COMM_WORLD
#endif

RETURN
END

SUBROUTINE GC_INIT_FINAL (ME,NPROC,COMM)

USE MPL, ONLY :                                                   &
    MPL_TAG_UB,                                                   &
    MPL_COMM_WORLD

USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__ME,                                                       &
    GC__NPROC,                                                    &
    GC__INITED

USE GCOM_MOD ,ONLY :                                              &
    GC_FORCE_BITREP,                                              &
    GC_ALLTOALL_VERSION,                                          &
    GC_OFF,                                                       &
    GC_ALLTOALL_ORIG

#if defined(MPI_SRC)
USE GC_GLOBALS_MOD, ONLY :                                        &
    GC__MY_MPI_COMM_WORLD,                                        &
    GC__MPI_MAXTAG
#endif


IMPLICIT NONE

#include "gc_kinds.h"

INTEGER (KIND=GC_INT_KIND) :: ME, NPROC
INTEGER (KIND=GC_INT_KIND) :: COMM ! IN : Global communicator to
                                   !       use for MPI

#if defined(MPI_SRC)
INTEGER (KIND=GC_INT_KIND)       :: INFO
LOGICAL (KIND=GC_LOG_KIND)       :: MPIFLAG
INTEGER (KIND=GC_INT_KIND), SAVE ::                               &
                            BSEND_BUFFER(MPI_BSEND_BUFFER_SIZE)
#endif

LOGICAL (KIND=GC_LOG_KIND),SAVE  :: INITED
INTEGER (KIND=GC_INT_KIND)       :: I
INTEGER (KIND=GC_INT_KIND)       :: IERR
DATA INITED/.FALSE./


IF (INITED) THEN
   RETURN
ELSE
   !INITED=.TRUE.
   CALL GC__INIT_INTERNALS(INITED)
#if defined(GCOM_TIMER)
   CALL GC__INIT_TIMER()
#endif
ENDIF

#if defined(MPI_SRC)
GC__MY_MPI_COMM_WORLD=COMM

CALL MPL_COMM_RANK(GC__MY_MPI_COMM_WORLD, ME, INFO)
CALL MPL_COMM_SIZE(GC__MY_MPI_COMM_WORLD, NPROC, INFO)

! Warning!! If using the mpi_64bit_fixes then this routine will
! only return the value of MPI_TAG_UB, no matter what is requested.
! This is hardcoded into the the C function in mpi_c_fix.c 

CALL MPL_COMM_GET_ATTR(MPL_COMM_WORLD,MPL_TAG_UB,                 &
                       GC__MPI_MAXTAG,MPIFLAG,INFO)

CALL MPL_BUFFER_ATTACH(BSEND_BUFFER,                              &
                      MPI_BSEND_BUFFER_SIZE*GC__ISIZE, IERR)

#endif

#if defined(SERIAL_SRC)
ME = 0
NPROC = 1
#endif

GC__NPROC = NPROC
GC__ME = ME

! Set default options
! Non-bit reproducible options by default
CALL GC_SETOPT(GC_FORCE_BITREP, GC_OFF, IERR)
! Original RALLTOALLE
CALL GC_SETOPT(GC_ALLTOALL_VERSION, GC_ALLTOALL_ORIG, IERR)

GC__INITED = .TRUE.

IF (GC__ME  ==  0) CALL GC__STAMP()

RETURN
END
