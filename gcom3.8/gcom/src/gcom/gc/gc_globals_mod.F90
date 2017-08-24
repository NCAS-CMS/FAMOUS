! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

MODULE GC_GLOBALS_MOD
!     *****************************************************************
!     * Purpose:
!     * Global variables - replacing many common blocks
!     *****************************************************************
#if defined(MPI_SRC)
USE MPL, ONLY :                                                   &
    MPL_ADDRESS_KIND
#endif

USE GCOM_MOD, ONLY :                                              &
    GC_INT_KIND,                                                  &
    GC_LOG_KIND

IMPLICIT NONE

INTEGER (KIND=GC_INT_KIND) :: GC__NPROC       ! Number of processors
INTEGER (KIND=GC_INT_KIND) :: GC__ME          ! My processor
LOGICAL (KIND=GC_LOG_KIND) :: GC__INITED      ! Has GCOM been initialised?

#if defined(MPI_SRC)
INTEGER (KIND=MPL_ADDRESS_KIND):: GC__MPI_MAXTAG        ! Maximum tag
INTEGER (KIND=GC_INT_KIND)     :: GC__MY_MPI_COMM_WORLD ! Communicator

INTEGER, PARAMETER             :: MAX_MPI_TYPES = 1024  ! Hard-wired
INTEGER (KIND=GC_INT_KIND)     ::                                 &
                             MPI_TYPES_ARRAY(5, MAX_MPI_TYPES)
INTEGER (KIND=GC_INT_KIND)     :: N_MPI_TYPES_DEFINED   ! How many defined
INTEGER (KIND=GC_INT_KIND)     :: MPI_TYPE_SEQ          
#endif

! Maximum number of rotations
INTEGER (KIND=GC_INT_KIND), PARAMETER :: MAX_ROTATE = 8192

! Maximum number of options
INTEGER (KIND=GC_INT_KIND), PARAMETER :: GC__MAX_OPTS = 2

! Options array
INTEGER (KIND=GC_INT_KIND)            :: GC__OPTIONS(GC__MAX_OPTS)

END MODULE GC_GLOBALS_MOD
