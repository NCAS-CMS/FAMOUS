! GC.limits
!   GCOM dimensioning. The user may change this by preprocessor
!   directives on the command line when compiling GCOM. The
!   default setting corresponds to:
#if defined(MPI_SRC)
!    -DMPI_BSEND_BUFFER_SIZE=160000
#if !defined(MPI_BSEND_BUFFER_SIZE)
#define MPI_BSEND_BUFFER_SIZE 160000
#endif
! Size of MPI Send buffer is MPI_BSEND_BUFFER_SIZE
#endif
!============================================================
