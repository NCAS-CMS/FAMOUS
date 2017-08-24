! GC_assorted
!   Collection of flags/constants used within GCOM that do not natuarally
!   fall within the scope of any other include files

#if !defined(GC__FORTERRUNIT)
#define GC__FORTERRUNIT *
#endif
! Error messages are written to unit GC__FORTERRUNIT

#if defined(MPI_SRC)
#if !defined(MPIABORT_ERRNO)
#define MPIABORT_ERRNO 9
#endif
! GC_ABORT will cause exit code MPIABORT_ERRNO to be output
#endif

!============================================================
