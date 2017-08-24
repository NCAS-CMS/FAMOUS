! GC.precision
!   Architecture dependent flags for word length
!   Cray defaults to 64bit REAL, 64bit INTEGER
!   Otherwise default is 32bit REAL, 32bit INTEGER

#include "gc_set_precision.h"

#if defined(PREC_32B)
! Using 32bit INTEGERs & REALs
#define GC_INT_TYPE "32bit INTEGERs"
#define GC_REAL_TYPE "32bit REALs"
#else
! Using 64bit INTEGERs & REALs
#define GC_INT_TYPE "64bit INTEGERs"
#define GC_REAL_TYPE "64bit REALs"
#endif

! And some (possibly outdated) derived quantities
#define GC__STREAM_PAD 1

!============================================================
