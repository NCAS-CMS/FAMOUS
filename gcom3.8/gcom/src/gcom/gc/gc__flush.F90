! *****************************COPYRIGHT*******************************
! (C) CROWN COPYRIGHT , MET OFFICE, ALL RIGHTS RESERVED.
! PLEASE REFER TO COPYRIGHT FILE IN TOP LEVEL GCOM DIRECTORY
!                 FOR FURTHER DETAILS
! *****************************COPYRIGHT*******************************

!=======================================================================
!  THIS IS AN INTERNAL ROUTINE TO BE USED WITHIN THE GC INTERFACE ONLY.
!  IT IS A WRAPPER FOR THE FLUSH CALL FOUND IN MANY COMPILERS.
!=======================================================================

SUBROUTINE GC__FLUSH(LUNIT)

! If not set we won't use any of the subroutine
#if defined(GC__FLUSHUNIT6)

!     Required if using the NAG compiler
#if defined(LINUX_NAG_COMPILER)
USE F90_UNIX_IO,ONLY:FLUSH
#endif

IMPLICIT NONE
#include "gc_kinds.h"

!     The subroutine's arguments, whatever the compiler
INTEGER (KIND=GC_INT_KIND), INTENT(IN)  :: LUNIT

INTEGER (KIND=GC_INT_KIND)              :: ICODE

!     If on NAG, X1, XD1 or XT3 require two 32 bit arguments to flush.
#if defined(LINUX_NAG_COMPILER) || defined(_X1) || defined(XD1) \
 || defined(XT3)
INTEGER (KIND=GC_INTEGER32) :: ICODE1
INTEGER (KIND=GC_INTEGER32) :: LUNIT1
LUNIT1 = LUNIT
CALL FLUSH(LUNIT1,ICODE1)
ICODE = ICODE1

!     If on the T3E we require two 64 bit arguments
#elif defined(T3E)
CALL FLUSH(LUNIT,ICODE)

!     All others use one 64 bit argument
#elif defined(IBM)
FLUSH(LUNIT)
#else
CALL FLUSH(LUNIT)
#endif

#endif
END SUBROUTINE GC__FLUSH
