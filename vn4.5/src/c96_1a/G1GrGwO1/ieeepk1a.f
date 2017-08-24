C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.
C
C Use, duplication or disclosure of this code is subject to the
C restrictions as set forth in the contract.
C
C                Meteorological Office
C                London Road
C                BRACKNELL
C                Berkshire UK
C                RG12 2SZ
C
C If no contract has been raised with this copy of the code, the use,
C duplication or disclosure of it is strictly prohibited.  Permission
C to do so must first be obtained in writing from the Head of Numerical
C Modelling at the above address.
C ******************************COPYRIGHT******************************
C

!LL  SUBROUTINE PACK21
!LL
!LL  Purpose: Packs IEEE 64-bit data into IEEE 32-bit data.
!LL
!LL  Original author: Bob Carruthers, Cray Research
!LL  Revised by     : Paul Burton
!LL
!LL  Model    Date      Modification history
!LL  version
!LL
!LL  4.3      18/03/97  Original version
!LL  4.5      07/09/98   Make pack21 only have four arguments,
!LL                      with pack21_stride having the previously
!LL                      optional fifth argument.
!LL                      Similarly for expand21.
!LL                      Portable error message removed.
!LL                                                   P.Burton

      SUBROUTINE PACK21(N,IN,OUT,NEXP)

! Compress input array 'in' from 64-bit into 'out' in
! 32 bit

      IMPLICIT NONE

!     Arguments:

      INTEGER
     &  N      ! IN: number of floating point words to convert
     &, NEXP   ! IN: present for compatibility - ignored

! The types of the argument are in Fortran 90 format.
! Fortran 77 format is included in comments after it, so can
! be selected by writing a simple modset

      REAL (KIND=8)
!     REAL*8  - Fortran 77 version
     &  IN(N)  ! IN: input array of 64 bit numbers

      REAL (KIND=4)
!     REAL*4  - Fortran 77 version
     &  OUT(N) ! OUT: output array of 32 bit numbers

!     Local variables

      INTEGER
     &  INC    ! increment=1 to be passed to PACK21_STRIDE

      PARAMETER
     & (INC=1)

      CALL PACK21_STRIDE(N,IN,OUT,NEXP,INC)

      RETURN

      END

      SUBROUTINE PACK21_STRIDE(N,IN,OUT,NEXP,INC)

! Compress input array 'in' from 64-bit into 'out' in
! 32 bit

      IMPLICIT NONE

!     Arguments:

      INTEGER
     &  N      ! IN: number of floating point words to convert
     &, NEXP   ! IN: present for compatibility - ignored
     &, INC    ! IN: stride through input array

! The types of the argument are in Fortran 90 format.
! Fortran 77 format is included in comments after it, so can
! be selected by writing a simple modset

      REAL (KIND=8)
!     REAL*8  - Fortran 77 version
     &  IN(N)  ! IN: input array of 64 bit numbers

      REAL (KIND=4)
!     REAL*4  - Fortran 77 version
     &  OUT(N) ! OUT: output array of 32 bit numbers

!     Local variables

      INTEGER
     &  I,J  ! loop indexes to output and input arrays


      J=1
      DO I=1,N
        OUT(I)=IN(J)
        J=J+INC
      ENDDO

      RETURN
      END


      SUBROUTINE EXPAND21(N,IN,OUT,NEXP)

! Expands input array 'in' from 32-bit into 'out' in
! 64 bit

      IMPLICIT NONE

!     Arguments:

      INTEGER
     &  N      ! IN: number of floating point words to convert
     &, NEXP   ! IN: present for compatibility - ignored

! The types of the argument are in Fortran 90 format.
! Fortran 77 format is included in comments after it, so can
! be selected by writing a simple modset

      REAL (KIND=4)
!     REAL*4  - Fortran 77 version
     &  IN(N)  ! IN: input array of 32 bit numbers

      REAL (KIND=8)
!     REAL*8  - Fortran 77 version
     &  OUT(N) ! OUT: output array of 64 bit numbers

!     Local variables

      INTEGER
     &  INC    ! increment=1 to be passed to EXPAND21_STRIDE

      PARAMETER
     & (INC=1)

      CALL EXPAND21_STRIDE(N,IN,OUT,NEXP,INC)

      RETURN

      END

      SUBROUTINE EXPAND21_STRIDE(N,IN,OUT,NEXP,INC)

! Compress input array 'in' from 32-bit into 'out' in
! 64 bit

      IMPLICIT NONE

!     Arguments:

      INTEGER
     &  N      ! IN: number of floating point words to convert
     &, NEXP   ! IN: present for compatibility - ignored
     &, INC    ! IN: stride through output array

! The types of the argument are in Fortran 90 format.
! Fortran 77 format is included in comments after it, so can
! be selected by writing a simple modset

      REAL (KIND=4)
!     REAL*4  - Fortran 77 version
     &  IN(N)  ! IN: input array of 32 bit numbers

      REAL (KIND=8)
!     REAL*8  - Fortran 77 version
     &  OUT(N) ! OUT: output array of 64 bit numbers

!     Local variables

      INTEGER
     &  I,J  ! loop indexes to input and output arrays


      J=1
      DO I=1,N
        OUT(J)=IN(I)
        J=J+INC
      ENDDO

      RETURN
      END
