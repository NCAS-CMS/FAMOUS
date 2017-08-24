C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1996, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!    SUBROUTINE ROOT_FRAC---------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE ROOT_FRAC (NSHYD,DZ,ROOTD,F_ROOT)

      IMPLICIT NONE
!
! Description:
!     Calculates the fraction of the total plant roots within each
!     soil layer.
!
!
! Documentation : UM Documentation Paper 25
!
! Current Code Owner : Peter Cox
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.4      9/97     New deck    Peter Cox
!  4.5      6/98     Exponential profile Peter Cox
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!

! Subroutine arguments:
!   Scalar arguments with intent(IN) :
      INTEGER
     & NSHYD                ! IN Number of soil moisture layers.

      REAL
     & DZ(NSHYD)            ! IN Soil layer thicknesses (m).
     &,ROOTD                ! IN Rootdepth (m).

!   Array arguments with intent(OUT) :
      REAL
     & F_ROOT(NSHYD)        ! OUT Fraction of roots in each soil
!                           !     layer.
! Local scalars:
      INTEGER
     & N                    ! WORK Loop counters

      REAL
     & FTOT                 ! WORK Normalisation factor.
     &,ZTOT                 ! WORK Total depth of soil (m).
     &,Z1,Z2                ! WORK Depth of the top and bottom of the
!                           !      soil layers (m).

! Local parameters:
      REAL
     & P                    ! WORK Power describing depth dependence
!                                  of the root density profile.
      PARAMETER (P=2.0)

      Z2=0.0
      ZTOT=0.0

      DO N=1,NSHYD
        Z1=Z2
        Z2=Z2+DZ(N)
        ZTOT=ZTOT+DZ(N)
        F_ROOT(N)=EXP(-P*Z1/ROOTD)-EXP(-P*Z2/ROOTD)
      ENDDO

      FTOT=1.0-EXP(-P*ZTOT/ROOTD)
      DO N=1,NSHYD
        F_ROOT(N)=F_ROOT(N)/FTOT
      ENDDO

      RETURN
      END
