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
!    SUBROUTINE HYD_CON------------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE HYD_CON (NPNTS,SOIL_PTS,SOIL_INDEX,B,KS,THETAK,K
C LOGICAL LTIMER
     +,LTIMER
     +)

      IMPLICIT NONE
!
! Description:
!     Calculates the hydraulic conductivity         (Cox, 6/95)
!
!
! Documentation : UM Documentation Paper 25
!
! Current Code Owner : David Gregory
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.1      6/96     New deck.   Peter Cox
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
     & NPNTS            ! IN points in grid
     &,SOIL_PTS         ! IN Number of soil points.

!   Array arguments with intent(IN) :
      INTEGER
     & SOIL_INDEX(NPNTS)! IN Array of soil points.

      REAL
     & B(NPNTS)         ! IN Exponent in conductivity and soil water
!                       !    suction fits.
     &,KS(NPNTS)        ! IN The saturated hydraulic conductivity (kg/m2
     &,THETAK(NPNTS)    ! IN Fractional saturation.
C
      LOGICAL LTIMER    ! Logical switch for TIMER diags

!   Array arguments with intent(OUT) :
      REAL
     & K(NPNTS)         ! OUT The hydraulic conductivity (kg/m2/s).

! Local scalars:
      INTEGER
     & I,J              ! WORK Loop counter.

      IF (LTIMER) THEN
        CALL TIMER('HYDCON  ',103)
      ENDIF

      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)

        IF (THETAK(I).GE.0.0.AND.THETAK(I).LT.1.0) THEN
          K(I)=KS(I)*THETAK(I)**(2*B(I)+3)
        ELSEIF (THETAK(I).LT.0.0) THEN
          K(I)=0.0
        ELSE
          K(I)=KS(I)
        ENDIF

      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('HYDCON  ',104)
      ENDIF

      RETURN
      END
