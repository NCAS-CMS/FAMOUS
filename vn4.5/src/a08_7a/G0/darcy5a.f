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
!    SUBROUTINE DARCY--------------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE DARCY (NPNTS,SOIL_PTS,SOIL_INDEX,B,KS,SATHH,
     &                  STHU1,DZ1,STHU2,DZ2,WFLUX
C LOGICAL LTIMER
     +,LTIMER
     +)

      IMPLICIT NONE
!
! Description:
!     Calculates the Darcian fluxes between adjacent soil layers.
!                                                     (Cox, 6/95)
!
!
! Documentation : UM Documentation Paper 25
!
! Current Code Owner : David Gregory
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.1      6/96     New deck.  Peter Cox
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!

! Global variables:

! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER
     & NPNTS                ! IN Number of gridpoints.
     &,SOIL_PTS             ! IN Number of soil points.


!   Array arguments with intent(IN) :
      INTEGER
     & SOIL_INDEX(NPNTS)    ! IN Array of soil points.

      REAL
     & B(NPNTS)             ! IN Clapp-Hornberger exponent.
     &,DZ1                  ! IN Thickness of the upper layer (m).
     &,DZ2                  ! IN Thickness of the lower layer (m).
     &,KS(NPNTS)            ! IN Saturated hydraulic conductivity
!                           !    (kg/m2/s).
     &,SATHH(NPNTS)         ! IN Saturated soil water pressure (m).
     &,STHU1(NPNTS)         ! IN Unfrozen soil moisture content of upper
!                           !    layer as a fraction of saturation.
!
     &,STHU2(NPNTS)         ! IN Unfrozen soil moisture content of lower
!                           !    layer as a fraction of saturation.
C
      LOGICAL LTIMER        ! Logical switch for TIMER diags


!   Array arguments with intent(OUT) :
      REAL
     & WFLUX(NPNTS)         ! OUT The flux of water between layers
!                           !     (kg/m2/s).

! Local scalars:
      INTEGER
     & I,J,N                ! WORK Loop counters.

! Local arrays:
      REAL
     & THETA(NPNTS,2)       ! WORK Fractional saturation of the upper
!                           !      and lower layer respectively.
     &,THETAK(NPNTS)        ! WORK Fractional saturation at the layer
!                           !      boundary.
     &,K(NPNTS)             ! WORK The hydraulic conductivity between
!                           !      layers (kg/m2/s).
     &,PSI(NPNTS,2)         ! WORK The soil water suction of the upper
!                           !      and lower layer respectively (m).
      IF (LTIMER) THEN
        CALL TIMER('DARCY   ',103)
      ENDIF

!-----------------------------------------------------------------------
! Calculate the fractional saturation of the layers
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        THETA(I,1)=STHU1(I)
        THETA(I,2)=STHU2(I)
      ENDDO

!-----------------------------------------------------------------------
! Calculate the soil water suction of the layers.
!-----------------------------------------------------------------------
      DO N=1,2
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          IF (THETA(I,N).LE.0.01) THEN  ! Prevent blow up for dry soil.
            PSI(I,N)=SATHH(I)/(0.01**B(I))
          ELSEIF (THETA(I,N).GT.0.01.AND.THETA(I,N).LE.1.0) THEN
            PSI(I,N)=SATHH(I)/(THETA(I,N)**B(I))
          ELSE
            PSI(I,N)=SATHH(I)
          ENDIF
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Estimate the fractional saturation at the layer boundary by
! interpolating the soil moisture.
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        THETAK(I)=(DZ2*THETA(I,1)+DZ1*THETA(I,2))/(DZ2+DZ1)
      ENDDO
!-----------------------------------------------------------------------
! Calculate the hydraulic conductivities for transport between layers.
!-----------------------------------------------------------------------
      CALL HYD_CON (NPNTS,SOIL_PTS,SOIL_INDEX,B,KS,THETAK,K,LTIMER)

!-----------------------------------------------------------------------
! Calculate the Darcian flux from the upper to the lower layer.
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        WFLUX(I)=K(I)*(2.0*(PSI(I,2)-PSI(I,1))/(DZ2+DZ1)+1)
      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('DARCY   ',104)
      ENDIF

      RETURN
      END
