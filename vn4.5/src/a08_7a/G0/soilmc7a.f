C *****************************COPYRIGHT*******************************
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
!    SUBROUTINE SOILMC-------------------------------------------------

      SUBROUTINE SOILMC ( NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX,
     &                    DZ,STHU,V_SAT,V_WILT,SMC )

      IMPLICIT NONE
!
! Description:
!     Diagnoses the soil moisture in a layer at the surface
!
      INTEGER
     & NPNTS                ! IN Number of gridpoints.
     &,NSHYD                ! IN Number of soil moisture levels.
     &,SOIL_PTS             ! IN Number of soil points.
     &,SOIL_INDEX(NPNTS)    ! IN Array of soil points.

      REAL
     & DZ(NSHYD)            ! IN Thicknesses of the soil layers (m).
     &,STHU(NPNTS,NSHYD)    ! IN Unfrozen soil moisture content of
!                           !    each layer as a frac. of saturation.
     &,V_SAT(NPNTS)         ! IN Volumetric soil moisture conc. at
!                           !    saturation (m3 H2O/m3 soil).
     &,V_WILT(NPNTS)        ! IN Volumetric soil moisture conc. below
!                           !    which stomata close (m3 H2O/m3 soil).

      REAL
     & SMC(NPNTS)           ! OUT Soil moisture (kg/m2).

      REAL
     & Z1,Z2                ! WORK Depth of the top and bottom of the
!                           !      soil layers (m).
     &,ZSMC                 ! WORK Depth of layer for soil moisture
!                           !      diagnostic (m).
      PARAMETER ( ZSMC = 1. )

      INTEGER
     & I,J,N                ! WORK Loop counters

C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------

      DO I=1,NPNTS
        SMC(I) = 0.
      ENDDO

      Z2 = 0.
      DO N=1,NSHYD
        Z1 = Z2
        Z2 = Z2 + DZ(N)
        IF ( Z2.LT.ZSMC ) THEN
          DO J=1,SOIL_PTS
            I = SOIL_INDEX(J)
            SMC(I) = SMC(I) + RHO_WATER * DZ(N) *
     &                               ( STHU(I,N)*V_SAT(I) - V_WILT(I) )
          ENDDO
        ELSEIF ( Z2.GE.ZSMC .AND. Z1.LT.ZSMC ) THEN
          DO J=1,SOIL_PTS
            I = SOIL_INDEX(J)
            SMC(I) = SMC(I) + RHO_WATER * ( Z2 - ZSMC ) *
     &                               ( STHU(I,N)*V_SAT(I) - V_WILT(I) )
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END
