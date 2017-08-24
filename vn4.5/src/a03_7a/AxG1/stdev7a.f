C *****************************COPYRIGHT******************************
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
!
!!!  SUBROUTINES STDEV1_SEA and STDEV1_LAND ----------------------------
!!!
!!!  Purpose: Calculate the standard deviations of layer 1 turbulent
!!!           fluctuations of temperature and humidity using approximate
!!!           formulae from first order closure.
!!!
!!!  -------------------------------------------------------------------
!

!!!  SUBROUTINE STDEV1_SEA ---------------------------------------------
!!!  Layer 1 standard deviations for sea and sea-ice
!!!  -------------------------------------------------------------------
      SUBROUTINE STDEV1_SEA (
     & P_POINTS,P_FIELD,P1,LAND_MASK,
     & BQ_1,BT_1,FQW_1,FTL_1,ICE_FRACT,RHOKM_1,RHOSTAR,VSHR,
     & Z0MSEA,Z0_ICE,Z1_TQ,
     & Q1_SD,T1_SD,LTIMER
     & )

      IMPLICIT NONE

      INTEGER
     & P_POINTS              ! IN Number of points to be processed.
     &,P_FIELD               ! IN Total number points.
     &,P1                    ! IN First point to be processed.

      LOGICAL
     & LTIMER                ! IN logical for TIMER
     &,LAND_MASK(P_FIELD)    ! IN .TRUE. for land

      REAL
     & BQ_1(P_FIELD)         ! IN Buoyancy parameter.
     &,BT_1(P_FIELD)         ! IN Buoyancy parameter.
     &,FQW_1(P_FIELD)        ! IN Surface flux of QW.
     &,FTL_1(P_FIELD)        ! IN Surface flux of TL.
     &,ICE_FRACT(P_FIELD)    ! IN Fraction of gridbox which is sea-ice.
     &,RHOKM_1(P_FIELD)      ! IN Surface momentum exchange coefficient.
     &,RHOSTAR(P_FIELD)      ! IN Surface air density.
     &,VSHR(P_FIELD)         ! IN Magnitude of surface-to-lowest-level
!                            !    wind shear.
     &,Z0MSEA(P_FIELD)       ! IN Sea roughness length.
     &,Z0_ICE(P_FIELD)       ! IN Sea-ice roughness length.
     &,Z1_TQ(P_FIELD)        ! IN Height of lowest tq level.

      REAL
     & Q1_SD(P_FIELD)        ! OUT Standard deviation of turbulent
!                            !     fluctuations of surface layer
!                            !     specific humidity (kg/kg).
     &,T1_SD(P_FIELD)        ! OUT Standard deviation of turbulent
!                            !     fluctuations of surface layer
!                            !     temperature (K).

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------

!  Workspace --------------------------------------------------------
      INTEGER
     & I                     ! Loop counter (horizontal field index).
      REAL
     & VS                    ! Surface layer friction velocity
     &,VSF1_CUBED            ! Cube of surface layer free convective
!                            ! scaling velocity
     &,WS1                   ! Turbulent velocity scale for surface
!                            ! layer
     &,Z0                    ! Roughness length

      IF (LTIMER) THEN
        CALL TIMER('STDEV1  ',3)
      ENDIF

      DO I=P1,P1+P_POINTS-1
        IF ( .NOT.LAND_MASK(I) ) THEN

          Z0 = Z0MSEA(I)
          IF ( ICE_FRACT(I) .GT. 0. ) Z0 = Z0_ICE(I)
          VS = SQRT ( RHOKM_1(I)/RHOSTAR(I) * VSHR(I) )
          VSF1_CUBED = 1.25*G*(Z1_TQ(I) + Z0) *
     &                ( BT_1(I)*FTL_1(I) + BQ_1(I)*FQW_1(I) )/RHOSTAR(I)
          IF ( VSF1_CUBED .GT. 0.0 ) THEN
            WS1 = ( VSF1_CUBED + VS*VS*VS ) ** (1.0/3.0)
            T1_SD(I) = MAX ( 0.0 , 1.93*FTL_1(I) / (RHOSTAR(I)*WS1) )
            Q1_SD(I) = MAX ( 0.0 , 1.93*FQW_1(I) / (RHOSTAR(I)*WS1) )
          ENDIF

        ENDIF
      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('STDEV1  ',4)
      ENDIF

      RETURN
      END

!!!  SUBROUTINE STDEV1_LAND --------------------------------------------
!!!  Layer 1 standard deviations for land tiles
!!!  -------------------------------------------------------------------
      SUBROUTINE STDEV1_LAND (
     & P_FIELD,LAND_FIELD,TILE_PTS,LAND_INDEX,TILE_INDEX,
     & BQ_1,BT_1,FQW_1,FTL_1,RHOKM_1,RHOSTAR,VSHR,Z0M,Z1_TQ,
     & Q1_SD,T1_SD,LTIMER
     & )

      IMPLICIT NONE

      INTEGER
     & P_FIELD               ! IN Total number of P-grid points.
     &,LAND_FIELD            ! IN Total number of land points.
     &,TILE_PTS              ! IN Number of tile points.
     &,LAND_INDEX(P_FIELD)   ! IN Index of land points.
     &,TILE_INDEX(LAND_FIELD)! IN Index of tile points.

      LOGICAL
     & LTIMER                ! IN logical for TIMER

      REAL
     & BQ_1(P_FIELD)         ! IN Buoyancy parameter.
     &,BT_1(P_FIELD)         ! IN Buoyancy parameter.
     &,FQW_1(LAND_FIELD)     ! IN Surface flux of QW.
     &,FTL_1(LAND_FIELD)     ! IN Surface flux of TL.
     &,RHOKM_1(LAND_FIELD)   ! IN Surface momentum exchange coefficient.
     &,RHOSTAR(P_FIELD)      ! IN Surface air density.
     &,VSHR(P_FIELD)         ! IN Magnitude of surface-to-lowest-level
!                            !    wind shear.
     &,Z0M(LAND_FIELD)       ! IN Roughness length for momentum.
     &,Z1_TQ(P_FIELD)        ! IN Height of lowest tq level.

      REAL
     & Q1_SD(P_FIELD)        ! INOUT Standard deviation of turbulent
!                            !       fluctuations of surface layer
!                            !       specific humidity (kg/kg).
     &,T1_SD(P_FIELD)        ! INOUT Standard deviation of turbulent
!                            !       fluctuations of surface layer
!                            !       temperature (K).

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------

!  Workspace --------------------------------------------------------
      INTEGER
     & I                     ! Horizontal field index.
     &,J                     ! Tile index.
     &,L                     ! Land field inde.
      REAL
     & VS                    ! Surface layer friction velocity
     &,VSF1_CUBED            ! Cube of surface layer free convective
!                            ! scaling velocity
     &,WS1                   ! Turbulent velocity scale for surface
!                            ! layer

      IF (LTIMER) THEN
        CALL TIMER('STDEV1  ',3)
      ENDIF

      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        I = LAND_INDEX(L)

        VS = SQRT ( RHOKM_1(L)/RHOSTAR(I) * VSHR(I) )
        VSF1_CUBED = 1.25*G*(Z1_TQ(I) + Z0M(L)) *
     &             ( BT_1(I)*FTL_1(L) + BQ_1(I)*FQW_1(L) )/RHOSTAR(I)
        IF ( VSF1_CUBED .GT. 0.0 ) THEN
          WS1 = ( VSF1_CUBED + VS*VS*VS ) ** (1.0/3.0)
          T1_SD(I) = MAX ( 0.0 , T1_SD(I),
     &                     1.93*FTL_1(L) / (RHOSTAR(I)*WS1) )
          Q1_SD(I) = MAX ( 0.0 , Q1_SD(I),
     &                     1.93*FQW_1(L) / (RHOSTAR(I)*WS1) )
        ELSE
          T1_SD(I) = MAX(T1_SD(I), 0.0)
          Q1_SD(I) = MAX(Q1_SD(I), 0.0)
        ENDIF

      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('STDEV1  ',4)
      ENDIF

      RETURN
      END
