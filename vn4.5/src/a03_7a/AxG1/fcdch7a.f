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
!!!   SUBROUTINES FCDCH_SEA AND FCDCH_LAND-----------------------------
!!!
!!!   Purpose: Calculate bulk transfer coefficients at one or more
!!!            gridpoints, according to formulae derived by R N B Smith
!!!            October 1989.
!!!
!!!   Model            Modification history:
!!!   version  Date
!!!
!!!   4.4      7/97   Split into separate land and sea routines for the
!!!                   MOSES II tile model (Richard Essery).
!!!
!!!   Programming standard: Unified Model Documentation Paper No 4,
!!!                        Version 2, dated 18/1/90.
!!!
!!!   System component covered: Part of P243.
!!!
!!!   Documentation: UM Documentation Paper No 24, section P243.
!!!                  See especially sub-section (iv).
!!!   -----------------------------------------------------------------

!     SUBROUTINE FCDCH_SEA---------------------------------------------
!
!     Transfer coefficients for sea, sea-ice and leads
!
!     -----------------------------------------------------------------
      SUBROUTINE FCDCH_SEA (P_POINTS,P_FIELD,P1,LAND_MASK,
     &                      RIB,Z0M,Z0H,Z0F,Z1_UV,Z1_TQ,
     &                      CD,CH,LTIMER)

      IMPLICIT NONE

      INTEGER
     & P_POINTS           ! IN Number of gridpoints treated.
     &,P_FIELD            ! IN Size of field on p-grid.
     &,P1                 ! IN First p-point to be treated.

      LOGICAL
     & LTIMER             ! IN logical for TIMER
     &,LAND_MASK(P_FIELD) ! IN .TRUE. for land; .FALSE. elsewhere.

      REAL
     & RIB(P_FIELD)       ! IN Bulk Richardson number.
     &,Z0M(P_FIELD)       ! IN Roughness length for momentum transport
     &,Z0H(P_FIELD)       ! IN Roughness length for heat and moisture
     &,Z0F(P_FIELD)       ! IN Roughness length for free-convective heat
!                         !    and moisture transport (m).
     &,Z1_UV(P_FIELD)     ! IN Height of lowest uv level (m).
     &,Z1_TQ(P_FIELD)     ! IN Height of lowest tq level (m).

      REAL
     & CD(P_FIELD)        ! OUT Surface drag coefficient including form
!                               drag.
     &,CH(P_FIELD)        ! OUT Bulk transfer coefficient for
!                               heat/moisture.

      EXTERNAL TIMER

!----------------------------------------------------------------------
!  Common and local physical constants
C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------


      REAL ALPHAR,HETGEN,CZ,DM
      PARAMETER (
     & ALPHAR=5.0   ! Tunable parameter in FM and FH calculation.
     &,HETGEN=0.0   ! Tunable parameter to represent 'the degree of
!                     heterogeneity' of the surface; must be > or = 0.0
!                     and < or = 1.0
     &,CZ=4.0       ! Tunable parameter in unstable Fh, Fm calculations,
!                     equal to (3h)**-1.5 in the documentation.
     &,DM=2.0       ! Tunable parameter in unstable Fm calculation.
     &)

!  Define local variables (more or less in order of first appearance).

      INTEGER I       ! Loop counter; horizontal field index
      REAL
     & KARMAN2        ! Square of von Karman's constant.
     &,ZETAM          ! See documentation for definition.
     &,ZETAH          ! See documentation for definition.
     &,CDN            ! CD for neutral conditions.
     &,CHN            ! CH for neutral conditions.
     &,PRANDTL        ! Prandtl number at neutrality.
     &,RFZ            ! Temporary in calculation of FM and FH.
     &,RIF            ! Flux Richardson number.
     &,AM             ! Temporary in calculation of FM and FH.
     &,AH             ! Temporary in calculation of FM and FH.
     &,BM             ! Temporary in calculation of FM and FH.
     &,BH             ! Temporary in calculation of FM and FH.
     &,FM             ! Stability factor for CD.
     &,FH             ! Stability factor for CH.

      IF (LTIMER) THEN
        CALL TIMER('FCDCH   ',3)
      ENDIF

      KARMAN2=VKMAN*VKMAN

      DO I=P1,P1+P_POINTS-1
        CD(I) = 0.
        CH(I) = 0.
        IF ( .NOT. LAND_MASK(I) ) THEN

!-----------------------------------------------------------------------
!! 1. Calculate neutral CD, CH.
!-----------------------------------------------------------------------
!  (A) Store ZETAM, ZETAH.
          ZETAM = LOG( (Z1_UV(I) + Z0M(I)) / Z0M(I) )
          ZETAH = LOG( (Z1_TQ(I) + Z0M(I)) / Z0H(I) )
!  (B) Calculate neutral CD, CH.  Eqns P243.40, P243.41
          CDN = KARMAN2 / ( ZETAM * ZETAM )
          CHN = KARMAN2 / ( ZETAH * ZETAM )
          PRANDTL = CDN / CHN
!  (C) Calculate temporary quantities.
          AM = 2.0 * ALPHAR / PRANDTL
          AH = AM

!-----------------------------------------------------------------------
!! 2. Calculate functions Fm, Fh.
!-----------------------------------------------------------------------
          RFZ=0.0
          BM=0.0
          BH=0.0
          RIF = RIB(I) / PRANDTL

!  Case 1: stable boundary layer (RIB > 0).
          IF (RIB(I) .GT. 0.0) THEN
            IF ( 1.0/RIF .GT. HETGEN*ALPHAR ) THEN
              FM = 1.0 - HETGEN * ALPHAR * RIF
              FM = ( FM * FM ) /
     &             ( 1.0 + 2.0 * (1.0-HETGEN) * ALPHAR * RIF )
              FH = FM
            ELSE
              FM = 0.0
              FH = 0.0
            ENDIF

!  Case 2: unstable boundary layer (RIB < or = 0).
          ELSE
!  (A) Store 1/Fz in RFZ.  Eqn P243.51, as approximated by P243.52.
            RFZ = CZ * SQRT ( Z1_UV(I) / Z0F(I) )
!  (B) Store BM and BH.
            BM = DM * AM * CDN * RFZ
            BH = AH * CHN * RFZ
!  (C) Finally calculate FM and FH.
            FM = 1.0 - AM * RIB(I) / ( 1.0 + BM * SQRT(-RIB(I)) )
            FH = 1.0 - AH * RIB(I) / ( 1.0 + BH * SQRT(-RIB(I)) )

          ENDIF

!-----------------------------------------------------------------------
!! 3. Calculate output coefficients.  Eqns P243.53, P243.54.
!-----------------------------------------------------------------------
          CD(I) = CDN * FM
          CH(I) = CHN * FH

        ENDIF ! Sea points

      ENDDO  ! POINTS

      IF (LTIMER) THEN
        CALL TIMER('FCDCH   ',4)
      ENDIF

      RETURN
      END

!     SUBROUTINE FCDCH_LAND---------------------------------------------
!
!     Transfer coefficients for snow, land ice and snow-free land tiles
!
!     ------------------------------------------------------------------
      SUBROUTINE FCDCH_LAND (
     & P_FIELD,LAND_FIELD,TILE_PTS,TILE_INDEX,LAND_INDEX,
     & RIB,WIND_PROFILE_FACTOR,Z0M,Z0H,Z0F,Z1_UV,Z1_TQ,
     & CD,CH,CD_STD,LTIMER
     & )

      IMPLICIT NONE

      INTEGER
     & P_FIELD            ! IN Size of field on p-grid.
     &,LAND_FIELD         ! IN Number of land points.
     &,TILE_PTS           ! IN Number of tile points.
     &,TILE_INDEX(LAND_FIELD)
!                         ! IN Index of tile points.
     &,LAND_INDEX(P_FIELD)! IN Index of land points.

      LOGICAL
     & LTIMER             ! IN Logical for TIMER.

      REAL
     & RIB(LAND_FIELD)    ! IN Bulk Richardson number.
     &,WIND_PROFILE_FACTOR(LAND_FIELD)
!                         ! IN for adjusting the surface transfer
!                         !    coefficients to remove form drag effects.
     &,Z0M(LAND_FIELD)    ! IN Roughness length for momentum transport
     &,Z0H(LAND_FIELD)    ! IN Roughness length for heat and moisture
     &,Z0F(LAND_FIELD)    ! IN Roughness length for free-convective heat
!                         !    and moisture transport (m).
     &,Z1_UV(P_FIELD)     ! IN Height of lowest uv level (m).
     &,Z1_TQ(P_FIELD)     ! IN Height of lowest tq level (m).

      REAL
     & CD(LAND_FIELD)     ! OUT Surface drag coefficient including form
!                         !     drag.
     &,CH(LAND_FIELD)     ! OUT Bulk transfer coefficient for
!                         !     heat/moisture.
     &,CD_STD(LAND_FIELD) ! OUT Surface drag coefficient excluding form
!                               drag.

      EXTERNAL TIMER

!----------------------------------------------------------------------
!  Common and local physical constants
C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------


      REAL ALPHAR,HETGEN,CZ,DM
      PARAMETER (
     & ALPHAR=5.0   ! Tunable parameter in FM and FH calculation.
     &,HETGEN=0.0   ! Tunable parameter to represent 'the degree of
!                     heterogeneity' of the surface; must be > or = 0.0
!                     and < or = 1.0
     &,CZ=4.0       ! Tunable parameter in unstable Fh, Fm calculations,
!                     equal to (3h)**-1.5 in the documentation.
     &,DM=2.0       ! Tunable parameter in unstable Fm calculation.
     &)

!  Define local variables (more or less in order of first appearance).

      INTEGER
     & I              ! Horizontal field index
     &,J              ! Tile field index
     &,L              ! Land field index

      REAL
     & KARMAN2        ! Square of von Karman's constant.
     &,ZETAM          ! See documentation for definition.
     &,ZETAH          ! See documentation for definition.
     &,CDN            ! CD for neutral conditions.
     &,CHN            ! CH for neutral conditions.
     &,CDN_STD        ! CD_STD for neutral conditions.
     &,PRANDTL        ! Prandtl number at neutrality.
     &,RFZ            ! Temporary in calculation of FM and FH.
     &,RIF            ! Flux Richardson number.
     &,AM             ! Temporary in calculation of FM and FH.
     &,AH             ! Temporary in calculation of FM and FH.
     &,BM             ! Temporary in calculation of FM and FH.
     &,BH             ! Temporary in calculation of FM and FH.
     &,BM_STD         ! Temporary in calculation of FM_STD.
     &,FM             ! Stability factor for CD.
     &,FH             ! Stability factor for CH.
     &,FM_STD         ! Stability factor for CD_STD.

      IF (LTIMER) THEN
        CALL TIMER('FCDCH   ',3)
      ENDIF

      KARMAN2=VKMAN*VKMAN

      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        I = LAND_INDEX(L)

!-----------------------------------------------------------------------
!! 1. Calculate neutral CD, CH.
!-----------------------------------------------------------------------
!  (A) Store ZETAM, ZETAH.
          ZETAM = LOG( (Z1_UV(I) + Z0M(L)) / Z0M(L) )
          ZETAH = LOG( (Z1_TQ(I) + Z0M(L)) / Z0H(L) )
!  (B) Calculate neutral CD, CH.  Eqns P243.40, P243.41
          CDN = KARMAN2 / ( ZETAM * ZETAM )
          CHN = KARMAN2 / ( ZETAH * ZETAM ) * WIND_PROFILE_FACTOR(L)
          CDN_STD = CDN * WIND_PROFILE_FACTOR(L) *
     &                    WIND_PROFILE_FACTOR(L)
          PRANDTL = CDN_STD / CHN
!  (C) Calculate temporary quantities.
          AM = 2.0 * ALPHAR / PRANDTL
          AH = AM

!-----------------------------------------------------------------------
!! 2. Calculate functions Fm, Fh.
!-----------------------------------------------------------------------
          RFZ=0.0
          BM=0.0
          BH=0.0
          BM_STD=0.0
          RIF = RIB(L) / PRANDTL

!  Case 1: stable boundary layer (RIB > 0).
          IF (RIB(L) .GT. 0.0) THEN
            IF ( 1.0/RIF .GT. HETGEN*ALPHAR ) THEN
              FM = 1.0 - HETGEN * ALPHAR * RIF
              FM = ( FM * FM ) /
     &             ( 1.0 + 2.0 * (1.0-HETGEN) * ALPHAR * RIF )
              FH = FM
              FM_STD = FM
            ELSE
              FM = 0.0
              FH = 0.0
              FM_STD = 0.0
            ENDIF

!  Case 2: unstable boundary layer (RIB < or = 0).
          ELSE
!  (A) Store 1/Fz in RFZ.  Eqn P243.51, as approximated by P243.52.
            RFZ = CZ * SQRT ( Z1_UV(I) / Z0F(L) )
!  (B) Store BM, BH and BM_STD.
            BM = DM * AM * CDN * RFZ
            BH = AH * CHN * RFZ
            BM_STD = DM * AM * CDN_STD * RFZ
!  (C) Finally calculate FM, FH and FM_STD.
            FM = 1.0 - AM * RIB(L) / ( 1.0 + BM * SQRT(-RIB(L)) )
            FH = 1.0 - AH * RIB(L) / ( 1.0 + BH * SQRT(-RIB(L)) )
            FM_STD = 1.0 - AM * RIB(L) /
     &              ( 1.0 + BM_STD * SQRT(-RIB(L)) )

          ENDIF

!-----------------------------------------------------------------------
!! 3. Calculate output coefficients.  Eqns P243.53, P243.54.
!-----------------------------------------------------------------------
          CD(L) = CDN * FM
          CH(L) = CHN * FH
          CD_STD(L) = CDN_STD * FM_STD

      ENDDO  ! POINTS

      IF (LTIMER) THEN
        CALL TIMER('FCDCH   ',4)
      ENDIF

      RETURN
      END

