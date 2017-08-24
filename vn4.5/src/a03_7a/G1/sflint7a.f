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
CLL  SUBROUTINES SFL_INT_LAND and SFL_INT_SEA--------------------------
CLL
CLL  Purpose: To calculate interpolation coefficients for 10m winds
CLL           and 1.5m temperature/specific humidity diagnostics
CLL           using a generalisation of the method of Geleyn (1988).
CLL
CLL  Suitable for single column use (via *IF definition IBM).
CLL
CLL  Model            Modification history:
CLL version  Date
CLL
CLL   3.4  18/10/94   *DECK inserted into UM version 3.4. S Jackson
CLL
CLL   4.0  30/12/94   Modified calculation of 10m wind interpolation
CLL                   factor when effective roughness length used;
CLL                   10m wind assumed to lie on "local" profile at
CLL                   height z0m+10 metres above the surface.
CLL                                                    R.N.B.Smith
CLL   4.2   Oct. 96   T3E migration - *DEF CRAY removed
CLL                                     S J Swarbrick
CLL
CLL   4.4   Jul. 97   Split into separate land and sea routines for
CLL                   the MOSES II tile model (Richard Essery)
CLL
CLL  Programming standard: Unified Model Documentation Paper No 4,
CLL                        Version 2, dated 18/1/90.
CLL
CLL  Logical component covered: Part of P243.
CLL
CLL  System Task:
CLL
CLL  External Documentation: UMDP No.24
CLL

CLL---------------------------------------------------------------------
CLL  SUBROUTINE SFL_INT_LAND
CLL---------------------------------------------------------------------

      SUBROUTINE SFL_INT_LAND (
     & P_FIELD,LAND_FIELD,TILE_PTS,TILE_INDEX,LAND_INDEX,
     & CD_STD,CD,CH,RIB,TILE_FRAC,WIND_PROFILE_FACTOR,
     & Z0M,Z0M_EFF,Z0H,Z0F,Z1,
     & SU10,SV10,ST1P5,SQ1P5,LTIMER,
     & CDR10M,CHR1P5M
     & )

      IMPLICIT NONE

      INTEGER
     & P_FIELD           ! IN Number of P-grid points.
     &,LAND_FIELD        ! IN Number of land points.
     &,TILE_PTS          ! IN Number of tile points.
     &,TILE_INDEX(LAND_FIELD)
!                        ! IN Index of tile points.
     &,LAND_INDEX(P_FIELD)
!                        ! IN Index of land points.

      REAL
     & CD_STD(LAND_FIELD)! IN Surface drag coefficient for shear stress
!                        !    only, i.e. without orographic part of drag
     &,CD(LAND_FIELD)    ! IN Effective surface drag coefficient,
!                        !    including orographic part of drag
     &,CH(LAND_FIELD)    ! IN Surface transfer coefficient for heat and
!                        !    moisture.
     &,RIB(LAND_FIELD)   ! IN Bulk Richardson number for
!                        !    lowest layer.
     &,TILE_FRAC(LAND_FIELD)
!                        ! IN Tile fraction.
     &,WIND_PROFILE_FACTOR(LAND_FIELD)
!                        ! IN Factor for converting effective friction
!                        !    velocity to local one.
     &,Z0M(LAND_FIELD)   ! IN Roughness length for momentum (m).
     &,Z0M_EFF(LAND_FIELD)
!                        ! IN Effective roughness length for
!                        !    momentum (m).
     &,Z0H(LAND_FIELD)   ! IN Roughness length for heat and
!                        !    moisture (m).
     &,Z0F(LAND_FIELD)   ! IN Roughness length in the free
!                        !    convective limit (m).
     &,Z1(P_FIELD)       ! IN Height of lowest model level (m).

      LOGICAL
     & SU10              ! IN 10m U-wind diagnostic flag
     &,SV10              ! IN 10m V-wind diagnostic flag
     &,ST1P5             ! IN screen temp diagnostic flag
     &,SQ1P5             ! IN screen specific humidity
!                        !    diagnostic flag
     &,LTIMER            ! IN TIMER diagnostics flag

      REAL
     & CDR10M(P_FIELD)   ! INOUT GBM interpolation coeff. for 10m wind
     &,CHR1P5M(LAND_FIELD)
!                        ! OUT Local interpolation coefficient for 1.5m
!                        !     temperature

      EXTERNAL TIMER

! Local and other symbolic constants :-
C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------


      REAL Z1P5M,Z10M
      PARAMETER (
     & Z1P5M = 1.5  ! for diagnosis of screen values of temperature
!                   ! and humidity (assumed to be at 1.5m).
     &,Z10M = 10.0  ! for diagnosis of 10m winds.
     & )

!  Define local storage.

!  (a) Local work arrays.

      REAL
     & Z1E(LAND_FIELD)  ! Level 1 height + Z0M_EFF
     &,SQRTCD(LAND_FIELD) ! Square root of CD

!  (b) Scalars.

      REAL
     & Z1S              ! Level 1 height + Z0M_EFF - Z0H
     &,Z1P5ME           ! Z1P5M + Z0H
     &,Z10ME            ! Z10M + Z0M
     &,SQRTCD_K         ! Temporary storage in calc of 1.5 amd 10m diags
     &,Z_OVER_Z1        ! Temporary storage in calc of 1.5 amd 10m diags
     &,CDNZ             ! Neutral drag coef. for momentum @ 10m
     &,CDNZ1            ! Neutral drag coef. for momentum @ level1
     &,CHNZ             ! Neutral drag coef. for heat/moisture @ 1.5m
     &,CHNZ1            ! Neutral drag coef. for heat/moisture @ level1
     &,CDTEMP1          ! Workspace in calc of interpolation coeffs.
     &,CDTEMP2          ! Workspace in calc of interpolation coeffs.
     &,CDTEMP3          ! Workspace in calc of interpolation coeffs.
     &,CD10             ! Local interpolation coeff. for 10m wind

      INTEGER
     & I       ! Loop counter (horizontal field index).
     &,J       ! Loop counter (tile point index).
     &,L       ! Loop counter (land point field index).

      IF (LTIMER) THEN
        CALL TIMER('SFL_INT   ',3)
      ENDIF

!-----------------------------------------------------------------------
!  If selected calculate interpolation coefficient for 10m winds.
!-----------------------------------------------------------------------

      IF(SU10.OR.SV10) THEN
        DO J=1,TILE_PTS
          L = TILE_INDEX(J)
          I = LAND_INDEX(L)

          Z1E(L) = Z1(I) + Z0M_EFF(L)
          Z10ME = Z10M + Z0M(L)
          CDNZ = LOG( Z10ME / Z0M(L) )
          CDNZ1 = LOG( Z1E(L) / Z0M_EFF(L) )
          SQRTCD(L) = SQRT(CD(L))
          SQRTCD_K = SQRTCD(L) / VKMAN
          Z_OVER_Z1 = Z10M  / Z1(I)

          IF (RIB(L).GE.0.0) THEN

! Stable case

            CD10 = Z_OVER_Z1 + SQRTCD_K *
     &                  (CDNZ - Z_OVER_Z1*CDNZ1)
          ELSE

! Unstable Case

            CDTEMP1 = EXP( 1.0 / SQRTCD_K )
            CDTEMP2 = Z1E(L) / Z0M_EFF(L)
            CDTEMP3 = LOG( ( ( Z1E(L) - Z0M(L) ) * CDTEMP1 -
     &                       ( Z0M_EFF(L) - Z0M(L) ) * CDTEMP2 ) /
     &                     ( ( Z1E(L) - Z10ME ) * CDTEMP1 -
     &                       ( Z0M_EFF(L) - Z10ME ) * CDTEMP2 ) )
            CD10 = SQRTCD_K * ( CDNZ + CDTEMP3 )
          ENDIF
          CDR10M(I) = CDR10M(I) +
     &                          TILE_FRAC(L)*CD10*WIND_PROFILE_FACTOR(L)
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
!   If selected calculate interpolation coefficient for 1.5m screen
!   temperature and specific humidity.
!-----------------------------------------------------------------------

      IF (ST1P5.OR.SQ1P5) THEN
        DO J=1,TILE_PTS
          L = TILE_INDEX(J)
          I = LAND_INDEX(L)

! variables to be used later
          Z1E(L) = Z1(I) + Z0M_EFF(L)
          Z1S = Z1E(L) - Z0H(L)
          Z1P5ME = Z1P5M + Z0H(L)
          CHNZ = LOG( Z1P5ME / Z0H(L) )
          CHNZ1 = LOG( Z1E(L) / Z0H(L) )
          SQRTCD(L) = SQRT(CD_STD(L))
          SQRTCD_K =0.0
          IF (SQRTCD(L) .GT. 0.0) SQRTCD_K = CH(L) / (SQRTCD(L) * VKMAN)
          Z_OVER_Z1 =  Z1P5M  / Z1S

! Stable case

          IF (RIB(L).GE.0.0) THEN
            CHR1P5M(L) = Z_OVER_Z1 + SQRTCD_K *
     &                  (CHNZ - Z_OVER_Z1 * CHNZ1)
          ELSE

! Unstable Case

            CDTEMP1 = EXP(1.0 / SQRTCD_K)
            CDTEMP2 = ( Z1P5M * Z1E(L) +
     &                  Z0H(L) * ( Z1S - Z1P5M ) * CDTEMP1 ) /
     &                ( Z1S * Z1P5ME )
            CHR1P5M(L) = 1.0 - SQRTCD_K * LOG ( CDTEMP2 )
          ENDIF

        ENDDO
      ENDIF

      IF (LTIMER) THEN
        CALL TIMER('SFL_INT ',4)
      ENDIF

      RETURN
      END

CLL---------------------------------------------------------------------
CLL  SUBROUTINE SFL_INT_SEA
CLL---------------------------------------------------------------------

      SUBROUTINE SFL_INT_SEA (
     & P_POINTS,P_FIELD,P1,CD,CH,RIB,Z0M,Z0H,Z0F,Z1,
     & LAND_MASK,SU10,SV10,ST1P5,SQ1P5,LTIMER,
     & CDR10M,CHR1P5M
     & )

      IMPLICIT NONE

      INTEGER
     & P_POINTS           ! IN Number of points to be processed.
     &,P_FIELD            ! IN Total number points.
     &,P1                 ! IN First point to be processed.

      REAL
     & CD(P_FIELD)        ! IN Effective surface drag coefficient,
!                         !    including orographic part of drag
     &,CH(P_FIELD)        ! IN Surface transfer coefficient for heat
!                         !    and moisture.
     &,RIB(P_FIELD)       ! IN Bulk Richardson number for
!                         !    lowest layer.
     &,Z0M(P_FIELD)       ! IN Roughness length for momentum (m).
     &,Z0H(P_FIELD)       ! IN Roughness length for heat and
!                         !    moisture (m).
     &,Z0F(P_FIELD)       ! IN Roughness length in the free
!                         !    convective limit (m).
     &,Z1(P_FIELD)        ! IN Height of lowest model level (m).

      LOGICAL
     & LAND_MASK(P_FIELD) ! IN T for land points, F otherwise.
     &,SU10               ! IN 10m U-wind diagnostic flag
     &,SV10               ! IN 10m V-wind diagnostic flag
     &,ST1P5              ! IN screen temp diagnostic flag
     &,SQ1P5              ! IN screen specific humidity
!                         !    diagnostic flag
     &,LTIMER             ! IN TIMER diagnostics flag

      REAL
     & CDR10M(P_FIELD)    ! OUT interpolation coefficicent for 10m wind
     &,CHR1P5M(P_FIELD)   ! OUT Interpolation coefficient for 1.5m
!                         !     temperature

      EXTERNAL TIMER

!    Local and other symbolic constants :-
C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------

      REAL Z1P5M,Z10M
      PARAMETER (
     & Z1P5M = 1.5  ! for diagnosis of screen values of temperature
!                   ! and humidity (assumed to be at 1.5m).
     &,Z10M = 10.0  ! for diagnosis of 10m winds.
     & )

!  Define local storage.

!  (a) Local work arrays.

      REAL
     & Z1E(P_FIELD)     ! Level 1 height + Z0M
     &,SQRTCD(P_FIELD)  ! Square root of CD

!  (b) Scalars.

      REAL
     & Z1S              ! Level 1 height + Z0M - Z0H
     &,Z1P5ME           ! Z1P5M + Z0H
     &,Z10ME            ! Z10M + Z0M
     &,SQRTCD_K         ! Temporary storage in calc of 1.5 amd 10m diags
     &,Z_OVER_Z1        ! Temporary storage in calc of 1.5 amd 10m diags
     &,CDNZ             ! Neutral drag coef. for momentum @ 10m
     &,CDNZ1            ! Neutral drag coef. for momentum @ level1
     &,CHNZ             ! Neutral drag coef. for heat/moisture @ 1.5m
     &,CHNZ1            ! Neutral drag coef. for heat/moisture @ level1
     &,CDTEMP1          ! Workspace in calc of interpolation coeffs.
     &,CDTEMP2          ! Workspace in calc of interpolation coeffs.
     &,CDTEMP3          ! Workspace in calc of interpolation coeffs.

      INTEGER
     & I       ! Loop counter (horizontal field index).

      IF (LTIMER) THEN
        CALL TIMER('SFL_INT   ',3)
      ENDIF

      DO I=P1,P1+P_POINTS-1
        SQRTCD(I) = SQRT(CD(I))
        Z1E(I) = Z1(I) + Z0M(I)
      ENDDO

!-----------------------------------------------------------------------
!   If selected calculate interpolation coefficient for 10m winds.
!-----------------------------------------------------------------------

      IF(SU10.OR.SV10) THEN
        DO I=P1,P1+P_POINTS-1
          IF ( .NOT. LAND_MASK(I) ) THEN

            Z10ME = Z10M + Z0M(I)
            CDNZ = LOG( Z10ME / Z0M(I) )
            CDNZ1 = LOG( Z1E(I) / Z0M(I) )
            SQRTCD_K = SQRTCD(I) / VKMAN
            Z_OVER_Z1 = Z10M  / Z1(I)

            IF (RIB(I).GE.0.0) THEN

! Stable case

              CDR10M(I) = Z_OVER_Z1 + SQRTCD_K *
     &                    (CDNZ - Z_OVER_Z1*CDNZ1)
            ELSE

! Unstable Case

              CDTEMP1 = EXP( 1.0 / SQRTCD_K )
              CDTEMP2 = Z1E(I) / Z0M(I)
              CDTEMP3 = LOG( Z1(I)* CDTEMP1 /
     &                       ((Z1(I) - Z10M)*CDTEMP1 + Z10ME*CDTEMP2) )
              CDR10M(I) = SQRTCD_K * ( CDNZ + CDTEMP3 )
            ENDIF

          ENDIF
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
!   If selected calculate interpolation coefficient for 1.5m screen
!   temperature and specific humidity.
!-----------------------------------------------------------------------

      IF (ST1P5.OR.SQ1P5) THEN
        DO I=P1,P1+P_POINTS-1
          IF ( .NOT. LAND_MASK(I) ) THEN

! variables to be used later
            Z1S = Z1E(I) - Z0H(I)
            Z1P5ME = Z1P5M + Z0H(I)
            CHNZ = LOG( Z1P5ME / Z0H(I) )
            CHNZ1 = LOG( Z1E(I) / Z0H(I) )
            SQRTCD(I) = SQRT(CD(I))
            SQRTCD_K =0.0
            IF (SQRTCD(I) .GT. 0.0)
     &        SQRTCD_K = CH(I) / (SQRTCD(I) * VKMAN)
            Z_OVER_Z1 =  Z1P5M  / Z1S

! Stable case

            IF (RIB(I).GE.0.0) THEN
              CHR1P5M(I) = Z_OVER_Z1 + SQRTCD_K *
     &                    (CHNZ - Z_OVER_Z1 * CHNZ1)
            ELSE

! Unstable Case

              CDTEMP1 = EXP(1.0 / SQRTCD_K)
              CDTEMP2 = ( Z1P5M * Z1E(I) +
     &                   Z0H(I) * ( Z1S - Z1P5M ) * CDTEMP1 ) /
     &                  ( Z1S * Z1P5ME )
              CHR1P5M(I) = 1.0 - SQRTCD_K * LOG ( CDTEMP2 )
            ENDIF

          ENDIF
        ENDDO
      ENDIF

      IF (LTIMER) THEN
        CALL TIMER('SFL_INT ',4)
      ENDIF

      RETURN
      END
