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
!!!  SUBROUTINES SF_RIB_LAND and SF_RIB_SEA ---------------------------
!!!
!!!  Purpose: Calculate bulk Richardson number for surface layer
!!!
!!!  SJ, RE       <- programmer of some or all of previous code changes
!!!
!!!  ------------------------------------------------------------------

!    SUBROUTINE SF_RIB_LAND--------------------------------------------
!
!    Calculate RIB for land tiles
!
!    ------------------------------------------------------------------
      SUBROUTINE SF_RIB_LAND (
     & P_FIELD,LAND_FIELD,TILE_PTS,LAND_INDEX,TILE_INDEX,
     & BQ_1,BT_1,QSTAR,QW_1,RESFT,TL_1,TSTAR,VSHR,Z0H,Z0M,Z1_TQ,Z1_UV,
     & RIB,LTIMER
     & )

      IMPLICIT NONE

      INTEGER
     & P_FIELD             ! IN Total number of P-grid points.
     &,LAND_FIELD          ! IN Total number of land points.
     &,TILE_PTS            ! IN Number of tile points.
     &,LAND_INDEX(P_FIELD) ! IN Index of land points.
     &,TILE_INDEX(LAND_FIELD)! IN Index of tile points.

      LOGICAL
     & LTIMER              ! IN logical for TIMER

      REAL
     & BQ_1(P_FIELD)       ! IN A buoyancy parameter for lowest atm
!                          !    level. ("beta-q twiddle").
     &,BT_1(P_FIELD)       ! IN A buoyancy parameter for lowest atm
!                          !    level. ("beta-T twiddle").
     &,QSTAR(LAND_FIELD)   ! IN Surface saturated sp humidity.
     &,QW_1(P_FIELD)       ! IN Total water content of lowest
!                          !    atmospheric layer (kg per kg air).
     &,RESFT(LAND_FIELD)   ! IN Total resistance factor.
     &,TL_1(P_FIELD)       ! IN Liquid/frozen water temperature for
!                          !    lowest atmospheric layer (K).
     &,TSTAR(LAND_FIELD)   ! IN Surface temperature (K).
     &,VSHR(P_FIELD)       ! IN Magnitude of surface-to-lowest-level
!                          !    wind shear.
     &,Z0H(LAND_FIELD)     ! IN Roughness length for heat and moisture m
     &,Z0M(LAND_FIELD)     ! IN Effective roughness length for momentum
     &,Z1_TQ(P_FIELD)      ! IN Height of lowest TQ level (m).
     &,Z1_UV(P_FIELD)      ! IN Height of lowest UV level (m).

      REAL
     & RIB(LAND_FIELD)     ! OUT Bulk Richardson number for lowest layer

!  Symbolic constants -----------------------------------------------

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_R_CP-------------------------------------
C R IS GAS CONSTANT FOR DRY AIR
C CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
C PREF IS REFERENCE SURFACE PRESSURE
      REAL R,CP,KAPPA,PREF  ! PREFTOKI

      PARAMETER(R=287.05,
     &          CP=1005.,
     &          KAPPA=R/CP,
     &          PREF=100000.)
C*----------------------------------------------------------------------


!  Workspace --------------------------------------------------------
      INTEGER
     & I                   ! Horizontal field index.
     &,J                   ! Tile field index.
     &,L                   ! Land field index.

      REAL
     & DQ(LAND_FIELD)      ! Sp humidity difference between surface
!                          ! and lowest atmospheric level (Q1 - Q*).
     &,DTEMP(LAND_FIELD)   ! Modified temperature difference between
!                            surface and lowest atmospheric level.

      IF (LTIMER) THEN
        CALL TIMER('SF_RIB  ',3)
      ENDIF

!-----------------------------------------------------------------------
!!  1 Calculate temperature (strictly, liquid/ice static energy) and
!!    humidity jumps across the surface layer.
!-----------------------------------------------------------------------
      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        I = LAND_INDEX(L)
        DTEMP(L) = TL_1(I) - TSTAR(L) + (G/CP)*(Z1_TQ(I)+Z0M(L)-Z0H(L))
!                                                             ! P243.118
        DQ(L) = QW_1(I) - QSTAR(L)                            ! P243.119
      ENDDO

!-----------------------------------------------------------------------
!!  2 Calculate bulk Richardson numbers for the surface layer.
!-----------------------------------------------------------------------
      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        I = LAND_INDEX(L)
        RIB(L) = G*Z1_UV(I)*(BT_1(I)*DTEMP(L) + BQ_1(I)*RESFT(L)*DQ(L))
     &             / ( VSHR(I)*VSHR(I) )                       ! P243.43
      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('SF_RIB  ',4)
      ENDIF

      RETURN
      END

!    SUBROUTINE SF_RIB_SEA---------------------------------------------
!
!    Calculate RIB for sea, sea-ice and sea-ice leads
!
!    ------------------------------------------------------------------
      SUBROUTINE SF_RIB_SEA (
     & P_POINTS,P_FIELD,P1,LAND_MASK,NSICE,SICE_INDEX,
     & BQ_1,BT_1,ICE_FRACT,QSTAR_ICE,QSTAR_SEA,QW_1,TL_1,TSTAR_ICE,
     & TSTAR_SEA,VSHR,Z0H_ICE,Z0H_SEA,Z0M_ICE,Z0M_SEA,Z1_TQ,Z1_UV,
     & RIB_SEA,RIB_ICE,LTIMER
     & )

      IMPLICIT NONE

      INTEGER
     & P_POINTS            ! IN Number of P-grid points to be processed.
     &,P_FIELD             ! IN Total number of P-grid points.
     &,P1                  ! IN First P-point to be processed.
     &,NSICE               ! IN Number of sea-ice points.
     &,SICE_INDEX(P_FIELD) ! IN Index of sea-ice points.

      LOGICAL
     & LTIMER              ! IN logical for TIMER
     &,LAND_MASK(P_FIELD)  ! IN .TRUE. for land; .FALSE. elsewhere. F60.

      REAL
     & BQ_1(P_FIELD)       ! IN A buoyancy parameter for lowest atm
!                          !    level. ("beta-q twiddle").
     &,BT_1(P_FIELD)       ! IN A buoyancy parameter for lowest atm
!                          !    level. ("beta-T twiddle").
     &,ICE_FRACT(P_FIELD)  ! IN Fraction of gridbox which is sea-ice.
     &,QSTAR_ICE(P_FIELD)  ! IN Surface saturated sp humidity over
!                          !    sea-ice.
     &,QSTAR_SEA(P_FIELD)  ! IN Surface saturated sp humidity over
!                          !    sea and sea-ice leads.
     &,QW_1(P_FIELD)       ! IN Total water content of lowest
!                          !    atmospheric layer (kg per kg air).
     &,TL_1(P_FIELD)       ! IN Liquid/frozen water temperature for
!                          !    lowest atmospheric layer (K).
     &,TSTAR_ICE(P_FIELD)  ! IN Surface temperature of sea-ice (K).
     &,TSTAR_SEA(P_FIELD)  ! IN Surface temperature of sea and sea-ice
!                          !    leads (K).
     &,VSHR(P_FIELD)       ! IN Magnitude of surface-to-lowest-level
!                          !    wind shear.
     &,Z0H_ICE(P_FIELD)    ! IN Roughness length for heat and moisture
!                          !    transport over sea-ice (m).
     &,Z0H_SEA(P_FIELD)    ! IN Roughness length for heat and moisture
!                          !    transport over sea or sea-ice leads (m).
     &,Z0M_ICE(P_FIELD)    ! IN Roughness length for momentum over
!                          !    sea-ice (m).
     &,Z0M_SEA(P_FIELD)    ! IN Roughness length for momentum over sea
!                          !    or sea-ice leads (m).
     &,Z1_TQ(P_FIELD)      ! IN Height of lowest TQ level (m).
     &,Z1_UV(P_FIELD)      ! IN Height of lowest UV level (m).

      REAL
     & RIB_SEA(P_FIELD)    ! OUT Bulk Richardson number for lowest layer
!                          !     over sea or sea-ice leads.
     &,RIB_ICE(P_FIELD)    ! OUT Bulk Richardson number for lowest layer
!                          !     over sea-ice.


!  Symbolic constants -----------------------------------------------

C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_R_CP-------------------------------------
C R IS GAS CONSTANT FOR DRY AIR
C CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
C PREF IS REFERENCE SURFACE PRESSURE
      REAL R,CP,KAPPA,PREF  ! PREFTOKI

      PARAMETER(R=287.05,
     &          CP=1005.,
     &          KAPPA=R/CP,
     &          PREF=100000.)
C*----------------------------------------------------------------------


!  Workspace --------------------------------------------------------
      INTEGER
     & I                   ! Horizontal field index.
     &,J                   !Sea-ice field index.
      REAL
     & DQ                  ! Sp humidity difference between surface
!                          ! and lowest atmospheric level (Q1 - Q*).
     &,DTEMP               ! Modified temperature difference between
!                          ! surface and lowest atmospheric level.

      IF (LTIMER) THEN
        CALL TIMER('SF_RIB  ',3)
      ENDIF

      DO I=P1,P1+P_POINTS-1
        IF ( .NOT.LAND_MASK(I) ) THEN
! Sea and sea-ice leads
          DTEMP = TL_1(I) - TSTAR_SEA(I)                   ! P243.118
     &                  + (G/CP)*(Z1_TQ(I) + Z0M_SEA(I) - Z0H_SEA(I))
          DQ = QW_1(I) - QSTAR_SEA(I)                      ! P243.119
          RIB_SEA(I) = G*Z1_UV(I)*( BT_1(I)*DTEMP + BQ_1(I)*DQ ) /
     &                                 ( VSHR(I)*VSHR(I) )
        ENDIF
      ENDDO

      DO J=1,NSICE
        I = SICE_INDEX(J)
! Sea-ice
        DTEMP = TL_1(I) - TSTAR_ICE(I)
     &                  + (G/CP)*(Z1_TQ(I) + Z0M_ICE(I) - Z0H_ICE(I))
        DQ = QW_1(I) - QSTAR_ICE(I)
        RIB_ICE(I) = G*Z1_UV(I)*( BT_1(I)*DTEMP + BQ_1(I)*DQ ) /
     &                                ( VSHR(I) * VSHR(I) )
      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('SF_RIB  ',4)
      ENDIF

      RETURN
      END
