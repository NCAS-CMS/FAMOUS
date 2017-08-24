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
!**********************************************************************
! Routine to calculate the gridbox mean land surface parameters from
! the areal fractions of the surface types and the structural
! properties of the plant functional types.
!
! Written by Peter Cox (June 1997)
!**********************************************************************
      SUBROUTINE SPARM (LAND_FIELD,LAND1,LAND_PTS,TILE_PTS,TILE_INDEX
     &,                 ALBSOIL,FRAC,HT,LAI
     &,                 ALBSNC,ALBSNF,CATCH_T,Z0,Z0_T)

      IMPLICIT NONE

      INTEGER
     + NNVG                       ! Number of non-vegetation surface
C                                 ! types.
     +,NPFT                       ! Number of plant functional types.
     +,NTYPE                      ! Number of surface types.
     +,SOIL                       ! Index of the surface type 'Soil'
      PARAMETER (NNVG=4, NPFT=5, NTYPE=9, SOIL=8)
C                                 ! Land surface types :
C                                 !     1 - Broadleaf Tree
C                                 !     2 - Needleleaf Tree
C                                 !     3 - C3 Grass
C                                 !     4 - C4 Grass
C                                 !     5 - Shrub
C                                 !     6 - Urban
C                                 !     7 - Water
C                                 !     8 - Soil
C                                 !     9 - Ice

      INTEGER
     & LAND_FIELD            ! IN Number of land points in whole grid.
     &,LAND1                 ! IN First land point to be processed.
     &,LAND_PTS              ! IN Number of land points to be processed.
     &,TILE_PTS(NTYPE)              ! IN Number of land points which
!                                   !    include the nth surface type.
     &,TILE_INDEX(LAND_FIELD,NTYPE) ! IN Indices of land points which
!                                   !    include the nth surface type.

      REAL
     & ALBSOIL(LAND_FIELD)        ! IN Soil albedo.
     &,FRAC(LAND_FIELD,NTYPE)     ! IN Fractional cover of each
!                                 !    surface type.
     &,HT(LAND_FIELD,NPFT)        ! IN Vegetation height (m).
     &,LAI(LAND_FIELD,NPFT)       ! IN Leaf area index.

      REAL
     & ALBSNC(LAND_FIELD)         ! OUT Snow-covered albedo.
     &,ALBSNF(LAND_FIELD)         ! OUT Snow-free albedo.
     &,CATCH_T(LAND_FIELD,NTYPE-1)! OUT Canopy capacity for each type
!                                 !     apart from ice (kg/m2).
     &,Z0(LAND_FIELD)             ! OUT Roughness length (m).
     &,Z0_T(LAND_FIELD,NTYPE)     ! OUT Roughness length for each
!                                 !     type (m).
      REAL
     & ALBSNC_T(LAND_FIELD,NTYPE) ! WORK Snow-covered albedo for each
!                                 !      type.
     &,ALBSNF_T(LAND_FIELD,NTYPE) ! WORK Snow-free albedo for each type.
     &,CATCH(LAND_FIELD)          ! WORK Canopy capacity (kg/m2).
     &,FZ0(LAND_FIELD)            ! WORK Aggregation function of Z0.

      INTEGER
     & J,L,N                      ! WORK Loop counters

!-----------------------------------------------------------------------
! Local parameters.
!-----------------------------------------------------------------------
      REAL
     & ALBSNCS                    ! Snow-covered albedo of bare soil.
      PARAMETER (ALBSNCS = 0.8)

C-----------------------------------------------------------------------
C Surface parameters for each Plant Functional Type
C-----------------------------------------------------------------------
      REAL
     + ALBSNC_MAX(NPFT)           ! Snow-covered albedo for large LAI.
     +,ALBSNC_MIN(NPFT)           ! Snow-covered albedo for zero LAI.   
     +,ALBSNF_MAX(NPFT)           ! Snow-free albedo for large LAI.
     +,DZ0V_DH(NPFT)              ! Rate of change of vegetation
C                                 ! roughness length with height.
     +,CATCH0(NPFT)               ! Minimum canopy capacity (kg/m2).
     +,DCATCH_DLAI(NPFT)          ! Rate of change of canopy capacity
C                                 ! with LAI.
     +,INFIL_F(NPFT)              ! Infiltration enhancement factor.
     +,KEXT(NPFT)                 ! Light extinction coefficient.
     +,ROOTD_FT(NPFT)             ! Rootdepth (m).
C----------------------------------------------------------------------
C                           BT    NT   C3G   C4G    S
C----------------------------------------------------------------------
      DATA ALBSNC_MAX  /  0.15, 0.15, 0.60, 0.60, 0.40 /                
      DATA ALBSNC_MIN  /  0.30, 0.30, 0.80, 0.80, 0.80 /                
      DATA ALBSNF_MAX  /  0.10, 0.10, 0.20, 0.20, 0.20 /
      DATA DZ0V_DH     /  0.05, 0.05, 0.10, 0.10, 0.10 /
      DATA CATCH0      /  0.50, 0.50, 0.50, 0.50, 0.50 /
      DATA DCATCH_DLAI /  0.05, 0.05, 0.05, 0.05, 0.05 /
      DATA INFIL_F     /  4.00, 4.00, 2.00, 2.00, 2.00 /
      DATA KEXT        /  0.50, 0.50, 0.50, 0.50, 0.50 /
      DATA ROOTD_FT    /  3.00, 1.00, 0.50, 0.50, 0.50 /                
      REAL
     + ALBSNC_NVG(NNVG)           ! Snow-covered albedo.
     +,ALBSNF_NVG(NNVG)           ! Snow-free albedo.
     +,CATCH_NVG(NNVG)            ! Canopy capacity (kg/m2).
     +,GS_NVG(NNVG)               ! Surface conductance (m/s).
     +,INFIL_NVG(NNVG)            ! Infiltration enhancement factor.
     +,ROOTD_NVG(NNVG)            ! Rootdepth (m).
     +,Z0_NVG(NNVG)               ! Roughness length (m).
C----------------------------------------------------------------------
C                         Urban  Water  Soil   Ice
C----------------------------------------------------------------------
      DATA ALBSNC_NVG  /  0.40,  0.80,  0.80,  0.80 /
      DATA ALBSNF_NVG  /  0.18,  0.06, -1.00,  0.75 /
      DATA CATCH_NVG   /  0.50,  0.00,  0.50,  0.00 /
      DATA GS_NVG      /  5E-3,   1E2,  1E-2,  0.00 /
      DATA INFIL_NVG   /  0.10,  1.00,  0.50,  0.00 /
      DATA ROOTD_NVG   /  0.50,  1.00,  0.10,  0.00 /
      DATA Z0_NVG      /  1.50,  3E-4,  3E-4,  1E-4 /
      REAL
     + LB                         ! Blending height (m).
      PARAMETER (LB = 550.0)

!----------------------------------------------------------------------
! Set parameters for vegetated surface types
!----------------------------------------------------------------------
      DO N=1,NPFT
        CALL PFT_SPARM (LAND_FIELD,N,TILE_INDEX(1,N),TILE_PTS(N)
     &,                 ALBSOIL,HT(1,N),LAI(1,N)
     &,                 ALBSNC_T(1,N),ALBSNF_T(1,N),CATCH_T(1,N)
     &,                 Z0_T(1,N))
      ENDDO

!----------------------------------------------------------------------
! Set parameters for non-vegetated surface types
!----------------------------------------------------------------------
      DO N=NPFT+1,NTYPE
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          ALBSNC_T(L,N) = ALBSNC_NVG(N-NPFT)
          ALBSNF_T(L,N) = ALBSNF_NVG(N-NPFT)
          IF ( ALBSNF_NVG(N-NPFT).LT.0. ) ALBSNF_T(L,N) = ALBSOIL(L)
          Z0_T(L,N) = Z0_NVG(N-NPFT)
        ENDDO
      ENDDO

      DO N=NPFT+1,NTYPE-1
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          CATCH_T(L,N) = CATCH_NVG(N-NPFT)
        ENDDO
      ENDDO

!----------------------------------------------------------------------
! Form area means
!----------------------------------------------------------------------
      DO L=1,LAND_FIELD
        ALBSNC(L) = 0.0
        ALBSNF(L) = 0.0
        CATCH(L) = 0.0
        FZ0(L) = 0.0
      ENDDO

      DO N=1,NTYPE
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          ALBSNC(L) = ALBSNC(L) + FRAC(L,N) * ALBSNC_T(L,N)
          ALBSNF(L) = ALBSNF(L) + FRAC(L,N) * ALBSNF_T(L,N)
          FZ0(L) = FZ0(L) + FRAC(L,N) / (LOG(LB / Z0_T(L,N)))**2
        ENDDO
      ENDDO

      DO N=1,NTYPE-1
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          CATCH(L) = CATCH(L) + FRAC(L,N) * CATCH_T(L,N)
        ENDDO
      ENDDO

!----------------------------------------------------------------------
! Calculate the effective roughness length
!----------------------------------------------------------------------
      DO L=LAND1,LAND1+LAND_PTS-1
        Z0(L) = LB * EXP(-SQRT(1. / FZ0(L)))
      ENDDO

      RETURN
      END
