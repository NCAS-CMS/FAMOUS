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
! Routine to calculate the land surface parameters of a given PFT from
! its areal fraction and structural properties.
C
! Written by Peter Cox (June 1997)
C**********************************************************************
      SUBROUTINE PFT_SPARM  (LAND_FIELD,N,TILE_INDEX,TILE_PTS
     &,                      ALBSOIL,HT,LAI
     &,                      ALBSNC_T,ALBSNF_T,CATCH_T
     &,                      Z0_T)


      IMPLICIT NONE

      INTEGER
     & LAND_FIELD                 ! IN Number of land points.
     &,N                          ! IN Plant functional type.
     &,TILE_PTS                   ! IN Number of land points which
!                                 !    include the surface type.
     &,TILE_INDEX(LAND_FIELD)     ! IN Indices of land points which
!                                 !    include the surface type.
     &,J,L                        ! WORK Loop counters.

      REAL
     & ALBSOIL(LAND_FIELD)        ! IN Soil albedo.
     &,HT(LAND_FIELD)             ! IN Vegetation height (m).
     &,LAI(LAND_FIELD)            ! IN Leaf area index.
     &,ALBSNC_T(LAND_FIELD)       ! OUT Snow-covered albedo.
     &,ALBSNF_T(LAND_FIELD)       ! OUT Snow-free albedo.
     &,CATCH_T(LAND_FIELD)        ! OUT Canopy capacity (kg/m2).
     &,Z0_T(LAND_FIELD)           ! OUT Roughness length (m).
     &,FLIT                       ! WORK Weighting factor for albedo.

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

      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        FLIT = 1.0 - EXP(-KEXT(N) * LAI(L))
        ALBSNC_T(L) = ALBSNC_MIN(N) * (1 - FLIT)
     &              + ALBSNC_MAX(N) * FLIT
        ALBSNF_T(L) = ALBSOIL(L) * (1 - FLIT) + ALBSNF_MAX(N) * FLIT
        Z0_T(L) = DZ0V_DH(N) * HT(L)
        CATCH_T(L) = CATCH0(N) + DCATCH_DLAI(N) * LAI(L)
      ENDDO


      RETURN
      END
