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
! Version 2A of vegetation section: models leaf phenology and vegetation
! competition
!
! Subroutine Interface:
      SUBROUTINE VEG_IC(P_FIELD,FIRST_POINT,LAST_POINT,LAND_FIELD
     &,                 LAND1,LAND_PTS,LAND_INDEX,P_ROWS,ROW_LENGTH
     &,                 EW_Halo,NS_Halo
     &,                 A_STEP,ASTEPS_SINCE_TRIFFID
     &,                 PHENOL_PERIOD,TRIFFID_PERIOD
     &,                 L_PHENOL,L_TRIFFID,L_TRIF_EQ
     &,                 ALB_SOIL,ATIMESTEP,FRAC_DISTURB
     &,                 G_LEAF_AC,G_LEAF_PHEN_AC,NPP_AC
     &,                 RESP_S_AC,RESP_W_AC
     &,                 CS,FRAC,LAI,HT
     &,                 ALBSNC,ALBSNF,CATCH_T,Z0_P,Z0_T
     &,                 C_VEG,CV,LIT_C,LIT_C_MN,G_LEAF_DAY,G_LEAF_PHEN
     &,                 LAI_PHEN,G_LEAF_DR_OUT,NPP_DR_OUT,RESP_W_DR_OUT
     &,                 RESP_S_DR_OUT
     &                  )


      IMPLICIT NONE
!
! Description:
!   Updates Leaf Area Index for Plant Functional Types (PFTs) and uses
!   this to derive new vegetation parameters for PFTs along with gridbox
!   mean values where appropriate.
!
! Method:
!   Calls PHENOL which models phenology and updates Leaf Area Index
!   (LAI), then passes new LAI into SPARM along with canopy height
!   and fractional cover of Plant Functional Types.  SPARM uses this to
!   derive the vegetation parameters for each PFT, and also derives
!   gridbox means where this is required.
!
! Current Code Owner:  Richard Betts
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.4    8/10/97   Original code.  Richard Betts
!   4.5    5/8/98    Pass info on grid and halo dimensions into VEG.
!                    Richard Betts
!   4.5   23/11/98   Output G_LEAF_DAY, G_LEAF_PHEN, LAI_PHEN, 
!                    G_LEAF_DR_OUT, NPP_DR_OUT, RESP_W_DR_OUT and 
!                    RESP_S_DR_OUT as diagnostics.  Richard Betts
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.

      INTEGER
     & P_FIELD                      ! IN Number of P-points.
     &,FIRST_POINT                  ! IN First P-point to be processed.
     &,LAST_POINT                   ! IN Last P-point to be processed.
     &,LAND_FIELD                   ! IN Number of land points.
     &,LAND1                        ! IN First land point to be processe
     &,LAND_PTS                     ! IN Number of land points.
     &,P_ROWS                       ! IN Number of rows on P grid.
     &,ROW_LENGTH                   ! IN Number of P points in a row.
     &,EW_Halo                      ! IN Halo size in the EW direction.
     &,NS_Halo                      ! IN Halo size in the NS direction.
     &,A_STEP                       ! IN Atmospheric timestep number.
     &,ASTEPS_SINCE_TRIFFID         ! INOUT Number of atmospheric
!                                   !       timesteps since last call
!                                   !       to TRIFFID.
     &,PHENOL_PERIOD                ! IN Phenology period (days).
     &,TRIFFID_PERIOD               ! IN TRIFFID period (days).

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
     & LAND_INDEX(LAND_FIELD)       ! IN I=LAND_INDEX(L) => the Ith
C                                   !    point in P_FIELD is the Lth
C                                   !    land point.
      LOGICAL
     & L_PHENOL                     ! IN .T. for interactive leaf
C                                   !    phenology.
     &,L_TRIFFID                    ! IN .T. for interactive vegetation.
     &,L_TRIF_EQ                    ! IN .T. for vegetation equilibrium.

      REAL
     & ALB_SOIL(LAND_FIELD)         ! IN snow-free albedo of soil.
     &,ATIMESTEP                    ! IN Atmospheric timestep (s).
     &,FRAC_DISTURB(LAND_FIELD)     ! IN Fraction of gridbox in which
!                                   !    vegetation is disturbed.
     &,G_LEAF_AC(LAND_FIELD,NPFT)   ! INOUT Accumulated leaf turnover
!                                   !       rate.
     &,G_LEAF_PHEN_AC(LAND_FIELD,NPFT)! INOUT Accumulated leaf turnover
C                                   !       rate including phenology.
     &,NPP_AC(LAND_FIELD,NPFT)      ! INOUT Accumulated NPP (kg C/m2).
     &,RESP_W_AC(LAND_FIELD,NPFT)   ! INOUT Accumulated wood respiration
C                                   !       (kg C/m2).
     &,RESP_S_AC(LAND_FIELD)        ! INOUT Accumulated soil respiration
C                                   !       (kg C/m2).
     &,CS(LAND_FIELD)               ! INOUT Soil carbon content
!                                   !       (kg C/m2).
     &,FRAC(LAND_FIELD,NTYPE)       ! INOUT Fractions of surface types.
     &,LAI(LAND_FIELD,NPFT)         ! INOUT LAI of plant functional
!                                   !       types.
     &,HT(LAND_FIELD,NPFT)          ! INOUT Height of plant functional
C                                   !       types (m).
     &,ALBSNC(LAND_FIELD)           ! OUT Snow-covered albedo.
     &,ALBSNF(LAND_FIELD)           ! OUT Snow-free albedo.
     &,CATCH_T(LAND_FIELD,NTYPE-1)  ! OUT Canopy capacity for each type
C                                   !     aside from ice (kg/m2).
     &,G_LEAF_DAY(LAND_FIELD,NPFT)  ! OUT Mean leaf turnover rate for
!                                   !     input to PHENOL (/360days).
     &,G_LEAF_PHEN(LAND_FIELD,NPFT) ! OUT Mean leaf turnover rate over
!                                   !     phenology period (/360days).
     &,G_LEAF_DR_OUT(LAND_FIELD,NPFT) ! OUT Mean leaf turnover rate for 
!                                   !       driving TRIFFID (/360days).
     &,LAI_PHEN(LAND_FIELD,NPFT)    ! OUT LAI of PFTs after phenology.
     &,NPP_DR_OUT(LAND_FIELD,NPFT)  ! OUT Mean NPP for driving TRIFFID 
!                                   !     (kg C/m2/360days).
     &,RESP_W_DR_OUT(LAND_FIELD,NPFT) ! OUT Mean wood respiration for
!                                   !       driving TRIFFID 
!                                   !       (kg C/m2/360days).
     &,RESP_S_DR_OUT(LAND_FIELD)    ! OUT Mean soil respiration for
!                                   !     driving TRIFFID 
!                                   !     (kg C/m2/360days).
     &,Z0_P(P_FIELD)                ! OUT Effective roughness length
C                                   !     on full grid (m).
     &,Z0_T(LAND_FIELD,NTYPE)       ! OUT Roughness length for each type
C                                   !     (m).
     &,C_VEG(LAND_FIELD,NPFT)       ! OUT Total carbon content of
C                                   !     the vegetation (kg C/m2).
     &,CV(LAND_FIELD)               ! OUT Gridbox mean vegetation
C                                   !     carbon (kg C/m2).
     &,LIT_C(LAND_FIELD,NPFT)       ! OUT Carbon Litter 
!                                   !     (kg C/m2/360days).
     &,LIT_C_MN(LAND_FIELD)         ! OUT Gridbox mean carbon litter
!                                   !     (kg C/m2/360days).


      CALL VEG(P_FIELD,FIRST_POINT,LAST_POINT,LAND_FIELD
     &,        LAND1,LAND_PTS,LAND_INDEX,P_ROWS,ROW_LENGTH
     &,        EW_Halo,NS_Halo
     &,        A_STEP,ASTEPS_SINCE_TRIFFID
     &,        PHENOL_PERIOD,TRIFFID_PERIOD
     &,        L_PHENOL,L_TRIFFID,L_TRIF_EQ
     &,        ALB_SOIL,ATIMESTEP,FRAC_DISTURB
     &,        G_LEAF_AC,G_LEAF_PHEN_AC,NPP_AC
     &,        RESP_S_AC,RESP_W_AC
     &,        CS,FRAC,LAI,HT
     &,        ALBSNC,ALBSNF,CATCH_T,Z0_P,Z0_T
     &,        C_VEG,CV,LIT_C,LIT_C_MN,G_LEAF_DAY,G_LEAF_PHEN
     &,        LAI_PHEN,G_LEAF_DR_OUT,NPP_DR_OUT,RESP_W_DR_OUT
     &,        RESP_S_DR_OUT
     &         )

      RETURN
      END
