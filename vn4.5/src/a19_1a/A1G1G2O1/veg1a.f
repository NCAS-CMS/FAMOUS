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
! Version 1A of vegetation section: models leaf phenology
!
! Subroutine Interface:
      SUBROUTINE VEG(P_FIELD,FIRST_POINT,LAST_POINT,LAND_FIELD
     &,              LAND1,LAND_PTS,LAND_INDEX,P_ROWS,ROW_LENGTH
     &,              EW_Halo,NS_Halo
     &,              A_STEP,PHENOL_PERIOD,L_PHENOL
     &,              ALB_SOIL,ATIMESTEP
     &,              G_LEAF_AC,FRAC,LAI,HT
     &,              ALBSNC,ALBSNF,CATCH_T,Z0_P,Z0_T
     &,              G_LEAF_DAY,G_LEAF_PHEN,LAI_PHEN
     &               )


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
!   4.5   16/09/98   Call SWAPB_LAND to update halo regions of input
!                    fields.   Richard Betts
!   4.5   23/11/98   Output G_LEAF_DAY, G_LEAF_PHEN and LAI_PHEN as
!                    diagnostics.  Richard Betts
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.


      INTEGER
     & P_FIELD               ! IN Number of P-points in whole grid.
     &,FIRST_POINT           ! IN First P-point to be processed.
     &,LAST_POINT            ! IN Number of P-points to be processed.
     &,LAND_FIELD            ! IN Number of land points.
     &,LAND1                 ! IN First land point to be processed.
     &,LAND_PTS              ! IN Number of land points to be processed.
     &,P_ROWS                       ! IN Number of rows on P grid.
     &,ROW_LENGTH                   ! IN Number of P points in a row.
     &,EW_Halo                      ! IN Halo size in the EW direction.
     &,NS_Halo                      ! IN Halo size in the NS direction.
     &,A_STEP                ! IN Atmospheric timestep number.
     &,PHENOL_PERIOD         ! IN Phenology period (days).

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
!                                   !    point in P_FIELD is the Lth
!                                   !    land point.

      INTEGER
     & I,J,L,N                      ! WORK loop counters.

      LOGICAL
     & L_PHENOL                     ! IN .T. for interactive leaf
!                                   !    phenology.
      REAL
     & ALB_SOIL(LAND_FIELD)         ! IN snow-free albedo of soil.
     &,ATIMESTEP                    ! IN Atmospheric timestep (s).
     &,G_LEAF_AC(LAND_FIELD,NPFT)   ! INOUT Accumulated leaf turnover
!                                   !       rate.
     &,FRAC(LAND_FIELD,NTYPE)       ! INOUT Fractions of surface types.
     &,LAI(LAND_FIELD,NPFT)         ! INOUT LAI of plant functional
!                                   !       types.
     &,HT(LAND_FIELD,NPFT)          ! INOUT Height of plant functional
!                                   !       types (m).
     &,ALBSNC(LAND_FIELD)           ! OUT Snow-covered albedo.
     &,ALBSNF(LAND_FIELD)           ! OUT Snow-free albedo.
     &,CATCH_T(LAND_FIELD,NTYPE-1)  ! OUT Canopy capacity for each type
!                                   !     aside from ice (kg/m2).
     &,LAI_PHEN(LAND_FIELD,NPFT)    ! OUT LAI of PFTs after phenology.
!                                   !     Required as separate variable
!                                   !     for top-level argument list
!                                   !     matching with VEG_IC2A.
     &,Z0_P(P_FIELD)                ! OUT Effective roughness length
!                                   !     on full grid (m).
     &,Z0_T(LAND_FIELD,NTYPE)       ! OUT Roughness length for each type
!                                   !     (m).

      INTEGER
     & NSTEP_PHEN                   ! WORK Number of atmospheric
!                                   !      timesteps between calls to
!                                   !      PHENOL.
     &,TILE_PTS(NTYPE)              ! WORK Number of land points which
!                                   !      include the nth surface type.
     &,TILE_INDEX(LAND_FIELD,NTYPE) ! WORK Indices of land points which
!                                   !      include the nth surface type.

      REAL
     & DTIME_PHEN                   ! WORK The phenology timestep (yr).
     &,G_LEAF_DAY(LAND_FIELD,NPFT)  ! WORK Mean leaf turnover rate for
!                                   !      input to PHENOL (/360days).
     &,G_LEAF_PHEN(LAND_FIELD,NPFT) ! WORK Mean leaf turnover rate over
!                                   !      phenology period (/360days).
     &,Z0(LAND_FIELD)               ! WORK Roughness length on
!                                   !      land points (m).

!-----------------------------------------------------------------------
! Initialisations
!-----------------------------------------------------------------------
      DO L=1,LAND_FIELD
        ALBSNC(L)=0.0
        ALBSNF(L)=0.0
        Z0(L)=0.0
      ENDDO

      DO N=1,NTYPE
        DO L=1,LAND_FIELD
          Z0_T(L,N)=0.0
        ENDDO
      ENDDO

      DO N=1,NTYPE-1
        DO L=1,LAND_FIELD
          CATCH_T(L,N)=0.0
        ENDDO
      ENDDO

      DO N=1,NPFT
        DO L=1,LAND_FIELD
          G_LEAF_PHEN(L,N)=0.0
          G_LEAF_DAY(L,N)=0.0
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Calculate the number of atmospheric timesteps between calls to PHENOL
! and TRIFFID.
!-----------------------------------------------------------------------
      NSTEP_PHEN=INT(86400.0*PHENOL_PERIOD/ATIMESTEP)

!-----------------------------------------------------------------------
! Update halos on input fields
!-----------------------------------------------------------------------
      CALL SWAPB_LAND(LAI,LAND_FIELD,P_FIELD,
     &                ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,
     &                NPFT,LAND_INDEX)

      CALL SWAPB_LAND(HT,LAND_FIELD,P_FIELD,
     &                ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,
     &                NPFT,LAND_INDEX)

      CALL SWAPB_LAND(G_LEAF_AC,LAND_FIELD,P_FIELD,
     &                ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,
     &                NPFT,LAND_INDEX)

!-----------------------------------------------------------------------
! Create the TILE_INDEX array of land points with each surface type
!-----------------------------------------------------------------------
      CALL TILEPTS(P_FIELD,LAND_FIELD,LAND1,LAND_PTS,
     &             FRAC,TILE_PTS,TILE_INDEX)

      IF (L_PHENOL .AND. MOD(A_STEP,NSTEP_PHEN).EQ.0) THEN

!-----------------------------------------------------------------------
! Calculate the phenology timestep in years.
!-----------------------------------------------------------------------
        DTIME_PHEN=FLOAT(PHENOL_PERIOD)/360.0


        DO N=1,NPFT

!-----------------------------------------------------------------------
! Calculate the mean turnover rate and update the leaf phenological
! state.
!-----------------------------------------------------------------------
          DO J=1,TILE_PTS(N)
            L=TILE_INDEX(J,N)
            G_LEAF_DAY(L,N)=G_LEAF_AC(L,N)/DTIME_PHEN
          ENDDO

          WRITE(6,*) 'Calling phenology'

          CALL PHENOL (LAND_FIELD,TILE_PTS(N),TILE_INDEX(1,N),N,
     &                 G_LEAF_DAY(1,N),HT(1,N),DTIME_PHEN,
     &                 G_LEAF_PHEN(1,N),LAI(1,N))

          WRITE(6,*) 'Phenology completed normally'

          DO L=1,LAND_FIELD
            LAI_PHEN(L,N)=LAI(L,N)
          ENDDO

!-----------------------------------------------------------------------
! Reset the accumulation over atmospheric model timesteps to zero.
!-----------------------------------------------------------------------
          DO L=1,LAND_FIELD
            G_LEAF_AC(L,N)=0.0
          ENDDO
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
! Calculate gridbox mean vegetation parameters from fractions of
! surface functional types
!-----------------------------------------------------------------------
      CALL SPARM (LAND_FIELD,LAND1,LAND_PTS,TILE_PTS,TILE_INDEX
     &,           ALB_SOIL,FRAC,HT,LAI
     &,           ALBSNC,ALBSNF,CATCH_T,Z0,Z0_T)

!-----------------------------------------------------------------------
! Copy Z0 from land field to full field
!-----------------------------------------------------------------------
      DO L=LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)
        Z0_P(I)=Z0(L)
      ENDDO

      RETURN
      END
