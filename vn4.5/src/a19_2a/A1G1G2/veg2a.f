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
      SUBROUTINE VEG(P_FIELD,FIRST_POINT,LAST_POINT,LAND_FIELD
     &,              LAND1,LAND_PTS,LAND_INDEX,P_ROWS,ROW_LENGTH
     &,              EW_Halo,NS_Halo
     &,              A_STEP,ASTEPS_SINCE_TRIFFID
     &,              PHENOL_PERIOD,TRIFFID_PERIOD
     &,              L_PHENOL,L_TRIFFID,L_TRIF_EQ
     &,              ALB_SOIL,ATIMESTEP,FRAC_DISTURB
     &,              G_LEAF_AC,G_LEAF_PHEN_AC,NPP_AC
     &,              RESP_S_AC,RESP_W_AC
     &,              CS,FRAC,LAI,HT
     &,              ALBSNC,ALBSNF,CATCH_T,Z0_P,Z0_T
     &,              C_VEG,CV,LIT_C,LIT_C_MN,G_LEAF_DAY,G_LEAF_PHEN
     &,              LAI_PHEN,G_LEAF_DR_OUT,NPP_DR_OUT,RESP_W_DR_OUT
     &,              RESP_S_DR_OUT
     &               )


      IMPLICIT NONE
!
! Description:
!   Updates Leaf Area Index for Plant Functional Types (PFTs) and uses
!   this to derive new vegetation parameters for PFTs along with gridbox
!   mean values where appropriate.
!
! Method:
!   Calls PHENOL which models phenolgy and updates Leaf Area Index
!   (LAI), then calls TRIFFID to update vegetation and soil fractions, 
!   LAI, canopy height, veg and soil carbon and carbon fluxes.  Passes
!   fractions, LAI and canopy height to SPARM which derives the 
!   vegetation parameters for each PFT and also the gridbox means where 
!   this is required.
!
! Current Code Owner:  Richard Betts
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.4    8/10/97   Original code.  Richard Betts
!   4.5   12/05/98   Find total fraction of gridbox covered by 
!                    vegetation or soil, use this to derive indices of 
!                    land points on which TRIFFID may operate, and pass 
!                    both to TRIFFID.  Initialise top and bottom rows 
!                    for all variables.  Richard Betts
!   4.5   30/06/98   Add second call to TILEPTS to update TILE_INDEX
!                    after TRIFFID.  Richard Betts
!   4.5    6/08/98   Call SWAPB_LAND to update halo regions of input
!                    fields.   Richard Betts
!   4.5   23/11/98   Output G_LEAF_DAY, G_LEAF_PHEN, LAI_PHEN, 
!                    G_LEAF_DR_OUT, NPP_DR_OUT, RESP_W_DR_OUT and 
!                    RESP_S_DR_OUT as diagnostics.  Richard Betts
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
     &,LAND_PTS              ! IN Number of land points.
     &,P_ROWS                ! IN Number of rows on P grid.
     &,ROW_LENGTH            ! IN Number of P points in a row.
     &,EW_Halo               ! IN Halo size in the EW direction.
     &,NS_Halo               ! IN Halo size in the NS direction.
     &,A_STEP                ! IN Atmospheric timestep number.
     &,ASTEPS_SINCE_TRIFFID  ! INOUT Number of atmosphere
C                                    timesteps since last call
C                                           to TRIFFID.
     &,PHENOL_PERIOD         ! IN Phenology period (days).
     &,TRIFFID_PERIOD        ! IN TRIFFID period (days).

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
C                                   !    P-point is the Lth land
C                                   !    point.

      INTEGER
     & I,J,K,L,N                    ! WORK loop counters.
     &,KITER                        ! WORK Number of TRIFFID iterations.

      LOGICAL
     & L_PHENOL                     ! IN .T. for interactive leaf
C                                   !    phenology.
     &,L_TRIFFID                    ! IN .T. for interactive vegetation.
     &,L_TRIF_EQ                    ! IN .T. for vegetation equilibrium.

      REAL
     & ALB_SOIL(LAND_FIELD)         ! IN snow-free albedo of soil.
     &,ATIMESTEP                    ! IN Atmospheric timestep (s).
     &,FRAC_DISTURB(LAND_FIELD)     ! IN Fraction of gridbox in which
C                                   !    vegetation is disturbed.
     &,G_LEAF_AC(LAND_FIELD,NPFT)   ! INOUT Accumulated leaf turnover
C                                   !       rate.
     &,G_LEAF_PHEN_AC(LAND_FIELD,NPFT)! INOUT Accumulated leaf turnover
C                                   !       rate including phenology.
     &,NPP_AC(LAND_FIELD,NPFT)      ! INOUT Accumulated NPP (kg C/m2).
     &,RESP_W_AC(LAND_FIELD,NPFT)   ! INOUT Accumulated wood respiration
C                                   !       (kg C/m2).
     &,RESP_S_AC(LAND_FIELD)        ! INOUT Accumulated soil respiration
C                                   !       (kg C/m2).
     &,CS(LAND_FIELD)               ! INOUT Soil carbon content
C                                   !       (kg C/m2).
     &,FRAC(LAND_FIELD,NTYPE)       ! INOUT Fractions of surface types.
     &,LAI(LAND_FIELD,NPFT)         ! INOUT LAI of plant functional
C                                   !       types.
     &,HT(LAND_FIELD,NPFT)          ! INOUT Height of plant functional
C                                   !       types (m).
     &,ALBSNC(LAND_FIELD)           ! OUT Snow-covered albedo.
     &,ALBSNF(LAND_FIELD)           ! OUT Snow-free albedo.
     &,CATCH_T(LAND_FIELD,NTYPE-1)  ! OUT Canopy capacity for each type
C                                   !     aside from ice (kg/m2).
     &,Z0_P(P_FIELD)                ! OUT Effective roughness length
C                                   !     on full grid (m).
     &,Z0_T(LAND_FIELD,NTYPE)       ! OUT Roughness length for each type
C                                   !     (m).
     &,C_VEG(LAND_FIELD,NPFT)       ! OUT Total carbon content of
C                                   !     the vegetation (kg C/m2).
     &,CV(LAND_FIELD)               ! OUT Gridbox mean vegetation
C                                   !     carbon (kg C/m2).
     &,G_LEAF_DAY(LAND_FIELD,NPFT)  ! OUT Mean leaf turnover rate for
!                                   !      input to PHENOL (/360days).
     &,G_LEAF_DR_OUT(LAND_FIELD,NPFT) ! OUT Mean leaf turnover rate for 
!                                   !       driving TRIFFID (/360days).
     &,LAI_PHEN(LAND_FIELD,NPFT)    ! OUT LAI of PFTs after phenology.
     &,LIT_C(LAND_FIELD,NPFT)       ! OUT Carbon Litter 
!                                   !     (kg C/m2/360days).
     &,LIT_C_MN(LAND_FIELD)         ! OUT Gridbox mean carbon litter
!                                   !     (kg C/m2/360days).
     &,NPP_DR_OUT(LAND_FIELD,NPFT)  ! OUT Mean NPP for driving TRIFFID 
!                                   !     (kg C/m2/360days).
     &,RESP_W_DR_OUT(LAND_FIELD,NPFT) ! OUT Mean wood respiration for
!                                   !       driving TRIFFID 
!                                   !       (kg C/m2/360days).
     &,RESP_S_DR_OUT(LAND_FIELD)    ! OUT Mean soil respiration for
!                                   !     driving TRIFFID 
!                                   !     (kg C/m2/360days).

      INTEGER
     & NSTEP_PHEN                   ! WORK Number of atmospheric
C                                   !      timesteps between calls to
C                                   !      PHENOL.
     &,NSTEP_TRIF                   ! WORK Number of atmospheric
C                                   !      timesteps between calls to
C                                   !      TRIFFID.
     &,TILE_PTS(NTYPE)              ! WORK Number of land points which
C                                   !      include the nth surface type.
     &,TILE_INDEX(LAND_FIELD,NTYPE) ! WORK Indices of land points which
C                                   !      include the nth surface type.
     &,TRIF_PTS                     ! WORK Number of points on which
!                                   !      TRIFFID may operate
     &,TRIF_INDEX(LAND_FIELD)       ! WORK Indices of land points on 
!                                   !      which TRIFFID may operate

      REAL
     & DTIME_PHEN                   ! WORK The phenology timestep (yr).
     &,FORW                         ! WORK Forward timestep weighting
C                                   !      for TRIFFID.
     &,GAMMA                        ! WORK Inverse TRIFFID timestep
!                                   !      (/360days).
     &,GAM_TRIF                     ! WORK Inverse TRIFFID coupling
!                                   !      timestep (/360days).
     &,G_ANTH(LAND_FIELD)           ! WORK Anthropogenic disturbance
!                                   !      rate (/360days).
     &,G_LEAF_PHEN(LAND_FIELD,NPFT) ! WORK Mean leaf turnover rate over
!                                   !      phenology period (/360days).
     &,G_LEAF_DR(LAND_FIELD,NPFT)   ! WORK Mean leaf turnover rate
!                                   !      for driving TRIFFID 
!                                   !      (/360days).
     &,NPP_DR(LAND_FIELD,NPFT)      ! WORK Mean NPP for driving
!                                   !      TRIFFID (kg C/m2/360days).
     &,RESP_W_DR(LAND_FIELD,NPFT)   ! WORK Mean wood respiration for
!                                   !      driving TRIFFID 
!                                   !      (kg C/m2/360days).
     &,RESP_S_DR(LAND_FIELD)        ! WORK Mean soil respiration for
!                                   !      driving TRIFFID 
!                                   !      (kg C/m2/360days).
     &,FRAC_VS(LAND_FIELD)          ! WORK Total fraction of gridbox 
!                                   !      covered by veg or soil.
     &,Z0(LAND_FIELD)               ! WORK Roughness length on
C                                   !      land points (m).
C-----------------------------------------------------------------------
C Local parameters
C-----------------------------------------------------------------------
      INTEGER
     + ITER_EQ                    ! Number of TRIFFID iterations for
C                                 ! gradient descent to equilibrium.
      REAL
     + GAMMA_EQ                   ! Inverse timestep for gradient
C                                 ! descent to equilibrium (/360days).
      PARAMETER(GAMMA_EQ = 1.0E-4, ITER_EQ = 10)

      REAL
     + FRAC_MIN                   ! Minimum ("seed") areal fraction.
      PARAMETER(FRAC_MIN = 0.01)
      REAL
     & G_ANTH0                      ! Anthropogenic disturbance rate
!                                   ! (/360days).
      PARAMETER (G_ANTH0=0.0)

C-----------------------------------------------------------------------
C Initialisations
C-----------------------------------------------------------------------
      DO N=1,NPFT
        DO L=1,LAND_FIELD
          G_LEAF_PHEN(L,N)=0.0
          G_LEAF_DAY(L,N)=0.0
          G_LEAF_DR(L,N)=0.0
          NPP_DR(L,N)=0.0
          RESP_W_DR(L,N)=0.0
          C_VEG(L,N)=0.0
          LIT_C(L,N)=0.0
        ENDDO
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

      DO L=1,LAND_FIELD
        ALBSNC(L)=0.0
        ALBSNF(L)=0.0
        G_ANTH(L)=0.0
        RESP_S_DR(L)=0.0
        Z0(L)=0.0
        CV(L)=0.0
        LIT_C_MN(L)=0.0
        FRAC_VS(L) = 0.0
      ENDDO

C-----------------------------------------------------------------------
C Calculate the number of atmospheric timesteps between calls to PHENOL
C and TRIFFID.
C-----------------------------------------------------------------------
      NSTEP_PHEN=INT(86400.0*PHENOL_PERIOD/ATIMESTEP)
      NSTEP_TRIF=INT(86400.0*TRIFFID_PERIOD/ATIMESTEP)

!-----------------------------------------------------------------------
! Update halos on input fields
!-----------------------------------------------------------------------
      CALL SWAPB_LAND(ALB_SOIL,LAND_FIELD,P_FIELD,
     &                ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,
     &                1,LAND_INDEX)

      CALL SWAPB_LAND(LAI,LAND_FIELD,P_FIELD,
     &                ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,
     &                NPFT,LAND_INDEX)

      CALL SWAPB_LAND(HT,LAND_FIELD,P_FIELD,
     &                ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,
     &                NPFT,LAND_INDEX)

      CALL SWAPB_LAND(G_LEAF_AC,LAND_FIELD,P_FIELD,
     &                ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,
     &                NPFT,LAND_INDEX)


      CALL SWAPB_LAND(FRAC_DISTURB,LAND_FIELD,P_FIELD,
     &                ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,
     &                1,LAND_INDEX)

      CALL SWAPB_LAND(G_LEAF_PHEN_AC,LAND_FIELD,P_FIELD,
     &                ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,
     &                NPFT,LAND_INDEX)

      CALL SWAPB_LAND(NPP_AC,LAND_FIELD,P_FIELD,
     &                ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,
     &                NPFT,LAND_INDEX)

      CALL SWAPB_LAND(RESP_W_AC,LAND_FIELD,P_FIELD,
     &                ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,
     &                NPFT,LAND_INDEX)

      CALL SWAPB_LAND(RESP_S_AC,LAND_FIELD,P_FIELD,
     &                ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,
     &                1,LAND_INDEX)

      CALL SWAPB_LAND(CS,LAND_FIELD,P_FIELD,
     &                ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,
     &                1,LAND_INDEX)

      CALL SWAPB_LAND(FRAC,LAND_FIELD,P_FIELD,
     &                ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,
     &                    NTYPE,LAND_INDEX)


!-----------------------------------------------------------------------
! Find total fraction of gridbox covered by vegetation and soil, and use
! this to set indices of land points on which TRIFFID may operate
!-----------------------------------------------------------------------
      TRIF_PTS = 0
      DO L=LAND1,LAND1+LAND_PTS-1
        DO N=1,NPFT
          FRAC_VS(L) = FRAC_VS(L) + FRAC(L,N)
        ENDDO
        N=SOIL
        FRAC_VS(L) = FRAC_VS(L) + FRAC(L,N)
        IF (FRAC_VS(L).GE.(NPFT*FRAC_MIN)) THEN
          TRIF_PTS = TRIF_PTS + 1
          TRIF_INDEX(TRIF_PTS) = L
        ENDIF
      ENDDO

C-----------------------------------------------------------------------
C Create the TILE_INDEX array of land points with each surface type
C-----------------------------------------------------------------------
      CALL TILEPTS(P_FIELD,LAND_FIELD,LAND1,LAND_PTS,
     &             FRAC,TILE_PTS,TILE_INDEX)

      IF (L_PHENOL .AND. MOD(A_STEP,NSTEP_PHEN).EQ.0) THEN

C-----------------------------------------------------------------------
C Calculate the phenology timestep in years.
C-----------------------------------------------------------------------
        DTIME_PHEN=FLOAT(PHENOL_PERIOD)/360.0

        DO N=1,NPFT

C-----------------------------------------------------------------------
C Calculate the mean turnover rate and update the leaf phenological
! state, and take copy of updated LAI field for output as diagnostic.
C-----------------------------------------------------------------------
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

C-----------------------------------------------------------------------
C Increment the leaf turnover rate for driving TRIFFID and reset the
C accumulation over atmospheric model timesteps to zero.
C-----------------------------------------------------------------------
          DO J=1,TILE_PTS(N)
            L=TILE_INDEX(J,N)
            G_LEAF_PHEN_AC(L,N)=G_LEAF_PHEN_AC(L,N)
     &                    +G_LEAF_PHEN(L,N)*DTIME_PHEN
          ENDDO

          DO L=1,LAND_FIELD
            G_LEAF_AC(L,N)=0.0
          ENDDO

        ENDDO
      ENDIF

C-----------------------------------------------------------------------
C Call TRIFFID vegetation model to update vegetation and terrestrial
C carbon storage.
C-----------------------------------------------------------------------
      IF (L_TRIFFID .AND.
     &   (ASTEPS_SINCE_TRIFFID.EQ.NSTEP_TRIF)) THEN

C-----------------------------------------------------------------------
C Calculate the TRIFFID inverse coupling timestep.
C-----------------------------------------------------------------------
        GAM_TRIF=360.0/FLOAT(TRIFFID_PERIOD)

C-----------------------------------------------------------------------
C Diagnose the mean fluxes over the coupling period.
C-----------------------------------------------------------------------
        DO L=LAND1,LAND1+LAND_PTS-1
          RESP_S_DR(L)=RESP_S_AC(L)*GAM_TRIF
        ENDDO

        DO N=1,NPFT
          DO J=1,TILE_PTS(N)
            L=TILE_INDEX(J,N)
            G_LEAF_DR(L,N)=G_LEAF_PHEN_AC(L,N)*GAM_TRIF
            NPP_DR(L,N)=NPP_AC(L,N)*GAM_TRIF
            RESP_W_DR(L,N)=RESP_W_AC(L,N)*GAM_TRIF
          ENDDO
        ENDDO

C-----------------------------------------------------------------------
C Diagnose the mean leaf turnover rates over the coupling period.
C-----------------------------------------------------------------------
        IF (L_PHENOL) THEN
          DO N=1,NPFT
            DO J=1,TILE_PTS(N)
              L=TILE_INDEX(J,N)
              G_LEAF_DR(L,N)=G_LEAF_PHEN_AC(L,N)*GAM_TRIF
            ENDDO
          ENDDO
        ELSE
          DO N=1,NPFT
            DO J=1,TILE_PTS(N)
              L=TILE_INDEX(J,N)
              G_LEAF_DR(L,N)=G_LEAF_AC(L,N)*GAM_TRIF
            ENDDO
          ENDDO
        ENDIF

C-----------------------------------------------------------------------
C Calculate the anthropogenic disturbance rate
C-----------------------------------------------------------------------
        DO L=LAND1,LAND1+LAND_PTS-1
          G_ANTH(L)=G_ANTH0*FRAC_DISTURB(L)
        ENDDO

!-----------------------------------------------------------------------
! Take copies of TRIFFID input variables for output as diagnostics.
!-----------------------------------------------------------------------
        DO N=1,NPFT
          DO L=1,LAND_FIELD
            G_LEAF_DR_OUT(L,N)=G_LEAF_DR(L,N)
            NPP_DR_OUT(L,N)=NPP_DR(L,N)
            RESP_W_DR_OUT(L,N)=RESP_W_DR(L,N)
          ENDDO
        ENDDO
        DO L=1,LAND_FIELD
          RESP_S_DR_OUT(L)=RESP_S_DR(L)
        ENDDO

C-----------------------------------------------------------------------
C Select timestep and forward timestep weighting parameters for
C equilibrium or dynamic vegetation and call TRIFFID.
C-----------------------------------------------------------------------
        IF (L_TRIF_EQ) THEN
          FORW=1.0
          GAMMA=GAMMA_EQ
          KITER=ITER_EQ
        ELSE
          FORW=0.5
          GAMMA=GAM_TRIF
          KITER=1
        ENDIF

        DO K=1,KITER

          WRITE(6,*) 'Calling TRIFFID'

          CALL TRIFFID (LAND_FIELD,TRIF_PTS,TRIF_INDEX,FORW,GAMMA
     &,                 FRAC_VS,G_ANTH,G_LEAF_DR,NPP_DR,RESP_S_DR
     &,                 RESP_W_DR,CS,FRAC,HT,LAI
     &,               C_VEG,CV,LIT_C,LIT_C_MN)

          WRITE(6,*) 'TRIFFID completed normally'
 
        ENDDO

C-----------------------------------------------------------------------
C Update TILE_INDEX for new surface type fractions.  
C-----------------------------------------------------------------------
        CALL TILEPTS(P_FIELD,LAND_FIELD,LAND1,LAND_PTS,
     &               FRAC,TILE_PTS,TILE_INDEX)

C-----------------------------------------------------------------------
C Reset the accumulations to zero.
C-----------------------------------------------------------------------
        DO L=LAND1,LAND1+LAND_PTS-1
          RESP_S_AC(L)=0.0
        ENDDO

        DO N=1,NPFT
          DO L=LAND1,LAND1+LAND_PTS-1
            NPP_AC(L,N)=0.0
            RESP_W_AC=0.0
          ENDDO
        ENDDO

        IF (L_PHENOL) THEN
          DO N=1,NPFT
            DO L=LAND1,LAND1+LAND_PTS-1
              G_LEAF_PHEN_AC(L,N)=0.0
            ENDDO
          ENDDO
        ELSE
          DO N=1,NPFT
            DO L=LAND1,LAND1+LAND_PTS-1
              G_LEAF_AC(L,N)=0.0
            ENDDO
          ENDDO
        ENDIF

        ASTEPS_SINCE_TRIFFID=0

      ENDIF

C-----------------------------------------------------------------------
C Calculate gridbox mean vegetation parameters from fractions of
C surface functional types
C-----------------------------------------------------------------------
      CALL SPARM (LAND_FIELD,LAND1,LAND_PTS,TILE_PTS,TILE_INDEX
     &,           ALB_SOIL,FRAC,HT,LAI
     &,           ALBSNC,ALBSNF,CATCH_T,Z0,Z0_T)

C-----------------------------------------------------------------------
C Copy Z0 from land field to full field
C-----------------------------------------------------------------------
      DO L=LAND1,LAND1+LAND_PTS-1
        I=LAND_INDEX(L)
        Z0_P(I)=Z0(L)
      ENDDO

      RETURN
      END
