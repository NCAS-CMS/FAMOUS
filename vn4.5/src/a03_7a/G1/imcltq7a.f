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
!!!  SUBROUTINE IM_CAL_TQ ----------------------------------------------
!!!
!!!  Purpose: Calculate increments for T and Q in the boundary layer,
!!!           using an implicit numerical scheme. The tridiagonal
!!!           matrices are inverted using simple Gaussian elimination.
!!!
!!!  Suitable for single column use; activate *IF definition IBM.
!!!
!!!  Model           Modification history from model version 4.1
!!! version  Date
!!!
!!!  SJ, RE  <- Programmers of some or all of previous code or changes
!!!
!!!  Version for MOSES II land surface tile model (RE 21/3/97)
!!!
!!!  Programming standard: UM Documentation Paper No 4, Version 2,
!!!                        dated 18/1/90
!!!
!!!  System component covered: P244
!!!
!!!  Project task: P24
!!!
!!!  Documentation: UM Documentation Paper No 24.
!!!
!!!---------------------------------------------------------------------

      SUBROUTINE IM_CAL_TQ (P_FIELD,P1,P_POINTS,BL_LEVELS,LAND_FIELD,
     &                      LAND_INDEX,NTYPE,TILE_INDEX,TILE_PTS,
     &                      LAND_MASK,LTIMER,
     &                      ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_SNOW,
     &                      DTL_NT,DQW_NT,DTRDZ,ICE_FRACT,
     &                      RDZ,RESFT,RHOKH,
     &                      RHOKH_1,RHOKH1_SICE,RHOKPM,RHOKPM_SICE,
     &                      TILE_FRAC,
     &                      FQW,FQW_ICE,FQW_TILE,E_SEA,
     &                      FTL,FTL_ICE,FTL_TILE,H_SEA,QW,TL)

      IMPLICIT NONE

      INTEGER
     & P_FIELD                     ! IN No. of points in P-grid.
     &,P1                          ! IN First point to be processed in
!                                  !    P-grid.
     &,P_POINTS                    ! IN Number of P-grid points to be
!                                  !    processed.
     &,BL_LEVELS                   ! IN No. of atmospheric levels for
!                                  !    which boundary layer fluxes are
!                                  !    calculated.
     &,LAND_FIELD                  ! IN Total number of land points.
     &,LAND_INDEX(P_FIELD)         ! IN Index of land points.
     &,NTYPE                       ! IN Number of land surface tiles.
     &,TILE_PTS(NTYPE)             ! IN Number of tiles.
     &,TILE_INDEX(LAND_FIELD,NTYPE)! IN Index of tile points.

      REAL
     & ALPHA1(LAND_FIELD,NTYPE)    ! IN Gradient of saturated specific
!                                  !    humidity with respect to
!                                  !    temperature between the bottom
!                                  !    model layer and tile surfaces.
     &,ALPHA1_SICE(P_FIELD)        ! IN ALPHA1 for sea-ice
     &,ASHTF(P_FIELD)              ! IN Coefficient to calculate surface
!                                  !    heat flux into soil or sea-ice
!                                  !    (W/m2/K).
     &,ASHTF_SNOW(P_FIELD)         ! IN Coefficient to calculate surface
!                                  !    heat flux into snow (W/m2/K).
     &,DTL_NT(P_FIELD,BL_LEVELS)   ! IN Non-turbulent increment for TL.
     &,DQW_NT(P_FIELD,BL_LEVELS)   ! IN Non-turbulent increment for QW.
     &,DTRDZ(P_FIELD,BL_LEVELS)    ! IN -g.dt/dp for model layers.
     &,ICE_FRACT(P_FIELD)          ! IN Fraction of grid-box which is
!                                  !    sea-ice (decimal fraction).
     &,RDZ(P_FIELD,BL_LEVELS)      ! IN 1./dz for model layers.
     &,RESFT(LAND_FIELD,NTYPE)     ! IN Total resistance factor
     &,RHOKH(P_FIELD,2:BL_LEVELS)  ! IN Exchange coeff for FTL above
!                                  !    surface.
     &,RHOKH_1(LAND_FIELD,NTYPE)   ! IN Land surface exchange coeff.
     &,RHOKH1_SICE(P_FIELD)        ! IN Sea and sea-ice surface exchange
!                                  !    coeff.
     &,RHOKPM(LAND_FIELD,NTYPE)    ! IN Land surface exchange coeff.
     &,RHOKPM_SICE(P_FIELD)        ! IN Sea-ice surface exchange coeff.
     &,TILE_FRAC(LAND_FIELD,NTYPE)
!                                  ! IN Tile fractions.

       LOGICAL
     & LAND_MASK(P_FIELD)          ! IN T for land, F elsewhere.
     &,LTIMER                      ! IN Logical for TIMER.

      REAL
     & FQW(P_FIELD,BL_LEVELS)      ! INOUT Flux of QW (ie., for surface,
!                                  !       total evaporation). Kg/sq m/s
     &,FQW_ICE(P_FIELD)            ! INOUT Surface flux of QW for
!                                  !       sea-ice.
     &,FQW_TILE(LAND_FIELD,NTYPE)  ! INOUT Surface flux of QW for land
!                                  !       tiles.
     &,FTL(P_FIELD,BL_LEVELS)      ! INOUT Flux of TL (ie., for surface,
!                                  !       H/Cp where H is sensible heat
!                                  !       in W per sq m).
     &,FTL_ICE(P_FIELD)            ! INOUT H/Cp for sea-ice.
     &,FTL_TILE(LAND_FIELD,NTYPE)  ! INOUT H/Cp for land tiles.
     &,E_SEA(P_FIELD)              ! INOUT Evaporation from sea times
!                                  !       leads fraction (kg/m2/s).
!                                  !       Zero over land.
     &,H_SEA(P_FIELD)              ! INOUT Surface sensible heat flux
!                                  !       over sea times leads fraction
!                                  !       Zero over land.
     &,QW(P_FIELD,BL_LEVELS)       ! INOUT Total water content (kg per
!                                  !       kg air).  From P243.
     &,TL(P_FIELD,BL_LEVELS)       ! INOUT Liquid/frozen water
!                                          temperature (K).  From P243.

!  External references :-
      EXTERNAL TIMER

!  Local and other symbolic constants :-
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
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

C*L------------------COMDECK C_GAMMA------------------------------------
C GAMMA holds the implicit weighting coeff. for up to 30 B.L. levels.
C It is only required for the the number of B.L. levels actually used,
C so it does not need to be set up to 30 when less BL levels are used.
      REAL GAMMA(30)       ! Max of 30 Boundary Layer levels assumed.
C
      DATA GAMMA / 2 * 2.0 , 1.5 , 27 * 1.0 /
C*----------------------------------------------------------------------
      REAL LS
      PARAMETER (
     & LS=LC+LF     ! Latent heat of sublimation (J per kg).
     &)
! Workspace :-
      REAL
     & AQ(P_FIELD,BL_LEVELS)   ! Matrix elements.
     &,AT(P_FIELD,BL_LEVELS)   ! Matrix elements.
     &,BQ1(P_FIELD)            ! Matrix element.
     &,BT1(P_FIELD)            ! Matrix element.
     &,CQ1(P_FIELD)            ! Matrix element.
     &,CT1(P_FIELD)            ! Matrix element.
     &,DQW(P_FIELD,BL_LEVELS)  ! Delta QW.
     &,DTL(P_FIELD,BL_LEVELS)  ! Delta TL.

!  Local scalars :-

      REAL
     & CQ       ! Matrix element.
     &,CT       ! Matrix element.
     &,RBQ      ! Reciprocal of BQ'.
     &,RBT      ! Reciprocal of BT'.

      INTEGER
     & BLM1     ! BL_LEVELS minus 1.
     &,I        ! Loop counter (horizontal field index).
     &,J        ! Loop counter (tile field index).
     &,K        ! Loop counter (vertical index).
     &,L        ! Loop counter (horizontal field index).
     &,N        ! Loop counter (tile index).

      IF (LTIMER) THEN
        CALL TIMER('IMCALTQ ',3)
      ENDIF

      BLM1 = BL_LEVELS-1

!-----------------------------------------------------------------------
!! 1.  Calculate those matrix and vector elements on the LHS of eqn
!!     P244.79 which are to do with implicit solution of the moisture
!!     transport problem - upward sweep through lower half of matrix.
!-----------------------------------------------------------------------
!! 1.1 Row of matrix for QW on top boundary layer level.
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1
        DQW(I,BL_LEVELS) = DTRDZ(I,BL_LEVELS) * FQW(I,BL_LEVELS)
     &                     +DQW_NT(I,BL_LEVELS)               ! P244.58
        AQ(I,BL_LEVELS) = - DTRDZ(I,BL_LEVELS)                ! P244.56
     &       * GAMMA(BL_LEVELS)*RHOKH(I,BL_LEVELS)*RDZ(I,BL_LEVELS)
        RBQ = 1.0 / ( 1.0 - AQ(I,BL_LEVELS) )   ! Reciprocal of P244.98
        DQW(I,BL_LEVELS) = RBQ * DQW(I,BL_LEVELS)             ! P244.99
        AQ(I,BL_LEVELS) = RBQ * AQ(I,BL_LEVELS)               ! P244.100
      ENDDO

!-----------------------------------------------------------------------
!! 1.2 Rows of matrix for QW on intermediate levels.
!-----------------------------------------------------------------------
      DO K=BLM1,2,-1
        DO I=P1,P1+P_POINTS-1
          DQW(I,K) = - DTRDZ(I,K) * ( FQW(I,K+1) - FQW(I,K) )
     &               + DQW_NT(I,K)                            ! P244.54
          CQ = - DTRDZ(I,K)*GAMMA(K+1)*RHOKH(I,K+1)*RDZ(I,K+1)! P244.52
          AQ(I,K) = - DTRDZ(I,K)*GAMMA(K)*RHOKH(I,K)*RDZ(I,K) ! P244.51
          RBQ = 1.0 / ( 1.0 - AQ(I,K) - CQ*(1.0 + AQ(I,K+1)) )
!                                             ... reciprocal of P244.101
          DQW(I,K) = RBQ * ( DQW(I,K) - CQ*DQW(I,K+1) )       ! P244.102
          AQ(I,K) = RBQ * AQ(I,K)                             ! P244.103
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!! 1.3 Row of matrix for QW on lowest level.
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1
        DQW(I,1) = - DTRDZ(I,1)*(FQW(I,2) - FQW(I,1)) + DQW_NT(I,1)
        AQ(I,1) = 0.
        CQ1(I) = - GAMMA(2)*DTRDZ(I,1)*RHOKH(I,2)*RDZ(I,2)
        BQ1(I) = 1. - (1. + AQ(I,2))*CQ1(I)
        IF ( .NOT. LAND_MASK(I) ) THEN
! Sea or sea-ice
          AQ(I,1) = - CP*GAMMA(1)*DTRDZ(I,1)*ICE_FRACT(I) *
     &                      ALPHA1_SICE(I)*RHOKH1_SICE(I)*RHOKPM_SICE(I)
          BQ1(I) = BQ1(I) + GAMMA(1)*DTRDZ(I,1) *
     &                      ( (1. - ICE_FRACT(I))*RHOKH1_SICE(I) +
     &        ICE_FRACT(I)*(ASHTF(I)+CP*RHOKH1_SICE(I))*RHOKPM_SICE(I) )
        ENDIF
      ENDDO

! Snow-free land tiles
      DO N=1,NTYPE-1
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          I = LAND_INDEX(L)
          AQ(I,1) = AQ(I,1) - CP*GAMMA(1)*DTRDZ(I,1)*TILE_FRAC(L,N) *
     &                   RESFT(L,N)*ALPHA1(L,N)*RHOKH_1(L,N)*RHOKPM(L,N)
          BQ1(I) = BQ1(I) + GAMMA(1)*DTRDZ(I,1)*TILE_FRAC(L,N) *
     &               RESFT(L,N)*RHOKPM(L,N)*(CP*RHOKH_1(L,N) + ASHTF(I))
        ENDDO
      ENDDO

! Snow and land-ice tile
      N = NTYPE
      DO J=1,TILE_PTS(N)
        L = TILE_INDEX(J,N)
        I = LAND_INDEX(L)
        AQ(I,1) = AQ(I,1) - CP*GAMMA(1)*DTRDZ(I,1)*TILE_FRAC(L,N) *
     &                              ALPHA1(L,N)*RHOKH_1(L,N)*RHOKPM(L,N)
        BQ1(I) = BQ1(I) + GAMMA(1)*DTRDZ(I,1)*TILE_FRAC(L,N) *
     &                     RHOKPM(L,N)*(CP*RHOKH_1(L,N) + ASHTF_SNOW(I))
      ENDDO

      DO I=P1,P1+P_POINTS-1
        DQW(I,1) = (DQW(I,1)  - CQ1(I)*DQW(I,2)) / BQ1(I)
        AQ(I,1) = AQ(I,1) / BQ1(I)
      ENDDO

!-----------------------------------------------------------------------
!! 2.  Calculate those matrix and vector elements on the LHS of eqn
!!     P244.79 which are to do with implicit solution of the heat
!!     transport problem - upward sweep through upper half of matrix.
!-----------------------------------------------------------------------
!! 2.1 Row of matrix for TL on lowest level.
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1
        DTL(I,1) = - DTRDZ(I,1)*(FTL(I,2) - FTL(I,1)) + DTL_NT(I,1)
        AT(I,1) = - DTRDZ(I,1)*GAMMA(2)*RHOKH(I,2)*RDZ(I,2)
        BT1(I) = 1. - AT(I,1)
        CT1(I) = 0.
        IF ( .NOT.LAND_MASK(I) ) THEN
! Sea or sea-ice
          CT1(I) = - LS*GAMMA(1)*DTRDZ(I,1)*ICE_FRACT(I) *
     &                                     RHOKH1_SICE(I)*RHOKPM_SICE(I)
          BT1(I) = BT1(I) - ALPHA1_SICE(I)*CT1(I) + GAMMA(1)*DTRDZ(I,1)
     &                        * ( ICE_FRACT(I)*ASHTF(I)*RHOKPM_SICE(I)
     &                            + (1. - ICE_FRACT(I))*RHOKH1_SICE(I) )
        ENDIF
      ENDDO

! Snow-free land tiles
      DO N=1,NTYPE-1
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          I = LAND_INDEX(L)
          CT1(I) = CT1(I) - LC*GAMMA(1)*DTRDZ(I,1)*TILE_FRAC(L,N) *
     &                               RESFT(L,N)*RHOKH_1(L,N)*RHOKPM(L,N)
          BT1(I) = BT1(I) + GAMMA(1)*DTRDZ(I,1)*TILE_FRAC(L,N) *
     &             RHOKPM(L,N)*( LC*ALPHA1(L,N)*RESFT(L,N)*RHOKH_1(L,N)
     &                                                      + ASHTF(I) )
        ENDDO
      ENDDO

! Snow and land-ice tile
      N = NTYPE
      DO J=1,TILE_PTS(N)
        L = TILE_INDEX(J,N)
        I = LAND_INDEX(L)
        CT1(I) = CT1(I) - LS*GAMMA(1)*DTRDZ(I,1)*TILE_FRAC(L,N) *
     &                                          RHOKH_1(L,N)*RHOKPM(L,N)
        BT1(I) = BT1(I) + GAMMA(1)*DTRDZ(I,1)*TILE_FRAC(L,N) *
     &                        RHOKPM(L,N)*( LS*ALPHA1(L,N)*RHOKH_1(L,N)
     &                                                 + ASHTF_SNOW(I) )
      ENDDO

      DO I=P1,P1+P_POINTS-1
        BT1(I) = BT1(I) - CT1(I)*AQ(I,1)
        DTL(I,1) = (DTL(I,1) - CT1(I)*DQW(I,1)) / BT1(I)
        AT(I,1) = AT(I,1) / BT1(I)
      ENDDO

!-----------------------------------------------------------------------
!! 2.2 Rows of matrix for TL on intermediate levels.
!-----------------------------------------------------------------------
      DO K=2,BLM1
        DO I=P1,P1+P_POINTS-1
          DTL(I,K) = - DTRDZ(I,K) * ( FTL(I,K+1) - FTL(I,K) )
     &               + DTL_NT(I,K)                            ! P244.38
          AT(I,K) = - DTRDZ(I,K)*GAMMA(K+1)*RHOKH(I,K+1)*RDZ(I,K+1)
!                                                             ! P244.35
          CT = - DTRDZ(I,K)*GAMMA(K)*RHOKH(I,K)*RDZ(I,K)      ! P244.36
          RBT = 1.0 / ( 1.0 - AT(I,K) - CT*(1.0 + AT(I,K-1)) )
!                                             ... Reciprocal of P244.113
          DTL(I,K) = RBT * ( DTL(I,K) - CT*DTL(I,K-1) )       ! P244.114
          AT(I,K) = RBT * AT(I,K)                             ! P244.115
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!! 2.3 Row of matrix for TL on top boundary layer level.
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1
        DTL(I,BL_LEVELS) = DTRDZ(I,BL_LEVELS) * FTL(I,BL_LEVELS)
     &                   + DTL_NT(I,BL_LEVELS)                ! P244.42
        CT = - DTRDZ(I,BL_LEVELS)*GAMMA(BL_LEVELS)*
     &         RHOKH(I,BL_LEVELS)*RDZ(I,BL_LEVELS)
        RBT = 1.0 / ( 1.0 - CT*(1.0 + AT(I,BLM1)) )
!                                             ... Reciprocal of P244.116
        DTL(I,BL_LEVELS) = RBT * ( DTL(I,BL_LEVELS) - CT*DTL(I,BLM1) )
!                                                             ! P244.117
      ENDDO

!-----------------------------------------------------------------------
!! 3.  Downward sweep through whole matrix.
!-----------------------------------------------------------------------
!! 3.1 Increment TL.
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1
        TL(I,BL_LEVELS) = TL(I,BL_LEVELS) + DTL(I,BL_LEVELS)  ! P244.127
      ENDDO
      DO K=BLM1,1,-1
        DO I=P1,P1+P_POINTS-1
          DTL(I,K) = DTL(I,K) - AT(I,K)*DTL(I,K+1)            ! P244.118
          TL(I,K) = TL(I,K) + DTL(I,K)                        ! P244.127
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!! 3.2 Increment QW.
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1
        DQW(I,1) = DQW(I,1) - AQ(I,1)*DTL(I,1)
        QW(I,1) = QW(I,1) + DQW(I,1)
      ENDDO
      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          DQW(I,K) = DQW(I,K) - AQ(I,K)*DQW(I,K-1)            ! P244.121
          QW(I,K) = QW(I,K) + DQW(I,K)                        ! P244.128
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!! 4.  Calculate final implicit fluxes of heat and moisture.
!-----------------------------------------------------------------------
!! 4.1 Surface fluxes.
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1
        FTL(I,1) = 0.
        FQW(I,1) = 0.
        IF ( .NOT.LAND_MASK(I) ) THEN
          IF ( ICE_FRACT(I) .GT. 0. ) THEN
! Sea-ice and leads
            FTL_ICE(I) = FTL_ICE(I) +
     &                           GAMMA(1)*ICE_FRACT(I)*RHOKPM_SICE(I) *
     &        ( LS*RHOKH1_SICE(I)*(DQW(I,1) - ALPHA1_SICE(I)*DTL(I,1))
     &                                             - ASHTF(I)*DTL(I,1) )
            FQW_ICE(I) = FQW_ICE(I) -
     &                           GAMMA(1)*ICE_FRACT(I)*RHOKPM_SICE(I) *
     &        ( CP*RHOKH1_SICE(I)*(DQW(I,1) - ALPHA1_SICE(I)*DTL(I,1))
     &                                             + ASHTF(I)*DQW(I,1) )
            H_SEA(I) = H_SEA(I) - (1.0 - ICE_FRACT(I)) * CP *
     &                                  GAMMA(1)*RHOKH1_SICE(I)*DTL(I,1)
            E_SEA(I) = E_SEA(I) - (1.0 - ICE_FRACT(I)) *
     &                                  GAMMA(1)*RHOKH1_SICE(I)*DQW(I,1)
            FTL(I,1) = FTL_ICE(I) + H_SEA(I)/CP
            FQW(I,1) = FQW_ICE(I) + E_SEA(I)
          ELSE
! Unfrozen sea
            H_SEA(I) = H_SEA(I) - CP*GAMMA(1)*RHOKH1_SICE(I)*DTL(I,1)
            E_SEA(I) = E_SEA(I) - GAMMA(1)*RHOKH1_SICE(I)*DQW(I,1)
            FTL(I,1) = H_SEA(I)/CP
            FQW(I,1) = E_SEA(I)
          ENDIF
        ENDIF
      ENDDO

! Snow-free land tiles
      DO N=1,NTYPE-1
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          I = LAND_INDEX(L)
          FTL_TILE(L,N) = FTL_TILE(L,N) + GAMMA(1)*RHOKPM(L,N) *
     &     ( LC*RESFT(L,N)*RHOKH_1(L,N)*(DQW(I,1)-ALPHA1(L,N)*DTL(I,1))
     &                                             - ASHTF(I)*DTL(I,1) )
          FQW_TILE(L,N) = FQW_TILE(L,N) - GAMMA(1)*RESFT(L,N) *
     &    RHOKPM(L,N)*( CP*RHOKH_1(L,N)*(DQW(I,1)-ALPHA1(L,N)*DTL(I,1))
     &                                             + ASHTF(I)*DQW(I,1) )
          FTL(I,1) = FTL(I,1) + TILE_FRAC(L,N)*FTL_TILE(L,N)
          FQW(I,1) = FQW(I,1) + TILE_FRAC(L,N)*FQW_TILE(L,N)
        ENDDO
      ENDDO

! Snow and land-ice tile
      N = NTYPE
      DO J=1,TILE_PTS(N)
        L = TILE_INDEX(J,N)
        I = LAND_INDEX(L)
        FTL_TILE(L,N) = FTL_TILE(L,N) + GAMMA(1)*RHOKPM(L,N) *
     &              ( LS*RHOKH_1(L,N)*(DQW(I,1) - ALPHA1(L,N)*DTL(I,1))
     &                                        - ASHTF_SNOW(I)*DTL(I,1) )
        FQW_TILE(L,N) = FQW_TILE(L,N) - GAMMA(1)*RHOKPM(L,N) *
     &              ( CP*RHOKH_1(L,N)*(DQW(I,1) - ALPHA1(L,N)*DTL(I,1))
     &                                        + ASHTF_SNOW(I)*DQW(I,1) )
        FTL(I,1) = FTL(I,1) + TILE_FRAC(L,N)*FTL_TILE(L,N)
        FQW(I,1) = FQW(I,1) + TILE_FRAC(L,N)*FQW_TILE(L,N)
      ENDDO

!-----------------------------------------------------------------------
!! 4.2 Fluxes at layer interfaces above the surface.
!-----------------------------------------------------------------------
      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          FTL(I,K)  = FTL(I,K) - GAMMA(K)*RHOKH(I,K)*RDZ(I,K) ! P244.33
     &                              * ( DTL(I,K) - DTL(I,K-1) )
          FQW(I,K)  = FQW(I,K) - GAMMA(K)*RHOKH(I,K)*RDZ(I,K) ! P244.44
     &                              * ( DQW(I,K) - DQW(I,K-1) )
        ENDDO
      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('IMCALTQ ',4)
      ENDIF

      RETURN
      END
