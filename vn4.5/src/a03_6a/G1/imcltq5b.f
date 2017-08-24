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
!!!  Purpose: Calculate increments for
!!!           T and Q in the boundary layer, using an
!!!           implicit numerical scheme.  The tridiagonal matrices are
!!!           inverted using simple Gaussian elimination.
!!!
!!!
!!!  Model           Modification history
!!! version  Date
CLL  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!!!
!!!  JJ  <- Programmers of some or all of previous code or changes
!!!
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
!!  Arguments :-
      SUBROUTINE IM_CAL_TQ (
     & P_FIELD,P1
     &,LAND_INDEX
     &,LAND_PTS,LAND1
     &,P_POINTS,BL_LEVELS,N_TYPES,TILE_FRAC
     &,ALPHA1_GB,ALPHA1,ASHTF
     &,DTRDZ,DTRDZ_RML,RHOKH,RDZ
     &,ICE_FRACT,LYING_SNOW,RADNET_C,RESFT,RHOKPM_TILE                  
     &,RHOKPM_POT_TILE
     &,TIMESTEP,LAND_MASK
     &,EPOT,EPOT_TILE
     &,FQW,FQW_TILE,FTL,FTL_TILE,E_SEA,H_SEA,DQW_NT,QW
     &,GAMMA,RHOKE,RHOKH_1,DTL_NT,TL
     &,SURF_HT_FLUX,NRML
     &,LTIMER
     &)

      IMPLICIT NONE

      LOGICAL LTIMER

      INTEGER
     & P_FIELD                     ! IN No. of points in P-grid.
     &,P1                          ! IN First point to be processed in
!                                       P-grid.
     &,LAND1                       ! IN First land point to be processed
     &,P_POINTS                    ! IN Number of P-grid points to be
!                                       processed.
     &,LAND_PTS                    ! IN Number of land points to be
!                                       processed.
     &,BL_LEVELS                   ! IN No. of atmospheric levels for
!                                       which boundary layer fluxes are
!                                       calculated.
     &,N_TYPES                     ! IN Number of land surface tiles
     &,LAND_INDEX(P_FIELD)         ! IN LAND_INDEX(I)=J => the Jth
!                                     point in P_FIELD is the Ith
!                                     land point.

      REAL
     & DTRDZ(P_FIELD,BL_LEVELS)    ! IN dz for bottom BL_LEVELS
     &,RDZ(P_FIELD,BL_LEVELS)      ! IN 1./dz
     &,ALPHA1_GB(P_FIELD)          ! IN Gradient of saturated specific
!                                       humidity with respect to
!                                       temperature between the bottom
!                                       model layer and the surface.
     &,ALPHA1(P_FIELD,N_TYPES)     ! IN Gradient of saturated specific
!                                       humidity with respect to
!                                       temperature between the bottom
!                                       model layer and the surface.
     &,ASHTF(P_FIELD)              ! IN Coefficient to calculate surface
!                                       heat flux into soil or sea-ice
!                                       (W/m2/K).
     &,RHOKH(P_FIELD,2:BL_LEVELS)  ! IN Exchange coeff for FTL above
!                                       surface.
     &,GAMMA(BL_LEVELS)            ! IN Implicit weighting.

      REAL                         ! Split to avoid > 19 continuations.
     & ICE_FRACT(P_FIELD)          ! IN Fraction of grid-box which is
!                                       sea-ice (decimal fraction).
     &,LYING_SNOW(P_FIELD)         ! IN Lying snow (kg per sq m ie "mm")
     &,RADNET_C(P_FIELD,N_TYPES)   ! IN Area weighted ice component
!                                       of surface net radiation.
!                                       Modified for thermal canopy
!                                       over land.
!                                       (+ve downwards, W per sq m)
     &,RESFT(P_FIELD,N_TYPES)      ! IN Total resistance factor
     &,RHOKPM(P_FIELD)             ! IN Surface exchange coeff.
     &,RHOKPM_TILE(P_FIELD,N_TYPES)! IN Surface exchange coeff for tiles
     &,RHOKPM_POT_TILE(P_FIELD,N_TYPES)
!                                    IN Surface exchange coeff. for
!                                       potential evaporation.
     &,TILE_FRAC(P_FIELD,N_TYPES)  ! IN Tile fraction
     &,TIMESTEP                    ! IN Timestep in seconds.
     &,DTL_NT(P_FIELD,BL_LEVELS)   ! IN Non-turbulent increment for TL.
     &,DQW_NT(P_FIELD,BL_LEVELS)   ! IN Non-turbulent increment for QW.

       LOGICAL
     & LAND_MASK(P_FIELD)          ! IN T for land, F elsewhere.

!  Next 5 arrays are all IN as "explicit" fluxes from P243 (SF_EXCH and
!  possibly KMKH), and OUT as "implicit" fluxes.

      REAL
     & FQW(P_FIELD,BL_LEVELS)      ! INOUT Flux of QW (ie., for surface,
!                                          total evaporation). Kg/sq m/s
     &,FQW_TILE(P_FIELD,N_TYPES)   ! INOUT Tile flux of QW. Kg/sq m/s
     &,EPOT(P_FIELD)               ! INOUT Potential evaporation rate.
     &,EPOT_TILE(P_FIELD,N_TYPES)  ! INOUT Tile potential evaporation
!                                          rate.
     &,FTL(P_FIELD,BL_LEVELS)      ! INOUT Flux of TL (ie., for surface,
!                                          H/Cp where H is sensible heat
!                                          in W per sq m).
     &,FTL_TILE(P_FIELD,N_TYPES)   ! INOUT Tile flux of TL
     &,E_SEA(P_FIELD)              ! INOUT Evaporation from sea times
!                                          leads fraction (kg/m2/s).
!                                          Zero over land.
     &,H_SEA(P_FIELD)              ! INOUT Surface sensible heat flux ov
!                                          sea times leads fraction (W/m
!                                          Zero over land.
     &,QW(P_FIELD,BL_LEVELS)       ! INOUT Total water content (kg per
!                                          kg air).  From P243.
     &,RHOKE(P_FIELD,N_TYPES)      ! IN    Surface exchange coeff. for
!                                          FQW.
!                                      OUT =RHOKE to satisfy STASH.
     &,RHOKH_1(P_FIELD)            ! IN    Surface exchange coeffs for
!                                          FTL,
!                                      OUT =RHOKH_1 to satisfy STASH.
     &,TL(P_FIELD,BL_LEVELS)       ! INOUT Liquid/frozen water
!                                          temperature (K).  From P243.

      REAL
     & SURF_HT_FLUX(P_FIELD,N_TYPES)! OUT Net downward heat flux at
!                                        surface over land or sea-ice
!                                        fraction of gridbox (W/m2).
     &,DTRDZ_RML(P_FIELD)          ! OUT dz for the rapidly mixing layer
!                                      (needed in P245).
      INTEGER
     & NRML(P_FIELD)               ! IN The number of model layers
!                                       in the unstable rapidly mixing
!                                       layer. Zero if surface layer
!                                       is stable.

!  External references :-
      EXTERNAL TIMER

!  Local and other symbolic constants :-
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
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

C*L-----------COMDECK C_SICEHC FOR SUBROUTINE IMPL_CAL----------
C AI  = reciprocal effective areal heat capacity of sea-ice,
C          ( 1 / (J per sq m per K)).
      REAL AI

      PARAMETER(AI  = 4.8E-6)
C*----------------------------------------------------------------------


      REAL LS
      PARAMETER (
     & LS=LC+LF     ! Latent heat of sublimation (J per kg).
     &)

! Workspace :-
! 6*BL_LEVELS + 2 blocks of real workspace are required.
      REAL
     & AQ_AM(P_FIELD,BL_LEVELS)    ! As AQ: "Q" elements of matrix in
!                                     eqn P244.79 (modified during
!                                     Gaussian elimination process).
!                                     As AM: elements of matrix in eqn
!                                     P244.80 (also get modified).
     &,AQ_RML(P_FIELD)             ! Matrix element for humidity in
!                                     rapidly mixing layer. Then briefly
!                                     used for DELTAP on the UV grid.
     &,AT_ATQ(P_FIELD,BL_LEVELS)   ! Elements in atmospheric T rows of
!                                     matrix in eqn P244.79 (modified
!                                     during Gaussian elimination).
     &,AT_RML(P_FIELD)             ! Matrix element for temperature in
!                                     rapidly mixing layer.
     &,BPM(P_FIELD,N_TYPES)        ! Used in calculating elements of
!                                     TL(1) and QW(1) rows of matrix.
     &,BPM_GB(P_FIELD)             ! Used in calculating elements of
!                                     TL(1) and QW(1) rows of matrix.
     &,DELTAP(P_FIELD,BL_LEVELS)   ! -g.dt/dp for the rapidly mixing
!                                      layer
     &,DELTAP_RML(P_FIELD)         ! -g.dt/dp for the rapidly mixing
!                                      layer
     &,DELTA_QW(P_FIELD)           ! Increment in QW.
     &,DELTA_TL(P_FIELD)           ! Increment in TL.
     &,DQW_RML(P_FIELD)            ! Rapidly mixing layer increment
!                                     to QW.
     &,DQW(P_FIELD,BL_LEVELS)      ! Delta QW elements of vector
!                                     on RHS, then LHS, of eqn P244.79.
     &,DTL(P_FIELD,BL_LEVELS)      ! Delta TL (for atmosphere)
!                                     elements of vector on RHS, then
!                                     LHS, of eqn P244.79.
     &,DTL_RML(P_FIELD)            ! Delta TL for rapidly mixing layer.
     &,FQW_ICE(P_FIELD)            ! "Explicit" surface flux of QW for
!                                     sea-ice fraction of gridsquare.
     &,FTL_ICE(P_FIELD)            ! "Explicit" surface flux of TL for
!                                     sea-ice fraction of gridsquare.
     &,GAMMA_RKE_DQ(P_FIELD)       ! Gamma*rhoke*dq
     &,LAT_HEAT(P_FIELD)           ! Latent heat of evaporation for
!                                     snow-free land or sublimation for
!                                     snow-covered land.
     &,RHOKE_PM_GB(P_FIELD)        ! tile avg of RHOKE*RHOKPM

!  Local scalars :-
      REAL
     & CTQ      ! Matrix element in P244.??, for local increments to rml
     &,CQ       ! Matrix element in "Q" row in eqn P244.79.
     &,CQ_RML   ! As above but for rapidly mixing layer increment.
     &,CT       ! Matrix element in "T" row in eqn P244.79.
     &,CT_RML   ! As above but for rapidly mixing layer increment.
     &,RBTQ     ! Reciprocal of B P244.??, for local increments to rml
     &,RBQ      ! Reciprocal of BQ(') (eqns P244.98, 101, 104).
     &,RBQ_RML  ! As above but for the rapidly mixing layer increment.
     &,RBT      ! Reciprocal of BT' (eqns P244.107, 110, 113).
     &,RBT_RML  ! As above but for the rapidly mixing layer increment.
     &,temp1    ! temporary
     &,temp2    ! temporary
     &,temp3    ! temporary 

      INTEGER
     & BLM1     ! BL_LEVELS minus 1.
     &,NRMLP1   ! NRML plus 1.
     &,I        ! Loop counter (horizontal field index).
     &,L        ! Loop counter (horizontal field index).
     &,ITILE    ! Loop counter (land surface tile index).
     &,K        ! Loop counter (vertical index).

!-----------------------------------------------------------------------
!!  0.  Check that the scalars input to define the grid are consistent.
!       See comments to routine SF_EXCH for details.
!-----------------------------------------------------------------------

      IF (LTIMER) THEN
        CALL TIMER('IMCALTQ ',3)
      ENDIF

      BLM1 = BL_LEVELS-1

!-----------------------------------------------------------------------
!! (A) Calculations on P-grid.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!! 1.  Calculate implicit T and Q increments due to local mixing within
!!     the rapidly mixing layer (where it exists).
!!     The surface fluxes FTL(I,1), FQW(I,1) are used for calculating
!!     the rapidly mixing layer (rml) increments but not here.
!!     Therefore the matrix equation we must solve to find the implicit
!!     T and Q increments due to local mixing within the rml does not
!!     have a "surface" row and we can solve for the T and Q increments
!!     for K = 1 to NRML simultaneously.
!-----------------------------------------------------------------------

      DO K = 1,BL_LEVELS
        DO I = P1,P1+P_POINTS-1
          DELTAP(I,K) = -G * TIMESTEP/DTRDZ(I,K)
        ENDDO ! Loop over p-points
      ENDDO ! Loop over bl-levels

!-----------------------------------------------------------------------
!! 1.1 Start 'upward sweep' with lowest model layer, which will be the
!!     bottom of the rapidly mixing layer (rml) if it exists.
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1

        RHOKE_PM_GB(I) = 0.0  ! set to zero
        BPM_GB(I) = 0.0

        LAT_HEAT(I) = LC
        IF (LAND_MASK(I)) THEN
          IF (LYING_SNOW(I).GT.0.0) LAT_HEAT(I) = LS
        ELSE
          IF (ICE_FRACT(I).GT.0.0) LAT_HEAT(I) = LS
        ENDIF
        IF (NRML(I) .GE. 2) THEN

!  "Explicit" increments due to local mixing within the rml.
!  P244.49/31 but surface flux used in rml increment calculations.

! Add non-turbulent increments here
          DQW(I,1) = -DTRDZ(I,1) * FQW(I,2) + DQW_NT(I,1)

          DTL(I,1) = -DTRDZ(I,1) * FTL(I,2) + DTL_NT(I,1)

!  Define matrix elements (CTQ always zero for this case).

          AT_ATQ(I,1) = -DTRDZ(I,1) * GAMMA(2)*RHOKH(I,2)*RDZ(I,2)
!                                                        ! P244.28
          RBTQ = 1.0 / ( 1.0 - AT_ATQ(I,1) ) ! Reciprocal of P244.110

!  Now start Gaussian elimination

          DQW(I,1) = RBTQ * DQW(I,1)                  ! P244.102
          DTL(I,1) = RBTQ * DTL(I,1)                  ! P244.111

          AT_ATQ(I,1) = RBTQ * AT_ATQ(I,1)                    ! P244.112

!  Start calculating DELTAP_RML. Mid-level depths added in 2.2 below.

          DELTAP_RML(I) = DELTAP(I,1)
        ELSE ! No rapidly mixing layer calculations
          DTRDZ_RML(I) = 1.E30
          DQW_RML(I) = 1.E30
          DTL_RML(I) = 1.E30
          AQ_RML(I) = 1.E30
          AT_RML(I) = 1.E30
        ENDIF
      ENDDO ! Loop over p_points

!-----------------------------------------------------------------------
!! 2.2 Continue upward sweep through middle of the rapidly mixing layer
!!     (if it exists) and to its top. NB NRML is always < temp= BLM1.
!-----------------------------------------------------------------------

      DO K=2,BLM1
        DO I=P1,P1+P_POINTS-1

!   If in the top rapidly mixing layer then do not include flux at its
!   top in the calculation ie FQW(I,NRML+1) and FTL(I,NRML+1) are not
!   included here; they are taken account of in the non-local mixing
!   through the "rapidly mixing layer".

          IF ( K .EQ. NRML(I) ) THEN

!   Add final DELTAP contribution to DELTAP_RML and then calculate
!   DTRDZ_RML.  Lower level contributions added in 2.1 and below.

            DELTAP_RML(I) = DELTAP_RML(I) + DELTAP(I,K)
            DTRDZ_RML(I) =-G * TIMESTEP / DELTAP_RML(I)


!  "Explicit" flux divergence across layer giving explicit
!  increment due to the local mixing at the top of rml.

            DQW(I,K) = DTRDZ(I,K) * FQW(I,K) + DQW_NT(I,K)

            DTL(I,K) = DTRDZ(I,K) * FTL(I,K) + DTL_NT(I,K)


!  Define matrix elements (A always zero for this case).

            CTQ = -DTRDZ(I,K) * GAMMA(k)*RHOKH(I,K)*RDZ(I,K)  ! P244.36
            RBTQ = 1.0 / ( 1.0  - CTQ*( 1.0 + AT_ATQ(I,K-1) ) )
!                                            ... Reciprocal of P244.113
!  Now start Gaussian elimination

            DQW(I,K) = RBTQ * ( DQW(I,K) - CTQ*DQW(I,K-1) )
            DTL(I,K) = RBTQ * ( DTL(I,K) - CTQ*DTL(I,K-1) )
          ELSEIF (K .LT. NRML(I)) THEN

!  Add layer depths to form total rml depth.

            DELTAP_RML(I) = DELTAP_RML(I) + DELTAP(I,K)

!  "Explicit" flux divergence across layer giving explicit
!  increment due to the local mixing.P244.54/38

            DQW(I,K) = -DTRDZ(I,K) * ( FQW(I,K+1) - FQW(I,K) )
     &                + DQW_NT(I,K)

            DTL(I,K) = -DTRDZ(I,K) * ( FTL(I,K+1) - FTL(I,K) )
     &                + DTL_NT(I,K)

!  Define matrix elements.

            AT_ATQ(I,K) = -DTRDZ(I,K) *
     &      GAMMA(K+1)*RHOKH(I,K+1)*RDZ(I,K+1)
            CTQ = -DTRDZ(I,K) * GAMMA(k)*RHOKH(I,K)*RDZ(I,K)! P244.36
            RBTQ = 1.0 / ( 1.0  - AT_ATQ(I,K)
     &                          - CTQ * ( 1.0 + AT_ATQ(I,K-1) ) )
!                                             ... Reciprocal of P244.113
!  Now start Gaussian elimination

            DQW(I,K) = RBTQ * ( DQW(I,K) - CTQ*DQW(I,K-1) )
            DTL(I,K) = RBTQ * ( DTL(I,K) - CTQ*DTL(I,K-1) )
            AT_ATQ(I,K) = RBTQ * AT_ATQ(I,K)                 ! P244.115
          ENDIF
        ENDDO  !p_points
      ENDDO !blm1

!-----------------------------------------------------------------------
!! 2.3 Downward sweep through matrix. Add implicit increments due to
!!     local mixing within the rapidly mixing layer.  Update fluxes of
!!     heat and moisture at the surface and the top-of-rml using
!!     local mixing increments for layers 1 and NRML respectively.
!-----------------------------------------------------------------------

      DO K=BLM1,1,-1
        DO I=P1,P1+P_POINTS-1
          IF ((NRML(I) .GE. 2) .AND. (K .EQ. NRML(I))) THEN
            QW(I,K) = QW(I,K) + DQW(I,K)                   ! P244.128
            TL(I,K) = TL(I,K) + DTL(I,K)                   ! P244.127
            FQW(I,K+1) = FQW(I,K+1)
     &                 + GAMMA(K+1)*RHOKH(I,K+1)*RDZ(I,K+1)*DQW(I,K)
            FTL(I,K+1) = FTL(I,K+1)
     &                + GAMMA(K+1)*RHOKH(I,K+1)*RDZ(I,K+1)*DTL(I,K)
          ELSEIF ((NRML(I) .GE. 2) .AND. (K .LT. NRML(I))) THEN
            DQW(I,K) = DQW(I,K)
     &                           - AT_ATQ(I,K)*DQW(I,K+1)  ! P244.???
            DTL(I,K) = DTL(I,K)
     &                           - AT_ATQ(I,K)*DTL(I,K+1)  ! P244.???
            QW(I,K) = QW(I,K) + DQW(I,K)                   ! P244.128
            TL(I,K) = TL(I,K) + DTL(I,K)                   ! P244.127
          ENDIF

          IF ((NRML(I) .GE. 2) .AND. (K .EQ. 1)) THEN

            IF (LAND_MASK(I)) THEN

              DO ITILE=1,N_TYPES
                GAMMA_RKE_DQ(I) = GAMMA(1)*RHOKE(I,ITILE) *
     &                     (DQW(I,1) - ALPHA1(I,ITILE)*DTL(I,1))

                FTL_TILE(I,ITILE) = FTL_TILE(I,ITILE) +
     &              RHOKPM_TILE(I,ITILE)*(LAT_HEAT(I)*
     &              GAMMA_RKE_DQ(I)- GAMMA(1)*ASHTF(I)*DTL(I,1))

                FQW_TILE(I,ITILE) = FQW_TILE(I,ITILE) -
     &                      RHOKPM_TILE(I,ITILE) * (CP*GAMMA_RKE_DQ(I) +
     &                      GAMMA(1)*RESFT(I,ITILE)*ASHTF(I)*DQW(I,1))
                EPOT_TILE(I,ITILE) = EPOT_TILE(I,ITILE) -
     &                      RHOKPM_POT_TILE(I,ITILE) * 
     &                      (CP*GAMMA_RKE_DQ(I) +
     &                      GAMMA(1)*ASHTF(I)*DQW(I,1))

              ENDDO ! ITILE

            ELSEIF (ICE_FRACT(I) .GT. 0.0) THEN

              GAMMA_RKE_DQ(I) = GAMMA(1)*RHOKE(I,1) *
     &                     (DQW(I,1) - ALPHA1(I,1)*DTL(I,1))

              FQW_ICE(I) = FQW(I,1) - E_SEA(I) - RHOKPM_TILE(I,1) *
     &             (CP*GAMMA_RKE_DQ(I) + GAMMA(1)*ASHTF(I)*DQW(I,1))
              FTL_ICE(I) = FTL(I,1) - H_SEA(I)/CP + RHOKPM_TILE(I,1) *
     &             (LS*GAMMA_RKE_DQ(I) - GAMMA(1)*ASHTF(I)*DTL(I,1))
              E_SEA(I) = E_SEA(I) - (1.0 - ICE_FRACT(I)) *
     &                                  GAMMA(1)*RHOKH_1(I)*DQW(I,1)
              H_SEA(I) = H_SEA(I) - (1.0 - ICE_FRACT(I)) * CP *
     &                                  GAMMA(1)*RHOKH_1(I)*DTL(I,1)
              FTL_TILE(I,1) = FTL_ICE(I) + H_SEA(I)/CP
              FQW_TILE(I,1) = FQW_ICE(I) + E_SEA(I)
              EPOT_TILE(I,1) = FQW_ICE(I) + E_SEA(I)

            ELSE ! ordinary sea point
              FTL_TILE(I,1) = FTL(I,1) -
     &                            GAMMA(1)*RHOKH_1(I)*DTL(I,1)
              FQW_TILE(I,1) = FQW(I,1) -
     &                            GAMMA(1)*RHOKH_1(I)*DQW(I,1)
              EPOT_TILE(I,1) = EPOT(I) -
     &                            GAMMA(1)*RHOKH_1(I)*DQW(I,1)


            ENDIF ! land/ice/sea

          ENDIF ! nrml >= 2 and k=1

! Reset level 1 fluxes to zero before reaveraging tiles

        ENDDO  ! p_points
      ENDDO  !blm,-1,1


      DO I=P1,P1+P_POINTS-1
        IF(LAND_MASK(I)) THEN
          FTL(I,1)=0.0
          FQW(I,1)=0.0
          EPOT(I)=0.0
        ELSE
          FTL(I,1)=FTL_TILE(I,1)
          FQW(I,1)=FQW_TILE(I,1)
          RHOKE_PM_GB(I)=RHOKE(I,1)*RHOKPM_TILE(I,1)
        ENDIF
      ENDDO


      DO ITILE=1,N_TYPES
        DO I=P1,P1+P_POINTS-1
          IF(LAND_MASK(I)) THEN
            FTL(I,1)=FTL(I,1)+FTL_TILE(I,ITILE)*TILE_FRAC(I,ITILE)
            FQW(I,1)=FQW(I,1)+FQW_TILE(I,ITILE)*TILE_FRAC(I,ITILE)
            EPOT(I)=EPOT(I)+EPOT_TILE(I,ITILE)*TILE_FRAC(I,ITILE)
            RHOKE_PM_GB(I) = RHOKE_PM_GB(I) + TILE_FRAC(I,ITILE)*
     &                       RHOKE(I,ITILE)*RHOKPM_TILE(I,ITILE)

          ELSE
            FTL_TILE(I,ITILE) = FTL(I,1)
            FQW_TILE(I,ITILE) = FQW(I,1)
            EPOT_TILE(I,ITILE) = EPOT(I)
          ENDIF
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!! 3.  Calculate those matrix and vector elements on the LHS of eqn
!!     P244.79 which are to do with implicit solution of the moisture
!!     transport problem at the surface, above the rml (if it exists)
!!     and between all levels if it does not.
!!     Begin with "upward sweep" through lower half of matrix).
!-----------------------------------------------------------------------
!! 3.1 Row of matrix applying to QW transport into top "boundary"
!!     layer of model atmosphere.
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1
! Include non-turbulent increments.
        DQW(I,BL_LEVELS) = DTRDZ(I,BL_LEVELS) * FQW(I,BL_LEVELS)
     &                     +DQW_NT(I,BL_LEVELS)
!                                                            ... P244.58

        AQ_AM(I,BL_LEVELS) = -DTRDZ(I,BL_LEVELS)               ! P244.56
     &       * GAMMA(BL_LEVELS)*RHOKH(I,BL_LEVELS)*
     &          RDZ(I,BL_LEVELS)

        RBQ = 1.0 / ( 1.0 - AQ_AM(I,BL_LEVELS) )    ! Reciprocal P244.98
        DQW(I,BL_LEVELS) = RBQ * DQW(I,BL_LEVELS)   ! P244.99
        AQ_AM(I,BL_LEVELS) = RBQ * AQ_AM(I,BL_LEVELS)         ! P244.100
      ENDDO

!-----------------------------------------------------------------------
!! 3.2 Rows of matrix applying to "middle of boundary layer" model
!!     layers, i.e. all but the topmost and bottom layers.
!-----------------------------------------------------------------------

      DO K=BLM1,2,-1
        DO I=P1,P1+P_POINTS-1

!  "Explicit" flux divergence across layer giving explicit QW increment.

          IF ( K .GT. NRML(I) ) THEN
            DQW(I,K) = -DTRDZ(I,K) * ( FQW(I,K+1) - FQW(I,K) )
     &                + DQW_NT(I,K)
!                                                            ... P244.54

            CQ = -DTRDZ(I,K) * GAMMA(K+1)*RHOKH(I,K+1)*RDZ(I,K+1)
!                                                              ! P244.52
            AQ_AM(I,K) = -DTRDZ(I,K) * GAMMA(K)*RHOKH(I,K)*RDZ(I,K)
!                                                              ! P244.51
            RBQ = 1.0 / ( 1.0 - AQ_AM(I,K) - CQ*( 1.0 + AQ_AM(I,K+1) ) )
!                       1                       2                    2 1
!                                             ... reciprocal of P244.101

! Now include non-turbulent increments.
            DQW(I,K) = RBQ * (DQW(I,K) - CQ*DQW(I,K+1) )      ! P244.102

            AQ_AM(I,K) = RBQ * AQ_AM(I,K)                     ! P244.103
          ENDIF
        ENDDO ! P_points
      ENDDO !blm1,2,-1

!-----------------------------------------------------------------------
!! 3.3 Bottom model layer QW row of matrix equation.
!-----------------------------------------------------------------------


      DO I=P1,P1+P_POINTS-1

        DO ITILE=1,N_TYPES
          IF (LAND_MASK(I)) THEN
            BPM(I,ITILE) = GAMMA(1)*ASHTF(I)*RHOKPM_TILE(I,ITILE)
          ELSE
            BPM(I,ITILE) = (1. - ICE_FRACT(I))*GAMMA(1)*RHOKH_1(I) +
     &                      GAMMA(1)*ASHTF(I)*RHOKPM_TILE(I,1)
          ENDIF

          BPM_GB(I)=BPM_GB(I) + BPM(I,ITILE) * TILE_FRAC(I,ITILE)

        ENDDO !n_types

        IF ( NRML(I) .GE. 2 ) THEN

!-----------------------------------------------------------------------
!! 3.3.1 Start calculating rapidly mixing layer increments.
!-----------------------------------------------------------------------

          NRMLP1 = NRML(I) + 1

!  "Explicit" QW increment for the rapidly mixing layer.

          DQW_RML(I) = -DTRDZ_RML(I) * ( FQW(I,NRMLP1) - FQW(I,1) )

!  Define coefficients A,B,C, for implicit calculations.

          AQ_RML(I) = - DTRDZ_RML(I)*GAMMA(1)*CP*RHOKE_PM_GB(I)

          CQ_RML = -DTRDZ_RML(I) * GAMMA(NRMLP1)
     &             *RHOKH(I,NRMLP1)*RDZ(I,NRMLP1)

          RBQ_RML=0.0

          IF (LAND_MASK(I)) THEN
            DO ITILE=1,N_TYPES

              RBQ_RML = RBQ_RML + TILE_FRAC(I,ITILE) *
     &              1./(1. - AQ_RML(I) - CQ_RML*(1. + AQ_AM(I,NRMLP1))
     &               + DTRDZ_RML(I)*RESFT(I,ITILE)*BPM(I,ITILE) )

            ENDDO ! n_types
          ELSE
            RBQ_RML = 1./(1. - AQ_RML(I)
     &                - CQ_RML*(1. + AQ_AM(I,NRMLP1))
     &                + DTRDZ_RML(I)*RESFT(I,1)*BPM_GB(I) )
          ENDIF

          DQW_RML(I) = RBQ_RML * ( DQW_RML(I)
     &                                 - CQ_RML * DQW(I,NRMLP1) )
          AQ_RML(I) = RBQ_RML * AQ_RML(I)
        ELSE

!  "Explicit" increment for QW(1) when there is no rapidly mixing
!  layer or when it is only one model layer in depth.

          DQW(I,1) = - DTRDZ(I,1)*(FQW(I,2) - FQW(I,1)) + DQW_NT(I,1)

          AQ_AM(I,1) = - DTRDZ(I,1)*GAMMA(1)*CP*RHOKE_PM_GB(I)

          CQ = - DTRDZ(I,1)*GAMMA(2)*RHOKH(I,2)*RDZ(I,2)

          RBQ=0.0

          DO ITILE=1,N_TYPES
            IF (LAND_MASK(I)) THEN
              RBQ = RBQ + TILE_FRAC(I,ITILE) *
     &                    1./( 1. - AQ_AM(I,1) - CQ*(1. + AQ_AM(I,2))
     &                    + DTRDZ(I,1)*RESFT(I,ITILE)*BPM(I,ITILE) )
            ELSEIF(ITILE .EQ. 1) THEN
              RBQ =  1./( 1. - AQ_AM(I,1) - CQ*(1. + AQ_AM(I,2))
     &                + DTRDZ(I,1)*RESFT(I,1)*BPM(I,1) )
            ENDIF
          ENDDO ! n_types

          DQW(I,1) = RBQ*(DQW(I,1)  - CQ*DQW(I,2))
          AQ_AM(I,1) = RBQ*AQ_AM(I,1)
        ENDIF
      ENDDO

!-----------------------------------------------------------------------
!! 4.  Calculate those matrix and vector elements on the LHS of eqn
!!     P244.79 which are to do with implicit solution of the heat
!!     transport problem (i.e. for ice/liquid
!!     water temperatures in atmospheric boundary layer), and begin the
!!     solution algorithm (perform "upward sweep" through upper half of
!!     matrix).
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!! 4.2 Lowest atmospheric layer TL row of matrix.
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1
      IF (NRML(I) .GE. 2) THEN
       NRMLP1 = NRML(I) + 1

!  "Explicit" rapidly mixing layer increment for TL.

          DTL_RML(I) = - DTRDZ_RML(I) * ( FTL(I,NRMLP1) - FTL(I,1) )
          AT_RML(I) = - DTRDZ_RML(I)*GAMMA(NRMLP1)*
     &                  RHOKH(I,NRMLP1)*RDZ(I,NRMLP1)

          CT_RML = - DTRDZ_RML(I)*GAMMA(1)*RHOKE_PM_GB(I)*LAT_HEAT(I)

          RBT_RML=0.0
          IF(LAND_MASK(I)) THEN
            DO ITILE=1,N_TYPES

              RBT_RML = RBT_RML + TILE_FRAC(I,ITILE) *
     &                1./(1. - AT_RML(I) - ALPHA1(I,ITILE)*CT_RML*
     &                 (1.+AQ_RML(I)) + DTRDZ_RML(I)*BPM(I,ITILE) )
            ENDDO
          ELSE
              RBT_RML = 1./(1. - AT_RML(I) - ALPHA1(I,1)*CT_RML*
     &                 (1.+AQ_RML(I)) + DTRDZ_RML(I)*BPM(I,1) )
          ENDIF

          DTL_RML(I) = RBT_RML*(DTL_RML(I) - CT_RML*DQW_RML(I))
          AT_RML(I) = RBT_RML * AT_RML(I)
        ELSE

!  "Explicit" increment to TL(1) when there is no rapidly mixing layer
!  or it does not extend beyond the bottom model layer.

          DTL(I,1) = - DTRDZ(I,1)*(FTL(I,2) - FTL(I,1))
     &               + DTL_NT(I,1)
          CT = - DTRDZ(I,1)*GAMMA(1)*RHOKE_PM_GB(I)*LAT_HEAT(I)

          AT_ATQ(I,1) = - DTRDZ(I,1)*GAMMA(2)*RHOKH(I,2)*RDZ(I,2)

          RBT=0.0

          DO ITILE=1,N_TYPES

            IF(LAND_MASK(I)) THEN
              RBT = RBT + TILE_FRAC(I,ITILE) *
     &                1./( 1. - AT_ATQ(I,1) - ALPHA1(I,ITILE)*CT*
     &                (1. + AQ_AM(I,1)) + DTRDZ(I,1)*BPM(I,ITILE) )
            ELSEIF(ITILE .EQ.1) THEN
              RBT = 1./( 1. - AT_ATQ(I,1) - ALPHA1(I,1)*CT*
     &                (1. + AQ_AM(I,1)) + DTRDZ(I,1)*BPM(I,1) )
            ENDIF
          ENDDO

          DTL(I,1) = RBT*(DTL(I,1) - CT*DQW(I,1))
          AT_ATQ(I,1) = RBT*AT_ATQ(I,1)
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
!! 4.3 Rows of matrix applying to TL transport into model layers in the
!!     "middle" of the "boundary" layer, i.e. all but the bottom layer
!!     and the top "boundary" layer.
!-----------------------------------------------------------------------
      DO K=2,BLM1
        DO I=P1,P1+P_POINTS-1

!   "Explicit" flux divergence across layer giving explicit TL increment
!   due to mixing above rml if it exists or for all levels if it does
!   not.

          NRMLP1 = NRML(I) + 1
          IF (K .GT. NRML(I)) THEN
            DTL(I,K) = -DTRDZ(I,K) * ( FTL(I,K+1) - FTL(I,K) )
     &                               +DTL_NT(I,K)
!                                                            ... P244.38

            AT_ATQ(I,K) = -DTRDZ(I,K)*GAMMA(K+1)*
     &      RHOKH(I,K+1)*RDZ(I,K+1)                            ! P244.35
            CT = -DTRDZ(I,K) * GAMMA(K)*RHOKH(I,K)*RDZ(I,K)  ! P244.36
            IF ((NRML(I) .GE. 2) .AND. (K .EQ. NRMLP1)) THEN
              RBT = 1.0 / ( 1.0 - AT_ATQ(I,K) - CT*( 1.0 + AT_RML(I) ) )
              DTL(I,K) = RBT * ( DTL(I,K) - CT*DTL_RML(I) )
            ELSE
              RBT = 1.0 / ( 1.0 - AT_ATQ(I,K)
     &                          - CT*( 1.0 + AT_ATQ(I,K-1) ) )
!                         1          2                     2 1
!                                             ... Reciprocal of P244.113

              DTL(I,K) = RBT * ( DTL(I,K) - CT*DTL(I,K-1) )
!                                                           ... P244.114
            ENDIF
            AT_ATQ(I,K) = RBT * AT_ATQ(I,K)                   ! P244.115
          ENDIF
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!! 4.4 Top "boundary" layer TL row of matrix.  TL for this layer can
!!     then be, and is, updated.
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1
        DTL(I,BL_LEVELS) = DTRDZ(I,BL_LEVELS) * FTL(I,BL_LEVELS)
     &                   + DTL_NT(I,BL_LEVELS)
!                                                            ... P244.42
        CT = -DTRDZ(I,BL_LEVELS) * GAMMA(BL_LEVELS)*
     &        RHOKH(I,BL_LEVELS)*RDZ(I,BL_LEVELS)

        IF (NRML(I) .EQ. BLM1) THEN
          RBT = 1.0 / ( 1.0 - CT*( 1.0 + AT_RML(I) ) )
          DTL(I,BL_LEVELS) = RBT * ( DTL(I,BL_LEVELS)
     &                                  - CT*DTL_RML(I) )
        ELSE
         RBT = 1.0 / ( 1.0 - CT*( 1.0 + AT_ATQ(I,BLM1) ) )
!                                             ... Reciprocal of P244.116

          DTL(I,BL_LEVELS) = RBT * ( DTL(I,BL_LEVELS)
     &              - CT*DTL(I,BLM1) )   ! P244.117
        ENDIF

        TL(I,BL_LEVELS) = TL(I,BL_LEVELS) + DTL(I,BL_LEVELS) ! P244.127

      ENDDO

!-----------------------------------------------------------------------
!! 5.  "Downward sweep" through whole matrix.  TL, QW and
!!     updated when the final implicit increments have been calculated.
!-----------------------------------------------------------------------
!! 5.1 Remaining TL rows of matrix and add implicit increments above
!!     the rml or at all levels if it is less than two layers thick.
!-----------------------------------------------------------------------

      DO K=BLM1,1,-1
        DO I=P1,P1+P_POINTS-1
          IF ( (K .GT. NRML(I)) .OR. (NRML(I) .LT. 2) ) THEN
            DTL(I,K) = DTL(I,K) - AT_ATQ(I,K)*DTL(I,K+1)
                                                              ! P244.118
            TL(I,K) = TL(I,K) + DTL(I,K)                   ! P244.127
          ENDIF
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!! 5.2 Rapidly mixing layer increments
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1
        IF ( NRML(I) .GE. 2 ) THEN
          NRMLP1 = NRML(I) + 1
          DTL_RML(I) = DTL_RML(I) - AT_RML(I) * DTL(I,NRMLP1)
          DQW_RML(I) = DQW_RML(I)
     &                      - AQ_RML(I) * ALPHA1_GB(I) * DTL_RML(I)
          TL(I,1) = TL(I,1) + DTL_RML(I)
          QW(I,1) = QW(I,1) + DQW_RML(I)
        ELSE

!-----------------------------------------------------------------------
!! 5.3 Lowest-level QW row of matrix; local mixing where there is no rml
!-----------------------------------------------------------------------

          DQW(I,1) = DQW(I,1) - ALPHA1_GB(I)*AQ_AM(I,1)*DTL(I,1)
          QW(I,1) = QW(I,1) + DQW(I,1)     ! P244.128
        ENDIF
      ENDDO

!-----------------------------------------------------------------------
!! 5.4 Remaining QW rows of matrix + updating of QW's.
!!     Add implicit increments due to mixing above rml or at all levels
!!     if there it does not exist.
!-----------------------------------------------------------------------

      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          IF (K .GT. NRML(I)) THEN

            IF ((NRML(I) .GE. 2) .AND. (K-1 .EQ. NRML(I))) THEN
              DQW(I,K) = DQW(I,K) - AQ_AM(I,K)*DQW_RML(I)
            ELSE
              DQW(I,K) = DQW(I,K) - AQ_AM(I,K)*DQW(I,K-1) ! P244.121
            ENDIF

            QW(I,K) = QW(I,K) + DQW(I,K)                  ! P244.128
          ELSE

!  Add the increments due to rapid mixing if in the rapidly mixing layer

            TL(I,K) = TL(I,K) + DTL_RML(I)
            QW(I,K) = QW(I,K) + DQW_RML(I)
          ENDIF

        ENDDO !p_points
      ENDDO !bl_levels

!-----------------------------------------------------------------------
!! 6.  Calculate final implicit fluxes of heat and moisture.
!-----------------------------------------------------------------------
!! 6.1 Surface fluxes for the 3 surface types: land, sea-ice, ordinary
!!     sea. Pass out the value of RHOKH(,1) in GAMMA*RHOKH_1 to satisfy
!!     STASH GAMMA*RHOKH_RDZ will contain precisely that on output.
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1
        IF ( NRML(I) .GE. 2 ) THEN
           DELTA_TL(I) =  DTL_RML(I)
           DELTA_QW(I) = DQW_RML(I)
        ELSE
           DELTA_TL(I) = DTL(I,1)
           DELTA_QW(I) = DQW(I,1)
        ENDIF

        GAMMA_RKE_DQ(I) = GAMMA(1)*RHOKE(I,1)*
     &    (DELTA_QW(I) - ALPHA1(I,1)*DELTA_TL(I))
      ENDDO !p_points


      DO ITILE=1,N_TYPES

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO L=LAND1,LAND1+LAND_PTS-1
           I = LAND_INDEX(L)

! overwrite gamma_rke_dq for land points only
            GAMMA_RKE_DQ(I) = GAMMA(1)*RHOKE(I,ITILE)*
     &                  (DELTA_QW(I) - ALPHA1(I,ITILE)*DELTA_TL(I))

            temp1= RHOKPM_TILE(I,ITILE)*( LAT_HEAT(I)*GAMMA_RKE_DQ(I) -
     &                                   GAMMA(1)*ASHTF(I)*DELTA_TL(I) )

            FTL_TILE(I,ITILE) = FTL_TILE(I,ITILE) + temp1
            FTL(I,1) = FTL(I,1) + temp1 * TILE_FRAC(I,ITILE)


            temp2= RHOKPM_TILE(I,ITILE)*( CP*GAMMA_RKE_DQ(I) +
     &             GAMMA(1)*RESFT(I,ITILE)*ASHTF(I)*DELTA_QW(I) )

            FQW_TILE(I,ITILE) = FQW_TILE(I,ITILE) - temp2
            FQW(I,1) = FQW(I,1) - temp2 * TILE_FRAC(I,ITILE)
            temp3= RHOKPM_TILE(I,ITILE)*( CP*GAMMA_RKE_DQ(I) +
     &             GAMMA(1)*ASHTF(I)*DELTA_QW(I) )

            EPOT_TILE(I,ITILE) = EPOT_TILE(I,ITILE) - temp3
            EPOT(I) = EPOT(I) - temp3 * TILE_FRAC(I,ITILE)

            SURF_HT_FLUX(I,ITILE) = RADNET_C(I,ITILE) -           
     &         LAT_HEAT(I)*FQW_TILE(I,ITILE) - CP*FTL_TILE(I,ITILE)



            IF(ITILE .EQ. 1) THEN
              FTL_ICE(I) = 0.0
              FQW_ICE(I) = 0.0
            ENDIF

          ENDDO ! land points
        ENDDO ! End of Tile loop


        DO I=P1,P1+P_POINTS-1

        IF (.NOT. LAND_MASK(I) .AND. ICE_FRACT(I).GT.0.0) THEN
          FQW_ICE(I) = FQW(I,1) - E_SEA(I) - RHOKPM_TILE(I,1) *
     &            (CP*GAMMA_RKE_DQ(I) + GAMMA(1)*ASHTF(I)*DELTA_QW(I))
          FTL_ICE(I) =  FTL(I,1) - H_SEA(I)/CP + RHOKPM_TILE(I,1) *
     &            (LS*GAMMA_RKE_DQ(I) - GAMMA(1)*ASHTF(I)*DELTA_TL(I))
          E_SEA(I) = E_SEA(I) - (1.0 - ICE_FRACT(I)) *
     &                               GAMMA(1)*RHOKE(I,1)*DELTA_QW(I)
          H_SEA(I) = H_SEA(I) - (1.0 - ICE_FRACT(I)) * CP *
     &                               GAMMA(1)*RHOKH_1(I)*DELTA_TL(I)
          FTL(I,1) = FTL_ICE(I) + H_SEA(I)/CP
          FQW(I,1) = FQW_ICE(I) + E_SEA(I)

          FTL_TILE(I,1)=FTL(I,1)
          FQW_TILE(I,1)=FQW(I,1)
          EPOT(I) = FQW_ICE(I) + E_SEA(I)
          EPOT_TILE(I,1)=EPOT(I)


          SURF_HT_FLUX(I,1) = RADNET_C(I,1) - LS*FQW_ICE(I) -           
     &                                    CP*FTL_ICE(I)

! ordinary sea point
        ELSEIF(.NOT. LAND_MASK(I) .AND. ICE_FRACT(I).EQ.0.0) THEN

          FTL(I,1) = FTL(I,1) - GAMMA(1)*RHOKH_1(I)*DELTA_TL(I)
          FQW(I,1) = FQW(I,1) - GAMMA(1)*RHOKH_1(I)*DELTA_QW(I)

          FTL_TILE(I,1)=FTL(I,1)
          FQW_TILE(I,1)=FQW(I,1)
          EPOT(I) = EPOT(I) - GAMMA(1)*RHOKH_1(I)*DELTA_QW(I)
          EPOT_TILE(I,1)=EPOT(I)

          E_SEA(I) = FQW(I,1)
          H_SEA(I) = FTL(I,1)*CP

          FTL_ICE(I) = 0.0
          FQW_ICE(I) = 0.0

          SURF_HT_FLUX(I,1) = RADNET_C(I,1) - LC*E_SEA(I) - H_SEA(I)    

        ENDIF
      ENDDO

C-----------------------------------------------------------------------
CL  Ensures that the potential evaporation rate equals the evaporation
CL  rate, when there is net condensation. Otherwise E/Ep could be
CL  <0 or >1 when the implicit adjustment is added.
C-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1      
        IF(FQW(I,1).LT.0.0)THEN
          EPOT(I)=FQW(I,1)
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
!! 6.2 Fluxes at layer interfaces above the surface.
!-----------------------------------------------------------------------
      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1

!  Calculate and store fluxes due to local mixing.
!  FTL(local mixing) stored in array AT,
!  FQW(local mixing) stored in array AQ_AM.

          NRMLP1 = NRML(I) + 1
          IF ((NRML(I) .GE. 2) .AND. (K .EQ. NRMLP1)) THEN
            AT_ATQ(I,K) = FTL(I,K) - GAMMA(K)*RHOKH(I,K)*RDZ(I,K)
     &                              * ( DTL(I,K) - DTL_RML(I) )
            AQ_AM(I,K) = FQW(I,K) - GAMMA(K)*RHOKH(I,K)*RDZ(I,K)
     &                              * ( DQW(I,K) - DQW_RML(I) )
          ELSE
            AT_ATQ(I,K) = FTL(I,K) - GAMMA(K)*RHOKH(I,K)*RDZ(I,K)
     &                              * ( DTL(I,K) - DTL(I,K-1) )
            AQ_AM(I,K) = FQW(I,K) - GAMMA(K)*RHOKH(I,K)*RDZ(I,K)
     &                              * ( DQW(I,K) - DQW(I,K-1) )
          ENDIF

!  Now calculate the implicit fluxes including both local mixing and
!  if appropriate also the fluxes due to rapid mixing through layers.

          IF ( K .EQ. 2 ) THEN
            IF ( NRML(I) .GE. 2 ) THEN
              FTL(I,K) = AT_ATQ(I,K)
     &                    + FTL(I,K-1) - DTL_RML(I) / DTRDZ(I,K-1)
              FQW(I,K) = AQ_AM(I,K)
     &                    + FQW(I,K-1) - DQW_RML(I) / DTRDZ(I,K-1)
            ELSE
              FTL(I,K) = AT_ATQ(I,K)
              FQW(I,K) = AQ_AM(I,K)
            ENDIF
          ELSEIF ( K .LE. NRML(I) ) THEN
            FTL(I,K) = AT_ATQ(I,K) - AT_ATQ(I,K-1)
     &                    + FTL(I,K-1) - DTL_RML(I) / DTRDZ(I,K-1)
            FQW(I,K) = AQ_AM(I,K) - AQ_AM(I,K-1)
     &                    + FQW(I,K-1) - DQW_RML(I) / DTRDZ(I,K-1)
          ELSE
            FTL(I,K) = AT_ATQ(I,K)
            FQW(I,K) = AQ_AM(I,K)
          ENDIF
        ENDDO ! p_points
      ENDDO ! bl_levels

      IF (LTIMER) THEN
        CALL TIMER('IMCALTQ ',4)
      ENDIF

      RETURN
      END
