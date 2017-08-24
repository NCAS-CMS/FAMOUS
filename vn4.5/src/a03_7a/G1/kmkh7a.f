C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1998, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!  SUBROUTINE KMKH---------------------------------------------------
!!!
!!!  Purpose: To calculate the turbulent mixing coefficients KM and KH
!!!
!!!  Suitable for single column use.
!!!
!!!  Model            Modification history
!!! version  date
!!!
!!!   4.5  13/05/98   New deck.  Richard Betts
!!!
!!!  Documentation: UMDP No.24
!!!
!!!---------------------------------------------------------------------

! Arguments :-

      SUBROUTINE KMKH (
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,TIMESTEP,P,CCA,BT,BQ,BF,CF,DZL,DTRDZ
     &,RDZ,U,V,FTL,FQW
     &,RHO,Z0M,ZLB,H_BLEND
     &,QW,QCF,RHOKM,RHO_KM,RHOKH,TL,ZH
     &,CCB,CCT,L_MOM
     &,NRML,L_BL_LSPICE
     &,LTIMER
     & )

      IMPLICIT NONE

      LOGICAL LTIMER             ! IN Flag for TIMER diagnostics

      LOGICAL
     & L_MOM       ! IN Switch for convective momentum transport.
     &,L_BL_LSPICE       ! IN
!                              TRUE  Use scientific treatment of mixed
!                                    phase precip scheme.
!                              FALSE Do not use mixed phase precip
!                                    considerations

      INTEGER
     & P_FIELD                ! IN No. of P-grid points in whole field
     &,P1                     ! IN First P-grid point to be processed.
     &,P_POINTS               ! IN No. of P-grid points to be
!                                  processed.
     &,BL_LEVELS              ! IN No. of atmospheric levels for
!                                  which boundary layer fluxes are
!                                  calculated.

      REAL
     & TIMESTEP               ! IN Timestep in seconds.
     &,BQ(P_FIELD,BL_LEVELS)  ! IN A buoyancy parameter (beta q tilde)
     &,BT(P_FIELD,BL_LEVELS)  ! IN A buoyancy parameter (beta T tilde)
     &,BF(P_FIELD,BL_LEVELS)  ! IN A buoyancy parameter (beta F tilde)
     &,P(P_FIELD,BL_LEVELS)   ! IN Pressure at pressure points.
     &,CCA(P_FIELD)           ! IN Convective Cloud Amount.
     &,CF(P_FIELD,BL_LEVELS)  ! IN Cloud fractions for boundary levs.
     &,DZL(P_FIELD,BL_LEVELS) ! IN Layer depths (m).  DZL(,K) is the
!                                  distance from layer boundary K-1/2
!                                  to layer boundary K+1/2.  For K=1
!                                  the lower boundary is the surface.
     &,DTRDZ(P_FIELD,BL_LEVELS)
!                             ! IN -g.dt/dp for model layers.

      REAL                    ! Split to keep continuation cards
     & RDZ(P_FIELD,BL_LEVELS) ! IN Reciprocal of distance between
!                                  hybrid levels (m-1).  1/RDZ(,K) is
!                                  the vertical distance from level
!                                  K-1 to level K, except that for
!                                  K=1 it is just the height of the
!                                  lowest atmospheric level.
     &,RHO(P_FIELD,BL_LEVELS) ! IN Density at theta_levels
     &,U(P_FIELD,BL_LEVELS)   ! IN Westerly wind component on P-grid
!                                  (m per s).
     &,V(P_FIELD,BL_LEVELS)   ! IN Southerly wind component on P-grid
!                                  (m per s).
     &,Z0M(P_FIELD)           ! IN Roughness length for momentum (m).
     &,ZLB(P_FIELD,BL_LEVELS) ! IN ZLB(,K) is height above surface of
!                                  lower boundary of layer K (m).
     &,H_BLEND(P_FIELD)       ! IN Blending height for effective
!                                  roughness scheme passed through
!                                  EX_COEF

      REAL                    ! Note: for all these INOUT arrays,
!                                     apart from ZH, level 1 is IN
!                                     (though not always used), and
!                                     the other levels are all OUT.
     & QW(P_FIELD,BL_LEVELS)  ! INOUT Total water content (kg per kg
!                                        air).
     &,QCF(P_FIELD,BL_LEVELS)  ! INOUT Ice water content (kg per kg
!                                        air).
     &,RHOKM(P_FIELD,BL_LEVELS)
!                             ! INOUT Layer K-1 - to - layer K
!                                     exchange coefficient for
!                                     momentum .
     &,TL(P_FIELD,BL_LEVELS)  ! INOUT Liquid/frozen water temperature
!                                     (K).
     &,ZH(P_FIELD)            ! INOUT Boundary layer height (m).

      REAL
     & RHO_KM(P_FIELD,2:BL_LEVELS)
!                             ! OUT RHO * KM before interpolation
!                                   to UV-grid.
     &,RHOKH(P_FIELD,BL_LEVELS)
!                             ! OUT Layer K-1 - to - layer K
!                                   exchange coefficient for FTL,
!                                   on P-grid.

      INTEGER
     & CCB(P_FIELD)           ! IN  Convective Cloud Base.
     &,CCT(P_FIELD)           ! IN  Convective Cloud Top.
     &,NRML(P_FIELD)          ! INOUT Number of model layers in the
!                                     unstable Rapidly Mixing Layer.

!!----------------------------------------------------------------------
!    External references :-
      EXTERNAL TIMER, EX_COEF

!!----------------------------------------------------------------------
!    Local and other symbolic constants :-
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

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_GAMMA------------------------------------
C GAMMA holds the implicit weighting coeff. for up to 30 B.L. levels.
C It is only required for the the number of B.L. levels actually used,
C so it does not need to be set up to 30 when less BL levels are used.
      REAL GAMMA(30)       ! Max of 30 Boundary Layer levels assumed.
C
      DATA GAMMA / 2 * 2.0 , 1.5 , 27 * 1.0 /
C*----------------------------------------------------------------------

      REAL ETAR,GRCP,LCRCP,LFRCP,LS,LSRCP
      PARAMETER (
     & ETAR=1.0/(1.0-EPSILON)   ! Used in buoyancy parameter BETAC.
     &,GRCP=G/CP                ! Used in DZTL, FTL calculations.
     &,LCRCP=LC/CP              ! Latent heat of condensation / CP.
     &,LFRCP=LF/CP              ! Latent heat of fusion / CP.
     &,LS=LC+LF                 ! Latent heat of sublimation.
     &,LSRCP=LS/CP              ! Latent heat of sublimation / CP.
     &)


!  Define local storage.

!  (a) Workspace. 6*BL_LEVELS-1 blocks of real workspace are required
!      plus 1 block of logical.


      LOGICAL
     & TOPBL(P_FIELD)            ! Flag set when top of boundary layer
!                                ! is reached.
      REAL
     & FQW(P_FIELD,BL_LEVELS)    ! "Explicit" flux of QW (i.e.
!                                  evaporation) from layer below,
!                                  on P-grid (kg per sq m per s).
     &,FTL(P_FIELD,BL_LEVELS)    ! "Explicit" flux of TL = H/CP
!                                  (sensible heat/CP) from layer
!                                   below, on P-grid.
     &,RI(P_FIELD,2:BL_LEVELS)   ! Richardson number for lower interface
!                                  of layer.
     &,RIM(P_FIELD,2:BL_LEVELS)  ! Modified Richardson number for lower
!                                  interface of layer.
     &,DTRDZ_RML (P_FIELD)       ! -g.dt/dp for the rapidly mixing layer
     &,DELTAP (P_FIELD,BL_LEVELS)! -g.dt/dp for the rapidly mixing layer
     &,DELTAP_RML (P_FIELD)      ! -g.dt/dp for the rapidly mixing layer

      INTEGER
     & NTML(P_FIELD)             ! No. of turbulently mixed model levels

!  (b) Scalars.

      REAL
     & DZB           ! Temporary in calculation of RI(,K).
     &,DZQW          ! Difference of QW between "current" and "lower"
!                       levels.
     &,DZTL          ! Liquid/ice static energy difference between
!                      adjacent levels.
     &,DZQI          ! Difference of QCF between "current" and "lower"
!                       levels
     &,DZU           ! Westerly wind shear between adjacent levels.
     &,DZV           ! Southerly wind shear between adjacent levels.
     &,DVMOD2        ! Square of modulus of wind shear between adjacent
!                      levels
     &,DTL_RML_P     ! Explicit TL increment for the rapidly mixing
!                      layer
     &,DQW_RML_P     ! Explicit QW increment for the rapidly mixing
!                      layer
     &,DTL_RMLP1_P   ! Explicit TL increment for the model layer above
!                      the rapidly mixing layer.
     &,DQW_RMLP1_P   ! Explicit QW increment for the model layer above
!                      the rapidly mixing layer.
     &,RIT           ! Temporary calculation of the modified Richardson
!                      number at the interface between the rapidly
!                      mixing layer and the model layer above it.

      INTEGER
     & I             ! Loop counter (horizontal field index).
     &,K             ! Loop counter (vertical level index).
     &,MBL           ! Maximum number of model layers allowed in the
!                      rapidly mixing layer; set to BL_LEVELS-1.
     &,NRMLP1        ! NRML + 1
     &,NRMLP2        ! NRML + 2
     &,IT_COUNTER    ! Iteration loop counter.

!-----------------------------------------------------------------------
!!  0.  Check that the scalars input to define the grid are consistent.
!!      See comments to routine SF_EXCH for details.
!-----------------------------------------------------------------------

      IF (LTIMER) THEN
        CALL TIMER('KMKH    ',3)
      ENDIF

!  Set MBL, "maximum number of boundary levels" for the purposes of
!  boundary layer height calculation.

      MBL = BL_LEVELS - 1

!-----------------------------------------------------------------------
!! 1 Initialise flag for having reached top of boundary layer
!!   and also the number of turbulently mixed layers
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1
        TOPBL(I) = .FALSE.
        NTML(I) = 1
      ENDDO

      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1

!-----------------------------------------------------------------------
!! 2.1 Calculate wind shear and other "jumps" between "current" and
!!     previous (lower) level.
!-----------------------------------------------------------------------

          DZU = U(I,K) - U(I,K-1)
          DZV = V(I,K) - V(I,K-1)
          DVMOD2 = MAX ( 1.0E-12 , DZU*DZU + DZV*DZV ) ! Used in P243.C1
          DZTL = TL(I,K) - TL(I,K-1) + GRCP/RDZ(I,K)   ! Used in P243.C2
          DZQW = QW(I,K) - QW(I,K-1)                   ! Used in P243.C2
          DZQI = QCF(I,K) - QCF(I,K-1)

          IF (L_BL_LSPICE) THEN
            DZB =  BT(I,K-1)*DZTL + BQ(I,K-1)*DZQW - BF(I,K-1)*DZQI
          ELSE
            DZB =  BT(I,K-1)*DZTL + BQ(I,K-1)*DZQW
          ENDIF

!-----------------------------------------------------------------------
!! 2.2 Calculate Richardson number Ri at the same interface.
!-----------------------------------------------------------------------

          RI(I,K) = G * DZB / ( RDZ(I,K) * DVMOD2 )

!-----------------------------------------------------------------------
!! 2.3 If either a stable layer (Ri > 1) or the maximum boundary layer
!!     height has been reached, set boundary layer height (ZH) to
!!     the height of the lower boundary of the current layer and set
!!     the number of rapidly mixing layers if the surface layer is
!!     unstable (as determined in subroutine SF_EXCH).
!-----------------------------------------------------------------------

          IF ( .NOT.TOPBL(I) .AND. (RI(I,K).GT.1.0 .OR. K.GT.MBL) ) THEN
            TOPBL(I) = .TRUE.
            ZH(I) = ZLB(I,K)
            NTML(I) = K-1
            IF ( NRML(I) .GT. 0 ) NRML(I) = K-1
          ENDIF
        ENDDO ! BL_LEVELS
      ENDDO ! P_POINTS

!-----------------------------------------------------------------------
!! 3.  Subroutine EX_COEF.
!-----------------------------------------------------------------------

      CALL EX_COEF (
     & P_FIELD,P1,P_POINTS,BL_LEVELS,
     & CCB,CCT,NTML,L_MOM,
     & CCA,DZL,RDZ,RI,U,V,
     & RHO,ZH,ZLB,Z0M,H_BLEND,
     & RHOKM,RHOKH,LTIMER
     &)

!-----------------------------------------------------------------------
!! 4.  Calculate "explicit" fluxes of TL and QW.
!-----------------------------------------------------------------------

      DO K=2,BL_LEVELS
!-----------------------------------------------------------------------
!! 4.1 "Explicit" fluxes of TL and QW, on P-grid.
!-----------------------------------------------------------------------
        DO I=P1,P1+P_POINTS-1
          FTL(I,K) = -RHOKH(I,K) *
     &     ( ( ( TL(I,K) - TL(I,K-1) ) * RDZ(I,K) ) + GRCP )
!          1 2 3                     3            2        1

          FQW(I,K)=-RHOKH(I,K) * ( QW(I,K) - QW(I,K-1) )*RDZ(I,K)

        ENDDO ! P_POINTS
      ENDDO ! BL_LEVELS

!-----------------------------------------------------------------------
!! 4.2 Use explicit fluxes to calculate a modified Richardson number
!!     at the top of the rapidly mixing layer (if it exists); if this
!!     indicates that the r.m.l. can deepen due to heat and/or moisture
!!     input from the surface then increase NRML(I).
!-----------------------------------------------------------------------

! Initialise height difference dz, for the rapidly mixing boundary layer

      DO I = P1,P1+P_POINTS-1
          DELTAP_RML(I) = 0.0
      ENDDO ! loop over P_POINTS.

! Calculate pressure differences, dp, and -g.dt/dp for model layers
! and dp for the rapidly mixing layer.

      DO K = 1,BL_LEVELS
        DO I = P1,P1+P_POINTS-1
          DELTAP(I,K) = -G * TIMESTEP/DTRDZ(I,K)
          IF (K .LE. NRML(I) )
     &        DELTAP_RML(I) = DELTAP_RML(I) + DELTAP(I,K)
          IF (K .GE. 2) RIM(I,K) = RI(I,K)
        ENDDO ! Loop over p-points
      ENDDO ! Loop over bl-levels

!-----------------------------------------------------------------------
!! 4.2.1 Iterate BL_LEVELS-2 times; this allows an initial r.m.l. of
!!       1 model layer to deepen to BL_LEVELS-1 model layers a layer at
!!       a time in each iteration.  Some of these iterations may be
!!       redundant if in fact the r.m.l. cannot deepen.
!-----------------------------------------------------------------------
      DO IT_COUNTER = 1,BL_LEVELS-2

        DO I = P1,P1+P_POINTS-1
!         !-------------------------------------------------------------
!         ! Only check to see if the r.m.l. can deepen if it exists
!         ! in the first place and has less than its maximum depth.
!         !-------------------------------------------------------------
          IF ((NRML(I) .GE. 1) .AND. (NRML(I) .LE. BL_LEVELS-2)) THEN
            NRMLP1 = NRML(I) + 1
            NRMLP2 = NRML(I) + 2
            DTRDZ_RML(I) =-G * TIMESTEP / DELTAP_RML(I)
!           !-----------------------------------------------------------
!           ! "Explicit" rapidly mixing layer increments to TL and QW.
!           !-----------------------------------------------------------
            DTL_RML_P = -DTRDZ_RML(I) * ( FTL(I,NRMLP1) - FTL(I,1) )
            DQW_RML_P = -DTRDZ_RML(I) * ( FQW(I,NRMLP1) - FQW(I,1) )
!           !-----------------------------------------------------------
!           ! "Explicit" increments to TL and QW in the model layer
!           ! above the rapidly mixing layer.
!           !-----------------------------------------------------------
            DTL_RMLP1_P = -DTRDZ(I,NRMLP1) *
     &                           ( FTL(I,NRMLP2) - FTL(I,NRMLP1) )
            DQW_RMLP1_P = -DTRDZ(I,NRMLP1) *
     &                           ( FQW(I,NRMLP2) - FQW(I,NRMLP1) )
            DZU = U(I,NRMLP1) - U(I,NRML(I))
            DZV = V(I,NRMLP1) - V(I,NRML(I))
            DVMOD2 = MAX (1.0E-12 , DZU*DZU + DZV*DZV)
!           !-----------------------------------------------------------
!           ! Calculate a modified Richardson number for the interface
!           ! between the rapidly mixing layer and the model layer above
!           !-----------------------------------------------------------
              RIT = RI(I,NRMLP1) + ( G/(RDZ(I,NRMLP1) * DVMOD2) ) *
     &              ( BT(I,NRML(I)) * ( DTL_RMLP1_P - DTL_RML_P)
     &              + BQ(I,NRML(I)) * ( DQW_RMLP1_P - DQW_RML_P ) )
            IF ( RIT .LE. 1.0 ) THEN
!             !---------------------------------------------------------
!             ! Deepen the rapidly mixing layer by one model layer
!             !---------------------------------------------------------
              NRML(I) = NRMLP1
              DELTAP_RML(I) = DELTAP_RML(I) + DELTAP(I,NRMLP1)
              ZH(I) = ZLB(I,NRMLP2)
              IF (RIT .LT. RI(I,NRMLP1)) RIM(I,NRMLP1) = RIT
            ENDIF
          ENDIF ! NRML(I) .GE. 1 and .LE. BL_LEVELS-2
        ENDDO ! Loop over p-points
      ENDDO! Loop over iterations


!-----------------------------------------------------------------------
!! 4.2.2 Adjust the Richardson number at the top of the rapidly mixing
!!       layer- the 'DZ' in the expression (P243.C1) is set to 100.0 m
!!       rather than the model level separation (if the latter is
!!       greater than 100 m).  This is a simple way of adjusting for
!!       inaccuracies in calculating gradients at inversions when the
!!       model's vertical resolution is coarse.
!-----------------------------------------------------------------------

      DO K = 2,BL_LEVELS
        DO I = P1,P1+P_POINTS-1
          IF ( (K-1 .EQ. NRML(I)) .AND. (100.0*RDZ(I,K) .LT. 1.0) )
     &       RIM(I,K) = 100.0*RDZ(I,K)*RIM(I,K)
        ENDDO ! Loop over p-points
      ENDDO ! Loop over bl-levels

!-----------------------------------------------------------------------
!! 5.  Subroutine EX_COEF using adjusted Richardson Number, RIM
!-----------------------------------------------------------------------

      CALL EX_COEF (
     & P_FIELD,P1,P_POINTS,BL_LEVELS,
     & CCB,CCT,NTML,L_MOM,
     & CCA,DZL,RDZ,RIM,U,V,
     & RHO,ZH,ZLB,Z0M,H_BLEND,
     & RHOKM,RHOKH,LTIMER
     &)

!-----------------------------------------------------------------------
!! 5.1 store RHO_KM on P-grid for otput.
!-----------------------------------------------------------------------
      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          RHO_KM(I,K) = RHOKM(I,K)
        ENDDO ! P_POINTS
      ENDDO ! BL_LEVELS

      IF (LTIMER) THEN
        CALL TIMER('KMKH    ',4)
      ENDIF

      RETURN
      END
