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
!!!  SUBROUTINE EXCF_NL ----------------------------------------------
!!!
!!!  Purpose: To calculate non-local exchange coefficients for
!!!           boundary layer subroutine KMKH.
!!!
!!!  Suitable for single column use (via *IF definition IBM).
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!  4.4    Feb 1997  Written by R N B Smith
!!!
!!!  Programming standard:
!!!
!!!  System component covered: Part of P243.
!!!
!!!  Project task:
!!!
!!!  Documentation: UMDP No.24
!!!
!!!---------------------------------------------------------------------
!
!!   Arguments :-
      SUBROUTINE EXCF_NL (
     & P_FIELD,P1,P_POINTS,BL_LEVELS,
     & NTML,CF,
     & RDZ,ZH,Z_UV,Z_TQ,RHO_UV,RHO_TQ,
     & Z0M,V_S,FB_SURF,DB_TOP,
     & BFLUX_SURF,BFLUX_SAT_SURF,ZETA_S,BT_TOP,BTT_TOP,
     & DF_TOP_OVER_CP,ZETA_R,ALPHA_R,BTC_TOP,                 
     & DB_TOP_CLD,CHI_S_TOP,ZC,
     & RHOKM,RHOKH,RHOKM_TOP,RHOKH_TOP,
     & ZHSC,DSCDEPTH,NTDSC,DB_DSCT,SVL,
     & BT_DSCT,BTT_DSCT,
     & DF_DSCT_OVER_CP,ZETA_R_DSC,ALPHA_R_DSC,BTC_DSCT,
     & DB_DSCT_CLD,CHI_S_DSCT,ZC_DSC,COUPLED,
     & LTIMER
     &)
!
      IMPLICIT NONE
!
      LOGICAL LTIMER
!
      INTEGER
     & P_FIELD     ! IN No. of P-grid points in whole field.
     &,P1          ! IN First P-grid point to be processed.
     &,P_POINTS    ! IN No. of P-grid points to be processed.
     &,BL_LEVELS   ! IN maximum number of boundary layer levels
!
      LOGICAL
     & COUPLED(P_FIELD)          ! IN  Flag to indicate Sc layer weakly
!                                !     (de)coupled to surface.
!
      INTEGER
     & NTML(P_FIELD)             ! IN  Number of turbulently mixed
!                                !     layers.
     &,NTDSC(P_FIELD)            ! IN  Top level of any decoupled
!                                !     turbulently mixed Sc layer.
!
      REAL
     & RDZ(P_FIELD,BL_LEVELS)    ! IN Reciprocal of distance between
!                                !    T,q-levels (m^-1). 1/RDZ(,K) is
!                                !    the vertical distance from level
!                                !    K-1 to level K, except that for
!                                !    K=1 it is just the height of the
!                                !    lowest atmospheric level.
     &,ZH(P_FIELD)               ! IN Boundary layer height (m).
     &,Z_UV(P_FIELD,BL_LEVELS)   ! IN For a vertically staggered grid
!                                !    with a u,v-level first above the
!                                !    surface, Z_UV(*,K) is the height
!                                !    of the k-th u,v-level (half level
!                                !    k-1/2) above the surface;
!                                !    for an unstaggered grid the
!                                !    heights of the half-levels
!                                !    0.5 to BL_LEVELS-0.5 should be
!                                !    input to elements 1 to BL_LEVELS.
!                                !    (1st value not used in either
!                                !     case.)
     &,Z_TQ(P_FIELD,BL_LEVELS)   ! IN For a vertically staggered grid
!                                !    with a u,v-level first above the
!                                !    surface, Z_TQ(*,K) is the height
!                                !    of the k-th T,q-level (full level
!                                !    k) above the surface;
!                                !    for an unstaggered grid the
!                                !    heights of the half levels
!                                !    1.5 to BL_LEVELS+0.5 should be
!                                !    input to elements 1 to BL_LEVELS.
!                                !    (value for BL_LEVELS not used
!                                !    in either case.)
!
     &,RHO_UV(P_FIELD,BL_LEVELS) ! IN For a vertically staggered grid
!                                !    with a u,v-level first above the
!                                !    surface, RHO_UV(*,K) is the
!                                !    density at the k-th u,v-level
!                                !    above the surface;
!                                !    for an unstaggered grid the
!                                !    densities at the layer interfaces
!                                !    (half-levels) 0.5 to BL_LEVELS-0.5
!                                !    should be input to elements 1 to
!                                !    BL_LEVELS.
!                                !    (1st value not used in either
!                                !    case.)
     &,RHO_TQ(P_FIELD,BL_LEVELS) ! IN For a vertically staggered grid
!                                !    with a u,v-level first above the
!                                !    surface, RHO_TQ(*,K) is the
!                                !    density of the k-th T,q-level
!                                !    above the surface;
!                                !    for an unstaggered grid the
!                                !    densities at the layer interfaces
!                                !    (half-levels) 1.5 to BL_LEVELS+0.5
!                                !    should be input to elements 1 to
!                                !    BL_LEVELS.
!                                !    (value for BL_LEVELS not used
!                                !    in either case.)
!
     &,Z0M(P_FIELD)              ! IN Roughness length for momentum (m).
     &,V_S(P_FIELD)              ! IN Surface friction velocity (m/s).
     &,FB_SURF(P_FIELD)          ! IN Buoyancy flux at the surface over
!                                !    density (m^2/s^3).
     &,BFLUX_SURF(P_FIELD)       ! IN Surface buoyancy flux (clear air
!                                !    term) (kg/m/s^3).   
     &,BFLUX_SAT_SURF(P_FIELD)   ! IN Surface buoyancy flux (cloudy air
!                                !    term) (kg/m/s^3).   
     &,DB_TOP(P_FIELD)           ! IN Buoyancy jump across top of b.l
!                                !    (m/s^2)
     &,DF_TOP_OVER_CP(P_FIELD)   ! IN Radiative flux change at cloud top
!                                !    divided by c_P (K.kg/m^2/s).
     &,BT_TOP(P_FIELD)           ! IN Buoyancy parameter at the top of
!                                !    the b.l. (m/s^2/K).
     &,BTT_TOP(P_FIELD)          ! IN In-cloud buoyancy parameter at
!                                !    the top of the b.l. (m/s^2/K).
     &,BTC_TOP(P_FIELD)          ! IN Cloud fraction weighted buoyancy
!                                !    parameter at the top of the b.l.
     &,DB_TOP_CLD(P_FIELD)       ! IN In-cloud buoyancy jump at the
!                                !    top of the b.l. (m/s^2).
     &,CHI_S_TOP(P_FIELD)        ! IN Mixing fraction of just saturated
!                                !    mixture at top of the b.l.
     &,ZETA_S(P_FIELD)           ! IN Non-cloudy fraction of mixing
!                                !    layer for surface forced
!                                !    entrainment term.
     &,ZETA_R(P_FIELD)           ! IN Non-cloudy fraction of mixing
!                                !    layer for cloud top radiative
!                                !    cooling entrainment term.
     &,ALPHA_R(P_FIELD)          ! IN Fraction of the cloud top
!                                !    radiative cooling assumed to
!                                !    act above the minimum turbulent
!                                !    flux level.
     &,ZC(P_FIELD)               ! IN Cloud depth (not cloud fraction
!                                !    weighted) (m).
!
     &,CF(P_FIELD,BL_LEVELS)     ! IN Cloud fraction.
      REAL
     & ZHSC(P_FIELD)             ! IN Cloud layer height (m).
     &,DSCDEPTH(P_FIELD)         ! IN Decoupled cloud layer depth (m).
     &,DB_DSCT(P_FIELD)          ! IN Buoyancy parameter at the top of
!                                !    the decoupled Sc layer (m/s^2/K).
     &,SVL(P_FIELD,BL_LEVELS)    ! IN Liquid/frozen water virtual
!                                !    static energy over CP.
     &,DF_DSCT_OVER_CP(P_FIELD)  ! IN Radiative flux change at DSC top
!                                !    divided by c_P (K.kg/m^2/s).
     &,BT_DSCT(P_FIELD)          ! IN Buoyancy parameter at the top of
!                                !    the decoupled Sc  (m/s^2/K).
     &,BTT_DSCT(P_FIELD)         ! IN In-cloud buoyancy parameter at
!                                !    the top of the decoupled Sc
!                                !    (m/s^2/K).
     &,BTC_DSCT(P_FIELD)         ! IN Cloud fraction weighted buoyancy
!                                !    parameter at the top of the 
!                                !    decoupled Sc layer.
     &,DB_DSCT_CLD(P_FIELD)      ! IN In-cloud buoyancy jump at the
!                                !    top of the decoupled SC (m/s^2).
     &,CHI_S_DSCT(P_FIELD)       ! IN Mixing fraction of just saturated
!                                !    mixture at the top of the
!                                !    decoupled Sc layer.
     &,ZETA_R_DSC(P_FIELD)       ! IN Non-cloudy fraction of decoupled 
!                                !    Sc layer for cloud top radiative
!                                !    cooling entrainment term.
     &,ALPHA_R_DSC(P_FIELD)      ! IN Fraction of the cloud top
!                                !    radiative cooling assumed to
!                                !    act above the minimum turbulent
!                                !    flux level for the decoupled Sc.
     &,ZC_DSC(P_FIELD)           ! IN Cloud depth (not cloud fraction
!                                !    weighted) of the decoupled Sc (m).
      REAL                                                              
     & RHOKM(P_FIELD,2:BL_LEVELS)! OUT Layer k-1 - to - layer k
!                                !     turbulent mixing coefficient
!                                !     for momentum (kg/m/s).
     &,RHOKH(P_FIELD,2:BL_LEVELS)! OUT Layer k-1 - to - layer k
!                                !     turbulent mixing coefficient
!                                !     for heat and moisture (kg/m/s).
     &,RHOKM_TOP(P_FIELD,2:BL_LEVELS)
!                                ! OUT Non-local turbulent mixing
!                                !     coefficient for top-down mixing
!                                !     of momentum.
     &,RHOKH_TOP(P_FIELD,2:BL_LEVELS)
!                                ! OUT Non-local turbulent mixing
!                                !     coefficient for top-down mixing
!                                !     of heat and moisture.
!*
!*L---------------------------------------------------------------------
      EXTERNAL TIMER
!*
!*L---------------------------------------------------------------------
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

C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------

      REAL A_ENT_1,A_ENT_2,A_ENT_SHR,C_T
      PARAMETER (
     & A_ENT_1=0.23             ! Entrainment parameter.
     &,A_ENT_2=0.056            ! Entrainment parameter.
     &,A_ENT_SHR=5.0            ! Entrainment parameter.
     &,C_T=1.0                  ! Parameter in Zilitinkevich term.
     &)
!*
!
!  Define local storage.
!
!  (a) Workspace.
!
      REAL
     & W_M_TOP(P_FIELD)         ! Turbulent velocity scale for momentum
!                               ! evaluated at the top of the b.l.
     &,W_H_TOP(P_FIELD)         ! Turbulent velocity scale for scalars
!                               ! evaluated at the top of the b.l.
     &,V_TOP(P_FIELD)           ! Velocity scale for top-down mixing.
     &,PRANDTL_TOP(P_FIELD)     ! Turbulent Prandtl number
!                               ! evaluated at the top of the b.l.
     &,KH_TOP_FACTOR(P_FIELD)   ! Factor to ensure K_H profile is
!                               ! continuous at z_h.
     &,KM_TOP_FACTOR(P_FIELD)   ! Factor to ensure K_M profile is
!                               ! continuous at z_h.
     &,KH_SCT_FACTOR(P_FIELD)   ! Factor to ensure top-down K_H profile
!                               ! is continuous at z_h.
     &,KM_SCT_FACTOR(P_FIELD)   ! Factor to ensure top-down K_M profile
!                               ! is continuous at z_h.
     &,KH_DSCT_FACTOR(P_FIELD)  ! Factor to ensure K_H profile is 
!                               ! continuous at ZHSC.
     &,KM_DSCT_FACTOR(P_FIELD)  ! Factor to ensure K_M profile is 
!                               ! continuous at ZHSC.
     &,V_TOP_DSC(P_FIELD)       ! Velocity scale for top-down convection
!                               ! in decoupled Sc layer.
     &,SVL_PLUME(P_FIELD)       ! s_VL of descending plume to find
!                               ! decoupled Sc layer base.
     &,SCDEPTH(P_FIELD)         ! Depth of top-driven mixing in
!                               ! well-mixed Sc layer.
!  (b) Scalars.
!
      REAL
     & PRANDTL    ! Turbulent Prandtl number.
     &,ZK_UV      ! Height above surface of u,v-level.
     &,ZK_TQ      ! Height above surface of T,q-level.
     &,ZH_M       ! Height taken as top of the turbulent mixing layer
!                 ! for momentum.
     &,W_S_CUBED_UV ! Cube of free-convective velocity scale (u,v-level)
     &,W_S_CUBED_TQ ! Cube of free-convective velocity scale (T,q-level)
     &,W_M_UV     ! Turbulent velocity scale for momentum (u,v-level).
     &,W_M_TQ     ! Turbulent velocity scale for momentum (T,q-level).
     &,W_H_UV     ! Turbulent velocity scale for scalars (u,v-level).
     &,SF_TERM    ! Surface buoyancy flux term for entrainment paramn.
     &,SF_SHEAR_TERM
!                 ! Surface shear term for the entrainment paramn.
     &,IR_TERM    ! Indirect radiative term for entrainment paramn.
     &,DR_TERM    ! Direct radiative term for entrainment paramn.
     &,EVAP_TERM  ! Evaporative term in entrainment parametrization.
     &,ZIL_CORR   ! Zilitinkevich correction term in entrn. paramn.
     &,ZETA_S_FAC ! Factor involving ZETA_S.   
     &,ZETA_R_SQ  ! ZETA_R squared.
     &,ZR         ! Ratio ZC/ZH.
     &,Z_PR       ! Height above surface layer.
     &,ZH_PR      ! Height of cloud layer top above surface layer.
     &,ZCML_BASE  ! Height of base of cloud mixed layer.
     &,RHOKH_ENT  ! entrainment eddy viscosity.
     &,FRAC_TOP   ! Fraction of turbulent mixing driven from the top.
     &,FACTOR     ! Temporary scalar.
     &,ENT_FACTOR ! Factor to weight entrainment by CF.
     &,DTDZ_PAR   ! Eddy turnover timescale.
!
      INTEGER
     & I       ! Loop counter (horizontal field index).
     &,K       ! Loop counter (vertical level index).
!
      LOGICAL
     & SCBASE(P_FIELD)     ! Flag to signal base of convective mixed
!                          ! layer reached.
!
      IF (LTIMER) THEN
        CALL TIMER('EXCF_NL  ',3)
      ENDIF
!
!-----------------------------------------------------------------------
!! 0.  Calculate top-of-b.l. velocity scales and Prandtl number.        
!-----------------------------------------------------------------------
!
      DO I = P1,P1-1+P_POINTS

        V_TOP(I) = 0.0
        V_TOP_DSC(I) = 0.0
        SCBASE(I)=.FALSE.
        SVL_PLUME(I) = 0.0

        IF (FB_SURF(I) .GE. 0.0) THEN
!
!         By definition the top of the b.l. is in the 'outer layer' so
!         the free-convective velocity scale cubed is
!
          IF (COUPLED(I)) THEN
            W_S_CUBED_UV = 0.25 * ZHSC(I) * FB_SURF(I)
          ELSE
          W_S_CUBED_UV = 0.25 * ZH(I) * FB_SURF(I)
          ENDIF
!
!         Turbulent velocity scale for momentum
!
          W_M_TOP(I) = (V_S(I)*V_S(I)*V_S(I) + W_S_CUBED_UV)**(1.0/3.0)
!
!         Turbulent Prandtl number and velocity scale for scalars
!
          PRANDTL_TOP(I) = 0.75 * ( V_S(I)*V_S(I)*V_S(I)*V_S(I) +
     &                        (4.0/25.0)*W_S_CUBED_UV*W_M_TOP(I) ) /
     &                            ( V_S(I)*V_S(I)*V_S(I)*V_S(I) +
     &                        (8.0/25.0)*W_S_CUBED_UV*W_M_TOP(I) )
          W_H_TOP(I) = W_M_TOP(I) / PRANDTL_TOP(I)
        ELSE
          W_M_TOP(I) = V_S(I)
          PRANDTL_TOP(I) = 0.75
          W_H_TOP(I) = W_M_TOP(I) / PRANDTL_TOP(I)
        ENDIF
      ENDDO
!
!-----------------------------------------------------------------------
!! 1.  Loop round levels; calculate the top-of-b.l. entrainment
!!     mixing coefficients.
!-----------------------------------------------------------------------
!
      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
!
!-----------------------------------------------------------------------
!! 1.2 Calculate top-of-b.l. entrainment mixing coefficients
!!     and store b.l. top quantities for later use.
!-----------------------------------------------------------------------
!!         First the top of the surface mixed layer (if not coupled)
!-----------------------------------------------------------------------
          IF ( (K .EQ. NTML(I)+1) .AND. (DB_TOP(I) .GT. 0.0) 
     &          .AND. .NOT.COUPLED(I) ) THEN
!           !-----------------------------------------------------------
!           ! Calculate the surface buoyancy flux term
!           !-----------------------------------------------------------
            ZETA_S_FAC = (1.0 - ZETA_S(I)) * (1.0 - ZETA_S(I))   
            SF_TERM = A_ENT_1 * MAX ( 0.0 ,   
     &                          ( (1.0 - ZETA_S_FAC) * BFLUX_SURF(I)
     &                             + ZETA_S_FAC * BFLUX_SAT_SURF(I) ) )
!           !-----------------------------------------------------------
!           ! Calculate the surface shear term
!           !-----------------------------------------------------------
            SF_SHEAR_TERM = A_ENT_SHR * V_S(I) * V_S(I) * V_S(I)
     &                        * RHO_UV(I,K) / ZH(I)
!           !-----------------------------------------------------------
!           ! Calculate the indirect radiative term
!           !-----------------------------------------------------------
            ZETA_R_SQ = ZETA_R(I)*ZETA_R(I)
            IR_TERM = (BT_TOP(I)*ZETA_R_SQ + BTT_TOP(I)*(1.0-ZETA_R_SQ))
     &                 * A_ENT_1 * DF_TOP_OVER_CP(I)
!           !-----------------------------------------------------------
!           ! Calculate the evaporative term
!           !-----------------------------------------------------------
            ZR = SQRT( ZC(I) / ZH(I) )
            EVAP_TERM = A_ENT_2 * RHO_UV(I,K)
     &                  * CHI_S_TOP(I) * CHI_S_TOP(I)
     &                  * ZR * ZR * ZR * DB_TOP_CLD(I)
     &                  * SQRT( ZH(I) * DB_TOP(I) )
!           !-----------------------------------------------------------
!           ! Calculate the direct radiative term
!           !-----------------------------------------------------------
            DR_TERM = BTC_TOP(I) * ALPHA_R(I) * DF_TOP_OVER_CP(I)
!           !-----------------------------------------------------------
!           ! Finally combine terms to calculate the entrainment
!           ! mixing coefficients
!           !-----------------------------------------------------------
            IF (CF(I,NTML(I)) .GE. 0.9) THEN
              ENT_FACTOR = 1.0
            ELSE
              ENT_FACTOR = EXP(-((0.90-CF(I,NTML(I)))**3.0)/0.075)
            ENDIF
!
            ZIL_CORR = C_T * ( (SF_TERM + SF_SHEAR_TERM + ENT_FACTOR *
     &                          (IR_TERM + EVAP_TERM) ) /
     &                   (RHO_UV(I,K) * SQRT(ZH(I))) )**(2.0/3.0)

            RHOKH_ENT = (SF_TERM + SF_SHEAR_TERM
     &         + ENT_FACTOR * (IR_TERM + EVAP_TERM
     &          + DR_TERM)) / ((DB_TOP(I) + ZIL_CORR) * RDZ(I,K) )

            V_TOP(I) = ( ENT_FACTOR * (IR_TERM + EVAP_TERM) * ZH(I) / 
     &                      (A_ENT_1*RHO_UV(I,K)) )**(1.0/3.0)
            FRAC_TOP = V_TOP(I) / ( V_TOP(I) + W_H_TOP(I) + 1.0E-14 )

            RHOKH(I,K) = RHOKH_ENT * ( 1.0 - FRAC_TOP )
            RHOKM(I,K) = PRANDTL_TOP(I) * RHOKH(I,K)
     &                   * RHO_TQ(I,K-1) / RHO_UV(I,K)

            RHOKH_TOP(I,K) = RHOKH_ENT * FRAC_TOP
            RHOKM_TOP(I,K) = 0.75 * RHOKH_TOP(I,K)
     &                       * RHO_TQ(I,K-1) / RHO_UV(I,K)
          ELSE
            RHOKH(I,K) = 0.0
            RHOKM(I,K) = 0.0
            RHOKH_TOP(I,K) = 0.0
            RHOKM_TOP(I,K) = 0.0
          ENDIF
!-----------------------------------------------------------------------
!!        Then the top of the decoupled Sc (if coupled use ZHSC
!!        length-scale).
!-----------------------------------------------------------------------
          IF ( (K .EQ. NTDSC(I)+1) .AND. (DB_DSCT(I) .GT. 0.0) ) THEN
             IF (COUPLED(I)) THEN
!              !--------------------------------------------------------
!              ! Calculate the surface buoyancy flux term
!              !--------------------------------------------------------
               ZETA_S_FAC = (1.0 - ZETA_S(I)) * (1.0 - ZETA_S(I))   
               SF_TERM = A_ENT_1 * MAX ( 0.0 ,   
     &                          ( (1.0 - ZETA_S_FAC) * BFLUX_SURF(I)
     &                             + ZETA_S_FAC * BFLUX_SAT_SURF(I) ) )
!              !--------------------------------------------------------
!              ! Calculate the surface shear term
!              !--------------------------------------------------------
               SF_SHEAR_TERM = A_ENT_SHR * V_S(I) * V_S(I) * V_S(I)
     &                           * RHO_UV(I,K) / DSCDEPTH(I)
             ELSE
               SF_TERM = 0.0
               SF_SHEAR_TERM = 0.0
             ENDIF
!           !-----------------------------------------------------------
!           ! Calculate the indirect radiative term
!           !-----------------------------------------------------------
            ZETA_R_SQ = ZETA_R_DSC(I)*ZETA_R_DSC(I)
            IR_TERM = (BT_DSCT(I)*ZETA_R_SQ+BTT_DSCT(I)*(1.0-ZETA_R_SQ))
     &                 * A_ENT_1 * DF_DSCT_OVER_CP(I)
!           !-----------------------------------------------------------
!           ! Calculate the evaporative term
!           !-----------------------------------------------------------
            ZR = SQRT( ZC_DSC(I) / DSCDEPTH(I) )
            EVAP_TERM = A_ENT_2*RHO_UV(I,K) *CHI_S_DSCT(I)*CHI_S_DSCT(I)
     &                  * ZR * ZR * ZR * DB_DSCT_CLD(I)
     &                  * SQRT( DSCDEPTH(I) * DB_DSCT(I) )
!           !-----------------------------------------------------------
!           ! Calculate the direct radiative term
!           !-----------------------------------------------------------
            DR_TERM = BTC_DSCT(I) * ALPHA_R_DSC(I) * DF_DSCT_OVER_CP(I)
!           !-----------------------------------------------------------
!           ! Finally combine terms to calculate the entrainment
!           ! mixing coefficients
!           !-----------------------------------------------------------

            IF (CF(I,NTDSC(I)).GE.0.9) THEN
              ENT_FACTOR = 1.0
            ELSE
              ENT_FACTOR = EXP(-((0.90-CF(I,NTDSC(I)))**3.0)/0.075)
            ENDIF

            V_TOP_DSC(I) =( ENT_FACTOR * (IR_TERM + EVAP_TERM) * 
     &                      DSCDEPTH(I) / 
     &                         (A_ENT_1*RHO_UV(I,K)) )**(1.0/3.0)
            ZIL_CORR = C_T * ( (SF_TERM + SF_SHEAR_TERM +
     &                 ENT_FACTOR * (IR_TERM + EVAP_TERM)) /
     &                 (RHO_UV(I,K) * SQRT(DSCDEPTH(I))) )**(2.0/3.0)
            RHOKH_TOP(I,K) = ( SF_TERM + SF_SHEAR_TERM
     &               + ENT_FACTOR * (IR_TERM + EVAP_TERM + DR_TERM) )
     &                / ((DB_DSCT(I) + ZIL_CORR) * RDZ(I,K) )
            RHOKM_TOP(I,K) = 0.75 * RHOKH_TOP(I,K) 
     &                   * RHO_TQ(I,K-1) / RHO_UV(I,K)
          ENDIF
        ENDDO                                                           
      ENDDO                                                             
!-----------------------------------------------------------------------
!! 1.3  Parcel descent to find base of decoupled Sc layer (if there is
!!      one) can now include the entrainment source and so give a more
!!      accurate value for DSCDEPTH.
!-----------------------------------------------------------------------
!       Take the current radiatively determined DSCDEPTH as a starting
!       point (although restricted to less than 1.2*ZC in case the 
!       presence of entrainment makes the parcel positively buoyant
!       below cloud base - assume the cloud is well mixed and allow
!       for an overshoot).  Assume that the cloudy air parcel (with 
!       s_VL taken from NTDSC assumed to be representative of the 
!       convectively mixed layer) mixes with entrained air and cools
!       radiatively on a timescale given by the eddy turnover time
!       ( = DSCDEPTH/V_TOP_DSC) and a length scale of DSCDEPTH so the
!       flux into the parcel is multiplied by DTDZ_PAR = 1.0/V_TOP_DSC.
!
!          NOTE (1): if V_TOP_DSC = 0, DSCDEPTH is unchanged from its
!                    value estimated in KMKH.
!          NOTE (2): the entrainment calculation could be repeated
!                    using this new value of DSCDEPTH but the dependence
!                    on DSCDEPTH is small so this isn't done.
!          NOTE (3): the radiative increment has already been added at
!                    this stage suggesting SVL_PARCEL may have double
!                    counted its radiative contribution (leading to less
!                    decoupling) but it makes little difference to the
!                    results.
!-----------------------------------------------------------------------
        DO I = P1,P1-1+P_POINTS
          IF ( V_TOP_DSC(I) .GT. 0.0 ) THEN
            K=NTDSC(I)
            IF ( NTDSC(I) .LE. 2 ) THEN
!
!             ! svl_plume not needed
!
              DSCDEPTH(I) = Z_UV(I,K+1)
              SCBASE(I) = .TRUE.
            ELSE
              DTDZ_PAR = 1.0 / V_TOP_DSC(I)
              SVL_PLUME(I) =            !  plume perturbation
     &          ( RDZ(I,K+1) * RHOKH_TOP(I,K+1) * (SVL(I,K+1)-SVL(I,K))
     &            - DF_DSCT_OVER_CP(I) ) * DTDZ_PAR / RHO_UV(I,K+1)
              IF ( Z_UV(I,K-1) .GE. ZHSC(I)-ZC_DSC(I) .AND.
     &             SVL(I,K-1) .LT. SVL(I,K) ) THEN
!
!               ! cloud layer is at least two layers thick and there is
!               ! an inversion near cloud top, take layer NTDSC-1
!               ! as representative
!
                SVL_PLUME(I) = SVL(I,K-1) + SVL_PLUME(I)
              ELSE
!
!               ! cloud layer less than two layers thick so take
!               ! layer NTDSC as representative
!
                SVL_PLUME(I) = SVL(I,K) + SVL_PLUME(I)
              ENDIF
            ENDIF
          ENDIF
        ENDDO
!-----------------------------------------------------------------------
!  Loop over levels to find where plume becomes positively buoyant
!-----------------------------------------------------------------------
        DO K=BL_LEVELS-1,1,-1
          DO I = P1,P1-1+P_POINTS
 
            IF ( V_TOP_DSC(I).GT.0.0 .AND. NTDSC(I).GT.2 ) THEN
              IF ( .NOT.SCBASE(I) .AND. 
!                   not yet found Sc base
     &              ( SVL_PLUME(I) .GT. SVL(I,K) 
!                     plume positively buoyant
     &                .OR. K .EQ. 1 ) )
!                          plume reached the surface
     &          THEN
                DSCDEPTH(I) = MAX( ZHSC(I)-Z_UV(I,K+1), 
     &                          MIN(1.2*ZC_DSC(I),DSCDEPTH(I)) )
!
!                 depth must be at least MIN( 1.2*ZC , 
!                                         radiatively-determined-depth )
!
                SCBASE(I)=.TRUE.
              ENDIF
            ENDIF
          ENDDO
        ENDDO
!-----------------------------------------------------------------------
!! 1.4  Identical plume descent to the above but to find depth of 
!       top-down mixing in well mixed Sc layer, SCDEPTH.
!-----------------------------------------------------------------------
        DO I = P1,P1-1+P_POINTS
          SCBASE(I) = .FALSE.
          SVL_PLUME(I) = 0.0
          SCDEPTH(I) = 0.0   ! default for V_TOP eq 0
        ENDDO
        DO I = P1,P1-1+P_POINTS
          IF ( V_TOP(I) .GT. 0.0 ) THEN
            K=NTML(I)
            IF ( NTML(I).LE.2 ) THEN
! 
!             ! svl_plume not needed
!
              SCDEPTH(I) = Z_UV(I,K+1)
              SCBASE(I) = .TRUE.
            ELSE
              DTDZ_PAR = 1.0 / V_TOP(I)
              SVL_PLUME(I) =            !  plume perturbation
     &                ( ( SVL(I,K+1)-SVL(I,K) ) * RDZ(I,K+1)
     &                  * ( RHOKH(I,K+1)+RHOKH_TOP(I,K+1) )
     &                 - DF_TOP_OVER_CP(I) ) * DTDZ_PAR / RHO_UV(I,K+1)
              IF ( Z_UV(I,K-1) .GE. ZH(I)-ZC(I) .AND.
     &             SVL(I,K-1) .LT. SVL(I,K) ) THEN
!
!               ! cloud layer is at least two layers thick and there is
!               ! an inversion near cloud top, take layer NTML-1 as
!               ! representative
!
                SVL_PLUME(I) = SVL(I,K-1) + SVL_PLUME(I)
              ELSE
!
!               ! cloud layer less than two layers thick so take layer
!               ! NTDSC as representative
!
                SVL_PLUME(I) = SVL(I,K) + SVL_PLUME(I)
              ENDIF
            ENDIF
          ENDIF
        ENDDO
!-----------------------------------------------------------------------
!  Loop over levels to find where plume becomes positively buoyant
!-----------------------------------------------------------------------
        DO K=BL_LEVELS-1,1,-1
          DO I = P1,P1-1+P_POINTS
 
            IF ( V_TOP(I).GT.0.0 .AND. NTML(I).GT.2 ) THEN
              IF ( .NOT.SCBASE(I) .AND.  
!                   not yet found scbase
     &             (SVL_PLUME(I) .GT. SVL(I,K) 
!                   plume positively buoyant
     &              .OR. K .EQ. 1 )  )
!                        reached the surface
     &          THEN
                SCDEPTH(I) = MAX( ZH(I)-Z_UV(I,K+1), 1.2*ZC(I) )
!                     !   depth must be at least 1.2*ZC 
                SCBASE(I)=.TRUE.
              ENDIF
            ENDIF
          ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!     Calculate factors required to ensure that the non-local turbulent
!     mixing coefficient profiles are continuous as the entrainment 
!     level is approached.
!-----------------------------------------------------------------------
!
      DO I = P1,P1-1+P_POINTS
        K=NTML(I)+1  
!       ! for cubic form of KH and KM:                                  
            KH_TOP_FACTOR(I) = MAX( 0.7 , 1.0 - SQRT( RHOKH(I,K) /
     &                 ( RHO_UV(I,K)*W_H_TOP(I)*VKMAN*ZH(I) ) ) )
            KM_TOP_FACTOR(I) = MAX( 0.7 , 1.0 - SQRT( RHOKM(I,K) /
     &                 ( RHO_TQ(I,K-1)*W_M_TOP(I)*VKMAN*ZH(I) ) ) )
!
!       ! for quadratic form of KH and KM:                              
!           KH_TOP_FACTOR(I) = MAX( 0.9 , 1.0 - RHOKH(I,K) /
!    &                   ( RHO_UV(I,K)*W_H_TOP(I)*VKMAN*ZH(I) ) )
!           KM_TOP_FACTOR(I) = MAX( 0.9 , 1.0 - RHOKM(I,K) /
!    &                   ( RHO_TQ(I,K-1)*W_M_TOP(I)*VKMAN*ZH(I) ) )
!
        Z_PR = MIN( 0.9*ZH(I) , SCDEPTH(I) )
        FACTOR = 0.85 * RHO_UV(I,K) * V_TOP(I) * VKMAN * Z_PR
        IF ( FACTOR .GT. 0.0) THEN
          KH_SCT_FACTOR(I) = 1.0 -
     &         ( (RHOKH_TOP(I,K)*RHOKH_TOP(I,K)) / (FACTOR*FACTOR) )
        ELSE
          KH_SCT_FACTOR(I) = 1.0
          ENDIF
        Z_PR = MIN( 0.9*Z_TQ(I,K-1) , SCDEPTH(I) )
        FACTOR = 0.85 * RHO_TQ(I,K-1) * V_TOP(I) * VKMAN * Z_PR * 0.75
        IF ( FACTOR .GT. 0.0) THEN
          KM_SCT_FACTOR(I) = 1.0 -
     &         ( (RHOKM_TOP(I,K)*RHOKM_TOP(I,K)) / (FACTOR*FACTOR) )
        ELSE
          KM_SCT_FACTOR(I) = 1.0
        ENDIF
!
!       !---------------------------------------------------------------
!       !  Set up factors to ensure K profile continuity at ZHSC;
!       !  no need to limit size of factor as precise shape of top-down 
!       !  mixing profile not important.
!       !---------------------------------------------------------------
!
       IF (NTDSC(I) .GT. 0) THEN
!       !-------------------------------------------------------------
!       ! Only calculate _DSCT_FACTORs when a decoupled stratocumulus 
!       ! layer exists, i.e. NTDSC > 0.
!       !-------------------------------------------------------------
        K=NTDSC(I)+1
        Z_PR = MAX(0.0 , MIN( ZHSC(I) - 0.1*ZH(I) , DSCDEPTH(I) ) )
        FACTOR = 0.85*RHO_UV(I,K)*V_TOP_DSC(I)*VKMAN*Z_PR
        IF ( FACTOR .GT. 0.0) THEN
          KH_DSCT_FACTOR(I) = 1.0 - 
     &          ( (RHOKH_TOP(I,K)*RHOKH_TOP(I,K)) / (FACTOR*FACTOR) )
        ELSE
          KH_DSCT_FACTOR(I) = 1.0
        ENDIF

        Z_PR = MAX(0.0, MIN( Z_TQ(I,K-1) - 0.1*Z_TQ(I,NTML(I)) ,
     &                                                DSCDEPTH(I) ) )
        FACTOR = 0.85*RHO_TQ(I,K-1)*V_TOP_DSC(I)*VKMAN*Z_PR*0.75
        IF ( FACTOR .GT. 0.0) THEN
          KM_DSCT_FACTOR(I) = 1.0 -
     &         ( (RHOKM_TOP(I,K)*RHOKM_TOP(I,K)) / (FACTOR*FACTOR) )
        ELSE
          KM_DSCT_FACTOR(I) = 1.0
        ENDIF
       ENDIF
        ENDDO
!
!-----------------------------------------------------------------------
!! 2.  Loop around levels again calculating height dependent turbulent
!!     transport coefficients within the mixing layer.
!-----------------------------------------------------------------------
!
      DO K=2,BL_LEVELS
        DO I = P1,P1-1+P_POINTS
!
!           Calculate the height of u,v-level above the surface
!
            ZK_UV = Z_UV(I,K) + Z0M(I)
!
!           Calculate the height of T,q-level above the surface
!
            ZK_TQ = Z_TQ(I,K-1) + Z0M(I)
!
!         !-------------------------------------------------------------
!         ! Calculate RHOK(H/M)_TOP, top-down turbulent mixing profiles 
!         ! for the surface mixed layer.
!         ! This is a variation on an up-side-down version of the cubic
!         ! surface-forced profiles below.  Implement between at least
!         ! the top of the `surface layer' (at Z=0.1*ZH) and ZH.
!         !-------------------------------------------------------------
!
          ZCML_BASE = MAX( 0.1*ZH(I) , ZH(I)-SCDEPTH(I) )
          IF ( ZK_UV .LT. ZH(I) .AND.
     &         ZK_UV .GT. ZCML_BASE ) THEN
            Z_PR = ZK_UV - ZCML_BASE
            ZH_PR = ZH(I) - ZCML_BASE
            RHOKH_TOP(I,K) = RHO_UV(I,K) * V_TOP(I) * 0.85 * VKMAN *
     &             SQRT( 1.0 - KH_SCT_FACTOR(I)*Z_PR/ZH_PR ) 
     &                                         * Z_PR * Z_PR / ZH_PR
          ENDIF
         IF (NTDSC(I) .GT. 0) THEN
!         !-------------------------------------------------------------
!         ! Only add contribution to top-down mixing coefficient
!         ! profiles for decoupled stratocumulus layers when
!         ! one exists, i.e. NTDSC > 0.
!         !-------------------------------------------------------------
          ZCML_BASE = MAX( 0.1*ZH(I) , ZHSC(I)-DSCDEPTH(I) )
          IF ( ZK_UV .LT. ZHSC(I) .AND.
     &            ZK_UV .GT. ZCML_BASE ) THEN
!           !-----------------------------------------------------------
!           ! Calculate RHOK(H/M)_TOP, top-down turbulent mixing 
!           ! profiles and add to any generated in the surface mixing
!           ! layer.
!           ! This is a variation on an up-side-down version of the
!           ! cubic surface-forced profiles above.  Implement between
!           ! at least the top of the `surface layer' (at Z=0.1*ZH) and
!           ! ZHSC.
!           !-----------------------------------------------------------
            Z_PR = ZK_UV - ZCML_BASE
            ZH_PR = ZHSC(I) - ZCML_BASE
            RHOKH_TOP(I,K) = RHOKH_TOP(I,K) + 
     &              RHO_UV(I,K)*V_TOP_DSC(I)*0.85*VKMAN*
     &              SQRT( 1.0 - KH_DSCT_FACTOR(I)*Z_PR/ZH_PR ) 
     &                                         * Z_PR * Z_PR / ZH_PR
          ENDIF
         ENDIF
!
          ZCML_BASE = MAX( 0.1*Z_TQ(I,NTML(I)) , 
     &                                  Z_TQ(I,NTML(I))-SCDEPTH(I) )
          IF ( ZK_TQ .LT. Z_TQ(I,NTML(I)) .AND.
     &         ZK_TQ .GT. ZCML_BASE ) THEN
            Z_PR = ZK_TQ - ZCML_BASE
            ZH_PR = Z_TQ(I,NTML(I)) - ZCML_BASE
            RHOKM_TOP(I,K) = 0.75 * RHO_TQ(I,K-1) * V_TOP(I) * 0.85 *
     &               VKMAN * SQRT( 1.0 - KM_SCT_FACTOR(I)*Z_PR/ZH_PR ) 
     &                                         * Z_PR * Z_PR / ZH_PR
          ENDIF
         IF (NTDSC(I) .GT. 0) THEN
!         !-------------------------------------------------------------
!         ! Only add contribution to top-down mixing coefficient
!         ! profiles for decoupled stratocumulus layers when
!         ! one exists, i.e. NTDSC > 0.
!         !-------------------------------------------------------------
          ZCML_BASE = MAX( 0.1*Z_TQ(I,NTML(I)) , 
     &                                  Z_TQ(I,NTDSC(I))-DSCDEPTH(I) )
          IF ( ZK_TQ .LT. Z_TQ(I,NTDSC(I)) .AND.
     &            ZK_TQ .GT. ZCML_BASE ) THEN
              Z_PR = ZK_TQ - ZCML_BASE
              ZH_PR = Z_TQ(I,NTDSC(I)) - ZCML_BASE
              RHOKM_TOP(I,K) = RHOKM_TOP(I,K) + 
     &           0.75*RHO_TQ(I,K-1)*V_TOP_DSC(I)*0.85*VKMAN*
     &           SQRT( 1.0 - KM_DSCT_FACTOR(I)*Z_PR/ZH_PR ) 
     &                                      * Z_PR * Z_PR / ZH_PR
          ENDIF
         ENDIF
!                                                                       
          IF (FB_SURF(I) .GE. 0.0) THEN
!           Calculate the free-convective scaling velocity at z(k)
!
            IF (ZK_UV .LE. 0.1*ZH(I)) THEN
!
!             Surface layer calculation
!
              W_S_CUBED_UV = 2.5 * ZK_UV * FB_SURF(I)
            ELSE
!
!             Outer layer calculation
!
              IF (COUPLED(I)) THEN  !  coupled and cloudy
                W_S_CUBED_UV = 0.25 * ZHSC(I) * FB_SURF(I)
              ELSE
              W_S_CUBED_UV = 0.25 * ZH(I) * FB_SURF(I)
            ENDIF
            ENDIF                                                       
            IF (ZK_TQ .LE. 0.1*Z_TQ(I,NTML(I))) THEN
!
!             Surface layer calculation
!
              W_S_CUBED_TQ = 2.5 * ZK_TQ * FB_SURF(I)
            ELSE
!
!             Outer layer calculation
!
              IF (COUPLED(I)) THEN  !  coupled and cloudy
                W_S_CUBED_TQ = 0.25 * ZHSC(I) * FB_SURF(I)
              ELSE
              W_S_CUBED_TQ = 0.25 * ZH(I) * FB_SURF(I)
              ENDIF
 
            ENDIF
!
!           Turbulent velocity scale for momentum
!
            W_M_UV = (V_S(I)*V_S(I)*V_S(I) + W_S_CUBED_UV)**(1.0/3.0)
!
            W_M_TQ = (V_S(I)*V_S(I)*V_S(I) + W_S_CUBED_TQ)**(1.0/3.0)
!
!           Turbulent Prandtl number and velocity scale for scalars
!
            PRANDTL = 0.75 * ( V_S(I)*V_S(I)*V_S(I)*V_S(I) +
     &                   (4.0/25.0)*W_S_CUBED_UV*W_M_UV ) /
     &                       ( V_S(I)*V_S(I)*V_S(I)*V_S(I) +
     &                   (8.0/25.0)*W_S_CUBED_UV*W_M_UV )
            W_H_UV = W_M_UV / PRANDTL
!
            IF ( K .LE. NTML(I) ) THEN
!             !---------------------------------------------------------
!             ! Calculate RHOKH(w_h,z/z_h)
!             !---------------------------------------------------------
!             N.B. ZH(I) = Z_UV(I,NTML(I)+1)
!
!               with cubic form
!
              RHOKH(I,K) = RHO_UV(I,K) * W_H_UV * VKMAN * ZK_UV *
     &                ( 1.0 - KH_TOP_FACTOR(I) * ( ZK_UV / ZH(I) ) ) *
     &                ( 1.0 - KH_TOP_FACTOR(I) * ( ZK_UV / ZH(I) ) )
!
!               with quadratic form
!
!             RHOKH(I,K) = RHO_UV(I,K) * W_H_UV * VKMAN * ZK_UV *
!    &                ( 1.0 - KH_TOP_FACTOR(I) * ( ZK_UV / ZH(I) ) )
!             !---------------------------------------------------------
!             ! Calculate RHOKM(w_m,z/z_h)
!             !---------------------------------------------------------
              ZH_M = Z_TQ(I,NTML(I))
!
!               with cubic form
!
              RHOKM(I,K) = RHO_TQ(I,K-1) * W_M_TQ * VKMAN * ZK_TQ *
     &                ( 1.0 - KM_TOP_FACTOR(I) * ( ZK_TQ / ZH_M ) ) *
     &                ( 1.0 - KM_TOP_FACTOR(I) * ( ZK_TQ / ZH_M ) )
!
!               with quadratic form
!
!             RHOKM(I,K) = RHO_TQ(I,K-1) * W_M_TQ * VKMAN * ZK_TQ *
!    &                ( 1.0 - KM_TOP_FACTOR(I) * ( ZK_TQ / ZH_M ) )

            ENDIF
          ENDIF
        ENDDO
      ENDDO
      IF (LTIMER) THEN
        CALL TIMER('EXCF_NL  ',4)
      ENDIF
      RETURN
      END
