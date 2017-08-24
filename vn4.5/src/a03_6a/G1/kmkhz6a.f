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
!  SUBROUTINE KMKHZ  --------------------------------------------------
!!!
!!!  Purpose: To calculate the turbulent mixing coefficients KM and KH
!!!
!!!
!!!  Model            Modification history
!!! version  date
!!!
!!!   4.4  05/09/97   New deck   R.N.B.Smith
!!!  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!!!
!!!  Programming standard:
!!!
!!!  System component covered: Part of P243.
!!!
!!!  Documentation: UMDP No.24
!!!
!!!---------------------------------------------------------------------

! Arguments :-

      SUBROUTINE KMKHZ (
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,P,P_HALF,T,Q,QCL,QCF,BT,BQ,CF,DZL
     &,RDZ,DELTAP,FTL,FQW
     &,Z0M,Z_FULL,Z_HALF,Z_UV,Z_TQ,V_S,FB_SURF
     &,QW,RHOKM,DB_SVL,RHOKH,TL,ZH,TV1_SD,T1_SD,Q1_SD,NTML
     &,GRAD_T_ADJ,GRAD_Q_ADJ
     &,BTM,BQM,DQSDT,BTM_CLD,BQM_CLD,A_QSM,A_DQSDTM,RHO_TQ,RHO_UV
     &,RAD_HR,RADHR_DIM1,CUMULUS,Z_LCL,RHOKM_TOP,RHOKH_TOP,ZHT  
     &,BL_TYPE_1,BL_TYPE_2,BL_TYPE_3,BL_TYPE_4,BL_TYPE_5,BL_TYPE_6
     &,UNSTABLE,NTDSC,DSC
     &,LTIMER
     & )

      IMPLICIT NONE

      LOGICAL LTIMER             ! IN Flag for TIMER diagnostics

      INTEGER
     & P_FIELD                ! IN No. of P-grid points in whole field
     &,P1                     ! IN First P-grid point to be processed.
     &,P_POINTS               ! IN No. of P-grid points to be
!                                  processed.
!                                  ( = P_POINTS = 1 for SCM.)
     &,BL_LEVELS              ! IN No. of atmospheric levels for
!                                  which boundary layer fluxes are
!                                  calculated.
     &,RADHR_DIM1             ! IN Dimension of Radiative heating
!                             !    rate (P_FIELD but used for
!                             !    dynamic allocation)

      REAL
     & BQ(P_FIELD,BL_LEVELS)  ! IN A buoyancy parameter for clear air
!                             !    on p,T,q-levels (full levels).
     &,BT(P_FIELD,BL_LEVELS)  ! IN A buoyancy parameter for clear air
!                             !    on p,T,q-levels (full levels).
     &,BQM(P_FIELD,BL_LEVELS) ! IN A buoyancy parameter for clear air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,BTM(P_FIELD,BL_LEVELS) ! IN A buoyancy parameter for clear air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,BQM_CLD(P_FIELD,BL_LEVELS)
!                             ! IN A buoyancy parameter for cloudy air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,BTM_CLD(P_FIELD,BL_LEVELS)
!                             ! IN A buoyancy parameter for cloudy air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,A_QSM(P_FIELD,BL_LEVELS)
!                             ! IN Saturated lapse rate factor
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,A_DQSDTM(P_FIELD,BL_LEVELS)
!                             ! IN Saturated lapse rate factor
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,P(P_FIELD,BL_LEVELS)   ! IN P(*,K) is pressure at full level k.
     &,P_HALF(P_FIELD,BL_LEVELS)
!                             ! IN P_HALF(*,K) is pressure at half
!                             !    level k-1/2.
     &,T(P_FIELD,BL_LEVELS)   ! IN Temperature (K).
     &,Q(P_FIELD,BL_LEVELS)   ! IN Sp humidity (kg water per kg air).
     &,QCF(P_FIELD,BL_LEVELS) ! IN Cloud ice (kg per kg air).
     &,QCL(P_FIELD,BL_LEVELS) ! IN Cloud liquid water (kg per kg air).
     &,CF(P_FIELD,BL_LEVELS)  ! IN Cloud fractions for boundary levs.
     &,DZL(P_FIELD,BL_LEVELS) ! IN Layer depths (m).  DZL(,K) is the
!                                  distance from layer boundary K-1/2
!                                  to layer boundary K+1/2.  For K=1
!                                  the lower boundary is the surface.

      REAL
     & RDZ(P_FIELD,BL_LEVELS) ! IN Reciprocal of distance between
!                                  full levels (m-1).  1/RDZ(,K) is
!                                  the vertical distance from level
!                                  K-1 to level K, except that for
!                                  K=1 it is the height of the
!                                  lowest atmospheric full level.
     &,Z0M(P_FIELD)           ! IN Roughness length for momentum (m).
     &,Z_FULL(P_FIELD,BL_LEVELS) ! IN Z_FULL(*,K) is the height of the
!                                !    k-th full level above the surface.
     &,Z_HALF(P_FIELD,BL_LEVELS) ! IN Z_HALF(*,K) is the height of level
!                                     k-1/2 above the surface (m).
     &,DELTAP(P_FIELD,BL_LEVELS) ! IN Pressure difference between
!                                !    half-levels (Pa).
     &,RAD_HR(RADHR_DIM1,BL_LEVELS)
!                                ! IN Radiative heating rate
!                                !    (Kelvins/s)
     &,V_S(P_FIELD)              ! IN Surface friction velocity (m/s)
     &,FB_SURF(P_FIELD)          ! IN Surface buoyancy flux over density
!                                     (m^2/s^3).
     &,FQW(P_FIELD,BL_LEVELS)    ! IN "Explicit" flux of QW (i.e.
!                                      evaporation) from layer below
!                                      on P-grid (kg per sq m per s).
     &,FTL(P_FIELD,BL_LEVELS)    ! IN "Explicit" flux of TL = H/CP
!                                     (sensible heat/CP) from layer
!                                     below, on P-grid.
     &,TV1_SD(P_FIELD)           ! IN Standard Deviation of level 1
!                                !    virtual temperature (K).
     &,T1_SD(P_FIELD)            ! IN Standard Deviation of level 1
!                                !    temperature (K).
     &,Q1_SD(P_FIELD)            ! IN Standard Deviation of level 1
!                                !    specific humidity (kg/kg).
     &,DQSDT(P_FIELD,BL_LEVELS)  ! IN Partial derivative of QSAT w.r.t.
!                                !    temperature.
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
     &,QW(P_FIELD,BL_LEVELS)     ! IN Total water content (kg per kg
!                                     air).
     &,TL(P_FIELD,BL_LEVELS)     ! IN Liquid/frozen water temperature
!                                     (K).
      REAL
     & ZH(P_FIELD)            ! INOUT Height of the top of the surface
!                             !       based turbulently mixed layer (m).

      REAL
     & DB_SVL(P_FIELD,2:BL_LEVELS)
!                             ! OUT Buoyancy jump across layer
!                                   interface based on SVL.
     &,GRAD_T_ADJ(P_FIELD)    ! OUT Temperature gradient adjustment
!                                   for non-local mixing in unstable
!                                   turbulent boundary layer.
     &,GRAD_Q_ADJ(P_FIELD)    ! OUT Humidity gradient adjustment
!                                   for non-local mixing in unstable
!                                   turbulent boundary layer.
     &,Z_LCL(P_FIELD)         ! OUT Height of lifting condensation
!                                   level.
     &,RHOKM(P_FIELD,2:BL_LEVELS)
!                             ! OUT Non-local turbulent mixing
!                             !     coefficient for momentum.   
     &,RHOKH(P_FIELD,2:BL_LEVELS)
!                             ! OUT Non-local turbulent mixing
!                             !     coefficient for scalars.
     &,RHOKM_TOP(P_FIELD,2:BL_LEVELS)
!                             ! OUT Top-down turbulent mixing
!                             !     coefficient for momentum.
     &,RHOKH_TOP(P_FIELD,2:BL_LEVELS)
!                             ! OUT Top-down turbulent mixing
!                             !     coefficient for scalars. 
     &,ZHT(P_FIELD)           ! OUT Height below which there may be
!                             !     turbulent mixing (m).
     &,BL_TYPE_1(P_FIELD)     ! OUT Indicator set to 1.0 if stable
!                             !     b.l. diagnosed, 0.0 otherwise.
     &,BL_TYPE_2(P_FIELD)     ! OUT Indicator set to 1.0 if Sc over
!                             !     stable surface layer diagnosed,
!                             !     0.0 otherwise.
     &,BL_TYPE_3(P_FIELD)     ! OUT Indicator set to 1.0 if well mixed
!                             !     b.l. diagnosed, 0.0 otherwise.
     &,BL_TYPE_4(P_FIELD)     ! OUT Indicator set to 1.0 if decoupled
!                             !     Sc layer (not over cumulus)
!                             !     diagnosed, 0.0 otherwise.
     &,BL_TYPE_5(P_FIELD)     ! OUT Indicator set to 1.0 if decoupled
!                             !     Sc layer over cumulus diagnosed,
!                             !     0.0 otherwise.
     &,BL_TYPE_6(P_FIELD)     ! OUT Indicator set to 1.0 if a cumulus
!                             !     capped b.l. diagnosed,
!                             !     0.0 otherwise.

      INTEGER
     & NTML(P_FIELD)          ! OUT Number of model levels in the
!                                   turbulently mixed layer.
     &,NTDSC(P_FIELD)         ! OUT Top level for turbulent mixing in
!                             !     cloud layer.

      LOGICAL
     & CUMULUS(P_FIELD)       ! OUT Flag for cumulus in the b.l.
     &,UNSTABLE(P_FIELD)      ! OUT Flag to indicate an unstable
!                                   boundary layer forced from
!                                   the surface.
     &,DSC(P_FIELD)           ! OUT Flag set if decoupled stratocumulus
!                             !     layer found.

                                                                        
!!----------------------------------------------------------------------
!    External references :-
      EXTERNAL TIMER, QSAT, EXCF_NL

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
     &,GRCP=G/CP                ! Adiabatic lapse rate.
     &,LCRCP=LC/CP              ! Latent heat of condensation / CP.
     &,LFRCP=LF/CP              ! Latent heat of fusion / CP.
     &,LS=LC+LF                 ! Latent heat of sublimation.
     &,LSRCP=LS/CP              ! Latent heat of sublimation / CP.
     &)
      REAL A_PLUME,B_PLUME,A_GRAD_ADJ,MAX_T_GRAD
     &    ,MIN_SVL_GRAD,SC_CFTOL,CT_RESID,LGF,FGF
      PARAMETER (
     & A_PLUME=0.4   
     &,B_PLUME=3.26
     &,A_GRAD_ADJ=3.26
     &,MAX_T_GRAD=1.0E-3
     &,MIN_SVL_GRAD=1.0E-3  ! Minimum SVL gradient in a mixed layer.
     &,SC_CFTOL=0.1         ! Cloud fraction required for a
                            ! stratocumulus layer to be diagnosed.
     &,CT_RESID=500.0       ! Approximate parcel cloud top residence
                            ! time (s).
     &,LGF=1.0              ! Adiabatic gradient factor for cloud
!                           ! liquid water.
     &,FGF=0.0)             ! Adiabatic gradient factor for cloud
!                           ! frozen water.

!
!  Define local storage.
!
!  (a) Workspace. 6*BL_LEVELS-1 blocks of real workspace are required
!      plus 1 block of logical.
!
!
      LOGICAL
     & TOPBL(P_FIELD)            ! Flag set when top of boundary layer
!                                ! is reached.
     &,CLOUD_BASE(P_FIELD)       ! Flag set when cloud base is reached.
!
     &,ABOVE_LCL(P_FIELD)        ! Flag set when parcel above LCL.
     &,COUPLED(P_FIELD)          ! Flag to indicate stratocumulus layer
!                                ! weakly coupled to surface (i.e.
!                                ! weakly decoupled).

      REAL
     & QS(P_FIELD)               ! Saturated sp humidity at pressure and
!                                ! temperature of sucessive levels.
     &,QCL_IC_BOT(P_FIELD,BL_LEVELS)  ! In-cloud liquid water content
!                                     ! at the bottom of the model layer
     &,QCF_IC_BOT(P_FIELD,BL_LEVELS)  ! In-cloud frozen water content
!                                     ! at the bottom of the model layer
     &,QCL_IC_TOP(P_FIELD,BL_LEVELS)  ! In-cloud liquid water content
!                                     ! at the top of the model layer.
     &,QCF_IC_TOP(P_FIELD,BL_LEVELS)  ! In-cloud frozen water content
!                                     ! at the top of the model layer.
     &,CFL(P_FIELD,BL_LEVELS)         ! Liquid cloud fraction.
     &,CFF(P_FIELD,BL_LEVELS)         ! Frozen cloud fraction.
     &,DQCLDZ(P_FIELD,BL_LEVELS)      ! Vertical gradient of in-cloud
!                                     ! liquid cloud water in a
!                                     ! well mixed layer.               
     &,DQCFDZ(P_FIELD,BL_LEVELS)      ! Vertical gradient of in-cloud
!                                     ! frozen cloud water in a
!                                     ! well-mixed layer.
     &,CHI_S(P_FIELD,2:BL_LEVELS)     ! Mixing fraction of just saturate
!                                     ! mixture.
     &,DB(P_FIELD,2:BL_LEVELS)        ! Buoyancy jump across layer
!                                     ! interface.
     &,BFLUX_SURF(P_FIELD)            ! Surface buoyancy flux using
!                                     ! clear air buoyancy parameters.  
     &,BFLUX_SAT_SURF(P_FIELD)        ! Surface buoyancy flux using 
!                                     ! cloudy buoyancy parameters.   
     &,DB_TOP(P_FIELD)                ! Buoyancy jump at the top of the
!                                     ! boundary layer.
     &,DB_CLD(P_FIELD,2:BL_LEVELS)    ! In-cloud buoyancy jump across
!                                     ! layer interface.
     &,DF_OVER_CP(P_FIELD,BL_LEVELS)  ! Radiative flux change over layer
!                                     ! divided by c_P.
     &,DF_TOP_OVER_CP(P_FIELD)        ! Radiative flux change at cloud t
!                                     ! divided by c_P.
     &,SVL_PLUME(P_FIELD)             ! Liquid/frozen water static      
!                                     ! energy over CP for a plume      
!                                     ! parcel rising without dilution  
!                                     ! from level 1.
     &,ENV_SVL_KM1(P_FIELD)           ! Density (virtual) static energy
!                                     ! over CP for last layer.
!                                     ! considered.
     &,PAR_SVL_KM1(P_FIELD)           ! Density (virtual) static energy
                                      ! over CP of parcel at last level
!                                     ! considered.
     &,Z_CLD(P_FIELD)                 ! Cloud fraction weighted
!                                     ! thickness of b.l. cloud.
     &,BT_TOP(P_FIELD)                ! Buoyancy parameter at the top of
!                                     ! the turbulently mixed layer.   
     &,BTT_TOP(P_FIELD)               ! In-cloud buoyancy parameter at
!                                     ! the top of the turbulently
!                                     ! mixed layer.  
     &,BTC_TOP(P_FIELD)               ! Cloud fraction weighted buoyancy
!                                     ! parameter at the top of the b.l.
     &,DB_TOP_CLD(P_FIELD)            ! In-cloud buoyancy jump at the
!                                     ! top of the b.l.
     &,CLD_FACTOR(P_FIELD)            ! Fraction of grid box potentially
!                                     ! giving evaporative entrainment.
     &,CHI_S_TOP(P_FIELD)             ! Mixing fraction of just saturate
!                                     ! mixture at top of the b.l.
     &,ZETA_S(P_FIELD)                ! Non-cloudy fraction of mixing
!                                     ! layer for surface forced
!                                     ! entrainment term.
     &,ZETA_R(P_FIELD)                ! Non-cloudy fraction of mixing
!                                     ! layer for cloud top radiative
!                                     ! cooling entrainment term.
     &,ALPHA_R(P_FIELD)               ! Fraction of the cloud top
!                                     ! radiative cooling assumed to
!                                     ! act above the minimum turbulent
!                                     ! flux level.
     &,ZC(P_FIELD)                    ! Cloud depth (not cloud fraction
!                                     ! weighted).
     &,T_LCL(P_FIELD)                 ! Temperature at liftng condensati
!                                     ! level.
     &,P_LCL(P_FIELD)                 ! Pressure at lifting condensation
!                                     ! level.
     &,ZHSC(P_FIELD)                  ! Height of cloud layer top above
!                                     ! surface (m).
     &,DSVL_TOP(P_FIELD)              ! s_VL jump at cloud layer top(K)
     &,DSCDEPTH(P_FIELD)              ! Depth of cloud layer (m).
     &,DB_DSCT(P_FIELD)               ! Buoyancy jump at the top of the
!                                     ! surface mixed layer.
     &,SVL(P_FIELD,BL_LEVELS)         ! Liquid/frozen water virtual
!                                     ! static energy over CP.
     &,DSVL_KP2(P_FIELD)              ! Liquid/frozen water virtual
!                                     ! static energy gradient over CP.
     &,DF_DSCT_OVER_CP(P_FIELD)       ! Radiative flux change at
!                                     ! decoupled stratocumulus top
!                                     ! divided by c_P.
     &,BT_DSCT(P_FIELD)               ! Buoyancy parameter at the top of
!                                     ! the decoupled Sc.
     &,BTT_DSCT(P_FIELD)              ! In-cloud buoyancy parameter at
!                                     ! the top of the decoupled Sc.
     &,BTC_DSCT(P_FIELD)              ! Cloud fraction weighted buoyancy
!                                     ! parameter at the top of the 
!                                     ! decoupled Sc.
     &,DB_DSCT_CLD(P_FIELD)           ! In-cloud buoyancy jump at the
!                                     ! top of the decoupled Sc.
     &,CHI_S_DSCT(P_FIELD)            ! Mixing fraction of just
!                                     ! saturated mixture at top of
!                                     ! the decoupled Sc.
     &,ZETA_R_DSC(P_FIELD)            ! Non-cloudy fraction of decoupled
!                                     ! Sc layer for cloud top radiative
!                                     ! cooling entrainment term.
     &,ALPHA_R_DSC(P_FIELD)           ! Fraction of the cloud top
!                                     ! radiative cooling assumed to
!                                     ! act above the minimum turbulent
!                                     ! flux level for the decoupled Sc.
     &,ZC_DSC(P_FIELD)                ! Cloud depth (not cloud fraction
!                                     ! weighted).
     &,Z_CLD_DSC(P_FIELD)             ! Cloud fraction weighted
!                                     ! thickness of the decoupled Sc
!                                     ! layer.
     &,CLD_FACTOR_DSC(P_FIELD)        ! Fraction of grid box potentially
!                                     ! giving evaporative entrainment.
     &,ZHPAR(P_FIELD)                 ! Temporary store of ZH before LCL
!                                     ! test.
!
      INTEGER
     & K_CLOUD_TOP(P_FIELD)      ! Level number of top of b.l. cloud.
     &,NLCL(P_FIELD)             ! No. of model layers below the lifting
!                                ! condensation level.
     &,NTPAR(P_FIELD)            ! Temporary location to save initial
!                                ! NTML before LCL test.

!
!  (b) Scalars.
!
      REAL
     & VIRT_FACTOR       ! Temporary in calculation of buoyancy
!                        ! parameters.
     &,DTLDZ             ! Vertical gradient of TL in a well-mixed
!                        ! layer.
     &,DQW               ! Total water content change across layer
!                        ! interface.
     &,DTL               ! Liquid/ice static energy change across
!                        ! layer interface.
     &,DQCL              ! Cloud liquid water change across layer
!                        ! interface.
     &,DQCF              ! Cloud frozen water change across layer
!                        ! interface.
     &,VAP_PRESS         ! Vapour pressure.
     &,GRAD_CLD          ! Humidity gradient in layer above LCL.
     &,GRAD_SUB          ! Humidity gradient in layer below LCL.
     &,Q_VAP_PARC        ! Vapour content of parcel.
     &,Q_LIQ_PARC        ! Condensed water content of parcel.
     &,T_PARC            ! Temperature of parcel.
     &,T_DENS_PARC       ! Density potential temperature of parcel.
     &,T_DENS_ENV        ! Density potential temperature of environment.
     &,DENV_BYDZ         ! Gradient of density potential
!                        ! temperature in the environment.         
     &,DPAR_BYDZ         ! Gradient of density potential 
!                        ! temperature of the parcel.
     &,D_SIEMS           ! Siems et al. (1990) cloud-top entrainment
!                        ! instability parameter.
     &,DSVL_KP1          ! s_VL gradient between previous layers.
!
      INTEGER
     & I       ! Loop counter (horizontal field index).
     &,J       ! Offset counter in certain I loops.
     &,K       ! Loop counter (vertical level index).
     &,KM1     ! K-1
     &,MBL     ! Maximum number of model layers allowed in the rapidly
!              ! mixing layer; set to BL_LEVELS-1.
     &,INV_DROP   ! Set to 1 if inversion spread over 2 levels.
     &,K_RAD_SMLT ! Highest SML level for radiative divergence calc.
     &,K_RAD_LIM   
!
!-----------------------------------------------------------------------
!!  0.  Check that the scalars input to define the grid are consistent.
!!      See comments to routine SF_EXCH for details.
!-----------------------------------------------------------------------

      IF (LTIMER) THEN
        CALL TIMER('KMKHZ   ',3)
      ENDIF

!  Set MBL, "maximum number of boundary levels" for the purposes of
!  boundary layer height calculation.

      MBL = BL_LEVELS - 1

!
!-----------------------------------------------------------------------
!! 0.1 Set TOPBL to .FALSE. and calculate boundary layer top using
!!     a non-local, plume method.    
! ----------------------------------------------------------------------
!
      DO I = P1,P1-1+P_POINTS
        SVL_PLUME(I) = TL(I,1) + GRCP * 0.5 * DZL(I,1) + A_PLUME +
     &                 MIN( MAX_T_GRAD * ZH(I) , B_PLUME * TV1_SD(I) )

!
!-----------------------------------------------------------------------
!! 0.2 Calculate temperature and pressure of lifting condensation level
!!     using approximations from Bolton (1980)
!-----------------------------------------------------------------------
!
        VAP_PRESS = Q(I,1) * P(I,1) / (100.0*EPSILON)
        IF (VAP_PRESS .GT .0.0) THEN
          T_LCL(I) =  2840.0/(ALOG(T(I,1))/KAPPA-ALOG(VAP_PRESS)-4.805)
     &                + 55.0
          P_LCL(I) = P(I,1) * (T_LCL(I)/T(I,1))**(1.0/KAPPA)
        ELSE
          P_LCL(I) = P_HALF(I,1)
        ENDIF
        PAR_SVL_KM1(I) = SVL_PLUME(I) * ( 1.0 + C_VIRTUAL*Q(I,1) -
     &                                          QCL(I,1) - QCF(I,1) )
        ENV_SVL_KM1(I) = TL(I,1) + GRCP * 0.5 * DZL(I,1)
        CUMULUS(I)= .FALSE.
        TOPBL(I) = .FALSE.
        NTML(I) = 1
        NLCL(I) = 1
        Z_LCL(I) = Z_HALF(I,1)
      ENDDO
!
!-----------------------------------------------------------------------
!! 0.3 Now compare plume s_VL with each model layer s_VL in turn to
!!     find the first time that plume has negative buoyancy.
!-----------------------------------------------------------------------
!
      DO  K = 2,BL_LEVELS
        CALL QSAT(QS(P1),T(P1,K),P(P1,K),P_POINTS)
        DO I = P1,P1-1+P_POINTS
          IF ( P_LCL(I) .LT. P_HALF(I,K) ) THEN
!           !-----------------------------------------------------------
!           ! Below the lifting condensation level
!           ! potential temperature is constant.
!           !-----------------------------------------------------------
            T_PARC = SVL_PLUME(I) - GRCP * Z_FULL(I,K)
            Q_VAP_PARC = Q(I,1)
            Q_LIQ_PARC = 0.0
            Z_LCL(I) = Z_HALF(I,K)
            NLCL(I) = K-1
            ABOVE_LCL(I) = .FALSE.
          ELSE
!           !-----------------------------------------------------------
!           ! Above the lifting condenstion level
!           ! calculate parcel potential temperature by linearising
!           ! q_sat about the environmental temperature.
!           !-----------------------------------------------------------
            T_PARC = ( SVL_PLUME(I) - GRCP * Z_FULL(I,K) 
     &                 + LCRCP * (QW(I,1) - QS(I) + DQSDT(I,K)*T(I,K)) 
     &               ) / (1.0 + LCRCP*DQSDT(I,K))
            IF (T_PARC .LT. TM) THEN
              T_PARC = ( SVL_PLUME(I) - GRCP * Z_FULL(I,K)
     &                   + LSRCP * (QW(I,1) - QS(I) + DQSDT(I,K)*T(I,K))
     &                 ) / (1.0 + LSRCP*DQSDT(I,K))
          ENDIF
            Q_VAP_PARC = QS(I) + DQSDT(I,K) * (T_PARC - T(I,K))
            Q_LIQ_PARC = QW(I,1) - Q_VAP_PARC
            ABOVE_LCL(I) = .TRUE.
          ENDIF                                                         
!
          T_DENS_PARC = T_PARC *
     &                  (1.0 + C_VIRTUAL*Q_VAP_PARC - Q_LIQ_PARC)
          T_DENS_ENV = T(I,K) *
     &                  (1.0 + C_VIRTUAL*Q(I,K) - QCL(I,K) - QCF(I,K))
!         !-------------------------------------------------------------
!         ! Find vertical gradients in parcel and environment SVL
!         ! (using values from level below (i.e. K-1)).
!         !-------------------------------------------------------------
          DPAR_BYDZ = (T_DENS_PARC + GRCP*Z_FULL(I,K) -
     &                 PAR_SVL_KM1(I)) / (Z_FULL(I,K) - Z_FULL(I,K-1))
          DENV_BYDZ = (T_DENS_ENV + GRCP*Z_FULL(I,K) -
     &                 ENV_SVL_KM1(I)) / (Z_FULL(I,K) - Z_FULL(I,K-1))
          IF ( .NOT.TOPBL(I) .AND.
     &         (  ( T_DENS_PARC-T_DENS_ENV .LE. 0.0 ) .OR.
!
!                      plume non buoyant
!
     &         (ABOVE_LCL(I) .AND. (DENV_BYDZ .GT. 1.25*DPAR_BYDZ)) .OR.
!
!                      or environmental virtual temperature gradient
!                      significantly larger than parcel gradient
!                      above lifting condensation level
!
     &           (K .GT. MBL)
!                      or max no. of model layers in b.l. reached       
     &         )
     &         ) THEN
!
            TOPBL(I) = .TRUE.
            ZH(I) = Z_HALF(I,K)
            NTML(I) = K-1
          ENDIF
          ENV_SVL_KM1(I) = T_DENS_ENV + GRCP*Z_FULL(I,K)
          PAR_SVL_KM1(I) = T_DENS_PARC + GRCP*Z_FULL(I,K)
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!! 0.3a Look for decoupled cloudy mixed layer above b.l. top at ZH:
!!      find cloud base above surface mixed layer inversion, 
!!        i.e. above NTML+1, 
!!      then cloud top (i.e. CF < SC_CFTOL) 
!!      and finally check that cloud is well mixed.  
!-----------------------------------------------------------------------
!      Initialise variables and calculate SVL: conserved variable used 
!      to test for well mixed layers.
!
      DO I = P1,P1-1+P_POINTS
        CLOUD_BASE(I) = .FALSE.
        DSC(I) = .FALSE.
        COUPLED(I) = .FALSE.
        ZHSC(I) = 0.0
        NTDSC(I) = 0
      ENDDO

      DO K = 1,BL_LEVELS
        DO I = P1,P1-1+P_POINTS
         SVL(I,K) = TL(I,K) * ( 1.0 + C_VIRTUAL*QW(I,K) ) +
     &           GRCP * Z_FULL(I,K) 
        ENDDO
      ENDDO
!
      DO K=3,MBL
        DO I = P1,P1-1+P_POINTS
!       !--------------------------------------------------------------
!       ! Find cloud base (where cloud here means CF > SC_CFTOL).
!       !--------------------------------------------------------------
          IF ( K .GT. NTML(I)+1 .AND. CF(I,K) .GT. SC_CFTOL 
     &                          .AND. .NOT.CLOUD_BASE(I)
!                                  not yet found cloud base
     &                          .AND. .NOT.DSC(I) ) THEN
!                                  not yet found a Sc layer
            CLOUD_BASE(I) = .TRUE.
          ENDIF
          IF ( CLOUD_BASE(I) .AND. .NOT.DSC(I) .AND. (
!                  found cloud base but not yet reached cloud top
     &                            (CF(I,K+1).LT.SC_CFTOL) )
!                  got to cloud top
     &       ) THEN
!           !-----------------------------------------------------------
!           ! Look to see if at least top of cloud is well mixed:
!           ! test SVL-gradient for top 2 pairs of levels, in case
!           ! cloud top extends into the inversion.
!           ! Parcel descent in Section 4.0 below will determine depth
!           ! of mixed layer.
!           !----------------------------------------------------------
            IF ( 2.0*(SVL(I,K)-SVL(I,K-1)) / 
     &           (DZL(I,K)+DZL(I,K-1)) .LT. MIN_SVL_GRAD ) THEN
              DSC(I) = .TRUE.
              NTDSC(I) = K
              ZHSC(I)  = Z_HALF(I,K+1)
            ELSEIF ( 2.0*(SVL(I,K-1)-SVL(I,K-2)) /
     &               (DZL(I,K-1)+DZL(I,K-2)) .LT. MIN_SVL_GRAD ) THEN
              DSC(I) = .TRUE.
              NTDSC(I) = K-1
              ZHSC(I)  = Z_HALF(I,K)
            ELSE
              CLOUD_BASE(I) = .FALSE.
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!! 0.4 Save parcel ascent top: this will be used to allow mixing and
!!     entrainment into decoupled Sc of single layer thickness when it
!!     occurs above Cu.
!-----------------------------------------------------------------------
      DO I = P1,P1-1+P_POINTS
        ZHPAR(I) = ZH(I)
        NTPAR(I) = NTML(I)
      ENDDO
! 
!-----------------------------------------------------------------------
!     Test height derived above against lifting condensation level
!-----------------------------------------------------------------------
      DO I = P1,P1-1+P_POINTS
!-----------------------------------------------------------------------
!     Check lifting condensation levels against height of parcel ascent,
!     if lifting condensation level lower than parcel ascent then decide
!     on type of cloudy layer. If lifting condensation level at or below
!     low grid point, assume fog layer and turbulent mixing. For
!     gradient tests assume any if LCL and top of parcel ascent is less
!     than two levels then stratocumulus.
!-----------------------------------------------------------------------
        IF ( (NTML(I)-NLCL(I) .GE. 2) .AND. (NLCL(I) .GT. 1) ) THEN
!-----------------------------------------------------------------------
!     Cloudy boundary layer, diagnose whether stratocumulus or cumulus.
!     For stratocumulus top of mixed layer = ZH
!     For cumulus top of mixed layer = ZLCL
!     Diagnosis is done by comparing gradients
!-----------------------------------------------------------------------
          IF (NTML(I) .EQ. MBL) THEN                                    
            CUMULUS(I) = .TRUE.
            ZH(I) = Z_LCL(I)
            NTML(I) = NLCL(I)
          ELSE
            GRAD_CLD =  ABS( QW(I,NTML(I)) - QW(I,NLCL(I)) ) /
     &                 ( Z_FULL(I,NTML(I)) - Z_FULL(I,NLCL(I)) )
            GRAD_SUB =  ABS( QW(I,NLCL(I)) - QW(I,1) ) /
     &                 ( Z_FULL(I,NLCL(I)) - Z_FULL(I,1) )
            IF (GRAD_CLD .GT. 1.10*GRAD_SUB) THEN
!-----------------------------------------------------------------------
!     Not well mixed, however it is possible that the depth of a well
!     mixed boundary layer has increased but not yet been mixed yet so
!     test gradient from next level down.
!-----------------------------------------------------------------------

              GRAD_CLD =  ABS( QW(I,NTML(I)-1) - QW(I,NLCL(I)) ) /
     &                   ( Z_FULL(I,NTML(I)-1) - Z_FULL(I,NLCL(I)) )

              IF ( GRAD_CLD .GT. 1.10*GRAD_SUB) THEN
!-----------------------------------------------------------------------
!      Diagnoss a cumulus layer and set mixed layer depth to Z_LCL
!-----------------------------------------------------------------------
                IF (P_LCL(I) .LT. (P(I,NLCL(I)))) THEN
!-----------------------------------------------------------------------
!      If LCL is diagnosed in the upper half of the layer set Z_LCL to
!      the height of the upper layer interface
!      (in code above LCL is always set to the lower interface).
!-----------------------------------------------------------------------
                  NLCL(I) = NLCL(I)+1
                  Z_LCL(I) = Z_HALF(I,NLCL(I)+1)
              ENDIF
                CUMULUS(I) = .TRUE.
                ZH(I) = Z_LCL(I)                                        
                NTML(I) = NLCL(I)                                       
            ENDIF
          ENDIF
        ENDIF
        ENDIF                                                           
      ENDDO
!
      CALL QSAT(QS(P1),T(P1,1),P(P1,1),P_POINTS)
!
!-----------------------------------------------------------------------
!!  0.4a Check whether the inversion at NTML is spread over 2 levels.   
!!       If it is, we want the lower as it may indicate a deepening 
!!       layer which will deepen too fast if we jump up a level too
!!       soon. Note that NTML does not change if it is at NLCL (ie. Cu).
!-----------------------------------------------------------------------
!
      DO I = P1,P1-1+P_POINTS
        IF ( NTML(I).GT.1 .AND. .NOT.CUMULUS(I) ) THEN
          K=NTML(I)
          IF ( 2.0*(SVL(I,K)-SVL(I,K-1)) /
     &         (DZL(I,K)+DZL(I,K-1)) .GT. MIN_SVL_GRAD ) THEN
            NTML(I) = K-1
            ZH(I)  = Z_HALF(I,K)
          ENDIF
        ENDIF
      ENDDO
!
!-----------------------------------------------------------------------
!! 0.4aa If no decoupled cloud layer above ZHPAR was found but the layer
!!       to ZHPAR is cloudy, declare this layer a decoupled cloud layer 
!!       and set ZHSC and NTDSC accordingly.  Note that if the 
!!       subsequent tests indicate a well mixed cloud layer to the
!!       surface, the decoupled layer will be switched off.
!-----------------------------------------------------------------------
!
      DO I = P1,P1-1+P_POINTS
        IF ( .NOT.DSC(I) .AND. CF(I,NTPAR(I)) .GT. SC_CFTOL ) THEN
          DSC(I)  = .TRUE.

          IF ( ZHPAR(I) .GT. ZH(I) .AND. NTPAR(I) .GT. 1 ) THEN 
!           ! May still need to lower this inversion, as at section 0.4a
            K=NTPAR(I)
            IF ( 2.0*(SVL(I,K)-SVL(I,K-1)) /
     &           (DZL(I,K)+DZL(I,K-1)) .GT. MIN_SVL_GRAD ) THEN
              NTPAR(I) = K - 1
              ZHPAR(I)  = Z_HALF(I,K)
            ENDIF
          ENDIF
          ZHSC(I) = ZHPAR(I)
          NTDSC(I)= NTPAR(I)
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
!! 0.4b Calculate the radiative flux changes across cloud top for the
!!     stratocumulus layer and thence a first guess for the top-down
!!     mixing depth of this layer, DSCDEPTH.  
!-----------------------------------------------------------------------
!     Initialise variables
!------------------------------
      DO I = P1,P1-1+P_POINTS
        K_CLOUD_TOP(I) = 0
        DF_DSCT_OVER_CP(I) = 0.0
      ENDDO

      DO K = 1,BL_LEVELS-1
        DO I = P1,P1-1+P_POINTS
          
          DF_OVER_CP(I,K) = DELTAP(I,K)/G * RAD_HR(I,K)
!         !-------------------------------------------------------------
!         ! Find the layer with the greatest radiative flux jump and 
!         ! assume that this is the top decoupled Sc layer.  If this is 
!         ! above the surface mixed layer (SML) (which implies strong
!         ! decoupling), limit the search to levels above the SML.
!         !-------------------------------------------------------------
          K_RAD_LIM = 0
          IF ( NTDSC(I) .GT. NTML(I) ) K_RAD_LIM = NTML(I)+1

          IF ( DSC(I) .AND. K .GT. K_RAD_LIM .AND. 
     &         K .LE. NTDSC(I)+2 .AND.
     &         DF_OVER_CP(I,K) .GT. DF_DSCT_OVER_CP(I) ) THEN
            K_CLOUD_TOP(I) = K
            DF_DSCT_OVER_CP(I) = DF_OVER_CP(I,K)
          ENDIF

        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     If cloud top extends into the level above NTDSC (so that the
!     radiative divergence maximum is at NTDSC+1) but there is  
!     additional `cloud top' cooling at NTDSC, add this on to 
!     DF_DSCT_OVER_CP.
!-----------------------------------------------------------------------
      DO I = P1,P1-1+P_POINTS
        IF ( DSC(I) .AND. K_CLOUD_TOP(I) .EQ. NTDSC(I)+1 ) THEN
          K=NTDSC(I)
          DF_OVER_CP(I,K) = DELTAP(I,K)/G * RAD_HR(I,K)
          IF (DF_OVER_CP(I,K) .GT. 0.0 )
     &      DF_DSCT_OVER_CP(I) = DF_DSCT_OVER_CP(I) + DF_OVER_CP(I,K)
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
!     Set DSCDEPTH, the depth of the DSC layer.  Note that this estimate
!     will be the length-scale used to calculate the entrainment rate
!     (although the dependence is only weak) but that a more accurate 
!     plume descent (dependent on w_e) is subsequently used to determine
!     the depth over which the top-down mixing profiles will be applied.
!     If DSC is FALSE, DSCDEPTH = 0.  The plume descent here uses 
!     a radiative perturbation to the cloud-top SVL based roughly on a 
!     typical cloud-top residence time.  If the plume does not sink and 
!     the cloud is decoupled from the surface, then it is assumed to be 
!     stable, i.e. St rather than Sc, and no mixing or entrainment is
!     applied to it.
!-----------------------------------------------------------------------
      DO I = P1,P1-1+P_POINTS
        DSCDEPTH(I) = 0.0
        IF ( DSC(I) ) THEN
          IF (K_CLOUD_TOP(I) .GT. 0) THEN
            SVL_PLUME(I) = SVL(I,NTDSC(I))
     &             + CT_RESID * MIN( RAD_HR(I,K_CLOUD_TOP(I)), 0.0)
          ELSE
            SVL_PLUME(I) = SVL(I,NTDSC(I))
          ENDIF

          IF ( K_CLOUD_TOP(I) .EQ. NTDSC(I)+1 ) 
     &      SVL_PLUME(I)=SVL_PLUME(I)
     &             + CT_RESID * MIN( RAD_HR(I,NTDSC(I)), 0.0)

        ELSE
          SVL_PLUME(I)=0.0
        ENDIF
      ENDDO

      DO K=BL_LEVELS-1,1,-1
        DO I = P1,P1-1+P_POINTS
          IF ( K.LT.NTDSC(I)       ! layer must at least 2 levels deep
     &             .AND. SVL_PLUME(I)-SVL(I,K) .LT. 0.01 ) THEN
            DSCDEPTH(I) = ZHSC(I) - Z_HALF(I,K)
          ENDIF
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!! 0.4c If the layer to ZH is cloud capped and there is no Cu or 
!!     decoupled Sc above (indicated by NTML=NTDSC), set spare logical 
!!     CLOUD_BASE to true and try looking for a weak, 
!!     possibly diurnally decoupling, inversion below ZH using local 
!!     gradients method.  If one is found, lower ZH to this decoupling 
!!     inversion leaving ZHSC at cloud-top.  
!!     If DSCDEPTH was found to be zero above, set DSCDEPTH=ZHSC-ZH.
!-----------------------------------------------------------------------
      K=1
      DO I = P1,P1-1+P_POINTS
        CLOUD_BASE(I) = .FALSE.
        IF ( NTML(I) .EQ. NTDSC(I) ) THEN
!-----------------------------------------------------------------------
!       Well mixed cloudy layer mixing up to ZH: reset decoupled layer 
!       parameters ready to look for decoupling inversion below ZH.
!-----------------------------------------------------------------------
          DSC(I)  = .FALSE.
          ZHSC(I) = 0.0
          NTDSC(I)= 0
          CLOUD_BASE(I) = .TRUE.
        ENDIF
        DSVL_KP2(I)= 2.0*(SVL(I,K+2) - SVL(I,K+1)) /
     &                   (DZL(I,K+2) + DZL(I,K+1))
      ENDDO

      DO K=2,MBL-1
        DO I = P1,P1-1+P_POINTS
!----------------------------------------------------------------------
!! 0.4d Local gradients method: find where s_VL gradient decreases from 
!!      a non-trivial positive value (ie. greater than MAX_T_GRAD).
!!      The inversion is then at the max s_VL-gradient height.
!----------------------------------------------------------------------
          DSVL_KP1 = DSVL_KP2(I)
          DSVL_KP2(I) = 2.0*(SVL(I,K+2) - SVL(I,K+1)) /
     &                      (DZL(I,K+2) + DZL(I,K+1))
          IF ( K.LT.NTML(I) .AND. .NOT.DSC(I) .AND. CLOUD_BASE(I) .AND. 
!                     Below apparently well-mixed cloud-layer top
     &      (DSVL_KP1 .GT. MAX_T_GRAD .AND. DSVL_KP2(I) .LT. DSVL_KP1)
!                     and gradient conditions satisfied
     &       ) THEN
            DSC(I) = .TRUE.
            NTDSC(I) = NTML(I)
            ZHSC(I)  = Z_HALF(I,NTDSC(I)+1)
            NTML(I) = K
            ZH(I)  = Z_HALF(I,NTML(I)+1)
!           !-----------------------------------------------------------
!           ! As in section 0.4c, check whether this ZH is at the base
!           ! of the inversion and if it is not, lower it.
!           !-----------------------------------------------------------
            IF ( 2.0*(SVL(I,K) - SVL(I,K-1)) /
     &               (DZL(I,K) + DZL(I,K-1)) .GT. MIN_SVL_GRAD ) THEN
              NTML(I) = K-1
              ZH(I) = Z_HALF(I,K)
            ENDIF
!           !-----------------------------------------------------------
!           ! Reset DSCDEPTH, if top-down mixing is only weak at present
!           !-----------------------------------------------------------
            IF ( DSCDEPTH(I) .EQ. 0.0 ) 
     &               DSCDEPTH(I) = ZHSC(I) - Z_HALF(I,NTDSC(I)-1)
          ENDIF
        ENDDO
      ENDDO
!----------------------------------------------------------------------
!!  0.4e Tidy up variables associated with decoupled layer.
!----------------------------------------------------------------------
      DO I = P1,P1-1+P_POINTS
        IF (CUMULUS(I) .AND. DSC(I)) 
     &    DSCDEPTH(I)=MIN( DSCDEPTH(I),ZHSC(I)-Z_HALF(I,NTML(I)+2) )
        IF ( DSCDEPTH(I) .EQ. 0.0 ) THEN
          IF (NTDSC(I) .EQ. NTPAR(I)) THEN
!           !----------------------------------------------------------
!           ! Indicates a Sc layer over Cu: force mixing over single
!           ! layer.
!           !----------------------------------------------------------
            DSCDEPTH(I) = DZL(I,NTDSC(I))
          ELSE
            DSC(I)=.FALSE.
            NTDSC(I)=0
            ZHSC(I)=0.0
            DF_DSCT_OVER_CP(I) = 0.0
          ENDIF
        ENDIF
!       !---------------------------------------------------------------
!       ! These next two tests should not be necessary but put in to be
!       ! fail-safe.
!       !---------------------------------------------------------------
        IF (.NOT.DSC(I)) DSCDEPTH(I) = 0.0
        IF (NTDSC(I) .EQ. NTML(I)) THEN
!         ! Clearly modelling the same layer 
          DSC(I) = .FALSE.
          NTDSC(I) = 0
          ZHSC(I) = 0.0
          DF_DSCT_OVER_CP(I) = 0.0
          DSCDEPTH(I) = 0.0
        ENDIF
      ENDDO
!
!----------------------------------------------------------------------
!     If decoupled cloud layer found test to see if it is, in fact, 
!     coupled to the surface mixed-layer: if SVL difference between 
!     NTML and NTDSC is less than 0.2K then assume coupled.  This will 
!     mean that the surface-driven entrainment term will be applied at 
!     ZHSC, no entrainment will be applied at ZH and ZHSC will be the 
!     lengthscale used in the entrainment inputs.
!----------------------------------------------------------------------
      DO I = P1,P1-1+P_POINTS
        COUPLED(I) = .FALSE.
        IF ( DSC(I) ) THEN
!         !------------------------------------------------------------
!         ! Note this IF test structure is required because if DSC is
!         ! false then NTDSC = 0 and cannot be used to index SVL.
!         !------------------------------------------------------------
          IF ( SVL(I,NTDSC(I)) - SVL(I,NTML(I)) .LT. 0.2 ) 
     &       COUPLED(I) = .TRUE.
        ENDIF
      ENDDO
!
      DO I = P1,P1+P_POINTS-1
!
!-----------------------------------------------------------------------
!! 0.5 Calculate the within-layer vertical gradients of cloud liquid
!!     and frozen water for the layer 1
!-----------------------------------------------------------------------
!
        VIRT_FACTOR = 1.0 + C_VIRTUAL*Q(I,1) - QCL(I,1) - QCF(I,1)
!
        GRAD_T_ADJ(I) = MIN( MAX_T_GRAD ,
     &                       A_GRAD_ADJ * T1_SD(I) / ZH(I) )
!        IF (T1_SD(I) .GT. 0.0) THEN
!          GRAD_Q_ADJ(I) = (Q1_SD(I) / T1_SD(I)) * GRAD_T_ADJ(I)
!        ELSE
          GRAD_Q_ADJ(I) = 0.0
!        ENDIF
        DTLDZ = -GRCP + GRAD_T_ADJ(I)
        DQCLDZ(I,1) = -LGF * ( DTLDZ*DQSDT(I,1) +
     &                   G*QS(I)/(R*T(I,1)*VIRT_FACTOR) )
     &                  / (1.0 + LCRCP*DQSDT(I,1))
        DQCFDZ(I,1) = -FGF * ( DTLDZ*DQSDT(I,1) +
     &                   G*QS(I)/(R*T(I,1)*VIRT_FACTOR) )
     &                  / (1.0 + LSRCP*DQSDT(I,1))
!
!-----------------------------------------------------------------------
!! 0.6 Calculate the cloud liquid and frozen water contents at the
!!     top and bottom of layer 1
!-----------------------------------------------------------------------
!
        IF ( QCL(I,1) + QCF(I,1) .GT. 0.0 ) THEN
          CFL(I,1) = CF(I,1) * QCL(I,1)/(QCL(I,1)+QCF(I,1))
          CFF(I,1) = CF(I,1) * QCF(I,1)/(QCL(I,1)+QCF(I,1))
        ELSE
          CFL(I,1) = 0.0
          CFF(I,1) = 0.0
        ENDIF
!
        IF (CFL(I,1) .GT. 0.0) THEN
          QCL_IC_TOP(I,1) = QCL(I,1) / CFL(I,1) +
     &                       0.5*DZL(I,1)*DQCLDZ(I,1)
        ELSE
          QCL_IC_TOP(I,1) = 0.0
        ENDIF
!
        IF (CFF(I,1) .GT. 0.0) THEN
          QCF_IC_TOP(I,1) = QCF(I,1) / CFF(I,1) +
     &                       0.5*DZL(I,1)*DQCFDZ(I,1)
        ELSE
          QCF_IC_TOP(I,1) = 0.0
        ENDIF
!
        QCL_IC_BOT(I,1) = 0.0
        QCF_IC_BOT(I,1) = 0.0
!
      ENDDO
!
!-----------------------------------------------------------------------
!! 1.  First loop round boundary layer levels.
!-----------------------------------------------------------------------
!
      DO K=2,BL_LEVELS
!
        CALL QSAT(QS(P1),T(P1,K),P(P1,K),P_POINTS)
!
        DO I=P1,P1+P_POINTS-1
!
!-----------------------------------------------------------------------
!! 1.4 Calculate the within-layer vertical gradients of cloud liquid
!!     and frozen water for the current layer
!-----------------------------------------------------------------------
!
          VIRT_FACTOR = 1.0 + C_VIRTUAL*Q(I,K) - QCL(I,K) - QCF(I,K)
!
          IF (K. LE. NTML(I)) THEN
            DTLDZ = -GRCP + GRAD_T_ADJ(I)
          ELSE
            DTLDZ = -GRCP
          ENDIF
!
          DQCLDZ(I,K) = -LGF * ( DTLDZ*DQSDT(I,K)
     &                   + G*QS(I)/(R*T(I,K)*VIRT_FACTOR) )
     &                    / ( 1.0 + LCRCP*DQSDT(I,K) )
          DQCFDZ(I,K) = -FGF * ( DTLDZ*DQSDT(I,K)
     &                   + G*QS(I)/(R*T(I,K)*VIRT_FACTOR) )
     &                    / ( 1.0 + LSRCP*DQSDT(I,K) )
!
!-----------------------------------------------------------------------
!! 1.5 Calculate the cloud liquid and frozen water contents at the
!!     top and bottom of the current layer
!-----------------------------------------------------------------------
!
          IF ( QCL(I,K) + QCF(I,K) .GT. 0.0 ) THEN
            CFL(I,K) = CF(I,K) * QCL(I,K)/( QCL(I,K) + QCF(I,K) )
            CFF(I,K) = CF(I,K) * QCF(I,K)/( QCL(I,K) + QCF(I,K) )
          ELSE
            CFL(I,K) = 0.0
            CFF(I,K) = 0.0
          ENDIF
!
          IF (CFL(I,K) .GT. 0.0) THEN
            QCL_IC_TOP(I,K) = QCL(I,K) / CFL(I,K) +
     &                         0.5*DZL(I,K)*DQCLDZ(I,K)
            QCL_IC_BOT(I,K) = MAX( 0.0 , QCL(I,K) / CFL(I,K) -
     &                                    0.5*DZL(I,K)*DQCLDZ(I,K) )
          ELSE
            QCL_IC_TOP(I,K) = 0.0
            QCL_IC_BOT(I,K) = 0.0
          ENDIF
!
          IF (CFF(I,K) .GT. 0.0) THEN
            QCF_IC_TOP(I,K) = QCF(I,K) / CFF(I,K) +
     &                         0.5*DZL(I,K)*DQCFDZ(I,K)
            QCF_IC_BOT(I,K) = MAX( 0.0 , QCF(I,K) / CFF(I,K) -
     &                                    0.5*DZL(I,K)*DQCFDZ(I,K) )
          ELSE
            QCF_IC_TOP(I,K) = 0.0
            QCF_IC_BOT(I,K) = 0.0
          ENDIF
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!! 2.  Second loop round boundary layer levels.
!-----------------------------------------------------------------------
!
      DO K=2,BL_LEVELS
        KM1 = K-1
        DO I=P1,P1+P_POINTS-1
!
!-----------------------------------------------------------------------
!! 2.1 Calculate the jumps of QW and TL across the layer interface
!!     at level k-1/2.
!-----------------------------------------------------------------------
!
          DQW = QW(I,K) - QW(I,KM1)                    ! Used in P243.C2
!-----------------------------------------------------------------------
!         ! TL_K_BOT   = TL(I,K)   - 0.5*DZL(I,K)  *(-GRCP)
!         ! TL_KM1_TOP = TL(I,KM1) + 0.5*DZL(I,KM1)*(-GRCP)
!         ! DTL = TL_K_BOT - TL_KM1_TOP   so therefore
!-----------------------------------------------------------------------
          DTL = TL(I,K) - TL(I,KM1) + GRCP/RDZ(I,K)    ! Used in P243.C2
          IF (K .LE. NTML(I)) THEN
            DTL = DTL - GRAD_T_ADJ(I)/RDZ(I,K)
            DQW = DQW - GRAD_Q_ADJ(I)/RDZ(I,K)
          ENDIF
!
          DQCL = CFL(I,K)*QCL_IC_BOT(I,K) -
     &             CFL(I,KM1)*QCL_IC_TOP(I,KM1)
          DQCF = CFF(I,K)*QCF_IC_BOT(I,K) -
     &             CFF(I,KM1)*QCF_IC_TOP(I,KM1)
!
!-----------------------------------------------------------------------
!! 2.3 Calculate the buoyancy jumps across the interface between layers
!!     k and k-1
!-----------------------------------------------------------------------
!
          DB(I,K) = G * ( BTM(I,KM1)*DTL + BQM(I,KM1)*DQW +
     &                    (LCRCP*BTM(I,KM1) - ETAR*BQM(I,KM1)) * DQCL +
     &                    (LSRCP*BTM(I,KM1) - ETAR*BQM(I,KM1)) * DQCF )
!
          DB_SVL(I,K) = G * ( BTM(I,KM1)*DTL + BQM(I,KM1)*DQW )
!
          DB_CLD(I,K) = G * ( BTM_CLD(I,KM1)*DTL + BQM_CLD(I,KM1)*DQW )
!
          CHI_S(I,K) = -QCL_IC_TOP(I,KM1) /
     &                  (A_QSM(I,KM1)*DQW + A_DQSDTM(I,KM1)*DTL)
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!! 3.0a Calculate surface buoyancy flux
!-----------------------------------------------------------------------
!
        DO I = P1,P1-1+P_POINTS

        BFLUX_SURF(I) = G * ( BTM(I,NTML(I))*FTL(I,1) + 
     &                        BQM(I,NTML(I))*FQW(I,1) )

        IF ( BFLUX_SURF(I) . GT. 0.0 ) THEN
          BFLUX_SAT_SURF(I) = G * ( BTM_CLD(I,NTML(I))*FTL(I,1) +
     &                              BQM_CLD(I,NTML(I))*FQW(I,1) )
          IF ( COUPLED(I) ) THEN
             BFLUX_SAT_SURF(I) = G * ( BTM_CLD(I,NTDSC(I))*FTL(I,1) +
     &                                 BQM_CLD(I,NTDSC(I))*FQW(I,1) )
          ENDIF
        ELSE
          BFLUX_SAT_SURF(I) = 0.0
          ENDIF

        UNSTABLE(I) = (BFLUX_SURF(I) .GT. 0.0) .AND. (NTML(I) .GT. 1)

        ENDDO
!-----------------------------------------------------------------------
!! 3.0aa Calculate Sc layer cloud depth (not cloud fraction weighted).
!!       (If DSC(I)=.FALSE. then NTDSC=0 and ZC_DSC remains equal to 0.)
!-----------------------------------------------------------------------
!
!     Initialise variables
!
      DO I = P1,P1-1+P_POINTS
        ZC_DSC(I) = 0.0
        CLOUD_BASE(I) = .TRUE.
      ENDDO
!
      DO K=BL_LEVELS,2,-1
        KM1=K-1
        DO I = P1,P1-1+P_POINTS
          IF ( (K .EQ. NTDSC(I)) .AND. (CF(I,K) .GT. SC_CFTOL) )
     &      CLOUD_BASE(I) = .FALSE.
          IF ( (K .LE. NTDSC(I)) .AND. (CF(I,K) .GT. SC_CFTOL) .AND.
     &        .NOT.CLOUD_BASE(I) ) THEN
            ZC_DSC(I) = ZC_DSC(I) + 0.5*DZL(I,K)
            IF (CF(I,KM1) .GT. SC_CFTOL) THEN
              ZC_DSC(I) = ZC_DSC(I) + 0.5*DZL(I,K)
        ELSE
              ZC_DSC(I) = ZC_DSC(I) +
     &                  MIN( 0.5*(DZL(I,K) + DZL(I,KM1)) * CFL(I,K) ,
     &                        QCL(I,K) / DQCLDZ(I,K) ) / CF(I,K)
              IF (DQCFDZ(I,K) .GT. 0.0) THEN
                ZC_DSC(I) = ZC_DSC(I) +
     &                   MIN( 0.5*(DZL(I,K) + DZL(I,KM1)) * CFF(I,K) ,
     &                        QCF(I,K) / DQCFDZ(I,K) ) / CF(I,K)
              ELSE
                ZC_DSC(I) = ZC_DSC(I) +
     &                0.5*(DZL(I,K) + DZL(I,KM1)) * CFF(I,K) / CF(I,K)
              ENDIF
              CLOUD_BASE(I) = .TRUE.
        ENDIF
            ENDIF
      ENDDO
      ENDDO

      DO I = P1,P1-1+P_POINTS
        IF ( DSC(I) .AND. .NOT.CLOUD_BASE(I)
     &              .AND. (CF(I,1) .GT. SC_CFTOL) ) THEN
          ZC_DSC(I) = ZC_DSC(I) + 0.5*DZL(I,1) +
     &               MIN( 0.5*DZL(I,1)*CFL(I,1) ,
     &                    QCL(I,1) / DQCLDZ(I,1) ) / CF(I,1)
          IF(DQCFDZ(I,1) .GT. 0.0) THEN
            ZC_DSC(I) = ZC_DSC(I) +
     &               MIN( 0.5*DZL(I,1)*CFF(I,1) ,
     &                    QCF(I,1) / DQCFDZ(I,1) ) / CF(I,1)
          ELSE
            ZC_DSC(I) = ZC_DSC(I) + 0.5*DZL(I,1)*CFF(I,1)/CF(I,1)
          ENDIF
          CLOUD_BASE(I) = .TRUE.
        ENDIF
      ENDDO
!     !-----------------------------------------------------------------
!     !  Check for cloud within inversion.
!     !-----------------------------------------------------------------
      DO I = P1,P1-1+P_POINTS
        IF ( DSC(I) .AND. (CF(I,NTDSC(I)+1) .GT. SC_CFTOL) ) THEN
            K=NTDSC(I)+1
            ZC_DSC(I) = ZC_DSC(I) + 0.5*DZL(I,K)
            IF (CF(I,K-1) .GT. SC_CFTOL) THEN
              ZC_DSC(I) = ZC_DSC(I) + 0.5*DZL(I,K)
            ELSE
              ZC_DSC(I) = ZC_DSC(I) +
     &                  MIN( 0.5*(DZL(I,K) + DZL(I,K-1)) * CFL(I,K) ,
     &                        QCL(I,K) / DQCLDZ(I,K) ) / CF(I,K)
              IF (DQCFDZ(I,K) .GT. 0.0) THEN
                ZC_DSC(I) = ZC_DSC(I) +
     &                   MIN( 0.5*(DZL(I,K) + DZL(I,K-1)) * CFF(I,K) ,
     &                        QCF(I,K) / DQCFDZ(I,K) ) / CF(I,K)
              ELSE
                ZC_DSC(I) = ZC_DSC(I) +
     &                0.5*(DZL(I,K) + DZL(I,K-1)) * CFF(I,K) / CF(I,K)
              ENDIF
            ENDIF
        ENDIF
      ENDDO
!     !-----------------------------------------------------------------
!     !  Layer cloud depth cannot be > the layer depth itself.
!     !-----------------------------------------------------------------
      DO I = P1,P1-1+P_POINTS
        ZC_DSC(I) = MIN( ZC_DSC(I), DSCDEPTH(I) )
      ENDDO
!-----------------------------------------------------------------------
!! 3.0b Calculate surface mixed layer cloud depth (not cloud fraction 
!!      weighted)
!-----------------------------------------------------------------------
!
      DO I = P1,P1-1+P_POINTS
        ZC(I) = 0.0
        CLOUD_BASE(I) = .TRUE.
      ENDDO
!
      DO K=BL_LEVELS-1,2,-1
        KM1=K-1
        DO I = P1,P1-1+P_POINTS
          IF ( (K .EQ. NTML(I)) .AND. (CF(I,K) .GT. 0.0) )
     &      CLOUD_BASE(I) = .FALSE.
          IF ( (K .LE. NTML(I)) .AND. (CF(I,K) .GT. 0.0) .AND.
     &        .NOT.CLOUD_BASE(I) ) THEN
            ZC(I) = ZC(I) + 0.5*DZL(I,K)
            IF (CF(I,KM1) .GT. 0.0) THEN
              ZC(I) = ZC(I) + 0.5*DZL(I,K)
            ELSE
              ZC(I) = ZC(I) +
     &                   MIN( 0.5*(DZL(I,K)+DZL(I,KM1))*CFL(I,K) ,
     &                        QCL(I,K) / DQCLDZ(I,K) ) / CF(I,K)
              IF (DQCFDZ(I,K) .GT. 0.0) THEN
                ZC(I) = ZC(I) +
     &                   MIN( 0.5*(DZL(I,K)+DZL(I,KM1))*CFF(I,K) ,
     &                        QCF(I,K) / DQCFDZ(I,K) ) / CF(I,K)
              ELSE
                ZC(I) = ZC(I) +
     &                   0.5*(DZL(I,K)+DZL(I,KM1))*CFF(I,K)/CF(I,K)
              ENDIF
              CLOUD_BASE(I) = .TRUE.
            ENDIF
            ENDIF
        ENDDO
      ENDDO
      DO I = P1,P1-1+P_POINTS
        IF ( (CF(I,1) .GT. 0.0) .AND. .NOT.CLOUD_BASE(I) ) THEN
          ZC(I) = ZC(I) + 0.5*DZL(I,1) +
     &               MIN( 0.5*DZL(I,1)*CFL(I,1) ,
     &                    QCL(I,1) / DQCLDZ(I,1) ) / CF(I,1)
          IF (DQCFDZ(I,1) .GT. 0.0) THEN
            ZC(I) = ZC(I) + 
     &               MIN( 0.5*DZL(I,1)*CFF(I,1) ,
     &                    QCF(I,1) / DQCFDZ(I,1) ) / CF(I,1)
          ELSE
            ZC(I) = ZC(I) + 0.5*DZL(I,1)*CFF(I,1) / CF(I,1)
          ENDIF
          CLOUD_BASE(I) = .TRUE.
        ENDIF
      ENDDO
!     !-----------------------------------------------------------------
!     !  Check for cloud within inversion.
!     !-----------------------------------------------------------------
      DO I = P1,P1-1+P_POINTS
        IF ( CF(I,NTML(I)+1) .GT. SC_CFTOL ) THEN
            K=NTML(I)+1
            ZC(I) = ZC(I) + 0.5*DZL(I,K)
            IF (CF(I,K-1) .GT. SC_CFTOL) THEN
              ZC(I) = ZC(I) + 0.5*DZL(I,K)
            ELSE
              ZC(I) = ZC(I) +
     &                  MIN( 0.5*(DZL(I,K) + DZL(I,K-1)) * CFL(I,K) ,
     &                        QCL(I,K) / DQCLDZ(I,K) ) / CF(I,K)
              IF (DQCFDZ(I,K) .GT. 0.0) THEN
                ZC(I) = ZC(I) +
     &                   MIN( 0.5*(DZL(I,K) + DZL(I,K-1)) * CFF(I,K) ,
     &                        QCF(I,K) / DQCFDZ(I,K) ) / CF(I,K)
              ELSE
                ZC(I) = ZC(I) +
     &                0.5*(DZL(I,K) + DZL(I,K-1)) * CFF(I,K) / CF(I,K)
              ENDIF
            ENDIF
        ENDIF
      ENDDO
!
!-----------------------------------------------------------------------
!! 3.1 Calculate inputs for the top of b.l. entrainment parametrization.
!-----------------------------------------------------------------------
!
      DO I = P1,P1-1+P_POINTS
       ZETA_R_DSC(I) = 0.0
       ALPHA_R_DSC(I) = 0.0
       CHI_S_DSCT(I) = 0.0
       CLD_FACTOR_DSC(I) = 0.0
       BT_DSCT(I) = 0.0
       BTT_DSCT(I) = 0.0
       BTC_DSCT(I) = 0.0
       DB_DSCT(I) = 0.0
       DB_DSCT_CLD(I) = 0.0
       CHI_S_TOP(I) = 0.0
       CLD_FACTOR(I) = 0.0
       BT_TOP(I) = 0.0
       BTT_TOP(I) = 0.0
       BTC_TOP(I) = 0.0
       DB_TOP(I) = 0.0
       DB_TOP_CLD(I) = 0.0    ! default required if COUPLED
       Z_CLD(I) = 0.0
       Z_CLD_DSC(I) = 0.0
      ENDDO
!
      DO K = 1,BL_LEVELS
      DO I = P1,P1-1+P_POINTS
          IF ( K .LE. NTML(I)+1 )  THEN
!           !---------------------------------------------------
!           ! Calculation of cloud fraction weighted
!           ! thickness of cloud in the surface mixed layer
!           !---------------------------------------------------
            Z_CLD(I) = Z_CLD(I) +
     &                 CF(I,K) * 0.5 * DZL(I,K) +
     &                  MIN( CFL(I,K) * 0.5 * DZL(I,K) ,
     &                          QCL(I,K) / DQCLDZ(I,K) ) 
            IF (DQCFDZ(I,K) .GT. 0.0) THEN
              Z_CLD(I) = Z_CLD(I) +
     &                  MIN( CFF(I,K) * 0.5 * DZL(I,K) ,
     &                          QCF(I,K) / DQCFDZ(I,K) )
            ELSE
              Z_CLD(I) = Z_CLD(I) +
     &                    CFF(I,K) * 0.5 * DZL(I,K)
            ENDIF
          ENDIF
!
          IF ( DSC(I) .AND. K.LE.NTDSC(I)+1 .AND. 
     &         ( COUPLED(I) .OR. 
     &               Z_HALF(I,K+1).GE.ZHSC(I)-ZC_DSC(I) ) ) THEN
!           !----------------------------------------------------
!           ! Calculation of cloud fraction weighted thickness of 
!           ! cloud in the DSC layer (or to the surface if COUPLED)
!           !----------------------------------------------------
            Z_CLD_DSC(I) = Z_CLD_DSC(I) +
     &                 CF(I,K) * 0.5 * DZL(I,K) +
     &                  MIN( CFL(I,K) * 0.5 * DZL(I,K) ,
     &                          QCL(I,K) / DQCLDZ(I,K) ) 
            IF (DQCFDZ(I,K) .GT. 0.0) THEN
              Z_CLD_DSC(I) = Z_CLD_DSC(I) +
     &                  MIN( CFF(I,K) * 0.5 * DZL(I,K) ,
     &                          QCF(I,K) / DQCFDZ(I,K) )
            ELSE
              Z_CLD_DSC(I) = Z_CLD_DSC(I) +
     &                        CFF(I,K) * 0.5 * DZL(I,K) 
            ENDIF
          ENDIF
!
          IF (K .EQ. NTML(I)  .AND. .NOT.COUPLED(I) ) THEN
!           !------------------------------------------------------
!           ! Calculation of SML inputs.  If COUPLED then these are 
!           ! not used (as no entrainment is then applied at ZH)
!           !------------------------------------------------------
            CHI_S_TOP(I) = MAX( 0.0, MIN( CHI_S(I,K+1), 1.) ) 
            CLD_FACTOR(I) = MAX( 0.0 , CF(I,K)-CF(I,K+1) )
            BT_TOP(I) = G * BTM(I,K)
            BTT_TOP(I) = G * BTM_CLD(I,K)
            BTC_TOP(I) = BTT_TOP(I)
            DB_TOP(I) = DB(I,K+1)
            DB_TOP_CLD(I) = DB_CLD(I,K+1)
          ENDIF
          IF (K .EQ. NTDSC(I)) THEN
!           !---------------------------------------------------
!           ! Calculation of DSC inputs 
!           ! (if DSC=.FALSE. then K never equals NTDSC(=0))
!           !---------------------------------------------------
            CHI_S_DSCT(I) = MAX( 0.0, MIN( CHI_S(I,K+1), 1.) )
            CLD_FACTOR_DSC(I) = MAX( 0.0 , CF(I,K)-CF(I,K+1) )
            BT_DSCT(I) = G * BTM(I,K)
            BTT_DSCT(I) = G * BTM_CLD(I,K)
            BTC_DSCT(I) = BTT_DSCT(I)
            DB_DSCT(I) = DB(I,K+1)
            DB_DSCT_CLD(I) = DB_CLD(I,K+1)
          ENDIF
        ENDDO
      ENDDO
!     !-----------------------------------------------------------------
!     ! Next those terms dependent on the presence of buoyancy reversal.
!     !"----------------------------------------------------------------
      DO I = P1,P1-1+P_POINTS
        Z_CLD(I) = MIN( Z_CLD(I), ZH(I) )
        Z_CLD_DSC(I) = MIN( Z_CLD_DSC(I), ZHSC(I) )
!       !---------------------------------------------------------------
!       ! First the surface mixed layer.
!       !---------------------------------------------------------------
        IF ( COUPLED(I) ) THEN
          ZETA_S(I) = 1.0 - Z_CLD_DSC(I) / ZHSC(I)
          ZETA_R(I) = 1.0 - ZC_DSC(I) / ZHSC(I)    
        ELSE
          ZETA_S(I) = 1.0 - Z_CLD(I) / ZH(I)
          ZETA_R(I) = 1.0 - ZC(I) / ZH(I)     
        ENDIF
!
        IF (DB_TOP_CLD(I) .GE. 0.0) THEN
!         !--------------------------------------------------
!         ! i.e. no buoyancy reversal (or default if COUPLED)
!         !--------------------------------------------------
          ALPHA_R(I) = 0.2
          DB_TOP_CLD(I) = 0.0
        ELSE
!         !----------------------------
!         ! IF (DB_TOP_CLD(I) .LT. 0.0)
!         ! i.e. buoyancy reversal
!         !----------------------------
          DB_TOP_CLD(I) = -DB_TOP_CLD(I) * CLD_FACTOR(I)
          D_SIEMS = MAX( 0.0, 
     &             CHI_S_TOP(I) * DB_TOP_CLD(I) / (DB_TOP(I)+1.E-14) )
          ZETA_R(I) = MIN( ZETA_R(I)+10.0*(1.0-ZETA_R(I))
     &                                    *D_SIEMS, 1.0 )
          ALPHA_R(I) = MIN( 0.2+10.0*(1.0-0.2)*D_SIEMS, 1.0 )
        ENDIF
!       !---------------------------------------------------------------
!       ! Now the decoupled Sc layer (DSC).
!       !---------------------------------------------------------------
        IF (DSC(I)) THEN
          IF ( COUPLED(I) ) THEN
            ZETA_R_DSC(I) = 1.0 - ZC_DSC(I) / ZHSC(I)     
          ELSE
            ZETA_R_DSC(I) = 1.0 - ZC_DSC(I) / DSCDEPTH(I)
          ENDIF

          IF (DB_DSCT_CLD(I) .GE. 0.0) THEN
!           !----------------------------
!           ! i.e. no buoyancy reversal
!           !----------------------------
            ALPHA_R_DSC(I) = 0.2
            DB_DSCT_CLD(I) = 0.0
          ELSE
!           !----------------------------
!           ! IF (DB_DSCT_CLD(I) .LT. 0.0)
!           ! i.e. buoyancy reversal
!           !----------------------------
            DB_DSCT_CLD(I) = -DB_DSCT_CLD(I) * CLD_FACTOR_DSC(I)
            D_SIEMS = MAX( 0.0, 
     &           CHI_S_DSCT(I) * DB_DSCT_CLD(I) / (DB_DSCT(I)+1.E-14) )
            ZETA_R_DSC(I) = MIN( ZETA_R_DSC(I)+10.0*(1.0-ZETA_R_DSC(I))
     &                                    *D_SIEMS, 1.0 )
            ALPHA_R_DSC(I) = MIN( 0.2+10.0*(1.0-0.2)*D_SIEMS, 1.0 )
          ENDIF
        ENDIF
      ENDDO
!
!-----------------------------------------------------------------------
!! 4. Calculate the radiative flux change across cloud top for mixed
!!    layer to ZH.  If there is a decoupled layer above, restrict 
!!    search for maximum divergence to below NTML+1.  This may 
!!    introduce errors if NTML changes during radiative timestep but 
!!    can't be helped.
!-----------------------------------------------------------------------
!     Initialise variables
!------------------------------
      DO I = P1,P1-1+P_POINTS
        K_CLOUD_TOP(I) = 0
        DF_TOP_OVER_CP(I) = 0.0
      ENDDO
!
      DO K = 1,BL_LEVELS-1
        DO I = P1,P1-1+P_POINTS
          
          DF_OVER_CP(I,K) = DELTAP(I,K)/G * RAD_HR(I,K)

          K_RAD_SMLT = NTML(I)+1
          IF ( .NOT.DSC(I) ) K_RAD_SMLT = NTML(I)+2
!         !-------------------------------------------------------------
!         ! Find the layer with the greatest radiative flux jump below
!         ! K_RAD_SMLT and assume that this is the top of the SML.  
!         !-------------------------------------------------------------
          IF (DF_OVER_CP(I,K) .GT. DF_TOP_OVER_CP(I)
     &                .AND. K .LE. K_RAD_SMLT ) THEN
            K_CLOUD_TOP(I) = K
            DF_TOP_OVER_CP(I) = DF_OVER_CP(I,K)
          ENDIF

        ENDDO
      ENDDO
!     !-----------------------------------------------------------------
!     !  If cloud top extends into the level above NTML (so that the
!     !  radiative divergence maximum is at NTML+1) add on RAD_HR from
!     !  NTML if there is additional cloud top cooling there.
!     !-----------------------------------------------------------------
        DO I = P1,P1-1+P_POINTS
          IF ( K_CLOUD_TOP(I) .EQ. NTML(I)+1 ) THEN
           K=NTML(I)
           DF_OVER_CP(I,K) = DELTAP(I,K)/G * RAD_HR(I,K)
           IF (DF_OVER_CP(I,K) .GT. 0.0 )
     &        DF_TOP_OVER_CP(I) = DF_TOP_OVER_CP(I) + DF_OVER_CP(I,K)
          ENDIF
        ENDDO
!
!-----------------------------------------------------------------------
!! 5.1 Subroutine EXCF_NL.
!-----------------------------------------------------------------------
!
      CALL EXCF_NL (
     & P_FIELD,P1,P_POINTS,BL_LEVELS,
     & NTML,CF,
     & RDZ,ZH,Z_UV,Z_TQ,RHO_UV,RHO_TQ,
     & Z0M,V_S,FB_SURF,DB_TOP,
     & BFLUX_SURF,BFLUX_SAT_SURF,ZETA_S,BT_TOP,BTT_TOP,                 
     & DF_TOP_OVER_CP,ZETA_R,ALPHA_R,BTC_TOP,
     & DB_TOP_CLD,CHI_S_TOP,ZC,
     & RHOKM(1,2),RHOKH(1,2),RHOKM_TOP(1,2),RHOKH_TOP(1,2),
     & ZHSC,DSCDEPTH,NTDSC,DB_DSCT,SVL,
     & BT_DSCT,BTT_DSCT,
     & DF_DSCT_OVER_CP,ZETA_R_DSC,ALPHA_R_DSC,BTC_DSCT,
     & DB_DSCT_CLD,CHI_S_DSCT,ZC_DSC,COUPLED,
     & LTIMER
     &)
!
!-----------------------------------------------------------------------
!! 5.2 Diagnose boundary layer type.
!      Six different types are considered:
!      1 - Stable b.l.
!      2 - Stratocumulus over a stable surface layer.
!      3 - Well mixed b.l. (possibly with stratocumulus).
!      4 - Decoupled stratocumulus (not over cumulus).
!      5 - Decoupled stratocumulus over cumulus.
!      6 - Cumulus capped b.l.
!-----------------------------------------------------------------------
!      First initialise the type variables and set the diagnostic ZHT.
!     
      DO I = P1,P1-1+P_POINTS
        BL_TYPE_1(I) = 0.0
        BL_TYPE_2(I) = 0.0
        BL_TYPE_3(I) = 0.0
        BL_TYPE_4(I) = 0.0
        BL_TYPE_5(I) = 0.0
        BL_TYPE_6(I) = 0.0
        ZHT(I) = MAX( ZH(I) , ZHSC(I) )
      ENDDO
      DO I = P1,P1-1+P_POINTS
        IF (.NOT.UNSTABLE(I) .AND. .NOT.DSC(I) .AND.
     &       .NOT.CUMULUS(I)) THEN
!         Stable b.l.
          BL_TYPE_1(I) = 1.0
        ELSEIF (.NOT.UNSTABLE(I) .AND. DSC(I) .AND.
     &       .NOT.CUMULUS(I)) THEN
!         Stratocumulus over a stable surface layer
          BL_TYPE_2(I) = 1.0
        ELSEIF (UNSTABLE(I) .AND. .NOT.CUMULUS(I) .AND.
     &       (.NOT.DSC(I) .OR. COUPLED(I))) THEN
!         Well mixed b.l. (possibly with stratocumulus)
          BL_TYPE_3(I) = 1.0
        ELSEIF (UNSTABLE(I) .AND. DSC(I) .AND. .NOT.CUMULUS(I)) THEN
!         Decoupled stratocumulus (not over cumulus)
          BL_TYPE_4(I) = 1.0
        ELSEIF (DSC(I) .AND. CUMULUS(I)) THEN
!         Decoupled stratocumulus over cumulus
          BL_TYPE_5(I) = 1.0
        ELSEIF (.NOT.DSC(I) .AND. CUMULUS(I)) THEN
!         Cumulus capped b.l.
          BL_TYPE_6(I) = 1.0
        ENDIF
      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('KMKHZ   ',4)
      ENDIF

      RETURN
      END
