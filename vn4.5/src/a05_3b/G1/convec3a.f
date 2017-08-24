C (c) CROWN COPYRIGHT 1995, METEOROLOGICAL OFFICE, All Rights Reserved.
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
CLL  SUBROUTINE CONVECT------------------------------------------------
CLL
CLL  PURPOSE : TOP LEVEL OF THE MASS FLUX CONVECTION SCHEME.
CLL            LOOPS ROUND MODEL LEVELS FORM SURFACE UPWARDS
CLL            A STABILITY TEST IS CARRIED OUT TO DETERMINE WHICH
CLL            POINTS ARE TOO STABLE FOR CONVECTION TO OCCUR
CLL            SUBROUTINE LIFTP AND CONVC2 ARE CALLED TO CALCULATE
CLL            THE PARCEL ASCENT
CLL            SUBROUTINE POUR IS CALLED TO CALCULATE THE EVAPORATION
CLL            OF FALLING PRECIPITATION
CLL            SUBROUTINE DD_CALL CALLS THE DOWNDRAUGHT CODE
CLL            SUBROUTINE CORNRG IS CALLED TO CONSERVE MOIST STATIC
CLL            ENERGY ONCE OTHER CALCULATIONS ARE COMPLETE
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  4.0   5/05/95  : New version (based on 2B) incorporating;
CLL                   Tracer transports
CLL                   Convective momentum transports with cloud pressure
CLL                   gradients and eddy flux formulation
CLL                   CAPE closure and CAPE diagnostic
CLL                   Diagnosis of deep/shallow/mid convection
CLL                   Pressure dependency of evaporation of
CLL                   precipitation
CLL   4.1  10/6/96  : Changed dimensions of momentum arrays to
CLL                   allow convection to be split into segments
CLL                   with the momentum transport scheme.
CLL                                    Pete Inness
CLL  4.1   25/03/96 : CAPE closure restructured to avoid increments
CLL                   from split final detrainment being included
CLL                   in the scheme if convection re-initiates
CLL                   from a level at which split final detrainment
CLL                   has already occurred.
CLL                                        P. Inness.
CLL
CLL  4.1   10/05/96   Include check to prevent negative tracer values
CLL                   (involves use or CRAY-specific MINVAL function)
CLL                                           M. Woodage,  D. Roberts
CLL   4.2    Oct. 96  T3E migration: *DEF CRAY removed
CLL                   (was used to switch on WHENIMD & MINVAL) 
CLL                                    S.J.Swarbrick
CLL
CLL  4.2   26/9/96  : Four new diagnostics added  -
CLL                   (i)  Gridbox mean conv. cloud water
CLL                   (ii) Gridbox mean conv. cloud liquid water path
CLL                   (iii)Cloud base pressure weighted by convective
CLL                        cloud amount (CCA)
CLL                   (iv) Cloud top pressure weighted by CCA
CLL                                                          J.Cairns
CLL  4.3  Feb. 97   T3E optimisation: introduce recip_pstar,
CLL                   eliminate copying into workspace arrays
CLL                   for CORENG call.         S.J.Swarbrick
CLL  4.3   03/02/97  (i) Pass logical switch L_XSCOMP down to ENVIRON.
CLL                  (ii) Put setting of bottom model layer parcel
CLL                       excess to standard deviation of turbulent
CLL                       fluctuations under control of logical
CLL                       switch L_SDXS (also passed down to ENVIRON).
CLL                                                   R.N.B.Smith
!LL  4.4  Oct 97    Add halo mask to stop redundant calculations
!LL                                               Alan Dickinson
CLL  4.4   29/08/97  Pass switch L_CCW down to CLOUD_W to determine if
CLL                  precip is included in water path and pass in switch
CLL                  L_3D_CCA to determine if a 3D conv cloud amount
CLL                  should be calculated in new subroutine CALC_3D_CCA.
CLL  4.4  26/11/97  Levels loop for DTRABYDT should be NLEV
CLL                 not TRLEV.  RTHBarnes.
!LL  4.5  5/6/98    Updraught factor and L_CLOUD_DEEP passed into 
!LL                 convection as part of anvil scheme. J.Gregory 
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!LL  4.5  20/02/98  Remove redundant code. A. Dickinson               
CLL  4.5  05/05/98  Add Fujitsu vectorization directives.
CLL  4.5  05/05/98  Use fortran 90 intrinsic TINY instead of
CLL                 1.0E-100 for safety_margin. RBarnes@ecmwf.int
CLL
CLL                                                                     
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  Dated 05/02/92
CLL
CLL LOGICAL COMPONENTS INCLUDED:
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE CONVECT(NP_FIELD,NPNTS,NLEV,NBL,TH,Q,PSTAR,BLAND,U,V,
     *                   TRACER,DTHBYDT,DQBYDT,DUBYDT,DVBYDT,RAIN,SNOW,
     *                   CCA,ICCB,ICCT,CCLWP,CCW,ICCBPxCCA,ICCTPxCCA,
     *                   GBMCCWP,GBMCCW,LCBASE,LCTOP,LCCA,
     *                   LCCLWP,CAPE_OUT,EXNER,AK,BK,
     *                   AKM12,BKM12,DELAK,DELBK,TIMESTEP,T1_SD,Q1_SD,
     &                   L_MOM,L_TRACER,L_CAPE,NTRA,TRLEV,L_XSCOMP,
     &                   L_SDXS,N_CCA_LEV,L_3D_CCA,L_CCW,MPARWTR
     &                  ,ANVIL_FACTOR ,TOWER_FACTOR
     &     ,l_halo
     &                   ,UD_FACTOR,L_CLOUD_DEEP
     &                   ,UP_FLUX,FLG_UP_FLX,DWN_FLUX,FLG_DWN_FLX,
     &                    ENTRAIN_UP,FLG_ENTR_UP,DETRAIN_UP,
     &                    FLG_DETR_UP,ENTRAIN_DWN,FLG_ENTR_DWN,
     &                    DETRAIN_DWN,FLG_DETR_DWN
     &                  )
!
      IMPLICIT NONE
C
C
C--------------------------------------------------------------------
C MODEL CONSTANTS
C--------------------------------------------------------------------
C
      REAL THPIXS_DEEP,    ! INITIAL EXCESS POTENTIAL TEMPERATURE (K)
     *     QPIXS_DEEP      ! AND MIXING RATIO (KG/KG) FOR DEEP
                           ! CONVECTION
C
      REAL THPIXS_SHALLOW, ! INITIAL EXCESS POTENTIAL TEMPERATURE (K)
     *     QPIXS_SHALLOW   ! AND MIXING RATIO (KG/KG) FOR SHALLOW
                           ! CONVECTION
C
      REAL THPIXS_MID,     ! INITIAL EXCESS POTENTIAL TEMPERATURE (K)
     *     QPIXS_MID       ! AND MIXING RATIO (KG/KG) FOR MID-LEVEL
                           ! CONVECTION
C
      PARAMETER (THPIXS_DEEP= 0.2, QPIXS_DEEP =0.0)
      PARAMETER (THPIXS_SHALLOW = 0.2, QPIXS_SHALLOW =0.0)
      PARAMETER (THPIXS_MID= 0.2, QPIXS_MID =0.0)
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
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

      REAL XSBMIN !  MINIMUM EXCESS BUOYANCY TO CONTINUE PARCEL ASCENT
                  !  (K)
      PARAMETER (XSBMIN = 0.2)
C
      REAL MPARB  !  MINIMUM (PARCEL BUOYANCY/LAYER THICKNESS) (K/PA)
      PARAMETER (MPARB = 1.0)
C
      REAL DELTHST  !  DIFFERENCE IN POTENTIAL TEMPERATURE BETWEEN
                    !  LEVELS ABOVE WHICH THE ATMOSPHERE IF ASSUMED
                    !  TO BE TOO STABLE TO CONVECT (K)
      PARAMETER (DELTHST = 1.5)
C
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
      REAL C_DEEP,    ! CONSTANTS USED TO DETERMINE INITIAL CONVECTIVE
     *     D_DEEP     ! MASS FLUX FROM PARCEL BUOYANCY FOR DEEP
                      ! CONVECTION
C
      REAL C_SHALLOW, ! CONSTANTS USED TO DETERMINE INITIAL CONVECTIVE
     *     D_SHALLOW  ! MASS FLUX FROM PARCEL BUOYANCY FOR SHALLOW
                      ! CONVECTION
C
      REAL C_MID,     ! CONSTANTS USED TO DETERMINE INITIAL CONVECTIVE
     *     D_MID      ! MASS FLUX FROM PARCEL BUOYANCY FOR MID-LEVEL
                      ! CONVECTION
C
      PARAMETER (C_DEEP = 5.17E-4, D_DEEP = 0.0)
      PARAMETER (C_SHALLOW = 5.17E-4, D_SHALLOW = 0.0)
      PARAMETER (C_MID = 5.17E-4, D_MID = 0.0)
      REAL AE1,AE2,       ! COEFFICIENTS USED IN CALCULATION
     *     ENTCOEF,       ! OF ENTRAINMENT RATE
     *     SH_FAC
C
      PARAMETER(AE1=1.0,AE2=1.5)
      PARAMETER (ENTCOEF = 3.0,SH_FAC=1.0)
C
      REAL CAPE_TS      !  TIMESCALE FOR DESTRUCTION OF CONVECTIVE
                        !  AVAILABLE POTENTIAL ENERGY BY CONVECTION
                        !  WHEN A CAPE CLOSURE TO THE CONVECTION
                        !  SCHEME IS EMPLOYED (S)
      PARAMETER ( CAPE_TS = 7200.0 )
      REAL QSTICE   !  APPROXIMATION TO SATURATION MIXING RATIO
                    !  AT TEMPERATURE AT WHICH LIQUID WATER TURNS TO
                    !  ICE (SEE COMDECK TICE) (KG/KG)
      PARAMETER (QSTICE = 3.5E-3)
C
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C
C---------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C---------------------------------------------------------------------
C
      INTEGER NP_FIELD            ! LENGTH OF DATA (ALSO USED TO
                                  ! SPECIFY STARTING POINT OF
                                  ! DATA PASSED IN)
C
      INTEGER NPNTS               ! IN FULL VECTOR LENGTH
C
      INTEGER NLEV                ! IN NUMBER OF MODEL LAYERS
C
      INTEGER NBL                 ! IN NUMBER OF BOUNDARY LAYER LEVELS
C
      INTEGER NCONV               ! NUMBER OF POINTS WHICH PASS
                                  ! INITIAL STABILITY TEST IN LAYER K
C
      INTEGER NINIT               ! NUMBER OF POINTS AT WHICH
                                  ! CONVECTION OCCURS IN LAYER K
C
      INTEGER NTERM               ! NUMBER OF CONVECTING POINTS IN
                                  ! LAYER K AT WHICH CONVECTION IS
                                  ! TERMINATING
C
      INTEGER NCNLV               ! NUMBER OF POINTS AT WHICH CONVECTION
                                  ! OCCURS AT SOME LAYER OF THE DOMAIN
C
      INTEGER NTRA                ! NUMBER OF TRACER FIELDS
C
      INTEGER TRLEV               ! NUMBER OF MODEL LEVELS ON WHICH
                                  ! TRACERS ARE INCLUDED
C
      INTEGER I,K,KC,KTRA,K_TEST, ! LOOP COUNTERS
     *        KT
C
      INTEGER N_CCA_LEV           ! Number of levels for conv cloud
!                                 ! amount: 1 for 2D, nlevs for 3D.
C
C---------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C---------------------------------------------------------------------
C
      LOGICAL BLAND(NP_FIELD)     ! IN LAND/SEA MASK
C
      LOGICAL L_TRACER            ! IN SWITCH FOR INCLUSION OF TRACERS
C
      LOGICAL L_MOM               ! IN SWITCH FOR INCLUSION OF
                                  !    MOMENTUM TRANSPORTS
C
      LOGICAL L_CAPE              ! IN SWITCH FOR USE OF CAPE CLOSURE
C
      LOGICAL L_XSCOMP            ! IN Switch for allowing compensating
                                  !    cooling and drying of the
                                  !    environment in initiating layer
C
      LOGICAL L_SDXS              ! IN Switch for allowing parcel excess
                                  !    to be set to s.d. of turbulent
                                  !    fluctuations in lowest model
                                  !    layer
C
      LOGICAL L_3D_CCA            ! IN Switch for conv cld amt varying
!                                 !    with height (3D), or not (2D)
      LOGICAL L_CCW               ! IN Switch for allowing precip
                                  !    before calculation of water
                                  !    path.
!
      LOGICAL L_CLOUD_DEEP        ! IN Switch for depth criterion for
!                                 !    anvil clouds.
!
      REAL PSTAR(NP_FIELD)        ! IN SURFACE PRESSURE (PA)
C
      REAL EXNER(NP_FIELD,NLEV+1) ! IN EXNER RATIO
C
      REAL AK(NLEV),              ! IN HYBRID CO-ORDINATE COEFFICIENTS
     *     BK(NLEV)               !    DEFINE PRESSURE AT MID-POINT
                                  !    OF LAYER K
C
      REAL AKM12(NLEV+1),         ! IN HYBRID CO-ORDINATE COEFFICIENTS
     *     BKM12(NLEV+1)          !    TO DEFINE PRESSURE AT
                                  !    LEVEL K-1/2
C
      REAL DELAK(NLEV),           ! IN DIFFERENCE IN HYBRID CO-ORDINATE
     *     DELBK(NLEV)            !    COEFFICIENTS ACROSS LAYER K
C
      REAL TIMESTEP               ! IN MODEL TIMESTEP (SECS)
C
      REAL T1_SD(NP_FIELD)        ! IN Standard deviation of turbulent
C                                 !    fluctuations of layer 1
C                                 !    temperature (K).
      REAL Q1_SD(NP_FIELD)        ! IN Standard deviation of turbulent
C                                 !    fluctuations of layer 1
C                                 !    humidity (kg/kg).
      REAL MPARWTR                ! IN Reservoir of conv cld water left
!                                 !    in a layer after conv. precip.
      REAL ANVIL_FACTOR           ! IN used in calculation of cld. amt.
     &    ,TOWER_FACTOR           !    on model levels if L_3D_CCA = .T.
!
      REAL UD_FACTOR              ! IN Updraught factor: used in conv.
!                                 !    cloud water path as seen by rad.
!                                 !    if L_CCW is true.
!
      LOGICAL l_halo(NP_FIELD)  ! Mask for halos
      LOGICAL FLG_UP_FLX          ! STASH FLAG FOR UPDRAUGHT MASS FLUX
C
      LOGICAL FLG_DWN_FLX         ! STASH FLAG FOR DOWNDRAGHT MASS FLUX
!
      LOGICAL FLG_ENTR_UP         ! STASH FLAG FOR UPDRAUGHT ENTRAINMENT
!
      LOGICAL FLG_ENTR_DWN        ! STASH FLAG FOR DOWNDRAUGHT ENTRAINMN
!
      LOGICAL FLG_DETR_UP         ! STASH FLAG FOR UPDRAUGHT DETRAINMENT
!
      LOGICAL FLG_DETR_DWN        ! STASH FLAG FOR DOWNDRAUGHT DETRAINMN
!                                                                       
C---------------------------------------------------------------------
C  VARIABLES WHICH ARE INPUT AND OUTPUT
C---------------------------------------------------------------------
C
      REAL TH(NP_FIELD,NLEV)      ! INOUT
                                  ! IN MODEL POTENTIAL TEMPERATURE (K)
                                  ! OUT MODEL POTENTIAL TEMPERATURE
                                  !     AFTER CONVECTION (K)
C
      REAL Q(NP_FIELD,NLEV)       ! INOUT
                                  ! IN MODEL MIXING RATIO (KG/KG)
                                  ! OUT MODEL MIXING RATIO AFTER
                                  !     AFTER CONVECTION (KG/KG)
C
      REAL U(NP_FIELD,NLEV)       ! INOUT
                                  ! IN MODEL U FIELD (M/S)
                                  ! OUT MODEL U FIELD AFTER CONVECTIVE
                                  !     MOMENTUM TRANSPORT (M/S)
C
      REAL V(NP_FIELD,NLEV)       ! INOUT
                                  ! IN MODEL V FIELD (M/S)
                                  ! OUT MODEL V FIELD AFTER CONVECTIVE
                                  !     MOMENTUM TRANSPORT (M/S)
C
      REAL TRACER(NP_FIELD,TRLEV, ! INOUT
     *            NTRA)           ! IN  MODEL TRACER FIELDS (KG/KG)
                                  ! OUT MODEL TRACER FIELDS AFTER
                                  !     CONVECTION (KG/KG)
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL DTHBYDT(NP_FIELD,NLEV) ! OUT INCREMENTS TO POTENTIAL
                                  !     TEMPERATURE DUE TO CONVECTION
                                  !     (K/S)
C
      REAL DQBYDT(NP_FIELD,NLEV)  ! OUT INCREMENTS TO MIXING RATIO
                                  !     DUE TO CONVECTION (KG/KG/S)
C
      REAL DUBYDT(NP_FIELD,NLEV)  ! OUT INCREMENTS TO U DUE TO
                                  !     CONVECTIVE MOMENTUM TRANSPORT
                                  !     (M/S**2)
C
      REAL DVBYDT(NP_FIELD,NLEV)  ! OUT INCREMENTS TO V DUE TO
                                  !     CONVECTIVE MOMENTUM TRANSPORT
                                  !     (M/S**2)
C
      REAL RAIN(NP_FIELD)         ! OUT SURFACE CONVECTIVE RAINFALL
                                  !     (KG/M**2/S)
C
      REAL SNOW(NP_FIELD)         ! OUT SURFACE CONVECTIVE SNOWFALL
                                  !     (KG/M**2/S)
C
      REAL CCA(NP_FIELD,N_CCA_LEV)! OUT CONVECTIVE CLOUD AMOUNT (%)
C
      INTEGER ICCB(NP_FIELD)      ! OUT CONVECTIVE CLOUD BASE LEVEL
C
      INTEGER ICCT(NP_FIELD)      ! OUT CONVECTIVE CLOUD TOP LEVEL
C
      REAL CCLWP(NP_FIELD)        ! OUT CONDENSED WATER PATH (KG/M**2)
C
      REAL CCW(NP_FIELD,NLEV)     ! OUT CONVECTIVE CLOUD LIQUID WATER
                                  ! (G/KG) ON MODEL LEVELS
C
      REAL ICCBPxCCA(NP_FIELD)    ! OUT CONV. CLD BASE PRESSURE x CCA
C
      REAL ICCTPxCCA(NP_FIELD)    ! OUT CONV. CLD TOP PRESSURE x CCA
C
      REAL GBMCCWP(NP_FIELD)      ! OUT GRIDBOX MEAN CCWP
C
      REAL GBMCCW(NP_FIELD,NLEV)  ! OUT GRIDBOX MEAN CCW
C
      REAL LCCA(NP_FIELD)         ! OUT LOWEST CONV.CLOUD AMOUNT (%)
C
      INTEGER LCBASE(NP_FIELD)    ! OUT LOWEST CONV.CLOUD BASE LEVEL
C
      INTEGER LCTOP(NP_FIELD)     ! OUT LOWEST CONV.CLOUD TOP LEVEL
C
      REAL LCCLWP(NP_FIELD)       ! OUT CONDENSED WATER PATH (KG/M**2)
                                  !     FOR LOWEST CONV.CLOUD
C
      REAL CAPE_OUT(NPNTS)        ! OUT SAVED VALUES OF CONVECTIVE
                                  !     AVAILABLE POTENTIAL ENERGY
                                  !     FOR DIAGNOSTIC OUTPUT
      REAL UP_FLUX(NP_FIELD,NLEV)     ! OUT UPDRAUGHT MASS FLUX
C
      REAL DWN_FLUX(NP_FIELD,NLEV)  ! OUT DOWNDRAUGHT MASS FLUX
!
      REAL ENTRAIN_UP(NP_FIELD,NLEV)  ! FRACTIONAL ENTRAINMENT RATE
                                      ! INTO UPDRAUGHTS
      REAL DETRAIN_UP(NP_FIELD,NLEV)  ! FRACTIONAL DETRAINMENT RATE 
                                      ! FROM UPDRAUGHTS
      REAL ENTRAIN_DWN(NP_FIELD,NLEV) ! FRACTIONAL ENTRAINMENT RATE
                                      ! INTO DOWNDRAUGHTS
      REAL DETRAIN_DWN(NP_FIELD,NLEV) ! FRACTIONAL DETRAINMENT RATE
                                      ! FROM DOWNDRAUGHTS      
                
                                                                 
C----------------------------------------------------------------------
C VARIABLES DEFINED LOCALLY
C
      REAL WORK(NPNTS,NLEV*2),    !  WORK SPACE
     *     WORK2(NPNTS,NLEV*2)
      LOGICAL BWORK(NPNTS,4),     ! WORK SPACE FOR 'BIT' MASKS
     *        BWORK2(NPNTS,4)
C
      REAL CAPE(NPNTS)            ! CONVECTIVE AVAILABLE POTENTIAL
                                  ! ENERGY (J/KG)
C
      REAL DCPBYDT(NPNTS)         ! RATE OF CHANGE OF CAPE
C
      REAL CAPE_C(NPNTS)          ! CAPE - COMPRESSED
C
      REAL DCPBYDT_C(NPNTS)       ! RATE OF CHANGE OF CAPE - COMPRESSED
C
      REAL DTHEF(NPNTS)           ! THETA INCREMENT FROM CONVECTION
                                  ! IN MODEL LEVEL AT WHICH SPLIT
                                  ! FINAL DETRAINMENT LAST OCCURRED
                                  ! (K/S)
C
      REAL DQF(NPNTS)             ! SPECIFIC HUMIDITY INCREMENT FROM
                                  ! CONVECTION IN MODEL LEVEL AT WHICH
                                  ! SPLIT FINAL DETRAINMENT LAST
                                  ! OCCURRED (KG/KG/S)
C
      REAL DUEF(NPNTS)            ! AS DTHEF BUT FOR U INCREMENTS (ms-2)
!
      REAL DVEF(NPNTS)            ! AS DTHEF BUT FOR V INCREMENTS (ms-2)
!                                                                       
      LOGICAL BCONV(NPNTS)        ! MASK FOR POINTS WHERE STABILITY
                                  ! LOW ENOUGH FOR CONVECTION
                                  ! TO OCCUR
C
      REAL QSE(NPNTS,NLEV)        ! SATURATION MIXING RATIO OF CLOUD
                                  ! ENVIRONMENT (KG/KG)
C
      REAL TT(NPNTS)              ! TEMPORARY STORE FOR TEMPERATURE
                                  ! IN CALCULATION OF SATURATION
                                  ! MIXING RATIO (K)
C
      REAL TTKM1(NPNTS)           ! TEMPORARY STORE FOR TEMPERATURE
                                  ! IN LAYER K-1 FOR USE IN FREEZING
                                  ! LEV. CALCULATION FOR ANVIL (K)
C
      REAL PT(NPNTS)              ! TEMPORARY STORE FOR PRESSURE
                                  ! IN CALCULATION OF SATURATION
                                  ! MIXING RATIO (PA)
C
      REAL CCA_2DC(NPNTS)         ! COMPRESSED VALUES OF 2D CCA
C
      INTEGER ICCBC(NPNTS)        ! COMPRESSED VALUES OF CCB
C
      INTEGER ICCTC(NPNTS)        ! COMPRESSED VALUES OF CCT
C
      REAL TCW(NPNTS)             ! TOTAL CONDENSED WATER (KG/M**2/S)
C
      REAL TCWC(NPNTS)            ! COMPRESSED VALUES OF TCW
C
      REAL CCLWPC(NPNTS)          ! COMPRESSED VALUE OF CCLWP
C
      REAL LCCAC(NPNTS)           ! COMPRESSED VALUES OF LCCA
C
      INTEGER LCBASEC(NPNTS)      ! COMPRESSED VALUES OF LCBASE
C
      INTEGER LCTOPC(NPNTS)       ! COMPRESSED VALUES OF LCTOP
C
      REAL LCCLWPC(NPNTS)         ! COMPRESSED VALUE OF LCCLWP
C
      REAL DQSTHK(NPNTS)          ! GRADIENT OF SATURATION MIXING
                                  ! RATIO OF CLOUD ENVIRONMENT WITH
                                  ! POTENTIAL TEMPERATURE IN LAYER K
                                  ! (KG/KG/K)
C
      REAL DQSTHKP1(NPNTS)        ! GRADIENT OF SATURATION MIXING
                                  ! RATIO OF CLOUD ENVIRONMENT WITH
                                  ! POTENTIAL TEMPERATURE IN LAYER K+1
                                  ! (KG/KG/K)
C
      REAL DTRABYDT(NPNTS,NLEV,   ! INCREMENT TO TRACER DUE TO
     *              NTRA)         ! CONVECTION (KG/KG/S)
C
      REAL PRECIP(NPNTS,NLEV)     ! AMOUNT OF PRECIPITATION
                                  ! FROM EACH LAYER (KG/M*:2/S)
C
      REAL THPI(NPNTS)            ! INITIAL PARCEL POTENTIAL TEMPERATURE
                                  ! (K)
C
      REAL QPI(NPNTS)             ! INITIAL PARCEL MIXING RATIO
                                  ! (KG/KG)
C
      REAL TRAPI(NPNTS,NTRA)      ! INITIAL PARCEL TRACER CONTENT
                                  ! (KG/KG)
C
      REAL THP(NPNTS,NLEV)        ! PARCEL POTENTIAL TEMPERATURE
                                  ! IN LAYER K (K)
C
      REAL QP(NPNTS,NLEV)         ! PARCEL MIXING RATIO IN LAYER K
                                  ! (KG/KG)
C
      REAL UP(NPNTS,NLEV)         ! PARCEL U IN LAYER K (M/S)
C
      REAL VP(NPNTS,NLEV)         ! PARCEL V IN LAYER K (M/S)
C
      REAL TRAP(NPNTS,NLEV,NTRA)  ! PARCEL TRACER CONTENT IN LAYER K
                                  ! (KG/KG)
C
      REAL XPK(NPNTS,NLEV)        ! PARCEL CLOUD WATER IN LAYER K
                                  ! (KG/KG)
C
      REAL FLX(NPNTS,NLEV)        ! PARCEL MASSFLUX IN LAYER K (PA/S)
C
      REAL FLX_INIT(NPNTS)        ! INITIAL MASSFLUX AT CLOUD BASE
                                  ! (PA/S)
C
      REAL FLX_INIT_NEW(NPNTS)    ! INITIAL MASSFLUX AT CLOUD BASE,
                                  ! SCALED TO DESTROY CAPE OVER
                                  ! GIVEN TIMESCALE (PA/S)
C
      REAL FLXMAX_INIT(NPNTS)     ! MAXIMUM POSSIBLE INITIAL MASSFLUX
                                  ! LIMITED TO THE MASS IN TH INITIAL
                                  ! CONVECTING LAYER (PA/S)
C
      INTEGER START_LEV(NPNTS)    ! LEVEL AT WHICH CONVECTION INITIATES
C
      INTEGER DET_LEV(NPNTS)      ! LEVEL AT WHICH SPLIT FINAL
                                  ! DETRAINMENT LAST OCCURRED
C
      LOGICAL BINIT(NPNTS)        ! MASK FOR POINTS WHERE CONVECTION
                                  ! IS OCCURING
C
      LOGICAL BTERM(NPNTS)        ! MASK FOR POINTS WHERE CONVECTION
                                  ! TERMINATES IN LAYER K+1
C
      LOGICAL BWATER(NPNTS,2:NLEV) ! MASK FOR POINTS AT WHICH
                                   ! PRECIPITATION IS LIQUID
C
      LOGICAL BGMK(NPNTS)         ! MASK FOR POINTS WHERE PARCEL IN
                                  ! LAYER K IS SATURATED
C
      LOGICAL BCNLV(NPNTS)        ! MASK FOR THOSE POINTS AT WHICH
                                  ! CONVECTION HAS OCCURED AT SOME
                                  ! LEVEL OF THE MODEL
C
      REAL DEPTH(NPNTS)           ! DEPTH OF CONVECTIVE CLOUD (M)
C
      REAL FLXMAXK(NPNTS)         ! MAXIMUM INITIL CONVECTIVE MASSFLUX
                                  ! (PA/S)
C
      REAL FLXMAX2(NPNTS)         ! MAXIMUM INITIL CONVECTIVE MASSFLUX
                                  ! (PA/S)
C
      REAL PK(NPNTS)              ! PRESSURE AT MID-POINT OF LAYER K
                                  ! (PA)
C
      REAL PKP1(NPNTS)            ! PRESSURE AT MID-POINT OF LAYER K+1
                                  ! (PA)
C
      REAL DELPK(NPNTS)           ! PRESSURE DIFFERENCE ACROSS LAYER K
                                  ! (PA)
C
      REAL DELPKP1(NPNTS)         ! PRESSURE DIFFERENCE ACROSS LAYER K+1
                                  ! (PA)
C
      REAL DELPKP12(NPNTS)        ! PRESSURE DIFFERENCE BETWEEN
                                  ! LEVELS K AND K+1 (PA)
C
      REAL EKP14(NPNTS),          ! ENTRAINMENT COEFFICIENTS AT LEVELS
     *     EKP34(NPNTS)           ! K+1/4 AND K+3/4 MULTIPLIED BY
                                  ! APPROPRIATE LAYER THICKNESS
C
      REAL AMDETK(NPNTS)          ! MIXING DETRAINMENT COEFFICIENT AT
                                  ! LEVEL K MULTIPLIED BY APPROPRIATE
                                  ! LAYER THICKNESS
C
      REAL DELTAK(NPNTS)          ! FORCED DETRAINMENT RATE
C                                                                       
      REAL EXK(NPNTS)             ! EXNER RATIO AT LEVEL K
C
      REAL EXKP1(NPNTS)           ! EXNER RATIO AT LEVEL K+1
C
      REAL DELEXKP1(NPNTS)        ! DIFFERENCE IN EXNER RATIO
                                  ! ACROSS LAYER K+1
C
      REAL EMINDS(NPNTS)          ! MINIMUM BUOYANCY FOR CONVECTION TO
                                  ! INITIATE FROM LAYER K
C
      INTEGER INDEX1(NPNTS),      ! INDEX FOR COMPRESS AND
     *        INDEX2(NPNTS),      ! EXPAND
     *        INDEX3(NPNTS),
     *        INDEX4(NPNTS)
C
      LOGICAL L_SHALLOW(NPNTS)    ! CONVECTION LIKELY TO BE SHALLOW
                                  ! IF SET TO TR
C
      LOGICAL L_SHALLOW_C(NPNTS), ! CONVECTION LIKELY TO BE SHALLOW
     *        L_SHALLOW_C2(NPNTS) ! IF SET TO TRUE -- COMPRESSED
C
      LOGICAL L_MID(NPNTS)        ! CONVECTION STARTS ABOVE BOUNDARY
                                  ! LAYER IF SET TO TRUE
C
      LOGICAL L_MID_C(NPNTS),     ! CONVECTION STARTS ABOVE BOUNDARY
     *        L_MID_C2(NPNTS)     ! LAYER IF SET TO TRUE -- COMPRESSED
C
      REAL TRAPK_C(NPNTS,NTRA),   ! PARCEL TRACER CONTENT IN LAYER K
     *     TRAPK_C2(NPNTS,NTRA)   ! - COMPRESSED (KG/KG)
C
      REAL TRAPKP1_C(NPNTS,NTRA), ! PARCEL TRACER CONTENT IN LAYER K+1
     *     TRAPKP1_C2(NPNTS,NTRA) ! - COMPRESSED (KG/KG)
C
      REAL TRAEK_C(NPNTS,NTRA),   ! TRACER CONTENT OF CLOUD ENVIRONMENT
     *     TRAEK_C2(NPNTS,NTRA)   ! IN LAYER K - COMPRESSED (KG/KG)
C
      REAL TRAEKP1_C(NPNTS,NTRA), ! TRACER CONTENT OF CLOUD ENVIRONMENT
     *     TRAEKP1_C2(NPNTS,NTRA) ! IN LAYER K+1 - COMPRESSED (KG/KG)
C
      REAL DTRAEK_C(NPNTS,NTRA)   ! INCREMENTS TO MODEL TRACER
                                  ! DUE TO CONVECTION AT LEVEL K
                                  ! - COMPRESSED (KG/KG/S)
C
      REAL DTRAEKP1_C(NPNTS,NTRA) ! INCREMENTS TO MODEL TRACER DUE TO
                                  ! CONVECTION IN LAYER K+1 -COMPRESSED
                                  ! (KG/KG/S)
C
      REAL EFLUX_U_UD(NPNTS),     ! VERTICAL EDDY FLUX OF MOMENTUM DUE
     *     EFLUX_V_UD(NPNTS)      ! TO UD AT TOP OF A LAYER
C
      REAL EFLUX_U_DD(NPNTS),     ! VERTICAL EDDY FLUX OF MOMENTUM DUE
     *     EFLUX_V_DD(NPNTS)      ! TO DD AT BOTTOM OF A LAYER
C
      REAL LIMITED_STEP(NPNTS),   ! Reduced step size for tracer mixing
     &     STEP_TEST1(NLEV),      ! Work array used in reducing step
     &     STEP_TEST2(NLEV)       !       "
      REAL REDUCTION_FACTOR(NPNTS,NTRA)    ! Diagnostic array for time-
!                                 ! step reduction factor for tracers
      REAL SAFETY_MARGIN          ! Small no. used in tracer step reducn
C
      PARAMETER (SAFETY_MARGIN = 1.0E-100 )
!
C
      INTEGER FREEZE_LEV(NPNTS)   ! FREEZING LEVEL
C
      REAL CCA_2D(NPNTS)          ! Conv cloud amount on a single
!                                 ! level, as calculated in CONRAD
C
C
      REAL FLX2                   ! TEMPORARY STORE FOR MASS FLUX
C
      REAL AEKP14,AEKP34          ! CONSTANTS USED IN CALCULATION
                                  ! OF ENTRAINMENT COEFFICIENTS
C
      REAL EL                     ! LATENT HEAT OF CONDENSATION
                                  ! USED IN UNDILUTE ASCENT CALCULATION
C
      REAL THVUNDI,THVEKP1        ! VIRTUAL TEMPERATURE OF UNDILUTE
                                  ! PARCEL AND ENVIRONMENT USED IN
                                  ! BUOYANCY CALCULATIONS FOR THE
                                  ! UNDILUTE ASCENT
C
      REAL C,D                    ! MASS FLUX PARAMETERS
C
      REAL recip_PSTAR(NP_FIELD)  ! Reciprocal of pstar array 
C
C----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C----------------------------------------------------------------------
C
      EXTERNAL QSAT,FLAG_WET,LIFT_PAR,CONVEC2,LAYER_CN,
     *         DQS_DTH,COR_ENGY,DD_CALL,CALC_3D_CCA
C

      REAL
     &    PU,PL,PM
C*L------------------ COMDECK P_EXNERC ---------------------------------
C statement function to define exner pressure at model levels
C from exner pressures at model half-levels
      REAL P_EXNER_C
      REAL R_P_EXNER_C
      REAL P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM
C consistent with geopotential see eqn 26 DOC PAPER 10
      P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)/
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )
      R_P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )/
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)

C*------------------- --------------------------------------------------


C*---------------------------------------------------------------------
C
CL
CL---------------------------------------------------------------------
CL CALCULATE AN ARRAY OF SATURATION MIXING RATIOS
CL FIRST CONVERT POTENTIAL TEMPERATURE TO TEMPERATURE AND CALCULATE
CL PRESSURE OF LAYER K
CL
CL SUBROUTINE QSAT
CL UM DOCUMENTATION PAPER P282
CL---------------------------------------------------------------------
CL
C  Calculate reciprocal of pstar
      DO I=1,NPNTS
        RECIP_PSTAR(I)=1./PSTAR(I)
      ENDDO
C
      DO 20 K=1,NLEV
       DO 25 I = 1,NPNTS
        TTKM1(I)=TT(I)
        PU=PSTAR(I)*BKM12(K+1) + AKM12(K+1)
        PL=PSTAR(I)*BKM12(K) + AKM12(K)
        TT(I) = TH(I,K)* P_EXNER_C(EXNER(I,K+1),EXNER(I,K),PU,PL,KAPPA)
        PT(I) = AK(K)+BK(K)*PSTAR(I)
        IF (TT(I).LT.TM) THEN
          IF (K.EQ.1) THEN
            FREEZE_LEV(I)=K
          ELSEIF(TTKM1(I).GE.TM) THEN
            FREEZE_LEV(I)=K
          ENDIF
        ENDIF
   25  CONTINUE
C
       CALL QSAT (QSE(1,K),TT,PT,NPNTS)
C
  20  CONTINUE
CL
CL---------------------------------------------------------------------
CL CALCULATE BIT VECTOR WHERE WATER WILL CONDENSE RATHER THAN ICE
CL SUBROUTINE FLAG_WET
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (2B)
CL---------------------------------------------------------------------
CL
      CALL FLAG_WET(BWATER,TH,EXNER,PSTAR,AKM12,BKM12,
     &                    NP_FIELD,NPNTS,NLEV)
C
C----------------------------------------------------------------------
C INITIALISE PRECIPITATION, DTH/DT, DQ/DT, CCW
C DU/DT, DV/DT AND TRACER INCREMENT ARRAYS
C----------------------------------------------------------------------
C
      DO K=1,NLEV
       DO I=1,NPNTS
        PRECIP(I,K) = 0.0
        CCW(I,K) = 0.0
        GBMCCW(I,K) = 0.0
        DTHBYDT(I,K) = 0.0
        DQBYDT(I,K) = 0.0
        IF(L_MOM)THEN
          DUBYDT(I,K) = 0.0
          DVBYDT(I,K) = 0.0
        END IF
       END DO
      END DO
      IF(L_TRACER)THEN
      DO KTRA=1,NTRA
        DO K=1,NLEV
          DO I=1,NPNTS
            DTRABYDT(I,K,KTRA) = 0.0
          END DO
        END DO
      END DO
      END IF
      DO K=1,N_CCA_LEV
        DO I=1,NPNTS
          CCA(I,K) = 0.0
        ENDDO
      ENDDO
C
      DO 50 I=1,NPNTS
C
C----------------------------------------------------------------------
C INITIALISE BIT VECTORS FOR POINTS WHICH ARE ALREADY CONVECTING
C AND FOR POINTS AT WHICH CONVECTION OCCURS AT SOME LEVEL OF
C THE ATMOSPHERE. ALSO SET BIT VECTORS FOR SHALLOW AND MID LEVEL
C CONVECTION TO FALSE AS DEEP CONVECTION IS ASSUMED UNTIL TEST
C ASCENT IS PERFORMED.
C----------------------------------------------------------------------
C
        BINIT(I) = .FALSE.
        BCNLV(I) = .FALSE.
        BTERM(I) = .FALSE.
        L_SHALLOW(I) = .FALSE.
        L_MID(I) = .FALSE.
C
C----------------------------------------------------------------------
C INITIALISE RADIATION DIAGNOSTICS
C----------------------------------------------------------------------
C
       CCA_2D(I) = 0.0
       ICCB(I) = 0
       ICCT(I) = 0
       TCW(I) = 0.0
       CCLWP(I) = 0.0
C
C---------------------------------------------------------------------
C INITIALISE GRIDBOX MEAN DIAGNOSTICS
C---------------------------------------------------------------------
C
       GBMCCWP(I) = 0.0
       ICCBPxCCA(I) = 0.0
       ICCTPxCCA(I) = 0.0
C
CL-------------------------------------------------------------------
CL INITIALISE DIAGNOSTICS FOR CLOSURE CALCULATION
CL-------------------------------------------------------------------
C
       FLX_INIT(I) = 0.0
       FLX_INIT_NEW(I) = 0.0
       CAPE(I) = 0.0
       CAPE_OUT(I) = 0.0
       DCPBYDT(I) = 0.0
       CAPE_C(I) = 0.0
       DCPBYDT_C(I) = 0.0
       START_LEV(I) = 0
       DELTAK(I)=0.0
       DET_LEV(I) = 0
       DTHEF(I) = 0.0
       DQF(I) = 0.0
       DUEF(I) = 0.0
       DVEF(I) = 0.0
C
C---------------------------------------------------------------------
C INITIALISE EDDY FLUX ARRAYS FOR UD AND DD
C--------------------------------------------------------------------
C
       EFLUX_U_UD(I) = 0.0
       EFLUX_V_UD(I) = 0.0
       EFLUX_U_DD(I) = 0.0
       EFLUX_V_DD(I) = 0.0
C
C---------------------------------------------------------------------
C INITIALISE SURFACE PRECIPITATION ARRAYS
C---------------------------------------------------------------------
C
       RAIN(I) = 0.0
  50   SNOW(I) = 0.0
CL
CL=====================================================================
CL MAIN LOOP OVER LEVELS - FROM SURFACE TO TOP
CL=====================================================================
CL
      DO 60 K=1,NLEV-1
CL
CL---------------------------------------------------------------------
CL CALCULATE LEVEL PRESSURES, EXNER RATIO FOR MID POINTS, ENTRAINMENT
CL RATES, DETRAINMENTS RATES AND PRESSURE DIFFERENCE ACROS  LAYERS AS
CL A FUNCTION OF GRID-POINT
CL
CL SUBROUTINE LAYER_CN
CL---------------------------------------------------------------------
CL
      CALL LAYER_CN(K,NP_FIELD,NPNTS,NLEV,EXNER,AK,BK,AKM12,BKM12,
     *              DELAK,DELBK,PSTAR,PK,PKP1,DELPK,DELPKP1,
     *              DELPKP12,EKP14,EKP34,AMDETK,EXK,EXKP1,
     *              DELEXKP1,recip_PSTAR)      
CL
CL---------------------------------------------------------------------
CL CALCULATE DQS/DTH FOR LAYERS K AND K+1
CL
CL SUBROUTINE DQS_DTH
CL---------------------------------------------------------------------
CL
      IF (K.EQ.1) THEN
       CALL DQS_DTH(DQSTHK,K,TH(1,K),QSE(1,K),EXK,NPNTS)
      ELSE
       DO 65 I=1,NPNTS
        DQSTHK(I) = DQSTHKP1(I)
  65   CONTINUE
      END IF
C
       CALL DQS_DTH(DQSTHKP1,K+1,TH(1,K+1),QSE(1,K+1),EXKP1,NPNTS)
C
      DO 70 I=1,NPNTS
C
C---------------------------------------------------------------------
C SET OTHER GIRD-POINT DEPENDENT CONSTANTS
C---------------------------------------------------------------------
C
C---------------------------------------------------------------------
C MAXIMUM INITIAL CONVECTIVE MASSFLUX
C---------------------------------------------------------------------
C
       FLXMAXK(I) = DELPK(I)/((1.0 + EKP14(I)) * TIMESTEP)
C
C---------------------------------------------------------------------
C MAXIMUM CONVECTIVE MASSFLUX AT MID-POINT OF LAYER 2
C---------------------------------------------------------------------
C
      IF (K.EQ.1) FLXMAX2(I) = (PSTAR(I)-PKP1(I)) / TIMESTEP
C
C---------------------------------------------------------------------
C MINIMUM BUOYANCY FOR CONVECTION TO START FROM LAYER K
C---------------------------------------------------------------------
C
       EMINDS(I) = MPARB*DELPKP12(I)*RECIP_PSTAR(I)
C
C----------------------------------------------------------------------
C SET BIT VECTOR FOR POINTS WHERE CONVECTION HAS OCCURRED AT SOME
C LEVEL OF THE ATMOSPHERE
C-----------------------------------------------------------------------
C
       BCNLV(I) =  BCNLV(I) .OR. BINIT(I)
CL
CL---------------------------------------------------------------------
CL SET INITIAL VALUES FOR POINTS NOT ALREADY INITIATED
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (3), EQUATION(17)
CL---------------------------------------------------------------------
CL
       IF (.NOT.BINIT(I)) THEN
C
        IF (K.LT.NBL) THEN
C
C----------------------------------------------------------------------
C SET TO DEEP CONVECTIVE VALUES - MODIFIED LATER IF SHALLOW CONVECTION
C IS TO DEVELOP
C----------------------------------------------------------------------
C
         L_SHALLOW(I) = .FALSE.
         IF ( L_SDXS .AND. K .EQ. 1 ) THEN   
           THPI(I) = TH(I,K) + MAX ( THPIXS_DEEP , T1_SD(I)/EXK(I) )
           THP(I,K) = TH(I,K) + MAX ( THPIXS_DEEP , T1_SD(I)/EXK(I) )
           QPI(I) = Q(I,K) + MAX ( QPIXS_DEEP , Q1_SD(I) )
           QP(I,K) = Q(I,K) + MAX ( QPIXS_DEEP , Q1_SD(I) )
         ELSE
           THPI(I) = TH(I,K) + THPIXS_DEEP
           THP(I,K) = TH(I,K) + THPIXS_DEEP
           QPI(I) = Q(I,K) + QPIXS_DEEP
           QP(I,K) = Q(I,K) + QPIXS_DEEP
         END IF
C
        ELSE           ! IF(K.GE.NBL)
C
C----------------------------------------------------------------------
C SET TO VALUES FOR MID-LEVEL CONVECTION
C----------------------------------------------------------------------
C
         L_MID(I) = .TRUE.
         THPI(I) = TH(I,K) + THPIXS_MID
         THP(I,K) = TH(I,K) + THPIXS_MID
         QPI(I) = Q(I,K) + QPIXS_MID
         QP(I,K) = Q(I,K) + QPIXS_MID
C
        END IF         ! IF(K.LT.NBL) END
C
        XPK(I,K) = 0.0
        FLX(I,K) = 0.0
        BGMK(I) = .FALSE.
        DEPTH(I) = 0.0
C
       END IF          ! IF(.NOT.BINIT(I)) END
CL
CL----------------------------------------------------------------------
CL FORM A BIT VECTOR OF POINTS FOR WHICH CONVECTION MAY BE POSSIBLE
CL FROM LAYER K TO K+1 EITHER BECAUSE STABILITY IS LOW ENOUGH
CL OR BECAUSE CONVECTION OCCURRING FROM LAYER K+1 TO K
CL THIS BIT VECTOR IS USED IN THE FIRST COMPRESS OF THE DATA
CL TO CALCULATE PARCEL BUOYANCY IN LAYER K+1
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION(3), EQUATION(16)
CL----------------------------------------------------------------------
CL
        BCONV(I) = BINIT(I) .OR.
     *           ((TH(I,K) - TH(I,K+1) + DELTHST
     *           + MAX(0.0,(Q(I,K)-QSE(I,K+1)))*(LC/(CP*EXKP1(I))))
     *           .GT. 0.)
        BCONV(I) = l_halo(I).AND.BCONV(I)
  70  CONTINUE
C
CL----------------------------------------------------------------------
CL READ INITIAL VALUES OF MOMENTUM AND TRACER INTO THE PARCEL
CL----------------------------------------------------------------------
CL
        IF(L_MOM)THEN
         DO I=1,NPNTS
         IF(.NOT.BINIT(I))THEN
          UP(I,K)=U(I,K)
          VP(I,K) = V(I,K)
         END IF
         END DO
        END IF
C
        IF(L_TRACER)THEN
C
        DO KTRA = 1,NTRA
          DO I = 1,NPNTS
          IF(.NOT.BINIT(I))THEN
           TRAPI(I,KTRA) = TRACER(I,K,KTRA)
           TRAP(I,K,KTRA) = TRAPI(I,KTRA)
          END IF
          END DO
        END DO
C
        END IF
C
CL
CL----------------------------------------------------------------------
CL COMPRESS DOWN POINTS ON THE BASIS OF BIT VECTOR BCONV
CL----------------------------------------------------------------------
CL
      NCONV = 0
      DO 75 I=1,NPNTS
        IF(BCONV(I))THEN
          NCONV = NCONV + 1
          INDEX1(NCONV) = I
        END IF
  75  CONTINUE
C
C----------------------------------------------------------------------
C  WORK SPACE USAGE FOR FIRST COMPRESS ON BASIS OF SIMPLE
C  STABILITY TEST (SECTION (3), EQN(16))
C
C  REFERENCES TO WORK AND BWORK REFER TO STARTING ADDRESS
C
C  LENGTH OF COMPRESSES DATA = NCONV
C
C  WORK(1,1)  = TH(#,K)
C  WORK(1,2)  = TH(#,K+1)
C  WORK(1,3)  = Q(#,K)
C  WORK(1,4)  = Q(#,K+1)
C  WORK(1,5)  = QSE(#,K+1)
C  WORK(1,6)  = DQSTHKP1(#)
C  WORK(1,7)  = THP(#,K)
C  WORK(1,8)  = QP(#,K)
C  WORK(1,9)  = PKP1(#)
C  WORK(1,10) = EXKP1(#)
C  WORK(1,11) = EKP14(#)
C  WORK(1,12) = EKP34(#)
C  WORK(1,13) = PARCEL POT. TEMPERATURE IN LAYER K+1
C  WORK(1,14) = PARCEL MIXING RATIO IN LAYER K+1
C  WORK(1,15) = EXCESS WATER VAPOUR IN PARCEL ABOVE
C               SATURATION AFTER DRY ASCENT
C  WORK(1,16) = PARCEL BUOYANCY IN LAYER K+1
C  WORK(1,17) = DELPKP12(#)
C  WORK(1,18) = PSTAR(#)
C  WORK(1,19) = FLX(#,K)
C  WORK(1,20) = EMINDS(#)
C  WORK(1,21) = FLXMAXK(#)
C  WORK(1,22) = FLXMAX2(#)
C  WORK(1,23) = U(#,K)
C  WORK(1,24) = U(#,K+1)
C  WORK(1,25) = V(#,K)
C  WORK(1,26) = V(#,K+1)
C  WORK(1,27) = UP(#,K)
C  WORK(1,28) = VP(#,K)
C  WORK(1,29) = PARCEL U IN LAYER K+1
C  WORK(1,30) = PARCEL V IN LAYER K+1
C
C  BWORK(1,1) = BWATER(INDEX1(I),K+1)
C  BWORK(1,2) = .TRUE. IF PARCEL SATURATED IN LAYER K+1
C  BWORK(1,3) = .TRUE. IF CONVECTION INITIATE FROM LAYER K+1
C  BWORK(1,4) = BINIT(INDEX1(I))
C----------------------------------------------------------------------
C
      IF (NCONV .NE. 0) THEN
        DO 80 I=1,NCONV
          WORK(I,1)  = TH(INDEX1(I),K)
          WORK(I,2)  = TH(INDEX1(I),K+1)
          WORK(I,3)  = Q(INDEX1(I),K)
          WORK(I,4)  = Q(INDEX1(I),K+1)
          WORK(I,5)  = QSE(INDEX1(I),K+1)
          WORK(I,6)  = DQSTHKP1(INDEX1(I))
          WORK(I,7)  = THP(INDEX1(I),K)
          WORK(I,8)  = QP(INDEX1(I),K)
          WORK(I,9)  = PKP1(INDEX1(I))
          WORK(I,10) = EXKP1(INDEX1(I))
          WORK(I,11) = EKP14(INDEX1(I))
          WORK(I,12) = EKP34(INDEX1(I))
          WORK(I,17) = DELPKP12(INDEX1(I))
          WORK(I,18) = PSTAR(INDEX1(I))
          WORK(I,19) = FLX(INDEX1(I),K)
          WORK(I,20) = EMINDS(INDEX1(I))
          WORK(I,21) = FLXMAXK(INDEX1(I))
          WORK(I,22) = FLXMAX2(INDEX1(I))
          BWORK(I,1) = BWATER(INDEX1(I),K+1)
          BWORK(I,4) = BINIT(INDEX1(I))
          L_SHALLOW_C(I) = L_SHALLOW(INDEX1(I))
          L_MID_C(I) = L_MID(INDEX1(I))
C
  80    CONTINUE
C
        IF(L_MOM)THEN
         DO I=1,NCONV
           WORK(I,23) = U(INDEX1(I),K)
           WORK(I,24) = U(INDEX1(I),K+1)
           WORK(I,25) = V(INDEX1(I),K)
           WORK(I,26) = V(INDEX1(I),K+1)
           WORK(I,27) = UP(INDEX1(I),K)
           WORK(I,28) = VP(INDEX1(I),K)
         END DO
        END IF
C
        IF(L_TRACER)THEN
C
        DO KTRA = 1,NTRA
          DO I=1,NCONV
           TRAEK_C(I,KTRA)   = TRACER(INDEX1(I),K,KTRA)
           TRAEKP1_C(I,KTRA) = TRACER(INDEX1(I),K+1,KTRA)
           TRAPK_C(I,KTRA)   = TRAP(INDEX1(I),K,KTRA)
          END DO
        END DO
C
        END IF
C
      IF ( K.LT.NBL) THEN
C
CL
CL--------------------------------------------------------------------
CL CARRY OUT TEST ASCENT TO ASCERTAIN WHETHER DEEP CONVECTION OR
CL SHALLOW CONVECTION IS POSSIBLE.
CL
CL UM DOCUMENTATION PAPER P27-3. SECTION 2.
CL
CL CALCULATION ONLY CARRIED OUT FOR CONVECTION INITIATING WITHIN THE
CL BOUNDARY LAYER
CL--------------------------------------------------------------------
CL
       DO K_TEST=K,NBL        ! LOOP OVER BOUNDARY LAYER LEVELS
C---------------------------------------------------------------------
C  SET COEFFICIENTS FOR CALCULATION OF ENTRAINMENT RATES
C---------------------------------------------------------------------
       IF(K_TEST.EQ.1)THEN
         AEKP14 = AE1
         AEKP34 = AE2
       ELSE
         AEKP14 = AE2
         AEKP34 = AE2
       END IF
C
C--------------------------------------------------------------------
C SET VALUES FOR TEST ASCENT
C--------------------------------------------------------------------
C
        IF ( K_TEST .EQ. K ) THEN
C
         DO I=1,NCONV         ! 1ST COMPRESS LOOP
          WORK2(I,1) = WORK(I,1)          ! THEK
          WORK2(I,2) = WORK(I,2)          ! THEKP1
          WORK2(I,3) = WORK(I,3)          ! QEK
          WORK2(I,4) = WORK(I,4)          ! QEKP1
          WORK2(I,5) = WORK(I,5)          ! QSEKP1
          WORK2(I,6) = WORK(I,6)          ! DQSEKP1
          WORK2(I,7) = WORK(I,7)          ! THPK
          WORK2(I,8) = WORK(I,8)          ! QPK
          WORK2(I,9) = WORK(I,9)          ! PKP1
          WORK2(I,10) = WORK(I,10)        ! EXKP1
          WORK2(I,11) = WORK(I,11)        ! EKP14
          WORK2(I,12) = WORK(I,12)        ! EKP34
          BWORK2(I,1) = BWORK(I,1)        ! BWATER KP1
          BWORK2(I,3) = .FALSE.           ! POINT WHERE CONVECTION
                                          ! HAS INITIATED FROM LAYER K
                                          ! OR ABOVE
          WORK2(I,20) = WORK(I,20)        ! EMINDS
C
         END DO              ! END OF 1ST COMPRESS LOOP
C
C
        ELSE                 ! IF(K_TEST.NE.K)
C
         DO I=1,NCONV        ! 2ND COMPRESS LOOP
C
          WORK2(I,1) = WORK2(I,2)                          ! THEK
          WORK2(I,2) = TH(INDEX1(I),K_TEST+1)              ! THEKP1
          WORK2(I,3) = WORK2(I,4)                          ! QEK
          WORK2(I,4) = Q(INDEX1(I),K_TEST+1)               ! QEKP1
          WORK2(I,5) = QSE(INDEX1(I),K_TEST+1)             ! QSEKP1
          WORK2(I,7) = WORK2(I,13)                         ! THPK
          WORK2(I,8) = WORK2(I,14)                         ! QPK
          WORK2(I,9) = AK(K_TEST+1) + BK(K_TEST+1)
     *                 *WORK(I,18)                         ! PKP1
          PU = WORK(I,18)*BKM12(K_TEST+2)+AKM12(K_TEST+2)
          PL = WORK(I,18)*BKM12(K_TEST+1)+AKM12(K_TEST+1)
          PM = WORK(I,18)*BK(K_TEST)+AK(K_TEST)
          WORK2(I,10) = P_EXNER_C(EXNER(INDEX1(I),K_TEST+2),
     *              EXNER(INDEX1(I),K_TEST+1),PU,PL,KAPPA) ! EXKP1
          WORK2(I,11) = ENTCOEF * AEKP14 * PM *
     *                (PM-AKM12(K_TEST+1)-BKM12(K_TEST+1)*
     *                WORK(I,18))/(WORK(I,18)*WORK(I,18))  ! EKP14
          WORK2(I,12) = ENTCOEF *AEKP34 * (AKM12(K_TEST+1)
     *                +BKM12(K_TEST+1)*WORK(I,18))*
     *                (AKM12(K_TEST+1)+BKM12(K_TEST+1)*
     *                WORK(I,18)-WORK2(I,9))/(WORK(I,18)*
     *                WORK(I,18))                          ! EKP34
          WORK2(I,20) = EMINDS(INDEX1(I))                  ! EMINDS
          BWORK2(I,1) = BWATER(INDEX1(I),K_TEST+1)         ! BWATER KP1
C
         END DO              ! END OF 2ND COMPRESS LOOP
C
         CALL DQS_DTH(WORK2(1,6),K_TEST+1,WORK2(1,2),WORK2(1,5),
     *                WORK2(1,10),NCONV)
C
        END IF               ! IF(K_TEST.EQ.K) END
C
C--------------------------------------------------------------------
C CARRY OUT TEST ASCENT
C L_TRACER AND L_MOM(THE LOGICAL SWITCHES FOR INCLUSION OF TRACERS AND
C MOMENTUM ARE SET TO .FALSE. IN THIS CALL SINCE THIS ASCENT IS PURELY
C TO DIAGNOSE THE DEPTH OF THE INITIATED CONVECTION. THUS NOT
C INCLUDING THE TRACERS AND WINDS SAVES CPU TIME AND MEMORY.
C--------------------------------------------------------------------
C
      CALL LIFT_PAR (NCONV,NPNTS,WORK2(1,13),WORK2(1,14),WORK2(1,15),
     *               BWORK2(1,2),BWORK2(1,1),WORK2(1,7),WORK2(1,8),
     *               WORK2(1,2),WORK2(1,4),WORK2(1,1),WORK2(1,3),
     *               WORK2(1,5),WORK2(1,6),WORK2(1,9),WORK2(1,10),
     *               WORK2(1,11),WORK2(1,12),.FALSE.,WORK2(1,29),
     *               WORK2(1,30),WORK2(1,27),WORK2(1,28),WORK2(1,23),
     *               WORK2(1,24),WORK2(1,25),WORK2(1,26),.FALSE.,NTRA,
     *               TRAPKP1_C2,TRAPK_C2,TRAEKP1_C2,TRAEK_C2,
     *               L_SHALLOW_C)
C
! Fujitsu vectorization directive
!OCL NOVREC
      DO I=1,NCONV           ! 1ST LOOP OVER CONVECTING POINTS
CL
CL---------------------------------------------------------------------
CL CALCULATE BUOYANCY OF PARCEL IN LAYER K+1
CL---------------------------------------------------------------------
CL
        WORK2(I,16) = WORK2(I,13)*(1.0 +
     *                            C_VIRTUAL * WORK2(I,14))
     *               - WORK2(I,2)*(1.0 +
     *                            C_VIRTUAL * WORK2(I,4))
C
C----------------------------------------------------------------------
C INITIATE CONVECTION WHERE BUOYANCY IS LARGE ENOUGH
C----------------------------------------------------------------------
C
        IF ( .NOT.BWORK2(I,3) .AND. .NOT.BWORK(I,4) )
     *                               BWORK2(I,3) = WORK2(I,16) .GT.
     *                               (WORK2(I,20)+XSBMIN)
C
C----------------------------------------------------------------------
C CHECK TO SEE IF CONVECTION INITIATING BETWEEN LAYERS K AND NBL
C REACHES ZERO BUOYANCY BEFORE NBL+1
C---------------------------------------------------------------------
C
        IF ( BWORK2(I,3) .AND. .NOT.BWORK(I,4) .AND.
     *     .NOT.L_SHALLOW_C(I).AND. WORK2(I,16) .LE. 0.0) THEN
         L_SHALLOW_C(I) = .TRUE.
         L_SHALLOW(INDEX1(I)) = L_SHALLOW_C(I)
C
C----------------------------------------------------------------------
C IF IN TOP 2 LAYERS OF BOUNDARY LAYER, CALCULATE THE POTENTIAL
C TEMPERATURE OF AN UNDILUTE PARCEL FROM THE INITIAL CONVECTIVE
C LEVEL, (MIMICKING CODE IN ROUTINE TERM_CON) AND RESET L_SHALLOW
C TO FALSE IF THIS PARCEL IS STILL BUOYANT.
C----------------------------------------------------------------------
C
         IF(K_TEST.EQ.NBL.OR.K_TEST.EQ.NBL-1)THEN
          IF(BWORK2(I,1))THEN
            EL=LC
          ELSE
            EL=LC+LF
          END IF
          THVUNDI=(THPI(INDEX1(I))+(EL/(WORK2(I,10)*CP))*(QPI(INDEX1(I))
     *            -WORK2(I,5))+((LC-EL)/(WORK2(I,10)*CP))*MAX(0.0,
     *            (QPI(INDEX1(I))-QSTICE)))*(1.0+C_VIRTUAL*WORK2(I,5))
          THVEKP1=(WORK2(I,2)*(1.0+C_VIRTUAL*WORK2(I,4))+XSBMIN)
          IF(THVUNDI.GT.THVEKP1)THEN
            L_SHALLOW_C(I)=.FALSE.
            L_SHALLOW(INDEX1(I))=L_SHALLOW_C(I)
          END IF
         END IF            ! IF(K_TEST.EQ.NBL.OR.K_TEST.EQ.NBL-1) END
         BWORK2(I,3) = .FALSE.
        END IF             ! IF(BWORK2(I,3.AND..NOT.BWORK(I,4)...) END
C
       END DO              ! END OF 1ST I LOOP OVER CONVECTIVE POINTS
C
      END DO               ! END OF LOOP OVER BOUNDARY LAYER LEVELS
C
C----------------------------------------------------------------------
C RESET INITIAL THETA AND Q OF THE PARCEL AND ENTRAINMENT/
C DETRAINMENT RATES IF SHALLOW CONVECTION
C---------------------------------------------------------------------
C
! Fujitsu vectorization directive
!OCL NOVREC
       DO I=1,NCONV        ! 2ND LOOP OVER CONVECTING POINTS
C
        IF ( L_SHALLOW_C(I) ) THEN
C
        IF (.NOT.BWORK(I,4)) THEN
C
         IF ( L_SDXS .AND. K .EQ. 1 ) THEN  
           WORK(I,7) = WORK(I,1) + MAX(THPIXS_SHALLOW,
     *                 T1_SD(INDEX1(I))/EXK(INDEX1(I)))
           THPI(INDEX1(I)) = WORK(I,1) + MAX(THPIXS_SHALLOW,
     *                       T1_SD(INDEX1(I))/EXK(INDEX1(I)))
           WORK(I,8) = WORK(I,3) + MAX(QPIXS_SHALLOW,Q1_SD(INDEX1(I)))
           QPI(INDEX1(I)) = WORK(I,3) + MAX(QPIXS_SHALLOW,
     *                      Q1_SD(INDEX1(I)))
         ELSE
           WORK(I,7) = WORK(I,1) + THPIXS_SHALLOW
           THPI(INDEX1(I)) = WORK(I,1) + THPIXS_SHALLOW
           WORK(I,8) = WORK(I,3) + QPIXS_SHALLOW
           QPI(INDEX1(I)) = WORK(I,3) + QPIXS_SHALLOW
         END IF
C
        END IF             ! IF(.NOT.BWORK(I,4)) END
C
         WORK(I,11) = WORK(I,11)*SH_FAC
         EKP14(INDEX1(I)) = WORK(I,11)
         WORK(I,12) = WORK(I,12)*SH_FAC
         EKP34(INDEX1(I)) = WORK(I,12)
         AMDETK(INDEX1(I)) = AMDETK(INDEX1(I))*SH_FAC
C
        END IF             ! IF(L_SHALLOW_C(I)) END
C
       END DO              ! END OF 2ND I LOOP OVER CONVECTING POINTS
C
      END IF               ! IF(K.LT.NBL) END
CL
CL---------------------------------------------------------------------
CL LIFT PARCEL FROM LAYER K TO K+1
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (3) AND (4)
CL---------------------------------------------------------------------
CL
      CALL LIFT_PAR (NCONV,NPNTS,WORK(1,13),WORK(1,14),WORK(1,15),
     *               BWORK(1,2),BWORK(1,1),WORK(1,7),WORK(1,8),
     *               WORK(1,2),WORK(1,4),WORK(1,1),WORK(1,3),
     *               WORK(1,5),WORK(1,6),WORK(1,9),
     *               WORK(1,10),WORK(1,11),WORK(1,12),L_MOM,
     *               WORK(1,29),WORK(1,30),WORK(1,27),WORK(1,28),
     *               WORK(1,23),WORK(1,24),WORK(1,25),WORK(1,26),
     *               L_TRACER,NTRA,TRAPKP1_C,TRAPK_C,TRAEKP1_C,
     *               TRAEK_C,L_SHALLOW_C)
C
      DO 110 I=1,NCONV
CL
CL---------------------------------------------------------------------
CL CALCULATE BUOYANCY OF PARCEL IN LAYER K+1
CL---------------------------------------------------------------------
CL
        WORK(I,16) = WORK(I,13)*(1.0 +
     *                            C_VIRTUAL * WORK(I,14))
     *               - WORK(I,2)*(1.0 +
     *                            C_VIRTUAL * WORK(I,4))
C
C----------------------------------------------------------------------
C INITIATE CONVECTION WHERE BUOYANCY IS LARGE ENOUGH
C----------------------------------------------------------------------
C
        BWORK(I,3) = .NOT.BWORK(I,4) .AND. WORK(I,16) .GT.
     *      (WORK(I,20)+ XSBMIN)
C
C----------------------------------------------------------------------
C CALCULATE INITIAL MASSFLUX FROM LAYER K
C----------------------------------------------------------------------
C
        IF ( BWORK(I,3) ) THEN
C
          IF(L_SHALLOW_C(I))THEN
            C=C_SHALLOW
            D=D_SHALLOW
          ELSEIF(L_MID_C(I))THEN
            C=C_MID
            D=D_MID
          ELSE
            C=C_DEEP
            D=D_DEEP
          END IF
C
           WORK(I,19) = 1.0E-3 * WORK(I,18) *
     1                        ( D + C * WORK(I,18) *
     2                 ((WORK(I,16) - XSBMIN) / WORK(I,17)))
C
        END IF
  110 CONTINUE
C
C----------------------------------------------------------------------
C LIMIT MASSFLUX IN LOWEST CONVECTING LAYER TO BE <= MASS OF LAYER
C OR
C IF K=1 ADJUST ENTRAINMENT RATE IN BOTTOM HALF OF LAYER 2 SO
C NOT TO AFFECT THE MASS FLUX AT MID-POINT OF LAYER 2
C----------------------------------------------------------------------
C
      IF ( K .EQ. 1 ) THEN
C
       DO I=1,NCONV
C
C--------------------------------------------------------------------
C CARRY OUT CALCULATION IF CONVECTION WAS INITIATED FROM LAYER 1
C--------------------------------------------------------------------
C
        IF ( BWORK(I,3) ) THEN
C
C--------------------------------------------------------------------
C CALCULATE MASS FLUX AT MID-POINT OF LAYER 2 USING STANDARD
C ENTRAINMENT RATES
C--------------------------------------------------------------------
C
         FLX2 = WORK(I,19) * (1.0 + WORK(I,11)) * (1.0 + WORK(I,12))
C
C--------------------------------------------------------------------
C IF MASS FLUX IN LAYER 2 EXCEEDS MASS OF LAYER THEN LIMIT MASS FLUX
C OVER A TIMESTEP TO MASS OF LAYER
C--------------------------------------------------------------------
C
        IF (WORK(I,19) .GT. WORK(I,21)) THEN
C
         WORK(I,19) = WORK(I,21)
C
C--------------------------------------------------------------------
C IF MASS FLUX AT MID-POINT OF LAYER 2 EXCEEDS THE MASS OF THE COLUMN
C DOWN TO THE SURFACE OVER THE TIMESTEP THEN LIMIT MASS FLUX
C--------------------------------------------------------------------
C
        IF ( FLX2 .GT. WORK(I,22)) FLX2 = WORK(I,22)
C
C--------------------------------------------------------------------
C ADJUST ENTRAINMENT RATE IN BOTTOM HALF OF LAYER 2
C--------------------------------------------------------------------
C
       WORK(I,12) = (FLX2/(WORK(I,19) * (1.0 + WORK(I,11)))) - 1.0
       END IF
C
       END IF
      END DO
C
C---------------------------------------------------------------------
C RECALCULATE ASCENT FROM LAYER 1 TO 2 USING ADJUSTED ENTRAINMENT RATE
C---------------------------------------------------------------------
C
      CALL LIFT_PAR (NCONV,NPNTS,WORK(1,13),WORK(1,14),WORK(1,15),
     *               BWORK(1,2),BWORK(1,1),WORK(1,7),WORK(1,8),
     *               WORK(1,2),WORK(1,4),WORK(1,1),WORK(1,3),
     *               WORK(1,5),WORK(1,6),WORK(1,9),
     *               WORK(1,10),WORK(1,11),WORK(1,12),L_MOM,
     *               WORK(1,29),WORK(1,30),WORK(1,27),WORK(1,28),
     *               WORK(1,23),WORK(1,24),WORK(1,25),WORK(1,26),
     *               L_TRACER,NTRA,TRAPKP1_C,TRAPK_C,TRAEKP1_C,
     *               TRAEK_C,L_SHALLOW_C)
C
       DO I=1,NCONV
C
        IF ( BWORK(I,3) ) THEN
CL
CL---------------------------------------------------------------------
CL RECALCULATE BUOYANCY OF PARCEL IN LAYER K+1
CL---------------------------------------------------------------------
CL
        WORK(I,16) = WORK(I,13)*(1.0 +
     *                            C_VIRTUAL * WORK(I,14))
     *               - WORK(I,2)*(1.0 +
     *                            C_VIRTUAL * WORK(I,4))
C
C----------------------------------------------------------------------
C RESET MASK TO INITIATE CONVECTION WHERE BUOYANCY IS LARGE ENOUGH
C----------------------------------------------------------------------
C
        BWORK(I,3) = .NOT.BWORK(I,4) .AND. WORK(I,16) .GT.
     *      (WORK(I,20)+ XSBMIN)
C
        BWORK(I,4) = BWORK(I,4) .OR. BWORK(I,3)
C
       END IF
C
       FLX(INDEX1(I),K) = WORK(I,19)
       IF(FLG_UP_FLX) UP_FLUX(INDEX1(I),K)=WORK(I,19)
C
      END DO
C
C----------------------------------------------------------------------
C END OF CALCULATION FOR LAYER 1
C----------------------------------------------------------------------
C
      ELSE
C
       DO I=1,NCONV
C
C----------------------------------------------------------------------
C IF MASS FLUX OUT OF THE INITIAL LAYER IS GREATER THAN THE MASS OF
C THE LAYER OVER THE TIMESTEP THEN LIMIT MASS FLUX TO MASSS OF LAYER
C----------------------------------------------------------------------
C
        IF (BWORK(I,3) .AND. WORK(I,19).GT.WORK(I,21))
     1                 WORK(I,19) = WORK(I,21)
C
        BWORK(I,4) = BWORK(I,4) .OR. BWORK(I,3)
C
        FLX(INDEX1(I),K) = WORK(I,19)
        IF(FLG_UP_FLX) UP_FLUX(INDEX1(I),K)=WORK(I,19)
C
       END DO
C
      END IF
C
CL
CL--------------------------------------------------------------------
CL ZERO MIXING DETRAINMENT RATE WHEN CONVECTION STARTS FROM LAYER K
CL STORE DIAGNOSTIC LINKED TO INITIAL CONVECTIVE MASSFLUX FOR
CL CALCULATION OF FINAL CLOSURE FOR DEEP CONVECTION.
CL--------------------------------------------------------------------
CL
      DO I=1,NCONV
       IF ( BWORK(I,3) )THEN
        AMDETK(INDEX1(I))=0.0
        FLX_INIT(INDEX1(I))=WORK(I,19)
        START_LEV(INDEX1(I))=K
        FLXMAX_INIT(INDEX1(I))=WORK(I,21)
       END IF
      END DO
CL
CL--------------------------------------------------------------------
CL COMPRESS DOWN THOSE POINTS WHICH ARE NOT BUOYANT IN LAYER K+1.
CL--------------------------------------------------------------------
CL
      NINIT = 0
      DO 115 I=1,NCONV
        IF(BWORK(I,4))THEN
          NINIT = NINIT + 1
          INDEX2(NINIT) = I
        END IF
  115 CONTINUE
C
C
C----------------------------------------------------------------------
C  WORK SPACE USAGE FOR SECOND COMPRESS ON BASIS OF WHETHER
C  PARCEL A PARCEL STARTING FROM LAYER K IS BUOYANT IN LAYER
C  K+1 OR IF CONVECTION ALREADY EXISTS IN LAYER K
C
C  REFERENCES TO WORK, WORK2, BWORK AND BWORK2
C  REFER TO STARTING ADDRESS
C
C  LENGTH OF COMPRESSES DATA = NINIT
C
C  WORK2 AND BWORK2 ARE COMPRESSED DOWN FROM COMPRESSED
C  ARRAYS STORED IN WORK AND BWORK AFTER FIST COMPRESS
C
C  WORK2(1,1)  = TH(#,K)
C  WORK2(1,2)  = TH(#,K+1)
C  WORK2(1,3)  = Q(#,K)
C  WORK2(1,4)  = Q(#,K+1)
C  WORK2(1,5)  = QSE(#,K+1)
C  WORK2(1,6)  = DQSTHKP1(#)
C  WORK2(1,7)  = THP(#,K)
C  WORK2(1,8)  = QP(#,K)
C  WORK2(1,9)  = PKP1(#)
C  WORK2(1,10) = EXKP1(#)
C  WORK2(1,11) = EKP14(#)
C  WORK2(1,12) = EKP34(#)
C  WORK2(1,13) = PARCEL POT. TEMPERATURE IN LAYER K+1
C  WORK2(1,14) = PARCEL MIXING RATIO IN LAYER K+1
C  WORK2(1,15) = EXCESS WATER VAPOUR IN PARCEL ABOVE
C               SATURATION AFTER DRY ASCENT
C  WORK2(1,16) = PARCEL BUOYANCY IN LAYER K+1
C  WORK2(1,17) = NOT USED IN THIS SECTION
C  WORK2(1,18) = PSTAR(#)
C  WORK2(1,19) = FLX(#,K)
C
C  BWORK2(1,1) = BWATER(INDEX1(I),K+1)
C  BWORK2(1,2) = .TRUE. IF PARCEL SATURATED IN LAYER K+1
C  BWORK2(1,3) = .TRUE. IF CONVECTION INITIATE FROM LAYER K+1
C  WORK2(1,23) = U(#,K)
C  WORK2(1,24) = U(#,K+1)
C  WORK2(1,25) = V(#,K)
C  WORK2(1,26) = V(#,K+1)
C  WORK2(1,27) = UP(#,K)
C  WORK2(1,28) = VP(#,K)
C  WORK2(1,29) = PARCEL U IN LAYER K+1
C  WORK2(1,30) = PARCEL V IN LAYER K+1
C
C  WORK AND BWORK NOW CONTAIN DATA COMPRESSED DOWN
C  FROM FULL LENGTH VECTORS
C
C  WORK(1,1) = not used in this section
C  WORK(1,2) = QSE(#,K)
C  WORK(1,3) = DQSTHK(#)
C  WORK(1,4) = THPI(#)
C  WORK(1,5) = QPI(#)
C  WORK(1,6) = XPK(#,K+1)
C  WORK(1,7) = not used in this section
C  WORK(1,8) = DEPTH(#)
C  WORK(1,9) = PRECIP(#,K+1)
C  WORK(1,10) = DTHBYDT(#,K)
C  WORK(1,11) = DQBYDT(#,K)
C  WORK(1,12) = DTHBYDT(#,K+1)
C  WORK(1,13) = DQBYDT(#,K+1)
C  WORK(1,14) = AMDETK(#)
C  WORK(1,15) = NOY USED IN THIS SECTION
C  WORK(1,16) = PK(#)
C  WORK(1,17) = EXK(#)
C  WORK(1,18) = DELEXKP1(#)
C  WORK(1,19) = DELPK(#)
C  WORK(1,20) = DELPKP1(#)
C  WORK(1,21) = CCW(#,K+1)
C  WORK(1,22) = T1_SD(#)
C  WORK(1,23) = Q1_SD(#)
C  WORK(1,24) = DUBYDT(#,K)
C  WORK(1,25) = DUBYDT(#,K+1)
C  WORK(1,26) = DVBYDT(#,K)
C  WORK(1,27) = DVBYDT(#,K+1)
C  WORK(1,28) = EFLUX_U_UD(#)
C  WORK(1,29) = EFLUX_V_UD(#)
C
C  BWORK(1,1) = BGMK(#)
C  BWORK(1,2) = BLAND(#)
C  BWORK(1,3) = BTERM(#)
C  BWORK(1,2) = BLAND(#)
C----------------------------------------------------------------------
C
      IF (NINIT .NE. 0) THEN
C
C-----------------------------------------------------------------------
C FIRST COMPRESS DOWN QUANTITIES FROM PREVIOUSLY COMPRESSED ARRAY
C-----------------------------------------------------------------------
C
        DO 120 I=1,NINIT
          WORK2(I,1)  = WORK(INDEX2(I),1)
          WORK2(I,2)  = WORK(INDEX2(I),2)
          WORK2(I,3)  = WORK(INDEX2(I),3)
          WORK2(I,4)  = WORK(INDEX2(I),4)
          WORK2(I,5)  = WORK(INDEX2(I),5)
          WORK2(I,6)  = WORK(INDEX2(I),6)
          WORK2(I,7)  = WORK(INDEX2(I),7)
          WORK2(I,8)  = WORK(INDEX2(I),8)
          WORK2(I,9)  = WORK(INDEX2(I),9)
          WORK2(I,10) = WORK(INDEX2(I),10)
          WORK2(I,11) = WORK(INDEX2(I),11)
          WORK2(I,12) = WORK(INDEX2(I),12)
          WORK2(I,13) = WORK(INDEX2(I),13)
          WORK2(I,14) = WORK(INDEX2(I),14)
          WORK2(I,15) = WORK(INDEX2(I),15)
          WORK2(I,16) = WORK(INDEX2(I),16)
          WORK2(I,17) = WORK(INDEX2(I),17)
          WORK2(I,18) = WORK(INDEX2(I),18)
          WORK2(I,19) = WORK(INDEX2(I),19)
          BWORK2(I,1) = BWORK(INDEX2(I),1)
          BWORK2(I,2) = BWORK(INDEX2(I),2)
          BWORK2(I,3) = BWORK(INDEX2(I),3)
          L_SHALLOW_C2(I) = L_SHALLOW_C(INDEX2(I))
          L_MID_C2(I) = L_MID_C(INDEX2(I))
  120   CONTINUE
C
        IF(L_MOM)THEN
         DO I=1,NINIT
          WORK2(I,23) = WORK(INDEX2(I),23)
          WORK2(I,24) = WORK(INDEX2(I),24)
          WORK2(I,25) = WORK(INDEX2(I),25)
          WORK2(I,26) = WORK(INDEX2(I),26)
          WORK2(I,27) = WORK(INDEX2(I),27)
          WORK2(I,28) = WORK(INDEX2(I),28)
          WORK2(I,29) = WORK(INDEX2(I),29)
          WORK2(I,30) = WORK(INDEX2(I),30)
         END DO
        END IF
C
        IF(L_TRACER)THEN
C
        DO KTRA=1,NTRA
          DO I=1,NINIT
            TRAEK_C2(I,KTRA)=TRAEK_C(INDEX2(I),KTRA)
            TRAEKP1_C2(I,KTRA)=TRAEKP1_C(INDEX2(I),KTRA)
            TRAPK_C2(I,KTRA)=TRAPK_C(INDEX2(I),KTRA)
            TRAPKP1_C2(I,KTRA)=TRAPKP1_C(INDEX2(I),KTRA)
          END DO
        END DO
C
        END IF
C----------------------------------------------------------------------
C COMPRESS DOWN REST OF DATA FROM FULL ARRAYS
C
C FIRST EXPAND BACK BWORK(1,2) (=BINIT) BACK TO FULL VECTORS
C----------------------------------------------------------------------
C
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO 130 I=1,NCONV
        BINIT(INDEX1(I)) = BWORK(I,4)
  130 CONTINUE
C
      NINIT = 0
      DO 135 I=1,NPNTS
        IF(BINIT(I))THEN
          NINIT = NINIT + 1
          INDEX3(NINIT) = I
        END IF
  135 CONTINUE
C
      DO 140 I=1,NINIT
        WORK(I,2) = QSE(INDEX3(I),K)
        WORK(I,3) = DQSTHK(INDEX3(I))
        WORK(I,4) = THPI(INDEX3(I))
        WORK(I,5) = QPI(INDEX3(I))
        WORK(I,6) = XPK(INDEX3(I),K)
        WORK(I,8) = DEPTH(INDEX3(I))
        CCA_2DC(I)    = CCA_2D(INDEX3(I))
        ICCBC(I)   = ICCB(INDEX3(I))
        ICCTC(I)   = ICCT(INDEX3(I))
        TCWC(I)    = TCW(INDEX3(I))
        CCLWPC(I)  = CCLWP(INDEX3(I))
        LCCAC(I)   = LCCA(INDEX3(I))   ! beware - LCCAC & LCBASEC
        LCBASEC(I) = LCBASE(INDEX3(I)) ! are IN/OUT to lower levels
        LCTOPC(I)  = LCTOP(INDEX3(I))
        LCCLWPC(I) = LCCLWP(INDEX3(I))
        BWORK(I,1) = BGMK(INDEX3(I))
        BWORK(I,2) = BLAND(INDEX3(I))
        WORK(I,10) = DTHBYDT(INDEX3(I),K)
        WORK(I,11) = DQBYDT(INDEX3(I),K)
        WORK(I,12) = DTHBYDT(INDEX3(I),K+1)
        WORK(I,13) = DQBYDT(INDEX3(I),K+1)
        WORK(I,14) = AMDETK(INDEX3(I))
        WORK(I,16) = PK(INDEX3(I))
        WORK(I,17) = EXK(INDEX3(I))
        WORK(I,18) = DELEXKP1(INDEX3(I))
        WORK(I,19) = DELPK(INDEX3(I))
        WORK(I,20) = DELPKP1(INDEX3(I))
        WORK(I,22) = T1_SD(INDEX3(I))
        WORK(I,23) = Q1_SD(INDEX3(I))
        CAPE_C(I)  = CAPE(INDEX3(I))
        DCPBYDT_C(I) = DCPBYDT(INDEX3(I))
C
        BWORK(I,4) = .TRUE.
  140 CONTINUE
C
      IF(L_MOM)THEN
       DO I=1,NINIT
        WORK(I,24) = DUBYDT(INDEX3(I),K)
        WORK(I,25) = DUBYDT(INDEX3(I),K+1)
        WORK(I,26) = DVBYDT(INDEX3(I),K)
        WORK(I,27) = DVBYDT(INDEX3(I),K+1)
        WORK(I,28) = EFLUX_U_UD(INDEX3(I))
        WORK(I,29) = EFLUX_V_UD(INDEX3(I))
       END DO
      END IF
C
        IF(L_TRACER)THEN
C
        DO KTRA=1,NTRA
          DO I=1,NINIT
            DTRAEK_C(I,KTRA) = DTRABYDT(INDEX3(I),K,KTRA)
            DTRAEKP1_C(I,KTRA) = DTRABYDT(INDEX3(I),K+1,KTRA)
          END DO
        END DO
C
        END IF
C
CL
CL----------------------------------------------------------------------
CL CALCULATE REST OF PARCEL ASCENT AND EFFECT OF CONVECTION
CL UPON THE LARGE-SCALE ATMOSPHERE
CL
CL SUBROUTINE CONVEC2
CL
CL UM DOCUMENTATION PAPER P27
CL SECTIONS (5),(6),(7),(8),(9),(10)
CL----------------------------------------------------------------------
CL
      CALL CONVEC2 (NINIT,NPNTS,NLEV,K,WORK2(1,1),WORK2(1,2),WORK2(1,3),
     *             WORK2(1,4),WORK2(1,5),WORK2(1,6),WORK2(1,18),
     *             WORK2(1,7),WORK2(1,8),WORK2(1,13),WORK2(1,14),
     *             WORK2(1,15),WORK2(1,16),WORK(1,2),WORK(1,3),
     *             WORK(1,4),WORK(1,5),WORK(1,6),WORK2(1,19),
     *             BWORK2(1,1),BWORK2(1,2),BWORK(1,1),BWORK2(1,3),
     *             BWORK(1,2),BWORK(1,3),WORK(1,8),WORK(1,9),
     *             WORK(1,10),WORK(1,11),WORK(1,12),WORK(1,13),
     *             BWORK(1,4),CCA_2DC,ICCBC,ICCTC,TCWC,
     *             WORK2(1,11),WORK2(1,12),WORK(1,14),
     *             WORK(1,16),WORK2(1,9),WORK(1,17),WORK2(1,10),
     *             WORK(1,18),WORK(1,19),WORK(1,20),
     *             CCLWPC,WORK(1,21),LCCAC,LCBASEC,LCTOPC,LCCLWPC,
     *             WORK(1,22),WORK(1,23),L_MOM,WORK2(1,23),WORK2(1,24),
     *             WORK2(1,25),WORK2(1,26),WORK2(1,27),WORK2(1,28),
     *             WORK2(1,29),WORK2(1,30),WORK(1,24),WORK(1,25),
     *             WORK(1,26),WORK(1,27),WORK(1,28),WORK(1,29),
     *             L_SHALLOW_C2,L_MID_C2,
     *             L_TRACER,NTRA,TRAEK_C2,TRAEKP1_C2,TRAPK_C2,
     *             TRAPKP1_C2,DTRAEK_C,DTRAEKP1_C,CAPE_C,DCPBYDT_C,     
     &             L_XSCOMP,L_SDXS,L_CCW,MPARWTR,UD_FACTOR,
     &             DELTAK)
CL
CL---------------------------------------------------------------------
CL EXPAND REQUIRED VECTORS BACK TO FULL FIELDS
CL----------------------------------------------------------------------
CL
      DO 145 I=1,NPNTS
        THP(I,K+1) = 0.0
        QP(I,K+1) = 0.0
        XPK(I,K+1) = 0.0
        FLX(I,K+1)= 0.0
        DEPTH(I) = 0.0
        PRECIP(I,K+1) = 0.0
        BGMK(I) = .FALSE.
        BTERM(I) = .FALSE.
        BINIT(I) = .FALSE.
  145 CONTINUE
C
      IF(L_MOM)THEN
       DO I=1,NPNTS
        UP(I,K+1) = 0.0
        VP(I,K+1) = 0.0
       END DO
      END IF
C
      IF(L_TRACER)THEN
C
      DO KTRA=1,NTRA
        DO I=1,NPNTS
          TRAP(I,K+1,KTRA) = 0.0
        END DO
      END DO
C
      END IF
C
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO 150 I=1,NINIT
        THP(INDEX3(I),K+1) = WORK2(I,7)
        QP(INDEX3(I),K+1) = WORK2(I,8)
        XPK(INDEX3(I),K+1) = WORK(I,6)
        FLX(INDEX3(I),K+1) = WORK2(I,19)
        DEPTH(INDEX3(I)) = WORK(I,8)
        PRECIP(INDEX3(I),K+1) = WORK(I,9)
        DTHBYDT(INDEX3(I),K) = WORK(I,10)
        DQBYDT(INDEX3(I),K) = WORK(I,11)
        DTHBYDT(INDEX3(I),K+1) = WORK(I,12)
        DQBYDT(INDEX3(I),K+1) = WORK(I,13)
        CCA_2D(INDEX3(I)) = CCA_2DC(I)
        ICCB(INDEX3(I)) = ICCBC(I)
        ICCT(INDEX3(I)) = ICCTC(I)
        TCW(INDEX3(I)) = TCWC(I)
        CCLWP(INDEX3(I)) = CCLWPC(I)
        LCCA(INDEX3(I)) = LCCAC(I)
        LCBASE(INDEX3(I)) = LCBASEC(I)
        LCTOP(INDEX3(I)) = LCTOPC(I)
        LCCLWP(INDEX3(I)) = LCCLWPC(I)
        CCW(INDEX3(I),K+1) = WORK(I,21)
        CAPE(INDEX3(I)) = CAPE_C(I)
        DCPBYDT(INDEX3(I)) = DCPBYDT_C(I)
C
        BGMK(INDEX3(I)) = BWORK(I,1)
        BTERM(INDEX3(I)) = BWORK(I,3)
        BINIT(INDEX3(I)) = BWORK(I,4)
        IF(FLG_UP_FLX) UP_FLUX(INDEX3(I),K+1)=WORK2(I,19)
        IF(FLG_ENTR_UP) ENTRAIN_UP(INDEX3(I),K)=(1.0-DELTAK(I))*
     &               (1.0-WORK(I,14))*(WORK2(I,11)+WORK2(I,12)*
     &               (1.0+WORK2(I,11)))*FLX(INDEX3(I),K)
        IF(FLG_DETR_UP) DETRAIN_UP(INDEX3(I),K)=-(WORK(I,14)+
     &                          DELTAK(I)*(1.0-WORK(I,14)))*
     &                          FLX(INDEX3(I),K)
        IF(BTERM(INDEX3(I))) THEN
!
! TERMINAL DETRAINMENT
!
         IF(FLG_ENTR_UP) ENTRAIN_UP(INDEX3(I),K+1)=0.0
         IF(FLG_DETR_UP) DETRAIN_UP(INDEX3(I),K+1)=-(1.0-DELTAK(I))*
     &                                              FLX(INDEX3(I),K) 
        ENDIF     
  150 CONTINUE
C
      IF(L_MOM)THEN
       DO I=1,NINIT
        UP(INDEX3(I),K+1) = WORK2(I,27)
        VP(INDEX3(I),K+1) = WORK2(I,28)
        DUBYDT(INDEX3(I),K) = WORK(I,24)
        DVBYDT(INDEX3(I),K) = WORK(I,26)
        DUBYDT(INDEX3(I),K+1) = WORK(I,25)
        DVBYDT(INDEX3(I),K+1) = WORK(I,27)
        EFLUX_U_UD(INDEX3(I)) = WORK(I,28)
        EFLUX_V_UD(INDEX3(I)) = WORK(I,29)
       END DO
      END IF
C
      IF(L_TRACER)THEN
C
      DO KTRA=1,NTRA
        DO I=1,NINIT
          TRAP(INDEX3(I),K+1,KTRA)=TRAPK_C2(I,KTRA)
          DTRABYDT(INDEX3(I),K,KTRA)=DTRAEK_C(I,KTRA)
          DTRABYDT(INDEX3(I),K+1,KTRA)=DTRAEKP1_C(I,KTRA)
        END DO
      END DO
C
      END IF
C
C
      END IF
C
      END IF
C
C-------------------------------------------------------------------
C ADJUSTMENT OF CLOSURE FOR DEEP CONVECTION
C
C UM DOCUMENTATION PAPER P27-3. SECTION 5.
C
C ADJUST INITIAL MASS FLUX SO THAT CAPE IS REMOVED BY CONVECTION
C OVER TIMESCALE CAPE_TS
C-------------------------------------------------------------------
C
C
      DO I=1,NPNTS
      IF(L_CAPE)THEN
       IF(.NOT.L_SHALLOW(I).AND.BTERM(I))THEN
        IF(DCPBYDT(I).GT.0.0)THEN
          FLX_INIT_NEW(I)=FLX_INIT(I)*CAPE(I)/(CAPE_TS*DCPBYDT(I))
          IF(FLX_INIT_NEW(I).GT.FLXMAX_INIT(I))THEN
            FLX_INIT_NEW(I)=FLXMAX_INIT(I)
          END IF
        END IF
       END IF
      END IF
       IF(BTERM(I))THEN
        CAPE_OUT(I)=CAPE(I)
        CAPE(I)=0.0
        DCPBYDT(I)=0.0
       END IF
      END DO
C
C---------------------------------------------------------------------
C RESCALE Q1, Q2 MASS FLUX AND PRECIP FOR DEEP CONVECTION
C---------------------------------------------------------------------
C
      IF(L_CAPE)THEN
      DO KT=1,K+1
      DO I=1,NPNTS
       IF(KT.GE.START_LEV(I).AND..NOT.L_SHALLOW(I).AND.BTERM(I).
     *  AND.FLX_INIT_NEW(I).GT.0.0)THEN
          IF(KT.EQ.DET_LEV(I))THEN
            DTHBYDT(I,KT)=(DTHBYDT(I,KT) - DTHEF(I))
     *                     *FLX_INIT_NEW(I)/FLX_INIT(I)
            DTHBYDT(I,KT) = DTHBYDT(I,KT) + DTHEF(I)
            DQBYDT(I,KT)=(DQBYDT(I,KT) - DQF(I))
     *                     *FLX_INIT_NEW(I)/FLX_INIT(I)
            DQBYDT(I,KT) = DQBYDT(I,KT) +DQF(I)
            IF(L_MOM) THEN
             DUBYDT(I,KT)=(DUBYDT(I,KT)-DUEF(I))*
     &                    FLX_INIT_NEW(I)/FLX_INIT(I)
             DUBYDT(I,KT)=DUBYDT(I,KT)+DUEF(I)
             DVBYDT(I,KT)=(DVBYDT(I,KT)-DVEF(I))*
     &                    FLX_INIT_NEW(I)/FLX_INIT(I)
             DVBYDT(I,KT)=DVBYDT(I,KT)+DVEF(I)
            ENDIF
          ELSE
            DTHBYDT(I,KT)=DTHBYDT(I,KT)*FLX_INIT_NEW(I)/FLX_INIT(I)
            DQBYDT(I,KT)=DQBYDT(I,KT)*FLX_INIT_NEW(I)/FLX_INIT(I)
            IF(L_MOM) THEN
             DUBYDT(I,KT)=DUBYDT(I,KT)*FLX_INIT_NEW(I)/FLX_INIT(I)
             DVBYDT(I,KT)=DVBYDT(I,KT)*FLX_INIT_NEW(I)/FLX_INIT(I)
            ENDIF
          END IF
            FLX(I,KT)=FLX(I,KT)*FLX_INIT_NEW(I)/FLX_INIT(I)
            IF(FLG_UP_FLX) UP_FLUX(I,KT)=FLX(I,KT)
            IF(FLG_ENTR_UP) ENTRAIN_UP(I,KT)=ENTRAIN_UP(I,KT)*
     &                                     FLX_INIT_NEW(I)/FLX_INIT(I)
            IF(FLG_DETR_UP) DETRAIN_UP(I,KT)=DETRAIN_UP(I,KT)*
     &                                     FLX_INIT_NEW(I)/FLX_INIT(I)
            PRECIP(I,KT)=PRECIP(I,KT)*FLX_INIT_NEW(I)/FLX_INIT(I)
       END IF
      END DO
      END DO
      DO I=1,NPNTS  
       IF(.NOT.L_SHALLOW(I).AND.BTERM(I).AND.
     &                                   FLX_INIT_NEW(I).GT.0.0)THEN  
        IF(CCA_2D(I).GT.2.0E-5) CCA_2D(I)=CCA_2D(I)+
     &                            0.06*LOG(FLX_INIT_NEW(I)/FLX_INIT(I)) 
       ENDIF
      END DO  
      END IF
CL
CL---------------------------------------------------------------------
CL DOWNDRAUGHT CALCULATION
CL
CL CARRIED OUT FOR THOSE CLOUD WHICH ARE TERMINATING
CL
CL SUBROUTINE DD_CALL
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (11)
CL---------------------------------------------------------------------
CL
C
      NTERM = 0
      DO 160 I=1,NPNTS
        IF (BTERM(I)) THEN
         NTERM = NTERM + 1
         DTHEF(I) = DTHBYDT(I,K+1)
         DQF(I)   = DQBYDT(I,K+1)
         IF(L_MOM) THEN
          DUEF(I)=DUBYDT(I,K+1)
          DVEF(I)=DVBYDT(I,K+1)
         ENDIF
 
         DET_LEV(I) = K+1
        END IF
  160 CONTINUE
C
      IF (NTERM .NE. 0) THEN
C
         CALL DD_CALL (NP_FIELD,NPNTS,K,THP(1,1),QP(1,1),TH(1,1),
     *                 Q(1,1),DTHBYDT(1,1),DQBYDT(1,1),FLX(1,1),
     *                 PSTAR,AK,BK,AKM12,BKM12,DELAK,DELBK,EXNER(1,1),
     *                 PRECIP(1,1),RAIN,SNOW,ICCB,ICCT,BWATER(1,2),
     *                 BTERM,BGMK,TIMESTEP,CCA_2D,NTERM,L_MOM,UP(1,1),
     *                 VP(1,1),U(1,1),V(1,1),DUBYDT(1,1),DVBYDT(1,1),
     *                 EFLUX_U_DD,EFLUX_V_DD,
     *                 L_TRACER,NTRA,TRAP,TRACER,DTRABYDT,NLEV,TRLEV,
     &                 recip_pstar,                                     
     &                 DWN_FLUX,FLG_DWN_FLX,ENTRAIN_DWN,
     &                 FLG_ENTR_DWN,DETRAIN_DWN,FLG_DETR_DWN)
 
C
C---------------------------------------------------------------------
C ZERO CONVECTION START LEVEL IF CONVECTION TERMINATES
C---------------------------------------------------------------------
C
      DO I=1,NPNTS
       IF(BTERM(I))THEN
       START_LEV(I)=0.0
       END IF
      END DO
C
C---------------------------------------------------------------------
C ADJUSTMENT TO CLOUD BASE, TOP AND AMOUNT
C
C IF CLOUD BASE AND TOP ARE EQUAL THEN ERRORS OCCUR IN RADIATION SCHEME
C
C ONLY OCCURS IF CONVECTION SATURATES UPON FORCED DETRAINMENT
C
C WHEN OCCURS ZERO CLOUD BASE, TOP AND AMOUNT
C
C---------------------------------------------------------------------
C
      DO I=1,NPNTS
        IF (BTERM(I) .AND. ICCB(I) .EQ. ICCT(I)) THEN
          ICCB(I) = 0.0
          ICCT(I) = 0.0
          CCA_2D(I) = 0.0
          TCW(I) = 0.0
          CCLWP(I) = 0.0
        END IF
        IF (BTERM(I) .AND. LCBASE(I) .EQ. LCTOP(I)) THEN
          LCBASE(I) = 0
          LCTOP(I) = 0
          LCCA(I) = 0.0
          LCCLWP(I) = 0.0
        END IF
      END DO
C
C---------------------------------------------------------------------
C RESET BTERM TO FALSE
C---------------------------------------------------------------------
C
      DO 200 I=1,NPNTS
  200  BTERM(I) = .FALSE.
C
      END IF
CL
CL=====================================================================
CL END OF MAIN LOOP
CL=====================================================================
CL
  60  CONTINUE
CL
CL---------------------------------------------------------------------
CL BALANCE ENERGY BUDGET BY APPLYING CORRECTION TO THE TEMPERATURES
CL
CL SUBROUTINE COR_ENGY
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (12)
CL---------------------------------------------------------------------
CL
      NCNLV = 0
      DO 210 I=1,NPNTS
        IF(BCNLV(I))THEN
          NCNLV = NCNLV + 1
          INDEX4(NCNLV) = I
        END IF
  210 CONTINUE
C
C
C----------------------------------------------------------------------
C WORK SPACE USAGE FOR ENERGY CORRECTION CALCULATION
C
C  REFERENCES TO WORK AND WORK2
C  REFER TO STARTING ADDRESS
C
C  LENGTH OF COMPRESSES DATA = NCNLV
C
C  WORK(1,1 TO NLEV)        = DTHBYDT(#,1 TO NLEV)
C  WORK(1,NLEV+1 TO 2*NLEV) = DQBYDT(#,1 TO NLEV)
C  WORK2(1,1 TO NLEV+1)     = EXNER(#,1 TO NLEV+1)
C  WORK2(1,NLEV+2)          = TH(#,1)
C  WORK2(1,NLEV+3)          = PSTAR(#)
C----------------------------------------------------------------------
C
      IF (NCNLV .NE. 0)THEN
C
        CALL COR_ENGY (NP_FIELD,NPNTS,NCNLV,NLEV,DTHBYDT,DQBYDT,SNOW,  
     *                EXNER,PSTAR,DELAK,DELBK,AKM12,BKM12,INDEX4)      
CL
      IF (L_3D_CCA) THEN
        CALL CALC_3D_CCA(NP_FIELD,NPNTS,NLEV,NBL,ANVIL_FACTOR
     &                  ,TOWER_FACTOR,AKM12,BKM12,ICCB,ICCT
     &                  ,FREEZE_LEV,PSTAR,CCA_2D,CCA,L_CLOUD_DEEP)
      ELSE
        DO I=1,NPNTS
          CCA(I,1)=CCA_2D(I)
        ENDDO
      ENDIF
CL---------------------------------------------------------------------
CL  UPDATE MODEL POTENTIAL TEMPERATURE, MIXING RATIO, U, V
CL  AND TRACER WITH INCREMENTS DUE TO CONVECTION
CL---------------------------------------------------------------------
CL
        DO 250 K=1,NLEV
         DO 250 I=1,NPNTS
           TH(I,K) = TH(I,K) + DTHBYDT(I,K) * TIMESTEP
           Q(I,K) = Q(I,K) + DQBYDT(I,K) * TIMESTEP
C
C---------------------------------------------------------------------
C   Calculation of gridbox mean CCW and CCWP, and CCA x conv. cloud
C   base and top pressure.
C---------------------------------------------------------------------
C
           IF (CCA_2D(I) .NE. 0.0) THEN
             IF (L_3D_CCA) THEN
               GBMCCW(I,K) = CCA(I,K) * CCW(I,K)
               DELPK(I) = -DELAK(K) - DELBK(K)*PSTAR(I)
               GBMCCWP(I) = GBMCCWP(I) + CCW(I,K)*DELPK(I)*CCA(I,K)/G
                 IF (K.EQ.NLEV) THEN
                   ICCBPxCCA(I) = CCA(I,ICCB(I)) *
     *                      (AK(ICCB(I)) + BK(ICCB(I)) * PSTAR(I))
                   ICCTPxCCA(I) = CCA(I,ICCT(I)-1) *
     *                      (AK(ICCT(I)) + BK(ICCT(I)) * PSTAR(I))
                 ENDIF
             ELSE
               GBMCCW(I,K)  = CCA_2D(I) * CCW(I,K)
               IF (K.EQ.NLEV) THEN
                 GBMCCWP(I)   = CCA_2D(I) * CCLWP(I)
                 ICCBPxCCA(I) = CCA_2D(I) *
     *                        (AK(ICCB(I)) + BK(ICCB(I)) * PSTAR(I))
                 ICCTPxCCA(I) = CCA_2D(I) *
     *                        (AK(ICCT(I)) + BK(ICCT(I)) * PSTAR(I))
               END IF
             ENDIF
           ENDIF
  250   CONTINUE
C
        IF(L_TRACER)THEN
C
!
!      BEFORE UPDATING THE TRACER FIELD, ADJUST THE TIMESTEP TO
!      PREVENT ANY NEGATIVE VALUES INVADING THE TRACER FIELDS.
!      NOTE THAT THE ADJUSTED TIMESTEP IS A FUNCTION OF GEOGRAPHICAL
!      LOCATION AND THE PARTICULAR TRACER.
!
        DO KTRA=1,NTRA
!
          DO I=1,NPNTS
!
            DO K=1,NLEV
!
              STEP_TEST2(K) = DTRABYDT(I,K,KTRA)
              STEP_TEST1(K) = ( 0.9999*ABS(TRACER(I,K,KTRA)) ) /
     &                  ( ABS(STEP_TEST2(K)) + SAFETY_MARGIN )
!
            END DO         ! END OF LEVEL (K) LOOP.
!
!
!     The following fragment of code provides a standard Fortran
!     alternative to the use of the Cray MINVAL function.
!
           LIMITED_STEP(I) = TIMESTEP
            DO K = 1,NLEV
              IF( STEP_TEST2(K) .LT. 0.0 ) THEN
               IF ( STEP_TEST1(K) .LT. LIMITED_STEP(I) ) THEN
               LIMITED_STEP(I) = STEP_TEST1(K)
               ENDIF
              ENDIF
            END DO
!
!     End of alternative to MINVAL.
!
!
!   Diagnose the factor by which the tstep has been multiplied
!
            REDUCTION_FACTOR(I,KTRA) = LIMITED_STEP(I)/TIMESTEP
!
          END DO           ! END OF LOCATION (I) LOOP.
!
!     Now update tracer field using  LIMITED STEP.
!     We can reverse order of  I and K loop now.
!
          DO K=1,NLEV
            DO I=1,NPNTS
              TRACER(I,K,KTRA) = TRACER(I,K,KTRA) + DTRABYDT(I,K,KTRA)
     &                           * LIMITED_STEP(I)
            END DO
          END DO
!
        END DO             ! END OF LOOP OVER TRACER TYPES (KTRA).
!
C
        END IF
C
      END IF
C
      RETURN
      END
C
