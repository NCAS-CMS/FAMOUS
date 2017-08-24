C ******************************COPYRIGHT******************************
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
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  3.2  8/07/93 : added convective cloud condensed water diagnostic
CLL               : P Inness
CLL  3.4  11/03/94  Add lowest conv.cloud diagnostics.  R.T.H.Barnes.
CLL
CLL  3.4  21/09/94  Standard deviation of surface layer turbulent
CLL                 fluctuations of temperature and humidity input
CLL                                                     R.N.B.Smith.
CLL
CLL
CLL  3.4  06/08/94  BTERM initialised to make downdraught scheme
CLL                 reproducible when call to CONVECT macrotasked.
CLL                 Authors: A.Dickinson, D.Salmond, Reviewer: R.Barnes
CLL
CLL  4.0  5/05/95   Added CAPE diagnostic. Pete Inness.
CLL
CLL  4.0  5/05/95   References to surface fluxes removed as not
CLL                 required at this version. Pete Inness.
CLL
CLL   4.2    Oct. 96  T3E migration: *DEF CRAY removed 
CLL                   (was used to switch on WHENIMD)
CLL                                    S.J.Swarbrick
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
!LL  4.4  Oct 97    Add halo mask to stop redundant calculations
!LL                                               Alan Dickinson
CLL  4.5  Jul. 98   Kill the IBM specific lines (JCThil)
!LL  4.5  20/02/98  Remove redundant code. A. Dickinson               
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
      SUBROUTINE CONVECT (NP_FIELD,NPNTS,NLEV,TH,Q,PSTAR,BLAND,DTHBYDT,
     &                    DQBYDT,RAIN,SNOW,CCA,ICCB,ICCT,CCLWP,
     &                    CCW,ICCBPxCCA,ICCTPxCCA,GBMCCWP,GBMCCW,
     &                    LCBASE,LCTOP,LCCA,LCCLWP,CAPE_OUT,
     &                    EXNER,AK,BK,AKM12,BKM12,DELAK,DELBK,TIMESTEP
     &                    ,l_halo
     &                    )
C
      IMPLICIT NONE
C
C
C--------------------------------------------------------------------
C MODEL CONSTANTS
C--------------------------------------------------------------------
C
      REAL THPIXS, QPIXS ! INITIAL EXCESS POTENTIAL TEMPERATURE (K)
                         ! AND MIXING RATIO (KG/KG)
       PARAMETER ( THPIXS = 0.2, QPIXS = 0.0 )
C
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
      REAL C,D  !  CONSTANTS USED TO DETERMINE THE INITIAL CONVECTIVE
                !  MASS FLUX FROM PARCEL BUOYANCY
                ! MI = 1.E-3 * (D + C*BUOYANCY/DELP)
       PARAMETER ( C = 5.17E-4, D = 0.0 )
C
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
      INTEGER I,K,KC              ! LOOP COUNTERS
C
C
C---------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C---------------------------------------------------------------------
C
      LOGICAL BLAND(NP_FIELD)     ! IN LAND/SEA MASK
C
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
      LOGICAL l_halo(NP_FIELD)  ! IN:  Mask for halos
C
C
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
      REAL RAIN(NP_FIELD)         ! OUT SURFACE CONVECTIVE RAINFALL
                                  !     (KG/M**2/S)
C
      REAL SNOW(NP_FIELD)         ! OUT SURFACE CONVECTIVE SNOWFALL
                                  !     (KG/M**2/S)
C
      REAL CCA(NP_FIELD)          ! OUT CONVECTIVE CLOUD AMOUNT (%)
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
      REAL CAPE_OUT(NPNTS)        ! OUT SAVED VALUES OF CONVECTIVE
                                  !     AVAILABLE POTENTIAL ENERGY
                                  !     FOR DIAGNOSTIC OUTPUT (J/KG)
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
C----------------------------------------------------------------------
C VARIABLES DEFINED LOCALLY
C
      REAL WORK(NPNTS,NLEV*2),   ! WORK SPACE
     *     WORK2(NPNTS,NLEV*2)
      LOGICAL BWORK(NPNTS,4),    ! WORK SPACE FOR 'BIT' MASKS
     *        BWORK2(NPNTS,4)
C
      REAL CAPE(NPNTS)            ! CONVECTIVE AVAILABLE POTENTIAL
                                  ! ENERGY (J/KG)
C
      REAL CAPE_C(NPNTS)          ! CAPE - COMPRESSED
C
C
      LOGICAL BCONV(NPNTS)       ! MASK FOR POINTS WHERE STABILITY
                                  ! LOW ENOUGH FOR CONVECTION
                                  ! TO OCCUR
C
      REAL QSE(NPNTS,NLEV)       ! SATURATION MIXING RATIO OF CLOUD
                                  ! ENVIRONMENT (KG/KG)
C
      REAL TT(NPNTS)             ! TEMPORARY STORE FOR TEMPERATURE
                                  ! IN CALCULATION OF SATURATION
                                  ! MIXING RATIO (K)
C
      REAL PT(NPNTS)             ! TEMPORARY STORE FOR PRESSURE
                                  ! IN CALCULATION OF SATURATION
                                  ! MIXING RATIO (PA)
C
      REAL CCAC(NPNTS)            ! COMPRESSED VALUES OF CCA
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
      REAL PRECIP(NPNTS,NLEV)     ! AMOUNT OF PRECIPITATION
                                  ! FROM EACH LAYER (KG/M*:2/S)
C
      REAL THPI(NPNTS)            ! INITIAL PARCEL POTENTIAL TEMPERATURE
                                  ! (K)
C
      REAL QPI(NPNTS)             ! INITIAL PARCEL MIXING RATIO
                                  ! (KG/KG)
C
      REAL THP(NPNTS,NLEV)        ! PARCEL POTENTIAL TEMPERATURE
                                  ! IN LAYER K (K)
C
      REAL QP(NPNTS,NLEV)         ! PARCEL MIXING RATIO IN LAYER K
                                  ! (KG/KG)
C
      REAL XPK(NPNTS)             ! PARCEL CLOUD WATER IN LAYER K
                                  ! (KG/KG)
C
      REAL FLX(NPNTS,NLEV)        ! PARCEL MASSFLUX IN LAYER K (PA/S)
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
      REAL EXK(NPNTS)             ! EXNER RATIO AT LEVEL K
C
      REAL EXKP1(NPNTS)           ! EXNER RATIO AT LEVEL K+1
C
      REAL DELEXKP1(NPNTS)        ! DIFFERENCE IN EXNER RATIO
                                  ! ACROSS LAYER K+1
C
      REAL EMINDS(NPNTS)          !
C
      INTEGER INDEX1(NPNTS),      ! INDEX FOR COMPRESS AND
     *        INDEX2(NPNTS),      ! EXPAND
     *        INDEX3(NPNTS),
     *        INDEX4(NPNTS)
C
C
      REAL recip_PSTAR(NP_FIELD)  ! Reciprocal of pstar array 
      REAL FLX2                   ! TEMPORARY STORE FOR MASS FLUX
C
C----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C----------------------------------------------------------------------
C
      EXTERNAL QSAT,FLAG_WET,LIFT_PAR,CONVEC2,LAYER_CN,
     *         DQS_DTH,COR_ENGY,DD_CALL
C

      REAL
     &    PU,PL
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
        PU=PSTAR(I)*BKM12(K+1) + AKM12(K+1)
        PL=PSTAR(I)*BKM12(K) + AKM12(K)
        TT(I) = TH(I,K)* P_EXNER_C(EXNER(I,K+1),EXNER(I,K),PU,PL,KAPPA)
        PT(I) = AK(K)+BK(K)*PSTAR(I)
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
C INITIALISE PRECIPITATION, DTH/DT, DQ/DT, CCW ARRAYS
C----------------------------------------------------------------------
C
      DO 40 K=1,NLEV
       DO 40 I=1,NPNTS
        PRECIP(I,K) = 0.0
        CCW(I,K) = 0.0
        GBMCCW(I,K) = 0.0
        DTHBYDT(I,K) = 0.0
   40   DQBYDT(I,K) = 0.0
      DO 50 I=1,NPNTS
C
C----------------------------------------------------------------------
C INITIALISE BIT VECTORS FOR POINTS WHICH ARE ALREADY CONVECTING
C AND FOR POINTS AT WHICH CONVECTION OCCURS AT SOME LEVEL OF
C THE ATMOSPHERE
C----------------------------------------------------------------------
C
        BINIT(I) = .FALSE.
        BCNLV(I) = .FALSE.
        BTERM(I) = .FALSE.
C
C----------------------------------------------------------------------
C INITIALISE RADIATION DIAGNOSTICS
C----------------------------------------------------------------------
C
       CCA(I) = 0.0
       ICCB(I) = 0
       ICCT(I) = 0
       TCW(I) = 0.0
       CCLWP(I) = 0.0
C
C--------------------------------------------------------------------
C INITIALISE GRIDBOX MEAN DIAGNOSTICS
C--------------------------------------------------------------------
C
       GBMCCWP(I) = 0.0
       ICCBPxCCA(I) = 0.0
       ICCTPxCCA(I) = 0.0
C
CL-------------------------------------------------------------------
CL INITIALISE CAPE DIAGNOSTIC
CL-------------------------------------------------------------------
C
       CAPE(I) = 0.0
       CAPE_OUT(I) = 0.0
       CAPE_C(I) = 0.0
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
       EMINDS(I) = MPARB*DELPKP12(I)*recip_pstar(I)
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
         THPI(I) = TH(I,K) + THPIXS
         QPI(I) = Q(I,K) + QPIXS
         THP(I,K) = TH(I,K) + THPIXS
         QP(I,K) = Q(I,K) + QPIXS
         XPK(I) = 0.0
         FLX(I,K) = 0.0
         BGMK(I) = .FALSE.
         DEPTH(I) = 0.0
       END IF
CL
CL----------------------------------------------------------------------
CL FORM A BIT VECTOR OF POINTS FOR WHICH CONVECTION MAY BE POSSIBLE
CL FROM LAYER K TO K-1 EITHER BECAUSE STABILITY IS LOW ENOUGH
CL OR BECAUSE CONVECTION OCCURRING FROM LAYER K+1 TO K
CL THIS BIT VECTOR IS USED IN THE FIRST COMPREE OF THE DATA
CL TO CALCULATE PARCEL BUOYANCY IN LAYER K-1
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
C
C
  80    CONTINUE
CL
CL---------------------------------------------------------------------
CL LIFT PARCEL FROM LAYER K TO K-1
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (3) AND (4)
CL---------------------------------------------------------------------
CL
      CALL LIFT_PAR (NCONV,WORK(1,13),WORK(1,14),WORK(1,15),
     *               BWORK(1,2),BWORK(1,1),WORK(1,7),WORK(1,8),
     *               WORK(1,2),WORK(1,4),WORK(1,1),WORK(1,3),
     *               WORK(1,5),WORK(1,6),WORK(1,9),
     *               WORK(1,10),WORK(1,11),WORK(1,12))
C
      DO 110 I=1,NCONV
CL
CL---------------------------------------------------------------------
CL CALCULATE BUOYANCY OF PARCEL IN LAYER K-1
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
      IF (BWORK(I,3))
     1     WORK(I,19) = 1.0E-3 * WORK(I,18) * (D + C * WORK(I,18) *
     2                    ((WORK(I,16) - XSBMIN) / WORK(I,17)))
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
      CALL LIFT_PAR (NCONV,WORK(1,13),WORK(1,14),WORK(1,15),
     *               BWORK(1,2),BWORK(1,1),WORK(1,7),WORK(1,8),
     *               WORK(1,2),WORK(1,4),WORK(1,1),WORK(1,3),
     *               WORK(1,5),WORK(1,6),WORK(1,9),
     *               WORK(1,10),WORK(1,11),WORK(1,12))
C
       DO I=1,NCONV
C
        IF ( BWORK(I,3) ) THEN
CL
CL---------------------------------------------------------------------
CL RECALCULATE BUOYANCY OF PARCEL IN LAYER K-1
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
C
       END DO
C
      END IF
C
CL
CL--------------------------------------------------------------------
CL ZERO MIXING DETRAINMENT RATE WHEN CONVECTION STARTS FROM LAYER K
CL--------------------------------------------------------------------
CL
      DO I=1,NCONV
       IF ( BWORK(I,3) ) AMDETK(INDEX1(I))=0.0
      END DO
CL
CL--------------------------------------------------------------------
CL COMPRESS DOWN THOSE POINTS WHICH ARE NOT BUOYANT IN LAYER K-1
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
C
C  WORK AND BWORK NOW CONTAIN DATA COMPRESSED DOWN
C  FROM FULL LENGTH VECTORS
C
C  WORK(1,1) = not used in this section
C  WORK(1,2) = QSE(#,K)
C  WORK(1,3) = DQSTHK(#)
C  WORK(1,4) = THPI(#)
C  WORK(1,5) = QPI(#)
C  WORK(1,6) = XPK(#)
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
  120   CONTINUE
C
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
        WORK(I,6) = XPK(INDEX3(I))
        WORK(I,8) = DEPTH(INDEX3(I))
        CCAC(I)    = CCA(INDEX3(I))
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
        CAPE_C(I)  = CAPE(INDEX3(I))
C
        BWORK(I,4) = .TRUE.
  140 CONTINUE
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
         CALL CONVEC2 (NINIT,NLEV,K,WORK2(1,1),WORK2(1,2),WORK2(1,3),
     *                WORK2(1,4),WORK2(1,5),WORK2(1,6),WORK2(1,18),
     *                WORK2(1,7),WORK2(1,8),WORK2(1,13),WORK2(1,14),
     *                WORK2(1,15),WORK2(1,16),WORK(1,2),WORK(1,3),
     *                WORK(1,4),WORK(1,5),WORK(1,6),WORK2(1,19),
     *                BWORK2(1,1),BWORK2(1,2),BWORK(1,1),BWORK2(1,3),
     *                BWORK(1,2),BWORK(1,3),WORK(1,8),WORK(1,9),
     *                WORK(1,10),WORK(1,11),WORK(1,12),WORK(1,13),
     *                BWORK(1,4),CCAC,ICCBC,ICCTC,TCWC,
     *                WORK2(1,11),WORK2(1,12),WORK(1,14),
     *                WORK(1,16),WORK2(1,9),WORK(1,17),WORK2(1,10),
     *                WORK(1,18),WORK(1,19),WORK(1,20),
     *             CCLWPC,WORK(1,21),LCCAC,LCBASEC,LCTOPC,LCCLWPC,
     *             CAPE_C)
CL
CL---------------------------------------------------------------------
CL EXPAND REQUIRED VECTORS BACK TO FULL FIELDS
CL----------------------------------------------------------------------
CL
      DO 145 I=1,NPNTS
        THP(I,K+1) = 0.0
        QP(I,K+1) = 0.0
        XPK(I) = 0.0
        FLX(I,K+1)= 0.0
        DEPTH(I) = 0.0
        PRECIP(I,K+1) = 0.0
        BGMK(I) = .FALSE.
        BTERM(I) = .FALSE.
        BINIT(I) = .FALSE.
  145 CONTINUE
C
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO 150 I=1,NINIT
        THP(INDEX3(I),K+1) = WORK2(I,7)
        QP(INDEX3(I),K+1) = WORK2(I,8)
        XPK(INDEX3(I)) = WORK(I,6)
        FLX(INDEX3(I),K+1) = WORK2(I,19)
        DEPTH(INDEX3(I)) = WORK(I,8)
        PRECIP(INDEX3(I),K+1) = WORK(I,9)
        DTHBYDT(INDEX3(I),K) = WORK(I,10)
        DQBYDT(INDEX3(I),K) = WORK(I,11)
        DTHBYDT(INDEX3(I),K+1) = WORK(I,12)
        DQBYDT(INDEX3(I),K+1) = WORK(I,13)
        CCA(INDEX3(I)) = CCAC(I)
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
C
        BGMK(INDEX3(I)) = BWORK(I,1)
        BTERM(INDEX3(I)) = BWORK(I,3)
        BINIT(INDEX3(I)) = BWORK(I,4)
  150 CONTINUE
C
       END IF
C
      END IF
C-------------------------------------------------------------------
C IF CONVECTION IS TERMINATING, READ VALUE OF CAPE INTO DIAGNOSTIC
C OUTPUT ARRAY AND RESET TO ZERO
C-------------------------------------------------------------------
C
      DO I=1,NPNTS
       IF(BTERM(I))THEN
        CAPE_OUT(I)=CAPE(I)
        CAPE(I)=0.0
       END IF
      END DO
C
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
        END IF
  160 CONTINUE
C
      IF (NTERM .NE. 0) THEN
C
         CALL DD_CALL (NP_FIELD,NPNTS,K,THP(1,1),QP(1,1),TH(1,1),Q(1,1),
     *                 DTHBYDT(1,1),DQBYDT(1,1),FLX(1,1),PSTAR,
     *                 AK,BK,AKM12,BKM12,DELAK,DELBK,EXNER(1,1),
     *                 PRECIP(1,1),RAIN,SNOW,ICCB,ICCT,BWATER(1,2),
     *                 BTERM,BGMK,TIMESTEP,CCA,NTERM,recip_pstar)     
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
          CCA(I) = 0.0
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
C
CL
CL---------------------------------------------------------------------
CL  UPDATE MODEL POTENTIAL TEMPERATURE AND MIXING RATIO
CL  WITH INCREMENTS DUE TO CONVECTION
CL---------------------------------------------------------------------
CL
        DO 250 K=1,NLEV
         DO 250 I=1,NPNTS
           TH(I,K) = TH(I,K) + DTHBYDT(I,K) * TIMESTEP
           Q(I,K) = Q(I,K) + DQBYDT(I,K) * TIMESTEP
CL
CL---------------------------------------------------------------------
CL CALCULATE GRIDBOX MEAN DIAGNOSTICS
CL---------------------------------------------------------------------
CL
           IF (CCA(I) .NE. 0.0) THEN
             GBMCCW(I,K)  = CCA(I) * CCW(I,K)
             IF (K.EQ.NLEV) THEN
               GBMCCWP(I)   = CCA(I) * CCLWP(I)
               ICCBPxCCA(I) = CCA(I) *
     *                      (AK(ICCB(I)) + BK(ICCB(I)) * PSTAR(I))
               ICCTPxCCA(I) = CCA(I) *
     *                      (AK(ICCT(I)) + BK(ICCT(I)) * PSTAR(I))
             END IF
           ENDIF
  250   CONTINUE
C
      END IF
C
      RETURN
      END
C
