!    Subroutine SULPHUR ---------------------------------------------
!
! Purpose: To perform oxidation chemistry of sulphur dioxide to 3 modes
!          of Sulphate aerosol (Aitken, accumulation and dissolved),
!          and dimethyl sulphide to sulphur dioxide and methyl sulphonic
!          acid.
!          There is also exchange between the 3 modes of slphateaerosol
!          due to nucleation, diffusion and evaporation processes.
!          Called by CHEM_CTL
!
! Current owners of code:                M J Woodage
!
! History:
! Version     Date     Comment
! -------     ----     -------
!  4.1      10/06/96   Original code       M J Woodage
!                      + Additional code   R Colvile, D Roberts
!
!  4.4      30/09/97    Partition dry oxidised SO2 so that some becomes
!                       Aitken mode, and some accumulation mode SO4.
!
!                       Replenish H2O2 from HO2 in daylight only.
!
!                       Evaporate dissolved SO4 circulating through
!                       non_cloudy part of grid box when cloud fraction
!                       less than 0.95.
!                                                 (M. Woodage)
!  4.5    1/06/98   Add code for oxidation by ozone using ammonia
!                   as a buffer, to wet oxidation section.     MJW
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
! System components covered:
!
! System task:
!
! Documentation: Not yet available
!
!-----------------------------------------------------------------
!
      SUBROUTINE SULPHUR(SO2,
     &             SO4_AIT,SO4_ACC,SO4_DIS,
     &            NH3,NH3_DEP,O3,L_SULPC_OZONE,
     &             DMS,H2O2_MXR,CLOUDF,
     &             NPNTS,QPNTS,NPFLD,TSTEP,
     &             QCL,QCF,LAND_MASK,COSZA2D,
     &             P,T,Q,
     &             OH_CONC,H2O2_LMT,HO2_CONC,
     &             FIRST_POINT,LAST_POINT,
     &             MSA,L_SULPC_DMS)
!
!
      IMPLICIT NONE
!
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

!
! Description:
!
!  Contains various parameters for the SW effective radius
!  parametrisation, defined for land and sea areas.
!
!  NTOT is the number concentration (m-3) of activated CCN;
!  KPARAM is the ratio of volume mean radius to effective radius;
!  DCONRE is the effective radius (m) for deep convective clouds.
!
! Current Code Owner: Andy Jones
!
! History:
!
! Version   Date     Comment
! -------   ----     -------
!    1     040894   Original code.    Andy Jones
!
      REAL NTOT_LAND,
     &     NTOT_SEA,
     &     KPARAM_LAND,
     &     KPARAM_SEA,
     &     DCONRE_LAND,
     &     DCONRE_SEA

      PARAMETER ( NTOT_LAND = 6.0E08,
     &            NTOT_SEA = 1.5E08,
     &            KPARAM_LAND = 0.67,
     &            KPARAM_SEA = 0.80,
     &            DCONRE_LAND = 9.5E-06,
     &            DCONRE_SEA = 13.5E-06 )
C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------
!-------------------COMDECK C_SULCHM--------------------------------
! Parameters for Sulphur Cycle Chemistry
      REAL
     &     EVAPTAU,       ! timescale for dissolved SO4 to evaporate
     &     NUCTAU,        ! timescale for accumulation mode particles
!                           to nucleate once they enter a cloud.
     &     DIFFUSE_AIT,   ! diffusion coefficient of Aitken particles
     &     K_SO2OH_HI,            ! high pressure reaction rate limit
     &     K_DMS_OH,              ! reaction rate for DMS+OH  cc/mcl/s
     &     BRAT_SO2,              ! branching ratio for SO2 in DMS oxidn
     &     BRAT_MSA,              ! branching ratio for MSA in DMS oxidn
     &     AVOGADRO,             ! no. of molecules in 1 mole
     &     RMM_H2O2,             ! relative molecular mass H2O2 kg/mole
     &     RMM_AIR,              ! relative molecular mass dry air
     &     RMM_W,                ! relative molecular mass water
     &     RELM_S_H2O2,          ! rel atomic mass sulphur/RMM_H2O2
     &     RELM_S_2N,         ! rel atomic mass Sulphur/2*Nitrogen
     &     PARH,                ! power of temp dependence of K_SO2OH_LO
     &     K1,                  ! parameters for calcn of K_SO2OH_LO
     &     T1,                  !
     &     FC,                  ! parameters for interpolation between
     &     FAC1,                !   LO and HI reaction rate limits
     &     K2,K3,K4,            ! parameters for calcn of K_HO2_HO2
     &     T2,T3,T4,            !
     &     CLOUDTAU,              ! air parcel lifetime in cloud
     &     CHEMTAU,               ! chem lifetime in cloud before oxidn
     &     O3_MIN,            ! min mmr of O3 required for oxidn
     &     THOLD                  ! threshold for cloud liquid water
!
!
      PARAMETER (
     &           EVAPTAU = 300.0,             ! secs  (=5 mins) 
     &             NUCTAU = 30.0,         ! secs
     &       DIFFUSE_AIT = 1.7134E-9,        ! sq m/s
     &        K_SO2OH_HI = 1.5E-12,    ! cc/mcl/s from STOCHEM model
     &           K_DMS_OH = 9.1E-12,      ! cc/mcl/s
     &          BRAT_SO2 = 0.9,   
     &           BRAT_MSA = 1.0-BRAT_SO2,
     &           AVOGADRO = 6.022E23,     ! per mole
     &           RMM_H2O2 = 3.40E-2,      ! kg/mole
     &            RMM_AIR = 2.896E-2,     ! kg/mole
     &              RMM_W = 1.8E-2,       ! kg/mole
     &        RELM_S_H2O2 = 3.206/3.40,
     &           RELM_S_2N = 3.206/2.80,
     &               PARH = 3.3,
     &                K1 = 3.0E-31,    ! (cc/mcl)2/s from STOCHEM
     &                 T1 = 300.0,        ! K
     &                FC = 0.6,        ! from STOCHEM model
     &              FAC1 = 1.0317,  ! 0.75-1.27*LOG10(FC) from STOCHEM
     &                 K2 = 2.3E-13,      ! cc/mcl/s
     &                 K3 = 1.9E-33,      ! (cc/mcl)2/s
     &                 K4 = 1.4E-21,      ! cc/mcl
     &                 T2 = 600.0,        ! K
     &                 T3 = 890.0,        ! K
     &                 T4 = 2200.0,       ! K
     &           CLOUDTAU = 1.08E4,       ! secs (=3 hours)
     &            CHEMTAU = 9.0E2,        ! secs (=15 mins)
     &              O3_MIN = 1.6E-8,    !(kg/kg, equiv. 10ppbv)
     &              THOLD = 1.0E-8        ! kg/kg
     &          )
!
      REAL RAD_AIT,         ! median radius of Aitken mode particles
     &     DIAM_AIT,        !   "    diameter    " 
     &     RAD_ACC,         ! median radius of acccumulation mode
     &     DIAM_ACC,        !   "    diameter    "
     &     CHI,             ! mole fraction of S in particle
     &     RHO_SO4,         ! density of  SO4 particle
     &     SIGMA,           ! standard devn of particle size distn
     &     E_PARM,          ! param relating size distns of Ait & Acc
     &     NUM_STAR         ! threshold concn of accu mode particles
                            !  below which PSI=1
!
      PARAMETER (
     &           RAD_AIT = 24.0E-9,          ! m
     &          DIAM_AIT = 2.0*RAD_AIT,    
     &           RAD_ACC = 95.0E-9,          ! m
     &          DIAM_ACC = 2.0*RAD_ACC,    
     &               CHI = 32.0/132.0,
     &           RHO_SO4 = 1769.0,            ! kg/m3
     &             SIGMA = 1.4,
     &            E_PARM = 0.9398, 
     &          NUM_STAR = 1.0E6             ! m-3
     &          )
!
!*---------------------------------------------------------------------
C*L------------------COMDECK C_PI---------------------------------------
CLL
CLL 4.0 19/09/95  New value for PI. Old value incorrect
CLL               from 12th decimal place. D. Robinson
CLL
      REAL PI,PI_OVER_180,RECIP_PI_OVER_180

      PARAMETER(
     & PI=3.14159265358979323846, ! Pi
     & PI_OVER_180 =PI/180.0,   ! Conversion factor degrees to radians
     & RECIP_PI_OVER_180 = 180.0/PI ! Conversion factor radians to
     &                              ! degrees
     & )
C*----------------------------------------------------------------------

!
      INTEGER
     &        NPNTS,            !IN no. of pts in 3_D array on P_LEVS
     &        QPNTS,            !IN no. of pts in 3_D array on Q_LEVS
     &        NPFLD,            !IN no. of pts in 2_D field
     &        FIRST_POINT,      !IN first point for calcns to be done
     &        LAST_POINT        !IN last  point for calcns to be done
!
      REAL TSTEP                ! IN chemistry tstep: may be less than
!                                     physics tstep
!
      REAL SO2(NPNTS),      !INOUT mass mix rat of sulphur in SO2 gas
     &     DMS(NPNTS),      !INOUT mass mix rat of S in DMS gas
     &     SO4_AIT(NPNTS),  !INOUT mass mix rat of S in SO4 AITKEN mode
     &     SO4_ACC(NPNTS),  !INOUT mass mix rat of S in SO4 ACCUM  mods
     &     SO4_DIS(NPNTS),  !INOUT mass mix rat of S in dissolved SO4
     &     NH3(NPNTS),        !INOUT mass mix rat of N in NH3
     &     H2O2_MXR(NPNTS)  !INOUT mass mix rat of hydrogen peroxide
!
      REAL CLOUDF(QPNTS),      !IN  cloud fraction (range 0 TO 1)
     &     QCL(QPNTS),         !IN  cloud liquid water (mmr)
     &     QCF(QPNTS),         !IN  cloud frozen water (mmr)
     &     COSZA2D(NPFLD),     !IN  cos(zenith angle) 2_D field
     &     P(NPNTS),           !IN  pressure at model level
     &     T(NPNTS),           !IN  temperature
     &     Q(QPNTS)            !IN  specific humidity (kg/kg)
!
      REAL OH_CONC(NPNTS),     !IN  OH concentration (3_D array)
     &     HO2_CONC(NPNTS),    !IN  HO2 concn (3_D array)
     &     H2O2_LMT(NPNTS)     !IN  H2O2 max limit field (3_D array)
      REAL O3(NPNTS)          !IN O3 mass mix rat field (3_D array)
      REAL NH3_DEP(NPNTS)     !OUT NH3 mmr depleted by buffering
      REAL MSA(NPNTS)          !OUT mass mix rat of S in MSA
!
      LOGICAL L_SULPC_DMS,        !IN, TRUE if DMS chemistry required
     &        L_SULPC_OZONE,     !IN, TRUE if O3 oxidn required
     &        LAND_MASK(NPFLD)    !IN, TRUE IF LAND, FALSE IF SEA.
!
! Local variables
!
      INTEGER I,J,K               ! loop variables
!
      REAL
     &     DRYRATE,               ! rate of dry oxidn SO2  in 1/s
     &     WETRATE,               ! rate of wet oxidn SO2  in 1/s
     &     DMSRATE                ! rate of dry oxidn DMS in 1/s
!
      REAL CLEARF(QPNTS),         ! clear air fraction (1-CLOUDF)
     &     DELTAS_DRY(NPNTS),     ! amount SO2 dry oxidised per ts
     &     DELTAS_WET(NPNTS),     ! amount SO2 wet oxidised per ts
     &     DELTAS_WET_O3(NPNTS), !amount of SO2 oxidised by O3 per ts 
     &     DELTAS_TOT(NPNTS),     ! total amount SO2 oxidised per ts
     &     DELTA_DMS(NPNTS),      ! amount DMS dry oxidised per ts
     &     DELTAS_EVAP(NPNTS), ! amount of dissolved SO4 released due
!                                to the evapn of cloud droplets per ts
     &     DELTAS_NUCL(NPNTS), ! amount of accu mode SO4 transfered to
!                                dissolved SO4 by nucleation per ts
     &     DELTAS_DIFFUSE(NPNTS), ! amount of Aitken mode SO4 transfd
!                                   to dissolved SO4 by diffusion per ts
     &     QCTOTAL(QPNTS)      ! total condensed water amount.(QCL+QCF)
!
      REAL
     &     EVAPTIME,        ! timescale for cloud droplets to evaporate
     &     NUCTIME,        ! timescale for particles to enter a cloud
!                            and nucleate.
     &     DIFFUSE_TAU,    ! diffusive lifetime of Aitken mode particles
!                            once they enter a cloud
     &     NLAND23,        ! NTOT_LAND**(-2/3)
     &     NSEA23,         ! NTOT_SEA**(-2/3)
     &     RHO_CUBEROOT,   ! cube root of density of water.
     &     DIFFUSE_TAURATIO, ! CLOUDTAU/DIFFUSE_TAU
     &     PROBDIFF_INV,   ! inverse of probability of a particle being
!                            diffusionallly scavenged in a cloud.
     &     PROBDIFF_FN1,   ! PROBDIFF_INV - 0.5
     &     PROBDIFF_FN2,   ! PROBDIFF_INV*EXP(DIFFUSE_TAURATIO*0.5)
     &     PROBDIFF_CLOUD, ! probability of an Aitken mode particle
!                            being in cloud at the start of a step.
     &     PROBDIFF_CLEAR, ! probability of an Aitken mode particle
!                            being in clear air at the start of a step.
     &     LAMBDA_AITKEN,  ! ratio of concentrations of Aitken mode
!                            particles in cloudy to clear air.
     &     DIFFRATE ! rate of diffusive capture of Aitken mode particles
!
      REAL AIR_CONC,              ! air concentration on molecules/cc
     &     W_CONC,                ! H2O concentration in molecules/cc
     &     K_SO2_OH,              ! reaction rate for SO2+OH  cc/mcl/s
     &     K_SO2OH_LO,            ! low_pressure reaction rate limit
     &     K_LIMRAT,              ! ratio  reaction rate limits (LO/HI)
     &     K_HO2_HO2,             ! reaction rate for HO2 to H2O2
     &     H2O2_CON_TO_MXR        ! conversion factor =
!                                 ! RMM_H2O2/AVAGADRO/AIR DENSITY
      REAL FAC2,FAC3              ! for interp between LO and HI K_SO2OH
!                                 !    reaction rate limits
!
      REAL TERM1,                 ! dummy variables to assist
     &     TERM2,                 !  wet oxidn calculation
     &     DENOM,                 !
     &     TERM3,                 ! dummy variables to assist
     &     TERM4,                 ! calculation of diffusive capture.
     &     TAURATIO,              ! CLOUDTAU/CHEMTAU
     &     HFTAURAT,              ! CLOUDTAU/2*CHEMTAU
     &     PROBOX_INV,            ! 1/probability of oxidn in cloud
     &     PROBOX_FN1,            ! PROBOX_INV - 0.5
     &     PROBOX_FN2,            ! PROBOX_INV * EXP(-HFTAURAT)
     &     PROB_CLOUD,            ! probability SO2 starts in cloud
     &     PROB_CLEAR,            ! probability SO2 starts in clear air
     &     LAMDA,                 ! ratio SO2 in cloud/clear air
     &     H2O2_REP               ! replenished H2O2
!
      REAL PSI(NPNTS)      ! fraction of dry oxidised SO2 becoming Ait
                           !  mode SO4 (remainder becomes accu mode)
      REAL MMR_STAR        ! threshold mmr of SO4_ACC below which PSI=1
      REAL RHO_AIR         ! air density kg/m3
      REAL CON_1,          ! constants required for PSI calcn
     &     CON_2
      REAL SCALE_FACTOR          ! to prevent negative SO2
!
!--------------------------------------------------------------------
! 0. Set up constants and initialise DELTA increments to 0
!--------------------------------------------------------------------
!
! The next few constants cannot be declared as PARAMETERs because
! they involve non-integer exponents.
!
      NLAND23 = 1.0/(NTOT_LAND**0.666667)
      NSEA23 = 1.0/(NTOT_SEA**0.666667)
      RHO_CUBEROOT = RHO_WATER**0.333333
!
!
! Calculate parameters which are same for each grid box (for wet oxidn)
!
      TAURATIO=CLOUDTAU/CHEMTAU
      HFTAURAT=TAURATIO/2.0
      PROBOX_INV=1.0/(1.0-EXP(-TAURATIO))
      PROBOX_FN1=PROBOX_INV-0.5
      PROBOX_FN2=PROBOX_INV*EXP(-HFTAURAT)
!
      DO J=1,NPNTS-NPFLD+1,NPFLD      ! J loops over all model points
        DO I=FIRST_POINT,LAST_POINT   ! I loop omits N AND S polar rows
          K=I+J-1
          DELTAS_DRY(K)  = 0.0
          DELTAS_WET(K)  = 0.0
          DELTAS_WET_O3(K) = 0.0
          NH3_DEP(K)  = 0.0
          DELTAS_TOT(K)  = 0.0
          DELTA_DMS(K)   = 0.0
          DELTAS_EVAP(K) = 0.0
          DELTAS_NUCL(K) = 0.0
          DELTAS_DIFFUSE(K) = 0.0
        END DO
      END DO
!
!---------------------------------------------------------------------
! 1.  This section calculates amount of SO2 dry oxidised using a
!     3_D OH concentration field to control the oxidn rate when
!     cos(zenith angle) is G.T. 10-6  (else rate is zero).
!  Reaction rates are dependent on temperature and pressure (R.Colvile)
!---------------------------------------------------------------------
!
      CON_1 = EXP(4.5*LOG(SIGMA)*LOG(SIGMA))
      CON_2 = 4.0*PI*RHO_SO4*CHI*NUM_STAR*(RAD_ACC)**3.0
      DO J=1,NPNTS-NPFLD+1,NPFLD      ! J loops over all model points
        DO I=FIRST_POINT,LAST_POINT   ! I loop omits N AND S polar rows
           K=I+J-1
!
         AIR_CONC = (AVOGADRO*P(K)) / (R*T(K)*RMM_AIR*1.0E6)  ! mcls/cc
         K_SO2OH_LO = AIR_CONC * K1 * (T1/T(K))**PARH
!
         K_LIMRAT = K_SO2OH_LO / K_SO2OH_HI
             FAC2 = K_SO2OH_LO /(1.0+ K_LIMRAT)
             FAC3 = 1.0 + (LOG10(K_LIMRAT)/FAC1)**2.0
!
         K_SO2_OH = FAC2 * FC**(1.0/FAC3)         ! Variable K_SO2_OH
!
!  Only calculate the dry oxidation rate in the daytime.
!
         IF (COSZA2D(I).GT.1.0E-06) THEN
           DRYRATE=K_SO2_OH * OH_CONC(K)
         ELSE                     ! Zero rate if nightime
           DRYRATE=0.0
         ENDIF
!
         DELTAS_DRY(K) = DRYRATE*TSTEP*SO2(K) ! Amnt SO2 dry oxidised
!  Calculate fraction becoming Aitken mode SO4
      RHO_AIR = P(K)/R*T(K)
      MMR_STAR = CON_1*CON_2/(3.0*RHO_AIR)
!
      IF (SO4_ACC(K).LT.MMR_STAR)   THEN
        PSI(K) = 1.0
      ELSE
        PSI(K) = RAD_ACC*E_PARM*SO4_AIT(K)/
     &          (RAD_ACC*E_PARM*SO4_AIT(K)+RAD_AIT*SO4_ACC(K))
      END IF
!
        END DO
      END DO
!
!--------------------------------------------------------------
! 2. This section calculates amount of SO2 wet oxidised using a
!    3_D H2O2 field to control the oxidn rate. As H2O2 is depleted, it
!    is replenished from the HO2 concn field (via a reaction rate), but
!    it may not exceed the H2O2 LIMIT field .
!    The reaction rate is a function of pressure, temp and humidity.
!    It uses D.Roberts' formula for the wet oxidn rate, and R.Colvile's
!    H2O2 and HO2 chemistry.
!    Depletion and replenishment of H2O2 is done at the end of the
!    routine, after scaling to avoid SO2 becoming negative.
!
!    If H2O2 is limiting, further oxidation by O3 occurs if there is
!    sufficient NH3 available to buffer the reaction; Choularton's
!    suggested procedure is:
!    a) Do oxidn of SO2 by H2O2   
!    b) Neutralise sulphate produced with NH3
!    c) If SO2 remains do oxidn by O3 until SO2 or NH3 required to 
!       neutralise it is exhausted.
!    
!--------------------------------------------------------------
!
      DO J=1,QPNTS-NPFLD+1,NPFLD     ! J loops over all model points
        DO I=FIRST_POINT,LAST_POINT  ! I loop omits N AND S polar rows
           K=I+J-1
!
          CLEARF(K)=1.0-CLOUDF(K)    ! CALCULATE CLEAR AIR FRACTION
!
        IF ( (QCL(K).GE.THOLD) .AND. (CLOUDF(K).GT.0.0) ) THEN
!
!   CALCULATE LAMDA
          TERM1=(CLEARF(K)*TAURATIO)**2
          TERM1=TERM1+2.0*TAURATIO*CLEARF(K)*(CLEARF(K)-CLOUDF(K))
          TERM1=SQRT(1.0+TERM1)
          TERM2=2.0*CLOUDF(K)-TAURATIO*CLEARF(K)-1.0
          TERM2=TERM2+TERM1
          LAMDA=TERM2/(2.0*CLOUDF(K))
!
!   CALCULATE PROB_CLEAR AND PROB_CLOUD
           DENOM=CLEARF(K)+CLOUDF(K)*LAMDA
         PROB_CLEAR=CLEARF(K)/DENOM
         PROB_CLOUD=CLOUDF(K)*LAMDA/DENOM
!
!   CALCULATE EXPECTED LIFETIME OF SO2 (WHICH IS 1/WETRATE)
           TERM1=PROBOX_FN1*PROB_CLEAR
           TERM2=PROBOX_FN2*PROB_CLOUD
           TERM2=(TERM1+TERM2)*CLEARF(K)/CLOUDF(K)
        DENOM=TERM2*CLOUDTAU+CHEMTAU
        WETRATE=1.0/DENOM
!
!!  RC'S H2O2 CHEMISTRY: OXIDATION OF SO2 TO SO4 IS CONTROLLED BY THE
!                        SMALLER OF THE SO2 AND H2O2 MIX RAT FIELDS
!
       IF (SO2(K).LE.(H2O2_MXR(K)*RELM_S_H2O2)) THEN
!
           DELTAS_WET(K) = (1.0-EXP(-WETRATE*TSTEP))*SO2(K)
!
       ELSE
!
       DELTAS_WET(K)=(1.0-EXP(-WETRATE*TSTEP))*H2O2_MXR(K)*RELM_S_H2O2
!
!   H2O2 field is limiting, so oxidise SO2 further with O3, using NH3
!   as buffer, if sufficient O3 and NH3 available. 
!
         IF ( L_SULPC_OZONE .AND. (O3(K).GE.O3_MIN) .AND.
     &           (NH3(K) .GT. DELTAS_WET(K)/RELM_S_2N) ) THEN
!
! O3 oxidation is controlled by smaller of SO2 and NH3 fields. 
! 2 atoms of N required for each S atom in ammonium sulphate.
! Sufficient NH3 must be available to neutralise all sulphate produced
! in grid box for O3 reaction to continue (otherwise PH too low)
! Note that all SO2 or NH3 is used in DELTAS_WET_O3 calcn; adjustment by
! SCALE_FACTOR at end of routine prevents removing too much.

           IF (SO2(K) .LE. (NH3(K)*RELM_S_2N))  THEN   ! SO2 limiting
! 
           DELTAS_WET_O3(K)=(1.0-EXP(-WETRATE*TSTEP))*SO2(K)
!
           ELSE                                        ! NH3 limiting
!
           DELTAS_WET_O3(K)=(1.0-EXP(-WETRATE*TSTEP))*NH3(K)*RELM_S_2N
!
           END IF
!
         END IF               ! End ozone oxidation block
!
       END IF                 ! End peroxide oxidation block
!
       END IF                 ! End cloud present condition
!
!
!
!
      END DO
      END DO
!
!------------------------------------------------------------------
! 3.  This section calculates amount of DMS dry oxidised to SO2
!     and MSA using a 3D OH concentration field to control the oxidn
!     rate when cos zenith angle is G.T. 10-6  (else rate is zero).
!     MSA is not saved, and K_DMS_OH is constant in this version.
!------------------------------------------------------------------
!
      IF (L_SULPC_DMS) THEN
!
      DO J=1,NPNTS-NPFLD+1,NPFLD      ! J loops over all model points
        DO I=FIRST_POINT,LAST_POINT   ! I loop omits N AND S polar rows
          K=I+J-1
!
!  Calculate DMS oxidation rate if daytime.
          IF (COSZA2D(I).GT.1.0E-06) THEN
            DMSRATE = K_DMS_OH * OH_CONC(K)
          ELSE                 ! ZERO RATE IF NIGHTIME
            DMSRATE=0.0
          ENDIF                   ! END COSZA2D IF
!
!  CALCULATE FRACTION OF DMS OXIDISED
!
          DELTA_DMS(K)=DMSRATE*TSTEP*DMS(K)    !*****DELTA_DMS*****
!
        END DO
      END DO
!
      END IF                   ! END L_SULPC_DMS IF

!---------------------------------------------------------------------
! 4. Release of aerosol from evaporating cloud droplets:
!    if no condensed water (liquid + ice) in grid box, release
!    dissolved sulphate as accumulation mode aerosol.
!     If cloud fraction less than 0.95, release some in clear  air.
!--------------------------------------------------------------------
!
      DO J=1,QPNTS-NPFLD+1,NPFLD      ! J loops over all model points
        DO I=FIRST_POINT,LAST_POINT   ! I loop omits N AND S polar rows
          K=I+J-1
          QCTOTAL(K) = QCL(K) + QCF(K)
          IF ( QCTOTAL(K) .LT. THOLD ) THEN
            DELTAS_EVAP(K) = SO4_DIS(K)
          ELSE IF  ( CLOUDF(K).LT.0.95 ) THEN  
            EVAPTIME = EVAPTAU + 0.5*CLOUDTAU   
            DELTAS_EVAP(K) = ( 1.0 - EXP(-TSTEP/EVAPTIME) )*SO4_DIS(K)
          ELSE
            DELTAS_EVAP(K) = 0.0
          ENDIF
        END DO
      END DO
!
!     Also release any dissolved aerosol in a non-wet level,
!     where it should not be.
!
      IF (QPNTS .LT. NPNTS) THEN
      DO J = QPNTS+1,NPNTS-NPFLD+1,NPFLD
        DO I=FIRST_POINT,LAST_POINT    ! I loop omits N AND S polar rows
          K=I+J-1
          DELTAS_EVAP(K) = SO4_DIS(K)
        END DO
      END DO
!
      ENDIF
!
!-------------------------------------------------------------------
! 5. Nucleation of accumulation mode aerosol forming dissolved SO4
!-------------------------------------------------------------------
!
!    THIS CODE ASSUMES THAT THE PARAMETER NUCTAU, WHICH IS THE
!    TIMESCALE FOR NUCLEATION ONCE A PARTICLE ENTERS A CLOUD, IS
!    VERY SHORT COMPARED WITH CLOUDTAU.
!
      DO J=1,QPNTS-NPFLD+1,NPFLD     ! J loops over all model points
        DO I=FIRST_POINT,LAST_POINT  ! I loop omits N AND S polar rows
           K=I+J-1
          IF ((QCTOTAL(K) .GE. THOLD) .AND. (CLOUDF(K).GT.0.0)) THEN
            NUCTIME=NUCTAU + ( (CLEARF(K)*CLOUDTAU)/(2.0*CLOUDF(K)) )
            DELTAS_NUCL(K) = ( 1.0 - EXP(-TSTEP/NUCTIME) )*SO4_ACC(K)
          ENDIF
        END DO
      END DO
!
!-------------------------------------------------------------------
! 6. Diffusional scavenging of Aitken mode SO4 forming dissolved SO4
!-------------------------------------------------------------------
!
!    THIS IS A MUCH SLOWER PROCESS THAN NUCLEATION AND THEREFORE WE
!    CANNOT MAKE THE SAME APPROXIMATIONS USED IN THAT CASE.
!
      DO J=1,QPNTS-NPFLD+1,NPFLD     ! J loops over all model points
        DO I=FIRST_POINT,LAST_POINT  ! I loop omits N AND S polar rows
           K=I+J-1
!
          IF ((QCTOTAL(K) .GE. THOLD) .AND. (CLOUDF(K).GT.0.0)) THEN
!
!    FIRST COMPUTE IN-CLOUD TIMESCALE FOR DIFFUSIONAL CAPTURE,
!    USING TOTAL CONDENSED WATER WITHIN THE CLOUDY PORTION OF THE BOX.
!    THE DIFFERENCE BETWEEN LIQUID WATER AND ICE IS NEGLECTED HERE.
!    THIS SHOULD BE IMPROVED ON EVENTUALLY.
!
            DIFFUSE_TAU =
     &      (RHO_CUBEROOT*(CLOUDF(K)/QCTOTAL(K))**0.333333)
     &               / ( 7.8*DIFFUSE_AIT )
!
!    USE DIFFERENT DROPLET DENSITIES OVER LAND AND SEA.
!    (AVOID THE TEMPTATION TO USE THE AEROSOL TO CONTROL THE DROPLET
!    CONCENTRATION, AS WE ARE ONLY SIMULATING ONE KIND OF AEROSOL.)
!
            IF( LAND_MASK(I) ) THEN
              DIFFUSE_TAU = DIFFUSE_TAU*NLAND23
            ELSE
              DIFFUSE_TAU = DIFFUSE_TAU*NSEA23
            ENDIF
!
            DIFFUSE_TAURATIO = CLOUDTAU/DIFFUSE_TAU
            PROBDIFF_INV = 1.0/( 1.0 - EXP(-DIFFUSE_TAURATIO) )
            PROBDIFF_FN1 = PROBDIFF_INV - 0.5
            PROBDIFF_FN2 = PROBDIFF_INV*EXP(-(0.5*DIFFUSE_TAURATIO) )
!
!     CALCULATE LAMBDA_AITKEN.
!
            TERM3 = (CLEARF(K)*DIFFUSE_TAURATIO)**2
            TERM3 = TERM3 + ( 2.0*DIFFUSE_TAURATIO
     &                       *CLEARF(K)*(CLEARF(K)-CLOUDF(K)) )
            TERM3 = SQRT(1.0+TERM3)
          TERM4=2.0*CLOUDF(K)-DIFFUSE_TAURATIO*CLEARF(K)-1.0
            TERM4 = TERM4 + TERM3
            LAMBDA_AITKEN = TERM4/(2.0*CLOUDF(K))
!
!   CALCULATE PROBDIFF_CLEAR AND PROBDIFF_CLOUD
!
            DENOM = CLEARF(K) + CLOUDF(K)*LAMBDA_AITKEN
            PROBDIFF_CLEAR = CLEARF(K)/DENOM
            PROBDIFF_CLOUD = (CLOUDF(K)*LAMBDA_AITKEN)/DENOM
!
!   CALCULATE EXPECTED LIFETIME OF AN AITKEN MODE PARTICLE WITH
!   RESPECT TO DIFFUSIVE CAPTURE BY CLOUD DROPLETS.
!
            TERM3 = PROBDIFF_FN1*PROBDIFF_CLEAR
            TERM4 = PROBDIFF_FN2*PROBDIFF_CLOUD
            TERM4 = (TERM3+TERM4)*CLEARF(K)/CLOUDF(K)
            DENOM = TERM4*CLOUDTAU + DIFFUSE_TAU
            DIFFRATE = 1.0/DENOM
!
!   NOW COMPUTE THE AMOUNT OF AITKEN MODE SO4 CONVERTED TO DISSOLVED
!   SO4 BY DIFFUSIVE CAPTURE DURING THE STEP.
!
            DELTAS_DIFFUSE(K) = (1.0-EXP(-DIFFRATE*TSTEP))*SO4_AIT(K)
!
          ENDIF
        END DO
      END DO
!
!-------------------------------------------------------------------
! 7. UPDATE SO2, SO4 modes and DMS, assuming:
!    Dry oxidation SO2 produces SO4 Aitken mode aerosol
!    Wet oxidation SO2 produces dissolved SO4
!    Dry oxidation DMS produces SO2 and MSA
!---------------------------------------------------------------------
!
! CHECK THAT AMOUNTS OF SO2 & DMS OXIDISED NOT GREATER THAN SO2 & DMS
! AVAILABLE AND SCALE INCREMENTS IF NECESSARY
!
      DO J=1,NPNTS-NPFLD+1,NPFLD     ! J loops over all model points
        DO I=FIRST_POINT,LAST_POINT  ! I loop omits N AND S polar rows
           K=I+J-1
          DELTAS_TOT(K) = DELTAS_DRY(K) + DELTAS_WET(K)
     &                              + DELTAS_WET_O3(K)  

!
          IF (DELTAS_TOT(K) .GT. SO2(K))  THEN
               SCALE_FACTOR = SO2(K)/DELTAS_TOT(K)
              DELTAS_DRY(K) = DELTAS_DRY(K) * SCALE_FACTOR
              DELTAS_WET(K) = DELTAS_WET(K) * SCALE_FACTOR
!
            IF (L_SULPC_OZONE) THEN    
              DELTAS_WET_O3(K) = DELTAS_WET_O3(K) * SCALE_FACTOR
            END IF
!
            DELTAS_TOT(K) = SO2(K)
          END IF
!
          IF (DELTA_DMS(K) .GT. DMS(K)) THEN
            DELTA_DMS(K) = DMS(K)
          END IF
!
!  UPDATE SO2, SO4 MODES AND DMS
!
          SO2(K) = SO2(K) + DELTA_DMS(K)*BRAT_SO2 - DELTAS_TOT(K)
!
      SO4_AIT(K)=SO4_AIT(K) + PSI(K)*DELTAS_DRY(K) - DELTAS_DIFFUSE(K)
!
          SO4_DIS(K) = SO4_DIS(K) + DELTAS_WET(K) - DELTAS_EVAP(K)
     &             + DELTAS_NUCL(K) + DELTAS_DIFFUSE(K)
     &                   + DELTAS_WET_O3(K)         
!
          DMS(K) = DMS(K) - DELTA_DMS(K)
!
          MSA(K) = DELTA_DMS(K)*BRAT_MSA
!
          SO4_ACC(K) = SO4_ACC(K) + DELTAS_EVAP(K) - DELTAS_NUCL(K)
     &             + (1.0-PSI(K))*DELTAS_DRY(K)
!
        END DO             ! END OF I LOOP
      END DO             ! END OF J LOOP
!
!
!--------------------------------------------------------------------
! 8.  DEPLETE  NH3 and H2O2, REPLENISH H2O2 from HO2
!--------------------------------------------------------------------
!
      DO J=1,QPNTS-NPFLD+1,NPFLD     ! J loops over wet levels    
      DO I=FIRST_POINT,LAST_POINT  ! I loop omits N AND S polar rows  
        K=I+J-1                                                      
!                                                                       
!!  Deplete H2O2                                                  
!
        IF (DELTAS_WET(K) .GT. 0.0)  THEN
          H2O2_MXR(K)=H2O2_MXR(K) - DELTAS_WET(K)/RELM_S_H2O2
        END IF
!                                                                       
!   Replenish H2O2_MXR field from HO2 field in dry part of grid box
!   (whether or not oxidn has occurred) in daylight only.
!
        IF ( (H2O2_MXR(K) .LT. H2O2_LMT(K)) .AND.
     &          (COSZA2D(I) .GT. 1.0E-06) )  THEN    
!                                                                       
!   Calculate replenishment in concn, then convert to mix ratio
!   Reaction rate K_HO2_HO2 is a fn of T, P, AND Q                  
!                                                                       
          AIR_CONC = (AVOGADRO*P(K)) / (R*T(K)*RMM_AIR*1.0E6) !MCLS/CC 
                                                                        
            W_CONC = Q(K) * AIR_CONC * RMM_AIR / RMM_W        !MCL/CC  
!                                                                       
          H2O2_CON_TO_MXR = RMM_H2O2 / (RMM_AIR * AIR_CONC)   !CC/MCL  
!                                                                       
          K_HO2_HO2 = ( K2*EXP(T2/T(K)) + AIR_CONC*K3*EXP(T3/T(K)) ) *
     &                (1.0 + W_CONC*K4*EXP(T4/T(K)) )                   
!
!
          H2O2_REP = HO2_CONC(K)*HO2_CONC(K)*K_HO2_HO2*TSTEP*CLEARF(K)
!                                                                       
          H2O2_REP = H2O2_CON_TO_MXR * H2O2_REP
!                                                                       
!  Increment H2O2_MXR up to H2O2_LMT value                     
!                                                                       
          H2O2_MXR(K) = H2O2_MXR(K) + H2O2_REP
!
          IF ( H2O2_MXR(K) .GT. H2O2_LMT(K) ) THEN
            H2O2_MXR(K)=H2O2_LMT(K)                                     
          END IF
!                                                                       
      END IF                     !End replenishment condition
!                                                                       
!! Deplete NH3 if ozone oxidation included and wet oxidn has occurred
!
      IF (L_SULPC_OZONE .AND. ((DELTAS_WET(K) .GT. 0.0) .OR.
     &                      (DELTAS_WET_O3(K) .GT. 0.0)) )   THEN
!
        NH3_DEP(K) = (DELTAS_WET(K)+DELTAS_WET_O3(K))/RELM_S_2N
!
        IF ( NH3_DEP(K) .GT. NH3(K) ) THEN            
          NH3_DEP(K) = NH3(K)
        END IF
!
        NH3(K) = NH3(K) - NH3_DEP(K)
!
      END IF       ! end L_SULPC_OZONE
!
      END DO                                                            
      END DO                                                            
!                                                                      
!-----------------------------------------------------------------------
!
!
      RETURN
      END
!
