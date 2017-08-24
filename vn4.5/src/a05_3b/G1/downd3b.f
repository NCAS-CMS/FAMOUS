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
CLL  SUBROUTINE DOWND--------------------------------------------------
CLL
CLL  PURPOSE : CALL DOWNDRAUGHT CALCULATION
CLL
CLL            CHANGE OF PHASE CALCULATION WHERE NO DOWNDRAUGHT OCCURS
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   4.3  03/02/97   As DOWND3A (vn4.2) except TIMESTEP passed to
CLL                   CHG_PHSE
CLL                                                   R.N.B.Smith
CLL   4.5    20/7/98 : Type of ICCB changed to INTEGER 
CLL                    Kill the IBM specific lines (JCThil)
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DOWND (NPNTS,NP_FULL,K,KCT,THDD_K,QDD_K,THE_K,THE_KM1,
     &                  QE_K,QE_KM1,DTHBYDT_K,DTHBYDT_KM1,DQBYDT_K,
     &                  DQBYDT_KM1,FLX_DD_K,P_KM1,DELPK,DELPKM1,EXK,
     &                  EXKM1,DELTD,DELQD,AMDETK,EKM14,EKM34,PRECIP_K,
     &                  RAIN,SNOW,ICCB,BWATER_K,BDD_START,
     &                  BDDWT_K,BDDWT_KM1,BDD_ON,RAIN_ENV,SNOW_ENV,
     &                  RAIN_DD,SNOW_DD,FLX_UD_K,TIMESTEP,CCA,NDDON_A,
     &                  LR_UD_REF,L_MOM,UDD_K,VDD_K,UE_K,VE_K,UE_KM1,
     &                  VE_KM1,DUBYDT_K,DVBYDT_K,DUBYDT_KM1,DVBYDT_KM1,
     &                  DELUD,DELVD,EFLUX_U_DD,EFLUX_V_DD,
     &                  L_TRACER,NTRA,TRADD_K,TRAE_K,
     &                  TRAE_KM1,DTRABYDT_K,DTRABYDT_KM1,DELTRAD)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
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
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
C
      INTEGER I,KTRA             ! LOOP COUNTERS
C
      INTEGER K                  ! IN PRESENT MODEL LAYER
C
      INTEGER NPNTS              ! IN NUMBER OF POINTS
C
      INTEGER NDDON,NDDON_A      ! NUMBER OF POINTS AT WHICH
                                 ! DOWNDRAUGHT DOES OCCUR
C
      INTEGER NP_FULL            ! IN FULL VECTOR LENGTH
C
      INTEGER NTRA               ! NUMBER OF TRACERS
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      INTEGER KCT                ! IN CONVECTIVE CLOUD TOP LAYER
C
      REAL THDD_K(NPNTS)         ! IN MODEL POTENTIAL TEMPERATURE
                                 !    OF DOWNDRAUGHT IN LAYER K (K)
C
      REAL QDD_K(NPNTS)          ! IN MIXING RATIO OF DOWNDRAUGHT IN
                                 !    LAYER K (KG/KG)
C
      REAL UDD_K(NPNTS)          ! IN U IN DOWNDRAUGHT IN LAYER K (M/S)
C
      REAL VDD_K(NPNTS)          ! IN V IN DOWNDRAUGHT IN LAYER K (M/S)
C
      REAL TRADD_K(NP_FULL,NTRA) ! IN TRACER CONTENT OF DOWNDRAUGHT
                                 !    IN LAYER K (KG/KG)
C
      REAL THE_K(NPNTS)          ! IN POTENTIAL TEMPERATURE OF
                                 !    ENVIRONMENT IN LAYER K (K)
C
      REAL THE_KM1(NPNTS)        ! IN POTENTIAL TEMPERATURE OF
                                 !    ENVIRONMENT IN LAYER K-1 (K)
C
      REAL QE_K(NPNTS)           ! IN MIXING RATIO OF ENVIRONMENT IN
                                 !    LAYER K (KG/KG)
C
      REAL QE_KM1(NPNTS)         ! IN MIXING RATIO OF ENVIRONMENT IN
                                 !    LAYER K-1 (KG/KG)
C
      REAL UE_K(NPNTS)           ! IN U OF ENVIRONMENT IN LAYER K (M/S)
C
      REAL UE_KM1(NPNTS)         ! IN U OF ENVIRONMENT IN LAYER K-1
                                 !    (M/S)
C
      REAL VE_K(NPNTS)           ! IN V OF ENVIRONMENT IN LAYER K (M/S)
C
      REAL VE_KM1(NPNTS)         ! IN V OF ENVIRONMENT IN LAYER K-1
                                 !    (M/S)
C
      REAL TRAE_K(NP_FULL,NTRA)  ! IN TRACER CONTENT OF ENVIRONMENT
                                 !    IN LAYER K (KG/KG)
C
      REAL TRAE_KM1(NP_FULL,NTRA)! IN TRACER CONTENT OF ENVIRONMENT
                                 !    IN LAYER K-1 (KG/KG)
C
      REAL FLX_DD_K(NPNTS)       ! IN DOWNDRAUGHT MASS FLUX OF LAYER K
                                 !    (PA/S)
C
      REAL P_KM1(NPNTS)          ! IN PRESSURE OF LAYER K-1 (PA)
C
      REAL DELPK(NPNTS)          ! IN PRESSURE DIFFERENCE ACROSS
                                 !    LAYER K (PA)
C
      REAL DELPKM1(NPNTS)        ! IN PRESSURE DIFFERENCE ACROSS
                                 !    LAYER K-1 (PA)
C
      REAL EXK(NPNTS)            ! IN EXNER RATIO FOR LAYER K
C
      REAL EXKM1(NPNTS)          ! IN EXNER RATIO FOR LAYER K-1
C
      REAL PRECIP_K(NPNTS)       ! IN PRECIPITATION ADDED WHEN
                                 !    DESCENDING FROM LAYER K TO K-1
                                 !    (KG/M**2/S)
C
      REAL AMDETK(NPNTS)         ! IN MIXING DETRAINMENT AT LEVEL K
                                 !    MULTIPLIED BY APPROPRIATE LAYER
                                 !    THICKNESS
C
      REAL EKM14(NPNTS)          ! IN EXNER RATIO AT LAYER K-1/4
C
      REAL EKM34(NPNTS)          ! IN EXNER RATIO AT LAYER K-3/4
C
      REAL DELTD(NPNTS)          ! IN COOLING NECESSARY TO
                                 !    ACHIEVE SATURATION (K)
C
      REAL DELQD(NPNTS)          ! IN MOISTENING NECESSARY TO
                                 !    ACHIEVE SATURATION (KG/KG)
C
      REAL DELUD(NPNTS)          ! IN CHANGE IN ENVIRONMENT U DUE TO
                                 !    DOWNDRAUGHT FORMATION (M/S)
C
      REAL DELVD(NPNTS)          ! IN CHANGE IN ENVIRONMENT V DUE TO
                                 !    DOWNDRAUGHT FORMATION (M/S)
C
      REAL DELTRAD(NP_FULL,NTRA) ! IN DEPLETION OF ENVIRONMENT TRACER
                                 !    DUE TO DOWNDRAUGHT FORMATION
C
      INTEGER ICCB(NPNTS)           ! IN CLOUD BASE LEVEL               
C
      LOGICAL BWATER_K(NPNTS)    ! IN MASK FOR THOSE POINTS AT WHICH
                                 !    CONDENSATE IS WATER IN LAYER K
C
      LOGICAL BDDWT_K(NPNTS)     ! IN MASK FOR THOSE POINTS IN
                                 !    DOWNDRAUGHT WHERE PRECIPITATION
                                 !    IS LIQUID IN LAYER K
C
      LOGICAL BDDWT_KM1(NPNTS)   ! IN MASK FOR THOSE POINTS IN
                                 !    DOWNDRAUGHT WHERE PRECIPITATION
                                 !    IS LIQUID IN LAYER K-1
C
      LOGICAL L_TRACER           ! IN SWITCH FOR INCLUSION OF TRACERS
C
      LOGICAL L_MOM              ! IN SWITCH FOR INCLUSION OF
                                 !    MOMENTUM TRANSPORTS
C
      REAL RAIN_ENV(NPNTS)       ! IN AMOUNT OF RAIN FALLING THROUGH
                                 !    THE ENVIRONMENT
C
      REAL SNOW_ENV(NPNTS)       ! IN AMOUNT OF SNOW FALLING THROUGH
                                 !    THE ENVIRONMENT
C
      REAL RAIN_DD(NPNTS)        ! IN AMOUNT OF RAIN FALLING THROUGH
                                 !    THE DOWNDRAUGHT
C
      REAL SNOW_DD(NPNTS)        ! IN AMOUNT OF SNOW FALLING THROUGH
                                 !    THE DOWNDRAUGHT
C
      REAL FLX_UD_K(NPNTS)       ! IN UPDRAUGHT MASSFLUX AT LAYER K
C
      REAL TIMESTEP              ! IN MODEL TIMESTEP (S)
C
      REAL CCA(NPNTS)            ! IN CONVECTIVE CLOUD AMOUNT
C
      REAL LR_UD_REF(NPNTS)      ! IN UD PPN MIXING RATION IN LOWEST
                                 !    PRECIPITATING LAYER IN UD
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      LOGICAL BDD_START(NPNTS)   ! INOUT
                                 ! IN  MASK FOR THOSE POINTS WHERE
                                 !     DOWNDRAUGHT MAY FORM IN LAYER K
                                 ! OUT MASK FOR THOSE POINTS WHERE
                                 !     DOWNDRAUGHT MAY FORM IN LAYER
                                 !     K-1
C
      REAL DTHBYDT_K(NPNTS)      ! INOUT
                                 ! IN  INCREMENT TO MODEL POTENTIAL
                                 !     TEMPERATURE OF LAYER K (K/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     POTENTIAL TEMPERATURE OF LAYER K
                                 !     (K/S)
C
      REAL DTHBYDT_KM1(NPNTS)    ! INOUT
                                 ! IN  INCREMENT TO MODEL POTENTIAL
                                 !     TEMPERATURE OF LAYER K-1 (K/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     POTENTIAL TEMPERATURE OF
                                 !     LAYER K-1 (K/S)
C
      REAL DQBYDT_K(NPNTS)       ! INOUT
                                 ! IN  INCREMENT TO MODEL MIXING
                                 !     RATIO OF LAYER K (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     MIXING RATIO OF LAYER K (KG/KG/S)
C
      REAL DQBYDT_KM1(NPNTS)     ! INOUT
                                 ! IN  INCREMENT TO MODEL MIXING
                                 !     RATIO OF LAYER K-1 (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     MIXING RATIO OF
                                 !     LAYER K-1 (KG/KG/S)
C
      REAL DUBYDT_K(NPNTS)       ! INOUT
                                 ! IN  INCREMENT TO MODEL U IN
                                 !     LAYER K (M/S**2)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     U IN LAYER K (M/S**2)
C
      REAL DUBYDT_KM1(NPNTS)     ! INOUT
                                 ! IN  INCREMENT TO MODEL U
                                 !     IN LAYER K-1 (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     U IN LAYER K-1 (M/S**2)
C
      REAL DVBYDT_K(NPNTS)       ! INOUT
                                 ! IN  INCREMENT TO MODEL V IN
                                 !     LAYER K (M/S**2)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     V IN LAYER K (M/S**2)
C
      REAL DVBYDT_KM1(NPNTS)     ! INOUT
                                 ! IN  INCREMENT TO MODEL V
                                 !     IN LAYER K-1 (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     V IN LAYER K-1 (M/S**2)
C
      REAL DTRABYDT_K(NP_FULL,   ! INOUT
     *                NTRA)      ! IN  INCREMENT TO MODEL TRACER OF
                                 !     LAYER K (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MODEL

C
      REAL DTRABYDT_KM1(NP_FULL, ! INOUT
     *                 NTRA)     ! IN  INCREMENT TO MODEL TRACER OF
                                 !     LAYER K-1 (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     TRACER IN LAYER K-1
                                 !     (KG/KG/S)
C
      REAL RAIN (NPNTS)          ! INOUT
                                 ! IN  INITIALISED RAINFALL (KG/M**2/S)
                                 ! OUT SURFACE RAINFALL (KG/M**2/S)
C
      REAL SNOW(NPNTS)           ! INOUT
                                 ! IN  INITIALISED SNOWFALL (KG/M**2/S)
                                 ! OUT SURFACE SNOWFALL (KG/M**2/S)
C
      LOGICAL BDD_ON(NPNTS)      ! INOUT
                                 ! IN  MASK FOR THOSE POINTS WHERE DD
                                 !     HAS CONTINUED FROM PREVIOUS LAYER
                                 ! OUT MASK FOR THOSE POINTS WHERE DD
                                 !     CONTINUES TO LAYER K-1
C
      REAL EFLUX_U_DD(NPNTS),    ! INOUT
     *     EFLUX_V_DD(NPNTS)     ! IN  EDDY FLUX OF MOMENTUM DUE TO DD
                                 !     AT TOP OF A LAYER
                                 ! OUT EDDY FLUX OF MOMENTUM DUE TO DD
                                 !     AT BOTTOM OF A LAYER
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
C
      REAL WORK(NDDON_A,38)      !  WORK SPACE
C
      LOGICAL BWORK(NDDON_A,5)   !  WORK SPACE FOR 'BIT' MASKS
C
      INTEGER INDEX1(NDDON_A)    !  INDEX FOR COMPRESS AND
C
      LOGICAL B_DD_END(NPNTS)    !  MASK FOR POINTS WHERE DOWNDRAUGHT
                                 ! HAS ENDED
C
      REAL TRADD_K_C(NDDON_A,    ! TRACER CONTENT IN DOWNDRAUGHT AT
     *               NTRA)       ! LAYER K - COMPRESSED (KG/KG)
C
      REAL TRAE_K_C(NDDON_A,NTRA)! TRACER CONTENT OF ENVIRONMENT AT
                                 ! LAYER K - COMPRESSED (KG/KG)
C
      REAL TRAE_KM1_C(NDDON_A,   ! TRACER CONTENT OF ENVIRONMENT IN
     *                NTRA)      ! LAYER K-1 - COMPRESSED (KG/KG)
C
      REAL DTRA_K_C(NDDON_A,NTRA)! INCREMENT TO MODEL TRACER IN LAYER
                                 ! K - COMPRESSED (KG/KG/S)
C
      REAL DTRA_KM1_C(NDDON_A,   ! INCREMENT TO MODEL TRACER IN LAYER
     *                NTRA)      ! K-1 - COMPRESSED (KG/KG/S)
C
      REAL DELTRAD_C(NDDON_A,    ! DEPLETION OF ENVIRONMENT TRACER
     *               NTRA)       ! DUE TO DOWNDRAUGHT FORMATION
                                 ! COMPRESSED
C
C
      REAL FACTOR                !  PROPORTION OF RAINFALL GOING INTO
                                 !  DOWNDRAUGHT FROM UD
C
      REAL FACTOR_ENV            !  PROPORTION OF RAINFALL GOING INTO
                                 !  DD FROM FALLING PPN
C
      REAL PPN_DD_REF            !  REFERENCE DD PPN MASS
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C-----------------------------------------------------------------------
C
      EXTERNAL CHG_PHSE, PEVP_BCB, DDRAUGHT
C
C-----------------------------------------------------------------------
C START OF MAIN LOOP
C   UPDATE PRECIPITATION AND CALCULATE MASK FOR WHERE PRECIPITATION
C   IS LIQUID
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
        B_DD_END(I) = .FALSE.
      END DO
C
      IF (K.EQ.KCT+1) THEN
        DO I=1,NPNTS
         RAIN_DD(I) = 0.0
         RAIN_ENV(I) = 0.0
         SNOW_DD(I) = 0.0
         SNOW_ENV(I) = 0.0
        END DO
      END IF
C
C----------------------------------------------------------------------
C INJECTION OF PRECIPITATION FROM UD AT LEVEL K
C----------------------------------------------------------------------
C
      DO I=1,NPNTS
       FACTOR= 0.0
       IF (BDD_ON(I) .AND. FLX_UD_K(I).GT.0.0) THEN
        FACTOR = G * FLX_DD_K(I)/FLX_UD_K(I)
        FACTOR = AMIN1(FACTOR,1.0)
       END IF
c
       IF (BWATER_K(I)) THEN
        RAIN_DD(I) = RAIN_DD(I) + PRECIP_K(I)*FACTOR
        RAIN_ENV(I) = RAIN_ENV(I) + PRECIP_K(I)*(1.0-FACTOR)
       ELSE
        SNOW_DD(I) = SNOW_DD(I) + PRECIP_K(I)*FACTOR
        SNOW_ENV(I) = SNOW_ENV(I) + PRECIP_K(I)*(1.0-FACTOR)
       END IF
c
      END DO
C
C----------------------------------------------------------------------
C INTERACTION OF DOWNDRAUGHT WITH RESERVE OF PRECIPITATION OUTSIDE
C DOWNDRAUGHT
C
C BASED UPON CONTINUITY OF PRECIPITATION MIXING RATIO WITHIN
C DOWNDRAUGHT - EITHER AFTER INJECTION OF RAIN FROM UD IN LEVEL
C K OR WITH PPN MIXING RATIO IN LOWEST PRECIPITATING LAYER
C
C IF DOWNDRAUGHT INCREASES IN MASS THEN WATER INJECTED
C IF DOWNDRAUGHT DECREASES IN MASS THEN WATER IS REMOVED
C
C----------------------------------------------------------------------
C
      DO I=1,NPNTS
C
       IF (BDD_ON(I)) THEN
C
        FACTOR_ENV = 0.0
        IF (PRECIP_K(I).GT.0.0) THEN
C
C---------------------------------------------------------------------
C CALCULATE NEW REFERENCE PPN MIXING RATIO
C DD PPN MIXING RATIO IN LAYER KM1 BASED ON CONTINUITY
C WITH THAT IN LAYER K
C---------------------------------------------------------------------
C
         LR_UD_REF(I) = G * PRECIP_K(I)/FLX_UD_K(I)
         PPN_DD_REF = RAIN_DD(I)+SNOW_DD(I)
        ELSE
C
C---------------------------------------------------------------------
C DD PPN MIXING RATIO IN LAYER KM1 BASED ON CONTINUITY
C WITH THAT IN LAST PRECIPITATING UD LAYER
C---------------------------------------------------------------------
C
         PPN_DD_REF = LR_UD_REF(I) * FLX_DD_K(I)
        END IF
C
C--------------------------------------------------------------------
C INJECT PPN INTO DD FROM PPN FALLING OUTSIDE OF THE DD
C--------------------------------------------------------------------
C
        IF ((RAIN_ENV(I) + SNOW_ENV(I)) .GT. 0.0) THEN
!-------Already inside IF ( BDD_ON(I)) block----------------------------
         FACTOR_ENV = ( (PPN_DD_REF * (1.0+EKM14(I))*
     *                    (1.0+EKM34(I))*(1.0-AMDETK(I))) -
     *                         (RAIN_DD(I)+SNOW_DD(I)) ) /
     *                          (RAIN_ENV(I)+SNOW_ENV(I))
         FACTOR_ENV = AMIN1(FACTOR_ENV,1.0)
         FACTOR_ENV = AMAX1(FACTOR_ENV,-1.0)
        END IF
C
        IF (FACTOR_ENV.GT.0.0) THEN
         RAIN_DD(I) = RAIN_DD(I) + RAIN_ENV(I)*FACTOR_ENV
         RAIN_ENV(I) = RAIN_ENV(I) * (1.0-FACTOR_ENV)
         SNOW_DD(I) = SNOW_DD(I) + SNOW_ENV(I)*FACTOR_ENV
         SNOW_ENV(I) = SNOW_ENV(I) * (1.0-FACTOR_ENV)
        ELSE
         RAIN_ENV(I) = RAIN_ENV(I) - RAIN_DD(I)*FACTOR_ENV
         RAIN_DD(I) = RAIN_DD(I) * (1.0+FACTOR_ENV)
         SNOW_ENV(I) = SNOW_ENV(I) - SNOW_DD(I)*FACTOR_ENV
         SNOW_DD(I) = SNOW_DD(I) * (1.0+FACTOR_ENV)
        END IF
C
       END IF
C
C--------------------------------------------------------------------
C ZERO PRECIPITATION RATE IN LAYER K
C--------------------------------------------------------------------
C
       PRECIP_K(I) = 0.0
C
      END DO
C
C
C-----------------------------------------------------------------------
C COMPRESS OUT ON BASIS OF BIT VECTOR BDDON - THOSE POINTS WITH A
C DOWNDRAUGHT
C-----------------------------------------------------------------------
C
      NDDON=0
C
      DO I=1,NPNTS
        IF (BDD_ON(I)) THEN
           NDDON = NDDON+1
           INDEX1(NDDON) = I
        END IF
      END DO
C
      IF (NDDON .NE. 0) THEN
         DO I=1,NDDON
          WORK(I,1) = THDD_K(INDEX1(I))
          WORK(I,2) = QDD_K(INDEX1(I))
          WORK(I,3) = THE_K(INDEX1(I))
          WORK(I,4) = THE_KM1(INDEX1(I))
          WORK(I,5) = QE_K(INDEX1(I))
          WORK(I,6) = QE_KM1(INDEX1(I))
          WORK(I,7) = DTHBYDT_K(INDEX1(I))
          WORK(I,8) = DTHBYDT_KM1(INDEX1(I))
          WORK(I,9) = DQBYDT_K(INDEX1(I))
          WORK(I,10) = DQBYDT_KM1(INDEX1(I))
          WORK(I,11) = FLX_DD_K(INDEX1(I))
          WORK(I,12) = P_KM1(INDEX1(I))
          WORK(I,13) = DELPK(INDEX1(I))
          WORK(I,14) = DELPKM1(INDEX1(I))
          WORK(I,15) = EXK(INDEX1(I))
          WORK(I,16) = EXKM1(INDEX1(I))
          WORK(I,17) = DELTD(INDEX1(I))
          WORK(I,18) = DELQD(INDEX1(I))
          WORK(I,19) = AMDETK(INDEX1(I))
          WORK(I,20) = EKM14(INDEX1(I))
          WORK(I,21) = EKM34(INDEX1(I))
          WORK(I,22) = RAIN_DD(INDEX1(I))
          WORK(I,23) = SNOW_DD(INDEX1(I))
          WORK(I,24) = CCA(INDEX1(I))
          BWORK(I,1) = BDD_START(INDEX1(I))
          BWORK(I,2) = BDDWT_K(INDEX1(I))
          BWORK(I,3) = BDDWT_KM1(INDEX1(I))
          BWORK(I,4) = BDD_ON(INDEX1(I))
          BWORK(I,5) = B_DD_END(INDEX1(I))
      END DO
C
      IF(L_MOM)THEN
        DO I=1,NDDON
          WORK(I,25) = UDD_K(INDEX1(I))
          WORK(I,26) = VDD_K(INDEX1(I))
          WORK(I,27) = UE_K(INDEX1(I))
          WORK(I,28) = VE_K(INDEX1(I))
          WORK(I,29) = UE_KM1(INDEX1(I))
          WORK(I,30) = VE_KM1(INDEX1(I))
          WORK(I,31) = DUBYDT_K(INDEX1(I))
          WORK(I,32) = DUBYDT_KM1(INDEX1(I))
          WORK(I,33) = DVBYDT_K(INDEX1(I))
          WORK(I,34) = DVBYDT_KM1(INDEX1(I))
          WORK(I,35) = DELUD(INDEX1(I))
          WORK(I,36) = DELVD(INDEX1(I))
          WORK(I,37) = EFLUX_U_DD(INDEX1(I))
          WORK(I,38) = EFLUX_V_DD(INDEX1(I))
        END DO
      END IF
C
      IF(L_TRACER)THEN
C
      DO KTRA=1,NTRA
        DO I=1,NDDON
          TRADD_K_C(I,KTRA) = TRADD_K(INDEX1(I),KTRA)
          TRAE_K_C(I,KTRA) = TRAE_K(INDEX1(I),KTRA)
          TRAE_KM1_C(I,KTRA) = TRAE_KM1(INDEX1(I),KTRA)
          DTRA_K_C(I,KTRA)  = DTRABYDT_K(INDEX1(I),KTRA)
          DTRA_KM1_C(I,KTRA) = DTRABYDT_KM1(INDEX1(I),KTRA)
          DELTRAD_C(I,KTRA) = DELTRAD(INDEX1(I),KTRA)
        END DO
      END DO
C
      END IF
C
C
C-----------------------------------------------------------------------
C START DOWNDRAUGHT CALCULATION
C-----------------------------------------------------------------------
C
C
         CALL DDRAUGHT (NDDON,NDDON_A,K,KCT,WORK(1,1),WORK(1,2),
     &                  WORK(1,3),WORK(1,4),WORK(1,5),WORK(1,6),
     &                  WORK(1,7),WORK(1,8),WORK(1,9),WORK(1,10),
     &                  WORK(1,11),WORK(1,12),WORK(1,13),WORK(1,14),
     &                  WORK(1,15),WORK(1,16),WORK(1,17),WORK(1,18),
     &                  WORK(1,19),WORK(1,20),WORK(1,21),WORK(1,22),
     &                  WORK(1,23),BWORK(1,1),BWORK(1,2),BWORK(1,3),
     &                  BWORK(1,4),BWORK(1,5),WORK(1,24),L_MOM,
     &                  WORK(1,25),WORK(1,26),WORK(1,27),WORK(1,28),
     &                  WORK(1,29),WORK(1,30),WORK(1,31),WORK(1,32),
     &                  WORK(1,33),WORK(1,34),WORK(1,35),WORK(1,36),
     &                  WORK(1,37),WORK(1,38),
     &                  L_TRACER,NTRA,TRADD_K_C,TRAE_K_C,TRAE_KM1_C,
     &                  DTRA_K_C,DTRA_KM1_C,DELTRAD_C)
C
C-----------------------------------------------------------------------
C EXPAND REQUIRED VECTORS BACK TO FULL FIELDS
C-----------------------------------------------------------------------
C
      DO I=1,NDDON
       THDD_K(INDEX1(I)) = WORK(I,1)
       QDD_K(INDEX1(I)) = WORK(I,2)
       DTHBYDT_K(INDEX1(I)) = WORK(I,7)
       DTHBYDT_KM1(INDEX1(I)) = WORK(I,8)
       DQBYDT_K(INDEX1(I)) = WORK(I,9)
       DQBYDT_KM1(INDEX1(I)) = WORK(I,10)
       FLX_DD_K(INDEX1(I)) = WORK(I,11)
       RAIN_DD(INDEX1(I)) = WORK(I,22)
       SNOW_DD(INDEX1(I)) = WORK(I,23)
       BDD_START(INDEX1(I)) = BWORK(I,1)
       BDDWT_K(INDEX1(I)) = BWORK(I,2)
       BDDWT_KM1(INDEX1(I)) = BWORK(I,3)
       BDD_ON(INDEX1(I)) = BWORK(I,4)
       B_DD_END(INDEX1(I)) = BWORK(I,5)
      END DO
C
      IF(L_MOM)THEN
       DO I=1,NDDON
        UDD_K(INDEX1(I)) = WORK(I,25)
        VDD_K(INDEX1(I)) = WORK(I,26)
        DUBYDT_K(INDEX1(I)) = WORK(I,31)
        DUBYDT_KM1(INDEX1(I)) = WORK(I,32)
        DVBYDT_K(INDEX1(I)) = WORK(I,33)
        DVBYDT_KM1(INDEX1(I)) = WORK(I,34)
        EFLUX_U_DD(INDEX1(I)) = WORK(I,37)
        EFLUX_V_DD(INDEX1(I)) = WORK(I,38)
       END DO
      END IF
C
      IF(L_TRACER)THEN
C
      DO KTRA=1,NTRA
        DO I=1,NDDON
          TRADD_K(INDEX1(I),KTRA) = TRADD_K_C(I,KTRA)
          DTRABYDT_K(INDEX1(I),KTRA) = DTRA_K_C(I,KTRA)
          DTRABYDT_KM1(INDEX1(I),KTRA) = DTRA_KM1_C(I,KTRA)
        END DO
      END DO
C
      END IF
C
      END IF
C
C-----------------------------------------------------------------------
C RESET PRECIPITATION FALLING THROUGH ENVIRONMENT IF DOWNDRAUGHT
C DID NOT FORM
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
        IF (.NOT.BDD_ON(I).AND..NOT.B_DD_END(I)) THEN
          RAIN_ENV(I) = RAIN_ENV(I)+RAIN_DD(I)
          SNOW_ENV(I) = SNOW_ENV(I)+SNOW_DD(I)
          RAIN_DD(I) = 0.0
          SNOW_DD(I) = 0.0
        END IF
      END DO
C
C-----------------------------------------------------------------------
C CARRY OUT CHANGE OF PHASE CALCULATION FOR PRECIPITATION FALLING
C THROUGH ENVIRONMENT
C-----------------------------------------------------------------------
C
         CALL CHG_PHSE (NPNTS,K,RAIN_ENV,SNOW_ENV,DTHBYDT_KM1,
     &                  EXK,EXKM1,DELPKM1,THE_K,THE_KM1,
     &                  TIMESTEP,CCA)
C
C-----------------------------------------------------------------------
C EVAPORATE RAIN FALLING THROUGH ENVIRONMENT IF LAYER K BELOW
C CLOUD BASE
C-----------------------------------------------------------------------
C
         CALL PEVP_BCB (NPNTS,K-1,ICCB,THE_KM1,P_KM1,QE_KM1,DELPKM1,
     &                  RAIN_ENV,SNOW_ENV,DTHBYDT_KM1,DQBYDT_KM1,
     &                  EXKM1,TIMESTEP,CCA)
C
C-----------------------------------------------------------------------
C RESET PRECIPITATION FALLING THROUGH ENVIRONMENT IF DOWNDRAUGHT
C TERMINATES
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
        IF (B_DD_END(I)) THEN
          RAIN_ENV(I) = RAIN_ENV(I)+RAIN_DD(I)
          SNOW_ENV(I) = SNOW_ENV(I)+SNOW_DD(I)
          RAIN_DD(I) = 0.0
          SNOW_DD(I) = 0.0
        END IF
      END DO
C
C-----------------------------------------------------------------------
C UPDATE RAIN AND SNOW
C-----------------------------------------------------------------------
C
       IF (K.EQ.2) THEN
         DO I=1,NPNTS
           RAIN(I) = RAIN(I)+RAIN_DD(I)+RAIN_ENV(I)
           SNOW(I) = SNOW(I)+SNOW_DD(I)+SNOW_ENV(I)
         END DO
       END IF
C
      RETURN
      END
C
