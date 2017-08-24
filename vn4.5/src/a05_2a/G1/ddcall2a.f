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
CLL  SUBROUTINE DD_CALL------------------------------------------------
CLL
CLL  PURPOSE : CALCULATE INITIAL DOWNDRAUGHT MASSFLUX
CLL
CLL            RESET EN/DETRAINMENT RATES FOR DOWNDRAUGHT
CLL
CLL            COMPRESS/EXPAND VARIABLES
CLL
CLL            INITIALISE DOWNDRAUGHT
CLL
CLL            CALL DOWNDRAUGHT ROUTINE
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.3   23/12/93 : DG060893 : CORRECTION TO REDUCE OVER PREDICTION
CLL                               OF CONVECTIVE SNOW; TO PASS ADDITIONAL
CLL                               DATA DOWN TO DOWN2A AND PREVENT DD
CLL                               FORMING BELOW UPDRAUGHT BASE
CLL   4.2    Oct. 96  T3E migration: *DEF CRAY removed
CLL                   (was used to switch on WHENIMD) 
CLL                                    S.J.Swarbrick
CLL   4.3    Feb. 97  T3E migration: pass recip_pstar to LAYER_DD :
CLL                    recip_pstar is compressed in the same way as
CLL                    pstar before being passed to LAYER_DD.
CLL                                    S.J.Swarbrick
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!LL  4.5   20/02/98  Remove redundant code. A. Dickinson
CLL
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DD_CALL (NP_FIELD,NPNTS,KCT,THP,QP,THE,QE,DTHBYDT,
     *                    DQBYDT,FLX,PSTAR,AK,BK,AKM12,BKM12,DELAK,
     *                    DELBK,EXNER,PRECIP,RAIN,SNOW,ICCB,ICCT,
     *                    BWATER,BTERM,BGMK,TIMESTEP,CCA,NTERM,
     *                    recip_pstar)    
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
C
      INTEGER I                  ! LOOP COUNTER
C
      INTEGER K                  ! PRESENT MODEL LAYER
C
      INTEGER NPNTS              ! IN NUMBER OF POINTS
C
      INTEGER NDD,NTERM          ! COMPRESSED VECTOR LENGTH FOR
                                 ! DOWNDRAUGHT CALCULATION
C
      INTEGER NP_FIELD           ! IN FULL VECTOR LENGTH
C
      INTEGER NDDON_TMP          ! NUMBER OF POINTS WITH ACTIVE
                                 ! DOWNDRAUGHT
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      INTEGER KCT                ! IN CONVECTIVE CLOUD TOP LAYER
C
      REAL AK(KCT+1)             ! IN ) HYBRID CO-ORDINATE VALUES AT
      REAL BK(KCT+1)             ! IN ) MID-LAYER OF LAYER K
C
      REAL AKM12(KCT+2)          ! IN ) HYBRID CO-ORDINATE VALUES AT
      REAL BKM12(KCT+2)          ! IN ) LOWER LAYER BOUNDARY OF LAYER K
C
      REAL DELAK(KCT+1)          ! IN ) HYBRID CO-ORDINATE VALUES FOR
      REAL DELBK(KCT+1)          ! IN ) THICKNESS OF LAYER K
C
      REAL EXNER(NP_FIELD,KCT+2) ! IN EXNER FUNCTION AT LAYER BOUNDARIES
                                 !    STARTING AT LEVEL K-1/2
C
      REAL THP(NPNTS,KCT+1)      ! IN POTENTIAL TEMPERATURE OF
                                 !    PARCEL (K)
C
      REAL QP(NPNTS,KCT+1)       ! IN MODEL MIXING RATIO (KG/KG)
C
      REAL THE(NP_FIELD,KCT+1)   ! IN MODEL ENVIRONMENTAL POTENTIAL
                                 !    TEMPERATURE (K)
C
      REAL QE(NP_FIELD,KCT+1)    ! IN ENVIRONMENT MIXING RATIO
                                 !    (KG/KG)
C
      REAL FLX(NPNTS,KCT+1)      ! IN CONVECTIVE MASSFLUX (PA/S)
C
      REAL PSTAR(NP_FIELD)       ! IN SURFACE PRESSURE (PA)
C
      REAL PRECIP(NPNTS,KCT+1)   ! IN PRECIPITATION ADDED WHEN
                                 !    DESCENDING FROM LAYER K TO K-1
                                 !    (KG/M**2/S)
C
      INTEGER ICCB(NP_FIELD)     ! IN CLOUD BASE LEVEL
C
      INTEGER ICCT(NP_FIELD)     ! IN CLOUD TOP LEVEL
C
      REAL CCA(NP_FIELD)         ! IN CONVECTIVE CLOUD AMOUNT
C
      LOGICAL BWATER(NPNTS,2:KCT+1)!IN  MASK FOR THOSE POINTS AT WHICH
                                   !     CONDENSATE IS WATER IN LAYER K
C
      LOGICAL BTERM(NPNTS)       ! IN MASK FOR THOSE POINTS WHERE
                                 !    UPDRAUGHT IS TERMINATING
C
      LOGICAL BGMK(NPNTS)        ! IN MASK FOR POINTS WHERE PARCEL IN
                                 !    LAYER K IS SATURATED
C
      REAL TIMESTEP
      REAL recip_PSTAR(NP_FIELD)! Reciprocal of pstar array 
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      REAL DTHBYDT(NP_FIELD,KCT+1) ! INOUT
                                 ! IN  INCREMENT TO MODEL POTENTIAL
                                 !     TEMPERATURE (K/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     POTENTIAL TEMPERATURE (K/S)
C
      REAL DQBYDT(NP_FIELD,KCT+1) ! INOUT
                                 ! IN  INCREMENT TO MODEL MIXING
                                 !     RATIO (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     MIXING RATIO (KG/KG/S)
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL RAIN(NP_FIELD)   ! OUT RAINFALL AT SURFACE (KG/M**2/S)
C
      REAL SNOW(NP_FIELD)   ! OUT SNOWFALL AT SURFACE (KG/M**2/S)
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
C
      REAL EXNER_KM12_C(NTERM)   ! COMPRESSED EXNER FUNCTION AT
                                 ! LAYER K
C
      REAL EXNER_KP12_C(NTERM)   ! COMPRESSED EXNER FUNCTION AT
                                 ! LAYER K+1
C
      REAL EXNER_KM32_C(NTERM)   ! COMPRESSED EXNER FUNCTION AT
                                 ! LAYER K-1
C
      REAL PK(NTERM)             ! PRESSURE OF LAYER K (PA)
C
      REAL P_KM1(NTERM)          ! PRESSURE OF LAYER K-1 (PA)
C
      REAL EXK(NTERM)            ! EXNER RATIO FOR LAYER K
C
      REAL EXKM1(NTERM)          ! EXNER RATIO FOR LAYER K-1
C
      REAL DELPK(NTERM)          ! PRESSURE DIFFERENCE ACROSS LAYER K
                                 ! (PA)
C
      REAL DELPKM1(NTERM)        ! PRESSURE DIFFERENCE ACROSS
                                 ! LAYER K-1 (PA)
C
      REAL AMDETK(NTERM)         ! MIXING DETRAINMENT AT LEVEL K
                                 ! MULTIPLIED BY APPROPRIATE LAYER
                                 ! THICKNESS
C
      REAL EKM12(NTERM)          ! EXNER RATIO AT LAYER K-1/2
C
      REAL EKM14(NTERM)          ! EXNER RATIO AT LAYER K-1/4
C
      REAL EKM34(NTERM)          ! EXNER RATIO AT LAYER K-3/4
C
      LOGICAL BWATER_K_C(NTERM)  ! COMPRESSED MASK FOR THOSE
                                 ! POINTS AT WHICH CONDENSATE
                                 ! IS WATER IN LAYER K
C
      REAL PRECIP_K_C(NTERM)     ! COMPRESSED PRECIPITATION
                                 ! ADDED WHEN DESCENDING FROM
                                 ! LAYER K TO K-1 (KG/M**2/S)
C
      REAL Q_K_C(NTERM)          ! COMPRESSED PARCEL MIXING RATIO
                                 ! OF LAYER K (KG/KG)
C
      REAL TH_K_C(NTERM)         ! COMPRESSED PARCEL POTENTIAL
                                 ! TEMPERATURE OF LAYER K (K)
C
      REAL PSTAR_C(NTERM)        ! COMPRESSED SURFACE PRESSURE (PA)
C
      REAL recip_PSTAR_C(NTERM)  ! Reciprocal of comp. pstar array
C                                                                       
      INTEGER ICCB_C(NTERM)      ! COMPRESSED CLOUD BASE LEVEL
C
      REAL DTHBYDT_K_C(NTERM)    ! COMPRESSED INCREMENT TO MODEL
                                 ! POTENTIAL TEMPERATURE OF LAYER K
                                 ! (K/S)
C
      REAL DTHBYDT_KM1_C(NTERM)  ! COMPRESSED INCREMENT TO MODEL
                                 ! POTENTIAL TEMPERATURE OF LAYER K-1
                                 ! (K/S)
C
      REAL DQBYDT_K_C(NTERM)     ! COMPRESSED INCREMENT TO MODEL
                                 ! MIXING RATIO OF LAYER K (KG/KG/S)
C
      REAL DQBYDT_KM1_C(NTERM)   ! COMPRESSED INCREMENT TO MODEL
                                 ! MIXING RATIO OF LAYER K-1 (KG/KG/S)
C
      REAL DELTD(NTERM)          ! COOLING NECESSARY TO
                                 ! ACHIEVE SATURATION (K)
C
      REAL DELQD(NTERM)          ! MOISTENING NECESSARY TO
                                 ! ACHIEVE SATURATION (KG/KG)
C
      REAL QDD_K(NTERM)          ! MIXING RATIO OF DOWNDRAUGHT IN
                                 ! LAYER K (KG/KG)
C
      REAL THDD_K(NTERM)         ! MODEL POTENTIAL TEMPERATURE
                                 ! OF DOWNDRAUGHT IN LAYER K (K)
C
      REAL FLX_DD_K(NPNTS)       ! DOWNDRAUGHT INITIAL MASS FLUX
                                 ! (PA/S)
C
      REAL FLX_DD_K_C(NTERM)     ! COMPRESSED DOWNDRAUGHT INITIAL
                                 ! MASS FLUX (PA/S)
C
      LOGICAL BDDI(NPNTS)        ! MASK FOR POINTS WHERE DOWNDRAUGHT
                                 ! MIGHT OCCUR
C
      LOGICAL BDDI_C(NTERM)      ! COMPRESSED MASK FOR POINTS WHERE
                                 ! DOWNDRAUGHT MAY INITIATE
C
      INTEGER INDEX1(NTERM)      ! INDEX FOR COMPRESS AND EXPAND
C
      REAL QE_K_C(NTERM)         ! COMPRESSED ENVIRONMENT MIXING
                                 ! RATIO OF LAYER K (KG/KG)
C
      REAL QE_KM1_C(NTERM)       ! COMPRESSED ENVIRONMENT MIXING
                                 ! RATIO OF LAYER K-1 (KG/KG)
C
      REAL THE_K_C(NTERM)        ! COMPRESSED POTENTIAL TEMPERATURE
                                 ! OF ENVIRONMENT IN LAYER K (K)
C
      REAL THE_KM1_C(NTERM)      ! COMPRESSED POTENTIAL TEMPERATURE
                                 ! OF ENVIRONMENT IN LAYER K-1 (K)
C
      REAL RAIN_C(NTERM)         ! COMPRESSED SURFACE RAINFALL
                                 ! (KG/M**2/S)
C
      REAL SNOW_C(NTERM)         ! COMPRESSED SURFACE SNOWFALL
                                 ! (KG/M**2/S)
C
      REAL FLX_UD_K_C(NTERM)     ! UPDRAUGHT MASS FLUX AT LAYER K
C
      REAL RAIN_ENV(NTERM)       ! AMOUNT OF RAINFALL PASSING THROUGH
                                 ! ENVIRONMENT (KG/M**2/S)
C
      REAL SNOW_ENV(NTERM)       ! AMOUNT OF SNOWFALL PASSING THROUGH
                                 ! ENVIRONMENT (KG/M**2/S)
C
      REAL RAIN_DD(NTERM)        ! AMOUNT OF RAINFALL PASSING THROUGH
                                 ! DOWNDRAUGHT (KG/M**2/S)
C
      REAL SNOW_DD(NTERM)        ! AMOUNT OF SNOWFALL PASSING THROUGH
                                 ! DOWNDRAUGHT (KG/M**2/S)
C
      LOGICAL BDD_START(NPNTS)   ! MASK FOR THOSE POINT WHERE
                                 ! DOWNDRAUGHT IS ABLE TO START
                                 ! FROM LEVEL K
C
      LOGICAL BDD_START_C(NTERM) ! COMPRESSED MASK FOR THOSE POINT
                                 ! WHERE DOWNDRAUGHT IS ABLE TO START
                                 ! FROM LEVEL K
C
      LOGICAL BDDWT_K(NPNTS)     ! MASK FOR POINTS IN DOWNDRAUGHT
                                 ! WHERE PPT IN LAYER K IS LIQUID
C
      LOGICAL BDDWT_K_C(NTERM)   ! COMPRESSED MASK FOR POINTS IN DD
                                 ! WHERE PPT IN LAYER K IS LIQUID
C
      LOGICAL BDDWT_KM1(NPNTS)   ! MASK FOR POINTS IN DOWNDRAUGHT
                                 ! WHERE PPT IN LAYER K-1 IS LIQUID
C
      LOGICAL BDDWT_KM1_C(NTERM) ! COMPRESSED MASK FOR POINTS IN DD
                                 ! WHERE PPT IN LAYER K-1 IS LIQUID
C
      LOGICAL BDD_ON(NPNTS)      ! MASK FOR THOSE POINTS WHERE DD
                                 ! CONTINUES FROM LAYER K+1
C
      LOGICAL BDD_ON_C(NTERM)    ! COMPRESSED MASK FOR POINTS WHERE DD
                                 ! CONTINUES FROM LAYER K+1
C
      INTEGER KMIN(NTERM)        ! FREEZING LEVEL WHERE ENTRAINMENT
                                 ! RATES ARE INCREASED
C
      REAL FLX_STRT(NPNTS)       ! MASSFLUX AT LEVEL WHERE DOWNDRAUGHT
                                 ! STARTS (PA/S)
C
      REAL FLX_STRT_C(NTERM)     ! COMPRESSED VALUE OF FLX_STRT
C
      REAL CCA_C(NTERM)          ! COMPRESSED CONVECTIVE CLOUD AMOUNT
C
      INTEGER INDEX2(NTERM)      ! INDEX OF WHERE ACTICE DOWNDRAUGHT
                                 ! OCCURS
C
      REAL LR_UD_REF(NTERM)      ! PRECIPITATION MIXING RATIO AT LOWEST
                                 ! PRECIPITATING LEVEL OF UD
C
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C-----------------------------------------------------------------------
C
      EXTERNAL FLX_INIT, LAYER_DD, DD_INIT, DOWND
C
C-----------------------------------------------------------------------
C CALCULATE INDEX FOR COMPRESS ON BASIS OF BTERM
C-----------------------------------------------------------------------
C
      NDD = 0
      DO I=1,NPNTS
       IF (BTERM(I)) THEN
          NDD = NDD+1
          INDEX1(NDD) = I
       END IF
      END DO
C
C----------------------------------------------------------------------
C INITIALISE LOGICAL ARRAYS AS FALSE
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       BDDI(I) = .FALSE.
       BDD_START(I) = .FALSE.
       BDDWT_K(I) = .FALSE.
       BDDWT_KM1(I) = .FALSE.
       BDD_ON(I) = .FALSE.
C
C-----------------------------------------------------------------------
C CALCULATE MASK FOR THOSE POINT WHERE DOWNDRAUGHT MIGHT OCCUR
C AND LEVEL AT WHICH IT MIGHT INITIATE
C-----------------------------------------------------------------------
C
        IF (KCT .GE. 4 .AND. BTERM(I) .AND. BGMK(I) .AND. (KCT-ICCB(I))
     &       .GT. 2)  BDDI(I) = .TRUE.
      END DO
C
C----------------------------------------------------------------------
C CALCULATE INITIAL DOWNDRAUGHT MASS FLUX
C-----------------------------------------------------------------------
C
      IF (KCT .GE. 4)
     *  CALL FLX_INIT (NPNTS,KCT,ICCB,ICCT,FLX,FLX_DD_K,BDDI,FLX_STRT)
C
C-----------------------------------------------------------------------
C COMPRESS ALL INPUT ARRAYS FOR THE DOWNDRAUGHT CALCULATION
C-----------------------------------------------------------------------
C
      DO 10 K = KCT+1,2,-1
C
         DO I=1,NDD
            TH_K_C(I) = THP(INDEX1(I),K)
            Q_K_C(I) = QP(INDEX1(I),K)
            THE_K_C(I) = THE(INDEX1(I),K)
            THE_KM1_C(I) = THE(INDEX1(I),K-1)
            QE_K_C(I) = QE(INDEX1(I),K)
            QE_KM1_C(I) = QE(INDEX1(I),K-1)
            DTHBYDT_K_C(I) = DTHBYDT(INDEX1(I),K)
            DTHBYDT_KM1_C(I) = DTHBYDT(INDEX1(I),K-1)
            DQBYDT_K_C(I) = DQBYDT(INDEX1(I),K)
            DQBYDT_KM1_C(I) = DQBYDT(INDEX1(I),K-1)
            EXNER_KM12_C(I) = EXNER(INDEX1(I),K)
            EXNER_KP12_C(I) = EXNER(INDEX1(I),K+1)
            EXNER_KM32_C(I) = EXNER(INDEX1(I),K-1)
            PRECIP_K_C(I) = PRECIP(INDEX1(I),K)
            FLX_UD_K_C(I) = FLX(INDEX1(I),K)
            BWATER_K_C(I) = BWATER(INDEX1(I),K)
         END DO
         IF (K.EQ.KCT+1) THEN
          DO I=1,NDD
            FLX_DD_K_C(I) = FLX_DD_K(INDEX1(I))
            FLX_STRT_C(I) = FLX_STRT(INDEX1(I))
            PSTAR_C(I) = PSTAR(INDEX1(I))
            recip_pstar_c(I)=recip_pstar(index1(I))
            ICCB_C(I) = ICCB(INDEX1(I))
            BDDI_C(I) = BDDI(INDEX1(I))
            BDD_START_C(I) = BDD_START(INDEX1(I))
            RAIN_C(I) = RAIN(INDEX1(I))
            SNOW_C(I) = SNOW(INDEX1(I))
            BDDWT_K_C(I) = BDDWT_K(INDEX1(I))
            BDDWT_KM1_C(I) = BDDWT_KM1(INDEX1(I))
            BDD_ON_C(I) = BDD_ON(INDEX1(I))
            CCA_C(I) = CCA(INDEX1(I))
            LR_UD_REF(I) = 0.0
          END DO
         END IF
C
C----------------------------------------------------------------------
C IF BELOW CONVECTIVE CLOUD BASE DOWNDRAUGHT NOT ALLOWED TO FORM
C----------------------------------------------------------------------
C
      DO I=1,NDD
       IF (K.LT.ICCB_C(I)) BDDI_C(I)=.FALSE.
      END DO
C
C-----------------------------------------------------------------------
C RESET EN/DETRAINMENT RATES FOR DOWNDRAUGHT
C-----------------------------------------------------------------------
C
      CALL LAYER_DD (NDD,K,KCT,THE_K_C,THE_KM1_C,FLX_STRT_C,AK,BK,
     *               AKM12,BKM12,DELAK,DELBK,EXNER_KM12_C,EXNER_KP12_C,
     *               EXNER_KM32_C,PSTAR_C,PK,P_KM1,DELPK,DELPKM1,EXK,
     *               EXKM1,AMDETK,EKM14,EKM34,KMIN,BDDI_C,
     *               recip_pstar_c)   
C
C-----------------------------------------------------------------------
C INITIALISE DOWNDRAUGHT
C DOWNDRAUGHT NOT ALLOWED TO FORM FROM CLOUD TOP LAYER (KCT+1)
C OR FROM BELOW CLOUD BASE
C-----------------------------------------------------------------------
C
      IF (KCT .GE. 4 .AND. K.LT.KCT+1)
     & CALL DD_INIT(NDD,TH_K_C,Q_K_C,THE_K_C,QE_K_C,PK,EXK,THDD_K,
     &              QDD_K,DELTD,DELQD,BDD_START_C,K,BDDI_C,BDD_ON_C)
C
C-----------------------------------------------------------------------
C UPDATE MASK FOR WHERE DOWNDRAUGHT OCCURS
C-----------------------------------------------------------------------
C
      DO I=1,NDD
        IF (BDD_START_C(I).OR.BDD_ON_C(I)) BDD_ON_C(I)=.TRUE.
      END DO
C
      NDDON_TMP = 0
      DO I=1,NDD
        IF (BDD_ON_C(I)) THEN
          NDDON_TMP = NDDON_TMP+1
        END IF
      END DO
C
C-----------------------------------------------------------------------
C CALL DOWNDRAUGHT ROUTINE
C-----------------------------------------------------------------------
C

      CALL DOWND(NDD,K,KCT,THDD_K,QDD_K,THE_K_C,THE_KM1_C,QE_K_C,
     &           QE_KM1_C,DTHBYDT_K_C,DTHBYDT_KM1_C,DQBYDT_K_C,
     &           DQBYDT_KM1_C,FLX_DD_K_C,P_KM1,DELPK,DELPKM1,EXK,
     &           EXKM1,DELTD,DELQD,AMDETK,EKM14,EKM34,PRECIP_K_C,
     &           RAIN_C,SNOW_C,ICCB_C,BWATER_K_C,BDD_START_C,
     &           BDDWT_K_C,BDDWT_KM1_C,BDD_ON_C,RAIN_ENV,SNOW_ENV,
     &           RAIN_DD,SNOW_DD,FLX_UD_K_C,TIMESTEP,CCA_C,NDDON_TMP,
     &           LR_UD_REF)
C
C-----------------------------------------------------------------------
C DECOMPRESS/EXPAND THOSE VARIABLES WHICH ARE TO BE OUTPUT
C-----------------------------------------------------------------------
C
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO I=1,NDD
         DTHBYDT(INDEX1(I),K) = DTHBYDT_K_C(I)
         DTHBYDT(INDEX1(I),K-1) = DTHBYDT_KM1_C(I)
         DQBYDT(INDEX1(I),K) = DQBYDT_K_C(I)
         DQBYDT(INDEX1(I),K-1) = DQBYDT_KM1_C(I)
         IF (K.EQ.2) THEN
          RAIN(INDEX1(I)) = RAIN_C(I)
          SNOW(INDEX1(I)) = SNOW_C(I)
         END IF
         PRECIP(INDEX1(I),K) = PRECIP_K_C(I)
        END DO
C
C----------------------------------------------------------------------
C   END OF MAIN K LOOP
C----------------------------------------------------------------------
C
 10   CONTINUE
C
      RETURN
      END
C
