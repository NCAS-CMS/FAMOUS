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
CLL  SUBROUTINE EVP----------------------------------------------------
CLL
CLL  PURPOSE : CALCULATES THE EVAPORATION OF PRECIPITATION
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 23/7/92
CLL
CLL  LOGICAL COMPONENTS COVERED:
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C      Vn.4.2   Oct. 96  T3E migration: *DEF CRAY removed, HF functions
C                         replaced.
C                                      S.J.Swarbrick
CLL  4.3    Feb. 97   T3E optimisation of powers & sqrt
CLL                                  D.Salmond & S.J.Swarbrick
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE EVP(NPNTS,PRECIP,TEVP,CCA,RHO,DELQ,DELPKM1,EVAP,
     &               BEVAP,IPHASE,AREA_FAC)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS USED IN THIS SUBROUTINE
C-----------------------------------------------------------------------
C
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
      REAL P_LQ1,P_LQ2   ! EXPONENTS USED IN CALCULATION OF
                         ! EVAPORATION OF LIQUID
      PARAMETER (P_LQ1 = 0.42, P_LQ2 = 0.6)
C
      REAL P_ICE1,P_ICE2 ! EXPONENTS USED IN CALCULATION OF
                         ! EVAPORATION OF ICE
      PARAMETER (P_ICE1 = 0.55, P_ICE2 = 0.76)
C
      REAL RHO_LQP1,RHO_LQP2, ! EXPONENTS AND CONSTANTS ASSOCIATED
     &     RHO_LQA,RHO_LQB    ! WITH DENSITY TERM IN EVAPORATION
                              ! OF LIQUID
      PARAMETER (RHO_LQP1=0.21, RHO_LQP2=0.55, RHO_LQA=67.08,
     &           RHO_LQB=541.16)
C
      REAL RHO_ICP1,RHO_ICP2, ! EXPONENTS AND CONSTANTS ASSOCIATED
     &     RHO_ICEA,RHO_ICEB  ! WITH DENSITY TERM IN EVAPORATION
                              ! OF ICE
      PARAMETER (RHO_ICP1=0.28, RHO_ICP2=0.63, RHO_ICEA=1569.52,
     &           RHO_ICEB=32069.02)
C
      REAL LQ_A,LQ_B,LQ_C ! CONSTANTS USED IN QUADRATIC FORMULA
                          ! FOR EVAPORATION OF LIQUID
      PARAMETER ( LQ_A = 2.008E-9, LQ_B = -1.385E-6, LQ_C = 2.424E-4)
C
      REAL ICE_A,ICE_B,ICE_C ! CONSTANTS USED IN QUADRATIC FORMULA
                             ! FOR EVAPORATION OF ICE
      PARAMETER ( ICE_A = -5.2E-9, ICE_B = 2.5332E-6, ICE_C = -2.911E-4)
C
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER I                 ! LOOP COUNTER
C
      INTEGER NPNTS             ! IN VECTOR LENGTH
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
C
      REAL DELQ(NPNTS)          ! IN CHANGE IN HUMIDITY MIXING
                                !    RATIO ACROSS LAYER K (KG/KG)
C
      REAL TEVP(NPNTS)          ! IN TEMPERATURE OF LAYER K (K)
C
      LOGICAL BEVAP(NPNTS)      ! IN MASK FOR POINTS WHERE EVAPORATION
                                !    TAKES PLACE
C
      REAL PRECIP(NPNTS)        ! IN AMOUNT OF PRECIPITATION(KG/M**2/S)
C
      REAL DELPKM1(NPNTS)       ! IN CHANGE IN PRESSURE ACROSS
                                !    LAYER K-1
C
      REAL CCA(NPNTS)           ! IN CONVECTIVE CLOUD AMOUNT
C
      REAL RHO(NPNTS)           ! IN DENSITY OF AIR
C
      INTEGER IPHASE            ! IN INDICATION FOR RAIN (1), OR
                                !    SNOW (2)
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL EVAP(NPNTS)   ! OUT EVAPORATION
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE LOCALLY DEFINED
C-----------------------------------------------------------------------
C
      REAL ECON          ! QUADRATIC TERM
C
      REAL C1            ! CONSTANT
C
      REAL C2            ! CONSTANT
C
      REAL SR_RHO        ! SQUARE ROOT OF DENSITY
C
      REAL LRATE         ! LOCAL RATE OF PRECIPITATION
C
      REAL CA            ! LOCAL CLOUD AREA
C
      REAL AREA_FAC      ! FRACTION OF CONVECTIVE CLOUD AMOUNT TO GIVE
                         ! LOCAL CLOUD AREA
      real work1(npnts),work2(npnts),work3(npnts),
     1     r_work1(npnts),r_work2(npnts),
     1     r_rho(npnts)
      integer kk
      real tl1,ti1
C
C-----------------------------------------------------------------------
C START OF ROUTINE
C-----------------------------------------------------------------------
C
      tl1=0.5*P_LQ1
      ti1=0.5*P_ICE1

      kk=0
      do i=1,npnts
      IF (BEVAP(I).AND.PRECIP(I) .GT. 0.0)THEN
      kk=kk+1
           CA = AREA_FAC*CCA(I)
           LRATE = PRECIP(I)/CA
      work1(kk)=LRATE*LRATE*RHO(I)
      work2(kk)=LRATE
      work3(kk)=RHO(I)
      ENDIF
      enddo

      IF (IPHASE.EQ.1) THEN        ! RAIN
C
      do i=1,kk
        r_work1(i)=work1(i)**tl1
        r_work2(i)=work2(i)**P_LQ2
        r_rho  (i)=work3(i)**RHO_LQP2
      end do
      kk=0
      DO I=1,NPNTS
       IF (BEVAP(I)) THEN
         IF (PRECIP(I) .GT. 0.0) THEN
           kk=kk+1
           ECON = ((LQ_A*TEVP(I)+LQ_B)*TEVP(I)+LQ_C)
           CA = AREA_FAC*CCA(I)
           LRATE = PRECIP(I)/CA
           C1 = RHO_LQA*CA*r_work1(kk)
           C2 = RHO_LQB*CA*r_work2(kk)*r_rho(kk)
           EVAP(I) = MIN(ECON*(C1+C2)*DELQ(I)*DELPKM1(I)/G,LRATE)
         ELSE
           EVAP(I) = 0.0
         END IF
       END IF
      END DO
C
      ELSE IF (IPHASE.EQ.2) THEN        ! SNOW
C
      do i=1,kk
        r_work1(i)=work1(i)**ti1
        r_work2(i)=work2(i)**P_ICE2
        r_rho  (i)=work3(i)**RHO_ICP2
      end do
C
      kk=0
      DO I=1,NPNTS
       IF (BEVAP(I)) THEN
         IF (PRECIP(I) .GT. 0.0) THEN
         kk=kk+1
           ECON = ((ICE_A*TEVP(I)+ICE_B)*TEVP(I)+ICE_C)
           CA = AREA_FAC*CCA(I)
           LRATE = PRECIP(I)/CA
           C1 = RHO_ICEA*CA*r_work1(kk)
           C2 = RHO_ICEB*CA*r_work2(kk)*r_rho(kk)
           EVAP(I) = MIN(ECON*(C1+C2)*DELQ(I)*DELPKM1(I)/G,LRATE)
         ELSE
           EVAP(I) = 0.0
         END IF
       END IF
      END DO
C
      ENDIF
C
      RETURN
      END
C
