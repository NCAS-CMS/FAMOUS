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
CLL  SUBROUTINE QCTROP--------------------------------------------------
CLL
CLL   Purpose: Calculates tropopause temperature and pressure.
CLL            (This version uses approximation that z is linear with
CLL            Exner function within a layer)
CLL            Quality control checks output and substitutes new
CLL            values where necessary
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.2  19/04/93  Code for new real missing data indicator (TCJ).
CLL   4.2    Oct. 96  T3E migration: *DEF CRAY code removed - was used
CLL                    to switch on HF functions.
CLL                                   S.J.Swarbrick
CLL   4.4  13/08/97  Limit thickness value to prevent LOG of 
CLL                  negative numbers. D. Robinson.
!LL   4.5  20/04/98  Implement START,END args so that duplicate
!LL                  calculations in the NS halos can be avoided.
!LL                  S.D.Mullerworth
CLL
CLL Programming standard :
CLL
CLL Logical components covered : D442
CLL
CLL Project task :
CLL
CLL  Documentation: The interpolation formulae are described in
CLL                 unified model on-line documentation paper S1.
CLL                 Quality control documentation not yet published
CLL
CLLEND----------------------------------------------------------------
C
C*L  Arguments:-------------------------------------------------------
      SUBROUTINE QCTROP
     &  (THETA,P_HALF,P_EXNER_HALF,ZH,TT,PT,ZT,POINTS,P_LEVELS
     &  ,P_Z,P_T,PSTAR,Q,Q_LEVELS,Z_REF,T_REF
     &  ,MIN_TROP_LEVEL,BOTTOM_QC_LEVEL,TOP_QC_LEVEL,AKH,BKH
     &  ,START,END)

      IMPLICIT NONE

      INTEGER
     * POINTS    !IN Number of points to be processed
     *,P_LEVELS  !IN Number of model levels
     *,Q_LEVELS  !IN Number of wet levels
     *,Z_REF     !IN Reference level for height interpolation
     *,T_REF     !IN Reference level for temperature interpolation
     *,MIN_TROP_LEVEL  !IN Level no. for lowest possible tropopause.
     *                 !   Set to first level above boundary layer
     *,BOTTOM_QC_LEVEL  !IN Lower level no. for quality control scheme
     *,TOP_QC_LEVEL   !IN Upper level no. for quality control scheme
     &,START,END !IN Range of points to calculate

      REAL
     * THETA(POINTS,P_LEVELS) !IN Potential temperature at full levels
     *,P_HALF(POINTS,P_LEVELS+1) !IN
     *,P_EXNER_HALF(POINTS,P_LEVELS+1) !IN Exner pressure at model
     *                                 !   half levels
     *,ZH(POINTS,P_LEVELS+1)  !IN Height of model half levels
     *,P_Z(POINTS)      !IN Pressure at reference level Z_REF
     *,P_T(POINTS)      !IN Pressure at reference level T_REF
     *,PSTAR(POINTS)   !IN Surface pressure
     *,Q(POINTS,Q_LEVELS)     !IN Specific humidity at full levels
     *,AKH(P_LEVELS+1)        !IN Hybrid Coords. A values on half levels
     *,BKH(P_LEVELS+1)        !IN Hybrid Coords. B values on half levels
     *,TT(POINTS)             !OUT Temperature of tropopause
     *,PT(POINTS)             !OUT Pressure of tropopause
     *,ZT(POINTS)             !OUT Height of tropopause

C Workspace usage:-----------------------------------------------------
       REAL LAPSE_RATE(POINTS,MIN_TROP_LEVEL:P_LEVELS)
       LOGICAL LTROP(POINTS)
C Workspace used for quality control:----------------------------------
      INTEGER INDEX(1:TOP_QC_LEVEL-BOTTOM_QC_LEVEL+1)
      REAL WEIGHT_REF(POINTS),PT_QC(POINTS),WEIGHT_LR(POINTS,2)
     *,WEIGHT_LR_INT(POINTS),PT_TT(POINTS),WEIGHT_TT_INT(POINTS)
     *,WEIGHT_TOT_INT(POINTS),P_UPPER(POINTS),P_LOWER(POINTS),SD(POINTS)
     *,PZ(POINTS,2),THICK(POINTS),ZZ(POINTS,2),TT_QC(POINTS),PZS(2)
     *,TH_REF(8),A_REF(8),G_REF(8),SD_REF(8)
C External subroutines called:-----------------------------------------
      EXTERNAL V_INT_Z, V_INT_T
C*---------------------------------------------------------------------
C Define local variables:----------------------------------------------
      INTEGER I,J,K,JP1
      REAL PJP1,PJ,PJM1  !  Pressures at half levels J+1/J/J-1
      REAL P_EXNER_FULL_J,P_EXNER_FULL_JM1
     *,DEL_EXNER_J,DEL_EXNER_JM1,TERM1,TERM2
     *,ZDT,P_EXNER_T
     *,ZDJM1,ZDJ

C*---------------------------------------------------------------------
C Define local variables for quality control:--------------------------
      INTEGER LOOP,JINT,LOOP_INDEX
      REAL FUNC_LR,RIP,P_FACTOR,PSIGWT,TUNER_LR,TUNER_TT
     *,SILLY,SDLIM,TRP_MAX,REAL_INCR
     *,LO_WT
C*---------------------------------------------------------------------
C Define data for quality control:-------------------------------------
      DATA PZS/100000.,50000./
      DATA TH_REF/4380.,4740.,4920.,5100.,5280.,5460.,5640.,5820./
      DATA A_REF/1852.7,1483.7,1457.7,821.62
     &          ,1995.7,2758.5,568.56,1852.7/
      DATA G_REF/-3.0325,-2.4334,-2.2868,-1.0497
     &          ,-3.2610,-4.6646,-0.7939,-3.0325/
      DATA SD_REF/74.728,46.730,38.964,33.777
     &           ,33.252,39.477,8.0700,74.728/
C----------------------------------------------------------------------
C Constants from comdecks:---------------------------------------------
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
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

C*L------------------COMDECK C_LAPSE ----------------------------------
      REAL LAPSE,LAPSE_TROP
      PARAMETER(LAPSE=0.0065)     !  NEAR SURFACE LAPSE RATE
      PARAMETER(LAPSE_TROP=0.002) !  TROPOPAUSE LAPSE RATE
C*----------------------------------------------------------------------
C----------------------------------------------------------------------
      REAL CP_OVER_G,ONE_OVER_KAPPA,P_EXNER_500,P_EXNER_50
      PARAMETER(CP_OVER_G=CP/G)
      PARAMETER(ONE_OVER_KAPPA=1./KAPPA)

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


CL 1. Set up local constants and initialise arrays

      P_EXNER_500=(500./1000.)**KAPPA
      P_EXNER_50=(50./1000.)**KAPPA

      DO I=START,END

C Initialise logical string which indicates whether tropopause found
C at a lower level: LTROP=T is not found; LTROP=F is found.

      LTROP(I)=.TRUE.

C Initialise tropopause height and pressure to missing data

      PT(I)=RMDI
      TT(I)=RMDI
      ZT(I)=RMDI

      ENDDO ! DO I=START,END

CL 2. Compute lapse rate between full levels: equation (3.16)


      DO 200 J=MIN_TROP_LEVEL,P_LEVELS
      DO 210 I=START,END

C Exner pressure at full levels
      PJP1 = AKH(J+1) + BKH(J+1)*PSTAR(I)
      PJ   = AKH(J)   + BKH(J)  *PSTAR(I)
      PJM1 = AKH(J-1) + BKH(J-1)*PSTAR(I)
      P_EXNER_FULL_J = P_EXNER_C
     +(P_EXNER_HALF(I,J+1),P_EXNER_HALF(I,J),PJP1,PJ,KAPPA)
      P_EXNER_FULL_JM1 = P_EXNER_C
     +(P_EXNER_HALF(I,J),P_EXNER_HALF(I,J-1),PJ,PJM1,KAPPA)

C Exner pressure difference across half layers
      DEL_EXNER_J=P_EXNER_HALF(I,J)-P_EXNER_FULL_J
      DEL_EXNER_JM1=P_EXNER_FULL_JM1-P_EXNER_HALF(I,J)

C Denominator
      TERM2=CP_OVER_G*(THETA(I,J-1)*DEL_EXNER_JM1
     *       +THETA(I,J)*DEL_EXNER_J)

C Numerator
      TERM1=THETA(I,J-1)*P_EXNER_FULL_JM1-THETA(I,J)*P_EXNER_FULL_J

C Lapse rate between level j-1 and j
      LAPSE_RATE(I,J)=TERM1/TERM2
210   CONTINUE
200   CONTINUE

CL 2.1. Quality control set up

C Set up constants
      TUNER_LR=1.0
      TUNER_TT=0.5
      SILLY=101320.0
      LOOP_INDEX=9

C Evaluate 1000-500 HPa thickness in metres
      DO 777 K=1,2
        DO 778 I=START,END
          PZ(I,K)=PZS(K)
  778   CONTINUE
        CALL V_INT_Z(PZ(1,K),P_Z,PSTAR,P_EXNER_HALF,
     &               THETA,Q,ZH,ZZ(1,K),POINTS,P_LEVELS,Q_LEVELS,
     &               Z_REF,AKH,BKH,START,END)
  777 CONTINUE
      DO 779 I=START,END
        THICK(I)=ZZ(I,2)-ZZ(I,1)
        IF (THICK(I).GT.6100.0) THEN
          write(6,*) 'QCTROP : Thickness = ',THICK(I),' Reset to 6100.'
          THICK(I)=6100.0
        ENDIF
  779 CONTINUE

C Pick up relevant constants from arrays
      LOOP=1
      DO 211 I=START,END
        SD(I)=SD_REF(LOOP)*100.0
        PT_TT(I)=(A_REF(LOOP)*100.0)+G_REF(LOOP)*10.0*THICK(I)
  211 CONTINUE

      DO 212 LOOP=2,8
      DO 213 I=START,END
        IF(THICK(I).GE.TH_REF(LOOP))THEN
          SD(I)=SD_REF(LOOP)*100.0
          PT_TT(I)=(A_REF(LOOP)*100.0)+G_REF(LOOP)*10.0*THICK(I)
        ENDIF
  213 CONTINUE
  212 CONTINUE

C Set up an array of alternative values using weighting functions
      DO 220 I=START,END
        FUNC_LR=(LAPSE_RATE(I,BOTTOM_QC_LEVEL)-LAPSE_TROP)*1000.0
        IF(FUNC_LR.LT.0.0)THEN
          WEIGHT_LR(I,1)=1.0
        ELSE
          WEIGHT_LR(I,1)=2.0/(EXP(FUNC_LR)+EXP(-FUNC_LR))
        ENDIF
  220 CONTINUE

C Set up default value for weighting function and trop pressure
      DO 230 I=START,END
        WEIGHT_REF(I)=WEIGHT_LR(I,1)          ! Must be set
        PT_QC(I)=SILLY               ! May not need to be set
  230 CONTINUE

C Evaluate weighting functions
      DO 240 J=BOTTOM_QC_LEVEL,TOP_QC_LEVEL
        IF(J.GT.BOTTOM_QC_LEVEL)THEN
          DO 250 I=START,END
            WEIGHT_LR(I,1)=WEIGHT_LR(I,2)
  250     CONTINUE
        ENDIF
        DO 260 I=START,END
          FUNC_LR=(LAPSE_RATE(I,J+1)-LAPSE_TROP)*1000.0
          IF(FUNC_LR.LT.0.0)THEN
            WEIGHT_LR(I,2)=1.0
          ELSE
            WEIGHT_LR(I,2)=2.0/(EXP(FUNC_LR)+EXP(-FUNC_LR))
          ENDIF
  260   CONTINUE

        DO 270 JINT=1,LOOP_INDEX
C Contribution from lapse rate weighting function
          DO 280 I=START,END
            REAL_INCR=(P_HALF(I,J)-P_HALF(I,J+1))/9.0 ! FLOAT(LOOP_INDEX
            P_FACTOR=(REAL_INCR*(JINT-1))/(P_HALF(I,J+1)-P_HALF(I,J))
            WEIGHT_LR_INT(I)=WEIGHT_LR(I,1)+
     &                      (WEIGHT_LR(I,1)-WEIGHT_LR(I,2))*P_FACTOR
            IF(WEIGHT_LR_INT(I).GT.1.0)THEN
              WEIGHT_LR_INT(I)=1.0
            ELSEIF(WEIGHT_LR_INT(I).LT.0.0)THEN
              WEIGHT_LR_INT(I)=0.0
            ENDIF
C Contribution from pseudo-sigma weighting
            PSIGWT=(P_HALF(I,J)-REAL_INCR*(JINT-1))*0.00001
            WEIGHT_LR_INT(I)=(WEIGHT_LR_INT(I)+PSIGWT)*TUNER_LR
C Contribution from thickness weighting function
            RIP=P_HALF(I,J)-REAL_INCR*(JINT-1)
            TERM1=-((RIP-PT_TT(I))*(RIP-PT_TT(I)))
     &            /(2.0*SD(I)*SD(I))
            WEIGHT_TT_INT(I)=EXP(TERM1)*TUNER_TT
            WEIGHT_TOT_INT(I)=WEIGHT_LR_INT(I)+WEIGHT_TT_INT(I)
C Update trop pressure where total weighting function is a maximum
            IF(WEIGHT_REF(I).LT.WEIGHT_TOT_INT(I))THEN
              WEIGHT_REF(I)=WEIGHT_TOT_INT(I)
              PT_QC(I)=RIP
            ENDIF
  280     CONTINUE
  270   CONTINUE
  240 CONTINUE
C Set a silly value if total weighting function is less than 1.0
C (This will nullify any attempt to provide a substitute value)
      DO 290 I=START,END
        IF(WEIGHT_REF(I).LT.1.0)THEN
          PT_QC(I)=SILLY
        ENDIF
  290 CONTINUE

C Obtain corresponding temperatures for alternative tropopause pressures
      CALL V_INT_T(TT_QC,PT_QC,P_T,PSTAR,P_EXNER_HALF
     &  ,THETA,POINTS,P_LEVELS,T_REF,AKH,BKH
     &  ,START,END)
C Estimate range of trop pressures allowed by quality control scheme
      SDLIM=3.0
      DO 2110 I=START,END
        P_UPPER(I)=PT_TT(I)+SDLIM*SD(I)
        P_LOWER(I)=EXP(2.0*ALOG(PT_TT(I))-ALOG(P_UPPER(I)))
 2110 CONTINUE

CL 3. Calculate tropopause temperature, height  and pressure

      LO_WT=0.5

      DO 300 J=MIN_TROP_LEVEL+1,P_LEVELS
      JP1=MIN(J+1,P_LEVELS)
      DO 310 I=START,END

C Exner pressure at full levels
      PJP1 = AKH(J+1) + BKH(J+1)*PSTAR(I)
      PJ   = AKH(J)   + BKH(J)  *PSTAR(I)
      PJM1 = AKH(J-1) + BKH(J-1)*PSTAR(I)
      P_EXNER_FULL_J = P_EXNER_C
     +(P_EXNER_HALF(I,J+1),P_EXNER_HALF(I,J),PJP1,PJ,KAPPA)
      P_EXNER_FULL_JM1 = P_EXNER_C
     +(P_EXNER_HALF(I,J),P_EXNER_HALF(I,J-1),PJ,PJM1,KAPPA)

C Criteria for layer containing tropopause
C (where 'layer' is interval between level j and level j-1)
      IF(P_EXNER_FULL_J.LT.P_EXNER_500.AND.
     * P_EXNER_FULL_JM1.GT.P_EXNER_50.AND.
     * ((LAPSE_RATE(I,J)*LO_WT)+(LAPSE_RATE(I,JP1)*(1.0-LO_WT))).LT.
     * LAPSE_TROP.AND.LTROP(I))THEN

C Reset logical string to say tropause now found
      LTROP(I)=.FALSE.

C Z(j-1)-Z(j-1/2); Z(j)-Z(j-1/2)
      ZDJM1=CP_OVER_G*THETA(I,J-1)*(P_EXNER_HALF(I,J)-P_EXNER_FULL_JM1)
      ZDJ=CP_OVER_G*THETA(I,J)*(P_EXNER_HALF(I,J)-P_EXNER_FULL_J)

C Z(tropopause) - Z(j-1/2): equation (3.19)
      ZDT=(THETA(I,J-1)*P_EXNER_FULL_JM1-THETA(I,J)*P_EXNER_FULL_J
     *+LAPSE_RATE(I,J-1)*ZDJM1
     *-LAPSE_RATE(I,JP1)*ZDJ)
     */MAX(1.E-6,(LAPSE_RATE(I,J-1)-LAPSE_RATE(I,JP1)))

C Ensure trop level doesn't undershoot Z(j-1) (cannot overshoot Z(j) )
      ZDT=MAX(ZDT,ZDJM1)

C Tropopause Height
      ZT(I)=ZDT+ZH(I,J)
      ZDT=MIN(ZDT,ZDJ)

C Tropopause temperature : equation (3.20)
      TT(I)=THETA(I,J)*P_EXNER_FULL_J
     *     -LAPSE_RATE(I,JP1)*(ZDT-ZDJ)

C Exner pressure of tropopause: equation (3.22)
      IF(ZDT.GT.0.0)THEN
        P_EXNER_T=P_EXNER_HALF(I,J)-G*ZDT/(CP*THETA(I,J))
      ELSE
        P_EXNER_T=P_EXNER_HALF(I,J)-G*ZDT/(CP*THETA(I,J-1))
      ENDIF

C Pressure of tropopause: equation (3.21)
      PT(I)=PREF*(P_EXNER_T)**ONE_OVER_KAPPA

      ENDIF

310   CONTINUE
300   CONTINUE

CL 4. Apply quality control to results

      DO 400 I=START,END
        IF(PT(I).LT.P_LOWER(I).AND.PT_QC(I).LT.P_UPPER(I).AND.
     &     PT(I).LT.PT_QC(I))THEN
          PT(I)=PT_QC(I)
          TT(I)=TT_QC(I)
        ENDIF
  400 CONTINUE

C-----------------------------------------------------------------------
C Set arbitrary max tropopause
C-----------------------------------------------------------------------
      TRP_MAX=10100.0
      DO I=START,END
        IF(PT(I).LT.TRP_MAX)THEN
          PT(I)=TRP_MAX
          TT(I)=199
          ZT(I)=16180
        ENDIF
      ENDDO

      RETURN
      END
