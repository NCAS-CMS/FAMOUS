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
CLL  SUBROUTINE TROP--------------------------------------------------
CLL
CLL  Purpose:  Calculates tropopause temperature and pressure.
CLL            (This version uses approximation that z is linear with
CLL            Exner function within a layer)
CLL            Note also routine TROPIN - any physical changes to one
CLL            should be considered for mirroring in the other.
CLL
CLL AD,DR,RS    <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.2  19/04/93  Code for new real missing data indicator (TCJ).
CLL
CLL   3.3  25/11/93  Ensure no overshoot for trop height  (C Wilson)
CLL   4.2  Nov. 96   T3E migration: Hf functions & *DEF CRAY replaced
CLL                  by T3E function rtor_v & *DEF T3E; code
CLL                  restructured appropriately.  S.J.Swarbrick
CLL   4.4  18/06/97  Initialise P_EXNER_T. Set PT to RMDI if
CLL                  Tropopause not found. S. Lorrimer
CLL
CLL  Logical components covered: P3
CLL
CLL  Project task:
CLL
CLL  Documentation: The interpolation formulae are described in
CLL                 unified model on-line documentation paper S1.
CLL
CLLEND----------------------------------------------------------------
C
C*L  Arguments:-------------------------------------------------------
      SUBROUTINE TROP
     *(PSTAR,THETA,P_EXNER_HALF,ZH,TT,PT,ZT,POINTS,P_LEVELS,
     * MIN_TROP_LEVEL,AKH,BKH)

      IMPLICIT NONE

      INTEGER
     * POINTS    !IN Number of points to be processed
     *,P_LEVELS  !IN Number of model levels
     *,MIN_TROP_LEVEL  !IN Level no. for lowest possible tropopause.
     *                 !   Set to first level above boundary layer

      REAL
     * PSTAR(POINTS)          !IN P Star
     *,THETA(POINTS,P_LEVELS) !IN Potential temperature at full levels
     *,P_EXNER_HALF(POINTS,P_LEVELS+1) !IN Exner pressure at model
     *                                 !   half levels
     *,ZH(POINTS,P_LEVELS+1)  !IN Height of model half levels
     *,AKH(P_LEVELS+1)        !IN Hybrid Coords. A values for half levs.
     *,BKH(P_LEVELS+1)        !IN Hybrid Coords. B values for half levs.
     *,TT(POINTS)             !OUT Temperature of tropopause
     *,PT(POINTS)             !OUT Pressure of tropopause
     *,ZT(POINTS)             !OUT Height of tropopause

C Workspace usage:-----------------------------------------------------
       REAL LAPSE_RATE(POINTS,MIN_TROP_LEVEL:P_LEVELS)
       LOGICAL LTROP(POINTS)
C External subroutines called:-----------------------------------------
C*---------------------------------------------------------------------
C Define local variables:----------------------------------------------
      INTEGER I,J,JP1
      REAL PJP1,PJ,PJM1  !  Pressures at half levels J+1/J/J-1
      REAL P_EXNER_FULL_J,P_EXNER_FULL_JM1
     *,DEL_EXNER_J,DEL_EXNER_JM1,TERM1,TERM2
     *,ZDT    
     *,ZDJM1,ZDJ
C     REAL P_EXNER_T 
      REAL P_EXNER_T(POINTS)
      REAL POWER(POINTS)
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

CL 1. Set up local constants and initialise arrays

      REAL CP_OVER_G,ONE_OVER_KAPPA,P_EXNER_500,P_EXNER_50
      PARAMETER(CP_OVER_G=CP/G)
      PARAMETER(ONE_OVER_KAPPA=1./KAPPA)

C----------------------------------------------------------------------

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


      P_EXNER_500=(500./1000.)**KAPPA
      P_EXNER_50=(50./1000.)**KAPPA

      DO 100 I=1,POINTS

C Initialise logical string which indicates whether tropopause found
C at a lower level: LTROP=T is not found; LTROP=F is found.

      LTROP(I)=.TRUE.

C Initialise tropopause height and pressure to missing data

      P_EXNER_T(I)=0
      TT(I)=RMDI
      ZT(I)=RMDI

100   CONTINUE

CL 2. Compute lapse rate between full levels: equation (3.16)


      DO 200 J=MIN_TROP_LEVEL,P_LEVELS
      DO 210 I=1,POINTS

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

CL 3. Calculate tropopause temperature, height  and pressure

      DO 300 J=MIN_TROP_LEVEL+1,P_LEVELS

C 'J+1' level for lapse rate test; allows J iteration up to P_LEVELS
      JP1=MIN(J+1,P_LEVELS)

      DO 310 I=1,POINTS

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
     * LAPSE_RATE(I,J).LT.LAPSE_TROP.AND.
     * LAPSE_RATE(I,JP1).LT.LAPSE_TROP.AND.
     * LTROP(I))THEN

C Reset logical string to say tropopause now found
      LTROP(I)=.FALSE.

C Z(j-1)-Z(j-1/2); Z(j)-Z(j-1/2)
      ZDJM1=CP_OVER_G*THETA(I,J-1)*(P_EXNER_HALF(I,J)-P_EXNER_FULL_JM1)
      ZDJ=CP_OVER_G*THETA(I,J)*(P_EXNER_HALF(I,J)-P_EXNER_FULL_J)

C Z(tropopause) - Z(j-1/2): equation (3.19)
      ZDT=(THETA(I,J-1)*P_EXNER_FULL_JM1-THETA(I,J)*P_EXNER_FULL_J
     *+LAPSE_RATE(I,J-1)*ZDJM1
     *-LAPSE_RATE(I,JP1)*ZDJ)
     */(LAPSE_RATE(I,J-1)-LAPSE_RATE(I,JP1))

C Ensure trop level doesn't undershoot Z(j-1)
      ZDT=MAX(ZDT,ZDJM1)

C Ensure trop level doesn't overshoot Z(j)
      ZDT=MIN(ZDT,ZDJ)

C Tropopause height
      ZT(I)=ZDT+ZH(I,J)

C Tropopause temperature : equation (3.20)
      TT(I)=THETA(I,J)*P_EXNER_FULL_J
     *     -LAPSE_RATE(I,JP1)*(ZDT-ZDJ)

C Exner pressure of tropopause: equation (3.22)
      IF(ZDT.GT.0.0)THEN
        P_EXNER_T(I)=P_EXNER_HALF(I,J)-G*ZDT/(CP*THETA(I,J))      
      ELSE
        P_EXNER_T(I)=P_EXNER_HALF(I,J)-G*ZDT/(CP*THETA(I,J-1))        
      ENDIF

      ENDIF

310   CONTINUE
300   CONTINUE

C Pressure of tropopause: equation (3.21)
      DO I=1,POINTS
        P_EXNER_T(I)=P_EXNER_T(I)**ONE_OVER_KAPPA
      END DO
      DO I=1,POINTS
        PT(I)=PREF*P_EXNER_T(I)
        IF (LTROP(I)) THEN 
          PT(I)=RMDI
        ENDIF
 
      END DO

      RETURN
      END
