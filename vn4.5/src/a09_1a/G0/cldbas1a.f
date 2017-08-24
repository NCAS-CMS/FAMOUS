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
CLL  SUBROUTINE CLOUD_COVER_BASE -------------------------------------
CLL
CLL     PURPOSE:
CLL Return an array holding the lowest cloud base height (Kft) for each
CLL cloud amount (in octas) requested, from model level cloud cover and
CLL convective cloud base level and cover. (Cloud cover is input as a
CLL fraction and cloud base is the height of the half level at the base
CLL of the layer with cloud)
CLL Modification: Also return fraction of air below 1000 ft asl
CLL containing cloud and height of base and top of lowest cloud
CLL layer asl. Lowest cloud layer is defined as lowest set of
CLL contiguous levels with cloud amount greater than threshold.
CLL Fraction set to zero if orography > 1000ft.
CLL     COMMENT:
CLL Since only have Q_LEVELS of CLOUD_FRACTION, do we need P_LEVELS?
CLL
CLL  Model            Modification history:
CLL version  Date
CLL  3.1  05/11/92  New deck author P.Smith
CLL  3.1  20/01/93  New deck - used as mods at 2.7 & 2.8/3.0.
CLL                 Interfacing done by R.T.H.Barnes.
CLL  3.2  02/07/93  Modification to add cloud fraction below 1000 ft
CLL                 and low cloud base and top.
CLL                 Author Pete Clark.
CLL  3.4  24/05/94  Modification to add Wet bulb freezing level height
CLL                 and wet bulb temperature.
CLL                 Author Steve Woltering.
CLL  3.4  06/07/94  Modification to calculate model level heights (not
CLL                 output) and total cloud top height.
CLL                 Author Steve Woltering.
!LL  4.3 26/02/97  Add first & last points to arg.list. RTHBarnes.
CLL
CLL  Programming standard: U M Doc. Paper No. 4
CLL
CLL  Logical components covered :
CLL
CLL  Project task:
CLL
CLL  External documentation  UMDP
CLL
CLLEND-------------------------------------------------------------
C
C*L  Arguments:-------------------------------------------------------
      SUBROUTINE CLOUD_COVER_BASE
     +        (TEMP,Q,P_STAR,P_EXNER_HALF,OROG             !INPUT
     +        ,CONV_CLD_COVER,CONV_BASE_LEV,CLOUD_COVER    !INPUT
     +        ,CONV_TOP_LEV                                !INPUT
     +        ,P_FIELD,P_LEVELS,Q_LEVELS                   !INPUT
     +        ,AK,BK,AKH,BKH,OCTAS,N_OCTAS                 !INPUT
     +        ,CLD_COVER_RQD,LOW_CLD_RQD                   !INPUT
     +        ,WBFL_RQD,WBTEMP_RQD                         !INPUT
     +        ,CLD_TOP_RQD                                 !INPUT
     +        ,CLOUD_BASE                                  !OUTPUT
     +        ,LOW_C_FRAC                                  !OUTPUT
     +        ,LOW_C_BASE                                  !OUTPUT
     +        ,LOW_C_TOP                                   !OUTPUT
     +        ,WBFLH                                       !OUTPUT
     +        ,TW                                          !OUTPUT
     +        ,CLOUD_TOP                                   !OUTPUT
     +        ,FIRST_POINT,LAST_POINT)
      IMPLICIT NONE
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
C*L------------------COMDECK C_KT_FT-----------------------------------
      REAL KT2MS,    ! Knots to m/s conversion
     &     FT2M      ! Feet to meters conversion

      PARAMETER(
     & KT2MS=1852.0/3600.0,
     & FT2M =0.3048)
C*----------------------------------------------------------------------
CLL==========================COMDECK C_LOWCLD==========================
CLL Description:
CLL   This COMDECK contains declarations for ceiling and threshold
CLL constants for lowest cloud layer diagnostics.
CLL
CLL  Model            Modification history:
CLL version  Date
CLL  3.2    19/07/93  First created by Pete Clark.
CLL
CLL Programming Standards U.M.D.P. No. 3 Version 5 Dated 08/12/92
CLL
CLL Logical component:
CLL
CLLEND-------------------------------------------------------------
C Define Parameters:
      REAL
     &   STR_CEIL            ! Max height asl for 'low' cloud (1000ft)
     &   ,CLOUD_THRESHOLD    ! Cloud fraction threshold for low cloud.
      PARAMETER (STR_CEIL=1000.0)
      PARAMETER (CLOUD_THRESHOLD=0.05)
C*--------------------------------------------------------------------
C input variables-----------------------------------------------------
C---------------------------------------------------------------------
      INTEGER
     *        P_FIELD                   ! IN NO. points in field.
     *       ,P_LEVELS                  ! IN NO. of model levels.
     *       ,Q_LEVELS                  ! IN NO. of model wet levels.
     *       ,N_OCTAS                   ! IN NO. of cloud cover limits
     *       ,CONV_BASE_LEV(P_FIELD)    ! IN level number conv base
     *       ,CONV_TOP_LEV(P_FIELD)     ! IN level number conv top
     *       ,FIRST_POINT,LAST_POINT    ! IN 1st & last pts for calc
      REAL
     *        TEMP(P_FIELD,P_LEVELS)         ! IN temp on model levels
     *       ,Q(P_FIELD,Q_LEVELS)            ! IN spec humidity array
     *       ,P_STAR(P_FIELD)                ! IN surface press. array
     *       ,P_EXNER_HALF(P_FIELD,P_LEVELS+1)! IN 1/2 lev exner press
     *       ,OROG(P_FIELD)                  ! IN model orography array
     *       ,CONV_CLD_COVER(P_FIELD)        ! IN conv cloud cover arr
     *       ,CLOUD_COVER(P_FIELD,Q_LEVELS)  ! IN cloud cover -mod levs
     *       ,AKH(P_LEVELS+1)                ! IN A 1/2 lev hybrid coord
     *       ,BKH(P_LEVELS+1)                ! IN B 1/2 lev hybrid coord
     *       ,AK(P_LEVELS)                   ! IN A lev hybrid coord
     *       ,BK(P_LEVELS)                   ! IN B lev hybrid coord
     *       ,OCTAS(N_OCTAS)                 ! IN cloud cover limits
      LOGICAL
     *        CLD_COVER_RQD      ! IN TRUE if cloud cover data required
     *       ,LOW_CLD_RQD        ! IN TRUE if low cloud data required
     *       ,WBFL_RQD           ! IN TRUE if wet bulb freezing lev rqd
     *       ,WBTEMP_RQD         ! IN TRUE if wet bulb temp required
     *       ,CLD_TOP_RQD        ! IN TRUE if cloud top height required
C*--------------------------------------------------------------------
C Output variables----------------------------------------------------
C---------------------------------------------------------------------
      REAL
     *        CLOUD_BASE(P_FIELD,N_OCTAS) ! OUT cloud bases for amnts.
     *       ,LOW_C_FRAC(P_FIELD)         ! OUT cloud amt below 1000 ft.
     *       ,LOW_C_BASE(P_FIELD)         ! OUT base of lowest cloud.
     *       ,LOW_C_TOP(P_FIELD)          ! OUT top of lowest cloud.
     *       ,TW(P_FIELD,Q_LEVELS)        ! OUT Wet bulb temp.
     *       ,WBFLH(P_FIELD)              ! OUT Wet bulb freezing lev ht
     *       ,CLOUD_TOP(P_FIELD)          ! Cloud top height.
C*--------------------------------------------------------------------
C External subroutines called-----------------------------------------
C---------------------------------------------------------------------
      EXTERNAL V_INT_ZH, TWBULB
C*--------------------------------------------------------------------
C Local varables:-----------------------------------------------------
C---------------------------------------------------------------------
      INTEGER
     ;       I                               ! LOOP p_fields
     ;      ,J                               ! LOOP p_levels - q_levels
     ;      ,L                               ! LOOP p_levels - q_levels
     ;      ,K                               ! LOOP n_octas
     ;      ,CLOUD_BASE_LEV(P_FIELD,N_OCTAS) ! level num of cloud base
     ;      ,CLOUD_TOP_LEV(P_FIELD)          ! Level of total cld top
      REAL
     ;       PHI_STAR(P_FIELD)               ! geopotential
     ;      ,CONV_AMNT(P_FIELD)              ! conv cloud amnt octas
     ;      ,CLOUD_AMNT(P_FIELD,Q_LEVELS)    ! cloud amnt octas
     ;      ,THETA(P_FIELD,P_LEVELS)         ! pot. temp model levels
     ;      ,MODEL_HALF_HT(P_FIELD,P_LEVELS+1) ! hts of model half levs
     ;      ,P_EXNER(P_FIELD,P_LEVELS)       ! Level exner press.
     ;      ,HEIGHT(P_FIELD,P_LEVELS)        ! Level heights ASL.
     ;      ,Z                               ! Level heights.
     ;      ,PU,PL                      ! Upper & lower 1/2 lev pressure
     ;      ,M_TO_KFT                        ! convert metres to kiloft
     ;      ,PT                              ! p thickness accumulator
     ;      ,FT                              ! cloud fract accumulator
     ;      ,DP                              ! Layer pressure thickness
     ;      ,H_ASL                           ! Layer base height asl
     ;      ,H_ASLN                          ! Layer top height asl
     ;      ,FR                              ! Layer fraction below ceil
     ;      ,CP_OVER_G                       ! Used in level hghts calc.
     ;      ,THRESH                          ! Threshold val for cld top
     ;      ,CLD                             ! Intermediate cloud variab
     ;      ,FRAC                            ! Interpolation fraction
     ;      ,STR_CEILM                       ! STR_CEIL in M
      PARAMETER ( CP_OVER_G = CP / G)
      PARAMETER ( THRESH = 0.0627)
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

C----------------------------------------------------------------------C
C     Set metres to KiloFT conversion                                  C
C----------------------------------------------------------------------C
      M_TO_KFT = (1./FT2M)*0.001
      STR_CEILM=STR_CEIL * FT2M
C----------------------------------------------------------------------C
C            FIND HEIGHTS of MODEL HALF LEVELS                         C
C----------------------------------------------------------------------C
C     GEOPOTENTIAL
      DO I=1,P_FIELD
C       PHI_STAR(I) = OROG(I) * G ! for ht above sea level
        PHI_STAR(I) = 0.0         ! for ht above model orography
      ENDDO
C     TEMP to THETA
      DO J=1,P_LEVELS
        DO I=1,P_FIELD
          PU = AKH(J+1)+BKH(J+1)*P_STAR(I)
          PL = AKH(J)+BKH(J)*P_STAR(I)
          P_EXNER(I,J) = P_EXNER_C( P_EXNER_HALF(I,J+1),
     &                 P_EXNER_HALF(I,J),PU,PL,KAPPA )
          THETA(I,J) = TEMP(I,J)/P_EXNER(I,J)
        ENDDO
      ENDDO
      CALL V_INT_ZH(P_EXNER_HALF,THETA,Q,PHI_STAR,MODEL_HALF_HT,
     *  P_FIELD,P_LEVELS,Q_LEVELS)
C----------------------------------------------------------------------C
C            FIND HEIGHTS OF MODEL LEVELS                              C
C----------------------------------------------------------------------C
      IF (WBFL_RQD .OR. CLD_TOP_RQD) THEN
        DO J=1,Q_LEVELS
          DO I=1,P_FIELD
            Z = MODEL_HALF_HT(I,J) + CP_OVER_G*
     &     (1.0+C_VIRTUAL*Q(I,J))*THETA(I,J)
     &    *(P_EXNER_HALF(I,J) - P_EXNER(I,J))
            HEIGHT(I,J) = Z + OROG(I)
          ENDDO
        ENDDO
        IF(P_LEVELS .GT. Q_LEVELS) THEN
          IF(P_LEVELS .GT. Q_LEVELS) THEN
            DO J=Q_LEVELS+1,P_LEVELS
              DO I=1,P_FIELD
                Z = MODEL_HALF_HT(I,J) + CP_OVER_G*
     &            THETA(I,J)
     &         *(P_EXNER_HALF(I,J) - P_EXNER(I,J))
               HEIGHT(I,J) = Z + OROG(I)
              ENDDO
            ENDDO
          END IF
        ENDIF
      ENDIF
C----------------------------------------------------------------------C
C  CALCULATE THE WET BULB TEMP AND/OR WET BULB FREEZING LEVEL          C
C----------------------------------------------------------------------C
      IF (WBTEMP_RQD .OR. WBFL_RQD) THEN
        CALL TWBULB(Q,P_STAR,TEMP,AK,BK,P_FIELD,P_LEVELS,Q_LEVELS,
     *   TW,FIRST_POINT,LAST_POINT)
        IF (WBFL_RQD) THEN
          DO I = 1,P_FIELD
            DO L = 1,Q_LEVELS
              IF (TW(I,L) .NE. RMDI) THEN
                IF (TW(I,L) .LE. ZERODEGC) THEN
                  IF(L .EQ. 1) THEN
                    WBFLH(I) = HEIGHT(I,L)
                  ELSE
                    FRAC = (ZERODEGC - TW(I,L-1))/(TW(I,L) - TW(I,L-1))
                    WBFLH(I) = HEIGHT(I,L)*FRAC + HEIGHT(I,L-1)
     &              *(1.0-FRAC)
                  ENDIF
                  GOTO 100
                ENDIF
              ELSE
                WBFLH(I) = RMDI
              ENDIF
            ENDDO
 100      CONTINUE
          ENDDO
        ENDIF
      ENDIF
      IF (CLD_TOP_RQD) THEN
C----------------------------------------------------------------------C
C                  INITIALISE OUTPUT ARRAYS                            C
C----------------------------------------------------------------------C
        DO I=1,P_FIELD
          CLOUD_TOP(I)=RMDI
          CLOUD_TOP_LEV(I)=IMDI
        ENDDO
C----------------------------------------------------------------------C
C    CALCULATE TOTAL CLOUD TOP LEVELS                                  C
C----------------------------------------------------------------------C
        DO I=1,P_FIELD
          DO J=Q_LEVELS,1,-1
            IF (J .GT. CONV_TOP_LEV(I) .OR. J .LT. CONV_BASE_LEV(I))
     *      THEN
              IF (CLOUD_COVER(I,J) .GE. THRESH) THEN
                CLOUD_TOP_LEV(I) = J
              ENDIF
            ELSE
              CLD = CONV_CLD_COVER(I) +
     *        (1 - CONV_CLD_COVER(I))*CLOUD_COVER(I,J)
              IF (CLD .GE. THRESH) THEN
                CLOUD_TOP_LEV(I) = J
              ENDIF
            ENDIF
            IF (CLOUD_TOP_LEV(I) .GT. 0) GOTO 1000
          ENDDO
 1000     CONTINUE
        ENDDO
C----------------------------------------------------------------------C
C    CALCULATE TOTAL CLOUD TOP HEIGHTS                                 C
C----------------------------------------------------------------------C
        DO I=1,P_FIELD
          IF (CLOUD_TOP_LEV(I) .GT. 0) THEN
            CLOUD_TOP(I) = HEIGHT(I,CLOUD_TOP_LEV(I)) * M_TO_KFT
          ENDIF
        ENDDO
      ENDIF
      IF(CLD_COVER_RQD) THEN
C----------------------------------------------------------------------C
C            INITIALISE OUTPUT ARRAY                                   C
C----------------------------------------------------------------------C
      DO J=1,N_OCTAS
        DO I=1,P_FIELD
          CLOUD_BASE(I,J) = RMDI
          CLOUD_BASE_LEV(I,J) = IMDI
        ENDDO
      ENDDO
C----------------------------------------------------------------------C
C            CONVERT CONV.CLOUD COVER TO OKTAS                         C
C----------------------------------------------------------------------C
      DO I=1,P_FIELD
        CONV_AMNT(I) = CONV_CLD_COVER(I) * 8.0
      ENDDO
C----------------------------------------------------------------------C
C            CONVERT CLOUD COVER TO OKTAS                              C
C----------------------------------------------------------------------C
      DO J=1,Q_LEVELS
        DO I=1,P_FIELD
          CLOUD_AMNT(I,J) = CLOUD_COVER(I,J) * 8.0
        ENDDO
      ENDDO
C----------------------------------------------------------------------C
C            SET CLOUD BASE TO MODEL LEVELS FOR CLOUD BANDS            C
C----------------------------------------------------------------------C
      DO K=1,N_OCTAS
        DO J=1,Q_LEVELS
          DO I=1,P_FIELD
            IF(CLOUD_AMNT(I,J).GT.OCTAS(K)) THEN
              IF(CLOUD_BASE_LEV(I,K).LT.0) THEN
                CLOUD_BASE_LEV(I,K) = J
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C----------------------------------------------------------------------C
C            COMPARE WITH CONVECTIVE CLOUD AND MODIFY IF NEEDED        C
C----------------------------------------------------------------------C
      DO K=1,N_OCTAS
        DO I=1,P_FIELD
          IF(CONV_AMNT(I).GT.OCTAS(K)) THEN
            IF(CLOUD_BASE_LEV(I,K).LT.0 .OR.
     *         CLOUD_BASE_LEV(I,K).GT.CONV_BASE_LEV(I)) THEN
              CLOUD_BASE_LEV(I,K) = CONV_BASE_LEV(I)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C----------------------------------------------------------------------C
C            CONVERT LEVEL NUMBERS TO HEIGHTS (M converted to Kft)     C
C----------------------------------------------------------------------C
      DO K=1,N_OCTAS
        DO I=1,P_FIELD
          IF(CLOUD_BASE_LEV(I,K).GT.0) THEN
            CLOUD_BASE(I,K) = MODEL_HALF_HT(I,CLOUD_BASE_LEV(I,K))
            CLOUD_BASE(I,K) = CLOUD_BASE(I,K) * M_TO_KFT
          ENDIF
        ENDDO
      ENDDO
      ENDIF
      IF(LOW_CLD_RQD) THEN
C----------------------------------------------------------------------C
C            FIND CLOUD FRACTION IN AIR < STR_CEIL, LOW_C_BASE         C
C            AND LOW_C_TOP                                             C
C----------------------------------------------------------------------C
        DO I=1,P_FIELD
          PT=0.0
          FT=0.0
          LOW_C_BASE(I)=RMDI ! Initialise output variables
          LOW_C_TOP(I)=RMDI
          LOW_C_FRAC(I)=RMDI
          H_ASLN=OROG(I)
          DO J=1,Q_LEVELS
            H_ASL=H_ASLN
            H_ASLN=MODEL_HALF_HT(I,J+1)+OROG(I)
C
C     Check if have not already found low cloud base.
            IF(LOW_C_BASE(I).EQ.RMDI) THEN
C     If not, and cloud cover in this layer > threshold
              IF(CLOUD_COVER(I,J).GE.CLOUD_THRESHOLD) THEN
C     Then call the bottom of this layer the base of low cloud.
                LOW_C_BASE(I) = H_ASL / FT2M
              ENDIF
            ENDIF
C
C     Check if already found low cloud base but not top.
            IF(LOW_C_BASE(I).NE.RMDI.AND.LOW_C_TOP(I).EQ.RMDI)
     +      THEN
C     If not, and cloud cover in this layer < threshold
              IF(CLOUD_COVER(I,J).LT.CLOUD_THRESHOLD) THEN
C     Then call the bottom of this layer the top of low cloud.
                LOW_C_TOP(I) = H_ASL / FT2M
              ENDIF
            ENDIF
C
C     If bottom of layer is below low cloud ceiling (1000ft)
            IF(H_ASL.LT.STR_CEILM) THEN
C     Calculate top and bottom layer pressures
              PU = AKH(J+1)+BKH(J+1)*P_STAR(I)
              PL = AKH(J)+BKH(J)*P_STAR(I)
C     And accumulate pressure thickness and pressure weighted cloud amt
              DP = PU - PL
C     If whole layer below ceiling, simply accumulate whole layer.
              IF(H_ASLN.LT.STR_CEILM) THEN
                PT = PT + DP
                FT = FT + DP * CLOUD_COVER(I,J)
              ELSE
C     Otherwise height interpolate.
                FR = (STR_CEILM - H_ASL) / (H_ASLN - H_ASL)
                PT = PT + DP * FR
                FT = FT + DP * CLOUD_COVER(I,J) * FR
C     And set result
                LOW_C_FRAC(I) = FT / PT
              ENDIF
            ENDIF
          ENDDO ! J over Q_LEVELS
        ENDDO ! I over P_FIELD
      ENDIF
C----------------------------------------------------------------------C
C                                                                      C
C----------------------------------------------------------------------C
      RETURN
      END
