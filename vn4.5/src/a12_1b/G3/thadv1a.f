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
CLL   SUBROUTINE TH_ADV -------------------------------------------
CLL
CLL   PURPOSE:  CALCULATES MASS-WEIGHTED INCREMENTS TO THETAL
CLL DUE TO ADVECTION  BY USING EQUATION (35) TO CALCULATE PROVISIONAL
CLL VALUES OF THETAL AT THE NEW TIME-LEVEL, AND THEN RECALCULATING THE
CLL ADVECTION TERMS ON THE RIGHT-HAND SIDE OF (35) USING THESE
CLL PROVISIONAL VALUES. THE FINAL INCREMENTS ARE CALCULATED AS IN
CLL EQUATION (40). THOSE REQUIRING FILTERING ARE FILTERED AND ALL THE
CLL INCREMENTS ARE ADDED ONTO THE FIELDS USING (40).  IF RUNNING A
CLL GLOBAL MODEL POLAR IS CALLED TO UPDATE POLAR VALUES.
CLL
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL   VERSION FOR CRAY Y-MP
CLL
CLL MM, DR      <- PROGRAMMER OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL   3.4    06/08/94 Code restructured to improve parallel efficiency
CLL                   on C90.
CLL                   Authors: A. Dickinson, D. Salmond
CLL                   Reviewer: M. Mawson
CLL
CLL   3.4   23/06/94  Argument LWHITBROM added and passed to ADV_P_GD
CLL                                                      S.J.Swarbrick
!     3.5    28/03/95 MPP code: Change updateable area and
!                     add boundary swaps.  P.Burton
CLL
!     4.0   25/4/95   Corrected index for calculation of South Pole
!                     values of THETAL_INCREMENT    P.Burton
CLL   4.0   14/02/95  Option to run with half_timestep at top level
CLL                   removed.  Author: T Davies,  Reviewer: M Mawson
!     4.1    29/04/96 Remove MPP code (new QTADV1C version for MPP)
!                     and add TYPFLDPT arguments       P.Burton
!LL   4.2    16/08/96  Add TYPFLDPT arguments to FILTER subroutine
!LL                    and make the FILTER_WAVE_NUMBER arrays
!LL                    globally sized                       P.Burton
!LL   4.2    30/10/96  Move declaration of TYPFLDPT variables to top of
!LL                    declarations.  P.Burton
!LL 4.3      24/04/97 Fix to 4th order calculations -
!LL                   Calculation of NUY via ISMIN   P.Burton
!LL  4.5  05/05/98  Recode -DEF,CRAY loops to find minimum of NUX/NUY
!LL                 to vectorize on Fujitsu VPP700. RBarnes@ecmwf.int
!LL  4.5  22/06/98  Fujitsu vectorization directives. R.Barnes.
!LL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD B.
CLL
CLL   SYSTEM COMPONENTS COVERED: P121
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION:       THE EQUATIONS USED ARE (35) AND (40)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER NO. 10
CLL                        M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLL
CLLEND-------------------------------------------------------------

C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE TH_ADV
     1              (THETAL,PSTAR_OLD,PSTAR,U_MEAN,V_MEAN,
     2              SEC_P_LATITUDE,ETADOT_MEAN,RS,DELTA_AK,DELTA_BK,
     3              LATITUDE_STEP_INVERSE,ADVECTION_TIMESTEP,NU_BASIC,
     4              LONGITUDE_STEP_INVERSE,NORTHERN_FILTERED_P_ROW,
     5              SOUTHERN_FILTERED_P_ROW,P_LEVELS,
     6              U_FIELD,P_FIELD,ROW_LENGTH,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     7              TRIGS,IFAX,FILTER_WAVE_NUMBER_P_ROWS,SEC_U_LATITUDE,
     8              AKH,BKH,QCL,QCF,P_EXNER,OMEGA,
     9              Q_LEVELS,AK,BK,L_SECOND,
     &              LWHITBROM)

      IMPLICIT NONE

! All TYPFLDPT arguments are intent IN
! Comdeck TYPFLDPT
! Variables which point to useful positions in a horizontal field

      INTEGER
     &  FIRST_ROW        ! First updatable row on field
     &, TOP_ROW_START    ! First point of north-pole (global) or
!                        ! Northern (LAM) row
     &, P_LAST_ROW       ! Last updatable row on pressure point field
     &, U_LAST_ROW       ! Last updatable row on wind point field
     &, P_BOT_ROW_START  ! First point of south-pole (global) or
!                        ! Southern (LAM) row on press-point field
     &, U_BOT_ROW_START  ! First point of south-pole (global) or
!                        ! Southern (LAM) row on wind-point field
     &, upd_P_ROWS       ! number of P_ROWS to be updated
     &, upd_U_ROWS       ! number of U_ROWS to be updated
     &, FIRST_FLD_PT     ! First point on field
     &, LAST_P_FLD_PT    ! Last point on pressure point field
     &, LAST_U_FLD_PT    ! Last point on wind point field
     &, FIRST_VALID_PT   ! first valid point of data on field
     &, LAST_P_VALID_PT  ! last valid point of data on field
     &, LAST_U_VALID_PT  ! last valid point of data on field
     &, VALID_P_ROWS     ! number of valid rows of P data
     &, VALID_U_ROWS     ! number of valid rows of U data
     &, START_POINT_NO_HALO
!                        ! first non-polar point of field (misses
!                        ! halo for MPP code)
     &, START_POINT_INC_HALO
!                        ! first non-polar point of field (includes
!                        ! halo for MPP code)
     &, END_P_POINT_NO_HALO
!                        ! last non-polar point of P field (misses
!                        ! halo for MPP code)
     &, END_P_POINT_INC_HALO
!                        ! last non-polar point of P field (includes
!                        ! halo for MPP code)
     &, END_U_POINT_NO_HALO
!                        ! last non-polar point of U field (misses
!                        ! halo for MPP code)
     &, END_U_POINT_INC_HALO
!                        ! last non-polar point of U field (includes
!                        ! halo for MPP code)
     &, FIRST_ROW_PT     ! first data point along a row
     &, LAST_ROW_PT      ! last data point along a row
     &, tot_P_ROWS         ! total number of P_ROWS on grid
     &, tot_U_ROWS         ! total number of U_ROWS on grid
     &, GLOBAL_ROW_LENGTH  ! length of a global row
     &, GLOBAL_P_FIELD     ! size of a global P field
     &, GLOBAL_U_FIELD     ! size of a global U field
!


! End of comdeck TYPFLDPT

      INTEGER
     *  P_FIELD            !IN DIMENSION OF FIELDS ON PRESSSURE GRID.
     *, U_FIELD            !IN DIMENSION OF FIELDS ON VELOCITY GRID
     *, P_LEVELS           !IN NUMBER OF PRESSURE LEVELS.
     *, Q_LEVELS           !IN NUMBER OF MOIST LEVELS.
     *, ROW_LENGTH         !IN NUMBER OF POINTS PER ROW
     *, NORTHERN_FILTERED_P_ROW !IN ROW ON WHICH FILTERING STOPS
     *, SOUTHERN_FILTERED_P_ROW !IN ROW ON WHICH FILTERING STARTS AGAIN.
     &, FILTER_WAVE_NUMBER_P_ROWS(GLOBAL_P_FIELD/GLOBAL_ROW_LENGTH)
     &       ! LAST WAVE NUMBER NOT TO BE CHOPPED
     *, IFAX(10)           !IN HOLDS FACTORS OF ROW_LENGTH USED BY
     *                     ! FILTERING.


C LOGICAL VARIABLE
      LOGICAL
     *  L_SECOND     ! SET TO TRUE IF NU_BASIC IS ZERO.
     & ,LWHITBROM

      REAL
     * THETAL(P_FIELD,P_LEVELS)  !INOUT THETAL FIELD
     *                           ! MASS-WEIGHTED ON OUTPUT.

      REAL
     * U_MEAN(U_FIELD,P_LEVELS) !IN AVERAGED MASS-WEIGHTED U VELOCITY
     *                          !   FROM ADJUSTMENT STEP
     *,V_MEAN(U_FIELD,P_LEVELS) !IN AVERAGED MASS-WEIGHTED V VELOCITY
     *                          !   * COS(LAT) FROM ADJUSTMENT STEP
     *,ETADOT_MEAN(P_FIELD,P_LEVELS)  !IN AVERAGED MASS-WEIGHTED
     *                          !VERTICAL VELOCITY FROM ADJUSTMENT STEP
     *,PSTAR(P_FIELD)           !IN PSTAR FIELD AT NEW TIME-LEVEL
     *,PSTAR_OLD(P_FIELD)       !IN PSTAR AT PREVIOUS TIME-LEVEL
     *,RS(P_FIELD,P_LEVELS)     !IN RS FIELD
     *,QCL(P_FIELD,Q_LEVELS)    !IN. PRIMARY ARRAY FOR QCL
     *,QCF(P_FIELD,Q_LEVELS)    !IN. PRIMARY ARRAY FOR QCF
     *,OMEGA(U_FIELD,P_LEVELS)  !IN. TRUE VERTICAL VELOCITY DP/DT
     *,P_EXNER(P_FIELD,P_LEVELS+1) !IN. PRIMARY ARRAY FOR EXNER FUNCTION

      REAL
     * DELTA_AK(P_LEVELS)      !IN LAYER THICKNESS
     *,DELTA_BK(P_LEVELS)      !IN LAYER THICKNESS
     *,AK(P_LEVELS)            !IN HYBRID CO-ORDINATE AT FULL LEVELS
     *,BK(P_LEVELS)            !IN HYBRID CO-ORDINATE AT FULL LEVELS
     *,AKH(P_LEVELS+1)         !IN HYBRID CO-ORDINATE AT HALF LEVELS
     *,BKH(P_LEVELS+1)         !IN HYBRID CO-ORDINATE AT HALF LEVELS
     *,SEC_P_LATITUDE(P_FIELD) !IN 1/COS(LAT) AT P POINTS (2-D ARRAY)
     *,SEC_U_LATITUDE(U_FIELD) !IN 1/COS(LAT) AT U POINTS (2-D ARRAY)
     *,LONGITUDE_STEP_INVERSE  !IN 1/(DELTA LAMDA)
     *,LATITUDE_STEP_INVERSE   !IN 1/(DELTA PHI)
     *,ADVECTION_TIMESTEP      !IN
     *,NU_BASIC                !IN STANDARD NU TERM FOR MODEL RUN.
     *,TRIGS(ROW_LENGTH)       !IN HOLDS TRIGONOMETRIC FUNCTIONS USED
     *                         ! IN FILTERING.

C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: 24 ARE REQUIRED.
      REAL
     &    OMEGA_P(P_FIELD,P_LEVELS)    ! HOLDS OMEGA AT P POINTS.

      REAL
     * THETAL_FIRST_INC(P_FIELD,P_LEVELS) ! HOLDS THETAL INCREMENT
     *                           ! RETURNED BY FIRST CALL TO ADV_P_GD
     *,THETAL_SECOND_INC(P_FIELD)! HOLDS THETAL INCREMENT
     *                           ! RETURNED BY SECOND CALL TO ADV_P_GD
     *,THETAL_PROV(P_FIELD,P_LEVELS) ! HOLDS PROVISIONAL VALUE OF
     *                           ! THETAL

      REAL
     * THETAL_INCREMENT(P_FIELD,P_LEVELS) !HOLDS INCREMENT TO THETAL
     *,ZERO(P_FIELD)             !A FIELD OF ZEROES USED WHERE VERTICAL
     *                           !VELOCITY IS ZERO.

      REAL
     * NUX(P_FIELD,P_LEVELS) !COURANT NBR DEPENDENT NU AT P PTS USED
     *                    ! IN EAST-WEST ADVECTION.
     *,NUY(P_FIELD,P_LEVELS) !COURANT NBR DEPENDENT NU AT P PTS USED
     *                    ! IN NORTH-SOUTH ADVECTION.
      REAL
     * BRSP(P_FIELD,P_LEVELS) !MASS TERM AT LEVEL K

C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES
      INTEGER
     *  P_POINTS_UPDATE    ! NUMBER OF P POINTS TO BE UPDATED.
     *                     !  = ROWS*ROWLENGTH
     *, U_POINTS_UPDATE    ! NUMBER OF U POINTS TO BE UPDATED.
     *                     !  = (ROWS-1)*ROWLENGTH
     *, P_POINTS_REQUIRED  ! NUMBER OF P POINTS AT WHICH VALUES ARE
     *                     ! NEEDED TO UPDATE AT P_POINTS_UPDATE
     *, U_POINTS_REQUIRED  ! NUMBER OF U POINTS AT WHICH VALUES ARE
     *                     ! NEEDED TO UPDATE AT U_POINTS_UPDATE
     *, START_U_REQUIRED   ! FIRST U POINT OF VALUES REQUIRED TO UPDATE
     *                     ! AT P POINTS UPDATE.
     *, END_U_REQUIRED     ! LAST U POINT OF REQUIRED VALUES.

C REAL SCALARS
      REAL
     & SCALAR1,SCALAR2,CONST1,LC_LF,TIMESTEP
     &,PK, PK1         ! Pressure at half levels k and k1 (k1=k-1)
     &,P_EXNER_FULL    ! Exner pressure at full model level

C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     &  I,I1,J,K1,IK,K
     *, FILTER_SPACE ! HORIZONTAL DIMENSION OF SPACE NEEDED IN FILTERING
     *               ! ROUTINE.

C*L   EXTERNAL SUBROUTINE CALLS:---------------------------------------
      EXTERNAL ADV_P_GD,POLAR,UV_TO_P,FILTER
C*---------------------------------------------------------------------
CL    CALL COMDECK TO GET PHYSICAL CONSTANTS USED.

CL    COMDECK C_THADV HOLDS PHYSICAL CONSTANTS REQUIRED BY ROUTINE
CL    TH_ADV.
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

C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
CL    END OF COMDECK C_THADV.

CL    MAXIMUM VECTOR LENGTH ASSUMED IS P_FIELD.
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE INCLUDING SUBROUTINE CALLS:
CL---------------------------------------------------------------------
CL

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


CL---------------------------------------------------------------------
CL    SECTION 1.     INITIALISATION
CL---------------------------------------------------------------------
C INCLUDE LOCAL CONSTANTS FROM GENERAL CONSTANTS BLOCK

      LC_LF = LC + LF
      P_POINTS_UPDATE   = upd_P_ROWS*ROW_LENGTH
      U_POINTS_UPDATE   = upd_U_ROWS*ROW_LENGTH
      P_POINTS_REQUIRED = (upd_P_ROWS+2)*ROW_LENGTH
      U_POINTS_REQUIRED = (upd_U_ROWS+2)*ROW_LENGTH
      START_U_REQUIRED  = START_POINT_NO_HALO-ROW_LENGTH
      END_U_REQUIRED    = END_U_POINT_NO_HALO+ROW_LENGTH


C *IF -DEF,NOWHBR replaced by LWHITBROM logical
      IF (LWHITBROM) THEN
CL    CALCULATE BRSP TERM AT LEVEL K

      K=1
! Loop over entire field
      DO I=FIRST_VALID_PT,LAST_P_VALID_PT
        BRSP(I,K)=(3.*RS(I,K)+RS(I,K+1))*(RS(I,K)-RS(I,K+1))
     *                *BKH(K+1)*.25*(PSTAR(I)-PSTAR_OLD(I))
      ENDDO
      K=P_LEVELS
! Loop over entire field
      DO I=FIRST_VALID_PT,LAST_P_VALID_PT
        BRSP(I,K)=-(3.*RS(I,K)+RS(I,K-1))*(RS(I,K)-RS(I,K-1))
     *                *BKH(K)*.25*(PSTAR(I)-PSTAR_OLD(I))
      ENDDO

      DO K=2,P_LEVELS -1
! Loop over entire field
        DO I=FIRST_VALID_PT,LAST_P_VALID_PT
          BRSP(I,K)=((3.*RS(I,K)+RS(I,K+1))*(RS(I,K)-RS(I,K+1))*BKH(K+1)
     *              *.25*(PSTAR(I)-PSTAR_OLD(I)))
     *              -((3.*RS(I,K)+RS(I,K-1))*(RS(I,K)-RS(I,K-1))*BKH(K)
     *              *.25*(PSTAR(I)-PSTAR_OLD(I)))
        ENDDO

      ENDDO
      END IF
C *ENDIF


! Loop over entire field
      DO I=FIRST_VALID_PT,LAST_P_VALID_PT
        ZERO(I) = 0.
      END DO

CL LOOP OVER P_LEVELS+1.
CL    ON 1 TO P_LEVELS PROVISIONAL VALUES OF THE FIELD ARE CALCULATED.
CL    ON 2 TO P_LEVELS+1 THE FINAL INCREMENTS ARE CALCULATED AND ADDED
CL    ON. THE REASON FOR THIS LOGIC IS THAT THE PROVISIONAL VALUE AT
CL    LEVEL K+1 IS NEEDED BEFORE THE FINAL INCREMENT AT LEVEL K CAN BE
CL    CALCULATED.

      DO K=1,P_LEVELS+1

        TIMESTEP = ADVECTION_TIMESTEP
CL IF NOT AT P_LEVELS+1 THEN
        IF(K.LE.P_LEVELS) THEN

CL---------------------------------------------------------------------
CL    SECTION 2.     CALCULATE COURANT NUMBER DEPENDENT NU IF IN
CL                   FORECAST MODE.
CL                   CALCULATE
CL                   PROVISIONAL VALUES OF THETAL AT NEW TIME-LEVEL.
CL---------------------------------------------------------------------

C ---------------------------------------------------------------------
CL    SECTION 2.1    SET NU TO NU_BASIC DEPENDENT ON MAX COURANT
CL                   NUMBER.
C ---------------------------------------------------------------------
CL    IF NU_BASIC IS ZERO THEN DO NOT BOTHER TO CALCULATE NU
          IF(.NOT.L_SECOND) THEN
CL    CALCULATE COURANT NUMBER
C NOTE: RS AND TRIG TERMS WILL BE INCLUDED AFTER INTERPOLATION TO P
C       GRID.
CL    CALL UV_TO_P TO MOVE MEAN VELOCITIES ONTO P GRID

          CALL UV_TO_P(U_MEAN(START_U_REQUIRED,K),
     *                 NUX(START_POINT_NO_HALO,K),U_POINTS_REQUIRED,
     *                 P_POINTS_UPDATE,ROW_LENGTH,upd_P_ROWS+1)

          CALL UV_TO_P(V_MEAN(START_U_REQUIRED,K),
     *                 NUY(START_POINT_NO_HALO,K),U_POINTS_REQUIRED,
     *                 P_POINTS_UPDATE,ROW_LENGTH,upd_P_ROWS+1)

CL    CALCULATE NU FROM COURANT NUMBER INCLUDING TRIG AND RS TERMS.
          DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
            NUX(I,K) = NUX(I,K)*LONGITUDE_STEP_INVERSE
            NUY(I,K) = NUY(I,K)*LATITUDE_STEP_INVERSE
            SCALAR1 = TIMESTEP/(RS(I,K)*
     *                RS(I,K)*(DELTA_AK(K)+DELTA_BK(K)*PSTAR_OLD(I)))
            SCALAR2 = SEC_P_LATITUDE(I)*SCALAR1
            SCALAR1 = SCALAR1*SCALAR1
            SCALAR2 = SCALAR2*SCALAR2
            NUX(I,K) = (1. - NUX(I,K)*NUX(I,K)*SCALAR2)*NU_BASIC
            NUY(I,K) = (1. - NUY(I,K)*NUY(I,K)*SCALAR1)*NU_BASIC
          END DO

C     SET NUX EQUAL TO MINIMUM VALUE ALONG EACH ROW
          DO J = 1,upd_P_ROWS
          I1 = START_POINT_NO_HALO+(J-1)*ROW_LENGTH
          SCALAR1 = NUX(I1,K)
          DO I=I1+1,I1+ROW_LENGTH-1
            IF(NUX(I,K).LT.SCALAR1) THEN
              SCALAR1 = NUX(I,K)
            END IF
          END DO
          IF(SCALAR1.LT.0.) SCALAR1=0.
          DO I=I1,I1+ROW_LENGTH-1
            NUX(I,K) = SCALAR1
          END DO
          END DO

C     SET NUY EQUAL TO MINIMUM VALUE ALONG EACH COLUMN
          DO J = 1,ROW_LENGTH
          I1 = START_POINT_NO_HALO+ J-1
          SCALAR1 = NUY(I1,K)
          DO I=I1+ROW_LENGTH,END_P_POINT_NO_HALO,ROW_LENGTH
            IF(NUY(I,K).LT.SCALAR1) THEN
              SCALAR1 = NUY(I,K)
            END IF
          END DO
          IF(SCALAR1.LT.0.) SCALAR1=0.
            DO I=I1,END_P_POINT_NO_HALO,ROW_LENGTH
              NUY(I,K) = SCALAR1
            END DO
          END DO
          END IF


C ---------------------------------------------------------------------
CL    SECTION 2.2    CALL ADV_P_GD TO OBTAIN FIRST INCREMENT DUE TO
CL                   ADVECTION.
C ---------------------------------------------------------------------

CL    CALL ADV_P_GD FOR THETAL.
          K1=K+1

          IF(K.EQ.P_LEVELS) THEN
C PASS ANY THETAL VALUES AS THOSE APPARENTLY AT LEVEL K+1 AS ETADOT
C IS SET TO ZERO BY USING ARRAY ZERO.
            K1 = K-1

          CALL ADV_P_GD(THETAL(1,K1),THETAL(1,K),THETAL(1,K1),
     *                  U_MEAN(1,K),V_MEAN(1,K),ETADOT_MEAN(1,K),ZERO,
     *                  SEC_P_LATITUDE,
     *                  THETAL_FIRST_INC(1,K),NUX(1,K),NUY(1,K),P_FIELD,
     *                  U_FIELD,ROW_LENGTH,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &                 TIMESTEP,LATITUDE_STEP_INVERSE,
     *                  LONGITUDE_STEP_INVERSE,SEC_U_LATITUDE,
     *                  BRSP(1,K),L_SECOND,LWHITBROM)
          ELSE IF(K.EQ.1)THEN

C PASS ANY THETAL VALUES FOR LEVEL K-1 AS ETADOT AT LEVEL 1
C IS SET TO ZERO BY USING ARRAY ZERO.
          CALL ADV_P_GD(THETAL(1,K1),THETAL(1,K),THETAL(1,K1),
     *                  U_MEAN(1,K),V_MEAN(1,K),ZERO,
     *                  ETADOT_MEAN(1,K1),
     *                  SEC_P_LATITUDE,THETAL_FIRST_INC(1,K),
     *                  NUX(1,K),NUY(1,K),
     *                  P_FIELD,U_FIELD,ROW_LENGTH,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &                  TIMESTEP,
     *                  LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     *                  SEC_U_LATITUDE,BRSP(1,K),L_SECOND,LWHITBROM)
          ELSE
          CALL ADV_P_GD(THETAL(1,K-1),THETAL(1,K),THETAL(1,K1),
     *                  U_MEAN(1,K),V_MEAN(1,K),ETADOT_MEAN(1,K),
     *                  ETADOT_MEAN(1,K1),
     *                  SEC_P_LATITUDE,THETAL_FIRST_INC(1,K),
     *                  NUX(1,K),NUY(1,K),
     *                  P_FIELD,U_FIELD,ROW_LENGTH,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &                  TIMESTEP,
     *                  LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     *                  SEC_U_LATITUDE,BRSP(1,K),L_SECOND,LWHITBROM)

          END IF


C ---------------------------------------------------------------------
CL    SECTION 2.3    REMOVE MASS-WEIGHTING FROM INCREMENT AND ADD ONTO
CL                   FIELD TO OBTAIN PROVISIONAL VALUE.
C ---------------------------------------------------------------------

! Loop over field, missing top and bottom rows.
          DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
            SCALAR1 = RS(I,K)*RS(I,K)
     *                      *(DELTA_AK(K)+DELTA_BK(K)*PSTAR_OLD(I))
            THETAL_FIRST_INC(I,K)=THETAL_FIRST_INC(I,K)/SCALAR1
            THETAL_PROV(I,K) = THETAL(I,K)- THETAL_FIRST_INC(I,K)
          END DO

CL   GLOBAL MODEL CALCULATE PROVISIONAL POLAR VALUE.
! Fujitsu vectorization directive
!OCL NOVREC
          DO I=1,ROW_LENGTH
C NORTH POLE.
            IK = P_FIELD - ROW_LENGTH + I
            THETAL_PROV(I,K) = THETAL(I,K)
            THETAL_FIRST_INC(I,K) = -THETAL_FIRST_INC(I,K)/
     *       (RS(I,K)*RS(I,K)*(DELTA_AK(K)+DELTA_BK(K)*PSTAR_OLD(I)))
C SOUTH POLE.
            THETAL_PROV(IK,K) = THETAL(IK,K)
            THETAL_FIRST_INC(IK,K) = -THETAL_FIRST_INC(IK,K)/
     *                        (RS(IK,K)*RS(IK,K)*
     *                        (DELTA_AK(K)+DELTA_BK(K)*PSTAR_OLD(IK)))
          END DO



        END IF
CL END CONDITIONAL ON LEVEL BEING LESS THAN P_LEVELS+1
      enddo

      CALL POLAR(THETAL_PROV,THETAL_FIRST_INC,THETAL_FIRST_INC,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &           P_FIELD,P_FIELD,P_FIELD,
     &           TOP_ROW_START,P_BOT_ROW_START,
     &           ROW_LENGTH,P_LEVELS)


! Set up OMEGA_P array
! (Was SECTION 3.2):
!    SECTION 3.2    INTERPOLATE OMEGA TO P GRID AND CALCULATE
!                   REMAINING TERM IN ADVECTION EQUATION.
!                   CALCULATE TOTAL MASS-WEIGHTED INCREMENT TO FIELD.

            DO K1=1,P_LEVELS
              CALL UV_TO_P(OMEGA(START_U_REQUIRED,K1),
     &                     OMEGA_P(START_POINT_NO_HALO,K1),
     &                     U_POINTS_REQUIRED,
     &                     P_POINTS_UPDATE,ROW_LENGTH,upd_P_ROWS+1)
! Fujitsu vectorization directive
!OCL NOVREC
              DO I=1,ROW_LENGTH
                OMEGA_P(I,K1)=0.
                OMEGA_P(P_FIELD+1-I,K1)=0.
              ENDDO
            ENDDO
              CALL POLAR(OMEGA_P,OMEGA_P,OMEGA_P,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &                   P_FIELD,P_FIELD,P_FIELD,
     &                   START_POINT_NO_HALO,
     &                   END_P_POINT_NO_HALO-ROW_LENGTH+1,
     &                   ROW_LENGTH,P_LEVELS)
CL BEGIN CONDITIONAL ON LEVEL BEING GREATER THAN 1
            DO K=1,P_LEVELS+1
        IF(K.GT.1) THEN
CL---------------------------------------------------------------------
CL    SECTION 3.     ALL WORK IN THIS SECTION PERFORMED AT LEVEL-1.
CL                   CALCULATE SECOND INCREMENT DUE TO ADVECTION.
CL                   CALCULATE TOTAL INCREMENT TO FIELD AND FILTER
CL                   WHERE NECESSARY THEN UPDATE FIELD.
CL                   THE POLAR INCREMENTS ARE THEN CALCULATED AND ADDED
CL                   ON BY CALLING POLAR.
CL---------------------------------------------------------------------
CL RESET TIMESTEP APPROPRIATE TO LEVEL

         TIMESTEP = ADVECTION_TIMESTEP

         CONST1 = R/(CP*CP)*TIMESTEP
C ---------------------------------------------------------------------
CL    SECTION 3.1    CALL ADV_P_GD TO OBTAIN SECOND INCREMENT DUE TO
CL                   ADVECTION.
C ---------------------------------------------------------------------

CL    CALL ADV_P_GD FOR THETAL.
          K1=K-1
C K1 HOLDS K-1.
          IF(K.GT.P_LEVELS) THEN
C THE ZERO VERTICAL FLUX AT THE TOP IS ENSURED BY PASSING ETADOT AS
C ZERO.

          CALL ADV_P_GD(THETAL_PROV(1,K-2),THETAL_PROV(1,K-1),
     *                  THETAL_PROV(1,K-2),
     *                  U_MEAN(1,K1),V_MEAN(1,K1),ETADOT_MEAN(1,K-1),
     *                  ZERO,SEC_P_LATITUDE,
     *                  THETAL_SECOND_INC,NUX(1,K-1),NUY(1,K-1),P_FIELD,
     *                  U_FIELD,ROW_LENGTH,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &                  TIMESTEP,LATITUDE_STEP_INVERSE,
     *                  LONGITUDE_STEP_INVERSE,SEC_U_LATITUDE,
     *                  BRSP(1,K-1),L_SECOND,LWHITBROM)

          ELSE IF(K.EQ.2) THEN
C THE ZERO VERTICAL FLUX AT THE BOTTOM IS ENSURED BY PASSING ETADOT AS
C ZERO.
          CALL ADV_P_GD(THETAL_PROV(1,K),THETAL_PROV(1,K-1),
     *                  THETAL_PROV(1,K),
     *                  U_MEAN(1,K1),V_MEAN(1,K1),ZERO,
     *                  ETADOT_MEAN(1,K),
     *                 SEC_P_LATITUDE,THETAL_SECOND_INC,
     *                 NUX(1,K-1),NUY(1,K-1),
     *                  P_FIELD,U_FIELD,ROW_LENGTH,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &                  TIMESTEP,
     *                  LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     *                  SEC_U_LATITUDE,
     *                  BRSP(1,K-1),L_SECOND,LWHITBROM)
          ELSE

          CALL ADV_P_GD(THETAL_PROV(1,K-2),THETAL_PROV(1,K-1),
     *                  THETAL_PROV(1,K),
     *                  U_MEAN(1,K1),V_MEAN(1,K1),ETADOT_MEAN(1,K-1),
     *                  ETADOT_MEAN(1,K),
     *                 SEC_P_LATITUDE,THETAL_SECOND_INC,
     *                 NUX(1,K-1),NUY(1,K-1),
     *                  P_FIELD,U_FIELD,ROW_LENGTH,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &                  TIMESTEP,
     *                  LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     *                  SEC_U_LATITUDE,
     *                  BRSP(1,K-1),L_SECOND,LWHITBROM)
          END IF


C TOTAL MASS-WEIGHTED HORIZONTAL AND VERTICAL INCREMENTS ARE CALCULATED
C SEPARATELY.

          IF(K.LT.Q_LEVELS+2) THEN
! Loop over field, missing top and bottom rows
            DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO

              PK  = AKH(K)  + BKH(K) *PSTAR(I)
              PK1 = AKH(K1) + BKH(K1)*PSTAR(I)  !  K1 = K-1
              P_EXNER_FULL = P_EXNER_C
     *        (P_EXNER(I,K),P_EXNER(I,K1),PK,PK1,KAPPA)

              THETAL_INCREMENT(I,K1) = .5*(THETAL_SECOND_INC(I) +
     *                       THETAL_FIRST_INC(I,K-1)*RS(I,K1)*RS(I,K1)
     *                      *(DELTA_AK(K1)+DELTA_BK(K1)*PSTAR(I)))
     *                      -(LC*QCL(I,K1)+LC_LF*QCF(I,K1))*CONST1*
     &                       OMEGA_P(I,K1)/((AK(K1)+BK(K1)*PSTAR(I))
     *                       *(P_EXNER_FULL))

          END DO
          ELSE
! Loop over field, missing top and bottom rows
            DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
              THETAL_INCREMENT(I,K1) = .5*(THETAL_SECOND_INC(I) +
     *                       THETAL_FIRST_INC(I,K-1)*RS(I,K1)*RS(I,K1)
     *                      *(DELTA_AK(K1)+DELTA_BK(K1)*PSTAR(I)))
          END DO
          END IF

C ---------------------------------------------------------------------
CL    SECTION 3.3    IF GLOBAL MODEL CALCULATE POLAR INCREMENTS.
CL                   IF LIMITED AREA MASS-WEIGHT BOUNDARIES.
C ---------------------------------------------------------------------

CL    GLOBAL MODEL CALCULATE POLAR INCREMENT.
CL    CALCULATE MERIDIONAL FLUX AROUND POLES BY ADDING THE TWO
CL    INCREMENTS AND ALSO MASS-WEIGHTING POLAR FIELDS.
C NEGATIVE SIGN BEFORE FIRST INCS IS DUE TO THEIR SIGN HAVING BEEN
C CHANGED PRIOR TO THE CALCULATION OF THE INTERMEDIATE VALUE.
          IF(K.LT.Q_LEVELS+2) THEN
! Fujitsu vectorization directive
!OCL NOVREC
          DO I=1,ROW_LENGTH
C NORTH POLE OR NORTHERN BOUNDARY.
            IK = P_FIELD - ROW_LENGTH + I
            SCALAR1 = RS(I,K1)*RS(I,K1)
     *                       *(DELTA_AK(K1)+DELTA_BK(K1)*PSTAR(I))

            PK  = AKH(K)  + BKH(K) *PSTAR(I)
            PK1 = AKH(K1) + BKH(K1)*PSTAR(I)  !  K1 = K-1
            P_EXNER_FULL = P_EXNER_C
     *      (P_EXNER(I,K),P_EXNER(I,K1),PK,PK1,KAPPA)

            THETAL_INCREMENT(I,K1) = -.5*(THETAL_SECOND_INC(I)
     *                       - THETAL_FIRST_INC(I,K-1)*SCALAR1)
     *                      +(LC*QCL(I,K1)+LC_LF*QCF(I,K1))*CONST1*
     &                     OMEGA_P(I,K1)/((AK(K1)+BK(K1)*PSTAR(I))
     *                       *P_EXNER_FULL)
            THETAL(I,K1) = THETAL(I,K1)*SCALAR1
C SOUTH POLE OR SOUTHERN BOUNDARY.
            SCALAR2 = RS(IK,K1)*RS(IK,K1)
     *                      *(DELTA_AK(K1)+DELTA_BK(K1)*PSTAR(IK))

            PK  = AKH(K)  + BKH(K) *PSTAR(IK)
            PK1 = AKH(K1) + BKH(K1)*PSTAR(IK)  !  K1 = K-1
            P_EXNER_FULL = P_EXNER_C
     &      (P_EXNER(IK,K),P_EXNER(IK,K1),PK,PK1,KAPPA)

            THETAL_INCREMENT(IK,K1) = -.5*(THETAL_SECOND_INC(IK)
     *                      - THETAL_FIRST_INC(IK,K-1)*SCALAR2)
     *                      +(LC*QCL(IK,K1)+LC_LF*QCF(IK,K1))*CONST1*
     *                     OMEGA_P(IK,K1)/((AK(K1)+BK(K1)*PSTAR(IK))
     *                       *P_EXNER_FULL)
            THETAL(IK,K1) = THETAL(IK,K1)*SCALAR2
          END DO
          ELSE
! Fujitsu vectorization directive
!OCL NOVREC
          DO I=1,ROW_LENGTH
C NORTH POLE OR NORTHERN BOUNDARY.
            IK = P_FIELD - ROW_LENGTH + I
            SCALAR1 = RS(I,K1)*RS(I,K1)
     *                       *(DELTA_AK(K1)+DELTA_BK(K1)*PSTAR(I))
            THETAL_INCREMENT(I,K1) = -.5*(THETAL_SECOND_INC(I)
     *                       - THETAL_FIRST_INC(I,K-1)*SCALAR1)
            THETAL(I,K1) = THETAL(I,K1)*SCALAR1
C SOUTH POLE OR SOUTHERN BOUNDARY.
            SCALAR2 = RS(IK,K1)*RS(IK,K1)
     *                      *(DELTA_AK(K1)+DELTA_BK(K1)*PSTAR(IK))
            THETAL_INCREMENT(IK,K1) = -.5*(THETAL_SECOND_INC(IK)
     *                      - THETAL_FIRST_INC(IK,K-1)*SCALAR2)
            THETAL(IK,K1) = THETAL(IK,K1)*SCALAR2
          END DO
      ENDIF
CL END CONDITIONAL LEVEL GREATER THAN ONE
        END IF

CL END LOOP OVER P_LEVELS+1
      enddo

CL---------------------------------------------------------------------
CL    SECTION 4      IF GLOBAL MODEL THEN FILTER INCREMENTS AND
CL                   UPDATE POLAR VALUES BY CALLING POLAR.
CL                   UPDATE ALL OTHER VALUES.
CL---------------------------------------------------------------------


C ---------------------------------------------------------------------
CL    SECTION 4.1    CALL FILTER TO DO FILTERING.
C ---------------------------------------------------------------------

C SET FILTER_SPACE WHICH IS ROW_LENGTH+2 TIMES THE NUMBER OF ROWS TO
C BE FILTERED.

      FILTER_SPACE = (ROW_LENGTH+2)*(NORTHERN_FILTERED_P_ROW-1+
     *                P_FIELD/ROW_LENGTH-SOUTHERN_FILTERED_P_ROW)
CL    CALL FILTER FOR THETAL INCREMENTS

      CALL FILTER(THETAL_INCREMENT,P_FIELD,P_LEVELS,
     &            FILTER_SPACE,ROW_LENGTH,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &            FILTER_WAVE_NUMBER_P_ROWS,TRIGS,IFAX,
     *            NORTHERN_FILTERED_P_ROW,SOUTHERN_FILTERED_P_ROW)

C ---------------------------------------------------------------------
CL    SECTION 4.2    CALL POLAR TO UPDATE POLAR VALUES
C ---------------------------------------------------------------------

      CALL POLAR(THETAL,THETAL_INCREMENT,THETAL_INCREMENT,
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &           P_FIELD,P_FIELD,P_FIELD,
     &           TOP_ROW_START,P_BOT_ROW_START,
     &           ROW_LENGTH,P_LEVELS)

C ---------------------------------------------------------------------
CL    SECTION 4.3    UPDATE ALL OTHER POINTS.
C   OUTPUT IS MASS-WEIGHTED.
C   INCREMENTS ARE ALREADY MASS-WEIGHTED
C ---------------------------------------------------------------------

      DO K=1,P_LEVELS
C UPDATE THETAL.
CFPP$ SELECT(CONCUR)
! Loop over field, missing top and bottom rows
        DO I= START_POINT_NO_HALO,END_P_POINT_NO_HALO
        THETAL(I,K)=THETAL(I,K)*RS(I,K)*RS(I,K)*
     &        (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))-THETAL_INCREMENT(I,K)
        END DO
      END DO

CL    END OF ROUTINE TH_ADV

      RETURN
      END
