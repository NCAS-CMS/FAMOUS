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
CLL   SUBROUTINE ADV_CTL --------------------------------------------
CLL
CLL   PURPOSE:   CALCULATES THE RIGHT-HAND SIDES OF EQUATIONS (40) TO
CLL              (42) REPRESENTING THE MASS WEIGHTED FIELDS AFTER
CLL              ADVECTION AND THE ADDITION OF THE CORIOLIS TERM DUE
CLL              TO VERTICAL MOTION. THE SPATIAL DIFFERENCING SCHEME
CLL              (35) TO (38) IS USED. ONE MORE PRESSURE ROW THAN
CLL              VELOCITY ROW IS UPDATED. DIVERGENCE DAMPS VELOCITY
CLL              FIELDS AS DESCRIBED IN SECTION 3.4 OF DOCUMENTATION
CLL              PAPER NO. 10
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL   VERSION FOR CRAY Y-MP
CLL
CLL   WRITTEN BY M.H MAWSON.
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.1     24/02/93  Tidy code to remove QA Fortran messages.
CLL   3.4     22/06/94  Argument LLINTS added and passed to UV_ADV
CLL                     Argument LWHITBROM added and passed to UV_ADV,
CLL                                                TH_ADV, QT_ADV
CLL                                                 S.J.Swarbrick
CLL                     Argument X_FIELD passed to UV_ADV to reduce
CLL                     memory use of new macrotasking. R.Rawlins
CLL  4.0   1/4/95 TRACER ADVECTION OF THETAL AND QT INCLUDED AS AN
CLL               OPTION UNDER THE CONTROL OF LOGICAL L_TRACER_THETAL_QT
CLL               L_TRACER_THETAL_QT IF SET TO TRUE.
CLL               CALLS TO TH_ADV AND QT_ADV ARE REPLACED BY
CLL               CALLS TO TRAC_ADV AND TRAC_VERT_ADV
CLL               L_HALF_TIMESTEP_TOP REPLACED BY L_TRACER_THETAL_QT.
CLL               AUTHOR: T. DAVIES,  REVIEWER: M. MAWSON
!     3.5    28/03/95 MPP code: Modify P_TO_UV calls and
!                     add halo updates          P.Burton
!     4.1    22/04/96 Added TYPFLDPT arguments to dynamics routines
!                     which allows many of the differences between
!                     MPP and "normal" code to be at top level
!                     P.Burton
!  4.2  20/08/96  MPP mods for tracer advection.  RTHBarnes.
!LL   4.2    16/08/96  Make the FILTER_WAVE_NUMBER arrays globally
!LL                    sized                               P.Burton
!LL  4.2  25/11/96  Corrections to allow LAM to run in MPP mode.
!LL                                                   RTHBarnes.
!LL  4.3  17/03/97  Make initialisation of OMEGA_P safe for MPP. RTHB.
C     vn4.3    Mar. 97   T3E migration : optimisation changes
C                                       D.Salmond
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD B.
CLL
CLL   LOGIACL COMPONENTS COVERED: P12
CLL
CLL   PROJECT TASK: P1
CLL
CLL   DOCUMENTATION:        THE EQUATIONS USED ARE (35) TO (46) AND
CLL                         SECTION 3.4 IN UNIFIED MODEL DOCUMENTATION
CLL                         NO. 10  M.J.P. CULLEN, T.DAVIES AND
CLL                         M.H. MAWSON VERSION 17, DATED 11/02/91.
CLLEND-------------------------------------------------------------
C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE ADV_CTL
     1              (THETAL,QT,PSTAR_OLD,PSTAR,U_MEAN,V_MEAN,U,V,
     &              COS_U_LATITUDE,COS_P_LATITUDE,
     2              SEC_P_LATITUDE,ETADOT_MEAN,RS,DELTA_AK,DELTA_BK,
     3              LATITUDE_STEP_INVERSE,ADVECTION_TIMESTEP,NU_BASIC,
     4              LONGITUDE_STEP_INVERSE,NORTHERN_FILTERED_P_ROW,
     5              SOUTHERN_FILTERED_P_ROW,Q_LEVELS,
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
     7              P_LEVELS,SEC_U_LATITUDE,F1,F2,AK,BK,KD,AKH,BKH,
     8              COS_U_LONGITUDE,SIN_U_LONGITUDE,TRIGS,IFAX,
     9              FILTER_WAVE_NUMBER_P_ROWS,OMEGA,QCL,QCF,P_EXNER,
     &              LLINTS,LWHITBROM,
     &              L_TRACER_THETAL_QT,NSWEEP,L_SUPERBEE)

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
     *, P_LEVELS           !IN    NUMBER OF PRESSURE LEVELS.
     *, Q_LEVELS           !IN    NUMBER OF MOIST LEVELS.
     *, ROW_LENGTH         !IN    NUMBER OF POINTS PER ROW
     *, NORTHERN_FILTERED_P_ROW !IN ROW ON WHICH FILTERING STOPS
     *                     ! MOVING TOWARDS THE EQUATOR.
     *, SOUTHERN_FILTERED_P_ROW !IN ROW ON WHICH FILTERING STARTS AGAIN.
     *                     ! MOVING TOWARDS SOUTHPOLE.
     *, IFAX(10)           !IN HOLDS FACTORS OF ROW_LENGTH USED BY
     *                     ! FILTERING.
     *, NSWEEP(P_FIELD/ROW_LENGTH,P_LEVELS) !IN
     *                  ! NUMBER OF EAST_WEST TIMESTEPS NEEDED FOR
     *                  ! EACH LATITUDE WHEN USING TRACER ADVECTION.
     *, FIRST_POINT     !
     *, POINTS          !

      INTEGER
     &  FILTER_WAVE_NUMBER_P_ROWS(GLOBAL_P_FIELD/GLOBAL_ROW_LENGTH)
!       LAST WAVE NUMBER NOT TO BE CHOPPED ON A P ROW
      LOGICAL
     &  L_SUPERBEE             ! FORM OF LIMITER USED IN TRACER
     &                         ! ADVECTION
     & ,L_TRACER_THETAL_QT     ! LOGICAL TRUE IF USING TRACER
     &                         ! ADVECTION FOR THETAL & QT
       INTEGER
     &  P_POINTS_UPDATE
     & ,START_U_REQUIRED
     & ,P_POINTS_REQUIRED
     & ,U_POINTS_REQUIRED

     & ,LLINTS              !Logical switch for linear TS
     & ,LWHITBROM           !Log swch for White & Bromley terms

      REAL
     * U_MEAN(U_FIELD,P_LEVELS)  !IN AVERAGED MASS-WEIGHTED U VELOCITY
     *                           !   FROM ADJUSTMENT STEP
     *,V_MEAN(U_FIELD,P_LEVELS)  !IN AVERAGED MASS-WEIGHTED V VELOCITY
     *                           !   * COS(LAT) FROM ADJUSTMENT STEP
     *,ETADOT_MEAN(P_FIELD,P_LEVELS)  !IN AVERAGED MASS-WEIGHTED
     *                          !VERTICAL VELOCITY FROM ADJUSTMENT STEP
     *,PSTAR(P_FIELD)            !IN PSTAR FIELD AT NEW TIME-LEVEL
     *,PSTAR_OLD(P_FIELD)        !IN PSTAR AT PREVIOUS TIME-LEVEL
     *,RS(P_FIELD,P_LEVELS)      !IN RS FIELD
     *,TRIGS(ROW_LENGTH)        !IN HOLDS TRIGONOMETRIC FUNCTIONS USED
     *                          ! IN FILTERING.
     *,QCL(P_FIELD,Q_LEVELS)    !IN. PRIMARY ARRAY FOR QCL.
     *,QCF(P_FIELD,Q_LEVELS)    !IN. PRIMARY ARRAY FOR QCF.
     *,P_EXNER(P_FIELD,P_LEVELS+1) !IN. PRIMARY ARRAY FOR P_EXNER.

      REAL
     * U(U_FIELD,P_LEVELS)       !INOUT U FIELD, MASS-WEIGHTED ON OUT.
     *,V(U_FIELD,P_LEVELS)       !INOUT V FIELD, MASS-WEIGHTED ON OUT.
     *,THETAL(P_FIELD,P_LEVELS)  !INOUT THETAL FIELD
     *,QT(P_FIELD,Q_LEVELS)      !INOUT QT FIELD.

      REAL
     * DELTA_AK(P_LEVELS)     !IN    LAYER THICKNESS
     *,DELTA_BK(P_LEVELS)     !IN    LAYER THICKNESS
     *,AK(P_LEVELS)           !IN    FIRST TERM IN HYBRID CO-ORDS.
     *,BK(P_LEVELS)           !IN    SECOND TERM IN HYBRID CO-ORDS.
     *,AKH(P_LEVELS+1)        !IN    AK AT HALF LEVELS
     *,BKH(P_LEVELS+1)        !IN    BK AT HALF LEVELS
     &,COS_P_LATITUDE(P_FIELD) !IN  COS_LAT AT P_POINTS (2D ARRAY)
     *,SEC_P_LATITUDE(P_FIELD) !IN  1/COS(LAT) AT P POINTS (2-D ARRAY)
     *,COS_U_LATITUDE(U_FIELD) !IN  COS(LAT) AT U POINTS (2-D ARRAY)
     *,SEC_U_LATITUDE(U_FIELD) !IN  1/COS(LAT) AT U POINTS (2-D ARRAY)
     *,COS_U_LONGITUDE(ROW_LENGTH) !IN COS(LONGITUDE) AT U-POINTS.
     *,SIN_U_LONGITUDE(ROW_LENGTH) !IN SIN(LONGITUDE) AT U-POINTS.
     *,LONGITUDE_STEP_INVERSE !IN 1/(DELTA LAMDA)
     *,LATITUDE_STEP_INVERSE  !IN 1/(DELTA PHI)
     *,ADVECTION_TIMESTEP     !IN
     *,NU_BASIC               !IN STANDARD NU TERM FOR MODEL RUN.
     *,F1(U_FIELD)            !IN A CORIOLIS TERM SEE DOCUMENTATION
     *,F2(U_FIELD)            !IN A CORIOLIS TERM SEE DOCUMENTATION
     *,KD(P_LEVELS)           !IN DIVERGENCE DAMPING COEFFICIENTS

      REAL
     * OMEGA(U_FIELD,P_LEVELS) !OUT TRUE VERTICAL VELOCITY
C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C  DEFINE LOCAL ARRAYS: 3 ARE REQUIRED
      REAL
     * WORK1(U_FIELD)         ! GENERAL WORKSPACE
     *,WORK2(P_FIELD)         ! GENERAL WORKSPACE
     &,   OMEGA_P(P_FIELD)    ! HOLDS OMEGA AT P POINTS.

C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES

C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     &  I,K,K1

      INTEGER X_FIELD  ! 1 IF 2ND ORDER ELSE U_FIELD IF 4TH ORDER

C   REAL SCALARS
      REAL
     &  CONST1,LC_LF,TIMESTEP
     &  ,PK, PK1       ! Pressure at half levels
     &  ,P_EXNER_FULL  ! Exner pressure at full model level

C LOGICAL VARIABLE
      LOGICAL
     * L_SECOND                ! SET TO TRUE IF NU_BASIC EQUAL TO ZERO
C*L   EXTERNAL SUBROUTINE CALLS:---------------------------------------
      EXTERNAL TH_ADV,QT_ADV,UV_ADV,P_TO_UV,DIV_DAMP
     &         ,TRAC_ADV,TRAC_VERT_ADV,UV_TO_P,POLAR
C*---------------------------------------------------------------------


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



CL    MAXIMUM VECTOR LENGTH ASSUMED IS U_FIELD.
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE INCLUDING SUBROUTINE CALLS:
CL---------------------------------------------------------------------
CL
C****************************************************************
C         INTEGERS AND VARIABLES NEEDED WHEN USING
C         TRACER ADVECTION OF THETAL & QT
C*****************************************************************
      IF(L_TRACER_THETAL_QT)THEN
       LC_LF=LC + LF
       P_POINTS_UPDATE=upd_P_ROWS*ROW_LENGTH
       START_U_REQUIRED  = START_POINT_NO_HALO-ROW_LENGTH
       P_POINTS_REQUIRED = (upd_P_ROWS+2)*ROW_LENGTH
       U_POINTS_REQUIRED = (upd_U_ROWS+2)*ROW_LENGTH
      ENDIF
CL---------------------------------------------------------------------
CL    SECTION 1.     INTERPOLATE FIELDS ONTO U GRID.
CL---------------------------------------------------------------------

      IF(NU_BASIC.EQ.0.) THEN
        L_SECOND=.TRUE.
      X_FIELD=1
      ELSE
        L_SECOND=.FALSE.
      X_FIELD=U_FIELD
      END IF

C----------------------------------------------------------------------
CL    SECTION 1.1    INTERPOLATE PSTAR ONTO U GRID.
C----------------------------------------------------------------------

      CALL P_TO_UV(PSTAR,WORK1,P_FIELD,U_FIELD,ROW_LENGTH,tot_P_ROWS)

C----------------------------------------------------------------------
CL    SECTION 1.2    INTERPOLATE PSTAR_OLD ONTO U GRID.
C----------------------------------------------------------------------

      CALL P_TO_UV(PSTAR_OLD,WORK2,P_FIELD,U_FIELD,ROW_LENGTH,
     &             tot_P_ROWS)


CL
CL---------------------------------------------------------------------
CL    SECTION 2.     CALL DIV_DAMP TO PERFORM DIVERGENCE DAMPING.
CL---------------------------------------------------------------------

C PSTAR_OLD ON U GRID IS HELD IN WORK2.

      CALL DIV_DAMP(U,V,RS,SEC_U_LATITUDE,WORK2,COS_U_LATITUDE,KD,
     *              LONGITUDE_STEP_INVERSE,LATITUDE_STEP_INVERSE,
     *              P_FIELD,U_FIELD,ROW_LENGTH,P_LEVELS,
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
     *              BKH,ADVECTION_TIMESTEP,DELTA_AK,DELTA_BK,
     *              COS_U_LONGITUDE,SIN_U_LONGITUDE,SEC_P_LATITUDE)

CL
CL---------------------------------------------------------------------
CL    SECTION 3.     CALL UV_ADV TO ADVECT U AND V.
CL---------------------------------------------------------------------

C PSTAR ON U GRID IS HELD IN WORK1.
C PSTAR_OLD ON U GRID IS HELD IN WORK2.

      CALL UV_ADV (U,V,WORK2,WORK1,U_MEAN,V_MEAN,
     *             SEC_U_LATITUDE,ETADOT_MEAN,RS,DELTA_AK,DELTA_BK,AK,
     *             BK,F1,F2,LATITUDE_STEP_INVERSE,ADVECTION_TIMESTEP,
     *             NU_BASIC,LONGITUDE_STEP_INVERSE,U_FIELD,P_FIELD,
     *             ROW_LENGTH,P_LEVELS,
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
     *             COS_U_LONGITUDE,SIN_U_LONGITUDE,SEC_P_LATITUDE,
     &             AKH,BKH,OMEGA,L_SECOND,LLINTS,
     &             LWHITBROM,X_FIELD)
CL
CL---------------------------------------------------------------------
CL    SECTION 4.     CALL TH_ADV TO ADVECT THETAL AND QT_ADV TO ADVECT
CL                   QT USING STANDARD HEUN ADVECTION.
CL                   IF USING TRACER ADVECTION FOR THETAL & QT
CL                   THEN CALL APPROPRIATE TRACER ROUTINES.
CL---------------------------------------------------------------------
      IF(.NOT.L_TRACER_THETAL_QT)THEN
CL---------------------------------------------------------------
C    SECTION 4.1  HEUN ADVVECTION SCHEME
C
CL----------------------------------------------------------------

      CALL TH_ADV (THETAL,PSTAR_OLD,PSTAR,U_MEAN,V_MEAN,SEC_P_LATITUDE,
     *            ETADOT_MEAN,RS,DELTA_AK,DELTA_BK,LATITUDE_STEP_INVERSE
     *            ,ADVECTION_TIMESTEP,NU_BASIC,LONGITUDE_STEP_INVERSE,
     *            NORTHERN_FILTERED_P_ROW,SOUTHERN_FILTERED_P_ROW,
     *            P_LEVELS,U_FIELD,P_FIELD,ROW_LENGTH,
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
     *            TRIGS,IFAX,FILTER_WAVE_NUMBER_P_ROWS,SEC_U_LATITUDE,
     *            AKH,BKH,QCL,QCF,P_EXNER,OMEGA,
     &            Q_LEVELS,AK,BK,L_SECOND,
     &            LWHITBROM)                    

      CALL QT_ADV (QT,PSTAR_OLD,PSTAR,U_MEAN,V_MEAN,SEC_P_LATITUDE,
     *            ETADOT_MEAN,RS,DELTA_AK,DELTA_BK,LATITUDE_STEP_INVERSE
     *            ,ADVECTION_TIMESTEP,NU_BASIC,LONGITUDE_STEP_INVERSE,
     *            NORTHERN_FILTERED_P_ROW,SOUTHERN_FILTERED_P_ROW,
     *            Q_LEVELS,P_LEVELS,U_FIELD,P_FIELD,ROW_LENGTH,
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
     *            TRIGS,IFAX,FILTER_WAVE_NUMBER_P_ROWS,
     &            SEC_U_LATITUDE,AKH,BKH,L_SECOND,
     &            LWHITBROM)

      ELSE
CL---------------------------------------------------------------
C    SECTION 4.2  TRACER ADVECTION OF THETAL AND QT
C
CL----------------------------------------------------------------
      DO K=1,P_LEVELS
        CALL TRAC_ADV(THETAL(1,K),NSWEEP(1,K),U_MEAN(1,K),V_MEAN(1,K),
     &                U_FIELD,P_FIELD,ADVECTION_TIMESTEP,ROW_LENGTH,
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
     &                SEC_P_LATITUDE,COS_P_LATITUDE,RS(1,K),  
     &                PSTAR_OLD,DELTA_AK(K),DELTA_BK(K),
     &                LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     &                L_SUPERBEE)
      END DO

C  Set temperature flux through lower boundary to zero
      DO I=1,P_FIELD
        WORK2(I)=0.
      END DO

      FIRST_POINT=ROW_LENGTH
      POINTS = P_FIELD - 2*ROW_LENGTH + 2

      TIMESTEP=ADVECTION_TIMESTEP
      CONST1=R/(CP*CP)*TIMESTEP
      CALL TRAC_VERT_ADV(THETAL,ETADOT_MEAN,PSTAR,P_FIELD,
     &                   TIMESTEP,1,P_LEVELS,FIRST_POINT,
     &                   POINTS,P_LEVELS,1,P_LEVELS,RS,AK,BK,DELTA_AK,
     &                   DELTA_BK,WORK2,L_TRACER_THETAL_QT,L_SUPERBEE)
C ---------------------------------------------------------------------
CL                   INTERPOLATE OMEGA TO P GRID AND CALCULATE
CL                   REMAINING TERM IN ADVECTION EQUATION.
CL                   CALCULATE TOTAL MASS-WEIGHTED INCREMENT TO FIELD.
C ---------------------------------------------------------------------

          DO 110 K=1,P_LEVELS

          CALL UV_TO_P(OMEGA(START_U_REQUIRED,K),
     &                 OMEGA_P(START_POINT_NO_HALO),U_POINTS_REQUIRED,
     &                 P_POINTS_UPDATE,ROW_LENGTH,upd_P_ROWS+1)

          DO I = FIRST_VALID_PT,FIRST_VALID_PT+ROW_LENGTH-1
            OMEGA_P(I)=0.
          END DO
          DO I = LAST_P_VALID_PT-ROW_LENGTH+1,LAST_P_VALID_PT
            OMEGA_P(I)=0.
          END DO

C SET UP POLAR VALUE OF OMEGA

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
     &               P_FIELD,P_FIELD,P_FIELD,
     &               START_POINT_NO_HALO,
     &               END_P_POINT_NO_HALO-ROW_LENGTH+1,
     &               ROW_LENGTH,1)
C TOTAL MASS-WEIGHTED HORIZONTAL AND VERTICAL INCREMENTS ARE CALCULATED
C SEPARATELY.

          IF(K.LT.Q_LEVELS+1) THEN
            DO  I = FIRST_POINT,FIRST_POINT+POINTS-1 

              PK  = AKH(K+1)+ BKH(K+1)*PSTAR(I)
              PK1 = AKH(K) + BKH(K)*PSTAR(I)
              P_EXNER_FULL = P_EXNER_C
     &        (P_EXNER(I,K+1),P_EXNER(I,K),PK,PK1,KAPPA)

              WORK2(I) =
     &                   -(LC*QCL(I,K)+LC_LF*QCF(I,K))*CONST1*
     &                    OMEGA_P(I)/((AK(K)+BK(K)*PSTAR(I))
     &                    *(P_EXNER_FULL)*
     &        RS(I,K)*RS(I,K)*(DELTA_AK(K)+DELTA_BK(K)*PSTAR(I)))
              THETAL(I,K) =THETAL(I,K)+WORK2(I)
            END DO
          END IF

CL END LOOP OVER P_LEVELS+1
 110  CONTINUE

C  Copy polar values along row
      DO K=1,P_LEVELS
        DO I=1,ROW_LENGTH-1
          THETAL(I,K) = THETAL(ROW_LENGTH,K)
          THETAL(P_FIELD+1-I,K) = THETAL(P_FIELD+1-ROW_LENGTH,K)
        END DO
      END DO
                                                                        
      DO K=1,Q_LEVELS
        CALL TRAC_ADV(QT(1,K),NSWEEP(1,K),U_MEAN(1,K),V_MEAN(1,K),
     &                U_FIELD,P_FIELD,ADVECTION_TIMESTEP,ROW_LENGTH,
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
     &                SEC_P_LATITUDE,COS_P_LATITUDE,RS(1,K), 
     &                PSTAR_OLD,DELTA_AK(K),DELTA_BK(K),
     &                LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     &                L_SUPERBEE)
      END DO

C  Set moisture flux through lower boundary to zero
      DO I=1,P_FIELD
        WORK2(I)=0.
      END DO

! Values of FIRST_POINT and POINTS
! should be unaltered from those set for Thetal

      CALL TRAC_VERT_ADV(QT,ETADOT_MEAN,PSTAR,P_FIELD,
     &                   TIMESTEP,1,Q_LEVELS,FIRST_POINT,
     &                   POINTS,P_LEVELS,1,Q_LEVELS,RS,AK,BK,DELTA_AK,
     &                   DELTA_BK,WORK2,L_TRACER_THETAL_QT,L_SUPERBEE)

C     END DO

C  Copy polar values along row
      DO K=1,Q_LEVELS
        DO I=1,ROW_LENGTH-1
          QT(I,K) = QT(ROW_LENGTH,K)
          QT(P_FIELD+1-I,K) = QT(P_FIELD+1-ROW_LENGTH,K)
        END DO
      END DO
      ENDIF ! L_TRACER_THETAL_QT  

CL END OF ROUTINE ADV_CTL

      RETURN
      END
