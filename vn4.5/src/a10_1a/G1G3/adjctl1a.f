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
CLL   SUBROUTINE ADJ_CTL ---------------------------------------------
CLL
CLL   PURPOSE:   INTEGRATES SURFACE PRESSURE, POTENTIAL TEMPERATURE,
CLL              AND HORIZONTAL WIND COMPONENTS THROUGH A SPECIFIED
CLL              NUMBER OF ADJUSTMENT STEPS. AT THE END OF THE ROUTINE
CLL              UPDATED VALUES OF ALL THESE FIELDS ALONG WITH THE
CLL              UPDATED EXNER PRESSURE ARE HELD IN THE ARGUMENTS.
CLL              FOURIER FILTERING IS PERFORMED UNDER THE
CLL              UPDATE IDENTIFIER 'GLOBAL'. ONE MORE PRESSURE ROW IS
CLL              UPDATED THAN VELOCITY ROW.
CLL              FIRST_ROW IS NORTHERNMOST PRESSURE ROW TO BE UPDATED.
CLL              FIRST_U_ROW UPDATED IS THE FIRST ONE TO THE SOUTH OF
CLL              THE FIRST P ROW.
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL   VERSION FOR CRAY Y-MP
CLL   WRITTEN BY M.H MAWSON.
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.1     24/02/93  Tidy code to remove QA Fortran messages.
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.    Author: Tracey Smith.
CLL   3.4    22/06/94 Arguments LLINTS, LWHITBROM added and passed to
CLL                                               VERT_VEL, UV_ADJ
CLL                                             S.J.Swarbrick
CLL
CLL   3.4    06/08/94 Micro tasking directives inserted to improve
CLL                   parallel efficiency on C90.
CLL                   Authors: A. Dickinson, D. Salmond
CLL                   Reviewer: M. Mawson
!     3.5    28/03/95 MPP code: Change updateable area, add halo
!                     updates.                          P.Burton
!     4.1    02/04/96 Added TYPFLDPT arguments to dynamics routines
!                     which allows many of the differences between
!                     MPP and "normal" code to be at top level
!                     Added LEVELS argument to POLAR_UV
!                     P.Burton
!LL   4.2    16/08/96  Add TYPFLDPT arguments to FILTER subroutine
!LL                    and make the FILTER_WAVE_NUMBER arrays
!LL                    globally sized                      P.Burton
!LL  4.2  25/11/96  Corrections to allow LAM to run in MPP mode.
!LL                                                   RTHBarnes.
!     4.2    Oct. 96  T3E migration: *DEF CRAY removed; HF functions
!                      replaced by T3E rtor_v funtion (*DEF T3E)
!                      code restructured appropriately.
!                                      S.J.Swarbrick
C     vn4.3    Mar. 97   T3E migration : optimisation changes
C                                       D.Salmond
!  4.5  23/10/98  Introduce Single Column Model. JC Thil
!  4.5  12/05/98  Replace **k by exp(k*log( )) for faster running
!                 on Fujitsu VPP700.  RBarnes@ecmwf.int
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL   STANDARD B, VERSION 2, DATED 18/01/90
CLL   SYSTEM COMPONENTS COVERED: P11
CLL   SYSTEM TASK: P1
CLL   DOCUMENTATION:  THE EQUATIONS USED ARE (23) TO (30)
CLL                   IN UNIFIED MODEL DOCUMENTATION PAPER NO. 10
CLL                   M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON,
CLLEND-------------------------------------------------------------

C*L   ARGUMENTS:---------------------------------------------------

      SUBROUTINE ADJ_CTL
     1  (U,V,THETA,Q,PSTAR,OROG_HEIGHT,RS,U_MEAN,V_MEAN,P_EXNER,
     2   ETADOT_MEAN,PSTAR_OLD,COS_P_LATITUDE,COS_U_LATITUDE,
     3   SEC_P_LATITUDE,SEC_U_LATITUDE,TAN_U_LATITUDE,F1,F2,F3,
     4   LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,AK,BK,DELTA_AK,
     5   DELTA_BK,THETA_REF,ADJUSTMENT_TIMESTEP,ADJUSTMENT_STEPS,
     6   NORTHERN_FILTERED_P_ROW,SOUTHERN_FILTERED_P_ROW,ROW_LENGTH,
     7   P_LEVELS,Q_LEVELS,
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
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
     7   P_FIELD,U_FIELD,AKH,BKH,
     8   AKH_TO_THE_KAPPA,BKH_TO_THE_KAPPA,AK_TO_THE_KAPPA,
     9   BK_TO_THE_KAPPA,COS_U_LONGITUDE,
     *   SIN_U_LONGITUDE,TRIGS,IFAX,FILTER_WAVE_NUMBER_P_ROWS,
     *   FILTER_WAVE_NUMBER_U_ROWS,ERROR_CODE,ERROR_MESSAGE,
     & L_NEG_PSTAR,PHI_OUT,L_PHI_OUT,ADJ_TIME_SMOOTHING_WEIGHT,
     & ADJ_TIME_SMOOTHING_COEFF,LLINTS,LWHITBROM)

      IMPLICIT NONE

      LOGICAL
     *  L_NEG_PSTAR    !IN SWITCH, IF TRUE THEN NEGATIVE PSTAR VALUES
     *                 ! WILL BE DETECTED AND OUTPUT.
     *, L_PHI_OUT      !IN. IF TRUE THEN PHI REQUIRED AS OUTPUT.
     *, LLINTS         !Logical switch for linear TS
     *, LWHITBROM      !Logical switch for White & Bromley terms

      INTEGER
     *  P_FIELD            !IN DIMENSION OF FIELDS ON PRESSSURE GRID.
     *, U_FIELD            !IN DIMENSION OF FIELDS ON VELOCITY GRID
     *, P_LEVELS           !IN NUMBER OF PRESSURE LEVELS TO BE UPDATED.
     *, Q_LEVELS           !IN NUMBER OF MOIST LEVELS TO BE UPDATED.
     *, ROW_LENGTH         !IN    NUMBER OF POINTS PER ROW
     *, ADJUSTMENT_STEPS   !IN NUMBER OF ADJUSTMENT STEPS
! All TYPFLDPT arguments are intent IN
! Comdeck TYPFLDPT
! Variables which point to useful positions in a horizontal field

      INTEGER
     &  FIRST_ROW        ! First updatable row on field
     &, TOP_ROW_START    ! First point of north-pole (global) or
!                        ! Northern (LAM) row
!                        ! for processors not at top of LPG, this
!                        ! is the first point of valid data
!                        ! (ie. Northern halo).
     &, P_LAST_ROW       ! Last updatable row on pressure point field
     &, U_LAST_ROW       ! Last updatable row on wind point field
     &, P_BOT_ROW_START  ! First point of south-pole (global) or
!                        ! Southern (LAM) row on press-point field
     &, U_BOT_ROW_START  ! First point of south-pole (global) or
!                        ! Southern (LAM) row on wind-point field
!                        ! for processors not at base of LPG, this
!                        ! is the start of the last row of valid data
!                        ! (ie. Southern halo).
     &, upd_P_ROWS       ! number of P_ROWS to be updated
     &, upd_U_ROWS       ! number of U_ROWS to be updated
     &, FIRST_FLD_PT     ! First point on field
     &, LAST_P_FLD_PT    ! Last point on pressure point field
     &, LAST_U_FLD_PT    ! Last point on wind point field
! For the last three variables, these indexes are the start points
! and end points of "local" data - ie. missing the top and bottom
! halo regions.
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
! For the last two variables, these indexes are the start and
! end points along a row of the "local" data - ie. missing out
! the east and west halos
     &, tot_P_ROWS         ! total number of P_ROWS on grid
     &, tot_U_ROWS         ! total number of U_ROWS on grid
     &, GLOBAL_ROW_LENGTH  ! length of a global row
     &, GLOBAL_P_FIELD     ! size of a global P field
     &, GLOBAL_U_FIELD     ! size of a global U field
!

     &, MY_PROC_ID         ! my processor id
     &, NP_PROC_ID         ! processor number of North Pole Processor
     &, SP_PROC_ID         ! processor number of South Pole Processor
     &, GC_ALL_GROUP       ! group id of group of all processors
     &, GC_ROW_GROUP       ! group id of group of all processors on this
!                          ! processor row
     &, GC_COL_GROUP       ! group id of group of all processors on this
!                          ! processor column
     &, N_PROCS            ! total number of processors

     &, EW_Halo            ! Halo size in the EW direction
     &, NS_Halo            ! Halo size in the NS direction

     &, halo_4th           ! halo size for 4th order calculations
     &, extra_EW_Halo      ! extra halo size required for 4th order
     &, extra_NS_Halo      ! extra halo size required for 4th order
     &, LOCAL_ROW_LENGTH   ! size of local row
     &, FIRST_GLOBAL_ROW_NUMBER
!                          ! First row number on Global Grid    

! Variables which indicate if special operations are required at the
! edges.
      LOGICAL
     &  at_top_of_LPG    ! Logical variables indicating if this
     &, at_right_of_LPG  ! processor is at the edge of the Logical
     &, at_base_of_LPG   ! Processor Grid and should process its edge
     &, at_left_of_LPG   ! data differently.

! End of comdeck TYPFLDPT

      INTEGER
     *  ERROR_CODE         !INOUT. 0 ON ENTRY. NON-ZERO ON OUT IF
     *                     ! ABNORMAL RESULT OBTAINED.

      CHARACTER*80 ERROR_MESSAGE

      INTEGER
     *  NORTHERN_FILTERED_P_ROW !IN P ROW ON WHICH FILTERING STOPS
     *                          ! MOVING TOWARDS EQUATOR
     *, SOUTHERN_FILTERED_P_ROW !IN P ROW ON WHICH FILTERING STARTS
     *                          ! AGAIN MOVING TOWARDS SOUTH POLE
     &, FILTER_WAVE_NUMBER_P_ROWS(GLOBAL_P_FIELD/GLOBAL_ROW_LENGTH)
     &               ! LAST WAVE NUMBER NOT TO BE CHOPPED ON A P ROW
     &, FILTER_WAVE_NUMBER_U_ROWS(GLOBAL_U_FIELD/GLOBAL_ROW_LENGTH)
     &               ! LAST WAVE NUMBER NOT TO BE CHOPPED ON A U ROW
     *, IFAX(10)           !IN HOLDS FACTORS OF ROW_LENGTH USED BY
     *                     ! FILTERING.
     *,ADJ_TIME_SMOOTHING_WEIGHT(ADJUSTMENT_STEPS) !IN COEFFICIENTS FOR
     *                         ! FINITE DIFFERENCE SMOOTHING DERIVATIVE

      REAL
     * U(U_FIELD,P_LEVELS)    !INOUT U FIELD
     *,V(U_FIELD,P_LEVELS)    !INOUT V FIELD
     *,THETA(P_FIELD,P_LEVELS)!INOUT THETA FIELD
     *,P_EXNER(P_FIELD,P_LEVELS+1)!INOUT EXNER PRESSURE FIELD.
     *,Q(P_FIELD,Q_LEVELS)    !INOUT Q FIELD
     *,PSTAR(P_FIELD)         !INOUT PSTAR FIELD

      REAL
     * U_MEAN(U_FIELD,P_LEVELS) !OUT HOLDS MASS-WEIGHTED U
     *                        !  AVERAGED OVER ADJUSTMENT STEPS.
     *,V_MEAN(U_FIELD,P_LEVELS) !OUT HOLDS MASS-WEIGHTED V*COS(PHI)
     *                        !  AVERAGED OVER ADJUSTMENT STEPS.
     *,ETADOT_MEAN(P_FIELD,P_LEVELS) !OUT HOLDS MASS-WEIGHTED VERTICAL
     *                        ! VELOCITY AVERAGED OVER ADJUSTMENT
     *                        ! STEPS.
     *,PSTAR_OLD(P_FIELD)     !OUT HOLDS VALUE OF PSTAR ON PREVIOUS
     *                        ! TIMESTEP
     *,RS(P_FIELD,P_LEVELS)   !OUT RS FIELD
     *,PHI_OUT(P_FIELD,P_LEVELS) !OUT. HOLDS PHI IF DIAGNOSTIC
     *                           !     REQUIRED.

      REAL
     * DELTA_AK(P_LEVELS)       !IN    LAYER THICKNESS
     *,DELTA_BK(P_LEVELS)       !IN    LAYER THICKNESS
     *,AK(P_LEVELS)             !IN    VALUE OF A AT P POINTS
     *,BK(P_LEVELS)             !IN    VALUE OF B AT P POINTS
     *,AK_TO_THE_KAPPA(P_LEVELS)!IN (A/100000)**(R/CP) AT FULL LEVELS
     *,BK_TO_THE_KAPPA(P_LEVELS)!IN (B/100000)**(R/CP) AT FULL LEVELS
     *,AKH(P_LEVELS+1)          !IN    VALUE OF A AT HALF LEVELS.
     *,BKH(P_LEVELS+1)          !IN    VALUE OF B AT HALF LEVELS.
     *,AKH_TO_THE_KAPPA(P_LEVELS+1)!IN (A/100000)**(R/CP)
     *                                     !AT HALF LEVELS
     *,BKH_TO_THE_KAPPA(P_LEVELS+1)!IN (B/100000)**(R/CP)
     *                                     !AT HALF LEVELS
     *,OROG_HEIGHT(P_FIELD)     !IN OROGRAPHIC HEIGHT.

      REAL
     * F1(U_FIELD)             !IN A CORIOLIS TERM SEE DOCUMENTATION
     *,F2(U_FIELD)             !IN A CORIOLIS TERM SEE DOCUMENTATION
     *,F3(U_FIELD)             !IN A CORIOLIS TERM SEE DOCUMENTATION
     *,COS_U_LATITUDE(U_FIELD) !IN    COS(LAT) AT U POINTS (2-D ARRAY)
     *,COS_P_LATITUDE(P_FIELD) !IN    COS(LAT) AT P POINTS (2-D ARRAY)
     *,SEC_U_LATITUDE(U_FIELD) !IN  1/COS(LAT) AT U POINTS (2-D ARRAY)
     *,SEC_P_LATITUDE(P_FIELD) !IN  1/COS(LAT) AT P POINTS (2-D ARRAY)
     *,TAN_U_LATITUDE(U_FIELD) !IN    TAN(LAT) AT U POINTS (2-D ARRAY)
     *,COS_U_LONGITUDE(ROW_LENGTH) !IN COS(LONGITUDE) AT U POINTS
     *,SIN_U_LONGITUDE(ROW_LENGTH) !IN SIN(LONGITUDE) AT U POINTS

      REAL
     * THETA_REF(P_LEVELS)    !IN REFERENCE THETA PROFILE
     *,LONGITUDE_STEP_INVERSE !IN 1/LONGITUDE INCREMENT IN RADIANS
     *,LATITUDE_STEP_INVERSE  !IN 1/LATITUDE INCREMENT IN RADIANS
     *,ADJUSTMENT_TIMESTEP    !IN
     &,ADJ_TIME_SMOOTHING_COEFF !IN COEFFICIENT. ZERO = NO SMOOTHING
     *,TRIGS(ROW_LENGTH)      !IN HOLDS TRIGONOMETRIC FUNCTIONS USED
     *                        ! IN FILTERING.

C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C    DEFINE LOCAL ARRAYS: 6 ARE REQUIRED IF TIME SMOOTHING
      REAL
     * RS_DELTAP(P_FIELD)   !HOLDS RS * VERTICAL PRESSURE DIFFERENCE
     *                      !AT P POINTS.
     *,DIVERGENCE_FUNCTIONS(P_FIELD,P_LEVELS) !WORKSPACE FOR HOLDING
     *                               !QUANTITIES INVOLVING DIVERGENCE
     *,RS_DELTAP_UV(U_FIELD) !HOLDS RS_DELTAP AT U POINTS.
     *,RECIP_RS_SQUARED_SURFACE(P_FIELD) !HOLDS 1/(RS*RS) CALCULATED AT
     *                                   ! MODEL SURFACE.
     &,U_SMOOTH(U_FIELD,P_LEVELS) ! IN ACCUMULATES U DURING ADJUSTMENT
     &,V_SMOOTH(U_FIELD,P_LEVELS) ! IN ACCUMULATES V DURING ADJUSTMENT

C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES
      INTEGER
     *  NORTHERN_FILTERED_U_ROW ! U ROW ON WHICH FITERING STOPS MOVING
     *                     ! TOWARDS EQUATOR.
     *, SOUTHERN_FILTERED_U_ROW ! U ROW ON WHICH FILTERING STARTS AGAIN
     *                     ! MOVING TOWARDS SOUTH POLE.

      INTEGER
     *  I
     *, K
     *, ADJ_STEP_NUMBER    ! USED TO HOLD THE NUMBER OF THE
     *                     ! ADJUSTMENT STEP BEING EXECUTED.
     *, FILTER_SPACE_U     ! HORIZONTAL DIMENSION OF SPACE NEEDED IN
     *                     ! FILTERING ROUTINE FOR U ROWS.
     *, FILTER_SPACE_P     ! HORIZONTAL DIMENSION OF SPACE NEEDED IN
     *                     ! FILTERING ROUTINE FOR P ROWS.

      REAL
     *  RECIP_RS_DELTAP    ! HOLDS 1./RS_DELTAP
     *, RECIP_PREF         ! 1/PREF
     *, RECIP_PREF_TO_THE_KAPPA ! 1/PREF ** KAPPA
     *, RECIP_ADJUSTMENT_STEPS
     *, SCALAR
C Local workspace arrays used in T3E restructured code
      REAL    EXNER_wk(LAST_P_VALID_PT-FIRST_VALID_PT+1)
! No. of inputs for T3E vector library function
      integer n_inputs

C*L   EXTERNAL SUBROUTINE CALLS:---------------------------------------
      EXTERNAL UV_ADJ,P_TO_UV,VERT_VEL,FILTER,POLAR_UV,
     *         P_TH_ADJ
C*---------------------------------------------------------------------
CL    CALL COMDECK TO OBTAIN CONSTANTS USED.

CLL COMDECK C_ADJCTL HOLDS CONSTANTS FOR ROUTINE ADJ_CTL
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

CL  END OF COMDECK C_ADJCTL

CL    MAXIMUM VECTOR LENGTH ASSUMED IS P_FIELD
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE INCLUDING SUBROUTINE CALLS:
CL
CL---------------------------------------------------------------------
CL    SECTION 1.     INITIALISATION
CL---------------------------------------------------------------------
C INCLUDE LOCAL CONSTANTS FROM GENERAL CONSTANTS BLOCK

      RECIP_PREF = 1./PREF
      RECIP_PREF_TO_THE_KAPPA = RECIP_PREF**KAPPA
      RECIP_ADJUSTMENT_STEPS = 1./ADJUSTMENT_STEPS

CL    IF GLOBAL THEN SET FILTERING INFORMATION.

      NORTHERN_FILTERED_U_ROW = NORTHERN_FILTERED_P_ROW
      SOUTHERN_FILTERED_U_ROW = SOUTHERN_FILTERED_P_ROW - 1

C SET FILTER_SPACE WHICH IS ROW_LENGTH+2 TIMES THE NUMBER OF ROWS TO
C BE FILTERED.

      FILTER_SPACE_U = (ROW_LENGTH+2)*(NORTHERN_FILTERED_U_ROW-1+
     *                U_FIELD/ROW_LENGTH-SOUTHERN_FILTERED_U_ROW)
      FILTER_SPACE_P = (ROW_LENGTH+2)*(NORTHERN_FILTERED_P_ROW-1+
     *                P_FIELD/ROW_LENGTH-SOUTHERN_FILTERED_P_ROW)

C SET U_MEAN, ETADOT_MEAN, AND V_MEAN  TO ZERO

      DO 102 K = 1,P_LEVELS
        DO 104 I = 1,U_FIELD
          U_MEAN(I,K) = 0.
          V_MEAN(I,K) = 0.
 104  CONTINUE
        DO 106 I = 1,P_FIELD
          ETADOT_MEAN(I,K) = 0.
          DIVERGENCE_FUNCTIONS(I,K) = 0.0
 106  CONTINUE
 102  CONTINUE

CL LOOP OVER NUMBER OF ADJUSTMENT STEPS.

      DO 110 ADJ_STEP_NUMBER = 1,ADJUSTMENT_STEPS

CL
CL---------------------------------------------------------------------
CL    SECTION 2.    CALL UV_ADJ TO ADJUST U AND V. ALSO RETURNS RS.
CL---------------------------------------------------------------------

        CALL UV_ADJ(U,V,THETA,Q,OROG_HEIGHT,PSTAR,F1,F2,
     *              F3,SEC_U_LATITUDE,TAN_U_LATITUDE,AK,BK,DELTA_AK,
     *              DELTA_BK,LATITUDE_STEP_INVERSE,ADJUSTMENT_TIMESTEP,
     *              LONGITUDE_STEP_INVERSE,RS,
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
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
     *              U_FIELD,P_FIELD,ROW_LENGTH,P_LEVELS,
     *              Q_LEVELS,ADJ_STEP_NUMBER,AKH,BKH,P_EXNER,
     *              ADJUSTMENT_STEPS,L_PHI_OUT,PHI_OUT,LLINTS,
     *              LWHITBROM)
CL
CL---------------------------------------------------------------------
CL    SECTION 3.    MASS-WEIGHTING OF U AND V.
CL---------------------------------------------------------------------

! QAN fix
          DO I=1,P_FIELD
            RS_DELTAP(I)=0.0
          ENDDO
CL LOOP OVER P_LEVELS

CMIC@ DO ALL SHARED(P_LEVELS, P_FIELD, U_FIELD, ROW_LENGTH, RS,
CMIC@1   DELTA_AK, DELTA_BK, PSTAR,  U, V) PRIVATE(RS_DELTAP_UV,
CMIC@2   RS_DELTAP, K, I)
CMIC@a   SHARED(FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
CMIC@b   P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
CMIC@c   FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
CMIC@d   FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
CMIC@e   VALID_P_ROWS, VALID_U_ROWS,
CMIC@f   START_POINT_NO_HALO, START_POINT_INC_HALO,
CMIC@g   END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
CMIC@h   END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
CMIC@i   FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
CMIC@j   GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD)
        DO 300 K = 1,P_LEVELS

CL    CALCULATE RS * DELTA P AT ALL POINTS

! loop over all points, including valid halos
          DO 310 I= FIRST_VALID_PT , LAST_P_VALID_PT
            RS_DELTAP(I) = RS(I,K)*(DELTA_AK(K) + DELTA_BK(K)*PSTAR(I))
 310      CONTINUE

CL    INTERPOLATE RS DELTAP ONTO U GRID

          CALL P_TO_UV(RS_DELTAP,RS_DELTAP_UV,P_FIELD,U_FIELD,
     &                 ROW_LENGTH,tot_P_ROWS)

CL    CALCULATE MASS WEIGHTED U AND V COS(PHI) AT ALL POINTS.

! loop over "local" points - not including top and bottom halos
          DO 320 I= FIRST_FLD_PT,LAST_U_FLD_PT
            U(I,K) = U(I,K)*RS_DELTAP_UV(I)
            V(I,K) = V(I,K)*RS_DELTAP_UV(I)
 320      CONTINUE

CL END LOOP OVER P_LEVELS

 300    CONTINUE

CL
CL---------------------------------------------------------------------
CL    SECTION 4. FILTER U,V AND DIVERGENCE FUNCTIONS IF GLOBAL MODEL.
CL---------------------------------------------------------------------


C----------------------------------------------------------------------
CL    SECTION 4.1 U_FIELD
C----------------------------------------------------------------------

CL    CALL FILTER FOR U

        CALL FILTER(U,U_FIELD,P_LEVELS,FILTER_SPACE_U,ROW_LENGTH,
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
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
     *              FILTER_WAVE_NUMBER_U_ROWS,TRIGS,IFAX,
     *              NORTHERN_FILTERED_U_ROW,SOUTHERN_FILTERED_U_ROW)

C----------------------------------------------------------------------
CL    SECTION 4.2 V_FIELD
C----------------------------------------------------------------------

CL    CALL FILTER FOR V

        CALL FILTER(V,U_FIELD,P_LEVELS,FILTER_SPACE_U,ROW_LENGTH,
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
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
     *              FILTER_WAVE_NUMBER_U_ROWS,TRIGS,IFAX,
     *              NORTHERN_FILTERED_U_ROW,SOUTHERN_FILTERED_U_ROW)

        CALL POLAR_UV(U,V,ROW_LENGTH,U_FIELD,P_LEVELS,
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
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
     &                COS_U_LONGITUDE,SIN_U_LONGITUDE)

        DO  K = 1,P_LEVELS
C     MULTIPLY V BY COS(PHI).

! loop over "local" points - not including top and bottom halos
          DO I= FIRST_FLD_PT,LAST_U_FLD_PT
            V(I,K) = V(I,K)* COS_U_LATITUDE(I)
          ENDDO
        ENDDO


! Do halo update for U and V
        CALL SWAPBOUNDS(U,ROW_LENGTH,tot_P_ROWS,
     &                  EW_Halo,NS_Halo,P_LEVELS)
        CALL SWAPBOUNDS(V,ROW_LENGTH,tot_P_ROWS,
     &                  EW_Halo,NS_Halo,P_LEVELS)
CL
CL---------------------------------------------------------------------
CL    SECTION 5. CALCULATE U_MEAN,V_MEAN AND ETA DOT.
CL---------------------------------------------------------------------

      IF(ADJ_TIME_SMOOTHING_COEFF.NE.0.0) THEN
        IF(ADJ_STEP_NUMBER.EQ.1) THEN
          DO 520 K=1,P_LEVELS
! loop over all points, including valid halos
            DO I=FIRST_VALID_PT,LAST_U_VALID_PT
              U_SMOOTH(I,K)=ADJ_TIME_SMOOTHING_WEIGHT(ADJ_STEP_NUMBER)
     &                      *U(I,K)
              V_SMOOTH(I,K)=ADJ_TIME_SMOOTHING_WEIGHT(ADJ_STEP_NUMBER)
     &                      *V(I,K)
            END DO
 520      CONTINUE
        ELSE IF(ADJ_STEP_NUMBER.EQ.ADJUSTMENT_STEPS) THEN
          DO 530 K=1,P_LEVELS
! loop over all points, including valid halos
            DO I=FIRST_VALID_PT,LAST_U_VALID_PT
              U(I,K) = U(I,K)+ADJ_TIME_SMOOTHING_COEFF
     &                 *(ADJ_TIME_SMOOTHING_WEIGHT(ADJ_STEP_NUMBER)
     &                 *U(I,K)+U_SMOOTH(I,K))
              V(I,K) = V(I,K)+ADJ_TIME_SMOOTHING_COEFF
     &                 *(ADJ_TIME_SMOOTHING_WEIGHT(ADJ_STEP_NUMBER)
     &                 *V(I,K)+V_SMOOTH(I,K))
            END DO
 530      CONTINUE
        ELSE
          DO 540 K=1,P_LEVELS
! loop over all points, including valid halos
            DO I=FIRST_VALID_PT,LAST_U_VALID_PT
              U_SMOOTH(I,K)=U_SMOOTH(I,K) +
     &                      ADJ_TIME_SMOOTHING_WEIGHT(ADJ_STEP_NUMBER) *
     &                      U(I,K)
              V_SMOOTH(I,K)=V_SMOOTH(I,K) +
     &                      ADJ_TIME_SMOOTHING_WEIGHT(ADJ_STEP_NUMBER) *
     &                      V(I,K)
            END DO
 540      CONTINUE
        END IF
      END IF

CL    CALCULATE U_MEAN AND V_MEAN AT ALL POINTS AND ALL LEVELS.

        DO 500 K = 1,P_LEVELS
! loop over all points, including valid halos
          DO 510 I = FIRST_VALID_PT,LAST_U_VALID_PT
            U_MEAN(I,K)= U_MEAN(I,K) + U(I,K) * RECIP_ADJUSTMENT_STEPS
            V_MEAN(I,K)= V_MEAN(I,K) + V(I,K) * RECIP_ADJUSTMENT_STEPS
 510      CONTINUE
 500    CONTINUE

CL    CALL VERT_VEL TO CALCULATE ETA DOT.
CL    BOTH ETA DOT FOR THIS ADJUSTMENT STEP AND THE AVERAGED VALUE
CL    ARE RETURNED.
CL    THE SUM OF THE DIVERGENCES ARE HELD AT LEVEL 1 IN THE ARRAY.

C ETA DOT FOR THIS ADJUSTMENT STEP IS RETURNED IN DIVERGENCE FUNCTIONS.

        CALL VERT_VEL(U,V,ETADOT_MEAN,SEC_P_LATITUDE,
     *                DIVERGENCE_FUNCTIONS,
     *                U_FIELD,P_FIELD,P_LEVELS,
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
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
     *                ROW_LENGTH,LATITUDE_STEP_INVERSE,
     *                LONGITUDE_STEP_INVERSE,ADJUSTMENT_STEPS,AKH,BKH,
     *                RS,ADJ_STEP_NUMBER,RECIP_RS_SQUARED_SURFACE,
     *                PSTAR,LLINTS,LWHITBROM)

! Update halos for DIVERGENCE_FUNCTIONS
        CALL SWAPBOUNDS(DIVERGENCE_FUNCTIONS,ROW_LENGTH,tot_P_ROWS,
     &                  EW_Halo,NS_Halo,P_LEVELS)

CL
CL---------------------------------------------------------------------
CL    SECTION 6. RECREATE U AND V FROM MASS-WEIGHTING U AND V COS(PHI).
CL---------------------------------------------------------------------

CL LOOP OVER P_LEVELS

        DO 600 K = 1,P_LEVELS

CL    CALCULATE RS* DELTA P AT ALL POINTS

! loop over all points, including valid halos
          DO 610 I= FIRST_VALID_PT,LAST_P_VALID_PT
            RS_DELTAP(I) = RS(I,K)*(DELTA_AK(K) + DELTA_BK(K)*PSTAR(I))
 610      CONTINUE

CL    INTERPOLATE RS DELTAP ONTO U GRID

          CALL P_TO_UV(RS_DELTAP,RS_DELTAP_UV,P_FIELD,U_FIELD,
     &                 ROW_LENGTH,tot_P_ROWS)

CL    RECREATE U AND V FORM MASS-WEIGHTED U AND V COS(PHI) AT ALL POINTS

! loop over "local" points - not including top and bottom halos
          DO 620 I= FIRST_FLD_PT,LAST_U_FLD_PT
            RECIP_RS_DELTAP = 1./ RS_DELTAP_UV(I)
            U(I,K) = U(I,K) * RECIP_RS_DELTAP
            V(I,K) = V(I,K) * RECIP_RS_DELTAP * SEC_U_LATITUDE(I)
 620      CONTINUE

CL END LOOP OVER P_LEVELS

 600    CONTINUE

! Do boundary swap for U and V
        CALL SWAPBOUNDS(U,ROW_LENGTH,tot_P_ROWS,
     &                  EW_Halo,NS_Halo,P_LEVELS)
        CALL SWAPBOUNDS(V,ROW_LENGTH,tot_P_ROWS,
     &                  EW_Halo,NS_Halo,P_LEVELS)
CL
CL---------------------------------------------------------------------
CL      SECTION 7. CALL P_TH_ADJ TO ADJUST P* AND THETA.
CL---------------------------------------------------------------------

        CALL P_TH_ADJ(PSTAR,PSTAR_OLD,THETA,THETA_REF,
     *                DIVERGENCE_FUNCTIONS,RS,DELTA_AK,DELTA_BK,
     *                P_FIELD,P_LEVELS,
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
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
     *                ADJ_STEP_NUMBER,ADJUSTMENT_TIMESTEP,
     *                ERROR_CODE,ERROR_MESSAGE,
     *                RECIP_RS_SQUARED_SURFACE,L_NEG_PSTAR)

        IF(ERROR_CODE.NE.0) RETURN
! Do boundary swap for PSTAR and THETA
        CALL SWAPBOUNDS(PSTAR,ROW_LENGTH,tot_P_ROWS,
     &                  EW_Halo,NS_Halo,1)
!        CALL SET_SIDES(PSTAR,P_FIELD,ROW_LENGTH,1,fld_type_p)
!        CALL SWAPBOUNDS(THETA,ROW_LENGTH,lasize(2),
!     &                  EW_Halo,NS_Halo,P_LEVELS)
CL
CL---------------------------------------------------------------------
CL      SECTION 8. CALCULATE P_EXNER FOR PRESSURE AT NEW TIME-LEVEL.
CL                 CALCULATION PERFORMED AT ALL HALF-LEVELS.
CL---------------------------------------------------------------------
C                                                                       
        DO 800 K=1,P_LEVELS+1

C CALCULATE EXNER AT LEVEL K - 1/2

          IF(BKH(K).EQ.0.) THEN
C IF A CONSTANT PRESSURE SURFACE SET EXNER TO HELD CONSTANT VALUE.
            DO 810 I= 1,P_FIELD
              P_EXNER(I,K) = AKH_TO_THE_KAPPA(K)
 810        CONTINUE

          ELSE IF (K.GT.1.AND.AKH(K).EQ.0.) THEN
C IF A SIGMA LEVEL THEN THE LEVEL BELOW WAS A SIGMA LEVEL AND
C EXNER CAN BE CALCULATED BY RESCALING THE VALUE AT THE LOWER LEVEL.

            SCALAR = BKH_TO_THE_KAPPA(K)/BKH_TO_THE_KAPPA(K-1)
! loop over all points, including valid halos
            DO 820 I=FIRST_VALID_PT,LAST_P_VALID_PT
              P_EXNER(I,K) = P_EXNER(I,K-1)* SCALAR
 820        CONTINUE
          ELSE
C CALCULATE EXNER AS ((A+B*PSTAR)/100000)**(R/CP)

! loop over all points, including valid halos

            DO I=FIRST_VALID_PT,LAST_P_VALID_PT
              EXNER_wk(I-FIRST_VALID_PT+1)=AKH(K)+BKH(K)*PSTAR(I)
            END DO
            DO I=1,LAST_P_VALID_PT-FIRST_VALID_PT+1
              EXNER_wk(I)=EXNER_wk(I)**KAPPA
            END DO
            DO 830 I=FIRST_VALID_PT,LAST_P_VALID_PT
              P_EXNER(I,K) = EXNER_wk(I-FIRST_VALID_PT+1)
     &                              * RECIP_PREF_TO_THE_KAPPA
 830        CONTINUE
          END IF

 800    CONTINUE


CL END OF LOOP OVER ADJUSTMENT STEPS

 110  CONTINUE
! Update halos for ETADOT_MEAN
      CALL SWAPBOUNDS(ETADOT_MEAN,ROW_LENGTH,tot_P_ROWS,
     &                EW_Halo,NS_Halo,P_LEVELS)

CL    END OF ROUTINE ADJ_CTL

      RETURN
      END

