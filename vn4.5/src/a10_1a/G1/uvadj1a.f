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
CLL   SUBROUTINE UV_ADJ ---------------------------------------------
CLL
CLL   PURPOSE:  CALCULATES AND ADDS INCREMENTS TO U AND V USING
CLL             EQUATIONS 23 TO 26.
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL
CLL   VERSION FOR CRAY Y-MP
CLL
CLL MM, DR      <- PROGRAMMER OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.1     24/02/93  Tidy code to remove QA Fortran messages.
CLL   3.4     23/06/94  Argument LLINTS added and passed to CALC_TS
CLL                     DEF NOWHBR replaced by LOGICAL LWHITBROM
CLL                                                S.J.Swarbrick
CLL
CLL   3.4    06/08/94 Micro tasking directives inserted and code
CLL                   restructured to improve parallel efficiency
CLL                   on C90.
CLL                   Authors: A. Dickinson, D. Salmond
CLL                   Reviewer: M. Mawson
!     3.5    28/03/95 MPP code: Change updateable area P.Burton
!     4.1    02/04/96 Added TYPFLDPT arguments to dynamics routines
!                     which allows many of the differences between
!                     MPP and "normal" code to be at top level
!                     P.Burton
!LL   4.2    25/10/96 Initialise RECIP_RS_UV before use  P.Burton
!LL  4.2  25/11/96  Corrections to allow LAM to run in MPP mode.
!LL                                                   RTHBarnes.
!LL   4.3    17/01/97 Initialise PHI_OUT diagnostic so that halos
!LL                   contain real data                  P.Burton
C     vn4.3    Mar. 97   T3E migration : optimisation changes
C                                       D.Salmond
!LL   4.4    10/10/97 Correct loop bounds for u_field array
!LL                                                   P.Burton
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD B. VERSION 2, DATED 18/01/90
CLL   SYSTEM COMPONENTS COVERED: P111
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION:       THE EQUATIONS USED ARE (23) TO (26)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER NO. 10
CLL                        M.J.P. CULLEN,T.DAVIES, AND M.H.MAWSON
CLLEND-------------------------------------------------------------

C
C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE UV_ADJ
     1              (U,V,THETA,Q,OROG_HEIGHT,PSTAR,F1,
     2              F2,F3,SEC_U_LATITUDE,TAN_U_LATITUDE,AK,BK,DELTA_AK,
     3              DELTA_BK,LATITUDE_STEP_INVERSE,ADJUSTMENT_TIMESTEP,
     4              LONGITUDE_STEP_INVERSE,RS,
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
     5              U_FIELD,P_FIELD,ROW_LENGTH,P_LEVELS,
     6              Q_LEVELS,CALL_NUMBER,AKH,BKH,P_EXNER,
     8              ADJUSTMENT_STEPS,L_PHI_OUT,PHI_OUT,LLINTS,
     9              LWHITBROM)

      IMPLICIT NONE

      LOGICAL
     * L_PHI_OUT              !IN. TRUE IF PHI OUTPUT REQUIRED AS
     *                        !    DIAGNOSTIC.
     *,LLINTS                 !Switch for linear TS calc in CALC_TS
     *,LWHITBROM              !Switch for White & Bromley terms

      INTEGER
     *  P_FIELD               !IN DIMENSION OF FIELDS ON PRESSSURE GRID
     *, U_FIELD               !IN DIMENSION OF FIELDS ON VELOCITY GRID
     *, P_LEVELS              !IN    NUMBER OF PRESSURE LEVELS.
     *, Q_LEVELS              !IN    NUMBER OF MOIST LEVELS.
     *, ROW_LENGTH            !IN    NUMBER OF POINTS PER ROW
     *, CALL_NUMBER           !IN ADJUSTMENT STEP NUMBER
     *, ADJUSTMENT_STEPS      !IN NUMBER OF ADJUSTMENT STEPS
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

      REAL
     * THETA(P_FIELD,P_LEVELS)!INOUT THETA FIELD
     *,Q(P_FIELD,Q_LEVELS)    !INOUT Q FIELD
     *,PSTAR(P_FIELD)         !INOUT PSTAR FIELD
     *,RS(P_FIELD,P_LEVELS)   !INOUT PRIMARY MODEL ARRAY FOR RS FIELD
     *,U(U_FIELD,P_LEVELS)    !INOUT U FIELD
     *,V(U_FIELD,P_LEVELS)    !INOUT V FIELD

      REAL
     * P_EXNER(P_FIELD,P_LEVELS+1) !IN HOLDS EXNER PRESSURE AT HALF
     *                             ! LEVELS
     *,OROG_HEIGHT(P_FIELD)        !IN OROGRAPHIC HEIGHT FIELD

      REAL
     * DELTA_AK(P_LEVELS)        !IN  LAYER THICKNESS
     *,DELTA_BK(P_LEVELS)        !IN  LAYER THICKNESS
     *,AK(P_LEVELS)              !IN  VALUE OF A AT P POINTS
     *,BK(P_LEVELS)              !IN  VALUE OF B AT P POINTS
     *,AKH(P_LEVELS+1)           !IN  VALUE OF A AT HALF LEVELS.
     *,BKH(P_LEVELS+1)           !IN  VALUE OF B AT HALF LEVELS.
     *,SEC_U_LATITUDE(U_FIELD)   !IN 1/COS(LAT) AT U POINTS (2-D ARRAY)
     *,TAN_U_LATITUDE(U_FIELD)   !IN TAN(LAT) AT U POINTS (2-D ARRAY)

      REAL
     * F1(U_FIELD)            !IN A CORIOLIS TERM SEE DOCUMENTATION
     *,F2(U_FIELD)            !IN A CORIOLIS TERM SEE DOCUMENTATION
     *,F3(U_FIELD)            !IN A CORIOLIS TERM SEE DOCUMENTATION
     *,LONGITUDE_STEP_INVERSE !IN 1/LONGITUDE INCREMENT
     *,LATITUDE_STEP_INVERSE  !IN 1/LATITUDE INCREMENT
     *,ADJUSTMENT_TIMESTEP    !IN

      REAL
     * PHI_OUT(P_FIELD,P_LEVELS) !OUT. PHI DIAGNOSTIC
       REAL RECIP

C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: 15 ARE REQUIRED

      REAL
     * DPHI_BY_DLATITUDE          !HOLDS HORIZONTAL PRESSURE GRADIENT   
     *                            !IN X-DIRECTION AT U POINTS
     *,DPHI_BY_DLONGITUDE         !HOLDS HORIZONTAL PRESSURE GRADIENT  
     *                            !IN Y-DIRECTION AT U POINTS
     *,P(P_FIELD)                 !HOLDS PRESSURE AT A MODEL LEVEL
     *,RECIP_RS_UV(U_FIELD,P_LEVELS)     !HOLDS 1/RS AT U POINTS
     *,PHI_FULL_LEVEL(P_FIELD)    !HOLDS GEOPOTENTIAL AT A FULL LEVEL
     *,PHI_HALF_LEVEL(P_FIELD,P_LEVELS)  !HOLDS GEOPOT AT A HALF LEVEL
     *,DELTA_P_P_EXNER_BY_DELTAP(P_FIELD) !

      REAL
     * THETAS(P_FIELD,P_LEVELS)   !HOLDS THETAV + MU*THETAS
     *,TS(P_FIELD)                !HOLDS STANDARD TEMPERATURE
     *,WORK_U(U_FIELD)            !GENERAL WORKSPACE FOR VARIABLES
     *                            !AT U POINTS
     *,WORK_P(P_FIELD)            !GENERAL WORKSPACE FOR VARIABLES
     *                            !AT P POINTS
     *,U_TEMP_R(U_FIELD),V_TEMP_R(U_FIELD)
     *,U_TEMP_L(U_FIELD),V_TEMP_L(U_FIELD)
      INTEGER IP,IJP

C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES
      INTEGER POINTS  ! Number of points with valid part of field

      INTEGER row_start_offset,row_end_offset
! offsets required to mark out the updatable area for LAM MPP code
       REAL
     *  HALF_ADJUSTMENT_TIMESTEP
     *, RECIP_G

C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     *  I,IJ,IK,K
C WORK-SPACE SCALARS
      REAL
     *  TEMP1,TEMP2
     * ,PKP1,PK               ! Pressures at half levels k+1 and k
     * ,c1,c2,WORK_V
C LOGICAL VARIABLE
      LOGICAL
     *  CONSTANT_PRESSURE     ! TRUE IF ON A CONSTANT PRESSURE SURFACE

C*L   EXTERNAL SUBROUTINE CALLS:---------------------------------------

      EXTERNAL P_TO_UV,POLAR_UV,UV_TO_P
     *         ,CALC_TS,CALC_RS
C*---------------------------------------------------------------------
CL CALL COMDECK TO OBTAIN CONSTANTS USED.

CLL COMDECK C_UVADJ HOLDS CONSTANTS FOR ROUTINE UV_ADJ.
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

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
C*----------------------------------------------------------------------

CL END OF COMDECK C_UVADJ

CL    MAXIMUM VECTOR LENGTH ASSUMED IS P_FIELD
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
CL    SECTION 1.    INITIALISATION
CL---------------------------------------------------------------------
C INCLUDE LOCAL CONSTANTS FROM GENERAL CONSTANTS BLOCK

      POINTS=LAST_P_VALID_PT-FIRST_VALID_PT+1
! Number of points to be processed by CALC_RS/TS. For non-MPP runs
! this is simply P_FIELD, for MPP, it is all the points, minus any
! unused halo areas (ie. the halo above North pole row, and beneath
! South pole row)

      HALF_ADJUSTMENT_TIMESTEP = ADJUSTMENT_TIMESTEP*.5
      RECIP_G = 1./G

CL    SET PHI_HALF_LEVEL FOR LEVEL 1/2 = OROG_HEIGHT * G
! loop over all points, including valid halos
      DO 100 I=FIRST_VALID_PT,LAST_P_VALID_PT
        PHI_HALF_LEVEL(I,1) = OROG_HEIGHT(I) * G
 100  CONTINUE

CL LOOP OVER ALL PRESSURE LEVELS.

      DO K=1,P_LEVELS

CL---------------------------------------------------------------------
CL   IF (.NOT.LWHITBROM) THEN
CL    SECTION 2.    STORE RADIUS OF EARTH IN HORIZONTAL FIELD.
CL   ELSE
CL    SECTION 2.    CALCULATE RS AT LEVEL K.
CL   END IF
CL---------------------------------------------------------------------

C----------------------------------------------------------------------
CL   IF (.NOT.LWHITBROM) THEN
CL    SECTION 2.1.  STORE RADIUS OF EARTH IN HORIZONTAL FIELD.
CL   ELSE
CL    SECTION 2.1.  CALL CALC_RS TO GET RS ON FIRST CALL ONLY.
CL                  ALSO RETURNS TS SAVING CALL TO CALC_TS IN 3.4
CL   END IF
C----------------------------------------------------------------------

      IF (.NOT.LWHITBROM) THEN
! loop over all points, including valid halos
        DO 210 I=1,P_FIELD
          RS(I,K) = A
 210    CONTINUE
        DO I=1,U_FIELD
          RECIP_RS_UV(I,K)=1.0
        ENDDO
      ELSE
        IF(CALL_NUMBER.EQ.1) THEN
! QAN fix
          DO I=1,P_FIELD
            RS(I,K)=1.0
          ENDDO
          DO I=1,U_FIELD
            RECIP_RS_UV(I,K)=1.0
          ENDDO
          IF(K.NE.1) THEN
            CALL CALC_RS(PSTAR(FIRST_VALID_PT),AK,BK,TS(FIRST_VALID_PT),
     &                   RS(FIRST_VALID_PT,K-1),
     &                   RS(FIRST_VALID_PT,K),
     &                   POINTS,K,P_LEVELS,LLINTS)
          ELSE
C IF LEVEL 1 CALC_RS NEEDS A DUMMY ARRAY IN PLACE OF RS( ,K-1)
            CALL CALC_RS(PSTAR(FIRST_VALID_PT),AK,BK,TS(FIRST_VALID_PT),
     &                   RS(FIRST_VALID_PT,K+1),
     &                   RS(FIRST_VALID_PT,K),
     &                   POINTS,K,P_LEVELS,LLINTS)
          END IF
        END IF
      ENDIF ! LWHITBROM

C----------------------------------------------------------------------
CL   IF (.NOT.LWHITBROM) THEN
CL    SECTION 2.2.  STORE 1./RADIUS OF EARTH IN HORIZONTAL FIELD.
CL   ELSE
CL    SECTION 2.2.  CALL P_TO_UV TO GET RS AT U POINTS.
CL   END IF
C----------------------------------------------------------------------

      IF (.NOT.LWHITBROM) THEN
! loop over all points, including valid halos
         DO 220 I=FIRST_VALID_PT,LAST_U_VALID_PT
            RECIP_RS_UV(I,K) = 1./A
 220     CONTINUE
      ELSE
C STORE RS AT U POINTS IN RECIP_RS_UV

         CALL P_TO_UV(RS(1,K),RECIP_RS_UV(1,K),P_FIELD,
     &                U_FIELD,ROW_LENGTH,tot_P_ROWS)
! loop over "local" points - not including top and bottom halos
        DO  I=FIRST_FLD_PT,LAST_U_FLD_PT
            RECIP_RS_UV(I,K) = 1./RECIP_RS_UV(I,K)
       ENDDO

      ENDIF
      ENDDO
      IF (LWHITBROM) THEN
        CALL SWAPBOUNDS(RECIP_RS_UV,ROW_LENGTH,tot_P_ROWS,
     &                  EW_Halo,NS_Halo,P_LEVELS)
      ENDIF
CL---------------------------------------------------------------------
CL    SECTION 3.    CALCULATE PHI AT LEVEL K-1/2, EXNER AT LEVEL K,
CL   IF (.NOT.LWHITBROM) THEN
CL                  AND THETAV.
CL   ELSE
CL                  AND THETAV + MU * THETAS AT LEVEL K.
CL   END IF
CL---------------------------------------------------------------------

cmic$ do all shared (adjustment_timestep, ak, akh, bk, bkh, c_virtual)
cmic$*       shared (cp, delta_ak, delta_bk)
cmic$*       shared (epsilon, f1, f2, f3, half_adjustment_timestep)
cmic$*       shared ( kappa)
cmic$*       shared (longitude_step_inverse, p_exner, p_field)
cmic$*       shared (p_levels, phi_half_level, phi_out)
cmic$*       shared (pstar, q, q_levels, r, recip_g, row_length)
cmic$*       shared (rs, sec_u_latitude, points)
cmic$*       shared (tan_u_latitude, theta, thetas, u, u_field, v)
cmic$*       shared (call_number, lwhitbrom, llints)
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
cmic$*       private (constant_pressure, delta_p_p_exner_by_deltap)
cmic$*       private (dphi_by_dlatitude, dphi_by_dlatitude_p)
cmic$*       private (dphi_by_dlongitude, dphi_by_dlongitude_p)
cmic$*       private(dppebd_by_dlatitude_p,dppebd_by_dlongitude_p,i,ij)
cmic$*       private (k, kappa_dum, p, p_exl_dum, p_exu_dum)
cmic$*       private (phi_full_level, pk, pkp1, pl_dum, pu_dum)
cmic$*       shared (recip_rs_uv)
cmic$*       private (temp1, temp2, ts, work_p, work_u)
cmic$*       shared (MU_NORTH_POLE, MU_SOUTH_POLE)
      DO K=1,P_LEVELS

C----------------------------------------------------------------------
CL    SECTION 3.2.  CALCULATE EXNER AT LEVEL K.
C----------------------------------------------------------------------
C STORE EXNER AT LEVEL K IN WORK_P
! loop over all points, including valid halos
          DO 320 I= FIRST_VALID_PT,LAST_P_VALID_PT
            PKP1 = AKH(K+1) + BKH(K+1)*PSTAR(I)
            PK   = AKH(K)   + BKH(K)  *PSTAR(I)
            WORK_P(I) = P_EXNER_C
     +      (P_EXNER(I,K+1),P_EXNER(I,K),PKP1,PK,KAPPA)
 320      CONTINUE

      IF (LWHITBROM) THEN
C----------------------------------------------------------------------
CL    SECTION 3.3.  CALCULATES PRESSURE AT LEVEL K NEEDED FOR CALL
CL                  TO CALC_TS. PERFORMED ONLY IF CALL_NUMBER > 1.
C----------------------------------------------------------------------

          IF(BK(K).EQ.0.) THEN
C SET CONSTANT_PRESSURE BEFORE CALL TO TS AND P AT START ADDRESS AS
C THIS IS ALL TS NEEDS IN THIS CASE.
            CONSTANT_PRESSURE = .TRUE.
            P(FIRST_VALID_PT) = AK(K)
          ELSE
C SET CONSTANT_PRESSURE BEFORE CALL TO TS AND P.
! loop over all points, including valid halos
            DO 330 I=FIRST_VALID_PT,LAST_P_VALID_PT
              P(I) = AK(K) + BK(K)*PSTAR(I)
 330        CONTINUE
            CONSTANT_PRESSURE = .FALSE.
          END IF

C----------------------------------------------------------------------
CL    SECTION 3.4.  CALL CALC_TS TO GET STANDARD TEMPERATURE.
CL                  ONLY CALLED IF CALL_NUMBER GREATER THAN 1
CL                  AS TS CALCULATED IN SECTION 2.1 ON CALL_NUMBER 1.
CL                  THEN CALCULATE THETAS BY DIVIDING BY EXNER.
C----------------------------------------------------------------------
C EXNER AT LEVEL K IS IN WORK_P

          CALL CALC_TS(P(FIRST_VALID_PT),TS(FIRST_VALID_PT),POINTS,
     &                 CONSTANT_PRESSURE,LLINTS)

C       Convert TS to THETAS
! loop over all valid points - including top and bottom halos
        DO 340 I=FIRST_VALID_PT,LAST_P_VALID_PT
          THETAS(I,K) = TS(I)/WORK_P(I)
 340    CONTINUE

C----------------------------------------------------------------------
CL    SECTION 3.5.  CALCULATE MU
CL                  CALCULATE 1/RS AT U POINTS.
C----------------------------------------------------------------------

C MU IS CALCULATED AT U POINTS AND HELD IN WORK_U
! QAN fix
        DO I=1,U_FIELD
          WORK_U(I)=0.0
        ENDDO
! loop over all points, including valid halos
        DO 350 I=FIRST_VALID_PT,LAST_U_VALID_PT
            WORK_U(I) = (F2(I)*U(I,K) - F1(I)*V(I,K) +
     *              (U(I,K)*U(I,K)+V(I,K)*V(I,K))*RECIP_RS_UV(I,K))*
     *              RECIP_G
 350    CONTINUE
C CALL UV_TO_P TO INTERPOLATE MU ONTO P-GRID HELD IN WORK_P

        CALL UV_TO_P(WORK_U(START_POINT_NO_HALO-ROW_LENGTH),
     &               WORK_P(START_POINT_NO_HALO),
     &               U_FIELD-(START_POINT_NO_HALO-ROW_LENGTH)+1,
     &               P_FIELD-START_POINT_NO_HALO+1,
     &               ROW_LENGTH,upd_P_ROWS+1)

! Set WORK at North and South edges to one row in
        IF (at_top_of_LPG) THEN
          DO I=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
            WORK_P(I) = WORK_P(I+ROW_LENGTH)
          ENDDO
        ENDIF
        IF (at_base_of_LPG) THEN
          DO I=P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
            WORK_P(I) = WORK_P(I-ROW_LENGTH)
          ENDDO
        ENDIF

C----------------------------------------------------------------------
CL    SECTION 3.6.  CALCULATE THETAV + MU * THETAS
C----------------------------------------------------------------------

        IF(K.LE.Q_LEVELS) THEN
! loop over all points - including top and bottom halos
        DO 360 I=FIRST_VALID_PT,LAST_P_VALID_PT
            THETAS(I,K) = THETA(I,K)*(1.+ C_VIRTUAL
     *                     *Q(I,K))+ WORK_P(I)*THETAS(I,K)
 360      CONTINUE
        ELSE
! loop over all points - including top and bottom halos
        DO 362 I=FIRST_VALID_PT,LAST_P_VALID_PT
            THETAS(I,K) = THETA(I,K) + WORK_P(I)*THETAS(I,K)
 362      CONTINUE
        END IF

      ELSE      !   LWHITBROM

C----------------------------------------------------------------------
CL    SECTION 3.3.  CALCULATE THETAV
C----------------------------------------------------------------------

        IF(K.LE.Q_LEVELS) THEN
! loop over all points, including valid halos
        DO 460 I=FIRST_VALID_PT,LAST_P_VALID_PT
            THETAS(I,K) = THETA(I,K)*(1.+ C_VIRTUAL
     *                     *Q(I,K))
 460      CONTINUE
        ELSE
! loop over all points, including valid halos
        DO 462 I=FIRST_VALID_PT,LAST_P_VALID_PT
            THETAS(I,K) = THETA(I,K)
 462      CONTINUE
        END IF

      END IF    !   LWHITBROM


      ENDDO

      IF (LWHITBROM) THEN
        CALL SWAPBOUNDS(THETAS,ROW_LENGTH,tot_P_ROWS,
     &                  EW_Halo,NS_Halo,P_LEVELS)
      ENDIF

        c1=.5*LONGITUDE_STEP_INVERSE*ADJUSTMENT_TIMESTEP
        c2=.5*LATITUDE_STEP_INVERSE*ADJUSTMENT_TIMESTEP

cmic$ do all shared (adjustment_timestep, ak, akh, bk, bkh, c_virtual)
cmic$*       shared (cp, delta_ak, delta_bk)
cmic$*       shared (epsilon, f1, f2, f3, half_adjustment_timestep)
cmic$*       shared ( kappa)
cmic$*       shared (l_phi_out, latitude_step_inverse)
cmic$*       shared (longitude_step_inverse, p_exner, p_field)
cmic$*       shared (p_levels, phi_half_level, phi_out)
cmic$*       shared (pstar, q, q_levels, r, recip_g, row_length)
cmic$*       shared (rs, sec_u_latitude)
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
cmic$*       shared (tan_u_latitude, theta, thetas, u, u_field, v)
cmic$*       private (constant_pressure, delta_p_p_exner_by_deltap)
cmic$*       private (dphi_by_dlatitude, dphi_by_dlatitude_p)
cmic$*       private (dphi_by_dlongitude, dphi_by_dlongitude_p)
cmic$*       private(dppebd_by_dlatitude_p,dppebd_by_dlongitude_p,i,ij)
cmic$*       private (k, kappa_dum, p, p_exl_dum, p_exu_dum)
cmic$*       private (phi_full_level, pk, pkp1, pl_dum, pu_dum)
cmic$*       shared (recip_rs_uv)
cmic$*       private (ik,temp1, temp2, ts, work_p, work_u)
cmic$*       shared (c1,c2)
cmic$*       private (work_v,u_temp,v_temp)
      DO 110 K=1,P_LEVELS

CL--------------------------------------------------------------------- 
CL    SECTION 4.    CALCULATE PHI AT LEVEL K, EQUATION (26).            
CL--------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
CL    SECTION 4.1.  CALCULATE PHI AT LEVEL K                            
C---------------------------------------------------------------------- 
        TEMP2 = 1./(KAPPA+1.)

      if(k.ne.p_levels.and.k.ne.1)then

cdir$ nosplit
! loop over all points, including valid halos
        DO I=FIRST_VALID_PT,LAST_P_VALID_PT                         
         PHI_HALF_LEVEL(I,K) = PHI_HALF_LEVEL(I,K)+PHI_HALF_LEVEL(I,K-1)
         PHI_HALF_LEVEL(I,K+1) = -CP*THETAS(I,K)*
     &                            (P_EXNER(I,K+1) - P_EXNER(I,K) )
          DELTA_P_P_EXNER_BY_DELTAP(I) = (P_EXNER(I,K+1)*
     *        (AKH(K+1)+BKH(K+1)*PSTAR(I)) -
     *        P_EXNER(I,K)*(AKH(K)+BKH(K)*PSTAR(I)))
     *        / (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))*TEMP2
         PHI_FULL_LEVEL(I) = PHI_HALF_LEVEL(I,K) + CP*THETAS(I,K)*     
     *              (P_EXNER(I,K) - DELTA_P_P_EXNER_BY_DELTAP(I))       
      ENDDO

      else if(k.eq.1)then

cdir$ nosplit
! loop over all points, including valid halos                           
        DO I=FIRST_VALID_PT,LAST_P_VALID_PT
         PHI_HALF_LEVEL(I,K+1) = -CP*THETAS(I,K)*
     &                            (P_EXNER(I,K+1) - P_EXNER(I,K) )
         DELTA_P_P_EXNER_BY_DELTAP(I) = (P_EXNER(I,K+1)*
     *        (AKH(K+1)+BKH(K+1)*PSTAR(I)) -
     *        P_EXNER(I,K)*(AKH(K)+BKH(K)*PSTAR(I)))
     *        / (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))*TEMP2
         PHI_FULL_LEVEL(I) = PHI_HALF_LEVEL(I,K) + CP*THETAS(I,K)*
     *              (P_EXNER(I,K) - DELTA_P_P_EXNER_BY_DELTAP(I))
      ENDDO

      else if(k.eq.p_levels)then

cdir$ nosplit
! loop over all points, including valid halos                           
        DO I=FIRST_VALID_PT,LAST_P_VALID_PT
         PHI_HALF_LEVEL(I,K) = PHI_HALF_LEVEL(I,K)+PHI_HALF_LEVEL(I,K-1)
         DELTA_P_P_EXNER_BY_DELTAP(I) = (P_EXNER(I,K+1)*
     *        (AKH(K+1)+BKH(K+1)*PSTAR(I)) -
     *        P_EXNER(I,K)*(AKH(K)+BKH(K)*PSTAR(I)))
     *        / (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))*TEMP2
C CALCULATE PHI AT LEVEL K
          PHI_FULL_LEVEL(I) = PHI_HALF_LEVEL(I,K) + CP*THETAS(I,K)*
     *              (P_EXNER(I,K) - DELTA_P_P_EXNER_BY_DELTAP(I))
        ENDDO

       endif
                                                                        
CL    COPY PHI_FULL_LEVEL INTO OUTPUT ARRAY IF DIAGNOSTIC REQUIRED.

        IF(L_PHI_OUT) THEN
! loop over all points, including valid halos
          DO I=FIRST_VALID_PT,LAST_P_VALID_PT
            PHI_OUT(I,K) = PHI_FULL_LEVEL(I)
          END DO
! Initialise whole array so there are no NaNs for STASH to fall
! over on
          DO I=1,FIRST_VALID_PT-1
            PHI_OUT(I,K)=PHI_OUT(FIRST_VALID_PT,K)
          ENDDO
          DO I=LAST_P_VALID_PT+1,P_FIELD
            PHI_OUT(I,K)=PHI_OUT(LAST_P_VALID_PT,K)
          ENDDO
        END IF

CL---------------------------------------------------------------------
CL    SECTION 5.    CALCULATE HORIZONTAL PRESSURE GRADIENTS.
CL                  THEN CALCULATE CORIOLIS TERM AND IMPLICITLY UPDATE
CL                  U AND V. EQUATIONS (23) TO (25).
CL---------------------------------------------------------------------
C----------------------------------------------------------------------
CL    SECTION 5.1.  CALCULATE HORIZONTAL PRESSURE GRADIENT,
CL                  D(PHI)/D(LONGITUDE).
C----------------------------------------------------------------------
C---------------------------------------------------------------------- 
CL    SECTION 5.2.  CALCULATE HORIZONTAL PRESSURE GRADIENT,             
CL                  D(PHI)/D(LATITUDE).                                 
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
CL    SECTION 5.3.  UPDATE U AND V USING IMPLICIT                       
CL                  TREATMENT OF CORIOLIS TERMS.                        
C---------------------------------------------------------------------- 

        IF (at_left_of_LPG) THEN
        DO I=START_POINT_NO_HALO + FIRST_ROW_PT-1,
     &       END_U_POINT_NO_HALO,ROW_LENGTH
        U_TEMP_L(I)=U(I,K)
        V_TEMP_L(I)=V(I,K)
        ENDDO
        ENDIF
        IF (at_right_of_LPG) THEN
        DO I=START_POINT_NO_HALO + LAST_ROW_PT-1,
     &       END_U_POINT_NO_HALO,ROW_LENGTH
        U_TEMP_R(I)=U(I,K)
        V_TEMP_R(I)=V(I,K)
        U_TEMP_R(I-1)=U(I-1,K)
        V_TEMP_R(I-1)=V(I-1,K)
          ENDDO
        ENDIF


C LOOP OVER ALL POINTS TO BE UPDATED.                                   

cdir$ nosplit
cdir$ nounroll

        DO 530 I=START_POINT_NO_HALO,END_U_POINT_NO_HALO-1

          TEMP1 = HALF_ADJUSTMENT_TIMESTEP*
     *            (F3(I)+U(I,K)*TAN_U_LATITUDE(I)*RECIP_RS_UV(I,K))
          TEMP2 = TEMP1 * TEMP1
          RECIP=1.0/(1.+TEMP2)

          IJ = I + ROW_LENGTH
          DPHI_BY_DLONGITUDE = c1*(
     *                  (PHI_FULL_LEVEL(I+1)-PHI_FULL_LEVEL(I))+
     *                  (PHI_FULL_LEVEL(IJ+1)-PHI_FULL_LEVEL(IJ))+
     *     .5*CP*(THETAS(I+1,K)+THETAS(I,K))
     *                  *(DELTA_P_P_EXNER_BY_DELTAP(I+1) -
     *                    DELTA_P_P_EXNER_BY_DELTAP(I))+
     *     .5*CP*(THETAS(IJ+1,K)+THETAS(IJ,K))
     *                  *(DELTA_P_P_EXNER_BY_DELTAP(IJ+1) -
     *                    DELTA_P_P_EXNER_BY_DELTAP(IJ)))*
     *     SEC_U_LATITUDE(I)*RECIP_RS_UV(I,K)
          DPHI_BY_DLATITUDE = c2*(
     *                  (PHI_FULL_LEVEL(I)-PHI_FULL_LEVEL(IJ))+
     *                  (PHI_FULL_LEVEL(I+1)-PHI_FULL_LEVEL(IJ+1))+
     *     .5*CP*(THETAS(I,K)+THETAS(IJ,K))
     *                 *(DELTA_P_P_EXNER_BY_DELTAP(I) -
     *                    DELTA_P_P_EXNER_BY_DELTAP(IJ))+
     *     .5*CP*(THETAS(I+1,K)+THETAS(IJ+1,K))
     *                  *(DELTA_P_P_EXNER_BY_DELTAP(I+1) -
     *                    DELTA_P_P_EXNER_BY_DELTAP(IJ+1)))*
     *                                RECIP_RS_UV(I,K)

C CALCULATE V AT NEW TIME LEVEL.
C WORK_V HOLDS V AT NEW TIME-LEVEL.

          WORK_V= (V(I,K)*(1.-TEMP2)
     *                - TEMP1*(2.*U(I,K)-DPHI_BY_DLONGITUDE)
     *                - DPHI_BY_DLATITUDE)*RECIP

C CALCULATE U AT NEW TIME-LEVEL.

          U(I,K) = U(I,K) + TEMP1*(V(I,K)+WORK_V) -
     *                 DPHI_BY_DLONGITUDE

C SET V EQUAL TO V AT NEW TIME-LEVEL.

          V(I,K) = WORK_V

 530    CONTINUE
                                                        

C END LOOP OVER ALL POINTS TO BE UPDATED.                               

        IF (at_left_of_LPG) THEN
        DO I=START_POINT_NO_HALO + FIRST_ROW_PT-1,
     &       END_U_POINT_NO_HALO,ROW_LENGTH
        U(I,K)=U_TEMP_L(I)
        V(I,K)=V_TEMP_L(I)
        ENDDO
        ENDIF
        IF (at_right_of_LPG) THEN
        DO I=START_POINT_NO_HALO + LAST_ROW_PT-1,
     &       END_U_POINT_NO_HALO,ROW_LENGTH
        U(I,K)=U_TEMP_R(I)
        V(I,K)=V_TEMP_R(I)
        U(I-1,K)=U_TEMP_R(I-1)
        V(I-1,K)=V_TEMP_R(I-1)
        ENDDO
        ENDIF

CL END LOOP OVER P_LEVELS
 110  CONTINUE

CL END OF ROUTINE UV_ADJ

      RETURN
      END
