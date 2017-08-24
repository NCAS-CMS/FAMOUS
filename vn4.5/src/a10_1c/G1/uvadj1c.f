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
CLL   SUBROUTINE UV_ADJ ---------------------------------------------
CLL
CLL   PURPOSE:  CALCULATES AND ADDS INCREMENTS TO U AND V USING
CLL             EQUATIONS 23 TO 26.
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL
CLL   WAS VERSION FOR CRAY Y-MP
CLL
CLL MM, DR      <- PROGRAMMER OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
!LL   4.4    01/08/97 New version optimised for T3E
!LL                   Not bit reproducible with UVADJ1A.
CLL   4.4    01/08/97 Optimisation for T3E removing unnecessary
CLL                   array initialisations and reworking loops
CLL                   for streams.
CLL                   Author: D.Salmond
CLL                   Reviewer: A. Dickinson
!LL   4.5    21/08/98  Comment out cdir$ cache_bypass directives due
!LL                    to t3e hardware error with new compiler.
!LL                    S.D.Mullerworth
CLL
CLL
CLL   PROGRAMMING STANDARD:
CLL
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
C DEFINE LOCAL ARRAYS:

      REAL
     * DPHI_BY_DLATITUDE(P_FIELD) !HOLDS HORIZONTAL PRESSURE GRADIENT
     *                            !IN X-DIRECTION AT U POINTS
     *,DPHI_BY_DLONGITUDE(P_FIELD)!HOLDS HORIZONTAL PRESSURE GRADIENT
     *                            !IN Y-DIRECTION AT U POINTS
     *,P(P_FIELD)                 !HOLDS PRESSURE AT A MODEL LEVEL
     *,RECIP_RS_UV(U_FIELD,P_LEVELS)     !HOLDS 1/RS AT U POINTS
     *,PHI_FULL_LEVEL(P_FIELD)    !HOLDS GEOPOTENTIAL AT A FULL LEVEL
     *,PHI_HALF_LEVEL(P_FIELD)  !HOLDS GEOPOT AT A HALF LEVEL
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
      INTEGER IP,IJP,J

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
     *  I,IJ,IK,K,II,IX1
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

! Initialise work arrays
! cdir$ cache_bypass WORK_U
        DO I=1,U_FIELD
          WORK_U(I)=0.0
        ENDDO
! cdir$ cache_bypass WORK_P
        DO I=1,P_FIELD
          WORK_P(I)=0.0
        ENDDO

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
          RECIP_RS_UV(I,K) = 1.0
 210    CONTINUE
      ELSE
        IF(CALL_NUMBER.EQ.1) THEN

! Initialise RS so that P_TO_UV works in MPP mode
          DO I=1,FIRST_VALID_PT-1
            RS(I,K)=1.0
          ENDDO
          DO I=FIRST_VALID_PT+POINTS-1,P_FIELD
            RS(I,K)=1.0
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

      DO K=1,P_LEVELS

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
            PKP1 = AKH(K+1) + BKH(K+1)*PSTAR(I)
            PK   = AKH(K)   + BKH(K)  *PSTAR(I)
            WORK_P(I) = R_P_EXNER_C
     +      (P_EXNER(I,K+1),P_EXNER(I,K),PKP1,PK,KAPPA)
          THETAS(I,K) = TS(I)*WORK_P(I)
 340    CONTINUE

C----------------------------------------------------------------------
CL    SECTION 3.5.  CALCULATE MU
CL                  CALCULATE 1/RS AT U POINTS.
C----------------------------------------------------------------------

C MU IS CALCULATED AT U POINTS AND HELD IN WORK_U
! loop over all points, including valid halos
        DO 350 I=FIRST_VALID_PT,LAST_U_VALID_PT
            RECIP_RS_UV(I,K)=1.0/RECIP_RS_UV(I,K)
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

      DO 110 K=1,P_LEVELS

CL---------------------------------------------------------------------
CL    SECTION 4.    CALCULATE PHI AT LEVEL K, EQUATION (26).
CL---------------------------------------------------------------------
C----------------------------------------------------------------------
CL    SECTION 4.1.  CALCULATE PHI AT LEVEL K
C----------------------------------------------------------------------
        TEMP2 = 1./(KAPPA+1.)

      if(k.ne.1)then

cdir$ nosplit
        DO I=FIRST_VALID_PT,LAST_P_VALID_PT
         PHI_HALF_LEVEL(I) = PHI_HALF_LEVEL(I)
     &                           -CP*THETAS(I,K-1)*
     &                            (P_EXNER(I,K) - P_EXNER(I,K-1) )
      ENDDO
cdir$ nosplit
        DO I=FIRST_VALID_PT,LAST_P_VALID_PT
         TEMP1= 1.0/(DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))
          DELTA_P_P_EXNER_BY_DELTAP(I) = (P_EXNER(I,K+1)*
     *        (AKH(K+1)+BKH(K+1)*PSTAR(I)) -
     *        P_EXNER(I,K)*(AKH(K)+BKH(K)*PSTAR(I)))
     *        *TEMP1*TEMP2
      ENDDO

cdir$ nosplit
        DO I=FIRST_VALID_PT,LAST_P_VALID_PT
         PHI_FULL_LEVEL(I) = PHI_HALF_LEVEL(I) + CP*THETAS(I,K)*
     *              (P_EXNER(I,K) - DELTA_P_P_EXNER_BY_DELTAP(I))
      ENDDO

      else if(k.eq.1)then

cdir$ nosplit
! loop over all points, including valid halos
        DO I=FIRST_VALID_PT,LAST_P_VALID_PT

        PHI_HALF_LEVEL(I) = OROG_HEIGHT(I) * G
      ENDDO
cdir$ nosplit
        DO I=FIRST_VALID_PT,LAST_P_VALID_PT
         TEMP1= 1.0/(DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))
         DELTA_P_P_EXNER_BY_DELTAP(I) = (P_EXNER(I,K+1)*
     *        (AKH(K+1)+BKH(K+1)*PSTAR(I)) -
     *        P_EXNER(I,K)*(AKH(K)+BKH(K)*PSTAR(I)))
     *        *TEMP1*TEMP2
      ENDDO
cdir$ nosplit
        DO I=FIRST_VALID_PT,LAST_P_VALID_PT
         PHI_FULL_LEVEL(I) = PHI_HALF_LEVEL(I) + CP*THETAS(I,K)*
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


CL Save East and West edges so that they are held constant in LAM mode

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


C----------------------------------------------------------------------
CL    SECTION 5.1.  CALCULATE HORIZONTAL PRESSURE GRADIENT,
CL                  D(PHI)/D(LONGITUDE).
C----------------------------------------------------------------------
C----------------------------------------------------------------------
CL    SECTION 5.2.  CALCULATE HORIZONTAL PRESSURE GRADIENT,
CL                  D(PHI)/D(LATITUDE).
C----------------------------------------------------------------------
C
C

C The following loop is unrolled to level 2 by hand since the compiler
C is not able to do this at present.
C
cdir$ nosplit
      IX1=IAND(MAX(END_U_POINT_NO_HALO-START_POINT_NO_HALO,0),1)
      IF (IX1 .EQ. 1)THEN
          I=START_POINT_NO_HALO
          IJ = I + ROW_LENGTH
          DPHI_BY_DLONGITUDE(i) = c1*SEC_U_LATITUDE(I)*(
     *                ((PHI_FULL_LEVEL(I+1)- PHI_FULL_LEVEL(I))+
     *                 (PHI_FULL_LEVEL(IJ+1)-PHI_FULL_LEVEL(IJ)))+
     *     .5*CP*(((THETAS(I+1,K)+THETAS(I,K))
     *                  *(DELTA_P_P_EXNER_BY_DELTAP(I+1) -
     *                    DELTA_P_P_EXNER_BY_DELTAP(I)))+
     *           ((THETAS(IJ+1,K)+THETAS(IJ,K))
     *                  *(DELTA_P_P_EXNER_BY_DELTAP(IJ+1) -
     *                    DELTA_P_P_EXNER_BY_DELTAP(IJ)))))
          DPHI_BY_DLATITUDE(i) = c2*(
     *                ((PHI_FULL_LEVEL(I)-PHI_FULL_LEVEL(IJ))+
     *                 (PHI_FULL_LEVEL(I+1)-PHI_FULL_LEVEL(IJ+1)))+
     *     .5*CP*(((THETAS(I,K)+THETAS(IJ,K))
     *                 *(DELTA_P_P_EXNER_BY_DELTAP(I) -
     *                    DELTA_P_P_EXNER_BY_DELTAP(IJ)))+
     *           ((THETAS(I+1,K)+THETAS(IJ+1,K))
     *                  *(DELTA_P_P_EXNER_BY_DELTAP(I+1) -
     *                    DELTA_P_P_EXNER_BY_DELTAP(IJ+1)))))
      ENDIF
        DO II=IX1 + START_POINT_NO_HALO,END_U_POINT_NO_HALO-1,2
          I = II
          IJ = I + ROW_LENGTH
          DPHI_BY_DLONGITUDE(i) = c1*SEC_U_LATITUDE(I)*(
     *                ((PHI_FULL_LEVEL(I+1)- PHI_FULL_LEVEL(I))+
     *                 (PHI_FULL_LEVEL(IJ+1)-PHI_FULL_LEVEL(IJ)))+
     *     .5*CP*(((THETAS(I+1,K)+THETAS(I,K))
     *                  *(DELTA_P_P_EXNER_BY_DELTAP(I+1) -
     *                    DELTA_P_P_EXNER_BY_DELTAP(I)))+
     *           ((THETAS(IJ+1,K)+THETAS(IJ,K))
     *                  *(DELTA_P_P_EXNER_BY_DELTAP(IJ+1) -
     *                    DELTA_P_P_EXNER_BY_DELTAP(IJ)))))
          DPHI_BY_DLATITUDE(i) = c2*(
     *                ((PHI_FULL_LEVEL(I)-PHI_FULL_LEVEL(IJ))+
     *                 (PHI_FULL_LEVEL(I+1)-PHI_FULL_LEVEL(IJ+1)))+
     *     .5*CP*(((THETAS(I,K)+THETAS(IJ,K))
     *                 *(DELTA_P_P_EXNER_BY_DELTAP(I) -
     *                    DELTA_P_P_EXNER_BY_DELTAP(IJ)))+
     *           ((THETAS(I+1,K)+THETAS(IJ+1,K))
     *                  *(DELTA_P_P_EXNER_BY_DELTAP(I+1) -
     *                    DELTA_P_P_EXNER_BY_DELTAP(IJ+1)))))
          I = II + 1
          IJ = I + ROW_LENGTH
          DPHI_BY_DLONGITUDE(i) = c1*SEC_U_LATITUDE(I)*(
     *                ((PHI_FULL_LEVEL(I+1)- PHI_FULL_LEVEL(I))+
     *                 (PHI_FULL_LEVEL(IJ+1)-PHI_FULL_LEVEL(IJ)))+
     *     .5*CP*(((THETAS(I+1,K)+THETAS(I,K))
     *                  *(DELTA_P_P_EXNER_BY_DELTAP(I+1) -
     *                    DELTA_P_P_EXNER_BY_DELTAP(I)))+
     *           ((THETAS(IJ+1,K)+THETAS(IJ,K))
     *                  *(DELTA_P_P_EXNER_BY_DELTAP(IJ+1) -
     *                    DELTA_P_P_EXNER_BY_DELTAP(IJ)))))
          DPHI_BY_DLATITUDE(i) = c2*(
     *                ((PHI_FULL_LEVEL(I)-PHI_FULL_LEVEL(IJ))+
     *                 (PHI_FULL_LEVEL(I+1)-PHI_FULL_LEVEL(IJ+1)))+
     *     .5*CP*(((THETAS(I,K)+THETAS(IJ,K))
     *                 *(DELTA_P_P_EXNER_BY_DELTAP(I) -
     *                    DELTA_P_P_EXNER_BY_DELTAP(IJ)))+
     *           ((THETAS(I+1,K)+THETAS(IJ+1,K))
     *                  *(DELTA_P_P_EXNER_BY_DELTAP(I+1) -
     *                    DELTA_P_P_EXNER_BY_DELTAP(IJ+1)))))
      ENDDO



C----------------------------------------------------------------------
CL    SECTION 5.3.  UPDATE U AND V USING IMPLICIT
CL                  TREATMENT OF CORIOLIS TERMS.
C----------------------------------------------------------------------
C This loop calculates the reciprocal on the previous pass
C in order to mask the cost of the divide.

cdir$ nosplit
! cdir$ cache_bypass f3
        DO I=START_POINT_NO_HALO,END_U_POINT_NO_HALO-1
          TEMP1 = HALF_ADJUSTMENT_TIMESTEP*
     *            (F3(I)+U(I,K)*TAN_U_LATITUDE(I)*RECIP_RS_UV(I,K))
          TEMP2 = TEMP1 * TEMP1
          RECIP=1.0/(1.+TEMP2)

          WORK_V= (V(I,K)*(1.-TEMP2)
     *                - TEMP1*(2.*U(I,K)-DPHI_BY_DLONGITUDE(I)
     *                                             *RECIP_RS_UV(I,K))
     *                - DPHI_BY_DLATITUDE(i)*RECIP_RS_UV(I,K))*RECIP
          U(I,K) = U(I,K) + TEMP1*(V(I,K)+WORK_V) -
     *                 DPHI_BY_DLONGITUDE(i)*RECIP_RS_UV(I,K)
          V(I,K) = WORK_V

      ENDDO


C Reset East West values of U and V with input values
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
