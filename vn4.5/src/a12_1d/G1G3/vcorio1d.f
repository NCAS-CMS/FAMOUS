C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1996, METEOROLOGICAL OFFICE, All Rights Reserved.
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
CLL   SUBROUTINE V_CORIOL -----------------------------------------
CLL
CLL   PURPOSE:   CALCULATES APPROXIMATE VERTICAL VELOCITY AS IN
CLL              EQUATION (46) AT A MODEL LEVEL.
CLL      NOT SUITABLE FOR SINGLE COLUMN USE.
CLL      VERSION FOR CRAY Y-MP
CLL
CLL   WRITTEN BY M.H MAWSON.
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 4.2:
CLL VERSION  DATE
!LL   4.2   25/10/96  New deck for HADCM2-specific section A12_1D,
!LL                   as VCORIO1A but with reintroduced errors in
!LL                   calculation of WP, loops 210,212,220,222 and 224.
!LL                   T.Johns
!LL   4.3   10/04/97  Updated in line with MPP optimisations.  T.Johns
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD B. VERSION 2, DATED 18/01/90
CLL
CLL   SYSTEM COMPONENTS COVERED: P124
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION:       THE EQUATIONS USED ARE (44) TO (46)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER
CLL                        NO. 10 M.J.P. CULLEN, T.DAVIES AND
CLL                        M.H.MAWSON, VERSION 10, DATED 10/09/90.
CLLEND-------------------------------------------------------------

C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE V_CORIOL
     1                   (ETADOT_MINUS,ETADOT_PLUS,PSTAR,PSTAR_OLD,
     2                   U,V,RS,SEC_U_LATITUDE,ADVECTION_TIMESTEP,AK,
     3                   BK,DELTA_AK,DELTA_BK,DELTA_AKH_MINUS,
     4                   DELTA_BKH_MINUS,DELTA_AKH_PLUS,DELTA_BKH_PLUS,
     5                   ROW_LENGTH,
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
     6                   LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     7                   WK,U_FIELD,OMEGA,LLINTS)

      IMPLICIT NONE

      INTEGER
     *  U_FIELD            !IN DIMENSION OF FIELDS ON VELOCITY GRID
     *, ROW_LENGTH         !IN NUMBER OF POINTS PER ROW

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
     * U(U_FIELD)               !IN AVERAGED MASS-WEIGHTED U VELOCITY
     *                          !   FROM ADJUSTMENT STEP HELD AT P
     *                          !   POINTS WITH FIRST POINT OF FIELD
     *                          !   BEING FIRST P POINT ON SECOND ROW
     *                          !   OF P-GRID.
     *,V(U_FIELD)               !IN AVERAGED MASS-WEIGHTED V VELOCITY
     *                          !   * COS(LAT) FROM ADJUSTMENT STEP
     *                          !   STORAGE AS FOR U_MEAN.
     *,ETADOT_PLUS(U_FIELD)     !IN AVERAGED MASS-WEIGHTED
     *                          !VERTICAL VELOCITY FROM ADJUSTMENT STEP
     *                          ! AT LEVEL K+1/2.
     *,ETADOT_MINUS(U_FIELD)    !IN AVERAGED MASS-WEIGHTED
     *                          !VERTICAL VELOCITY FROM ADJUSTMENT STEP
     *                          ! AT LEVEL K-1/2.

      REAL
     * PSTAR(U_FIELD)           !IN PSTAR FIELD AT NEW TIME-LEVEL ON
     *                          ! U GRID.
     *,PSTAR_OLD(U_FIELD)       !INPSTAR AT PREVIOUS TIME-LEVEL ON
     *                          ! U GRID.
     *,RS(U_FIELD)              !IN RS FIELD ON U GRID.
     *,SEC_U_LATITUDE(U_FIELD)  !IN  1/COS(LAT) AT U POINTS (2-D ARRAY)
     *,LONGITUDE_STEP_INVERSE   !IN 1/LONGITUDE STEP
     *,LATITUDE_STEP_INVERSE    !IN 1/LATITUDE STEP
     *,ADVECTION_TIMESTEP       !IN

      REAL
     * WK(U_FIELD)              !OUT WK AS IN EQUATION (46).
     *,OMEGA(U_FIELD)           !OUT. HOLDS VERTICAL VELOCITY, OMEGA.


      REAL
     * AK                       !IN FIRST TERM IN HYBRID CO-ORDS.
     *,BK                       !IN SECOND TERM IN HYBRID CO-ORDS.
     *,DELTA_AK                 !IN LAYER THICKNESS
     *,DELTA_BK                 !IN LAYER THICKNESS
     *,DELTA_AKH_MINUS          !IN LAYER THICKNESS  AK(K) - AK(K-1)
     *,DELTA_BKH_MINUS          !IN LAYER THICKNESS  BK(K) - BK(K-1)
     *,DELTA_AKH_PLUS           !IN LAYER THICKNESS  AK(K+1) - AK(K)
     *,DELTA_BKH_PLUS           !IN LAYER THICKNESS  BK(K+1) - BK(K)
C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: 5 ARE REQUIRED

      REAL
     * DP_BY_DT(U_FIELD)
     *,WP(U_FIELD)
     *,WORK1(U_FIELD)
     *,WORK2(U_FIELD)
     *,TS(U_FIELD)

C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES
      REAL
     *  SCALAR

C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     *  I, POINTS

C LOGICAL VARIABLES
      LOGICAL
     * CONSTANT_PRESSURE
     *,LLINTS               ! Switch for linear TS calc in CALC_TS

C*L   EXTERNAL SUBROUTINE CALLS:---------------------------------------
      EXTERNAL CALC_TS
C*---------------------------------------------------------------------
CL    CALL COMDECK TO OBTAIN CONSTANTS USED.

CLL   COMDECK C_VCORI HOLDS CONSTANTS FOR ROUTINE V_CORIOL.
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

CL    END OF COMDECK C_VCORI

CL    MAXIMUM VECTOR LENGTH ASSUMED IS END_U_UPDATE+ROW_LENGTH+1-
CL                                   START_U_UPDATE
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE INCLUDING SUBROUTINE CALLS:
CL---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL    SECTION 1.     CALCULATE DP/DT
CL---------------------------------------------------------------------
!! First update halos of u and v fields since the indexing errors in 
!! later calculations are otherwise not using the correct data  
      CALL SWAPBOUNDS(u,local_row_length,tot_u_rows,ew_halo,ns_halo,1)
      CALL SWAPBOUNDS(v,local_row_length,tot_u_rows,ew_halo,ns_halo,1)

      IF(BK.EQ.0.) THEN
C A CONSTANT PRESSURE LEVEL SO DP/DT IS ZERO.
        CONSTANT_PRESSURE = .TRUE.
! Loop over U field missing top and bottom rows and halos
        DO 100 I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
          DP_BY_DT(I) = 0.
 100    CONTINUE
      ELSE
C CALCULATE DP/DT.
        CONSTANT_PRESSURE = .FALSE.
        SCALAR = BK/ADVECTION_TIMESTEP
! Loop over U field missing top and bottom rows and halos
        DO 110 I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
          DP_BY_DT(I) = (DELTA_AK+DELTA_BK*PSTAR_OLD(I))*RS(I)*RS(I)
     *                  *(PSTAR(I)-PSTAR_OLD(I))*SCALAR
 110    CONTINUE
      END IF

CL---------------------------------------------------------------------
CL    SECTION 2.     CALCULATE U.GRAD P
CL---------------------------------------------------------------------

C----------------------------------------------------------------------
CL    SECTION 2.1    CALCULATE U DP/D(LONGITUDE)
C----------------------------------------------------------------------

C CALCULATE U DP/D(LONGITUDE) BETWEEN P POINTS
! Loop over U field missing top and bottom rows and halos and
! last point (includes HADCM2-specific error)
      DO 210 I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
        WORK1(I) = .5*(U(I)+U(I-ROW_LENGTH))*(PSTAR(I+1)-PSTAR(I))*
     *             LONGITUDE_STEP_INVERSE*BK
 210  CONTINUE

      WORK1(END_U_POINT_NO_HALO) = 0.0
! MPP Code : No need to do recalculations of end points because cyclic
! boundary conditions means that halos do this for us automatically


C CALCULATE U DP/D(LONGITUDE) AT P POINTS

! Loop over U field missing top and bottom rows and halos and
! first point
      DO 214 I=START_POINT_NO_HALO+1,END_U_POINT_NO_HALO
        WP(I) = .5*(WORK1(I)+WORK1(I-1))
 214  CONTINUE

C----------------------------------------------------------------------
CL    SECTION 2.2    CALCULATE V DP/D(LATITUDE) AND HENCE U.GRAD P
C----------------------------------------------------------------------

C CALCULATE V DP/D(LATITUDE) BETWEEN P POINTS.

! Loop over U field missing bottom row, last point and top and
! bottom halos (includes HADCM2-specific error)
      DO 220 I=START_POINT_NO_HALO-ROW_LENGTH+1,END_U_POINT_NO_HALO
        WORK2(I) = .5*(V(I)+V(I-1))*(PSTAR(I)-PSTAR(I+ROW_LENGTH))
     *             *LATITUDE_STEP_INVERSE*BK
 220  CONTINUE


      WORK2(END_U_POINT_NO_HALO)=WORK2(END_U_POINT_NO_HALO-1)
      WP(START_POINT_NO_HALO)=WP(START_POINT_NO_HALO+1)
! MPP Code : No need to do recalculations of end points because cyclic
! boundary conditions means that halos do this for us automatically



C CALCULATE U.GRAD P

! Loop over field, missing top and bottom rows and halos
      DO 224 I=START_POINT_NO_HALO+1,END_U_POINT_NO_HALO
        WP(I)=(WP(I)+.5*(WORK2(I)+WORK2(I-ROW_LENGTH)))
     *         *SEC_U_LATITUDE(I)
 224  CONTINUE

CL---------------------------------------------------------------------
CL    SECTION 3.     CALL CALC_TS TO GET TS.
CL---------------------------------------------------------------------

C STORE PRESSURE IN WORK.
! Loop over field, missing top and bottom rows and halos
      DO 300 I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
        WORK2(I) = AK + BK*PSTAR(I)
 300  CONTINUE

C CALCULATE NUMBER OF POINTS CALC_TS TO BE CALLED FOR.
      POINTS = END_U_POINT_NO_HALO-START_POINT_NO_HALO+1

C TS IS RETURNED IN WORK.

      CALL CALC_TS(WORK2(START_POINT_NO_HALO),TS(START_POINT_NO_HALO),
     &             POINTS,CONSTANT_PRESSURE,LLINTS)

CL---------------------------------------------------------------------
CL    SECTION 4.     CALCULATE WK AS IN EQUATION (43).
CL---------------------------------------------------------------------

! Loop over field, missing top and bottom rows and halos
      DO 400 I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
           OMEGA(I)= WP(I)+DP_BY_DT(I)+.5*(ETADOT_PLUS(I)*
     *          (DELTA_AKH_PLUS+DELTA_BKH_PLUS*PSTAR(I))+ETADOT_MINUS(I)
     *          *(DELTA_AKH_MINUS+DELTA_BKH_MINUS*PSTAR(I)))
           WK(I) = -R*TS(I)*OMEGA(I)/(G*WORK2(I))
 400  CONTINUE


CL    END OF ROUTINE V_CORIOL

      RETURN
      END

