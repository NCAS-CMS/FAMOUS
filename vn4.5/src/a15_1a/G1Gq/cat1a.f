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
CLL  SUBROUTINE CAT-----------------------------------------------------
CLL
CLL  PURPOSE:   Calculates clear air turbulence
CLL  Tested under compiler CFT77
CLL  Tested under OS version 5.1
CLL
CLL  Author J.T.Heming            Date: 22/02/91
CLL D.Forrester <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   4.2    Oct. 96  T3E migration: *DEF CRAY removed
CLL                                   S.J.Swarbrick
!LL   4.5    Apr. 98  Add start-end arguments to V_INT calls.
!LL                   Removed unnecessary SWAPBOUNDS calls.
!LL                   WARNING: CAT does not initialise NS MPP halos of
!LL                   CAT_PROB. S.D.Mullerworth
CLL   4.5   05/02/98  Alter CAT probability calculation.
CLL                    Clare Bysouth and Tim Westmacott.
CLL
CLL  Logical component number: D413
CLL
CLL  Project task:
CLL
CLL  Programming standard: U M DOC  Paper NO. 4,
CLL
CLL  External documentation:
CLL
CLLEND------------------------------------------------------------------
C
C*L  ARGUMENTS:---------------------------------------------------------
      SUBROUTINE CAT(
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
! MPP For relating local array to global array
     &      GLSIZE,
C data in
     & U,V,PUV,PSTAR,PRESS_REQD,MAX_WIND_P,
C data out
     & CAT_PROB,
C constants in
     & P_FIELD,U_FIELD,P_LEVELS,ROW_LENGTH,P_ROWS,SEC_U_LATITUDE,
     & AK,BK,EW_SPACE,NS_SPACE)
C*L
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      EXTERNAL  V_INT,P_TO_UV,H_GRAD,SPLINE,EVAL_SP
C*----------------------------------------------------------------------
      INTEGER
     * P_FIELD                       ! IN  NO OF POINTS ON P/T GRID
     *,U_FIELD                       ! IN  NO OF POINTS ON U/V GRID
     *,P_LEVELS                      ! IN  NO OF MODEL LEVELS
     *,P_ROWS                        ! IN  NO OF ROWS ON P/T GRID
     *,ROW_LENGTH                    ! IN  NO OF COLUMNS
     &,GLSIZE                        ! IN  GLOBAL NO OF ROWS

C-----------------------------------------------------------------------
      REAL
     * U(U_FIELD,P_LEVELS)           ! IN  U FIELD   AT FULL LEVELS
     *,V(U_FIELD,P_LEVELS)           ! IN  V FIELD   AT FULL LEVELS
     *,PRESS_REQD(U_FIELD)           ! IN/OUT  PRESSURE AT WHICH TO
     &                               ! CALCULATE CAT (modified by
     &                               ! routine)
     *,MAX_WIND_P(U_FIELD)           ! IN PRESSURE OF MAXIMUM WIND
     *,CAT_PROB(U_FIELD)             ! OUT CAT PROBABILITY
     *                               !
     *,PUV(U_FIELD,P_LEVELS)         ! IN  PRESS FIELD AT U/V POINTS
     *,PSTAR(P_FIELD)                ! IN  SURFACE PRESSURE FIELD
     *,SEC_U_LATITUDE(U_FIELD)       ! IN  1/COS(LAT) AT U/V POINTS
     *,AK(P_LEVELS)                  ! IN  A ARRAY AT FULL LEVELS
     *,BK(P_LEVELS)                  ! IN  B ARRAY AT FULL LEVELS
     *,EW_SPACE                      ! IN  LONGITUDE GRID SPACING
     *,NS_SPACE                      ! IN  LATITUDE GRID SPACING
     *                               !     BOTH THE ABOVE IN DEGREES
C*
C*L
C-----------------------------------------------------------------------
C Local Variables
C-----------------------------------------------------------------------
      INTEGER
     * I,K                        ! LOOP COUNTERS
     & ,GLOBAL_ROW_NO             ! To store global no of local row
     & ,S_HEM_TOP                 ! Global row no of first s.hem row
C-----------------------------------------------------------------------
      REAL
     * HORIZ_SHEAR(U_FIELD)       ! HORIZONTAL WIND SHEAR
     *,VERT_SHEAR(U_FIELD)        ! VERTICAL WIND SHEAR
     *,ETA(P_LEVELS)              ! ETA VALUES
     *,TEMP_REQD                  ! ICAO TEMPERATURE OF PRESS_REQD
     *,WORK1(U_FIELD,P_LEVELS)    ! 3-DIMENSIONAL WORKSPACE
     *,WORK2(U_FIELD,P_LEVELS)    ! 3-DIMENSIONAL WORKSPACE
     *,WORK3(U_FIELD,P_LEVELS)    ! 3-DIMENSIONAL WORKSPACE
     *,WORK4(U_FIELD)             ! 2-DIMENSIONAL WORKSPACE
     *,WORK5(U_FIELD)             ! 2-DIMENSIONAL WORKSPACE
     *,ALPHA                      ! WORK VARIABLE
     *,PLEV(17)                   ! ARRAY OF STANDARD PRESSURES
     *,TEMP(17)                   ! ARRAY OF ICAO TEMPERATURES
     *,EVAL(8)                    ! ARRAY OF CAT PREDICTION INDEX VALUES
     *,CATVAL(8)                  ! ARRAY OF CAT PROBABILITY VALUES
C-----------------------------------------------------------------------
      LOGICAL
     * FOUND(U_FIELD)             ! USED IN INTERPOLATIONS

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
C-----------------------------------------------------------------------
CL    Constants required     G=acceleration due to gravity
CL                           R=gas constant
CL                           PREF=reference pressure
C-----------------------------------------------------------------------
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

C-----------------------------------------------------------------------
C*L   DATA FOR LOCAL ARRAYS
C-----------------------------------------------------------------------
      DATA PLEV/100000.,95000.,85000.,70000.,50000.,40000.,30000.,
     + 25000.,20000.,15000.,10000.,7000.,5000.,3000.,2000.,1000.,999./
      DATA TEMP/287.,285.,279.,269.,252.,241.,229.,
     + 221.,217.,217.,217.,217.,217.,220.,223.,228.,230./
      DATA EVAL/5.0,7.5,10.0,15.0,20.0,25.0,65.0,99.0/
      DATA CATVAL/0.0,1.75,2.0,2.5,3.0,3.5,7.5,7.5/
C*----------------------------------------------------------------------
CL    Interpolate the U field onto the required level
C     Stored in WORK5. WORK1 and WORK3 are dummy arguments
C-----------------------------------------------------------------------
      CALL V_INT(PUV,PRESS_REQD,U,WORK5,
     &  U_FIELD,P_LEVELS,WORK1,WORK3,.FALSE.
     &  ,FIRST_VALID_PT,LAST_U_VALID_PT)
C-----------------------------------------------------------------------
CL    Interpolate the V field onto the required level
C     Stored in WORK4. WORK1 and WORK3 are dummy arguments
C-----------------------------------------------------------------------
      CALL V_INT(PUV,PRESS_REQD,V,WORK4,
     &  U_FIELD,P_LEVELS,WORK1,WORK3,.FALSE.
     &  ,FIRST_VALID_PT,LAST_U_VALID_PT)
C-----------------------------------------------------------------------
CL    For calculation of vertical shear - if PRESS_REQD < 30mb above or
CL    below MAX_WIND_P set PRESS_REQD to 30mb below or above its present
CL    value where the vertical shear has a more representative value.
C-----------------------------------------------------------------------
      DO I=FIRST_FLD_PT,LAST_U_FLD_PT
        IF((MAX_WIND_P(I)-PRESS_REQD(I)).LT.3000.0.AND.
     &     (MAX_WIND_P(I)-PRESS_REQD(I)).GE.0.0)THEN
          PRESS_REQD(I)=PRESS_REQD(I)-3000.0
        ELSEIF((PRESS_REQD(I)-MAX_WIND_P(I)).LT.3000.0.AND.
     &         (PRESS_REQD(I)-MAX_WIND_P(I)).GT.0.0)THEN
          PRESS_REQD(I)=PRESS_REQD(I)+3000.0
        ENDIF
      ENDDO
C=======================================================================
C
C=======================================================================
CL 1. CALCULATION OF HORIZONTAL WIND SHEAR AT REQUIRED P-LEVEL
C-----------------------------------------------------------------------
C
C               (    dU       dU       dV       dV)
C  Horizontal = (U*V*-- - U*U*-- + V*V*-- - U*V*--) / (U*U+V*V)
C  wind shear   (    dX       dY       dX       dY)
C
C-----------------------------------------------------------------------
CL    Calculation of dU/dX,dU/dY,dV/dX & dV/dY
C-----------------------------------------------------------------------
C     WORK1(I,1) holds dU/dX         WORK2(I,1) holds dU/dY
C     WORK1(I,2) holds dV/dX         WORK2(I,2) holds dV/dY
C-----------------------------------------------------------------------
      CALL H_GRAD(
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
     &     WORK5,U_FIELD,SEC_U_LATITUDE,ROW_LENGTH,EW_SPACE,NS_SPACE,
     &     WORK1(1,1),WORK2(1,1))
C
      CALL H_GRAD(
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
     &     WORK4,U_FIELD,SEC_U_LATITUDE,ROW_LENGTH,EW_SPACE,NS_SPACE,
     &     WORK1(1,2),WORK2(1,2))

C-----------------------------------------------------------------------
CL    Calculation of HORIZ_SHEAR in m/s per 100km (i.e. per s x10E-5)
C-----------------------------------------------------------------------

      DO I=FIRST_FLD_PT,LAST_U_FLD_PT
        HORIZ_SHEAR(I)=((WORK5(I)*WORK4(I))*WORK1(I,1)-
     &    (WORK5(I)*WORK5(I))*WORK2(I,1))+
     &    ((WORK4(I)*WORK4(I))*WORK1(I,2)-(WORK5(I)*WORK4(I))
     &    *WORK2(I,2))
        HORIZ_SHEAR(I)=(HORIZ_SHEAR(I)/(WORK5(I)*WORK5(I)
     &    +WORK4(I)*WORK4(I)))*100000.0
      ENDDO
C
CL    Change sign of horizontal shear in s.hem of global model
C
C=======================================================================
C
C=======================================================================
CL 2. CALCULATION OF VERTICAL WIND SHEAR AT REQUIRED P-LEVEL
C-----------------------------------------------------------------------
C
C               (dU   dU   dV   dV)       dP
C  Vertical   = (-- * -- + -- * --)**0.5 *--
C  wind shear   (dP   dP   dP   dP)       dH
C
C-----------------------------------------------------------------------
CL    Calculate true A value (AK/PREF)- stored in WORK3
C-----------------------------------------------------------------------
      DO K=1,P_LEVELS
        ALPHA=AK(K)/PREF
        DO I=1,U_FIELD
          WORK3(I,K)=ALPHA
        ENDDO
      ENDDO
C-----------------------------------------------------------------------
CL    Calculate ETA values
C-----------------------------------------------------------------------
      DO K=1,P_LEVELS
        ETA(K)=WORK3(1,K)+BK(K)
      ENDDO
C-----------------------------------------------------------------------
CL    Calculation of dAK/dETA at ETA levels using cubic spline
C-----------------------------------------------------------------------
C     WORK3 holds AK values,WORK1 holds dAK/dETA values
C     WORK4 & 5 are dummy arguments
C-----------------------------------------------------------------------
      CALL SPLINE(ETA,WORK3,U_FIELD,P_LEVELS,WORK2)
      DO K=1,P_LEVELS
        CALL EVAL_SP(ETA,WORK3,WORK2,U_FIELD,P_LEVELS,ETA(K),WORK4,
     *       WORK1(1,K),WORK5)
      ENDDO
C-----------------------------------------------------------------------
CL    Calculation of dP/dETA at ETA levels
C-----------------------------------------------------------------------
C     P=(AK*PREF)+(BK*PSTAR)        ETA=AK+BK
C =>  dP/dETA = (PREF*dAK/dETA) + PSTAR *(1-dAK/dETA)
C
C-----------------------------------------------------------------------
CL    Interpolate PSTAR onto U/V points - held in WORK5
C-----------------------------------------------------------------------
      CALL P_TO_UV(PSTAR,WORK5,P_FIELD,U_FIELD,ROW_LENGTH,P_ROWS)
C-----------------------------------------------------------------------
CL    dP/dETA held in WORK1
C-----*********************---------------------------------------------
      DO K=1,P_LEVELS
        DO I=FIRST_FLD_PT,LAST_U_FLD_PT
          WORK1(I,K)=PREF*WORK1(I,K)+WORK5(I)*(1.0-WORK1(I,K))
        ENDDO
      ENDDO
C=======================================================================
C
C=======================================================================
CL    Calculation of dU/dETA at ETA levels using a cubic spline
C-----------------------------------------------------------------------
C     WORK2 holds dU/dETA values. WORK4 & 5 are dummy arguments.
C-----------------------------------------------------------------------
      CALL SPLINE(ETA,U,U_FIELD,P_LEVELS,WORK3)
      DO K=1,P_LEVELS
        CALL EVAL_SP(ETA,U,WORK3,U_FIELD,P_LEVELS,ETA(K),WORK4,
     *       WORK2(1,K),WORK5)
      ENDDO
C-----------------------------------------------------------------------
CL    Calculation of dU/dP at ETA levels
C-----------------------------------------------------------------------
C     dU/dP = (dU/dETA) / (dP/dETA)
C-----------------------------------------------------------------------
CL    dU/dP at ETA levels held in WORK2
C-----*********************************---------------------------------
      DO K=1,P_LEVELS
        DO I=FIRST_FLD_PT,LAST_U_FLD_PT
          WORK2(I,K)=WORK2(I,K)/WORK1(I,K)
        ENDDO
      ENDDO
C-----------------------------------------------------------------------
CL    Calculation of dU/dP at required P-level - held in WORK5
C--------------------*****************************************----------
      CALL V_INT(PUV,PRESS_REQD,WORK2,WORK5,
     & U_FIELD,P_LEVELS,WORK1,WORK4,.FALSE.,FIRST_FLD_PT,LAST_U_FLD_PT)
C=======================================================================
C
C=======================================================================
CL    Calculation of dV/dETA at ETA levels using a cubic spline
C-----------------------------------------------------------------------
C     WORK2 holds dV/dETA values. WORK4 & VERT_SHEAR are dummy arguments
C-----------------------------------------------------------------------
      CALL SPLINE(ETA,V,U_FIELD,P_LEVELS,WORK3)
      DO K=1,P_LEVELS
        CALL EVAL_SP(ETA,V,WORK3,U_FIELD,P_LEVELS,ETA(K),WORK4,
     *       WORK2(1,K),VERT_SHEAR)
      ENDDO
C-----------------------------------------------------------------------
CL    Calculation of dV/dP at ETA levels
C-----------------------------------------------------------------------
C     dV/dP = (dV/dETA) / (dP/dETA)
C-----------------------------------------------------------------------
CL    dV/dP at ETA levels held in WORK2
C-----*********************************---------------------------------
      DO K=1,P_LEVELS
        DO I=FIRST_FLD_PT,LAST_U_FLD_PT
          WORK2(I,K)=WORK2(I,K)/WORK1(I,K)
        ENDDO
      ENDDO
C-----------------------------------------------------------------------
CL    Calculation of dV/dP at required p-level - held in WORK4
C--------------------*****************************************----------
      CALL V_INT(PUV,PRESS_REQD,WORK2,WORK4,
     & U_FIELD,P_LEVELS,WORK1,WORK3,.FALSE.,FIRST_FLD_PT,LAST_U_FLD_PT)
C=======================================================================
C
C=======================================================================
CL    Calculation of dP/dH at required P level using ICAO standard atmos
C-----------------------------------------------------------------------
CL    Reset pressures greater than 1000mb or less than 10mb
C-----------------------------------------------------------------------
      DO I=FIRST_FLD_PT,LAST_U_FLD_PT
        IF(PRESS_REQD(I).GT.100000.) PRESS_REQD(I)=100000.
        IF(PRESS_REQD(I).LT.1000.) PRESS_REQD(I)=1000.
        FOUND(I)=.FALSE.
      ENDDO
C-----------------------------------------------------------------------
CL    Find the first value of PLEV which is less than PRESS_REQD
C-----------------------------------------------------------------------
CL *** Following loop labelled to workaround fmp mistranslation
C-----------------------------------------------------------------------
      DO 50 K=1,16
        DO I=FIRST_FLD_PT,LAST_U_FLD_PT
          IF((PLEV(K+1).LT.PRESS_REQD(I)).AND.(.NOT.FOUND(I))) THEN
            FOUND(I)=.TRUE.
C-----------------------------------------------------------------------
CL    Calculate the ICAO temperature by interpolation
CL    ALPHA=Interpolation weight
C-----------------------------------------------------------------------
            ALPHA=ALOG(PRESS_REQD(I)/PLEV(K))
     *       /ALOG(PLEV(K+1)/PLEV(K))
            TEMP_REQD=ALPHA*TEMP(K+1)+(1.-ALPHA)*TEMP(K)
C-----------------------------------------------------------------------
CL    Calculation of  dP/dH
C     P=Rho*R*T             dP/dH=-Rho*g
C     => dP/dH=-g*P/(R*T)
CL    dP/dH - held in ALPHA
C-----*********************---------------------------------------------
CL  Sign of ALPHA changed
            ALPHA=+G*PRESS_REQD(I)/(R*TEMP_REQD)
C=======================================================================
C
C=======================================================================
CL    Calculation of VERT_SHEAR in m/s per km (i.e. per s x10E-3)
C-----------------------------------------------------------------------
            VERT_SHEAR(I)=((WORK5(I)*WORK5(I)+WORK4(I)*WORK4(I))**0.5)
     &      *ALPHA*1000.0
CL    Modify VERT_SHEAR if too large
              IF(VERT_SHEAR(I).GT.6.0) THEN
                VERT_SHEAR(I)=0.67*VERT_SHEAR(I)+1.0
              ELSEIF(VERT_SHEAR(I).GT.1.0) THEN
                VERT_SHEAR(I)=0.80*VERT_SHEAR(I)+0.2
              ENDIF
          ENDIF
        ENDDO
 50   CONTINUE
C=======================================================================
C
C=======================================================================
CL 3. CALCULATE CAT PREDICTOR INDEX E - held in WORK1
C-----------------------------------------------------------------------
      DO I=FIRST_FLD_PT,LAST_U_FLD_PT
        WORK1(I,1)=1.25*HORIZ_SHEAR(I)+0.25*VERT_SHEAR(I)*VERT_SHEAR(I)
     &      +10.5
        FOUND(I)=.FALSE.
      ENDDO
C-----------------------------------------------------------------------
CL    If E is less than 5.0 CAT_PROB is set to 0.0
CL    If E is greater than 65.0 CAT_PROB is set to 7.5
C-----------------------------------------------------------------------
      DO I=FIRST_FLD_PT,LAST_U_FLD_PT
        IF(WORK1(I,1).LE.5.0) THEN
          CAT_PROB(I)=0.0
          FOUND(I)=.TRUE.
        ELSEIF(WORK1(I,1).GE.65.0)THEN
          CAT_PROB(I)=7.5
          FOUND(I)=.TRUE.
        ENDIF
      ENDDO
C-----------------------------------------------------------------------
CL    Find the first value of EVAL which is greater than E
C-----------------------------------------------------------------------
CL *** Following loop labelled to workaround fmp mistranslation
C-----------------------------------------------------------------------
      DO 100 K=1,7
        DO I=FIRST_FLD_PT,LAST_U_FLD_PT
          IF((EVAL(K+1).GT.WORK1(I,1)).AND.(.NOT.FOUND(I))) THEN
            FOUND(I)=.TRUE.
C-----------------------------------------------------------------------
CL    Calculate the CAT probability by interpolation
CL    ALPHA=Interpolation weight
C-----------------------------------------------------------------------
            ALPHA=(WORK1(I,1)-EVAL(K))/(EVAL(K+1)-EVAL(K))
            CAT_PROB(I)=ALPHA*CATVAL(K+1)+(1.-ALPHA)*CATVAL(K)
          ENDIF
        ENDDO
 100  CONTINUE
C=======================================================================
C     END OF SUBROUTINE CAT
C=======================================================================
      RETURN
      END
C=======================================================================
