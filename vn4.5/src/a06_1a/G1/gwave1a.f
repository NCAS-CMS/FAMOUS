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
CLL  SUBROUTINE G_WAVE-------------------------------------------
CLL
CLL  PURPOSE:   1) INTERPOLATE WINDS TO P/THETA POINTS
CLL             2) GATHER DATA FOR LAND POINTS ONLY
CLL             3) CALL SURFACE STRESS ROUTINE
CLL             4) CALL VERTICAL STRESS PROFILE ROUTINE TO CALCULATE
CLL                DRAG AT EACH LEVEL
CLL             5) INTERPOLATE ACCELERATION TO WIND POINTS AND UPDATE
CLL                WINDS
CLL  SUITABLE FOR SINGLE COLUMN USE, WITH CALLS TO: UV_TO_P REMOVED
CLL                                                 P_TO_UV REMOVED
!                                                (SCMA on)  
CLL  SUITABLE FOR ROTATED GRIDS
CLL
CLL  ORIGINAL VERSION FOR CRAY Y-MP
CLL  WRITTEN BY C. WILSON
CLL  FURTHER ALTERATIONS MAY BE REQUIRED FOR AUTOTASKING EFFICIENCY
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.3   25/10/93  Removal of DIAG06 directive. New arguments to
CLL                   dimension diagnostic arrays. D. Robinson.
CLL   3.4   11/05/94  Argument LFROUDE added and passed to GW_SURF
CLL                   DEF GWLINP replaced by LOGICAL LGWLINP
CLL                                                S.J.Swarbrick
!LL   4.1   31/05/96  Added MPP code    P.Burton
CLL
CLL   4.4   19/09/97  Remove *IF -DEF,CRAY compile options. S.Webster 
CLL  4.5    Jul. 98  Kill the IBM specific lines.
CLL                  Replace IBM with SCMA  (JCThil)
CLL                                                                     
CLL   4.5   17/03/97   Correct MPP GWD diagnostic bug. S.Webster
CLL
CLL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL  VERSION 1, DATED 12/09/89
CLL
CLL  SYSTEM TASK: CONTROL PART OF P22
CLL
CLL  DOCUMENTATION:
CLL
CLLEND-------------------------------------------------------------

C
C*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE G_WAVE
     1  (PSTAR,PEXNER,THETA,U,V,P_FIELD,U_FIELD,
     2   ROWS_P,ROW_LENGTH,START_LEVEL,LEVELS,
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
     3   AK,BK,AKH,BKH,DELTA_AK,DELTA_BK,SD_OROG_LAND,
     4   LAND_INDEX,LAND_POINTS, TIMESTEP,KAY,
     5   STRESS_UD,LEN_STRESS_UD,STRESS_UD_ON,U_LIST,LAND_POINTS_UD,
     6   STRESS_VD,LEN_STRESS_VD,STRESS_VD_ON,V_LIST,LAND_POINTS_VD,
     7  IRET,LFROUDE,LGWLINP)

      IMPLICIT NONE

      INTEGER
     *  P_FIELD            !IN    1ST DIMENSION OF FIELD OF PSTAR
     *, U_FIELD            !IN    1ST DIMENSION OF FIELD OF U,V
     *, ROWS_P             !IN    NUMBER OF ROWS of P grid
     *, ROW_LENGTH         !IN    NUMBER OF POINTS PER ROW
     *, START_LEVEL        !IN    START OF WAVE-BREAKING TEST
     *, LEVELS             !IN    NUMBER OF MODEL LEVELS
     *, LAND_POINTS        !IN    NUMBER OF LAND POINTS
     *, LAND_INDEX((ROWS_P)*ROW_LENGTH) ! INDEX FOR LAND POINTS
     *, LEN_STRESS_UD      !IN    ) Dimension of diagnostic arrays
     *, LEN_STRESS_VD      !IN    ) for GW stress - u and v
     *, LAND_POINTS_UD     !IN    ) No of land points in diagnostic
     *, LAND_POINTS_VD     !IN    ) arrays for GW stress - u and v
     *, IRET               ! RETURN CODE      :    IRET=0   NORMAL EXIT
C                          ! RETURN CODE      :    IRET=1   ?????

! All TYPFLDPT variables are intent IN
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
     * PSTAR(P_FIELD)         !IN    PRIMARY MODEL ARRAY FOR PSTAR FIELD
     *,PEXNER(P_FIELD,LEVELS+1) !IN    ARRAY FOR EXNER PRESSURE FIELD
     *,THETA(P_FIELD,LEVELS)  !IN    PRIMARY MODEL ARRAY FOR THETA FIELD
     *,U(U_FIELD,LEVELS)      !INOUT PRIMARY MODEL ARRAY FOR U FIELD
     *,V(U_FIELD,LEVELS)      !INOUT PRIMARY MODEL ARRAY FOR V FIELD
C            AK,BK  DEFINE HYBRID VERTICAL COORDINATES P=A+BP*,
C       DELTA_AK,DELTA_BK  DEFINE LAYER PRESSURE THICKNESS PD=AD+BDP*,

      REAL
     * DELTA_AK(LEVELS)       !IN    LAYER THICKNESS
     *,DELTA_BK(LEVELS)       !IN    LAYER THICKNESS
     *,AK (LEVELS)            !IN    VALUE AT LAYER CENTRE
     *,BK (LEVELS)            !IN    VALUE AT LAYER CENTRE
     *,AKH(LEVELS+1)          !IN    VALUE AT LAYER BOUNDARY
     *,BKH(LEVELS+1)          !IN    VALUE AT LAYER BOUNDARY
     *,SD_OROG_LAND(LAND_POINTS),  !IN STANDARD DEVIATION OF OROGRAPHY
     * TIMESTEP               !IN    TIMESTEP
     *,KAY                    !IN    surface stress constant ( m-1)
     *,STRESS_UD(LEN_STRESS_UD,*) !U STRESS DIAGNOSTIC
     *,STRESS_VD(LEN_STRESS_VD,*) !V STRESS DIAGNOSTIC

C WARNING: Storage will only be assigned by the calling routine for
C          for the number of levels required.

      LOGICAL
     * STRESS_UD_ON           !U stress diagnostic switch
     *,STRESS_VD_ON           !V stress diagnostic switch
     *,U_LIST(LEVELS+1),      ! Lists of levels for which stresses
     * V_LIST(LEVELS+1)       ! required.
     *,LFROUDE,LGWLINP        ! Logical switches
C*---------------------------------------------------------------------

C*L  WORKSPACE USAGE:-------------------------------------------------
C   DEFINE LOCAL WORKSPACE ARRAYS:
C   4 REAL ARRAYS AT FULL FIELD LENGTH REQUIRED
C   6*LEVELS+5 REAL ARRAYS OF LAND_POINTS LENGTH REQUIRED

C   2*LEVELS REAL ARRAYS OF LAND_POINTS LENGTH REQUIRED FOR DIAGNOSTICS


      REAL
     * WORK(P_FIELD,4)               ! GENERAL PURPOSE WORK
     *,UP_LAND(LAND_POINTS,LEVELS)   ! INTERPOLATED U COMPONENT ON PGRID
     *,VP_LAND(LAND_POINTS,LEVELS)   ! INTERPOLATED U COMPONENT ON PGRID
     *,THETA_LAND(LAND_POINTS,LEVELS)! THETA LAND POINTS
     *,S_STRESS(LAND_POINTS)         ! 'SURFACE' STRESS LAND POINTS
     *,SIN_A(LAND_POINTS)        ! SIN ('SURFACE' WIND ANGLE FROM NORTH)
     *,COS_A(LAND_POINTS)        ! COS ('SURFACE' WIND ANGLE FROM NORTH)
     *,PSTAR_LAND(LAND_POINTS)       ! PSTAR LAND POINTS
     *,PEXNER_LAND(LAND_POINTS,LEVELS+1)  ! PEXNER LAND POINTS
     *,DU_DT(LAND_POINTS,LEVELS)     ! U-ACCELERATION
     *,DV_DT(LAND_POINTS,LEVELS)     ! V-ACCELERATION

      REAL
     * STRESS_UD_LAND(LAND_POINTS_UD,LEVELS+1) !U STRESS DIAGNOSTIC
     *,STRESS_VD_LAND(LAND_POINTS_VD,LEVELS+1) !V STRESS DIAGNOSTIC


C*---------------------------------------------------------------------
C*L EXTERNAL SUBROUTINES CALLED---------------------------------------
      EXTERNAL GW_SURF,GW_LIN_P,GW_RICH         
     &,P_TO_UV,UV_TO_P
C*------------------------------------------------------------------
CL  MAXIMUM VECTOR LENGTH ASSUMED IS (ROWS_P+1) * ROWLENGTH
CL---------------------------------------------------------------------
C----------------------------------------------------------------------
C    DEFINE LOCAL VARIABLES
      INTEGER
     *  P_POINTS      !     NUMBER OF P POINTS NEEDED
     *, U_POINTS_1    !     No. U points used to interpolate to P-grid
     *, U_POINTS      !     NUMBER OF U POINTS UPDATED
     *, START_U       !     Start position of U points updated
     *, START_U1      !     Start position of diagnostics updated
C
      INTEGER   I,IW,K,     ! LOOP COUNTERS IN ROUTINE
     *          KOUT_U,KOUT_V

C-------------------------------------------------------------------
CL    INTERNAL STRUCTURE INCLUDING SUBROUTINE CALLS:
CL    1.     INITIALISATION
C--------------------------

C------------------------------------------------------------------
CL    1.1  SET UP DIMENSIONS
C------------------------------------------------------------------

        U_POINTS_1    = (ROWS_P+1)*ROW_LENGTH
        U_POINTS      = (ROWS_P-1)*ROW_LENGTH
        P_POINTS      =  ROWS_P*ROW_LENGTH
        START_U       =  ROW_LENGTH
C
C  Three separate cases for START_U1. These arise because of different
C  row offsets in the call to uv_to_p  (GWAVE3A.253,256)
C
        IF ( at_top_of_lpg ) THEN
          START_U1       = 2*ROW_LENGTH
        ELSE
          START_U1       = 0
        ENDIF

C------------------------------------------------------------------
CL    1.2 INTERPOLATE WINDS TO P/THETA-GRID
C------------------------------------------------------------------
      DO K=1,LEVELS


        CALL UV_TO_P(U(1,K),WORK(1,1),U_POINTS_1,P_POINTS,
     *   ROW_LENGTH,ROWS_P+1)
        CALL UV_TO_P(V(1,K),WORK(1,2),U_POINTS_1,P_POINTS,
     *   ROW_LENGTH,ROWS_P+1)

! Correct halos of interpolated U/V
!        CALL SWAPBOUNDS(WORK(1,1),LOCAL_ROW_LENGTH,ROWS_P+1,
!     &                  EW_Halo,NS_Halo,1) ! U field
!        CALL SWAPBOUNDS(WORK(1,2),LOCAL_ROW_LENGTH,ROWS_P+1,
!     &                  EW_Halo,NS_Halo,1) ! V field

C------------------------------------------------------------------
CL    1.3  GATHER WINDS AT LAND POINTS
C------------------------------------------------------------------

        DO I=1,LAND_POINTS
         UP_LAND(I,K) =WORK(LAND_INDEX(I),1)
         VP_LAND(I,K) =WORK(LAND_INDEX(I),2)
        END DO


      END DO

C------------------------------------------------------------------
CL    1.4  GATHER PSTAR,PEXNER,THETA,SD_OROG AT LAND POINTS
C------------------------------------------------------------------

      DO I=1,LAND_POINTS
        PSTAR_LAND(I) = PSTAR(LAND_INDEX(I))
      END DO

CL *** Following loop labelled to workaround fmp mistranslation

CFPP$ SELECT(CONCUR)
      DO 140 K=1,LEVELS
        DO I=1,LAND_POINTS
          PEXNER_LAND(I,K) = PEXNER(LAND_INDEX(I),K)
          THETA_LAND(I,K) = THETA(LAND_INDEX(I),K)
        END DO
 140  CONTINUE

      DO I=1,LAND_POINTS
        PEXNER_LAND(I,LEVELS+1) =PEXNER(LAND_INDEX(I),LEVELS+1)
      END DO

C------------------------------------------------------------------
CL    2. CALCULATE 'SURFACE' STRESS,CALL GW_SURF
C------------------------------------------------------------------

      CALL GW_SURF(PSTAR_LAND,PEXNER_LAND,THETA_LAND,UP_LAND,VP_LAND,
     *             SD_OROG_LAND,S_STRESS,LEVELS,LAND_POINTS,
     *             AK,BK,AKH,BKH,KAY,SIN_A,COS_A,LFROUDE)

C------------------------------------------------------------------
CL    3. CALCULATE STRESS PROFILE AND ACCELERATIONS,
CL       CALL GW_RICH  OR GW_LIN_P
C------------------------------------------------------------------

      IF (LGWLINP) THEN
      CALL GW_LIN_P(PSTAR_LAND,PEXNER_LAND,THETA_LAND,UP_LAND,VP_LAND,
     1             S_STRESS,START_LEVEL,LEVELS,LAND_POINTS,
     2             AKH,BKH,DELTA_AK,DELTA_BK,SIN_A,COS_A,
     3             DU_DT,DV_DT,
     4             STRESS_UD_LAND,LAND_POINTS_UD,STRESS_UD_ON,
     5             STRESS_VD_LAND,LAND_POINTS_VD,STRESS_VD_ON)
      ELSE
      CALL GW_RICH(PSTAR_LAND,PEXNER_LAND,THETA_LAND,UP_LAND,VP_LAND,
     1             S_STRESS,START_LEVEL,LEVELS,LAND_POINTS,
     2             AKH,BKH,DELTA_AK,DELTA_BK,KAY,SIN_A,COS_A,
     3             DU_DT,DV_DT,
     4             STRESS_UD_LAND,LAND_POINTS_UD,STRESS_UD_ON,
     5             STRESS_VD_LAND,LAND_POINTS_VD,STRESS_VD_ON)
      END IF

C------------------------------------------------------------------
CL    4. SCATTER ACCELERATIONS TO FULL AREA, INTERPOLATE TO UV-GRID
CL       AND UPDATE WINDS
C------------------------------------------------------------------

      DO I=1,P_FIELD
       DO IW=1,4
        WORK(I,IW) = 0.0
       END DO
      END DO

      DO K=1,LEVELS

! Fujitsu vectorization directive
!OCL NOVREC
CDIR$ IVDEP
        DO I=1,LAND_POINTS
          WORK(LAND_INDEX(I),1)= DU_DT(I,K)
          WORK(LAND_INDEX(I),2)= DV_DT(I,K)
        END DO

        CALL P_TO_UV(WORK(1,1),WORK(1,3),P_POINTS,U_POINTS,
     *               ROW_LENGTH,ROWS_P)
        CALL P_TO_UV(WORK(1,2),WORK(1,4),P_POINTS,U_POINTS,
     *               ROW_LENGTH,ROWS_P)

        DO I=1,U_POINTS
          U(START_U+I,K) = U(START_U+I,K) + TIMESTEP*WORK(I,3)
          V(START_U+I,K) = V(START_U+I,K) + TIMESTEP*WORK(I,4)
        END DO

      END DO

      IF (STRESS_UD_ON .OR. STRESS_VD_ON) THEN

        KOUT_U=0
        KOUT_V=0

        DO K=START_LEVEL,LEVELS+1

           IF (STRESS_UD_ON) THEN
             IF (U_LIST(K)) THEN
               KOUT_U=KOUT_U+1
! Fujitsu vectorization directive
!OCL NOVREC
CDIR$ IVDEP
               DO I=1,LAND_POINTS
                 WORK(LAND_INDEX(I),1)=STRESS_UD_LAND(I,K)
               END DO
               CALL P_TO_UV (WORK(1,1),STRESS_UD(START_U1+1,KOUT_U),
     *                       P_POINTS,U_POINTS,ROW_LENGTH,ROWS_P)
             ENDIF
           ENDIF

           IF (STRESS_VD_ON) THEN
             IF (V_LIST(K)) THEN
               KOUT_V=KOUT_V+1
! Fujitsu vectorization directive
!OCL NOVREC
CDIR$ IVDEP
               DO I=1,LAND_POINTS
                 WORK(LAND_INDEX(I),2)=STRESS_VD_LAND(I,K)
               END DO
               CALL P_TO_UV (WORK(1,2),STRESS_VD(START_U1+1,KOUT_V),
     *                       P_POINTS,U_POINTS,ROW_LENGTH,ROWS_P)
             ENDIF
           ENDIF

        ENDDO

      ENDIF

      IRET=0

      RETURN
      END

