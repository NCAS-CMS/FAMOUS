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
!+ Interfaces GWAV_CTL to versions 1A & 2A of gravity wave drag scheme.
!
      SUBROUTINE GWAV_INTCTL
     1  (PSTAR,PEXNER,THETA,Q,U,V,P_FIELD,U_FIELD,
     2  ROWS_P,ROW_LENGTH,START_LEVEL,LEVELS,Q_LEVELS,
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
     3  AK,BK,AKH,BKH,DELTA_AK,DELTA_BK,SD_OROG_LAND,
     4  OROG_GRAD_XX_LAND,OROG_GRAD_XY_LAND,OROG_GRAD_YY_LAND,
     5  LAND_INDEX,LAND_POINTS,TIMESTEP,KAY,KAY_LEE,  
     6  STRESS_UD,LEN_STRESS_UD,STRESS_UD_ON,U_LIST1,POINTS_STRESS_UD,
     7  STRESS_VD,LEN_STRESS_VD,STRESS_VD_ON,V_LIST1,POINTS_STRESS_VD,
     8  DU_DT_SATN,LEN_DU_DT_SATN,DU_DT_SATN_ON,
     9  U_LIST2,POINTS_DU_DT_SATN,
     &  DV_DT_SATN,LEN_DV_DT_SATN,DV_DT_SATN_ON,
     &  V_LIST2,POINTS_DV_DT_SATN,
     &  DU_DT_JUMP,LEN_DU_DT_JUMP,DU_DT_JUMP_ON,
     &  U_LIST3,POINTS_DU_DT_JUMP,
     &  DV_DT_JUMP,LEN_DV_DT_JUMP,DV_DT_JUMP_ON,
     &  V_LIST3,POINTS_DV_DT_JUMP,
     &  DU_DT_LEE,LEN_DU_DT_LEE,DU_DT_LEE_ON,U_LIST4,POINTS_DU_DT_LEE,
     &  DV_DT_LEE,LEN_DV_DT_LEE,DV_DT_LEE_ON,V_LIST4,POINTS_DV_DT_LEE,
     &  TRANS_D,LEN_TRANS_D,TRANS_D_ON,POINTS_TRANS_D,IRET,
     &  LFROUDE,LGWLINP)

      IMPLICIT NONE
!
! Description:
! Is called with all possible arguments for all versions of gravity
! wave drag scheme, and according to *IF DEF (A06_1A or A06_2A in this
! case) calls corresponding actual G_WAVE routine with only those
! arguments it requires.
!
! Current Code Owner: R.T.H.Barnes.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  3.4   22/11/94   Original code. R.T.H.Barnes.
!LL  4.3  13/02/97  Stop main code being called if no land points
!LL                                                 P.Burton
!  4.3   22/04/97   KAY_LEE passed in from namelist. S.Webster
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
! System component covered:
! System Task:
!
! Global variables (*CALLed COMDECKs etc...):

! Subroutine arguments

      INTEGER
     &  P_FIELD            !IN    1ST DIMENSION OF FIELD OF PSTAR
     &, U_FIELD            !IN    1ST DIMENSION OF FIELD OF U,V
     &, ROWS_P             !IN    NUMBER OF ROWS of P grid
     &, ROW_LENGTH         !IN    NUMBER OF POINTS PER ROW
     &, START_LEVEL        !IN    START OF WAVE-BREAKING TEST
     &, LEVELS             !IN    NUMBER OF MODEL LEVELS
     &, Q_LEVELS           !IN    NUMBER OF WET LEVELS
     &, LAND_POINTS        !IN    NUMBER OF LAND POINTS
     &, LAND_INDEX((ROWS_P)*ROW_LENGTH) ! INDEX FOR LAND POINTS
     &, IRET               ! RETURN CODE      :    IRET=0   NORMAL EXIT
!                          ! RETURN CODE      :    IRET=1   ?????
     &, LEN_STRESS_UD      !IN    ) Dimension of diagnostic arrays
     &, LEN_STRESS_VD      !IN    ) for GW stress - u and v
     &, POINTS_STRESS_UD   !IN    ) No of land points in diagnostic
     &, POINTS_STRESS_VD   !IN    ) arrays for GW stress - u and v
     &, LEN_DU_DT_SATN     !IN    ) Dimension of diagnostic arrays
     &, LEN_DV_DT_SATN     !IN    ) for GW satn - du and dv
     &, POINTS_DU_DT_SATN  !IN    ) No of land points in diagnostic
     &, POINTS_DV_DT_SATN  !IN    ) arrays for GW satn - du and dv
     &, LEN_DU_DT_JUMP     !IN    ) Dimension of diagnostic arrays
     &, LEN_DV_DT_JUMP     !IN    ) for GW satn - du and dv
     &, POINTS_DU_DT_JUMP  !IN    ) No of land points in diagnostic
     &, POINTS_DV_DT_JUMP  !IN    ) arrays for GW satn - du and dv
     &, LEN_DU_DT_LEE      !IN    ) Dimension of diagnostic arrays
     &, LEN_DV_DT_LEE      !IN    ) for GW lee - du and dv
     &, POINTS_DU_DT_LEE   !IN    ) No of land points in diagnostic
     &, POINTS_DV_DT_LEE   !IN    ) arrays for GW lee - du and dv
     &, LEN_TRANS_D        !IN    Dimension of diag for trans. coeff.
     &, POINTS_TRANS_D     !IN    No. of land points for trans. coeff.
C
C
! All TYPFLDPT variables are Intent IN
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
     & PSTAR(P_FIELD)         !IN    PRIMARY MODEL ARRAY FOR PSTAR FIELD
     &,PEXNER(P_FIELD,LEVELS+1) !IN    ARRAY FOR EXNER PRESSURE FIELD
     &,THETA(P_FIELD,LEVELS)  !IN    PRIMARY MODEL ARRAY FOR THETA FIELD
     &,Q(P_FIELD,Q_LEVELS)    !IN    SPECIFIC HUMIDITY AT FULL LEVELS
     &,U(U_FIELD,LEVELS)      !INOUT PRIMARY MODEL ARRAY FOR U FIELD
     &,V(U_FIELD,LEVELS)      !INOUT PRIMARY MODEL ARRAY FOR V FIELD
C            AK,BK  DEFINE HYBRID VERTICAL COORDINATES P=A+BP*,
C       DELTA_AK,DELTA_BK  DEFINE LAYER PRESSURE THICKNESS PD=AD+BDP*,

      REAL
     & DELTA_AK(LEVELS)       !IN    LAYER THICKNESS
     &,DELTA_BK(LEVELS)       !IN    LAYER THICKNESS
     &,AK (LEVELS)            !IN    VALUE AT LAYER CENTRE
     &,BK (LEVELS)            !IN    VALUE AT LAYER CENTRE
     &,AKH(LEVELS+1)          !IN    VALUE AT LAYER BOUNDARY
     &,BKH(LEVELS+1)          !IN    VALUE AT LAYER BOUNDARY
     &,SD_OROG_LAND(LAND_POINTS)  !IN STANDARD DEVIATION OF OROGRAPHY
     &,OROG_GRAD_XX_LAND(LAND_POINTS)
!                             !IN    DH/DX SQUARED GRADIENT OROGRAPHY
     &,OROG_GRAD_XY_LAND(LAND_POINTS)
!                             !IN   (DH/DX)(DH/DY) GRADIENT OROGRAPHY
     &,OROG_GRAD_YY_LAND(LAND_POINTS)
!                             !IN    DH/DY SQUARED GRADIENT OROGRAPHY
     &,TIMESTEP               !IN    TIMESTEP
     &,KAY                    !IN    surface stress constant ( m-1)
     &,KAY_LEE                !IN    TRAPPED LEE WAVE CONSTANT          

     &,STRESS_UD(LEN_STRESS_UD,*)    !U STRESS DIAGNOSTIC
     &,STRESS_VD(LEN_STRESS_VD,*)    !V STRESS DIAGNOSTIC
     &,DU_DT_SATN(LEN_DU_DT_SATN,*)  !U ACCELN DIAGNOSTIC  (SATURATION)
     &,DV_DT_SATN(LEN_DV_DT_SATN,*)  !V ACCELN DIAGNOSTIC  (SATURATION)
     &,DU_DT_JUMP(LEN_DU_DT_JUMP,*)  !U ACCELN DIAG  (HYDRAULIC JUMP)
     &,DV_DT_JUMP(LEN_DV_DT_JUMP,*)  !V ACCELN DIAG  (HYDRAULIC JUMP)
     &,DU_DT_LEE(LEN_DU_DT_LEE,*)    !U ACCELN DIAG  (TRAPPED LEE WAVE)
     &,DV_DT_LEE(LEN_DV_DT_LEE,*)    !V ACCELN DIAG  (TRAPPED LEE WAVE)
     &,TRANS_D(LEN_TRANS_D)          !TRANSMITTION COEFF DIAGN

C WARNING: Storage will only be assigned by the calling routine for
C          for the number of levels required.

      LOGICAL
     & STRESS_UD_ON           !U stress diagnostic switch
     &,STRESS_VD_ON           !V stress diagnostic switch
     &,U_LIST1(LEVELS+1)      ! Lists of levels for which stresses
     &,V_LIST1(LEVELS+1)      ! required.
     &,DU_DT_SATN_ON          !U accel (saturation) diagnostic switch
     &,DV_DT_SATN_ON          !V accel (saturation) diagnostic switch
     &,U_LIST2(LEVELS)        ! Lists of levels for which accelerations
     &,V_LIST2(LEVELS)        ! required.
     &,DU_DT_JUMP_ON          !U accel (hydr jump) diagnostic switch
     &,DV_DT_JUMP_ON          !V accel (hydr jump) diagnostic switch
     &,U_LIST3(LEVELS)        ! Lists of levels for which accelerations
     &,V_LIST3(LEVELS)        ! required.
     &,DU_DT_LEE_ON           !U accel (lee wave) diagnostic switch
     &,DV_DT_LEE_ON           !V accel (lee wave) diagnostic switch
     &,U_LIST4(LEVELS)        ! Lists of levels for which accelerations
     &,V_LIST4(LEVELS)        ! required.
     &,TRANS_D_ON             !Transmittion coefficient diag switch
     &,LFROUDE                !Switch for Froude number limiting
     &,LGWLINP                !Switch for linear stress profile

! Local dynamic arrays:

! Function & Subroutine calls:
      EXTERNAL G_WAVE

! Local parameters:

! Local scalars:

C-------------------------------------------------------------------
CL    1.   Call Gravity Wave Drag Scheme versions 1A or 2A.
C-------------------------------------------------------------------

      IF (LAND_POINTS .NE. 0) THEN
      CALL G_WAVE
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
     5   STRESS_UD,LEN_STRESS_UD,STRESS_UD_ON,U_LIST1,POINTS_STRESS_UD,
     6   STRESS_VD,LEN_STRESS_VD,STRESS_VD_ON,V_LIST1,POINTS_STRESS_VD,
     7   IRET,LFROUDE,LGWLINP)
      ENDIF ! if land_points .ne. 0

      RETURN
      END

