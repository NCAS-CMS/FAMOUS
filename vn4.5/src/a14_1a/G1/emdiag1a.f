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
CLL  SUBROUTINE ENG_MASS_DIAG------------------------------------------
CLL
CLL  PURPOSE : PART OF ENERGY CORRECTION SUITE OF ROUTINES
CLL            - TO GLOBALLY INTERGATE TOTAL ENERGY AMD MASS OF
CLL              THE ATMOSPHERE
CLL
CLL  NOT SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY D.GREGORY FEBRUARY 1991
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.4   26/05/94  Argument LLINTS added and passed to CALC_RS
CLL                   DEF NOWHBR replaced by LOGICAL LWHITBROM
CLL                                                  S.J.Swarbrick
!     4.1   24/11/95  Changed interface to ENERGY/MASS_SUM to make
!                     suitable for MPP use and added TYPFLDPT
!                     arguments.                          P.Burton
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  SYSTEM TASK : P##
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P###
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE ENG_MASS_DIAG (TL,U,V,AREA_P,AREA_UV,P_FIELD,
     &                          U_FIELD,ROW_LENGTH,ROWS,
     2                          DELTA_AK,DELTA_BK,AK,BK,TOT_ENERGY,
     3                          TOT_MASS_P,PART_MASS_P,P_LEVELS,PSTAR,
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
     4                          LLINTS,LWHITBROM)
C
      IMPLICIT NONE
      LOGICAL  LLINTS,LWHITBROM
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

C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
C*----------------------------------------------------------------------

C
C----------------------------------------------------------------------
C VECTOR LENGTHS
C----------------------------------------------------------------------
C
C
      INTEGER P_FIELD          ! IN VECTOR LENGTH OF VARIABLES ON
                               !    P GRID
C
      INTEGER U_FIELD          ! IN VECTOR LENGTH OF VARIABLES ON
                               !    UV GRID
C
C
      INTEGER ROW_LENGTH       ! IN NUMBER OF POINTS PER ROW
C
      INTEGER ROWS             ! IN NUMBER OF ROWS IN P GRID
C
      INTEGER P_LEVELS         ! IN NUMBER OF LEVELS IN VERTICAL

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

C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL TL(P_FIELD,P_LEVELS)         !IN TEMPERATURE
C
      REAL U(U_FIELD,P_LEVELS)          !IN COMPONENT OF WIND
C
      REAL V(U_FIELD,P_LEVELS)          !IN COMPONENT OF WIND
C
      REAL AREA_P(P_FIELD)              !IN AREA OF CELLS IN P GRID
C
      REAL AREA_UV(U_FIELD)             !IN AREA OF CELLS IN UV GRID
C
      REAL DELTA_AK(P_LEVELS)           ! IN |THICKNESS OF LAYERS IN
C                                            |
      REAL DELTA_BK(P_LEVELS)           ! IN |ETA CO-ORDINATES
C
      REAL AK(P_LEVELS)                 ! IN |ETA CO-ORDINATES OF
C                                            |
      REAL BK(P_LEVELS)                 ! IN |MID-LAYER POINTS
C
      REAL PSTAR(P_FIELD)               !IN PRESSURE AT SURFACE
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE IN AND OUT
C----------------------------------------------------------------------
C
      REAL TOT_ENERGY             !   TOTAL ENERGY OF ATMOSPHERE
C
      REAL TOT_MASS_P             !   TOTAL MASS OF ATMOSPHERE
C
      REAL PART_MASS_P            !   PARTIAL MASS OF ATMOSPHERE
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C----------------------------------------------------------------------
C
      REAL PSTAR_DELBK(P_FIELD)    ! PRESSURE_AT_SURFACE*DELTA_BK
C
      REAL DELP_P(P_FIELD)         ! MASS ELEMENTS ON P GRID
C
      REAL DELP_UV(U_FIELD)        ! MASS ELEMENTS ON UV GRID
C
      REAL RS_P_K(P_FIELD)         ! RADII ON P GRID
C
      REAL RS_UV_K(U_FIELD)        ! RADII ON UV GRID
C
      REAL WORK(P_FIELD)           ! DUMMY VARIABLE
C
      REAL TS(P_FIELD)             ! OUTPUT FROM SUBROUTINE CALC_RS
C
C
C----------------------------------------------------------------------
C INTERNAL LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER I                ! LOOP COUNTER
C
      INTEGER K                ! LOOP COUNTER
      INTEGER POINTS  ! Number of points for CALC_RS to process
C
C----------------------------------------------------------------------
C EXTERNAL SUBROUTINE CALLS  -  P_TO_UV,CALC_RS,ENERGY_SUM,MASS_SUM
C----------------------------------------------------------------------
C
C*---------------------------------------------------------------------
C
      POINTS=LAST_P_VALID_PT-FIRST_VALID_PT+1
! Number of points to be processed by CALC_RS. For non-MPP runs this
! is simply P_FIELD, for MPP, it is all the points, minus any
! unused halo areas (ie. the halo above North pole row, and beneath
! South pole row)
C----------------------------------------------------------------------
C ZERO MASS OF ATMOSPHERE
C----------------------------------------------------------------------
C
      TOT_MASS_P = 0.0
      PART_MASS_P = 0.0
! QAN fix
! Zero DELP_P and RS_P_Karray
      DO I=1,P_FIELD
        DELP_P(I)=0.0
        RS_P_K(I)=0.0
      ENDDO
C
C----------------------------------------------------------------------
C ZERO ENERGY OF ATMOSPHERE
C----------------------------------------------------------------------
C
      TOT_ENERGY = 0.0
C
C======================================================================
C MAIN LOOP OVER VERTICAL LEVELS
C======================================================================
C
      DO K=1,P_LEVELS
C
C----------------------------------------------------------------------
C CALCULATE MASS OF LEVEL K AT EACH GRID POINT AND ALSO
C P*DELTA_BK AT EACH GRID POINT
C----------------------------------------------------------------------
C
C
! Loop over all points, including halos
       DO I=FIRST_VALID_PT,LAST_P_VALID_PT
        PSTAR_DELBK(I) = -DELTA_BK(K)*PSTAR(I)
        DELP_P(I) = -DELTA_AK(K) + PSTAR_DELBK(I)
       END DO
C
C----------------------------------------------------------------------
C INTERPOLATE DELP_P TO UV GRID
C----------------------------------------------------------------------
C
       CALL P_TO_UV (DELP_P,DELP_UV,P_FIELD,U_FIELD,ROW_LENGTH,ROWS)
C
C----------------------------------------------------------------------
C CALCULATE RADIUS OF SPHERE AT LEVEL K
C----------------------------------------------------------------------
C
      IF (.NOT.LWHITBROM) THEN
C
       DO I=1,P_FIELD
        RS_P_K(I) = A
       END DO
C
       DO I=1,U_FIELD
        RS_UV_K(I) = A
       END DO
C
      ELSE
C
       CALL CALC_RS(PSTAR(FIRST_VALID_PT),AK,BK,TS(FIRST_VALID_PT),
     &              WORK(FIRST_VALID_PT),RS_P_K(FIRST_VALID_PT),
     &              POINTS,K,P_LEVELS,LLINTS)
C
C
C----------------------------------------------------------------------
C INTERPLOATE RADIUS OF SPHERE AT LEVEL K TO UV GRID
C----------------------------------------------------------------------
C
       CALL P_TO_UV (RS_P_K,RS_UV_K,P_FIELD,U_FIELD,ROW_LENGTH,ROWS)
C
      END IF
C
C----------------------------------------------------------------------
C SUM CP*TL OVER GLOBE FOR LEVEL K AND ADD TO TOTAL ENERGY SUM
C----------------------------------------------------------------------
C
       CALL ENERGY_SUM (TL(1,K),START_POINT_NO_HALO,
     &                  END_P_POINT_NO_HALO,P_FIELD,
     &                  AREA_P,DELP_P,RS_P_K,CP,TOT_ENERGY)
C
C
C----------------------------------------------------------------------
C SUM 0.5*U*U OVER GLOBE FOR LEVEL K AND ADD TO TOTAL ENERGY SUM
C----------------------------------------------------------------------
C
! Loop over all points except North and South Halos
       DO I=FIRST_FLD_PT,LAST_U_FLD_PT
        WORK(I) = U(I,K)*U(I,K)
       END DO
C
       CALL ENERGY_SUM (WORK,START_POINT_NO_HALO,
     &                  END_U_POINT_NO_HALO,U_FIELD,
     &                  AREA_UV,DELP_UV,RS_UV_K,0.5,TOT_ENERGY)
C
C
C----------------------------------------------------------------------
C SUM 0.5*V*V OVER GLOBE FOR LEVEL K AND ADD TO TOTAL ENERGY SUM
C----------------------------------------------------------------------
C
! Loop over all points except North and South Halos
       DO I=FIRST_FLD_PT,LAST_U_FLD_PT
        WORK(I) = V(I,K)*V(I,K)
       END DO
C
       CALL ENERGY_SUM (WORK,START_POINT_NO_HALO,
     &                  END_U_POINT_NO_HALO,U_FIELD,
     &                  AREA_UV,DELP_UV,RS_UV_K,0.5,TOT_ENERGY)
C
C
C----------------------------------------------------------------------
C SUM MASS OF LEVEL K OVER GLOBE AND ADD TO TOTAL ATMOSPHERIC MASS
C ON THE P_GRID
C----------------------------------------------------------------------
C
       CALL MASS_SUM (DELP_P,RS_P_K,AREA_P,
     &                START_POINT_NO_HALO,END_P_POINT_NO_HALO,
     &                P_FIELD,TOT_MASS_P)
C
C
C----------------------------------------------------------------------
C SUM PSTAR*DELBK FOR LEVEL K OVER THE GLOBE ON THE P GRID
C----------------------------------------------------------------------
C
       CALL MASS_SUM (PSTAR_DELBK,RS_P_K,AREA_P,
     &                START_POINT_NO_HALO,END_P_POINT_NO_HALO,
     &                P_FIELD,PART_MASS_P)
C
      IF (LWHITBROM) THEN
C
C----------------------------------------------------------------------
C STORE RADIUS OF SPHERE AT LEVEL K INTO WORK TO ALLOW CALCULATION
C OF RADIUS AT LEVEL K+1
C----------------------------------------------------------------------
C
       DO I=1,P_FIELD
        WORK(I) = RS_P_K(I)
       END DO
C
      END IF
C
C======================================================================
C END OF MAIN LOOP OVER LEVELS
C======================================================================
C
      END DO
C
      RETURN
      END
