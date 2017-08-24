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
CLL   SUBROUTINE TH_Q_DIF -----------------------------------------
CLL
CLL   PURPOSE:  CALCULATES DIFFUSION INCREMENTS FOR THETAL OR QT
CLL              IF STEEP SLOPE THEN EFFECTIVE DIFFUSION IS ZERO.
CLL
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL
CLL   WAS VERSION FOR CRAY Y-MP
CLL
CLL MM, TJ      <- PROGRAMMER OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
!LL   4.4   11/08/97  New version optimised for T3E.
!LL                   Not bit-reproducible with THQDIF1A.
!     4.4    17/07/97 SCALAR calculated using SEC_P_LATITUDE at both
!                     poles for non MPP code to enable bit comparison
!                     with MPP code.   I Edmond
CLL   4.4    25/07/97 Calling sequence changed from once per diffusion
CLL                   sweep per level to once per dynamics sweep, in
CLL                   order to improve MPP scalability.
CLL                   A. Dickinson
CLL
CLL
CLL   PROGRAMMING STANDARD:
CLL
CLL   SYSTEM COMPONENTS COVERED: P131
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION:       THE EQUATION USED IS (47)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER
CLL                        NO. 10 M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLL                        VERSION 16, DATED 09/01/91.
CLLEND-------------------------------------------------------------

C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE TH_Q_DIF
     1                  (TH_Q,RECIP_RS_SQUARED_DELTAP,
     2                   SEC_P_LATITUDE,ROW_LENGTH,
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
     5                   LEVEL_BASE,LEVELS,
     6                   KEXP_K1,ADVECTION_TIMESTEP,
     4                   P_FIELD,U_FIELD,
     5                   DIFFUSION_EW,DIFFUSION_NS)


      IMPLICIT NONE

      INTEGER
     *  U_FIELD            !IN DIMENSION OF FIELDS ON VELOCITY GRID
     *, P_FIELD            !IN DIMENSION OF FIELDS ON PRESSURE GRID
     *, ROW_LENGTH         !IN NUMBER OF POINTS PER ROW
     *, LEVELS             !IN NUMBER OF LEVELS
     *, LEVEL_BASE         !IN BOTTOM LEVEL AT WHICH DIFFUSION APPLIES
     *, KEXP_K1(LEVELS)    !IN ORDER OF DIFFUSION

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
     *  TH_Q(P_FIELD,LEVELS)    !IN. THETAL OR QT FIELD.
     *, RECIP_RS_SQUARED_DELTAP(P_FIELD,LEVELS) !IN  1/RS**2*DELTAP
     *,ADVECTION_TIMESTEP    !IN

      REAL
     * DIFFUSION_EW(P_FIELD,LEVELS) !IN EAST-WEST EFFECTIVE DIFFUSION
     *,DIFFUSION_NS(P_FIELD,LEVELS) !IN NORTH-SOUTH EFFECTIVE DIFFUSION
     *,SEC_P_LATITUDE(P_FIELD)         !IN 1/COS(LAT) AT P POINTS

C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: 5 ARE REQUIRED

      REAL
     * FIELD1(P_FIELD)      ! GENERAL WORKSPACE
     *,FIELD2(P_FIELD)      ! GENERAL WORKSPACE
     *,FIELD3(P_FIELD)      ! GENERAL WORKSPACE
     *,FIELD(P_FIELD)       ! GENERAL WORKSPACE
     *,FIELD_INC(P_FIELD)   ! GENERAL WORKSPACE
     *,NP_FLUX(ROW_LENGTH)  ! HOLDS NORTH POLAR FLUX
     *,SP_FLUX(ROW_LENGTH)  ! HOLDS SOUTH POLAR FLUX
C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES

C LOCAL REALS.
      REAL
     *  SCALAR
C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     *  I,J,IJ,K,JJ

C*L   EXTERNAL SUBROUTINE CALLS:---------------------------------------
      EXTERNAL
     *  POLAR
C*---------------------------------------------------------------------
CL    MAXIMUM VECTOR LENGTH ASSUMED IS END_P_UPDATE-START_P_UPDATE+1
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL

      DO K=LEVEL_BASE,LEVELS

        DO  I=FIRST_VALID_PT,LAST_P_VALID_PT
          FIELD(I) = TH_Q(I,K)
        END DO

C LOOP THROUGH CODE KEXP_K1 TIMES. THE ORDER OF THE DIFFUSION SCHEME IS
C DEL TO THE POWER 2*KEXP_K1.
        DO  JJ=1,KEXP_K1(K)

CL   ZERO INCREMENTS FOR FIRST AND LAST ROW
CL   OVERWRITTEN BY POLAR IN GLOBAL MODELS
          IF (at_top_of_LPG) THEN
            DO I=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
              FIELD_INC(I)=0.0
            ENDDO
          ENDIF
          IF (at_base_of_LPG) THEN
            DO I=P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
              FIELD_INC(I)=0.0
            ENDDO
          ENDIF


CL---------------------------------------------------------------------
CL    SECTION 1.    DELTA LAMBDA TERMS
CL---------------------------------------------------------------------
C----------------------------------------------------------------------
CL    SECTION 1.1    CALCULATE DELTAPHILAMBDA*1/(DELTALAMBDA)SQUARED
C----------------------------------------------------------------------


      DO  I=START_POINT_NO_HALO,END_P_POINT_NO_HALO-1
       FIELD1(I)=FIELD(I+1)-FIELD(I)
      END DO


! Set last point of field
      FIELD1(END_P_POINT_NO_HALO)=FIELD1(END_P_POINT_NO_HALO-1)

C----------------------------------------------------------------------
CL    SECTION 1.2    CALCULATE DELTA LAMBDA TERM
C----------------------------------------------------------------------

      DO I= START_POINT_NO_HALO+1,END_P_POINT_NO_HALO
       FIELD2(I)=(DIFFUSION_EW(I,K)*FIELD1(I)-
     &            DIFFUSION_EW(I-1,K)*FIELD1(I-1))*
     &            SEC_P_LATITUDE(I)
      END DO


      FIELD2(START_POINT_NO_HALO)=FIELD2(START_POINT_NO_HALO+1)

C----------------------------------------------------------------------
CL    SECTION 2    CALCULATE PHI DIRECTION TERM AND ADD
CL                   ONTO FIRST TERM TO GET TOTAL CORRECTION.
C----------------------------------------------------------------------

C   CALCULATE DELTA PHI TERMS

      DO  I=START_POINT_NO_HALO-ROW_LENGTH,END_P_POINT_NO_HALO
       FIELD1(I)=FIELD(I)-FIELD(I+ROW_LENGTH)
      END DO

C----------------------------------------------------------------------
CL    SECTION 2.3  CALCULATE DELTAPHI TERM AND ADD ONTO DELTALAMBDA TERM
C----------------------------------------------------------------------
cdir$ nosplit
      DO  I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
       FIELD_INC(I)= (FIELD2(I)+
     &            DIFFUSION_NS(I-ROW_LENGTH,K)*FIELD1(I-ROW_LENGTH)-
     &            DIFFUSION_NS(I,K)*FIELD1(I))*SEC_P_LATITUDE(I)
      END DO

C----------------------------------------------------------------------
CL    SECTION 3  CALCULATE DIFFUSION AT POLES
C----------------------------------------------------------------------

CL    LIMITED AREA MODEL ZEROES DEL-SQUARED ON BOUNDARIES.
      IF (at_left_of_LPG) THEN
        DO I=START_POINT_NO_HALO+FIRST_ROW_PT-1,
     &       END_P_POINT_NO_HALO,ROW_LENGTH
          FIELD_INC(I)=0.0
        ENDDO
      ENDIF

      IF (at_right_of_LPG) THEN
        DO I=START_POINT_NO_HALO+LAST_ROW_PT-1,
     &       END_P_POINT_NO_HALO,ROW_LENGTH
          FIELD_INC(I)=0.0
        ENDDO
      ENDIF

C DE-MASS-WEIGHT INCREMENT AND COPY INTO FIELD1
         DO I=FIRST_FLD_PT,LAST_P_FLD_PT
            FIELD(I) = FIELD_INC(I)*RECIP_RS_SQUARED_DELTAP(I,K)
         END DO

         if(jj.ne.KEXP_K1(K))then
         CALL SWAPBOUNDS(FIELD,ROW_LENGTH,tot_P_ROWS,
     &                   EW_Halo,NS_Halo,1)
         endif

C  END OF DIFFUSION SWEEPS
      END DO

CL ADD FINAL INCREMENT ONTO THETAL FIELD.
        SCALAR = (-1)**KEXP_K1(K)
        DO  I=FIRST_VALID_PT,LAST_P_VALID_PT
          TH_Q(I,K) = TH_Q(I,K) - FIELD(I) * ADVECTION_TIMESTEP
     &                   *SCALAR
        END DO

CL END LOOP OVER P_LEVELS FOR THETAL
      END DO
         CALL SWAPBOUNDS
     1  (TH_Q,ROW_LENGTH,tot_P_ROWS,
     &                   EW_Halo,NS_Halo,LEVELS)
CL    END OF ROUTINE TH_Q_DIF

      RETURN
      END
