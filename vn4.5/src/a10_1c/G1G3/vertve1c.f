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
CLL   SUBROUTINE VERT_VEL -------------------------------------------
CLL
CLL   PURPOSE:  CALCULATES DIVERGENCE FROM MASS-WEIGHTED HORIZONTAL
CLL             VELOCITY COMPONENTS USING EQUATION (30).
CLL             THEN DERIVES MASS-WEIGHTED VERTICAL VELOCITY
CLL             FIELD, EQUATION (29).
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL   WAS VERSION FOR CRAY Y-MP
CLL   WRITTEN BY M.H MAWSON.
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL
!LL   4.4   11/08/97  New version optimised for T3E
!LL                   Not bit reproducible with VERTVE1A.
CLL   4.4   11/08/97  Initialisation of ETADOT_MEAN to zero built into
CLL                   calculation of first adjustment timestep value.
CLL                   A. Dickinson
CLL
!     4.4    17/07/97 SCALAR calculated using SEC_P_LATITUDE at both
!                     poles for non MPP code to enable bit comparison
!                     with MPP code.   I Edmond
CLL  4.5  19/12/97  Set North- & Southmost rows of ETADOT_MEAN to
CLL                 zero for limited area runs.  RTHBarnes.
!     4.5    30/04/98 Loop merging for T3E optimisation for MES
!                     D.Salmond
!
CLL   PROGRAMMING STANDARD:
CLL
CLL   SYSTEM COMPONENTS COVERED: P112
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION:       THE EQUATIONS USED ARE (29) AND (30)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER NO. 10
CLL                        M.J.P. CULLEN,T.DAVIES AND M.H. MAWSON,
CLL                        VERSION 17 DATED 11/02/91.
CLLEND-------------------------------------------------------------

C*L   ARGUMENTS:---------------------------------------------------

      SUBROUTINE VERT_VEL
     1                  (U,V,ETADOT_MEAN,SEC_P_LATITUDE,SUM_DIVERGENCE,
     2                   U_FIELD,P_FIELD,P_LEVELS,
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
     3                   ROW_LENGTH,LATITUDE_STEP_INVERSE,
     4                   LONGITUDE_STEP_INVERSE,ADJUSTMENT_STEPS,AKH,
     5                   BKH,RS,CALL_NUMBER,RECIP_RS_SQUARED_SURFACE,
     6                   PSTAR,LLINTS,LWHITBROM)

      IMPLICIT NONE
      LOGICAL  LLINTS, LWHITBROM

      INTEGER
     *  ROW_LENGTH         !IN    NUMBER OF POINTS PER ROW
     *, P_LEVELS           !IN    NUMBER OF PRESSURE LEVELS OF DATA
     *, P_FIELD            !IN    NUMBER OF POINTS IN PRESSURE FIELD.
     *, U_FIELD            !IN    NUMBER OF POINTS IN VELOCITY FIELD.
     *, ADJUSTMENT_STEPS   !IN    HOLDS NUMBER OF ADJUSTMENT STEPS.
     *, CALL_NUMBER        !IN    CURRENT ADJUSTMENT STEP NUMBER
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
     * U(U_FIELD,P_LEVELS)    !IN. MASS WEIGHTED U VELOCITY.
     *,V(U_FIELD,P_LEVELS)    !IN. MASS WEIGHTED V VELOCITY*
     *                        ! COS(LATITUDE)
     *,SEC_P_LATITUDE(P_FIELD)!IN  1/COS(LAT) AT P POINTS
     *,LONGITUDE_STEP_INVERSE !IN 1/LONGITUDE INCREMENT
     *,LATITUDE_STEP_INVERSE  !IN 1/LATITUDE INCREMENT
     *,BKH(P_LEVELS+1)        !IN. HOLDS COEFFICIENT WHICH
     *                        ! MULTIPLIES PSTAR IN HYBRID CO-ORDS
     *                        ! AT LEVELS K-1/2
     *,AKH(P_LEVELS+1)        !IN. HOLDS FIRST COEFFICIENT
     *                        ! IN HYBRID CO-ORDS AT LEVELS K-1/2
     *,RS(P_FIELD,P_LEVELS)   !IN. RADIUS OF EARTH AT P POINTS.
     *,PSTAR(P_FIELD)         !IN. SURFACE PRESSURE AT P POINTS.

      REAL
     * SUM_DIVERGENCE(P_FIELD,P_LEVELS) !OUT. HOLDS MASS
     *                              ! WEIGHTED VERTICAL VELOCITY.

      REAL
     * ETADOT_MEAN(P_FIELD,P_LEVELS) !INOUT. HOLDS ACCUMULATED MASS-
     *                              ! WEIGHTED VERTICAL VELOCITY DIVIDED
     *                              ! BY NUMBER OF ADJUSTMENT_STEPS.
     *,RECIP_RS_SQUARED_SURFACE(P_FIELD) !INOUT. HOLDS 1./(RS*RS) AT
     *                              ! MODEL SURFACE. SET ON FIRST CALL
     *                              ! AND HELD CONSTANT FOR ALL
     *                              ! SUBSEQUENT ONES.
C*---------------------------------------------------------------------
! Parameters for MPP code

C*L   3 LOCAL ARRAYS NEEDED. -----------------------------------------

      REAL
     *  DU_DLONGITUDE(P_FIELD)
     *, DV_DLATITUDE(P_FIELD)
     *, DV_DLATITUDE2(U_FIELD)
      REAL
     & sum_tmp(ROW_LENGTH-2*EW_Halo,P_LEVELS),
     & sums(P_LEVELS)
C*---------------------------------------------------------------------

C DEFINE COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     *  I,J,K
     *  ,START_POINT,END_POINT
     *,LEVEL
     &, POINTS
      INTEGER info
C DEFINE LOCAL SCALARS
      REAL
     * RECIP_ADJUSTMENT_STEPS
      REAL
     *  SCALAR

      REAL SUM_N,SUM_S

C*---------------------------------------------------------------------
C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
C*----------------------------------------------------------------------

C*L   EXTERNAL SUBROUTINE CALLS:- ( IF LWHITBROM ) ---------------
      EXTERNAL CALC_RS
C*---------------------------------------------------------------------

CL    MAXIMUM VECTOR LENGTH ASSUMED IS ROWS*ROW_LENGTH.
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
! All references to poles in the comments, apply equally to the
! Northern and Southern rows when used in LAM configuration.
! References to halos apply only to MPP code.
CL---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL    SECTION 1. CALCULATE DIVERGENCE AS IN EQUATION (30).
CL---------------------------------------------------------------------

      POINTS=LAST_P_VALID_PT-FIRST_VALID_PT+1
! Number of points to be processed by CALC_RS. For non-MPP runs this
! is simply P_FIELD, for MPP, it is all the points, minus any
! unused halo areas (ie. the halo above North pole row, and beneath
! South pole row)

C LOOP OVER LEVELS
      DO 100 K=1,P_LEVELS

C CALCULATE DU/D(LAMDA)
! Loop over all points except South Pole, missing halos
        DO 110 I=START_POINT_NO_HALO - ROW_LENGTH +1,
     &           END_P_POINT_NO_HALO
          DU_DLONGITUDE(I) = LONGITUDE_STEP_INVERSE*(U(I,K)-U(I-1,K))
 110    CONTINUE

C CALCULATE DV/D(PHI)
! Loop over all non-polar points, missing halos
        DO 120 I=START_POINT_NO_HALO , END_P_POINT_NO_HALO
          DV_DLATITUDE(I) = LATITUDE_STEP_INVERSE*(V(I-ROW_LENGTH,K)
     *                         -V(I,K))
 120    CONTINUE

C CALCULATE AVERAGE OF DV_DLATITUDE
! Loop over all non-polar points, missing halos and first point
        DO 130 I=START_POINT_NO_HALO+1 , END_P_POINT_NO_HALO
          DV_DLATITUDE2(I) = DV_DLATITUDE(I) + DV_DLATITUDE(I-1)
 130    CONTINUE

! Set the first element of arrays where loops have skipped:
! Loop 110:
        DU_DLONGITUDE(START_POINT_NO_HALO - ROW_LENGTH)=
     &    DU_DLONGITUDE(START_POINT_NO_HALO - ROW_LENGTH + 1)
! Loop 130:
        DV_DLATITUDE2(START_POINT_NO_HALO)=
     &    DV_DLATITUDE2(START_POINT_NO_HALO+1)

C CALCULATE DIVERGENCES.

! Loop over all non-polar points, missing halos
        DO 150 J=START_POINT_NO_HALO , END_P_POINT_NO_HALO
          SUM_DIVERGENCE(J,K)= SEC_P_LATITUDE(J)*.5*(DU_DLONGITUDE(J)
     *                           + DU_DLONGITUDE(J-ROW_LENGTH)
     *                           + DV_DLATITUDE2(J))
 150    CONTINUE


C END LOOP OVER LEVELS
 100  CONTINUE
! Calculate divergence at poles by summing DV/D(LAT) around polar
! circle and averaging.

!      START_POINT=TOP_ROW_START+LAST_ROW_PT-1   ! Last point of NP row
!      END_POINT= P_BOT_ROW_START+FIRST_ROW_PT-1 ! First point of SP row
       START_POINT=START_POINT_NO_HALO-1
       END_POINT=END_P_POINT_NO_HALO+1


! New start and end points to include one point of each pole

! Do sum across North Pole
      IF (at_top_of_LPG) THEN
        SCALAR = SEC_P_LATITUDE(TOP_ROW_START)*LATITUDE_STEP_INVERSE /
     &           GLOBAL_ROW_LENGTH
        DO K=1,P_LEVELS
! Copy up the items to be summed into a temporary array
! Loop over all the non-halo points on a row
          DO I=FIRST_ROW_PT,LAST_ROW_PT
            sum_tmp(I-FIRST_ROW_PT+1,K)=
     &         -V(TOP_ROW_START+I-1,K)*SCALAR
          ENDDO
        ENDDO

! And perform the sum
        CALL GCG_RVECSUMR(ROW_LENGTH-2*EW_Halo , ROW_LENGTH-2*EW_Halo,
     &                    1,P_LEVELS,sum_tmp,GC_ROW_GROUP,
     &                    info,sums)

! And store the result back
        DO K=1,P_LEVELS
          SUM_DIVERGENCE(START_POINT,K)=sums(K)
        ENDDO
      ELSE  ! If this processor not at top of LPG
        START_POINT=START_POINT_NO_HALO  ! no North Pole point here
      ENDIF

! And sum across South Pole
      IF (at_base_of_LPG) THEN
        SCALAR=SEC_P_LATITUDE(P_BOT_ROW_START)*LATITUDE_STEP_INVERSE/
     &         GLOBAL_ROW_LENGTH

        DO K=1,P_LEVELS
! Copy up the items to be summed into a temporary array
! Loop over all the non-halo points on a row
          DO I=FIRST_ROW_PT,LAST_ROW_PT
            sum_tmp(I-FIRST_ROW_PT+1,K)=
     &        V(U_BOT_ROW_START+I-1,K)*SCALAR
          ENDDO
        ENDDO

! And perform the sum
        CALL GCG_RVECSUMR(ROW_LENGTH-2*EW_Halo , ROW_LENGTH-2*EW_Halo,
     &                    1,P_LEVELS,sum_tmp,GC_ROW_GROUP,
     &                    info,sums)

! And store the result back
        DO K=1,P_LEVELS
          SUM_DIVERGENCE(END_POINT,K)=sums(K)
        ENDDO
      ELSE  ! If this processor not at the bottom of LPG
        END_POINT=END_P_POINT_NO_HALO  ! no South Pole point here
      ENDIF

CL
CL---------------------------------------------------------------------
CL    SECTION 2. CALCULATE VERTICAL VELOCITY. EQUATION (29).
CL---------------------------------------------------------------------


C ---------------------------------------------------------------------
CL    SECTION 2.1 SUM DIVERGENCES THROUGHOUT ATMOSPHERE.
C ---------------------------------------------------------------------

C BY CODING THE SUMMATION AS FOLLOWS THE VALUES PUT INTO EACH LEVEL
C OF SUM_DIVERGENCE ARE THE ONES NEEDED FOR THE SECOND SUMMATION TERM
C IN EQUATION 29, WHILE THE TOTAL SUM IS HELD IN SUM_DIVERGENCE( ,1)

      DO 210 K=P_LEVELS-1,1,-1
        DO 212 I=START_POINT,END_POINT
          SUM_DIVERGENCE(I,K)= SUM_DIVERGENCE(I,K)+SUM_DIVERGENCE(I,K+1)
 212    CONTINUE
 210  CONTINUE

C ---------------------------------------------------------------------
CL    SECTION 2.2 CALCULATE MASS-WEIGHTED VERTICAL VELOCITY.
CL                CALCULATE 1/(RS*RS) IF THIS IS CALL NUMBER ONE.
C ---------------------------------------------------------------------

      IF(CALL_NUMBER.EQ.1) THEN
CL    CALCULATE 1/(RS*RS) AT MODEL SURFACE

      IF (.NOT.LWHITBROM) THEN

! loop over all points, including valid halos
        DO 220 I=FIRST_VALID_PT,LAST_P_VALID_PT
          RECIP_RS_SQUARED_SURFACE(I) = 1./(A*A)
 220    CONTINUE

      ELSE

        LEVEL=1
C DV_DLATITUDE,DU_DLONGITUDE ARE DUMMY ARRAYS REQUIRED BY CALC_RS AND
C THE CONTENTS TRANSFERED TO AND RETURNED FROM IT ARE IRRELEVANT.
        CALL CALC_RS(PSTAR(FIRST_VALID_PT),AKH,BKH,
     &               DV_DLATITUDE(FIRST_VALID_PT),
     &               DU_DLONGITUDE(FIRST_VALID_PT),
     *               RECIP_RS_SQUARED_SURFACE(FIRST_VALID_PT),
     &               POINTS,LEVEL,P_LEVELS,LLINTS)
! loop over all points, including valid halos
        DO 320 I=FIRST_VALID_PT,LAST_P_VALID_PT
          RECIP_RS_SQUARED_SURFACE(I)= 1./(RECIP_RS_SQUARED_SURFACE(I)*
     *                                 RECIP_RS_SQUARED_SURFACE(I))
 320    CONTINUE

      END IF      !     LWHITBROM

      END IF

C DP/D(PSTAR) IS NOTHING MORE THAN THE BK COEFFICENT.

      DO 222 K= P_LEVELS,2,-1
CFPP$ SELECT(CONCUR)
        DO 224 I=START_POINT,END_POINT
          SUM_DIVERGENCE(I,K)= SUM_DIVERGENCE(I,K) - BKH(K)
     *                          * SUM_DIVERGENCE(I,1)*RS(I,K)*RS(I,K)*
     *                          RECIP_RS_SQUARED_SURFACE(I)
 224    CONTINUE
 222  CONTINUE

C ---------------------------------------------------------------------
CL    SECTION 2.3 ACCUMULATE MASS-WEIGHTED VERTICAL VELOCITY DIVIDED
CL                BY NUMBER OF ADJUSTMENT TIMESTEPS.
C ---------------------------------------------------------------------

      RECIP_ADJUSTMENT_STEPS = 1./ ADJUSTMENT_STEPS
      DO K= 1,P_LEVELS

      IF (CALL_NUMBER.EQ.1)THEN

CFPP$ SELECT(CONCUR)
        DO  I=START_POINT,END_POINT
          ETADOT_MEAN(I,K)= SUM_DIVERGENCE(I,K)
     *                          * RECIP_ADJUSTMENT_STEPS
        ENDDO

      ELSE

CFPP$ SELECT(CONCUR)
        DO  I=START_POINT,END_POINT
          ETADOT_MEAN(I,K)= ETADOT_MEAN(I,K) + SUM_DIVERGENCE(I,K)
     *                          * RECIP_ADJUSTMENT_STEPS
        ENDDO

      ENDIF
      ENDDO

C IF GLOBAL MODEL SET ALL POINTS AT POLES TO THE UNIQUE VALUE.
      DO 240 K=1,P_LEVELS
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        IF (at_top_of_LPG) THEN
! Loop over North Pole points, missing out the last (START_POINT)
! point.
          DO I=TOP_ROW_START+FIRST_ROW_PT-1,
     &         TOP_ROW_START+LAST_ROW_PT-1
            ETADOT_MEAN(I,K)=ETADOT_MEAN(START_POINT,K)
            SUM_DIVERGENCE(I,K)=SUM_DIVERGENCE(START_POINT,K)
          ENDDO
        ENDIF  ! at North Pole

        IF (at_base_of_LPG) THEN
! Loop over South Pole points, missing out the first (END_POINT)
! point.
          DO I=P_BOT_ROW_START+FIRST_ROW_PT-1,
     &         P_BOT_ROW_START+LAST_ROW_PT-1
            ETADOT_MEAN(I,K)=ETADOT_MEAN(END_POINT,K)
            SUM_DIVERGENCE(I,K)=SUM_DIVERGENCE(END_POINT,K)
          ENDDO
        ENDIF  ! at South Pole
 240  CONTINUE

CL    END OF ROUTINE VERT_VEL

      RETURN
      END
