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
CLL  SUBROUTINE FILT_FLD   -----------------------------------------
CLL
CLL  PURPOSE: FOURIER DAMPS POTENTIAL TEMPERATURE FIELD.
CLL           SETS SURFACE PRESSURE,POTENTIAL TEMPERATURE AND
CLL           MOISTURE VARIABLES AT POLES TO THE MEAN VALUE OF THE
CLL           SURROUNDING ROWS.
CLL  NOT SUITABLE FOR I.B.M USE.
CLL
CLL  WRITTEN BY M.H MAWSON.
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.1     24/02/93  Tidy code to remove QA Fortran messages.
CLL   3.4     26/05/94  Argument LLINTS added and passed to CALC_RS
CLL                                                      S.J.Swarbrick
!LL   4.2     16/08/96  Added TYPFLDPT arguments and made
!LL                     FILTER_WAVE_NUMBER_P_ROWS globally sized.
!LL                     Add TYPFLDPT  args to FILTER.
!LL                                                        P.Burton
!LL   4.3     11/03/97  Added MPP code to for zonal sums.
!LL                     (MPP Non bit-reproducible for different
!LL                      numbers of processors)    P.Burton
!LL   4.5     28/10/98  Introduce Single Column Model. JC Thil
!LL 
CLL
CLL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL  SYSTEM COMPONENTS COVERED:  P142, P196.
CLL  SYSTEM TASK: P1
CLL  DOCUMENTATION:        SEE UNIFIED MODEL DOCUMENTATION PAPER
CLL                        NO. 10 M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLL                        FOR DETAILS OF FOURIER DAMPING.
CLLEND-------------------------------------------------------------

C
C*L  ARGUMENTS:---------------------------------------------------

      SUBROUTINE FILT_FLD
     1                   (P_FIELD,P_LEVELS,Q_LEVELS,ROW_LENGTH,
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
     2                    PSTAR,THETA,Q,QCL,QCF,
     3                    IFAX,TRIGS,FILTER_WAVE_NUMBER_P_ROWS,
     5                    NORTHERN_FILTERED_P_ROW,
     6                    SOUTHERN_FILTERED_P_ROW,
     7                    AK,BK,DELTA_AK,DELTA_BK,COS_P_LATITUDE,
     8                    RS_SQUARED_DELTAP,LATITUDE_STEP_INVERSE,
     9                    LLINTS)

      IMPLICIT NONE

! All FLDPTR arguments are intent IN
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

      LOGICAL  LLINTS  ! Arg passed to CALC_RS

      INTEGER
     1 P_FIELD,      !IN. NUMBER OF PRESSURE POINTS.
     2 ROW_LENGTH,   !IN. NUMBER OF POINTS ON A ROW.
     3 P_LEVELS,     !IN. NUMBER OF MODEL LEVELS.
     4 Q_LEVELS      !IN. NUMBER OF MOIST MODEL LEVELS.

      REAL
     1 PSTAR(P_FIELD),  !INOUT. PRIMARY ARRAY FOR SURFACE PRESSURE
     2 THETA(P_FIELD,P_LEVELS),!INOUT.PRIMARY ARRAY FOR POT. TEMP.
     3 Q(P_FIELD,Q_LEVELS), !INOUT. PRIMARY ARRAY FOR MOISTURE.
     4 QCL(P_FIELD,Q_LEVELS), !INOUT. PRIMARY ARRAY FOR CLOUD LIQUID
     4                        !       WATER.
     5 QCF(P_FIELD,Q_LEVELS)  !INOUT. PRIMARY ARRAY FOR CLOUD FROZEN
     5                        !       WATER.

      INTEGER
     *  NORTHERN_FILTERED_P_ROW !IN P ROW ON WHICH FILTERING STOPS
     *                          ! MOVING TOWARDS EQUATOR
     *, SOUTHERN_FILTERED_P_ROW !IN P ROW ON WHICH FILTERING STARTS
     *                          ! AGAIN MOVING TOWARDS SOUTH POLE
     &, FILTER_WAVE_NUMBER_P_ROWS(GLOBAL_P_FIELD/GLOBAL_ROW_LENGTH)
     &               ! LAST WAVE NUMBER NOT TO BE DAMPED ON A P ROW
     *, IFAX(10)           !IN HOLDS FACTORS OF ROW_LENGTH USED BY
     *                     ! FILTERING.

      REAL
     * TRIGS(ROW_LENGTH)      !IN HOLDS TRIGONOMETRIC FUNCTIONS USED
     *                        ! IN FILTERING.
     *,COS_P_LATITUDE(P_FIELD)!IN HOLDS COSINES OF LATITUDE AT P POINTS
     *,LATITUDE_STEP_INVERSE  !IN 1./(LATITUDE STEP IN RADIANS)
     *,AK(P_LEVELS)           !IN A PART OF ETA CO-ORDINATE
     *,BK(P_LEVELS)           !IN B PART OF ETA CO-ORDINATE
     *,DELTA_AK(P_LEVELS)     !IN LAYER THICKNESS OF A PART OF ETA
     *,DELTA_BK(P_LEVELS)     !IN LAYER THICKNESS OF B PART OF ETA
     *,RS_SQUARED_DELTAP(P_FIELD,P_LEVELS) !IN SPACE USED TO PUT
     *                                     ! MASS FIELD IN.

C*---------------------------------------------------------------------

C*L  DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C   10 LOCAL ARRAYS REQUIRED

      REAL
     & MEAN_MW_THETA(tot_P_ROWS,P_LEVELS)
     &,MEAN_MW_THETA_NEW(tot_P_ROWS,P_LEVELS)
     &,MEAN_MASS(tot_P_ROWS,P_LEVELS)
     &,MEAN_MW_NP_THETA_NEW(P_LEVELS)
     &,MEAN_MW_SP_THETA_NEW(P_LEVELS)
     &,MEAN_MW_NP_Q_NEW(Q_LEVELS)
     &,MEAN_MW_SP_Q_NEW(Q_LEVELS)
     &,MEAN_MW_NP_QCL_NEW(Q_LEVELS)
     &,MEAN_MW_SP_QCL_NEW(Q_LEVELS)
     &,MEAN_MW_NP_QCF_NEW(Q_LEVELS)
     &,MEAN_MW_SP_QCF_NEW(Q_LEVELS)
     &,MEAN_MASS_NP(P_LEVELS)
     &,MEAN_MASS_SP(P_LEVELS)
     &,NP_THETA(P_LEVELS)
     &,SP_THETA(P_LEVELS)
     &,NP_Q(Q_LEVELS)
     &,SP_Q(Q_LEVELS)
     &,NP_QCL(Q_LEVELS)
     &,SP_QCL(Q_LEVELS)
     &,NP_QCF(Q_LEVELS)
     &,SP_QCF(Q_LEVELS)
     &,MEAN_MW_NP_THETA(P_LEVELS)
     &,MEAN_MW_SP_THETA(P_LEVELS)
     &,MEAN_MW_NP_Q(Q_LEVELS)
     &,MEAN_MW_SP_Q(Q_LEVELS)
     &,MEAN_MW_NP_QCL(Q_LEVELS)
     &,MEAN_MW_SP_QCL(Q_LEVELS)
     &,MEAN_MW_NP_QCF(Q_LEVELS)
     &,MEAN_MW_SP_QCF(Q_LEVELS)
     &,WORK1(P_FIELD)
C*---------------------------------------------------------------------

C DEFINE LOCAL VARIABLES

      INTEGER
     *  FILTER_SPACE_P     ! HORIZONTAL DIMENSION OF SPACE NEEDED IN
     *                     ! FILTERING ROUTINE FOR P ROWS.
     &,  POINTS  ! number of updatable points
     &,  NORTH_FIRST_ROW,NORTH_LAST_ROW  ! limits for loop over
     &,  SOUTH_FIRST_ROW,SOUTH_LAST_ROW  ! filterable rows

      INTEGER
     1 I,K,J
     &, info  ! return code for communcations

      REAL
     & INCREMENT
     &,MEAN_RADIUS_NP
     &,MEAN_RADIUS_SP

      REAL
     & POLAR_COSINE

      REAL
     * NP_PSTAR,
     * SP_PSTAR
C ---------------------------------------------------------------------

C*L   EXTERNAL SUBROUTINE CALLS:---------------------------------------
      EXTERNAL FILTER,CALC_RS

C*---------------------------------------------------------------------

CL  MAXIMUM VECTOR LENGTH ASSUMED IS ROW_LENGTH
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL    SECTION 1.    INITIALISE CONSTANTS.
CL---------------------------------------------------------------------

C SET FILTER_SPACE WHICH IS ROW_LENGTH+2 TIMES THE NUMBER OF ROWS TO
C BE FILTERED.

      FILTER_SPACE_P = (ROW_LENGTH+2)*(NORTHERN_FILTERED_P_ROW-1+
     *                P_FIELD/ROW_LENGTH-SOUTHERN_FILTERED_P_ROW)

      POINTS=LAST_P_VALID_PT-FIRST_VALID_PT+1

! Set FIRST_ROW and LAST_ROW variables to point to the sections
! of the field being updated by FILTER      
! For the MPP code we must convert from global row numbers to local
! row numbers.
      NORTH_FIRST_ROW=FIRST_ROW
      NORTH_LAST_ROW=NORTHERN_FILTERED_P_ROW-FIRST_GLOBAL_ROW_NUMBER+
     &               NS_Halo+1  ! gives local row number
      IF (NORTH_LAST_ROW .GT. (tot_P_ROWS-NS_Halo))
     &  NORTH_LAST_ROW=tot_P_ROWS-NS_Halo
     
      SOUTH_FIRST_ROW=SOUTHERN_FILTERED_P_ROW-FIRST_GLOBAL_ROW_NUMBER+
     &                NS_Halo+1  ! gives local row number
      IF (SOUTH_FIRST_ROW .LT. (NS_Halo+1))
     &  SOUTH_FIRST_ROW=NS_Halo+1
      SOUTH_LAST_ROW=P_LAST_ROW
      
      

      POLAR_COSINE = 0.125/LATITUDE_STEP_INVERSE

CL
CL---------------------------------------------------------------------
CL    SECTION 2.    CALCULATE RS SQUARED AND MASS-WEIGHTED THETA ON
CL                  EACH ROW AT EACH LEVEL.
CL---------------------------------------------------------------------

CL    CALL CALC_RS TO GET RS FOR LEVEL 1.
C RS IS RETURNED IN RS_SQUARED_DELTAP( ,1)
C TS IS RETURNED IN WORK1, RS AT LEVEL K-1 IS INPUT IN
C RS_SQUARED_DELTAP( ,2) AS AT K-1= 0 THE INPUT IS NOT USED BY CALC_RS.

      CALL CALC_RS(PSTAR(FIRST_VALID_PT),AK,BK,
     &             WORK1(FIRST_VALID_PT),
     &             RS_SQUARED_DELTAP(FIRST_VALID_PT,2),
     &             RS_SQUARED_DELTAP(FIRST_VALID_PT,1),
     &             POINTS,1,P_LEVELS,LLINTS)

CL LOOP FROM 2 TO P_LEVELS
      DO K= 2,P_LEVELS

CL    CALL CALC_RS TO GET RS FOR LEVEL K.
C RS IS RETURNED IN RS_SQUARED_DELTAP(1,K)
C TS IS RETURNED IN WORK1, RS AT LEVEL K-1 IS INPUT AS
C RS_SQUARED_DELTAP(K-1).

        I=K
      CALL CALC_RS(PSTAR(FIRST_VALID_PT),AK,BK,
     &             WORK1(FIRST_VALID_PT),
     &             RS_SQUARED_DELTAP(FIRST_VALID_PT,K-1),
     &             RS_SQUARED_DELTAP(FIRST_VALID_PT,K),
     &             POINTS,I,P_LEVELS,LLINTS)

      END DO

CL END LOOP FROM 2 TO P_LEVELS.

CL FORM RS SQUARED * DELTA P * COSINE OF LATITUDE
CL AND ZONAL MEAN MASS-WEIGHTED THETA

      DO K=1,P_LEVELS
        DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
          RS_SQUARED_DELTAP(I,K) = RS_SQUARED_DELTAP(I,K)*
     &                             RS_SQUARED_DELTAP(I,K)*
     &                             (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))
     &                             *COS_P_LATITUDE(I)
        END DO
C SET POLAR VALUES.
C THE CORRECT COSINE VALUE IS DELTA_PHI/8
        IF (at_top_of_LPG) THEN
          DO I=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
            RS_SQUARED_DELTAP(I,K)=RS_SQUARED_DELTAP(I,K)*
     &        RS_SQUARED_DELTAP(I,K)*
     &        (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))*
     &        POLAR_COSINE
          ENDDO
        ENDIF
        
        IF (at_base_of_LPG) THEN
          DO I=P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
             RS_SQUARED_DELTAP(I,K)=RS_SQUARED_DELTAP(I,K)*
     &                             RS_SQUARED_DELTAP(I,K)*
     &                             (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))*
     &                             POLAR_COSINE
!         DONE
          ENDDO
        ENDIF
        DO J=1,tot_P_ROWS
          MEAN_MW_THETA(J,K)=0.0
          MEAN_MW_THETA_NEW(J,K)=0.0
          MEAN_MASS(J,K)=0.0
        ENDDO
        
        DO J=NORTH_FIRST_ROW,NORTH_LAST_ROW  
!       loop over rows to be filtered

          DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &         (J-1)*ROW_LENGTH+LAST_ROW_PT
            MEAN_MW_THETA(J,K) = MEAN_MW_THETA(J,K) + THETA(I,K)*
     &                           RS_SQUARED_DELTAP(I,K)
          END DO
        END DO
        DO J=SOUTH_FIRST_ROW,SOUTH_LAST_ROW  
!       loop over rows to be filtered

          DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &         (J-1)*ROW_LENGTH+LAST_ROW_PT
            MEAN_MW_THETA(J,K) = MEAN_MW_THETA(J,K) + THETA(I,K)*
     &                           RS_SQUARED_DELTAP(I,K)
          END DO
        END DO
      END DO

! So far MEAN_MW_THETA contains only the sum along my local part
! of the row. We must now do a sum, so that it contains the full
! sum for the entire global row
! NB : Since the partial sums on each processor will be different
! depending on the number of processors in the EW direction, the
! total sum will also be non-reproducible if the number of EW
! processors change.

      CALL GCG_RSUM(tot_P_ROWS*P_LEVELS,GC_ROW_GROUP,info,
     &              MEAN_MW_THETA)
     


CL
CL---------------------------------------------------------------------
CL    SECTION 3.    FILTER THETA FIELD.
CL---------------------------------------------------------------------

CL    CALL FILTER FOR THETA

        CALL FILTER(THETA,P_FIELD,P_LEVELS,FILTER_SPACE_P,ROW_LENGTH,
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
     *              FILTER_WAVE_NUMBER_P_ROWS,TRIGS,IFAX,
     *              NORTHERN_FILTERED_P_ROW,SOUTHERN_FILTERED_P_ROW)

CL
CL---------------------------------------------------------------------
CL    SECTION 4.    CALCULATE MASS-WEIGHTED THETA AFTER FILTERING.
CL                  CALCULATE CHANGE DUE TO FILTERING AND ADD AN
CL                  INCREMENT TO EACH POINT TO RETAIN CONSERVATION.
CL---------------------------------------------------------------------

CL CALCULATE ZONAL MEAN MASS-WEIGHTED THETA
CL CALCULATE INCREMENT NEEDED TO EACH THETA VALUE TO ENSURE
CL CONSERVATION.

      DO K=1,P_LEVELS
        DO J=NORTH_FIRST_ROW,NORTH_LAST_ROW
          DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &         (J-1)*ROW_LENGTH+LAST_ROW_PT
            MEAN_MW_THETA_NEW(J,K) = MEAN_MW_THETA_NEW(J,K) + 
     &                               THETA(I,K)*RS_SQUARED_DELTAP(I,K)
            MEAN_MASS(J,K)=MEAN_MASS(J,K) + RS_SQUARED_DELTAP(I,K)
          ENDDO
        ENDDO
        
        DO J=SOUTH_FIRST_ROW,SOUTH_LAST_ROW
          DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &         (J-1)*ROW_LENGTH+LAST_ROW_PT
            MEAN_MW_THETA_NEW(J,K) = MEAN_MW_THETA_NEW(J,K) + 
     &                               THETA(I,K)*RS_SQUARED_DELTAP(I,K)
            MEAN_MASS(J,K)=MEAN_MASS(J,K) + RS_SQUARED_DELTAP(I,K)
          ENDDO
        ENDDO
      ENDDO
        
! Do sum along rows for MEAN_MW_THETA_NEW and MEAN_MASS as before

      CALL GCG_RSUM(tot_P_ROWS*P_LEVELS,GC_ROW_GROUP,info,
     &              MEAN_MW_THETA_NEW)
      CALL GCG_RSUM(tot_P_ROWS*P_LEVELS,GC_ROW_GROUP,info,
     &              MEAN_MASS)

      DO K=1,P_LEVELS
        DO J=NORTH_FIRST_ROW,NORTH_LAST_ROW
          INCREMENT=(MEAN_MW_THETA_NEW(J,K)-MEAN_MW_THETA(J,K))/
     &              MEAN_MASS(J,K)
          DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &         (J-1)*ROW_LENGTH+LAST_ROW_PT
            THETA(I,K) = THETA(I,K) - INCREMENT
          ENDDO
        ENDDO
        
        DO J=SOUTH_FIRST_ROW,SOUTH_LAST_ROW
          INCREMENT=(MEAN_MW_THETA_NEW(J,K)-MEAN_MW_THETA(J,K))/
     &              MEAN_MASS(J,K)
          DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &         (J-1)*ROW_LENGTH+LAST_ROW_PT
            THETA(I,K) = THETA(I,K) - INCREMENT
          ENDDO
        ENDDO
        
      ENDDO


CL
CL---------------------------------------------------------------------
CL    SECTION 5.    SET THETA,Q,QCL,QCF AND PSTAR AT POLES TO MEAN OF
CL                  SURROUNDING ROW IN A CONSERVATIVE WAY.
CL---------------------------------------------------------------------

C ---------------------------------------------------------------------
CL    SECTION 5.1   CALCULATE MEAN MASS-WEIGHTED VALUES OF FIELDS
CL                  AROUND POLES.
C ---------------------------------------------------------------------

C CALCULATE MEAN MASS-WEIGHTED VALUES OF ALL FIELDS AROUND POLAR CAPS
C REMOVE DELTA P FROM RS_SQUARED FIELD.
! and calculate mean of pstar in row adjacent to pole
      DO K=1,P_LEVELS
        MEAN_MW_NP_THETA(K) = 0.0
        MEAN_MW_SP_THETA(K) = 0.0
        IF (K .LE. Q_LEVELS) THEN
          MEAN_MW_NP_Q(K) = 0.0
          MEAN_MW_SP_Q(K) = 0.0
          MEAN_MW_NP_QCL(K) = 0.0
          MEAN_MW_SP_QCL(K) = 0.0
          MEAN_MW_NP_QCF(K) = 0.0
          MEAN_MW_SP_QCF(K) = 0.0
        ENDIF
        
        IF (at_top_of_LPG) THEN
          IF (K .LE. Q_LEVELS) THEN
            DO J=FIRST_ROW-1,FIRST_ROW  ! NP and adjacent row
              DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &             (J-1)*ROW_LENGTH+LAST_ROW_PT
                MEAN_MW_NP_Q(K) = MEAN_MW_NP_Q(K) + Q(I,K)*
     &                              RS_SQUARED_DELTAP(I,K)
                MEAN_MW_NP_QCL(K) = MEAN_MW_NP_QCL(K) + QCL(I,K)*
     &                              RS_SQUARED_DELTAP(I,K)
                MEAN_MW_NP_QCF(K) = MEAN_MW_NP_QCF(K) + QCF(I,K)*
     &                              RS_SQUARED_DELTAP(I,K)
              ENDDO
            ENDDO
          ENDIF
          DO J=FIRST_ROW-1,FIRST_ROW  ! NP and adjacent row
            DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &           (J-1)*ROW_LENGTH+LAST_ROW_PT     
              MEAN_MW_NP_THETA(K) = MEAN_MW_NP_THETA(K) + THETA(I,K)*
     &                              RS_SQUARED_DELTAP(I,K)
              RS_SQUARED_DELTAP(I,K) = RS_SQUARED_DELTAP(I,K)/
     &                             (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))
            ENDDO
          ENDDO
          IF (K .EQ. 1) THEN
            MEAN_RADIUS_NP = 0.0
            DO J=FIRST_ROW-1,FIRST_ROW  ! NP and adjacent row
              DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &             (J-1)*ROW_LENGTH+LAST_ROW_PT
                MEAN_RADIUS_NP = MEAN_RADIUS_NP+RS_SQUARED_DELTAP(I,1)
              ENDDO
            ENDDO
            NP_PSTAR=0.0
            DO I=TOP_ROW_START+FIRST_ROW_PT-1,
     &           TOP_ROW_START+LAST_ROW_PT-1
              NP_PSTAR=NP_PSTAR+PSTAR(I+ROW_LENGTH)
            ENDDO
          ENDIF
                
            
          
        ENDIF
        
        IF (at_base_of_LPG) THEN
          IF (K .LE. Q_LEVELS) THEN
            DO J=P_LAST_ROW,P_LAST_ROW+1  ! SP and adjacent row
              DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &             (J-1)*ROW_LENGTH+LAST_ROW_PT
                MEAN_MW_SP_Q(K) = MEAN_MW_SP_Q(K) + Q(I,K)*
     &                              RS_SQUARED_DELTAP(I,K)
                MEAN_MW_SP_QCL(K) = MEAN_MW_SP_QCL(K) + QCL(I,K)*
     &                              RS_SQUARED_DELTAP(I,K)
                MEAN_MW_SP_QCF(K) = MEAN_MW_SP_QCF(K) + QCF(I,K)*
     &                              RS_SQUARED_DELTAP(I,K)
              ENDDO
            ENDDO
          ENDIF
          DO J=P_LAST_ROW,P_LAST_ROW+1  ! SP and adjacent row
            DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &           (J-1)*ROW_LENGTH+LAST_ROW_PT     
              MEAN_MW_SP_THETA(K) = MEAN_MW_SP_THETA(K) + THETA(I,K)*
     &                              RS_SQUARED_DELTAP(I,K)
              RS_SQUARED_DELTAP(I,K) = RS_SQUARED_DELTAP(I,K)/
     &                              (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))
            ENDDO
          ENDDO
          IF (K .EQ. 1) THEN
            MEAN_RADIUS_SP = 0.0
            DO J=P_LAST_ROW,P_LAST_ROW+1  ! NP and adjacent row
              DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &             (J-1)*ROW_LENGTH+LAST_ROW_PT
                MEAN_RADIUS_SP = MEAN_RADIUS_SP+RS_SQUARED_DELTAP(I,1)
              ENDDO
            ENDDO
            SP_PSTAR=0.0
            DO I=P_BOT_ROW_START+FIRST_ROW_PT-1,
     &           P_BOT_ROW_START+LAST_ROW_PT-1
              SP_PSTAR=SP_PSTAR+PSTAR(I-ROW_LENGTH)
            ENDDO
          ENDIF
                    
        ENDIF
      ENDDO ! K : loop over levels
      
! Need to sum the partial sums for the polar rows
! Once again, these sums will give different answers if the number of
! processors in the EW direction changes

      IF (at_top_of_LPG) THEN
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,MEAN_MW_NP_Q)
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,MEAN_MW_NP_QCL)
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,MEAN_MW_NP_QCF)
        CALL GCG_RSUM(P_LEVELS,GC_ROW_GROUP,info,MEAN_MW_NP_THETA)
        CALL GCG_RSUM(1,GC_ROW_GROUP,info,MEAN_RADIUS_NP)
        CALL GCG_RSUM(1,GC_ROW_GROUP,info,NP_PSTAR)        
      ENDIF
      IF (at_base_of_LPG) THEN
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,MEAN_MW_SP_Q)
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,MEAN_MW_SP_QCL)
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,MEAN_MW_SP_QCF)
        CALL GCG_RSUM(P_LEVELS,GC_ROW_GROUP,info,MEAN_MW_SP_THETA)
        CALL GCG_RSUM(1,GC_ROW_GROUP,info,MEAN_RADIUS_SP)
        CALL GCG_RSUM(1,GC_ROW_GROUP,info,SP_PSTAR)
      ENDIF
      

C ---------------------------------------------------------------------
CL    SECTION 5.2   CORRECT PSTAR VALUES.
C ---------------------------------------------------------------------


      IF (at_top_of_LPG) THEN
        NP_PSTAR = NP_PSTAR / GLOBAL_ROW_LENGTH
        IF (MY_PROC_ID .EQ. 0) THEN
          INCREMENT=GLOBAL_ROW_LENGTH*
     &      RS_SQUARED_DELTAP(TOP_ROW_START+FIRST_ROW_PT-1,1)*
     &      (NP_PSTAR-PSTAR(TOP_ROW_START+FIRST_ROW_PT-1))/
     &      MEAN_RADIUS_NP

        ENDIF
           
! We want all processors in polar row to have same value of
! INCREMENT as has been calculated by PE 0
        CALL GCG_RBCAST(101,1,0,GC_ROW_GROUP,info,INCREMENT)

        DO I=TOP_ROW_START+FIRST_ROW_PT-1,
     &       TOP_ROW_START+LAST_ROW_PT-1
          PSTAR(I)=NP_PSTAR - INCREMENT  
          PSTAR(I+ROW_LENGTH)=PSTAR(I+ROW_LENGTH) - INCREMENT
        ENDDO

      ENDIF
      
      IF (at_base_of_LPG) THEN
        SP_PSTAR = SP_PSTAR / GLOBAL_ROW_LENGTH
        IF (MY_PROC_ID .EQ. N_PROCS-1) THEN
          INCREMENT=GLOBAL_ROW_LENGTH*
     &      RS_SQUARED_DELTAP(P_BOT_ROW_START+LAST_ROW_PT-1,1)*
     &      (SP_PSTAR-PSTAR(P_BOT_ROW_START+LAST_ROW_PT-1))/
     &      MEAN_RADIUS_SP

        ENDIF
           
! We want all processors in polar row to have same value of
! INCREMENT as has been calculated by PE 0
        CALL GCG_RBCAST(101,1,N_PROCS-1,GC_ROW_GROUP,info,INCREMENT)

        DO I=P_BOT_ROW_START+FIRST_ROW_PT-1,
     &       P_BOT_ROW_START+LAST_ROW_PT-1
          PSTAR(I)=SP_PSTAR - INCREMENT  
          PSTAR(I-ROW_LENGTH)=PSTAR(I-ROW_LENGTH) - INCREMENT
        ENDDO

      ENDIF


C ---------------------------------------------------------------------
CL    SECTION 5.3   CORRECT VALUES OF OTHER FIELDS.
C ---------------------------------------------------------------------


      DO K=1,P_LEVELS
      
        IF (at_top_of_LPG) THEN
          DO J=FIRST_ROW-1,FIRST_ROW  ! NP and adjacent row
            DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &           (J-1)*ROW_LENGTH+LAST_ROW_PT     
              RS_SQUARED_DELTAP(I,K) = RS_SQUARED_DELTAP(I,K)*
     &                             (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))
            ENDDO
          ENDDO
          
          NP_THETA(K)=0.0
          
          DO I=TOP_ROW_START+FIRST_ROW_PT-1,
     &         TOP_ROW_START+LAST_ROW_PT-1
            NP_THETA(K)=NP_THETA(K) + THETA(I+ROW_LENGTH,K)
          ENDDO
          
          IF (K .LE. Q_LEVELS) THEN
          
            NP_Q(K)=0.0
            NP_QCL(K)=0.0
            NP_QCF(K)=0.0
            
            DO I=TOP_ROW_START+FIRST_ROW_PT-1,
     &           TOP_ROW_START+LAST_ROW_PT-1
              NP_Q(K)=NP_Q(K)+Q(I+ROW_LENGTH,K)
              NP_QCL(K)=NP_QCL(K)+QCL(I+ROW_LENGTH,K)
              NP_QCF(K)=NP_QCF(K)+QCF(I+ROW_LENGTH,K)
            ENDDO
          ENDIF
        ENDIF
        
        IF (at_base_of_LPG) THEN
          DO J=P_LAST_ROW,P_LAST_ROW+1  ! SP and adjacent row
            DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &           (J-1)*ROW_LENGTH+LAST_ROW_PT     
              RS_SQUARED_DELTAP(I,K) = RS_SQUARED_DELTAP(I,K)*
     &                             (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))
            ENDDO
          ENDDO
          
          SP_THETA(K)=0.0
          
          DO I=P_BOT_ROW_START+FIRST_ROW_PT-1,
     &         P_BOT_ROW_START+LAST_ROW_PT-1
            SP_THETA(K)=SP_THETA(K) + THETA(I-ROW_LENGTH,K)
          ENDDO
          
          IF (K .LE. Q_LEVELS) THEN
          
            SP_Q(K)=0.0
            SP_QCL(K)=0.0
            SP_QCF(K)=0.0
            
            DO I=P_BOT_ROW_START+FIRST_ROW_PT-1,
     &           P_BOT_ROW_START+LAST_ROW_PT-1
              SP_Q(K)=SP_Q(K)+Q(I-ROW_LENGTH,K)
              SP_QCL(K)=SP_QCL(K)+QCL(I-ROW_LENGTH,K)
              SP_QCF(K)=SP_QCF(K)+QCF(I-ROW_LENGTH,K)
            ENDDO
          ENDIF
        ENDIF
      ENDDO      
               
! Need to sum the partial sums for the polar rows
! Once again, these sums will give different answers if the number of
! processors in the EW direction changes
      IF (at_top_of_LPG) THEN
        CALL GCG_RSUM(P_LEVELS,GC_ROW_GROUP,info,NP_THETA)
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,NP_Q)
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,NP_QCL)
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,NP_QCF)
      ENDIF
      IF (at_base_of_LPG) THEN
        CALL GCG_RSUM(P_LEVELS,GC_ROW_GROUP,info,SP_THETA)
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,SP_Q)
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,SP_QCL)
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,SP_QCF)
      ENDIF

      DO K=1,P_LEVELS
      
        IF (at_top_of_LPG) THEN
          NP_THETA(K)=NP_THETA(K)/GLOBAL_ROW_LENGTH
          
          DO I=TOP_ROW_START+FIRST_ROW_PT-1,
     &         TOP_ROW_START+LAST_ROW_PT-1
            THETA(I,K)=NP_THETA(K)
          ENDDO
          
          MEAN_MW_NP_THETA_NEW(K)=0.0
          MEAN_MASS_NP(K)=0.0
          DO J=FIRST_ROW-1,FIRST_ROW  ! NP and adjacent row
            DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &           (J-1)*ROW_LENGTH+LAST_ROW_PT
              MEAN_MW_NP_THETA_NEW(K)=MEAN_MW_NP_THETA_NEW(K)+
     &                                THETA(I,K)*
     &                                RS_SQUARED_DELTAP(I,K)
              MEAN_MASS_NP(K)=MEAN_MASS_NP(K)+RS_SQUARED_DELTAP(I,K)
            ENDDO
          ENDDO
          
          IF (K .LE. Q_LEVELS) THEN                             
                    
            NP_Q(K)=NP_Q(K)/GLOBAL_ROW_LENGTH
            NP_QCL(K)=NP_QCL(K)/GLOBAL_ROW_LENGTH
            NP_QCF(K)=NP_QCF(K)/GLOBAL_ROW_LENGTH
            
            DO I=TOP_ROW_START+FIRST_ROW_PT-1,
     &           TOP_ROW_START+LAST_ROW_PT-1
              Q(I,K)=NP_Q(K)
              QCL(I,K)=NP_QCL(K)
              QCF(I,K)=NP_QCF(K)
            ENDDO
            
            MEAN_MW_NP_Q_NEW(K)=0.0
            MEAN_MW_NP_QCL_NEW(K)=0.0
            MEAN_MW_NP_QCF_NEW(K)=0.0
            DO J=FIRST_ROW-1,FIRST_ROW  ! NP and adjacent row
              DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &             (J-1)*ROW_LENGTH+LAST_ROW_PT
                MEAN_MW_NP_Q_NEW(K)=MEAN_MW_NP_Q_NEW(K)+
     &                                  Q(I,K)*RS_SQUARED_DELTAP(I,K)
                MEAN_MW_NP_QCL_NEW(K)=MEAN_MW_NP_QCL_NEW(K)+
     &                                  QCL(I,K)*RS_SQUARED_DELTAP(I,K)
                MEAN_MW_NP_QCF_NEW(K)=MEAN_MW_NP_QCF_NEW(K)+
     &                                  QCF(I,K)*RS_SQUARED_DELTAP(I,K)
              ENDDO
            ENDDO
          ENDIF  ! is this a wet level
        ENDIF  ! at_top_of_LPG
        
        IF (at_base_of_LPG) THEN
          SP_THETA(K)=SP_THETA(K)/GLOBAL_ROW_LENGTH
          
          DO I=P_BOT_ROW_START+FIRST_ROW_PT-1,
     &         P_BOT_ROW_START+LAST_ROW_PT-1
            THETA(I,K)=SP_THETA(K)
          ENDDO
          
          MEAN_MW_SP_THETA_NEW(K)=0.0
          MEAN_MASS_SP(K)=0.0
          DO J=P_LAST_ROW,P_LAST_ROW+1  ! SP and adjacent row
            DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &           (J-1)*ROW_LENGTH+LAST_ROW_PT
              MEAN_MW_SP_THETA_NEW(K)=MEAN_MW_SP_THETA_NEW(K)+
     &                                THETA(I,K)*
     &                                RS_SQUARED_DELTAP(I,K)
              MEAN_MASS_SP(K)=MEAN_MASS_SP(K)+RS_SQUARED_DELTAP(I,K)
            ENDDO
          ENDDO
          
          IF (K .LE. Q_LEVELS) THEN                             
                    
            SP_Q(K)=SP_Q(K)/GLOBAL_ROW_LENGTH
            SP_QCL(K)=SP_QCL(K)/GLOBAL_ROW_LENGTH
            SP_QCF(K)=SP_QCF(K)/GLOBAL_ROW_LENGTH
            
            DO I=P_BOT_ROW_START+FIRST_ROW_PT-1,
     &         P_BOT_ROW_START+LAST_ROW_PT-1
              Q(I,K)=SP_Q(K)
              QCL(I,K)=SP_QCL(K)
              QCF(I,K)=SP_QCF(K)
            ENDDO
            
            MEAN_MW_SP_Q_NEW(K)=0.0
            MEAN_MW_SP_QCL_NEW(K)=0.0
            MEAN_MW_SP_QCF_NEW(K)=0.0
            DO J=P_LAST_ROW,P_LAST_ROW+1  ! SP and adjacent row
              DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &             (J-1)*ROW_LENGTH+LAST_ROW_PT
                MEAN_MW_SP_Q_NEW(K)=MEAN_MW_SP_Q_NEW(K)+
     &                                  Q(I,K)*RS_SQUARED_DELTAP(I,K)
                MEAN_MW_SP_QCL_NEW(K)=MEAN_MW_SP_QCL_NEW(K)+
     &                                  QCL(I,K)*RS_SQUARED_DELTAP(I,K)
                MEAN_MW_SP_QCF_NEW(K)=MEAN_MW_SP_QCF_NEW(K)+
     &                                  QCF(I,K)*RS_SQUARED_DELTAP(I,K)
              ENDDO
            ENDDO
          ENDIF  ! is this a wet level
        ENDIF  ! at_base_of_LPG
      ENDDO ! K: loop over levels
               
! Need to sum the partial sums for the polar rows
! Once again, these sums will give different answers if the number of
! processors in the EW direction changes
      IF (at_top_of_LPG) THEN
        CALL GCG_RSUM(P_LEVELS,GC_ROW_GROUP,info,MEAN_MW_NP_THETA_NEW)
        CALL GCG_RSUM(P_LEVELS,GC_ROW_GROUP,info,MEAN_MASS_NP)
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,MEAN_MW_NP_Q_NEW)
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,MEAN_MW_NP_QCL_NEW)
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,MEAN_MW_NP_QCF_NEW)
      ENDIF
      IF (at_base_of_LPG) THEN
        CALL GCG_RSUM(P_LEVELS,GC_ROW_GROUP,info,MEAN_MW_SP_THETA_NEW)
        CALL GCG_RSUM(P_LEVELS,GC_ROW_GROUP,info,MEAN_MASS_SP)
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,MEAN_MW_SP_Q_NEW)
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,MEAN_MW_SP_QCL_NEW)
        CALL GCG_RSUM(Q_LEVELS,GC_ROW_GROUP,info,MEAN_MW_SP_QCF_NEW)
      ENDIF

      DO K=1,P_LEVELS
      
        IF (at_top_of_LPG) THEN
          MEAN_MW_NP_THETA_NEW(K) = (MEAN_MW_NP_THETA_NEW(K) -
     &      MEAN_MW_NP_THETA(K)) /MEAN_MASS_NP(K)

          DO J=FIRST_ROW-1,FIRST_ROW  ! NP and adjacent row
            DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &           (J-1)*ROW_LENGTH+LAST_ROW_PT
              THETA(I,K) = THETA(I,K) - MEAN_MW_NP_THETA_NEW(K)
            ENDDO
          ENDDO
          
          IF (K .LE. Q_LEVELS) THEN
            MEAN_MW_NP_Q_NEW(K) = (MEAN_MW_NP_Q_NEW(K) -
     &        MEAN_MW_NP_Q(K)) /MEAN_MASS_NP(K)
            MEAN_MW_NP_QCL_NEW(K) = (MEAN_MW_NP_QCL_NEW(K) -
     &        MEAN_MW_NP_QCL(K)) /MEAN_MASS_NP(K)
            MEAN_MW_NP_QCF_NEW(K) = (MEAN_MW_NP_QCF_NEW(K) -
     &        MEAN_MW_NP_QCF(K)) /MEAN_MASS_NP(K)

            DO J=FIRST_ROW-1,FIRST_ROW  ! NP and adjacent row
              DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &             (J-1)*ROW_LENGTH+LAST_ROW_PT
                Q(I,K) = Q(I,K) - MEAN_MW_NP_Q_NEW(K)
                QCL(I,K) = QCL(I,K) - MEAN_MW_NP_QCL_NEW(K)
                QCF(I,K) = QCF(I,K) - MEAN_MW_NP_QCF_NEW(K)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
        
        IF (at_base_of_LPG) THEN
          MEAN_MW_SP_THETA_NEW(K) = (MEAN_MW_SP_THETA_NEW(K) -
     &      MEAN_MW_SP_THETA(K)) /MEAN_MASS_SP(K)
          DO J=P_LAST_ROW,P_LAST_ROW+1  ! SP and adjacent row
            DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &           (J-1)*ROW_LENGTH+LAST_ROW_PT
              THETA(I,K) = THETA(I,K) - MEAN_MW_SP_THETA_NEW(K)
            ENDDO
          ENDDO
          
          IF (K .LE. Q_LEVELS) THEN
            MEAN_MW_SP_Q_NEW(K) = (MEAN_MW_SP_Q_NEW(K) -
     &        MEAN_MW_SP_Q(K)) /MEAN_MASS_SP(K)
            MEAN_MW_SP_QCL_NEW(K) = (MEAN_MW_SP_QCL_NEW(K) -
     &        MEAN_MW_SP_QCL(K)) /MEAN_MASS_SP(K)
            MEAN_MW_SP_QCF_NEW(K) = (MEAN_MW_SP_QCF_NEW(K) -
     &        MEAN_MW_SP_QCF(K)) /MEAN_MASS_SP(K)

            DO J=P_LAST_ROW,P_LAST_ROW+1  ! SP and adjacent row
              DO I=(J-1)*ROW_LENGTH+FIRST_ROW_PT,
     &             (J-1)*ROW_LENGTH+LAST_ROW_PT
                Q(I,K) = Q(I,K) - MEAN_MW_SP_Q_NEW(K)
                QCL(I,K) = QCL(I,K) - MEAN_MW_SP_QCL_NEW(K)
                QCF(I,K) = QCF(I,K) - MEAN_MW_SP_QCF_NEW(K)
              ENDDO
            ENDDO
          ENDIF  !  is this a wet level
        ENDIF  !  at_base_of_LPG
      ENDDO  !  K : loop over levels

CL    END OF ROUTINE FILT_FLD

      RETURN
      END
