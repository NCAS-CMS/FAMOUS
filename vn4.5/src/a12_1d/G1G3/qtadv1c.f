C *****************************COPYRIGHT******************************
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
CLL   SUBROUTINE QT_ADV -------------------------------------------
CLL
CLL   PURPOSE:  CALCULATES MASS-WEIGHTED INCREMENTS TO QT
CLL             DUE TO ADVECTION  BY USING EQUATION (36)
CLL             TO CALCULATE PROVISIONAL VALUES OF QT AT
CLL             THE NEW TIME-LEVEL, AND THEN RECALCULATING THE
CLL             ADVECTION TERMS ON THE RIGHT-HAND SIDE OF  (36)
CLL             USING THESE PROVISIONAL VALUES. THE FINAL INCREMENTS ARE
CLL             CALCULATED AS IN EQUATION (40). THOSE REQUIRING
CLL             FILTERING ARE FILTERED, THE INCREMENTS
CLL             ARE ADDED ONTO THE FIELDS USING (40).
CLL             IF RUNNING A GLOBAL MODEL POLAR IS CALLED
CLL             TO UPDATE POLAR VALUES.
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL   VERSION FOR CRAY Y-MP
CLL
CLL   WRITTEN BY M.H MAWSON.
CLL   MPP CODE ADDED BY P.BURTON
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 4.1:
CLL VERSION  DATE
CLL 4.1      29/11/95 New version of routine specifically for MPP
CLL                   Fourth order MPP code by Roar Skalin
CLL                                                      P.Burton
!LL   4.2    16/08/96  Add TYPFLDPT arguments to FILTER subroutine
!LL                    and make the FILTER_WAVE_NUMBER arrays
!LL                    globally sized.                   P.Burton
!LL   4.2    10/01/97  Initialise unprocessed points in QT_PROV.
!LL                    D. Robinson.
!LL 4.3      24/04/97 Fixes to 4th order calculations   P.Burton
C     vn4.3    Mar. 97   T3E migration : optimisation changes
C                                       D.Salmond
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD B.
CLL
CLL   SYSTEM COMPONENTS COVERED: P121
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION:       THE EQUATIONS USED ARE (36) AND (40)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER NO. 10
CLL                        M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLLEND-------------------------------------------------------------

C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE QT_ADV
     1              (QT,PSTAR_OLD,PSTAR,U_MEAN,V_MEAN,
     2              SEC_P_LATITUDE,ETADOT_MEAN,RS,DELTA_AK,DELTA_BK,
     3              LATITUDE_STEP_INVERSE,ADVECTION_TIMESTEP,NU_BASIC,
     4              LONGITUDE_STEP_INVERSE,NORTHERN_FILTERED_P_ROW,
     5              SOUTHERN_FILTERED_P_ROW,Q_LEVELS,P_LEVELS,
     6              U_FIELD,P_FIELD,ROW_LENGTH,
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
     7              TRIGS,IFAX,FILTER_WAVE_NUMBER_P_ROWS,SEC_U_LATITUDE,
     8              AKH,BKH,L_SECOND,
     9              extended_address,
     &              LWHITBROM)                         

      IMPLICIT NONE

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

      INTEGER
     *  P_FIELD            !IN DIMENSION OF FIELDS ON PRESSSURE GRID.
     *, U_FIELD            !IN DIMENSION OF FIELDS ON VELOCITY GRID
     *, P_LEVELS           !IN NUMBER OF PRESSURE LEVELS.
     *, Q_LEVELS           !IN NUMBER OF MOIST LEVELS.
     *, ROW_LENGTH         !IN NUMBER OF POINTS PER ROW
     *, NORTHERN_FILTERED_P_ROW !IN ROW ON WHICH FILTERING STOPS
     *                          ! MOVING TOWARDS EQUATOR
     *, SOUTHERN_FILTERED_P_ROW !IN ROW ON WHICH FILTERING STARTS AGAIN
     *                          ! MOVING TOWARDS SOUTH POLE
     &, FILTER_WAVE_NUMBER_P_ROWS(GLOBAL_P_FIELD/GLOBAL_ROW_LENGTH)
     &       ! LAST WAVE NUMBER NOT TO BE CHOPPED
     *, IFAX(10)           !IN HOLDS FACTORS OF ROW_LENGTH USED BY
     *                     ! FILTERING.

C LOGICAL VARIABLE
      LOGICAL
     *  L_SECOND     ! SET TO TRUE IF NU_BASIC IS ZERO.
     & ,LWHITBROM    ! LOGICAL SWITCH FOR WHITE & BROMLEY
      INTEGER extended_address(P_FIELD)

      REAL
     * U_MEAN(U_FIELD,P_LEVELS) !IN AVERAGED MASS-WEIGHTED U VELOCITY
     *                          !   FROM ADJUSTMENT STEP
     *,V_MEAN(U_FIELD,P_LEVELS) !IN AVERAGED MASS-WEIGHTED V VELOCITY
     *                          !   * COS(LAT) FROM ADJUSTMENT STEP
     *,ETADOT_MEAN(P_FIELD,P_LEVELS)  !IN AVERAGED MASS-WEIGHTED
     *                          !VERTICAL VELOCITY FROM ADJUSTMENT STEP
     *,PSTAR(P_FIELD)           !IN PSTAR FIELD AT NEW TIME-LEVEL
     *,PSTAR_OLD(P_FIELD)       !IN PSTAR AT PREVIOUS TIME-LEVEL
     *,RS(P_FIELD,P_LEVELS)     !IN RS FIELD

      REAL
     * DELTA_AK(P_LEVELS)      !IN    LAYER THICKNESS
     *,DELTA_BK(P_LEVELS)      !IN    LAYER THICKNESS
     *,AKH(P_LEVELS+1)         !IN HYBRID CO-ORDINATE AT HALF LEVELS
     *,BKH(P_LEVELS+1)         !IN HYBRID CO-ORDINATE AT HALF LEVELS
     *,SEC_P_LATITUDE(P_FIELD) !IN  1/COS(LAT) AT P POINTS (2-D ARRAY)
     *,SEC_U_LATITUDE(U_FIELD) !IN  1/COS(LAT) AT U POINTS (2-D ARRAY)
     *,LONGITUDE_STEP_INVERSE  !IN 1/(DELTA LAMDA)
     *,LATITUDE_STEP_INVERSE   !IN 1/(DELTA PHI)
     *,ADVECTION_TIMESTEP      !IN
     *,NU_BASIC                !IN STANDARD NU TERM FOR MODEL RUN.
     *,TRIGS(ROW_LENGTH)       !IN HOLDS TRIGONOMETRIC FUNCTIONS USED
     *                         ! IN FILTERING.

      REAL
     * QT(P_FIELD,Q_LEVELS)    !INOUT QT FIELD.
     *                         ! MASS-WEIGHTED ON OUTPUT.

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: 23 ARE REQUIRED

      REAL
     * QT_FIRST_INC(P_FIELD,Q_LEVELS) ! HOLDS QT INCREMENT
     *                       ! RETURNED BY FIRST CALL TO ADV_P_GD
     *,QT_SECOND_INC(P_FIELD)! HOLDS QT INCREMENT
     *                       !RETURNED BY SECOND CALL TO ADV_P_GD
     *,QT_PROV(P_FIELD,Q_LEVELS)      ! HOLDS PROVISIONAL VALUE OF QT


      REAL
     * NUX(P_FIELD,Q_LEVELS)     ! COURANT NUMBER DEPENDENT NU AT P PTS
     *                   ! IN EAST-WEST ADVECTION.
     *,NUY(P_FIELD,Q_LEVELS)     ! COURANT NUMBER DEPENDENT NU AT P PTS
     *                   ! IN NORTH-SOUTH ADVECTION.

      REAL NUX_MIN(upd_P_ROWS),  ! minimum value of NUX along a row
     &     NUY_MIN(ROW_LENGTH)   ! min of NUY along a column

      REAL
     & ZERO(P_FIELD)              ! ARRAY OF ZEROES.
     *,QT_INCREMENT(P_FIELD,Q_LEVELS)

      REAL
     * BRSP(P_FIELD,Q_LEVELS) !MASS TERM AT LEVEL K

! Work space required to allow the use of Fourth Order Advection
! U/V_MEAN_COPY and Q_COPY arrays are defined with an extra halo
! this is required for the bigger stencil of the 4th order operator.
      REAL  U_MEAN_COPY((ROW_LENGTH+2*extra_EW_Halo)*
     &                  (tot_U_ROWS+2*extra_NS_Halo),P_LEVELS),
     &  !    Copy of U_MEAN with extra halo space for 4th order
     &      V_MEAN_COPY((ROW_LENGTH+2*extra_EW_Halo)*
     &                  (tot_U_ROWS+2*extra_NS_Halo),P_LEVELS),
     &  !    Copy of V_MEAN with extra halo space for 4th order
     &      Q_COPY((ROW_LENGTH+2*extra_EW_Halo)*
     &             (tot_P_ROWS+2*extra_NS_Halo),Q_LEVELS)
     &  !    Copy of QT with extra halo space for 4th order

      INTEGER  extended_P_FIELD,
     &         extended_U_FIELD
!  These are the sizes of the arrays with the extra halos

C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES
      INTEGER
     *  P_POINTS_UPDATE    ! NUMBER OF P POINTS TO BE UPDATED.
     *                     !  = ROWS*ROWLENGTH
     *, U_POINTS_UPDATE    ! NUMBER OF U POINTS TO BE UPDATED.
     *                     !  = (ROWS-1)*ROWLENGTH
     *, P_POINTS_REQUIRED  ! NUMBER OF P POINTS AT WHICH VALUES ARE
     *                     ! NEEDED TO UPDATE AT P_POINTS_UPDATE
     *, U_POINTS_REQUIRED  ! NUMBER OF U POINTS AT WHICH VALUES ARE
     *                     ! NEEDED TO UPDATE AT U_POINTS_UPDATE
     *, START_U_REQUIRED   ! FIRST U POINT OF VALUES REQUIRED TO UPDATE
     *                     ! AT P POINTS UPDATE.
     *, END_U_REQUIRED     ! LAST U POINT OF REQUIRED VALUES.

      INTEGER info  ! return code from comms

C REAL SCALARS
      REAL
     & SCALAR1,SCALAR2,TIMESTEP

C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     &  I,J,K1,IK,K
     *, FILTER_SPACE ! HORIZONTAL DIMENSION OF SPACE NEEDED IN FILTERING
     *               ! ROUTINE.
     &, I_start,I_end


C*L   EXTERNAL SUBROUTINE CALLS:---------------------------------------
      EXTERNAL ADV_P_GD,POLAR,UV_TO_P,FILTER
C*---------------------------------------------------------------------

CL    MAXIMUM VECTOR LENGTH ASSUMED IS P_FIELD.
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE INCLUDING SUBROUTINE CALLS:
CL---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL    SECTION 1.     INITIALISATION
CL---------------------------------------------------------------------
C INCLUDE LOCAL CONSTANTS FROM GENERAL CONSTANTS BLOCK

      P_POINTS_UPDATE   = upd_P_ROWS*ROW_LENGTH
      U_POINTS_UPDATE   = upd_U_ROWS*ROW_LENGTH
      P_POINTS_REQUIRED = (upd_P_ROWS+2)*ROW_LENGTH
      U_POINTS_REQUIRED = (upd_U_ROWS+2)*ROW_LENGTH
      START_U_REQUIRED  = START_POINT_NO_HALO-ROW_LENGTH
      END_U_REQUIRED    = END_U_POINT_NO_HALO+ROW_LENGTH

C *IF -DEF,NOWHBR replaced by LWHITBROM logical
      IF (LWHITBROM) THEN
CL    CALCULATE BRSP TERM AT LEVEL K

      K=1
! Loop over entire field
      DO I=FIRST_VALID_PT,LAST_P_VALID_PT
        BRSP(I,K)=(3.*RS(I,K)+RS(I,K+1))*(RS(I,K)-RS(I,K+1))
     *                *BKH(K+1)*.25*(PSTAR(I)-PSTAR_OLD(I))
      ENDDO
      K=Q_LEVELS
! Loop over entire field
      DO I=FIRST_VALID_PT,LAST_P_VALID_PT
        BRSP(I,K)=-(3.*RS(I,K)+RS(I,K-1))*(RS(I,K)-RS(I,K-1))
     *                *BKH(K)*.25*(PSTAR(I)-PSTAR_OLD(I))
      ENDDO

      DO K=2,Q_LEVELS -1
! Loop over entire field
        DO I=FIRST_VALID_PT,LAST_P_VALID_PT
          BRSP(I,K)=((3.*RS(I,K)+RS(I,K+1))*(RS(I,K)-RS(I,K+1))*BKH(K+1)
     *              *.25*(PSTAR(I)-PSTAR_OLD(I)))
     *              -((3.*RS(I,K)+RS(I,K-1))*(RS(I,K)-RS(I,K-1))*BKH(K)
     *              *.25*(PSTAR(I)-PSTAR_OLD(I)))
        ENDDO

      ENDDO
      END IF
C *ENDIF

! Loop over entire field
      DO I=FIRST_VALID_PT,LAST_P_VALID_PT
        ZERO(I) = 0.
      ENDDO

! In order to use the same call to adv_p_gd for both the second and
! fourth order advection, U/V_MEAN are copied into _COPY arrays.
! In the case of second order advection some of the work space is
! wasted as there is more halo than we need.

! Calculate the size of the extended arrays which contain an
! extra halo:
      extended_U_FIELD=(ROW_LENGTH+2*extra_EW_Halo)*
     &                 (tot_U_ROWS+2*extra_NS_Halo)
      extended_P_FIELD=(ROW_LENGTH+2*extra_EW_Halo)*
     &                 (tot_P_ROWS+2*extra_NS_Halo)

      IF (L_SECOND) THEN

! Copy U/V_MEAN to U/V_MEAN_COPY with the same sized halos
        CALL COPY_FIELD(U_MEAN,U_MEAN_COPY,
     &                  U_FIELD,extended_U_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  EW_Halo,NS_Halo,
     &                  .FALSE.)
        CALL COPY_FIELD(V_MEAN,V_MEAN_COPY,
     &                  U_FIELD,extended_U_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  EW_Halo,NS_Halo,
     &                  .FALSE.)

      ELSE  ! if its fourth order:

        CALL COPY_FIELD(U_MEAN,U_MEAN_COPY,
     &                  U_FIELD,extended_U_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  halo_4th,halo_4th,
     &                  .TRUE.)
        CALL COPY_FIELD(V_MEAN,V_MEAN_COPY,
     &                  U_FIELD,extended_U_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  halo_4th,halo_4th,
     &                  .TRUE.)
        CALL COPY_FIELD(QT,Q_COPY,
     &                  P_FIELD,extended_P_FIELD,
     &                  ROW_LENGTH,tot_P_ROWS,Q_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  halo_4th,halo_4th,
     &                  .TRUE.)

       ENDIF ! IF (L_SECOND)

CL LOOP OVER Q_LEVELS+1.
CL    ON 1 TO Q_LEVELS PROVISIONAL VALUES OF THE FIELD ARE CALCULATED.
CL    ON 2 TO Q_LEVELS+1 THE FINAL INCREMENTS ARE CALCULATED AND ADDED
CL    ON. THE REASON FOR THIS LOGIC IS THAT THE PROVISIONAL VALUE AT
CL    LEVEL K+1 IS NEEDED BEFORE THE FINAL INCREMENT AT LEVEL K CAN BE
CL    CALCULATED.

      DO K=1,Q_LEVELS+1

        TIMESTEP = ADVECTION_TIMESTEP

CL IF NOT AT Q_LEVELS+1 THEN
        IF(K.LE.Q_LEVELS) THEN

CL---------------------------------------------------------------------
CL    SECTION 2.     CALCULATE COURANT NUMBER DEPENDENT NU IF IN
CL                   FORECAST MODE. CALCULATE PROVISIONAL VALUES OF
CL                   QT AT NEW TIME-LEVEL.
CL---------------------------------------------------------------------

C ---------------------------------------------------------------------
CL    SECTION 2.1    SET NU TO NU_BASIC DEPENDENT ON MAX COURANT
CL                   NUMBER.
C ---------------------------------------------------------------------
CL    IF NU_BASIC NOT SET TO ZERO
          IF(.NOT.L_SECOND) THEN
CL    THEN SET NU DEPENDING ON NU_BASIC AND MAX
CL    COURANT NUMBER.
CL    CALCULATE COURANT NUMBER
C NOTE: RS AND TRIG TERMS WILL BE INCLUDED AFTER INTERPOLATION TO P
C       GRID.
CL    CALL UV_TO_P TO MOVE MEAN VELOCITIES ONTO P GRID

          CALL UV_TO_P(U_MEAN(START_U_REQUIRED,K),
     *                 NUX(START_POINT_NO_HALO,K),U_POINTS_REQUIRED,
     *                 P_POINTS_UPDATE,ROW_LENGTH,upd_P_ROWS+1)

          CALL UV_TO_P(V_MEAN(START_U_REQUIRED,K),
     *                 NUY(START_POINT_NO_HALO,K),U_POINTS_REQUIRED,
     *                 P_POINTS_UPDATE,ROW_LENGTH,upd_P_ROWS+1)

CL    CALCULATE NU FROM COURANT NUMBER INCLUDING TRIG AND RS TERMS.
          DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
            NUX(I,K) = NUX(I,K)*LONGITUDE_STEP_INVERSE
            NUY(I,K) = NUY(I,K)*LATITUDE_STEP_INVERSE
            SCALAR1 = TIMESTEP/(RS(I,K)*
     *                RS(I,K)*(DELTA_AK(K)+DELTA_BK(K)*PSTAR_OLD(I)))
            SCALAR2 = SEC_P_LATITUDE(I)*SCALAR1
            SCALAR1 = SCALAR1*SCALAR1
            SCALAR2 = SCALAR2*SCALAR2
            NUX(I,K) = (1. - NUX(I,K)*NUX(I,K)*SCALAR2)*NU_BASIC
            NUY(I,K) = (1. - NUY(I,K)*NUY(I,K)*SCALAR1)*NU_BASIC
          ENDDO

! Set NUX equal to minimum value along each row

          DO J=FIRST_ROW,FIRST_ROW+upd_P_ROWS-1
            I_start=(J-1)*ROW_LENGTH+FIRST_ROW_PT ! start and end of
            I_end=(J-1)*ROW_LENGTH+LAST_ROW_PT    ! row, missing halos.

! Calculate minimum along this row
            SCALAR1=NUX(I_start,K)
            DO I=I_start+1,I_end
              IF (NUX(I,K) .LT. SCALAR1) SCALAR1=NUX(I,K)
            ENDDO
            NUX_MIN(J-FIRST_ROW+1)=SCALAR1
! The indexing of NUX_MIN goes from 1..ROWS
          ENDDO ! J : loop over rows

! So far we have only calculated the minimum along our local
! part of the row. Now we must find the minimum of all the
! local minimums along the row
          CALL GCG_RMIN(upd_P_ROWS,GC_ROW_GROUP,info,NUX_MIN)

! and now set all values of NUX to the minimum along the row
          DO J=FIRST_ROW,FIRST_ROW+upd_P_ROWS-1
            IF (NUX_MIN(J-FIRST_ROW+1) .LT. 0.0)
     &        NUX_MIN(J-FIRST_ROW+1)=0.0

            I_start=(J-1)*ROW_LENGTH+1  ! beginning and
            I_end=J*ROW_LENGTH          ! end of row

            DO I=I_start,I_end
              NUX(I,K)=NUX_MIN(J-FIRST_ROW+1)
            ENDDO

          ENDDO ! J : loop over rows

! Set NUY equal to minimum value along each column

          DO J=FIRST_ROW_PT,LAST_ROW_PT
            I_start=(FIRST_ROW-1)*ROW_LENGTH+J
! I_start points to the beginning of column J

! Calculate the minimum along this column
            I_end=I_start+(upd_P_ROWS-1)*ROW_LENGTH
! I_end points to the end of column J
            SCALAR1=NUY(I_start,K)
            DO I=I_start+ROW_LENGTH,I_end,ROW_LENGTH
              IF (NUY(I,K) .LT. SCALAR1) SCALAR1=NUY(I,K)
            ENDDO
            NUY_MIN(J)=SCALAR1

          ENDDO ! J : loop over columns
! Once again, this is only the minimum along our local part
! of each column. We must now find the miniumum of all the local
! minimums along the column
          CALL GCG_RMIN(ROW_LENGTH-2*EW_Halo,GC_COL_GROUP,info,
     &                  NUY_MIN(EW_Halo+1))

! and now set all values of NUY to the minimum along the column
          DO J=FIRST_ROW_PT,LAST_ROW_PT
            IF (NUY_MIN(J) .LT. 0.0) NUY_MIN(J)=0.0

            I_start=(FIRST_ROW-1)*ROW_LENGTH+J
            I_end=I_start+(upd_P_ROWS-1)*ROW_LENGTH

            DO I=I_start,I_end,ROW_LENGTH
              NUY(I,K)=NUY_MIN(J)
            ENDDO

          ENDDO ! J : loop over columns

        ENDIF  ! IF its fourth order advection
CL
C ---------------------------------------------------------------------
CL    SECTION 2.3    CALL ADV_P_GD TO OBTAIN FIRST INCREMENT DUE TO
CL                   ADVECTION.
C ---------------------------------------------------------------------

CL    CALL ADV_P_GD FOR QT.
          K1=K+1

          IF(K.EQ.Q_LEVELS) THEN
          K1=K-1
          CALL ADV_P_GD(QT(1,K1),QT(1,K),QT(1,K1),
     &                  U_MEAN_COPY(1,K),V_MEAN_COPY(1,K),
     &                  ETADOT_MEAN(1,K),ZERO,SEC_P_LATITUDE,
     *                  QT_FIRST_INC(1,K),NUX(1,K),NUY(1,K),P_FIELD,
     *                  U_FIELD,ROW_LENGTH,
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
     &                  TIMESTEP,LATITUDE_STEP_INVERSE,
     *                  LONGITUDE_STEP_INVERSE,
     &                  SEC_U_LATITUDE,BRSP(1,K),L_SECOND,LWHITBROM,
     &                  Q_COPY(1,K),extended_P_FIELD,extended_U_FIELD,
     &                  extended_address)  
          ELSE IF(K.EQ.1)THEN

C PASS ANY QT VALUES FOR LEVEL K-1 AS ETADOT AT LEVEL 1
C IS SET TO ZERO BY USING ARRAY ZERO.

          CALL ADV_P_GD(QT(1,K1),QT(1,K),QT(1,K1),
     &                  U_MEAN_COPY(1,K),V_MEAN_COPY(1,K),
     &                  ZERO,ETADOT_MEAN(1,K1),
     *                  SEC_P_LATITUDE,QT_FIRST_INC(1,K),
     *                  NUX(1,K),NUY(1,K),
     *                  P_FIELD,U_FIELD,ROW_LENGTH,
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
     &                  TIMESTEP,
     *                  LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     &                  SEC_U_LATITUDE,BRSP(1,K),L_SECOND,LWHITBROM,
     &                  Q_COPY(1,K),extended_P_FIELD,extended_U_FIELD,
     &                  extended_address)  
          ELSE
          CALL ADV_P_GD(QT(1,K-1),QT(1,K),QT(1,K1),
     &                  U_MEAN_COPY(1,K),V_MEAN_COPY(1,K),
     &                  ETADOT_MEAN(1,K),ETADOT_MEAN(1,K1),
     *                  SEC_P_LATITUDE,QT_FIRST_INC(1,K),
     *                  NUX(1,K),NUY(1,K),
     *                  P_FIELD,U_FIELD,ROW_LENGTH,
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
     &                  TIMESTEP,
     *                  LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     &                  SEC_U_LATITUDE,BRSP(1,K),L_SECOND,LWHITBROM,
     &                  Q_COPY(1,K),extended_P_FIELD,extended_U_FIELD,
     &                  extended_address)  
          END IF

C ---------------------------------------------------------------------
CL    SECTION 2.4    REMOVE MASS-WEIGHTING FROM INCREMENT AND ADD ONTO
CL                   FIELD TO OBTAIN PROVISIONAL VALUE.
C ---------------------------------------------------------------------

          DO I=1,START_POINT_NO_HALO-1
            QT_PROV(I,K) = 0.0
          ENDDO

          DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
            SCALAR1 = RS(I,K)*RS(I,K)
     *                      *(DELTA_AK(K)+DELTA_BK(K)*PSTAR_OLD(I))
            QT_FIRST_INC(I,K) = QT_FIRST_INC(I,K)/SCALAR1
            QT_PROV(I,K) = QT(I,K)-QT_FIRST_INC(I,K)
          ENDDO

          DO I=END_P_POINT_NO_HALO+1,P_FIELD
            QT_PROV(I,K) = 0.0
          ENDDO

          IF (at_top_of_LPG) THEN
! Do North Pole
            DO I=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
              QT_PROV(I,K) = QT(I,K)
              QT_FIRST_INC(I,K) = -QT_FIRST_INC(I,K)/(RS(I,K)*RS(I,K)
     &                       *(DELTA_AK(K)+DELTA_BK(K)*PSTAR_OLD(I)))
            ENDDO
          ENDIF

          IF (at_base_of_LPG) THEN
! Do South Pole
            DO I=P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
              QT_PROV(I,K) = QT(I,K)
              QT_FIRST_INC(I,K) = -QT_FIRST_INC(I,K)/(RS(I,K)*RS(I,K)
     &                       *(DELTA_AK(K)+DELTA_BK(K)*PSTAR_OLD(I)))
            ENDDO
          ENDIF


        END IF
CL END CONDITIONAL ON LEVEL BEING LESS THAN Q_LEVELS+1
      ENDDO

CL    CALL POLAR TO OBTAIN PROVISIONAL VALUE.
      CALL POLAR(QT_PROV,QT_FIRST_INC,QT_FIRST_INC,
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
     &           P_FIELD,P_FIELD,P_FIELD,
     &           TOP_ROW_START,P_BOT_ROW_START,
     &           ROW_LENGTH,Q_LEVELS)

      IF (L_SECOND) THEN

! Do a halo update on QT_PROV
        CALL SWAPBOUNDS(QT_PROV,ROW_LENGTH,tot_P_ROWS,
     &                  EW_Halo,NS_Halo,Q_LEVELS)
!        CALL SET_SIDES(QT_PROV,P_FIELD,ROW_LENGTH,
!     &                 Q_LEVELS,fld_type_p)

      ELSE  ! fourth order advection

! Copy QT_PROV into Q_COPY which has double halos for fourth
! order advection, and do swap to fill these halos
        CALL COPY_FIELD(QT_PROV,Q_COPY,
     &                  P_FIELD,extended_P_FIELD,
     &                  ROW_LENGTH,tot_P_ROWS,Q_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  halo_4th,halo_4th,
     &                  .TRUE.)

      ENDIF

CL BEGIN CONDITIONAL ON LEVEL BEING GREATER THAN 1
cmic$ do parallel
      DO K=1,Q_LEVELS+1
        IF(K.GT.1) THEN
CL---------------------------------------------------------------------
CL    SECTION 3.     ALL WORK IN THIS SECTION PERFORMED AT LEVEL-1.
CL                   CALCULATE SECOND INCREMENT DUE TO ADVECTION.
CL                   CALCULATE TOTAL INCREMENT TO FIELD AND FILTER
CL                   WHERE NECESSARY THEN UPDATE FIELD.
CL                   THE POLAR INCREMENTS ARE THEN CALCULATED AND ADDED
CL                   ON BY CALLING POLAR.
CL---------------------------------------------------------------------

          TIMESTEP = ADVECTION_TIMESTEP

C ---------------------------------------------------------------------
CL    SECTION 3.1    CALL ADV_P_GD TO OBTAIN SECOND INCREMENT DUE TO
CL                   ADVECTION.
C ---------------------------------------------------------------------

          K1=K-1
C K1 HOLDS K-1.

          IF(K.GT.Q_LEVELS) THEN
C THE ZERO VERTICAL FLUX AT THE TOP IS ENSURED BY PASSING ETADOT AS
C ZERO.

          CALL ADV_P_GD(QT_PROV(1,K-2),QT_PROV(1,K-1),
     *                  QT_PROV(1,K-2),
     &                  U_MEAN_COPY(1,K1),V_MEAN_COPY(1,K1),
     &                  ETADOT_MEAN(1,K-1),
     *                  ZERO,SEC_P_LATITUDE,
     *                  QT_SECOND_INC,NUX(1,K-1),NUY(1,K-1),P_FIELD,
     *                  U_FIELD,ROW_LENGTH,
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
     &                  TIMESTEP,LATITUDE_STEP_INVERSE,
     *                  LONGITUDE_STEP_INVERSE,SEC_U_LATITUDE,
     &                  BRSP(1,K-1),L_SECOND,LWHITBROM,
     &                  Q_COPY(1,K-1),extended_P_FIELD,extended_U_FIELD,
     &                  extended_address)

          ELSE IF(K.EQ.2) THEN

C THE ZERO VERTICAL FLUX AT THE BOTTOM IS ENSURED BY PASSING ETADOT AS
C ZERO.

          CALL ADV_P_GD(QT_PROV(1,K),QT_PROV(1,K-1),
     *                  QT_PROV(1,K),
     &                  U_MEAN_COPY(1,K1),V_MEAN_COPY(1,K1),ZERO,
     *                  ETADOT_MEAN(1,K),
     *                 SEC_P_LATITUDE,QT_SECOND_INC,
     *                 NUX(1,K-1),NUY(1,K-1),
     *                  P_FIELD,U_FIELD,ROW_LENGTH,
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
     &                  TIMESTEP,
     *                  LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     *                  SEC_U_LATITUDE,
     &                  BRSP(1,K-1),L_SECOND,LWHITBROM,
     &                  Q_COPY(1,K-1),extended_P_FIELD,extended_U_FIELD,
     &                  extended_address)
          ELSE

          CALL ADV_P_GD(QT_PROV(1,K-2),QT_PROV(1,K-1),
     *                  QT_PROV(1,K),
     &                  U_MEAN_COPY(1,K1),V_MEAN_COPY(1,K1),
     &                  ETADOT_MEAN(1,K-1),
     *                  ETADOT_MEAN(1,K),
     *                 SEC_P_LATITUDE,QT_SECOND_INC,
     *                 NUX(1,K-1),NUY(1,K-1),
     *                  P_FIELD,U_FIELD,ROW_LENGTH,
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
     &                  TIMESTEP,
     *                  LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     *                  SEC_U_LATITUDE,
     &                  BRSP(1,K-1),L_SECOND,LWHITBROM,
     &                  Q_COPY(1,K-1),extended_P_FIELD,extended_U_FIELD,
     &                  extended_address)

          END IF

C ---------------------------------------------------------------------
CL    SECTION 3.2    CALCULATE TOTAL MASS-WEIGHTED INCREMENT TO FIELD.
C ---------------------------------------------------------------------

C TOTAL MASS-WEIGHTED INCREMENT IS CALCULATED AND THEN STORED IN
C QT_INCREMENT.

          DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
            QT_INCREMENT(I,K1) = .5*(QT_SECOND_INC(I) +
     *                      QT_FIRST_INC(I,K-1)*RS(I,K1)*RS(I,K1)
     *                      *(DELTA_AK(K1)+DELTA_BK(K1)*PSTAR(I)))
          ENDDO

C ---------------------------------------------------------------------
CL    SECTION 3.3    IF GLOBAL MODEL CALCULATE POLAR INCREMENTS.
CL                   IF LIMITED AREA MASS-WEIGHT BOUNDARY VALUES.
C ---------------------------------------------------------------------

CL    GLOBAL MODEL SO CALCULATE POLAR INCREMENT.
CL    CALCULATE MERIDIONAL FLUX AROUND POLES BY ADDING THE TWO
CL    INCREMENTS AND ALSO MASS-WEIGHTING POLAR FIELDS.
C NEGATIVE SIGN BEFORE FIRST INCS IS DUE TO THEIR SIGN BEING CHANGED
C PRIOR TO THE INTERMEDIATE VALUE BEING CALCULATED.
          IF (at_top_of_LPG) THEN
! Do Northen boundary/pole
            DO I=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
              SCALAR1 = RS(I,K1)*RS(I,K1)*
     &                  (DELTA_AK(K1)+DELTA_BK(K1)*PSTAR(I))
              QT_INCREMENT(I,K1) = -.5*(QT_SECOND_INC(I) -
     &                             QT_FIRST_INC(I,K-1)*SCALAR1)
              QT(I,K1)=QT(I,K1)*SCALAR1
            ENDDO
          ENDIF ! (attop)

          IF (at_base_of_LPG) THEN
! Do Southern boundary/pole
            DO IK=P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
              SCALAR2 = RS(IK,K1)*RS(IK,K1)*
     &                  (DELTA_AK(K1)+DELTA_BK(K1)*PSTAR(IK))
              QT_INCREMENT(IK,K1) = -.5*(QT_SECOND_INC(IK) -
     &                              QT_FIRST_INC(IK,K-1)*SCALAR2)
              QT(IK,K1)     = QT(IK,K1)*SCALAR2
            ENDDO
          ENDIF ! (atbase)

CL END CONDITIONAL LEVEL GREATER THAN ONE
        END IF

CL END LOOP OVER Q_LEVELS+1
      ENDDO

CL---------------------------------------------------------------------
CL    SECTION 4      IF GLOBAL MODEL THEN FILTER INCREMENTS AND
CL                   UPDATE POLAR VALUES BY CALLING POLAR.
CL                   UPDATE ALL OTHER VALUES.
CL---------------------------------------------------------------------

C ---------------------------------------------------------------------
CL    SECTION 4.1    CALL FILTER TO DO FILTERING.
C ---------------------------------------------------------------------

C SET FILTER_SPACE WHICH IS ROW_LENGTH+2 TIMES THE NUMBER OF ROWS TO
C BE FILTERED.

      FILTER_SPACE = (ROW_LENGTH+2)*(NORTHERN_FILTERED_P_ROW-1+
     *                tot_P_ROWS-SOUTHERN_FILTERED_P_ROW)
CL    CALL FILTER FOR QT INCREMENTS

      CALL FILTER(QT_INCREMENT,P_FIELD,Q_LEVELS,
     &            FILTER_SPACE,ROW_LENGTH,
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
     &            FILTER_WAVE_NUMBER_P_ROWS,TRIGS,IFAX,
     *            NORTHERN_FILTERED_P_ROW,SOUTHERN_FILTERED_P_ROW)

C ---------------------------------------------------------------------
CL    SECTION 4.2    CALL POLAR TO UPDATE POLAR VALUES
C ---------------------------------------------------------------------

      CALL POLAR(QT,QT_INCREMENT,QT_INCREMENT,
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
     &           P_FIELD,P_FIELD,P_FIELD,
     &           TOP_ROW_START,P_BOT_ROW_START,
     &           ROW_LENGTH,Q_LEVELS)

C ---------------------------------------------------------------------
CL    SECTION 4.3    UPDATE ALL OTHER POINTS.
C ---------------------------------------------------------------------

      DO K=1,Q_LEVELS
C UPDATE QT.
CFPP$ SELECT(CONCUR)
        DO I= START_POINT_NO_HALO,END_P_POINT_NO_HALO
          QT(I,K)=QT(I,K)*RS(I,K)*RS(I,K)*
     &            (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))-QT_INCREMENT(I,K)
        ENDDO
      ENDDO

CL    END OF ROUTINE QT_ADV

      RETURN
      END
