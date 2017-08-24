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
!+ Interface routine for QTPOS allowing both normal and MPP code to use
!+ current algorithms
!
! Subroutine Interface
      SUBROUTINE QT_POS_CTL
     &                     (QT,RS_SQUARED_DELTAP,ROW_LENGTH,
     &                      P_FIELD,Q_LEVELS,
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
     &                      ERROR_CODE,ERROR_MESSAGE,
     &                      COS_P_LATITUDE,SEC_P_LATITUDE,
     &                      L_NEG_QT,L_QT_POS_LOCAL,DT,SF_QTFIX,QTFIX)

      IMPLICIT NONE

!
! Description:
! QT fields are mass weighted, and then
! if MPP *DEF is set, the field is gathered a level to a processor, and
! QT_POS is called by each processor with its levels. Otherwise (ie. non
! MPP), this routine just passed the mass-weighted QT field straight
! through to QTPOS. Finally, the mass weighting is removed, and
! the QT_FIX diagnostic is calculated (if required)
!
! Current code owner (QT_POS_CTL only) : Paul Burton
!
! History
!  Model    Date      Modification history from model version 4.1
!    4.1    13/05/96  New subroutine added to allow QT_POS to be used
!                     by both normal and MPP code.   P.Burton
!    4.3    24/04/97  Improved error trapping for MPP code
!                     Set QTFIX to safe values         P.Burton
!    4.4    29/07/97  Optimize communications for T3E.  P.Burton
!    4.5    29/09/98  T3E only: remove negative zeros from QT
!                                                    P.Burton
!
! Subroutine Arguments:

      INTEGER
     &  P_FIELD                  ! IN : Size of fields on P grid
     &, ROW_LENGTH               ! IN : number of points per row
     &, Q_LEVELS                 ! IN : number of moist levels
!
     &, ERROR_CODE               ! INOUT : 0 on entry
!                                !         1 on exit if error detected

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

      CHARACTER*(80)
     &  ERROR_MESSAGE            ! OUT : Error message

      REAL
     &  QT(P_FIELD,Q_LEVELS)     ! INOUT : QT field
     &, QTFIX(P_FIELD,Q_LEVELS) ! OUT : QT fix diagnostic from
!                                !       resetting and redistributing
!                                !       negative QT values
     &, RS_SQUARED_DELTAP(P_FIELD,Q_LEVELS)
!                                ! IN : RS*RS*DELTA_P
     &, COS_P_LATITUDE(P_FIELD)  ! IN : COS(LATITUDE) at P points
     &, SEC_P_LATITUDE(P_FIELD)  ! IN : SEC(LATITUDE) at P points
     &, DT                       ! IN : dynamics timestep

      LOGICAL
     &  L_NEG_QT                 ! IN : Switch to continue in event of
!                                !      excessive negative QT
!                                !      (non-conservative)
     &, L_QT_POS_LOCAL           ! IN : perform local correction
!                                !      (method 2)
     &, SF_QTFIX                 ! IN : STASHflag for output of QT fix
!                                !      diagnostic

! Local variables

      INTEGER
     &  K,I                      ! loop variables
     &, I_start                  ! start address for QTPOS calcs
     &, I_end                    ! end address of QTPOS calcs
     &, I_start_local            ! local start address for QTPOS calcs
     &, I_end_local              ! local end address for QTPOS calcs
     &, info                     ! return code from comms.
     &, MAP(Q_LEVELS)            ! processor number for level
     &, N_LEVS_ON_PROC(0:N_PROCS-1) ! number of levels on each processor
     &, int_FAILURE(Q_LEVELS)    ! version of LOGICAL FAILURE array,
!                                ! but using integers

      REAL
     &  FACTOR(Q_LEVELS)         ! scaling factor for revised QT
     &, local_FACTOR((Q_LEVELS/N_PROCS)+1)  ! scaling factors for the
!                                          ! levels I process
     &,  global_QT_data(GLOBAL_P_FIELD,(Q_LEVELS/N_PROCS)+1)
! Array containing the QT levels to be processed on this processor
     &,  NP_VALUE(Q_LEVELS,2)    ! value of QT + QTFIX at North Pole
     &,  SP_VALUE(Q_LEVELS,2)    ! value of QT + QTFIX at South Pole

      LOGICAL
     &  FAILURE(Q_LEVELS)  ! Indicates if a particular level has
!                          ! failed to be QT_POSd.
     &, local_FAILURE((Q_LEVELS/N_PROCS)+1)  ! FAILURE array for
!                                           ! levels I process

C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------

!---------------------------------------------------------------------
! NOTE : The MPP code assumes, that for global models, the value of
!        latitude remains constant over a row (in particular the polar
!        rows). This implies the code will not give the correct answers
!        if a rotated global grid is used.

      I_start=1
      I_end=GLOBAL_P_FIELD
      I_start_local=FIRST_FLD_PT
      I_end_local=LAST_P_FLD_PT

      DO K=1,Q_LEVELS
        DO I=I_start_local,I_end_local
          QT(I,K)=QT(I,K)*RS_SQUARED_DELTAP(I,K)*COS_P_LATITUDE(I)
        ENDDO

        IF (SF_QTFIX) THEN
          DO I=FIRST_VALID_PT,LAST_P_VALID_PT
            QTFIX(I,K)=QT(I,K)
          ENDDO
          DO I=1,FIRST_VALID_PT-1
            QTFIX(I,K)=QT(FIRST_VALID_PT,K)
          ENDDO
          DO I=LAST_P_VALID_PT+1,P_FIELD
            QTFIX(I,K)=QT(LAST_P_VALID_PT,K)
          ENDDO
        ENDIF

! Caclulate the mapping of which processor each level will go to
        MAP(K)=MOD((K-1),N_PROCS)  ! assumes first PE is PE 0

      ENDDO

! Distribute QT over the processors
      DO K=0,N_PROCS-1
        N_LEVS_ON_PROC(K)=0
      ENDDO

      DO K=1,Q_LEVELS
        N_LEVS_ON_PROC(MAP(K))=N_LEVS_ON_PROC(MAP(K))+1

        CALL GATHER_FIELD(QT(1,K),
     &                    global_QT_data(1,N_LEVS_ON_PROC(MAP(K))),
     &                    ROW_LENGTH,tot_P_ROWS,
     &                    GLOBAL_ROW_LENGTH,
     &                    GLOBAL_P_FIELD/GLOBAL_ROW_LENGTH,
     &                    MAP(K),GC_ALL_GROUP,info)

      ENDDO

! and call QT_POS with the whole levels on this processor

      DO K=1,N_LEVS_ON_PROC(MY_PROC_ID)
        local_FAILURE(K)=.FALSE.
      ENDDO

      CALL QT_POS(global_QT_data,global_ROW_LENGTH,global_P_FIELD,
     &            I_start,I_end,
     &            N_LEVS_ON_PROC(MY_PROC_ID),
     &            ERROR_CODE,ERROR_MESSAGE,
     &            local_FACTOR,local_FAILURE,
     &            L_NEG_QT,L_QT_POS_LOCAL)

! Set ERROR_CODE and ERROR_MESSAGE if any pe has set ERROR_CODE
      CALL GC_IMAX(1,n_procs,info,ERROR_CODE)

      IF (ERROR_CODE .NE. 0) THEN
        ERROR_MESSAGE=
     & 'QT_POS   : MASS-WEIGHTED QT SUMMED OVER A LEVEL WAS NEGATIVE.'
      ENDIF

! and gather back levels to original processors

      DO K=0,N_PROCS-1
        N_LEVS_ON_PROC(K)=0
      ENDDO

      DO K=1,Q_LEVELS
        N_LEVS_ON_PROC(MAP(K))=N_LEVS_ON_PROC(MAP(K))+1

! collect the FAILURE and FACTOR array elements for this level

        IF (MY_PROC_ID .EQ. MAP(K)) THEN
!         I processed this level
          IF (local_FAILURE(N_LEVS_ON_PROC(MY_PROC_ID))) THEN
            int_FAILURE(K)=1
          ELSE
            int_FAILURE(K)=0
          ENDIF
          FACTOR(K)=local_FACTOR(N_LEVS_ON_PROC(MY_PROC_ID))
        ELSE  ! I didn't process this level
          int_FAILURE(K)=0
          FACTOR(K)=0.0
        ENDIF

! and scatter this level back to it's original home

        CALL SCATTER_FIELD(QT(1,K),
     &                    global_QT_data(1,N_LEVS_ON_PROC(MAP(K))),
     &                    ROW_LENGTH,tot_P_ROWS,
     &                    GLOBAL_ROW_LENGTH,
     &                    GLOBAL_P_FIELD/GLOBAL_ROW_LENGTH,
     &                    MAP(K),GC_ALL_GROUP,info)

      ENDDO

! and distribute the FAILURE and FACTOR arrays

      CALL GC_ISUM(Q_LEVELS,N_PROCS,info,int_FAILURE)
      CALL GC_RSUM(Q_LEVELS,N_PROCS,info,FACTOR)

      DO K=1,Q_LEVELS
        IF (int_FAILURE(K) .EQ. 0) THEN
          FAILURE(K)=.FALSE.
        ELSE
          FAILURE(K)=.TRUE.
        ENDIF
      ENDDO


! Now remove the mass weighting, and calculate QT_FIX diagnostic

      DO K=1,Q_LEVELS
        IF (FAILURE(K)) THEN
          IF (SF_QTFIX) THEN
            DO I=I_start_local,I_end_local
              QTFIX(I,K) = (QT(I,K) * FACTOR(K) - QTFIX(I,K))*
     &                     SEC_P_LATITUDE(I)/(A*A*G*DT)
            ENDDO
          ENDIF
          DO I=I_start_local,I_end_local
            QT(I,K)=QT(I,K)*FACTOR(K)/RS_SQUARED_DELTAP(I,K)*
     &              SEC_P_LATITUDE(I)
          ENDDO
        ELSE  ! FAILURE(K) .EQ. .FALSE.
          IF (SF_QTFIX) THEN
            DO I=I_start_local,I_end_local
              QTFIX(I,K) = (QT(I,K) - QTFIX(I,K)) *
     &                     SEC_P_LATITUDE(I)/(A*A*G*DT)
            ENDDO
          ENDIF
          DO I=I_start_local,I_end_local
            QT(I,K)=QT(I,K)/RS_SQUARED_DELTAP(I,K)*
     &              SEC_P_LATITUDE(I)
          ENDDO
        ENDIF

      ENDDO ! K : loop over levels


! Bring halos up to date
      CALL SWAPBOUNDS(QT,ROW_LENGTH,tot_P_ROWS,
     &                EW_Halo,NS_Halo,Q_LEVELS)

      RETURN
      END

CLL   SUBROUTINE QT_POS -------------------------------------------
CLL
CLL   PURPOSE:   REMOVES NEGATIVE VALUES OF QT.  THERE ARE TWO ALTERNAT-
CLL              IVE METHODS: METHOD 1 IS THE ORIGINAL SCHEME WHILE
CLL              IF L_QT_POS_LOCAL = .TRUE. METHOD 2 IS USED.
CLL              METHOD 2 IS DESIGNED TO ELIMINATE THE
CLL              VERY SLOW CLIMATE DRIFT IN QT IN UPPER MODEL LEVELS.
CLL              IT IS UNNECESSARILY COMPLICATED (AND EXPENSIVE) FOR
CLL              FORECAST USE.
CLL
CLL        METHOD 1:  ACCUMULATES TOTALS OF
CLL              OF NEGATIVE AND POSITIVE VALUES ON A LEVEL, ZEROING
CLL              ALL NEGATIVE POINTS AND PROPORTIONALLY REMOVING THE
CLL              SUM OF THE NEGATIVES FROM ALL POSITIVE POINTS.
CLL        METHOD 2:  THE FOLLOWING METHOD IS APPLIED AT EACH LAYER:
CLL            STEP 1 - LOOP ROUND ALL POINTS. IF QT IS NEGATIVE,
CLL              (a):
CLL              SUM QT ROUND NEIGHBOURING POINTS WITHIN ONE ROW AND
CLL              ONE COLUMN OF POINT IN CURRENT LAYER.  IF THIS VALUE
CLL              IS POSITIVE AND LARGE ENOUGH,THE CENTRE NEGATIVE VALUE
CLL              IS SET EQUAL TO ZERO AND THE NEIGHBOURS ARE SCALED
CLL              TO CONSERVE MASS WEIGHTED QT. IF (a) DOES NOT WORK
CLL              THEN (b):
CLL              EXTEND SUMMATION WITH A LARGER SEARCH RADIUS AND REPEAT
CLL              PROCESS. (SEARCH RADIUS = 4 FOUND TO BE SUFFICIENT)
CLL            STEP 2 - IF GLOBAL MODEL, CHECK FOR NEGATIVE QT AT POLES
CLL              AND CORRECT USING ALL VALUES IN NEIGHBOURING ROW.
CLL            STEP 3 - IF ANY NEGATIVE VALUES REMAIN (VERY UNCOMMON)
CLL              PERFORM GLOBAL CORRECTION AS IN METHOD 1.
CLL
CLL              IN BOTH METHODS IF THERE ARE INSUFFICIENT
CLL              POSITIVE VALUES WITHIN THE WHOLE LAYER AN ERROR
CLL              MESSAGE IS PRODUCED.  PROGRAM CONTINUES IF L_NEG_Q
CLL              IS TRUE IN WHICH CASE QT IS NOT CONSERVED.
CLL   N.B.       DEFINITION OF RS SQUARED DELTAP MEANS THAT POSITIVE
CLL              QT VALUES HAVE A NEGATIVE SIGN!
CLL
CLL
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL   VERSION FOR CRAY Y-MP
CLL
CLL  MM, SB     <- PROGRAMMER OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL   3.4    02/08/94 Add QTFIX diagnostic.  Tim Johns
CLL   3.4     7/08/94 Method 2 added to perform local corrections.
CLL                   Author Chris Hall.
!LL   4.0     9/02/95 : Alter code for local correction option to help
!LL                     speed it up. Still scalar in parts. R A Stratton
CLL   4.0     1/02/95 Values of polar weighting changed for consistency
CLL                   with other parts of dynamics.  Chris Hall
!LL   4.1    14/05/96 Moved some code up to QT_POS_CTL, and changed
!LL                   arguments/comments.    P.Burton
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD A. VERSION 2, DATED 18/01/90
CLL
CLL   SYSTEM COMPONENTS COVERED: P 191
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION:       SECTION 3.7
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER
CLL                        NO. 10 M.J.P.CULLEN,T.DAVIES AND M.H.MAWSON
CLL                        VERSION 10, DATED 10/09/90.
CLLEND-------------------------------------------------------------

C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE QT_POS
     1                 (QT,ROW_LENGTH,P_FIELD,ISTART,IEND,Q_LEVELS,
     2                  ERROR_CODE,ERROR_MESSAGE,FACTOR,FAILURE,
     3                  L_NEG_QT,L_QT_POS_LOCAL)

      IMPLICIT NONE

      INTEGER
     &  ROW_LENGTH         ! IN : number of points per row
     &, P_FIELD            ! IN : size of horizontal P field
     &, ISTART             ! IN : first point to process
     &, IEND               ! IN : last point to process
     &, Q_LEVELS           ! IN : number of (moist) levels to process

      LOGICAL
     &  L_NEG_QT           ! IN : switch to continue in event of
!                          !      excessive negative QTs (non-conserv.)
     &, L_QT_POS_LOCAL     ! IN : perform local corrections (method 2)

      REAL
     &  QT(P_FIELD,Q_LEVELS) ! IN/OUT : QT field to update

      INTEGER
     &  ERROR_CODE         ! OUT : 1 on exit if error detected

      CHARACTER*(80)
     &  ERROR_MESSAGE      ! OUT : Error message if ERROR_CODE not zero

      REAL
     &  FACTOR(Q_LEVELS)   ! OUT : Scaling factor for QT

      LOGICAL
     &  FAILURE(Q_LEVELS)  ! OUT : Indicates a level needs rescaling

C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: NONE ARE REQUIRED
C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES

      INTEGER
     *  I,K,II,JJ,NN,J   ! Loop variables
     *, IN,IS,JE,JW        ! Row/column indices to north/south/east/west
     *, POINTER            ! Array pointer
     *, ROW                ! Row number counting from 0
     *, NPT                ! Point in row counting from 1
     *, ROWS               ! Number of rows
     *, MAX_SEARCH         ! Maximum search radius (points)
     *, NO_NEG             ! Count of number negative in step 3 & 1
     *, NO_NEG2            ! Count of number negative
     *, NINDEX(P_FIELD)    ! locations of negative q
      REAL
     *  SUM_NEIGHBOURS     ! Sum of QT around neighbours
     *, SUM_POSITIVE       ! Global sum of positive QT
     *, SUM_NEGATIVE       ! Global sum of negative QT
     *, TEMP1,TEMP2        ! Temporary sums

C*L   NO EXTERNAL SUBROUTINE CALLS:------------------------------------
C*---------------------------------------------------------------------
CL    MAXIMUM VECTOR LENGTH ASSUMED IS P_FIELD.
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL
      ROWS=P_FIELD/ROW_LENGTH
      MAX_SEARCH=4

CL LOOP OVER LEVELS
      DO K= 1,Q_LEVELS
        FAILURE(K)=.FALSE.
        FACTOR(K)=0.0
        IF (L_QT_POS_LOCAL) THEN
CL---------------------------------------------------------------------
CL    METHOD 2
CL---------------------------------------------------------------------
          NO_NEG=0

          DO I=1,P_FIELD
            IF(QT(I,K).GT.0.0) THEN
              NO_NEG=NO_NEG+1
              NINDEX(NO_NEG)=I
            ENDIF
          ENDDO

          DO J=1,NO_NEG
            I=NINDEX(J)
CL
CL---------------------------------------------------------------------
CL    SECTION 2.    STEP 1. WHERE QT IS NEGATIVE LOOP ROUND
CL              8 NEIGHBOURS TO PERFORM LOCAL CORRECTION.  THIS OCCURS
CL              AT APPROXIMATELY 2% OF THE POINTS.  THE FOLLOWING CODE
CL              IN SECTIONS 2 & 3 IS RATHER LONG-WINDED, BUT CPU IS
CL              SIGNIFICANTLY REDUCED OVER A SIMPLER APPROACH.
CL---------------------------------------------------------------------
C   If QT is negative find 8 neighbours ...
              ROW=(I-1)/ROW_LENGTH
              NPT=I-ROW*ROW_LENGTH
              JE=1
              JW=-1
C   Abandon search for neighbours on LAM boundary
              IF(NPT.EQ.ROW_LENGTH.OR.NPT.EQ.1) THEN
                FAILURE(K)=.TRUE.
                GOTO 200
              END IF
C   ... and sum neighbouring 8 values
              SUM_NEIGHBOURS = QT(I+JE,K) + QT(I+JW,K)
              IF(ROW.GT.0) THEN
                II=I-ROW_LENGTH
                SUM_NEIGHBOURS = SUM_NEIGHBOURS +
     *                           QT(II,K) + QT(II+JE,K) + QT(II+JW,K)
              END IF
              IF(ROW.LT.ROWS-1) THEN
                II=I+ROW_LENGTH
                SUM_NEIGHBOURS = SUM_NEIGHBOURS +
     *                           QT(II,K) + QT(II+JE,K) + QT(II+JW,K)
              END IF
C
C   If possible set QT=0 at point with negative value and adjust
C   neighbouring values by FACTOR to conserve QT
              IF(SUM_NEIGHBOURS.LT.0.0) THEN
                FACTOR(K)=1.0+QT(I,K)/SUM_NEIGHBOURS
                IF(FACTOR(K).GE.0.0) THEN
                  QT(I,K)=0.0
                  QT(I+JE,K) = QT(I+JE,K)*FACTOR(K)
                  QT(I+JW,K) = QT(I+JW,K)*FACTOR(K)
                  IF(ROW.GT.0) THEN
                    II=I-ROW_LENGTH
                    QT(II,K) = QT(II,K)*FACTOR(K)
                    QT(II+JE,K) = QT(II+JE,K)*FACTOR(K)
                    QT(II+JW,K) = QT(II+JW,K)*FACTOR(K)
                  END IF
                  IF(ROW.LT.ROWS-1) THEN
                    II=I+ROW_LENGTH
                    QT(II,K) = QT(II,K)*FACTOR(K)
                    QT(II+JE,K) = QT(II+JE,K)*FACTOR(K)
                    QT(II+JW,K) = QT(II+JW,K)*FACTOR(K)
                  END IF
                END IF
              END IF
  200       CONTINUE       ! LAM only
          END DO       ! loop over negative grid points

! recalculate number of negative gridpoints
          NO_NEG2=0
          DO J=1,NO_NEG
            I=NINDEX(J)
            IF(QT(I,K).GT.0.0) THEN
              NO_NEG2=NO_NEG2+1
              NINDEX(NO_NEG2)=I
            ENDIF
          ENDDO
!         WRITE(6,*)' Step 1 level ',k,' no neg ',no_NEG
!         WRITE(6,*)' Step 2 level ',k,' no neg ',no_NEG2
CL
CL---------------------------------------------------------------------
CL    SECTION 3.    STEP 2. WHERE QT IS NEGATIVE AND STEP 1 CORRECTIONS
CL              HAVE FAILED, ATTEMPT LOCAL CORRECTION BY EXTENDING
CL              SEARCH TO 24, 48 OR 80 SURROUNDING POINTS (SEARCH RADIUS
CL              EQUALS 2, 3 OR 4), OR FURTHER.  THIS OCCURS ABOUT
CL              0.1% OF THE TIME AT CLIMATE RESOLUTION AND MUCH MORE
CL              FREQUENTLY AT FORECAST RESOLUTION.  ALTHOUGH SMALL IN
CL              NUMBER, LOCAL CORRECTIONS AT THESE POINTS ARE REQUIRED
CL              TO SLOW DOWN THE LOSS OF STRATOSPHERIC MOISTURE THAT
CL              OCCURS ON CLIMATE TIME SCALES.
CL---------------------------------------------------------------------
         IF (NO_NEG2.GT.0) THEN         ! only do this if required
           DO J=1,NO_NEG2
            I=NINDEX(J)
              ROW=(I-1)/ROW_LENGTH
              NPT=I-ROW*ROW_LENGTH
C   Abandon search for neighbours on LAM boundary
              IF(NPT.EQ.ROW_LENGTH.OR.NPT.EQ.1) THEN
                FAILURE(K)=.TRUE.
                GOTO 300
              END IF
              DO NN=2,MAX_SEARCH
                JW=NPT-NN
                JE=NPT+NN
                IF(JE.GT.ROW_LENGTH) JE=ROW_LENGTH
                IF(JW.LT.1) JW=1
                IN=ROW-NN
                IS=ROW+NN
                IF(IN.LT.0) IN=0
                IF(IS.GE.ROWS) IS=ROWS-1
C   ... and sum QT at neighbouring points
                SUM_NEIGHBOURS=-QT(I,K)
                DO II=IN,IS
                  DO JJ=JW,JE
                    POINTER = II*ROW_LENGTH + JJ
                    SUM_NEIGHBOURS = SUM_NEIGHBOURS + QT(POINTER,K)
                  END DO
                END DO
C
C   If possible set QT=0 at point with negative value and adjust
C   neighbouring values by FACTOR to conserve QT
                IF(SUM_NEIGHBOURS.LT.0.0) THEN
                  FACTOR(K)=1.0+QT(I,K)/SUM_NEIGHBOURS
                  IF(FACTOR(K).GE.0.0) THEN
                    QT(I,K)=0.0
                    DO II=IN,IS
                      DO JJ=JW,JE
                        POINTER = II*ROW_LENGTH + JJ
                        QT(POINTER,K) = QT(POINTER,K)*FACTOR(K)
                      END DO
                    END DO
                    GOTO 300    ! Successfull correction performed
                  END IF
                END IF
              END DO     ! End loop on search radius (NN)
C
C  No correction possible within maximum search radius (=4)
              FAILURE(K)=.TRUE.
  300       CONTINUE
          END DO         ! End loop on points (I)
         END IF          ! end step 2

        END IF          !  End main section for method 1
CL---------------------------------------------------------------------
CL    SECTION 5.     METHOD 1.  ALSO METHOD 2 STEP 3.
CL              REMOVE REMAINING NEGATIVE QT BY SUMMING
CL              NEGATIVE AND POSITIVE VALUES SEPARATELY ON THE
CL              LEVEL AND PERFORMING A GLOBAL CORRECTION.
CL---------------------------------------------------------------------
        IF(FAILURE(K).OR..NOT.L_QT_POS_LOCAL) THEN
          NO_NEG=0            !zero count of -ve points
          SUM_POSITIVE = 0.0
          SUM_NEGATIVE = 0.0
          DO I=ISTART,IEND
            IF(QT(I,K).GT.0.)THEN
              NO_NEG=NO_NEG+1           !count number of -ve points
              SUM_NEGATIVE = SUM_NEGATIVE + QT(I,K)
              QT(I,K) = 0.0
            ELSE
              SUM_POSITIVE = SUM_POSITIVE + QT(I,K)
            END IF
          END DO
!         WRITE(6,*)' GLOBAL correction level ',k,' no neg ',no_NEG

CL If a negative value is found re-scale all positive points
CL
          IF (SUM_NEGATIVE.GT.0.0) THEN
            FAILURE(K)=.TRUE.
            FACTOR(K) = 1. + SUM_NEGATIVE/SUM_POSITIVE
            IF(L_NEG_QT.AND.FACTOR(K).LT.0)THEN
CL Skip abort if L_NEG_QT is true
              WRITE(6,*)' QT_POS : Mass weighted QT summed over level ',
     *            K,' was negative.  WARNING: QT not conserved'
              WRITE(6,*) NO_NEG,' points were -ve and the scaling ',
     *           'factor has been reset to 1'
              FACTOR(K) = 1.0
            ENDIF
            IF(FACTOR(K).LT.0.0) THEN
              ERROR_CODE = 1
              ERROR_MESSAGE=
     *  'QT_POS   : MASS-WEIGHTED QT SUMMED OVER A LEVEL WAS NEGATIVE.'
              RETURN
            END IF
          ELSE
            FAILURE(K)=.FALSE.
          END IF
        END IF

CL END LOOP OVER LEVELS
      END DO


CL    END OF ROUTINE QT_POS

      RETURN
      END
