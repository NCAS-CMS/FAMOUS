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
CLL  SUBROUTINE WINDMAX-------------------------------------------------
CLL
CLL  PURPOSE: CALCULATES THE MAXIMUM WIND SPEED THAT LIES BETWEEN
CLL           70000 AND 10000 Pa.(ALSO CONSTRAINED TO A PARTICULAR SET
CLL           OF MODEL LEVELS BY THE ORDER OF THE MATRIX POLYNOMIAL)
CLL  A QUARTIC CURVE IS FITTED TO THE ETA-LEVEL WIND SPEEDS,TREATED
CLL  AS LAYER MEANS CENTRED ON THE ETA-LEVELS,IN THE VICINITY
CLL  OF THE MAXIMUM ETA-LEVEL WIND.
CLL  THE MAXIMUM THAT THIS CURVE ATTAINS IN THE LAYER CENTRED ON THE
CLL  APPROPRIATE ETA-LEVEL GIVES THE MAXIMUM WIND SPEED AND LEVEL.
CLL  MAXIMUM WIND DIRECTION IS FOUND BY LINEAR INTERPOLATION (IN LOG(P))
CLL  FROM WIND DIRECTIONS AT SURROUNDING ETA-LEVELS.
CLL  THE PRESSURE OF THE MAXIMUM WIND IS CALCULATED FROM THE MAXIMUM
CLL  WIND ETA-LEVEL.
CLL
CLL  NOT SUITABLE FOR SINGLE COLUMN USE
CLL
CLL J.HEMING    <- PROGRAMMER OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL   4.3    17/03/96 MPP SWAPBOUNDS added. And fixed a bug where
CLL                   one loop had inappropriately been split into 2
CLL                   S.D. Mullerworth
!LL   4.5    17/04/98 Removed SWAPBOUNDS and ignore NS halos in loops
!LL                   WARNING - NS Halos need to be initialised by
!LL                   calling routine.
!LL                   Includes minor correction to (7) that changes
!LL                   bit comparison of diagnostic. S.D.Mullerworth
CLL                   
CLL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL  VERSION 2, DATED 18/01/90
CLL
CLL  SYSTEM TASK: D4
CLL
CLL  SYSTEM COMPONENT: D412
CLL
CLL  DOCUMENTATION:  UMDP NO.80
CLL
CLLEND------------------------------------------------------------------
C
C*L ARGUMENTS:----------------------------------------------------------
      SUBROUTINE WINDMAX
     1 (
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
     &  PSTAR,U,V,
     2  U_ROWS,P_ROWS,ROW_LENGTH,P_LEVELS,P_FIELD,U_FIELD,
     3  AK,BK,AKH,BKH,ETA_MATRIX_INV,MATRIX_P_O,
     4  U_MAXWIND,V_MAXWIND,PRESSURE_MAXWIND)
C
      IMPLICIT NONE
C
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
     *  U_ROWS       !IN   NUMBER OF ROWS FOR (U,V) FIELD
     *, P_ROWS       !IN   NUMBER OF ROWS FOR (P,T) FIELD
     *, ROW_LENGTH   !IN   NO OF POINTS PER ROW
     *, P_LEVELS     !IN   NO OF PRESSURE LEVELS
     *, P_FIELD      !IN   FIRST DIMENSION OF FIELD OF PSTAR
     *, U_FIELD      !IN   FIRST DIMENSION OF (U,V) FIELD
     *, MATRIX_P_O   !IN   ORDER OF POLYNOMIAL USED IN ETA MATRIX
C
      REAL
     *  PSTAR(P_FIELD)   !IN  PRESSURE AT GROUND LEVEL AT (P,T) POINTS
     *, U(U_FIELD,P_LEVELS) !IN  EASTERLY WIND COMPONENT IN M/S
     *, V(U_FIELD,P_LEVELS) !IN  NORTHERLY WIND COMPONENT IN M/S
C
C       AK,BK DEFINE HYBRID VERTICAL COORDINATES P=A+BP*
      REAL
     *  AK(P_LEVELS)     !IN  VALUE AT LAYER CENTRE
     *, BK(P_LEVELS)     !IN  VALUE AT LAYER CENTRE
     *, AKH(P_LEVELS+1)  !IN  LAYER THICKNESS
     *, BKH(P_LEVELS+1)  !IN  LAYER THICKNESS
     *, ETA_MATRIX_INV(MATRIX_P_O,MATRIX_P_O,P_LEVELS)
     *                   !IN INVERSE MATRIX OF ETA_HALF VALUES
C
      REAL
     *  U_MAXWIND(U_FIELD) !OUT EASTERLY COMPONENT OF MAXWIND IN M/S
     *, V_MAXWIND(U_FIELD) !OUT NORTHERLY COMPONENT OF MAXWIND IN M/S
     *, PRESSURE_MAXWIND(U_FIELD)!OUT PRESSURE AT LEVEL OF MAXIMUM WIND
C*----------------------------------------------------------------------
C
C*L WORKSPACE USAGE-----------------------------------------------------
C*----------------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED-----------------------------------------
      EXTERNAL P_TO_UV
C*----------------------------------------------------------------------
C
C*L------------------COMDECK C_PI---------------------------------------
CLL
CLL 4.0 19/09/95  New value for PI. Old value incorrect
CLL               from 12th decimal place. D. Robinson
CLL
      REAL PI,PI_OVER_180,RECIP_PI_OVER_180

      PARAMETER(
     & PI=3.14159265358979323846, ! Pi
     & PI_OVER_180 =PI/180.0,   ! Conversion factor degrees to radians
     & RECIP_PI_OVER_180 = 180.0/PI ! Conversion factor radians to
     &                              ! degrees
     & )
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

C
C-----------------------------------------------------------------------
C DEFINE LOCAL CONSTANTS
C-----------------------------------------------------------------------
      INTEGER
     *  N21      ! =(MATRIX_P_O+1)/2
C
      REAL
     *  TWOPI    ! = 2 * PI
     *, TLEV     ! TOP LEVEL IN Pa AT WHICH MAXWIND ALLOWED
     *, BLEV     ! BOTTOM LEVEL IN Pa AT WHICH MAXWIND ALLOWED
C
      PARAMETER(TWOPI=2*PI,
     *          TLEV=10000.0,BLEV=70000.0)
C-----------------------------------------------------------------------
C DEFINE LOCAL VARIABLES
C-----------------------------------------------------------------------
      INTEGER
     *  KL
     *, I,J,K        !  LOOP COUNTERS
     *, NO           !  GRIDPOINT NUMBER
     *, NI           !  = NO + ROW_LENGTH
     *, POINTS       !  NO OF POINTS BETWEEN TWO LEVELS WHICH
     *               !  ARE CHECKED FOR MAXIMUM WIND SPEED.
     *, LEVEL_MAXWIND(U_FIELD)! LEVEL AT WHICH MAXIMUM WIND OCCURS
C
      REAL
     *  ETA_LEVEL(P_LEVELS) ! ETA LEVELS
     *, ETA_HALF(P_LEVELS+1)! ETA HALF-LEVELS
     *, ETA_JMH             ! ETA OF J-1/2
     *, ETA_JPH             ! ETA OF J+1/2
     *, PSTAR_BAR(U_FIELD)  ! P* AT (U,V) POINTS
     *, WINDSPEED(U_FIELD,P_LEVELS) ! WINDSPEED AT GRIDPOINTS
     *, DIRECTION(U_FIELD,P_LEVELS) ! WIND DIRECTION AT GRIDPOINTS
     *, A(U_FIELD,MATRIX_P_O)      ! COEFFS OF QUARTIC POLYNOMIAL
     *, DIRECTION_MAXWIND(U_FIELD) ! DIRECTION OF MAXIMUM WIND
     *, SPEED_MAXWIND(U_FIELD)    ! SPEED OF MAXIMUM WIND
     *, ETA_MAXWIND(U_FIELD)      ! VALUE OF ETA AT MAXWIND LEVEL
C
      REAL
     *  INVPTS          ! 1/POINTS
     *, X
     *, XINC            ! INCREMENT IN X
     *, F(U_FIELD)      ! A(1)+...A(N)*ETA**(N-1)
     *, DIFF   ! DIFFERENCE IN DIRECTION BETWEEN ETA LEVELS EACH
     *         ! SIDE OF MAXIMUM WIND LEVEL
     *, ETA_J           ! ETA_LEVEL(J)    ) WHERE J=LEVEL_MAXWIND
     *, ETA_JP1         ! ETA_LEVEL(J+1)  )
     *, MULT
     *, P_JMH           ! P AT K=J-1/2 )  P=A(K)+B(K)*PSTAR_BAR
     *, P_JPH           ! P AT K=J+1/2 )
C
      LOGICAL
     *  S1         ! )
     *, S2         ! )  3 LOGICAL STATEMENTS
     *, S3         ! )
C-----------------------------------------------------------------------
CL  1. CALCULATE VALUES OF ETA ON FULL AND HALF LEVELS
C-----------------------------------------------------------------------
      DO K=1,P_LEVELS
        ETA_LEVEL(K)=AK(K)/PREF+BK(K)
        ETA_HALF(K)=AKH(K)/PREF+BKH(K)
      ENDDO
      ETA_HALF(P_LEVELS+1)=AKH(P_LEVELS+1)/PREF+BKH(P_LEVELS+1)
C-----------------------------------------------------------------------
CL  2. P* IS INTERPOLATED FROM THETA POINTS TO (U,V) POINTS
C-----------------------------------------------------------------------
      CALL P_TO_UV(PSTAR,PSTAR_BAR,P_FIELD,U_FIELD,ROW_LENGTH,P_ROWS)
C-----------------------------------------------------------------------
CL  3. WINDSPEED AND WIND DIRECTION ARE CALCULATED AT EACH GRIDPOINT
CL     AT EACH LEVEL
C-----------------------------------------------------------------------
      DO 6 K=1,P_LEVELS
        DO NO=FIRST_FLD_PT,LAST_U_FLD_PT
          WINDSPEED(NO,K)=SQRT(U(NO,K)**2+V(NO,K)**2)
          IF((U(NO,K).EQ.0.).AND.(V(NO,K).EQ.0.))THEN
            DIRECTION(NO,K)=0.
          ELSE
          DIRECTION(NO,K)=ATAN2(U(NO,K),V(NO,K))
          ENDIF
C-----------------------------------------------------------------------
C                - DIRECTION IN RADIANS FROM N
C                   (+VE BETWEEN 0 AND PI, -VE BETWEEN -PI AND 0)
C-----------------------------------------------------------------------
        ENDDO
 6    CONTINUE
C-----------------------------------------------------------------------
CL  4. FIND THE MAXIMUM ETA-LEVEL WIND SPEED BETWEEN 100mb and 700mb
CL     AND WHICH ALSO LIES BETWEEN LEVELS N21 AND P_LEVELS-N21+1
C-----------------------------------------------------------------------
      N21=(MATRIX_P_O+1)/2
      DO NO=FIRST_FLD_PT,LAST_U_FLD_PT
        LEVEL_MAXWIND(NO)=1
        SPEED_MAXWIND(NO)=0.
      ENDDO
      DO 20 I=N21,P_LEVELS-N21+1
        DO NO=FIRST_FLD_PT,LAST_U_FLD_PT
C-----------------------------------------------------------------------
C   CHECK THAT ETA LEVEL IS BELOW TLEV AND ABOVE BLEV
C-----------------------------------------------------------------------
          S1=WINDSPEED(NO,I).GE.SPEED_MAXWIND(NO)
          S2=TLEV.LT.(AK(I+1)+BK(I+1)*PSTAR_BAR(NO))
          S3=BLEV.GE.(AKH(I)+BKH(I)*PSTAR_BAR(NO))
          IF ((S1.AND.S2).AND.S3) THEN
C-----------------------------------------------------------------------
C   STORE NEW MAXIMUM SPEED
C-----------------------------------------------------------------------
            SPEED_MAXWIND(NO)=WINDSPEED(NO,I)
C-----------------------------------------------------------------------
C   STORE NEW MAXIMUM DIRECTION
C-----------------------------------------------------------------------
            DIRECTION_MAXWIND(NO)=DIRECTION(NO,I)
C-----------------------------------------------------------------------
C   STORE NEW MAXIMUM ETA
C-----------------------------------------------------------------------
            ETA_MAXWIND(NO)=ETA_LEVEL(I)
C-----------------------------------------------------------------------
C   STORE LEVEL
C-----------------------------------------------------------------------
            LEVEL_MAXWIND(NO)=I
          ENDIF
        ENDDO
  20  CONTINUE
      DO NO=FIRST_FLD_PT,LAST_U_FLD_PT
        IF (LEVEL_MAXWIND(NO).EQ.1) THEN
          DIRECTION_MAXWIND(NO)=0.0
          ETA_MAXWIND(NO)=0.0
        ENDIF
      ENDDO
C------------------------------------------------------------------
CL  5. FIND THE COEFFICIENTS A OF
CL          F(ETA)=A(1)+A(2)*ETA+.......+A(N)*ETA**(N-1)
C------------------------------------------------------------------
      DO 206 J=1,MATRIX_P_O
        DO NO=FIRST_FLD_PT,LAST_U_FLD_PT
          A(NO,J)=0.0
        ENDDO
 206  CONTINUE
      DO 207 J=1,MATRIX_P_O
        DO I=1,MATRIX_P_O
          DO NO=FIRST_FLD_PT,LAST_U_FLD_PT
            A(NO,J)=A(NO,J)+ETA_MATRIX_INV(I,J,LEVEL_MAXWIND(NO))
     *      *WINDSPEED(NO,LEVEL_MAXWIND(NO)-N21+I)
          ENDDO
        ENDDO
 207  CONTINUE
C-----------------------------------------------------------------------
CL  6. CALCULATE FROM THE CURVE THE SPEED OF THE MAXIMUM WIND AND
CL     THE VALUE OF ETA AT WHICH IT OCCURS
C-----------------------------------------------------------------------
C   FIND THE POINT WHICH GIVES THE MAXIMUM WINDSPEED BETWEEN
C   HALF LEVELS EITHER SIDE OF ETA(LEVEL_MAXWIND)
C   (THIS MAY OR MAY NOT BE A TURNING POINT OF THE CURVE)
C-----------------------------------------------------------------------
      POINTS=20
      INVPTS=1./POINTS
      DO NO=FIRST_FLD_PT,LAST_U_FLD_PT
        F(NO)=A(NO,MATRIX_P_O)
      ENDDO
      DO 208 I=MATRIX_P_O-1,1,-1
        DO NO =FIRST_FLD_PT,LAST_U_FLD_PT
          J=LEVEL_MAXWIND(NO)
          IF (J.GE.2) THEN
            X=ETA_HALF(J+1)
C-----------------------------------------------------------------------
C   CALCULATE THE WIND AT ETA(J+.5) USING A POLYNOMIAL
C   F=A(1)+A(2)*X+...+A(N)*X**(N-1)
C-----------------------------------------------------------------------
            F(NO)=F(NO)*X+A(NO,I)
          ENDIF
        ENDDO
 208  CONTINUE
      DO NO=FIRST_FLD_PT,LAST_U_FLD_PT
        J=LEVEL_MAXWIND(NO)
        IF (J.GE.2) THEN
          X=ETA_HALF(J+1)
          IF (SPEED_MAXWIND(NO).LE.F(NO)) THEN
            SPEED_MAXWIND(NO)=F(NO)
            ETA_MAXWIND(NO)=X
          ENDIF
        ENDIF
      ENDDO
C-----------------------------------------------------------------------
C   INCREMENT FROM ETA(J+.5) TO ETA(J-.5).
C   CALCULATE THE WIND AT EACH POINT AND SET IT AS THE SPEED OF THE
C   MAXIMUM WIND IF IT IS A LOCAL MAXIMUM
C-----------------------------------------------------------------------
      DO K=1,POINTS
        DO NO=FIRST_FLD_PT,LAST_U_FLD_PT
          F(NO)=A(NO,MATRIX_P_O)
        ENDDO
        DO I=MATRIX_P_O-1,1,-1
          DO NO=FIRST_FLD_PT,LAST_U_FLD_PT
            J=LEVEL_MAXWIND(NO)
            IF (J.GE.2) THEN
              XINC=(ETA_HALF(J)-ETA_HALF(J+1))*INVPTS
              X=ETA_HALF(J+1)+(XINC*K)
              F(NO)=F(NO)*X+A(NO,I)
            ENDIF
          ENDDO
        ENDDO
        DO NO=FIRST_FLD_PT,LAST_U_FLD_PT
          J=LEVEL_MAXWIND(NO)
          IF (J.GE.2) THEN
            XINC=(ETA_HALF(J)-ETA_HALF(J+1))*INVPTS
            X=ETA_HALF(J+1)+(XINC*K)
            IF (F(NO).GT.SPEED_MAXWIND(NO)) THEN
              SPEED_MAXWIND(NO)=F(NO)
              ETA_MAXWIND(NO)=X
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C-----------------------------------------------------------------------
CL  7. FIND THE DIRECTION OF THE MAXIMUM WIND BY LINEAR
CL     INTERPOLATION IN LOG(ETA) FROM SURROUNDING LEVELS
C-----------------------------------------------------------------------
C   MAKE SURE THAT THE MAXIMUM WIND IS IN THE LAYER ABOVE
C   LEVEL_MAXWIND. NOTE HERE J=LEVEL_MAXWIND(NO)
C-----------------------------------------------------------------------
      DO NO=FIRST_FLD_PT,LAST_U_FLD_PT
        J=LEVEL_MAXWIND(NO)
        IF (J.GE.2) THEN
          IF (ETA_LEVEL(J).LT.ETA_MAXWIND(NO)) J=J-1
          DIFF=DIRECTION(NO,J+1)-DIRECTION(NO,J)
C-----------------------------------------------------------------------
C PUT DIFF = ( D(J+1)-D(J)         IF ^D(J+1)-D(J)^.LT.PI
C            ( TWOPI-^D(J+1)-D(J)^ IF ^D(J+1)-D(J)^.GT.PI
C            (                               AND D(J+1).LT.D(J)
C            ( D(J+1)-D(J)-TWOPI   IF ^D(J+1)-D(J)^.GT.PI
C            (                               AND D(J+1).GT.D(J)
C    WHERE D=DIRECTION
C-----------------------------------------------------------------------
          IF (ABS(DIFF).GT.PI) THEN
            IF (DIFF.LT.0) THEN
              DIFF=TWOPI+DIFF
            ELSE
              DIFF=DIFF-TWOPI
            ENDIF
          ENDIF
C-----------------------------------------------------------------------
C   DIRECTION OF MAX WIND
C     =DIRECTION+DIFF*(LOG(ETA_MAXWIND/ETA(J))/LOG(ETA(J+1)/ETA(J)))
C-----------------------------------------------------------------------
          ETA_J=ETA_LEVEL(J)
          ETA_JP1=ETA_LEVEL(J+1)
          MULT=LOG(ETA_MAXWIND(NO)/ETA_J)/LOG(ETA_JP1/ETA_J)
          DIRECTION_MAXWIND(NO)=DIRECTION(NO,J)+DIFF*MULT
C-----------------------------------------------------------------------
C   OBTAIN A VALUE FOR DIRECTION_MAXWIND LYING BETWEEN 0 AND TWOPI
C-----------------------------------------------------------------------
          DIRECTION_MAXWIND(NO)=DIRECTION_MAXWIND(NO)+PI
          S1=TWOPI.LT.DIRECTION_MAXWIND(NO)
          S2=DIRECTION_MAXWIND(NO).LT.0.0
          IF (S1) DIRECTION_MAXWIND(NO)=DIRECTION_MAXWIND(NO)-TWOPI
          IF (S2) DIRECTION_MAXWIND(NO)=DIRECTION_MAXWIND(NO)+TWOPI
          IF (TWOPI.LT.DIRECTION_MAXWIND(NO).OR.
     &      DIRECTION_MAXWIND(NO).LT.0.0) DIRECTION_MAXWIND(NO)=0.0
        ENDIF
C-----------------------------------------------------------------------
CL  8. CONVERT THE SPEED OF MAXIMUM WIND AND DIRECTION TO U AND V
CL     COMPONENTS AND CONVERT THE ETA LEVEL TO PRESSURE LEVEL
C-----------------------------------------------------------------------
C  NOTE HERE J=LEVEL_MAXWIND(NO)
C-----------------------------------------------------------------------
        U_MAXWIND(NO)=SPEED_MAXWIND(NO)*SIN(DIRECTION_MAXWIND(NO)+PI)
        V_MAXWIND(NO)=SPEED_MAXWIND(NO)*COS(DIRECTION_MAXWIND(NO)+PI)
        ETA_JPH=ETA_HALF(J+1)
        ETA_JMH=ETA_HALF(J)
        P_JPH=AKH(J+1)+BKH(J+1)*PSTAR_BAR(NO)
        P_JMH=AKH(J)+BKH(J)*PSTAR_BAR(NO)
        MULT=(P_JMH-P_JPH)/(ETA_JMH-ETA_JPH)
        PRESSURE_MAXWIND(NO)=P_JPH+(ETA_MAXWIND(NO)-ETA_JPH)*MULT
      ENDDO
C=======================================================================
C END OF WINDMAX
C=======================================================================
      RETURN
      END
C=======================================================================
