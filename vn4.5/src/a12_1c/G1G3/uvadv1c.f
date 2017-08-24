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
CLL   SUBROUTINE UV_ADV -------------------------------------------
CLL
CLL                   PURPOSE:
CLL  CALCULATES MASS-WEIGHTED INCREMENTS TO U AND V DUE TO
CLL  ADVECTION  BY USING EQUATIONS (37) AND (38) TO CALCULATE
CLL  PROVISIONAL VALUES OF U AND V AT THE NEW TIME-LEVEL, AND THEN
CLL  RECALCULATING THE ADVECTION TERMS ON THE RIGHT-HAND SIDE OF (41)
CLL  AND (42) USING THESE PROVISIONAL VALUES.  THE CORIOLIS TERMS
CLL  ASSOCIATED WITH THE VERTICAL VELOCITY ARE CALCULATED AND INCLUDED
CLL  IN THE INCREMENTS.  THE FINAL INCREMENTS ARE CALCULATED AS IN
CLL  EQUATIONS (41) AND (42). IF RUNNING A GLOBAL MODEL POLAR_UV IS
CLL  CALLED TO UPDATE POLAR VALUES.
CLL
CLL                          CHANGES INCLUDE:-
CLL U_MEAN AND V_MEAN FIELDS NOT OVER-WRITTEN WHEN INTERPOLATION TO
CLL U_GRID PERFORMED. ETADOT AND RS FIELDS INTERPOLATED TO U_GRID INSIDE
CLL THIS ROUTINE INSTEAD OF INSIDE ADV_CTL. THIS COSTS 8 EXTRA
CLL HORIZONTAL FIELDS BUT ALLOWS ROUTINE TO BE CALLED BEFORE TH_ADV SO
CLL THAT OMEGA CALCULATED HERE CAN BE USED INSIDE TH_ADV TO CALCULATE
CLL EXTRA THERMODYNAMIC TERM.
CLL
CLL INCLUSION OF L_SECOND TO CHOOSE CHEAPER SECOND ORDER ADVECTION
CLL SCHEME ALONG WITH REMOVAL OF CODE PREVIOUSLY UNDER *DEF FORECAST.
CLL CODE INCLUDED TO ALLOW HALF-TIMESTEP TO BE USED AT TOP LEVEL.
CLL
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL   VERSION FOR CRAY Y-MP
CLL
CLL   WRITTEN  M.H MAWSON.
CLL   MPP CODE ADDED BY P.BURTON
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 4.1:
CLL VERSION  DATE
CLL 4.1      08/12/95 New version of routine specifically for MPP
CLL                   P.Burton
!LL 4.2      10/01/97 Initialise unprocessed points in U_PROV
!LL                   and V_PROV. D. Robinson.
!LL 4.3      24/04/97 Fixes to 4th order calculations   P.Burton
C     vn4.3    Mar. 97   T3E migration : optimisation changes
C                                       D.Salmond

CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD B.
CLL
CLL   SYSTEM COMPONENTS COVERED: P122
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION:       THE EQUATIONS USED ARE (37-38) AND (41-42)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER NO. 10
CLL                        M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLL
CLLEND-------------------------------------------------------------

C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE UV_ADV
     &              (U,V,PSTAR_OLD,PSTAR,U_MEAN,V_MEAN,SEC_U_LATITUDE,
     &              ETADOT_MEAN,RS,DELTA_AK,DELTA_BK,AK,BK,F1,F2,
     &              LATITUDE_STEP_INVERSE,ADVECTION_TIMESTEP,NU_BASIC,
     &              LONGITUDE_STEP_INVERSE,U_FIELD,P_FIELD,
     &              ROW_LENGTH,P_LEVELS,
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
     &              COS_U_LONGITUDE,SIN_U_LONGITUDE,SEC_P_LATITUDE,
     &              AKH,BKH,OMEGA,L_SECOND,LLINTS,
     &              extended_address,
     &              LWHITBROM,X_FIELD)

      IMPLICIT NONE

      INTEGER
     &  P_FIELD            !IN DIMENSION OF FIELDS ON PRESSSURE GRID.
     &, U_FIELD            !IN DIMENSION OF FIELDS ON VELOCITY GRID
     &, X_FIELD            !IN 1 IF 2ND ORDER ELSE U_FIELD
     &, P_LEVELS           !IN NUMBER OF PRESSURE LEVELS.
     &, ROW_LENGTH         !IN NUMBER OF POINTS PER ROW

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

C LOGICAL VARIABLE
      LOGICAL
     &  L_SECOND     ! SET TO TRUE IF NU_BASIC IS ZERO.
     & ,LLINTS              ! Switch for linear TS calc in CALC_TS
     & ,LWHITBROM           ! Switch for White & Bromley terms
      INTEGER extended_address(P_FIELD)

      REAL
     & U_MEAN(U_FIELD,P_LEVELS) !IN AVERAGED MASS-WEIGHTED U VELOCITY
     &                          !   FROM ADJUSTMENT STEP HELD AT U
     &                          !   POINTS.
     &,V_MEAN(U_FIELD,P_LEVELS) !IN AVERAGED MASS-WEIGHTED V VELOCITY
     &                          !   * COS(LAT) FROM ADJUSTMENT STEP
     &,ETADOT_MEAN(P_FIELD,P_LEVELS)  !IN AVERAGED MASS-WEIGHTED
     &                          !VERTICAL VELOCITY FROM ADJUSTMENT STEP

      REAL
     & U(U_FIELD,P_LEVELS)      !INOUT IN U FIELD,
     &                          !  OUT MASS-WEIGHTED U FIELD.
     &,V(U_FIELD,P_LEVELS)      !INOUT IN V FIELD,
     &                          !  OUT MASS-WEIGHTED V FIELD.

      REAL
     & PSTAR(U_FIELD)           !IN PSTAR FIELD AT NEW TIME-LEVEL ON
     &                          ! U GRID.
     &,PSTAR_OLD(U_FIELD)       !IN PSTAR AT PREVIOUS TIME-LEVEL ON
     &                          ! U GRID.
     &,RS(P_FIELD,P_LEVELS)     !IN RS FIELD.
     &,AK(P_LEVELS)             !IN FIRST TERM IN HYBRID CO-ORDS.
     &,BK(P_LEVELS)             !IN SECOND TERM IN HYBRID CO-ORDS.
     &,DELTA_AK(P_LEVELS)       !IN LAYER THICKNESS
     &,DELTA_BK(P_LEVELS)       !IN LAYER THICKNESS
     &,AKH(P_LEVELS+1)          !IN HYBRID CO-ORDINATE AT HALF LEVELS
     &,BKH(P_LEVELS+1)          !IN HYBRID CO-ORDINATE AT HALF LEVELS
     &,SEC_U_LATITUDE(U_FIELD)  !IN 1/COS(LAT) AT U POINTS (2-D ARRAY)
     &,SEC_P_LATITUDE(U_FIELD)  !IN 1/COS(LAT) AT P POINTS (2-D ARRAY)
     &,SIN_U_LONGITUDE(ROW_LENGTH)  !IN SIN(LONGITUDE) AT U POINTS.
     &,COS_U_LONGITUDE(ROW_LENGTH)  !IN COS(LONGITUDE) AT U POINTS.

      REAL
     & LONGITUDE_STEP_INVERSE   !IN 1/(DELTA LAMDA)
     &,LATITUDE_STEP_INVERSE    !IN 1/(DELTA PHI)
     &,ADVECTION_TIMESTEP       !IN
     &,NU_BASIC                 !IN STANDARD NU TERM FOR MODEL RUN.
     &,F1(U_FIELD)              !IN A CORIOLIS TERM (SEE DOCUMENTATION)
     &,F2(U_FIELD)              !IN A CORIOLIS TERM (SEE DOCUMENTATION)

      REAL
     & OMEGA(U_FIELD,P_LEVELS) !OUT TRUE VERTICAL VELOCITY

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: 35 ARE REQUIRED
      REAL
     & RS_U(U_FIELD,P_LEVELS)          ! RS AT U POINTS FOR CURRENT LEVE
     &,ETADOT_U(U_FIELD,P_LEVELS+1)  ! ETADOT AT U POINTS FOR CURRENT LE
     &,U_MEAN_P(U_FIELD,P_LEVELS) ! U MEAN AT P POINTS FOR CURRENT LEVEL
     &                   !   WITH FIRST POINT OF FIELD NOW
     &                   !   BEING FIRST P POINT ON SECOND ROW
     &                   !   OF P-GRID.
     &,V_MEAN_P(U_FIELD,P_LEVELS) ! V MEAN AT P POINTS FOR CURRENT LEVEL
     &                   !   WITH FIRST POINT OF FIELD NOW
     &                   !   BEING FIRST P POINT ON SECOND ROW
     &                   !   OF P-GRID.

      REAL
     & U_FIRST_INC(U_FIELD)       ! HOLDS U INCREMENT
     &                            !RETURNED BY FIRST CALL TO ADV_U_GD
     &,U_SECOND_INC(U_FIELD)      ! HOLDS U INCREMENT
     &                            !RETURNED BY SECOND CALL TO ADV_U_GD
     &,U_PROV(U_FIELD,P_LEVELS)            ! HOLDS PROVISIONAL VALUE OF

      REAL
     & V_FIRST_INC(U_FIELD)       ! HOLDS V INCREMENT
     &                            !RETURNED BY FIRST CALL TO ADV_U_GD
     &,V_SECOND_INC(U_FIELD)      ! HOLDS V INCREMENT
     &                            !RETURNED BY SECOND CALL TO ADV_U_GD
     &,V_PROV(U_FIELD,P_LEVELS)            ! HOLDS PROVISIONAL VALUE OF

C NP DENOTES NORTH POLE, SP DENOTES SOUTH POLE.
C POLAR INCREMENT ARRAYS ARE NOT USED IN LIMITED AREA MODEL BUT TO
C REMOVE THEM WOULD LEAD TO MODIFYING THE NUMBER OF VARIABLES
C PASSED TO ADV_U_GD. THE RETENTION OF THESE ARRAYS ADDS ONLY
C 12*ROW_LENGTH TO THE SPACE USED AND NOTHING TO THE CALCULATION
C TIME AS ALL USES OF THEM IN CALCULATION ARE CONTROLLED BY *IF'S.

      REAL
     & NUX(X_FIELD,P_LEVELS)      ! COURANT NUMBER DEPENDENT NU AT U POI
     &                    ! USED IN EAST-WEST ADVECTION.
     &,NUY(X_FIELD,P_LEVELS)      ! COURANT NUMBER DEPENDENT NU AT U POI
     &                    ! USED IN NORTH-SOUTH ADVECTION.

      REAL NUX_MIN(upd_P_ROWS),  ! minimum value of NUX along a row
     &     NUY_MIN(ROW_LENGTH)  ! min of NUY along a column

      REAL
     & DELTA_AKH(P_LEVELS+1)     ! LAYER THICKNESS  AK(K) - AK(K-1)
     &,DELTA_BKH(P_LEVELS+1)     ! LAYER THICKNESS  BK(K) - BK(K-1)
     &,WK(U_FIELD)               ! WK AS IN EQUATION (46).

! Work space required to allow the use of Fourth Order Advection
! U/V_MEAN_P_COPY and U/V_COPY arrays are defined with an extra halo
! this is required for the bigger stencil of the 4th order operator.

      REAL U_MEAN_P_COPY((ROW_LENGTH+2*extra_EW_Halo)*
     &                   (tot_U_ROWS+2*extra_NS_Halo),P_LEVELS),
     &  !    Copy of U_MEAN with extra halo space for 4th order
     &     V_MEAN_P_COPY((ROW_LENGTH+2*extra_EW_Halo)*
     &                   (tot_U_ROWS+2*extra_NS_Halo),P_LEVELS),
     &  !    Copy of V_MEAN with extra halo space for 4th order
     &     U_COPY((ROW_LENGTH+2*extra_EW_Halo)*
     &            (tot_U_ROWS+2*extra_NS_Halo),P_LEVELS),
     &  !    Copy of U with extra halo space for 4th order
     &     V_COPY((ROW_LENGTH+2*extra_EW_Halo)*
     &            (tot_U_ROWS+2*extra_NS_Halo),P_LEVELS)
     &  !    Copy of V with extra halo space for 4th order

      INTEGER  extended_P_FIELD,
     &         extended_U_FIELD
!  These are the sizes of the arrays with the extra halos

C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES
      INTEGER
     &  U_POINTS_UPDATE    ! NUMBER OF U POINTS TO BE UPDATED.
     &                     !  = (ROWS-1)*ROWLENGTH

C REAL SCALARS
      REAL
     & SCALAR1,SCALAR2,SCALAR3,SCALAR4,TIMESTEP

C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     &  I,J,KP,KM,IK,K,IL
      INTEGER I_start,I_end
      INTEGER info  ! return code from comms

C*L   EXTERNAL SUBROUTINE CALLS:---------------------------------------
      EXTERNAL ADV_U_GD,POLAR_UV,V_CORIOL,UV_TO_P,P_TO_UV
C*---------------------------------------------------------------------

CL    MAXIMUM VECTOR LENGTH ASSUMED IS (ROWS+1) * ROWLENGTH
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE INCLUDING SUBROUTINE CALLS:
CL---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL    SECTION 1.     INITIALISATION
CL---------------------------------------------------------------------
C INCLUDE LOCAL CONSTANTS FROM GENERAL CONSTANTS BLOCK

      U_POINTS_UPDATE = upd_U_ROWS*ROW_LENGTH

!QAN fix for RS_U
      DO K=1,P_LEVELS
        DO I=1,U_FIELD
          RS_U(I,K)=0.0
        ENDDO
      ENDDO

      DO K=1,P_LEVELS
CL    INTERPOLATE RS ONTO U GRID.
!          CALL P_TO_UV(RS(1,K),RS_U(1,K),P_FIELD,U_FIELD,ROW_LENGTH,
!     &                 tot_P_ROWS)
        CALL P_TO_UV(RS(FIRST_VALID_PT,K),RS_U(FIRST_VALID_PT,K),
     &               P_FIELD-FIRST_VALID_PT+1,U_FIELD-FIRST_VALID_PT+1,
     &               ROW_LENGTH,VALID_P_ROWS)
      ENDDO

CL    INTERPOLATE ETADOT ONTO U GRID AND INCLUDE BOTTOM AND TOP
CL    BOUNDARY CONDITION

      DO K =2, P_LEVELS
!        CALL P_TO_UV(ETADOT_MEAN(1,K),ETADOT_U(1,K),P_FIELD,U_FIELD,
!     &                 ROW_LENGTH,tot_P_ROWS)
        CALL P_TO_UV(ETADOT_MEAN(FIRST_VALID_PT,K),
     &               ETADOT_U(FIRST_VALID_PT,K),
     &               P_FIELD-FIRST_VALID_PT+1,U_FIELD-FIRST_VALID_PT+1,
     &               ROW_LENGTH,VALID_P_ROWS)
      END DO
      DO I = FIRST_VALID_PT,LAST_U_VALID_PT
        ETADOT_U(I,1) = 0.0
        ETADOT_U(I,P_LEVELS+1) = 0.0
      END DO

      IF (LWHITBROM) THEN

!        CALL FILL_HALOS(RS_U,U_FIELD,ROW_LENGTH,P_LEVELS,fld_type_u)

CL    CALCULATE BRSP TERM AT LEVEL K
C STORE IN OMEGA TO SAVE WORKSPACE

        K=1
        DO I=FIRST_VALID_PT,LAST_U_VALID_PT
          OMEGA(I,K)=(3.*RS_U(I,K)+RS_U(I,K+1))*(RS_U(I,K)-RS_U(I,K+1))
     &                  *BKH(K+1)*.25*(PSTAR(I)-PSTAR_OLD(I))
        ENDDO
        K=P_LEVELS
        DO I=FIRST_VALID_PT,LAST_U_VALID_PT
          OMEGA(I,K)=-(3.*RS_U(I,K)+RS_U(I,K-1))*(RS_U(I,K)-RS_U(I,K-1))
     &                   *BKH(K)*.25*(PSTAR(I)-PSTAR_OLD(I))
        ENDDO

        DO K=2,P_LEVELS -1
          DO I=FIRST_VALID_PT,LAST_U_VALID_PT
            OMEGA(I,K)=((3.*RS_U(I,K)+RS_U(I,K+1))
     &                *(RS_U(I,K)-RS_U(I,K+1))*BKH(K+1)
     &                *.25*(PSTAR(I)-PSTAR_OLD(I)))
     &                -((3.*RS_U(I,K)+RS_U(I,K-1))
     &                *(RS_U(I,K)-RS_U(I,K-1))*BKH(K)
     &                *.25*(PSTAR(I)-PSTAR_OLD(I)))
          ENDDO

        ENDDO
      ENDIF

! Precalculate U_MEAN and V_MEAN interpolated onto P grid - since it
! requires a call to SWAPBOUNDS, if we do it outside the main loop
! over levels, we can do just one call rather than a seperate call
! for each level (inefficient)

      DO K=1,P_LEVELS
! QAN fix
      DO I=1,U_FIELD
        U_MEAN_P(I,K)=0.0
        V_MEAN_P(I,K)=0.0
      ENDDO
        CALL UV_TO_P(U_MEAN(FIRST_VALID_PT,K),
     &               U_MEAN_P(FIRST_VALID_PT,K),
     &               U_FIELD-FIRST_VALID_PT+1,
     &               U_FIELD-FIRST_VALID_PT+1,
     &               ROW_LENGTH,upd_U_ROWS+2)
!     &               ROW_LENGTH,upd_P_ROWS+1)
!        CALL UV_TO_P(U_MEAN(1,K),U_MEAN_P(1,K),U_FIELD,U_FIELD,
!     &               ROW_LENGTH, tot_P_ROWS)

        CALL UV_TO_P(V_MEAN(FIRST_VALID_PT,K),
     &               V_MEAN_P(FIRST_VALID_PT,K),
     &               U_FIELD-FIRST_VALID_PT+1,
     &               U_FIELD-FIRST_VALID_PT+1,
     &               ROW_LENGTH,upd_U_ROWS+2)
!     &               ROW_LENGTH,upd_P_ROWS+1)
!        CALL UV_TO_P(V_MEAN(1,K),V_MEAN_P(1,K),U_FIELD,U_FIELD,
!     &               ROW_LENGTH, tot_P_ROWS)
      ENDDO

! This seems to be what was in the code originally but I can't
! understand. There seems no need since U_MEAN and V_MEAN are swapped
! in ADJ_CTL. But I would have thought the _P versions would need
! to be swapped as they are used in the advection routine. But there
! is no sign of a swap in the original code. Some experiment is
! required methinks!

!      CALL SWAPBOUNDS(U_MEAN,ROW_LENGTH,lasize(2),Offx,Offy,P_LEVELS)
!      CALL SWAPBOUNDS(V_MEAN,ROW_LENGTH,lasize(2),Offx,Offy,P_LEVELS)

CFPP$ NOCONCUR
      DO I=2,P_LEVELS
        DELTA_AKH(I) = AK(I) - AK(I-1)
        DELTA_BKH(I) = BK(I) - BK(I-1)
      ENDDO
C THESE ZERO VALUES SAVE HAVING TO PASS THE ZERO VERTICAL VELOCITIES
C ON LOWER AND UPPER BOUNDARIES TO V_CORIOL AS THE ZERO VELOCITIES ARE
C NOT HELD. (SEE CALL TO V_CORIOL IN SECTION 3.3)
      DELTA_AKH(1) = 0.0
      DELTA_BKH(1) = 0.0
      DELTA_AKH(P_LEVELS+1) = 0.0
      DELTA_BKH(P_LEVELS+1) = 0.0

! In order to use the same call to adv_u_gd for both the second and
! fourth order advection, U/V_MEAN_P are copied into _COPY arrays.
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
        CALL COPY_FIELD(U_MEAN_P,U_MEAN_P_COPY,
     &                  U_FIELD,extended_U_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  EW_Halo,NS_Halo,
     &                  .FALSE.)
        CALL COPY_FIELD(V_MEAN_P,V_MEAN_P_COPY,
     &                  U_FIELD,extended_U_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  EW_Halo,NS_Halo,
     &                  .FALSE.)

      ELSE  ! if its fourth order:

        CALL COPY_FIELD(U_MEAN_P,U_MEAN_P_COPY,
     &                  U_FIELD,extended_U_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  halo_4th,halo_4th,
     &                  .TRUE.)
        CALL COPY_FIELD(V_MEAN_P,V_MEAN_P_COPY,
     &                  U_FIELD,extended_U_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  halo_4th,halo_4th,
     &                  .TRUE.)
        CALL COPY_FIELD(U,U_COPY,
     &                  U_FIELD,extended_U_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  halo_4th,halo_4th,
     &                  .TRUE.)
        CALL COPY_FIELD(V,V_COPY,
     &                  U_FIELD,extended_U_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  halo_4th,halo_4th,
     &                  .TRUE.)

       ENDIF ! IF (L_SECOND)

CL---------------------------------------------------------------------
CL    SECTION 2.     ADVECTION OF U AND V.
CL                   SECTION 2 WILL CALCULATE PROVISIONAL VALUES OF
CL                   U AND V. SECTION 3 WILL CALCULATE FINAL VALUES.
CL---------------------------------------------------------------------

CL LOOP OVER P_LEVELS.

      DO K=1,P_LEVELS
CL SET TIMESTEP APPROPRIATE TO LEVEL

        TIMESTEP = ADVECTION_TIMESTEP

C ---------------------------------------------------------------------
CL    SECTION 2.1    SET NU DEPENDENT ON NU_BASIC AND MAX COURANT
CL                   NUMBER.
C ---------------------------------------------------------------------
CL IF NU_BASIC NOT EQUAL TO ZERO.
          IF(.NOT.L_SECOND) THEN
CL    THEN SET NU DEPENDENT ON NU_BASIC AND MAX
CL    COURANT NUMBER.
CL CALCULATE COURANT NUMBER SQUARED.
          DO I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
            SCALAR1 = U_MEAN_P(I,K)*LONGITUDE_STEP_INVERSE
            SCALAR2 = V_MEAN_P(I,K)*LATITUDE_STEP_INVERSE
            SCALAR3 = TIMESTEP/
     &                (RS_U(I,K)*RS_U(I,K)*(DELTA_AK(K)+DELTA_BK(K)*
     &                PSTAR_OLD(I)))
            SCALAR4 = SEC_U_LATITUDE(I)*SCALAR3
            SCALAR1 = SCALAR1*SCALAR1
            SCALAR2 = SCALAR2*SCALAR2
            SCALAR3 = SCALAR3*SCALAR3
            SCALAR4 = SCALAR4*SCALAR4
CL    CALCULATE NU PARAMETER.

            NUX(I,K) = (1.- SCALAR4*SCALAR1)*NU_BASIC
            NUY(I,K) = (1.- SCALAR3*SCALAR2)*NU_BASIC
          ENDDO

! Set NUX equal to minimum value along each row


          DO J=FIRST_ROW,FIRST_ROW+upd_U_ROWS-1
            I_start=(J-1)*ROW_LENGTH+FIRST_ROW_PT ! start and end of row
            I_end=(J-1)*ROW_LENGTH+LAST_ROW_PT    ! missing out halos
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
          CALL GCG_RMIN(upd_U_ROWS,GC_ROW_GROUP,info,NUX_MIN)

! and now set all values of NUX to the minimum along the row
          DO J=FIRST_ROW,FIRST_ROW+upd_U_ROWS-1
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
            I_end=I_start+(upd_U_ROWS-1)*ROW_LENGTH
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
            I_end=I_start+(upd_U_ROWS-1)*ROW_LENGTH

            DO I=I_start,I_end,ROW_LENGTH
              NUY(I,K)=NUY_MIN(J)
            ENDDO

          ENDDO ! J : loop over columns

        ENDIF  ! IF its fourth order advection

C ---------------------------------------------------------------------
CL    SECTION 2.3    CALL ADV_U_GD TO OBTAIN FIRST INCREMENT DUE TO
CL                   ADVECTION.
C ---------------------------------------------------------------------

          KP=K+1
          KM=K-1
          IF (K .EQ. P_LEVELS) THEN
            KP = K
          END IF
          IF (K .EQ. 1) THEN
            KM = K
          END IF

C BRSP IS CURRENTLY HELD IN OMEGA


          CALL ADV_U_GD(U(1,KM),U(1,K),U(1,KP),
     &                    U_MEAN_P_COPY(1,K),V_MEAN_P_COPY(1,K),
     &                    ETADOT_U(1,K),ETADOT_U(1,K+1),
     &                    SEC_U_LATITUDE,U_FIRST_INC,
     &                    NUX(1,K),NUY(1,K),U_FIELD,
     &                    ROW_LENGTH,
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
     &                    TIMESTEP,LATITUDE_STEP_INVERSE,
     &                    LONGITUDE_STEP_INVERSE,SEC_P_LATITUDE,
     &                    OMEGA(1,K),L_SECOND,LWHITBROM,
     &                    U_COPY(1,K),extended_U_FIELD,
     &                    extended_address)                 

CL    CALL ADV_U_GD FOR V.
          CALL ADV_U_GD(V(1,KM),V(1,K),V(1,KP),
     &                    U_MEAN_P_COPY(1,K),V_MEAN_P_COPY(1,K),
     &                    ETADOT_U(1,K),ETADOT_U(1,K+1),
     &                    SEC_U_LATITUDE,V_FIRST_INC,
     &                    NUX(1,K),NUY(1,K),U_FIELD,
     &                    ROW_LENGTH,
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
     &                    TIMESTEP,LATITUDE_STEP_INVERSE,
     &                    LONGITUDE_STEP_INVERSE,SEC_P_LATITUDE,
     &                    OMEGA(1,K),L_SECOND,LWHITBROM,
     &                    V_COPY(1,K),extended_U_FIELD,
     &                    extended_address)                 

C ---------------------------------------------------------------------
CL    SECTION 2.4    REMOVE MASS-WEIGHTING FROM INCREMENT AND ADD ONTO
CL                   FIELD TO OBTAIN INTERMEDIATE VALUE.
C ---------------------------------------------------------------------

          DO I=1,START_POINT_NO_HALO-1
            U_PROV(I,K) = 0.0
            V_PROV(I,K) = 0.0
          ENDDO

          DO I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
            SCALAR1 = 1./(RS_U(I,K)*RS_U(I,K)
     &                    *(DELTA_AK(K)+DELTA_BK(K)*PSTAR_OLD(I)))
            U_PROV(I,K) = U(I,K)- U_FIRST_INC(I)*SCALAR1
            V_PROV(I,K) = V(I,K)-V_FIRST_INC(I)*SCALAR1
          ENDDO

          DO I=END_U_POINT_NO_HALO+1,U_FIELD
            U_PROV(I,K) = 0.0
            V_PROV(I,K) = 0.0
          ENDDO


      ENDDO

!    IF GLOBAL MODEL CALCULATE PROVISIONAL POLAR VALUES.
!    CALL POLAR_UV TO FORM PROVISIONAL VALUES.

      CALL POLAR_UV(U_PROV,V_PROV,ROW_LENGTH,
     &              U_FIELD,P_LEVELS,
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
     &              COS_U_LONGITUDE,SIN_U_LONGITUDE)

!      CALL SET_SIDES(U_PROV,P_FIELD,ROW_LENGTH,P_LEVELS,fld_type_u)
!      CALL SET_SIDES(V_PROV,P_FIELD,ROW_LENGTH,P_LEVELS,fld_type_u)

      IF (L_SECOND) THEN

! Swap boundaries of U_PROV and V_PROV
        CALL SWAPBOUNDS(U_PROV,ROW_LENGTH,tot_U_ROWS,
     &                  EW_Halo,NS_Halo,P_LEVELS)
        CALL SWAPBOUNDS(V_PROV,ROW_LENGTH,tot_U_ROWS,
     &                  EW_Halo,NS_Halo,P_LEVELS)
!        CALL SET_SIDES(U_PROV,P_FIELD,ROW_LENGTH,P_LEVELS,fld_type_u)
!        CALL SET_SIDES(V_PROV,P_FIELD,ROW_LENGTH,P_LEVELS,fld_type_u)

      ELSE ! fourth order advection

! Copy U/V_PROV into U/V_COPY which have double halos for fourth
! order advection, and do swap to fill these halos
        CALL COPY_FIELD(U_PROV,U_COPY,
     &                  U_FIELD,extended_U_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  halo_4th,halo_4th,
     &                  .TRUE.)

        CALL COPY_FIELD(V_PROV,V_COPY,
     &                  U_FIELD,extended_U_FIELD,
     &                  ROW_LENGTH,tot_U_ROWS,P_LEVELS,
     &                  EW_Halo,NS_Halo,
     &                  halo_4th,halo_4th,
     &                  .TRUE.)

      ENDIF

      DO K=1,P_LEVELS
CL---------------------------------------------------------------------
CL    SECTION 3.     Second advection step.
CL---------------------------------------------------------------------

          TIMESTEP = ADVECTION_TIMESTEP
C ---------------------------------------------------------------------
CL    SECTION 3.1    CALL ADV_U_GD TO OBTAIN SECOND INCREMENT DUE TO
CL                   ADVECTION.
C ---------------------------------------------------------------------

          KP=K+1
          KM=K-1
          IF (K .EQ. P_LEVELS) THEN
            KP = K
          END IF
          IF (K .EQ. 1) THEN
            KM = K
          END IF

CL    CALL ADV_U_GD FOR U.

C BRSP IS CURRENTLY HELD IN OMEGA

          CALL ADV_U_GD(U_PROV(1,KM),U_PROV(1,K),U_PROV(1,KP),
     &                  U_MEAN_P_COPY(1,K),V_MEAN_P_COPY(1,K),
     &                  ETADOT_U(1,K),ETADOT_U(1,K+1),SEC_U_LATITUDE,
     &                  U_SECOND_INC,NUX(1,K),NUY(1,K),U_FIELD,
     &                  ROW_LENGTH,
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
     &                  LONGITUDE_STEP_INVERSE,SEC_P_LATITUDE,
     &                  OMEGA(1,K),L_SECOND,LWHITBROM,
     &                  U_COPY(1,K),extended_U_FIELD,
     &                  extended_address)                   

CL    CALL ADV_U_GD FOR V.
          CALL ADV_U_GD(V_PROV(1,KM),V_PROV(1,K),V_PROV(1,KP),
     &                  U_MEAN_P_COPY(1,K),V_MEAN_P_COPY(1,K),
     &                  ETADOT_U(1,K),ETADOT_U(1,K+1),SEC_U_LATITUDE,
     &                  V_SECOND_INC,NUX(1,K),NUY(1,K),U_FIELD,
     &                  ROW_LENGTH,
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
     &                  LONGITUDE_STEP_INVERSE,SEC_P_LATITUDE,
     &                  OMEGA(1,K),L_SECOND,LWHITBROM,
     &                  V_COPY(1,K),extended_U_FIELD,
     &                  extended_address)                   

C ---------------------------------------------------------------------
CL    SECTION 3.2    CALL V_CORIOL TO OBTAIN WK AS IN EQUATION (46).
C ---------------------------------------------------------------------


          CALL V_CORIOL(ETADOT_U(1,K),ETADOT_U(1,K+1),PSTAR,
     &              PSTAR_OLD,U_MEAN_P(1,K),V_MEAN_P(1,K),RS_U(1,K),
     &              SEC_U_LATITUDE,TIMESTEP,AK(K),BK(K),
     &              DELTA_AK(K),DELTA_BK(K),DELTA_AKH(K),
     &              DELTA_BKH(K),DELTA_AKH(K+1),DELTA_BKH(K+1),
     &              ROW_LENGTH,
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
     &              LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     &              WK,U_FIELD,OMEGA(1,K),LLINTS)

C ---------------------------------------------------------------------
CL    SECTION 3.3    CALCULATE TOTAL MASS-WEIGHTED INCREMENT TO FIELD
CL                   INCLUDING CORIOLIS TERM AND ADD ONTO MASS-WEIGHTED
CL                   FIELD.
CL                   IF GLOBAL CALL POLAR_UV TO UPDATE POLAR VALUES.
CL                   IF LIMITED AREA MASS-WEIGHT BOUNDARY VALUES.
C ---------------------------------------------------------------------

      DO I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
        SCALAR1=RS_U(I,K)*RS_U(I,K)*
     &                           (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I))
        U_SECOND_INC(I)=U_SECOND_INC(I)/SCALAR1
        V_SECOND_INC(I)=V_SECOND_INC(I)/SCALAR1
        WK(I)=WK(I)/SCALAR1
      END DO
CL    TOTAL MASS-WEIGHTED INCREMENT IS CALCULATED INCLUDING VERTICAL
CL    CORIOIS TERM AND ADDED ONTO MASS-WEIGHTED FIELD.

        DO I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
          SCALAR3 = 1./RS_U(I,K)
          U(I,K)=0.5 * (U(I,K)-U_SECOND_INC(I)+U_PROV(I,K))
          V(I,K)=0.5 * (V(I,K)-V_SECOND_INC(I)+V_PROV(I,K))
        ENDDO
        IF (LWHITBROM) THEN
          DO I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
            SCALAR3 = 1.0/RS_U(I,K)
            U(I,K) = U(I,K) -(F2(I) + U(I,K)*SCALAR3)*WK(I)*TIMESTEP
            V(I,K) = V(I,K) +(F1(I) - V(I,K)*SCALAR3)*WK(I)*TIMESTEP
          ENDDO
        ENDIF

CL    SET POLAR VALUES FOR OMEGA
      IF (at_top_of_LPG) THEN
        DO I=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
          OMEGA(I,K)=OMEGA(I+ROW_LENGTH,K)
        ENDDO
      ENDIF
      IF (at_base_of_LPG) THEN
        DO I=U_BOT_ROW_START,U_BOT_ROW_START+ROW_LENGTH-1
          OMEGA(I,K)=OMEGA(I-ROW_LENGTH,K)
        ENDDO
      ENDIF

CL END LOOP OVER P_LEVELS
      ENDDO

!    UPDATE POLAR VALUES BY CALLING POLAR_UV.

      CALL POLAR_UV(U,V,ROW_LENGTH,
     &              U_FIELD,P_LEVELS,
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
     &              COS_U_LONGITUDE,SIN_U_LONGITUDE)

CL MASS WEIGHT THE OUTPUT FIELDS

      DO K=1,P_LEVELS
        DO I=FIRST_FLD_PT,LAST_U_FLD_PT
          U(I,K)=U(I,K)*RS_U(I,K)*RS_U(I,K)*(DELTA_AK(K)+
     &           DELTA_BK(K)*PSTAR(I))
          V(I,K)=V(I,K)*RS_U(I,K)*RS_U(I,K)*(DELTA_AK(K)+
     &           DELTA_BK(K)*PSTAR(I))
        ENDDO
      ENDDO

CL    END OF ROUTINE UV_ADV

      RETURN
      END
