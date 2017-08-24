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
CLL
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL   VERSION FOR CRAY Y-MP
CLL
CLL   WRITTEN  M.H MAWSON.
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.4    06/08/94 New release 1B with error in second order term
CLL                   corrected in addition to faster multi-tasked
CLL                   code achieved by inserting micro tasking
CLL                   directives and code restructuring
CLL                   to improve parallel efficiency on C90.
CLL                   X_FIELD passed as argument to reduce memory
CLL                   usage when 2nd order advection used.
CLL                   Authors: A. Dickinson, D. Salmond
CLL                   Reviewer: M. Mawson
CLL   3.4   28/10/94  Argument LLINTS added and passed to V_CORIOL
CLL                   Argument LWHITBROM added and passed to ADV_U_GD
CLL                            R.T.H.Barnes pp. S.J.Swarbrick
!     3.5    28/03/95 MPP code: Change updateable area and
!                     add boundary swaps.  P.Burton
CLL
CLL   4.0    14/02/95 Option to run with half_timestep at top level
CLL                   removed. Author: T.Davies,  Reviewer: M. Mawson
!     4.1    29/04/96 Remove MPP code (new QTADV1C version for MPP)
!                     and add TYPFLDPT arguments       P.Burton
!LL 4.3      24/04/97 Fix to 4th order calculations -
!LL                   Calculation of NUY via ISMIN   P.Burton
!LL  4.5  05/05/98  Recode -DEF,CRAY loops to find minimum of NUX/NUY
!LL                 to vectorize on Fujitsu VPP700. Also improve
!LL                 efficiency in section 3.3.      RBarnes@ecmwf.int
!LL
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
     &              COS_U_LONGITUDE,SIN_U_LONGITUDE,SEC_P_LATITUDE,
     &              AKH,BKH,OMEGA,L_SECOND,LLINTS,
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
     &, P_LAST_ROW       ! Last updatable row on pressure point field
     &, U_LAST_ROW       ! Last updatable row on wind point field
     &, P_BOT_ROW_START  ! First point of south-pole (global) or
!                        ! Southern (LAM) row on press-point field
     &, U_BOT_ROW_START  ! First point of south-pole (global) or
!                        ! Southern (LAM) row on wind-point field
     &, upd_P_ROWS       ! number of P_ROWS to be updated
     &, upd_U_ROWS       ! number of U_ROWS to be updated
     &, FIRST_FLD_PT     ! First point on field
     &, LAST_P_FLD_PT    ! Last point on pressure point field
     &, LAST_U_FLD_PT    ! Last point on wind point field
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
     &, tot_P_ROWS         ! total number of P_ROWS on grid
     &, tot_U_ROWS         ! total number of U_ROWS on grid
     &, GLOBAL_ROW_LENGTH  ! length of a global row
     &, GLOBAL_P_FIELD     ! size of a global P field
     &, GLOBAL_U_FIELD     ! size of a global U field
!


! End of comdeck TYPFLDPT

C LOGICAL VARIABLE
      LOGICAL
     &  L_SECOND     ! SET TO TRUE IF NU_BASIC IS ZERO.
     & ,LLINTS              ! Switch for linear TS calc in CALC_TS
     & ,LWHITBROM           ! Switch for White & Bromley terms

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
C*---------------------------------------------------------------------

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

      REAL
     & DELTA_AKH(P_LEVELS+1)     ! LAYER THICKNESS  AK(K) - AK(K-1)
     &,DELTA_BKH(P_LEVELS+1)     ! LAYER THICKNESS  BK(K) - BK(K-1)
     &,WK(U_FIELD)               ! WK AS IN EQUATION (46).

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
     &  I,I1,J,KP,KM,IK,K

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
      DO K=1,P_LEVELS
CL    INTERPOLATE RS ONTO U GRID.
          CALL P_TO_UV(RS(1,K),RS_U(1,K),P_FIELD,U_FIELD,ROW_LENGTH,
     &                 tot_P_ROWS)
      ENDDO

CL    INTERPOLATE ETADOT ONTO U GRID AND INCLUDE BOTTOM AND TOP
CL    BOUNDARY CONDITION

      DO K =2, P_LEVELS
        CALL P_TO_UV(ETADOT_MEAN(1,K),ETADOT_U(1,K),P_FIELD,U_FIELD,
     &                 ROW_LENGTH,tot_P_ROWS)
      END DO
! Loop over field
      DO I = FIRST_VALID_PT,LAST_U_VALID_PT
        ETADOT_U(I,1) = 0.
        ETADOT_U(I,P_LEVELS+1) = 0.
      END DO

      IF (LWHITBROM) THEN
CL    CALCULATE BRSP TERM AT LEVEL K
C STORE IN OMEGA TO SAVE WORKSPACE

      K=1
! Loop over field
      DO I=FIRST_VALID_PT,LAST_U_VALID_PT
        OMEGA(I,K)=(3.*RS_U(I,K)+RS_U(I,K+1))*(RS_U(I,K)-RS_U(I,K+1))
     &                *BKH(K+1)*.25*(PSTAR(I)-PSTAR_OLD(I))
      ENDDO
      K=P_LEVELS
! Loop over field
      DO I=FIRST_VALID_PT,LAST_U_VALID_PT
        OMEGA(I,K)=-(3.*RS_U(I,K)+RS_U(I,K-1))*(RS_U(I,K)-RS_U(I,K-1))
     &                *BKH(K)*.25*(PSTAR(I)-PSTAR_OLD(I))
      ENDDO

      DO K=2,P_LEVELS -1
! Loop over field
        DO I=FIRST_VALID_PT,LAST_U_VALID_PT
          OMEGA(I,K)=((3.*RS_U(I,K)+RS_U(I,K+1))
     &              *(RS_U(I,K)-RS_U(I,K+1))*BKH(K+1)
     &              *.25*(PSTAR(I)-PSTAR_OLD(I)))
     &              -((3.*RS_U(I,K)+RS_U(I,K-1))
     &              *(RS_U(I,K)-RS_U(I,K-1))*BKH(K)
     &              *.25*(PSTAR(I)-PSTAR_OLD(I)))
        ENDDO

      ENDDO
      END IF

CFPP$ NOCONCUR
      DO I=2,P_LEVELS
        DELTA_AKH(I) = AK(I) - AK(I-1)
        DELTA_BKH(I) = BK(I) - BK(I-1)
      END DO
C THESE ZERO VALUES SAVE HAVING TO PASS THE ZERO VERTICAL VELOCITIES
C ON LOWER AND UPPER BOUNDARIES TO V_CORIOL AS THE ZERO VELOCITIES ARE
C NOT HELD. (SEE CALL TO V_CORIOL IN SECTION 3.3)
      DELTA_AKH(1) = 0.
      DELTA_BKH(1) = 0.
      DELTA_AKH(P_LEVELS+1) = 0.
      DELTA_BKH(P_LEVELS+1) = 0.

CL---------------------------------------------------------------------
CL    SECTION 2.     ADVECTION OF U AND V.
CL                   SECTION 2 WILL CALCULATE PROVISIONAL VALUES OF
CL                   U AND V. SECTION 3 WILL CALCULATE FINAL VALUES.
CL---------------------------------------------------------------------

CL LOOP OVER P_LEVELS.
cmic$ parallel shared (advection_timestep, akh, bkh)
cmic$*     shared(cos_u_longitude,sin_u_longitude)
cmic$*     shared(longitude_step_inverse,latitude_step_inverse)
cmic$*     shared(f1,f2,omega,delta_akh,delta_bkh,ak,bk)
cmic$*     shared (delta_ak, delta_bk)
cmic$*     shared (etadot_u, l_second, lwhitbrom, llints)
cmic$*     shared (nu_basic, nux, nuy)
cmic$*     shared (p_field, pstar)
cmic$*     shared (pstar_old, p_levels)
cmic$*     shared (row_length,rs_u, sec_p_latitude, sec_u_latitude)
CMIC@a   SHARED(FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
CMIC@b   P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
CMIC@c   FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
CMIC@d   FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
CMIC@e   VALID_P_ROWS, VALID_U_ROWS,
CMIC@f   START_POINT_NO_HALO, START_POINT_INC_HALO,
CMIC@g   END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
CMIC@h   END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
CMIC@i   FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
CMIC@j   GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD)
cmic$*     shared (u,v, u_field, u_mean)
cmic$*     shared (v_mean)
cmic$*     private (u_first_inc,v_first_inc)
cmic$*     shared (u_prov,v_prov)
cmic$*     shared (u_mean_p,v_mean_p)
cmic$*     private (const1, i, i1,wk)
cmic$*     private (ik, j, k, km, kp, kappa_dum )
cmic$*     private (omega_p, p_exl_dum, p_exner_full, p_exu_dum, pk)
cmic$*     private (pk1, pl_dum, pu_dum, scalar1, scalar2)
cmic$*     private (scalar3,scalar4)
cmic$*     private (u_second_inc, v_second_inc, timestep)
cmic$ do parallel

      DO K=1,P_LEVELS

        TIMESTEP = ADVECTION_TIMESTEP

CL---------------------------------------------------------------------
CL    SECTION 2.0    INTERPOLATE U_MEAN AND V_MEAN TO P GRID.
CL                   INTERPOLATE RS AND ETADOT TO U GRID.
CL---------------------------------------------------------------------

CL    INTERPOLATE U_MEAN ONTO P GRID.

        CALL UV_TO_P(U_MEAN(1,K),U_MEAN_P(1,K),U_FIELD,U_FIELD,
     &               ROW_LENGTH,upd_U_ROWS+2)

CL    INTERPOLATE V_MEAN ONTO P GRID.

        CALL UV_TO_P(V_MEAN(1,K),V_MEAN_P(1,K),U_FIELD,U_FIELD,
     &               ROW_LENGTH,upd_U_ROWS+2)

C ---------------------------------------------------------------------
CL    SECTION 2.1    SET NU DEPENDENT ON NU_BASIC AND MAX COURANT
CL                   NUMBER.
C ---------------------------------------------------------------------
CL IF NU_BASIC NOT EQUAL TO ZERO.
          IF(.NOT.L_SECOND) THEN
CL    THEN SET NU DEPENDENT ON NU_BASIC AND MAX
CL    COURANT NUMBER.
CL CALCULATE COURANT NUMBER SQUARED.
! Loop over field missing top and bottom rows
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
          END DO
C     SET NUX EQUAL TO MINIMUM ALONG EACH ROW
           DO J=1,upd_U_ROWS
          I1 = START_POINT_NO_HALO + (J-1)*ROW_LENGTH
          SCALAR1 = NUX(I1,K)
          DO I=I1+1,I1+ROW_LENGTH-1
            IF(NUX(I,K).LT.SCALAR1) THEN
              SCALAR1 = NUX(I,K)
            END IF
          END DO
          IF(SCALAR1.LT.0.) SCALAR1 = 0.
          DO I=I1,I1+ROW_LENGTH-1
            NUX(I,K) = SCALAR1
          END DO
          END DO

C     SET NUY EQUAL TO MINIMUM ALONG EACH COLUMN
          DO J=1,ROW_LENGTH
          I1 = START_POINT_NO_HALO+ J-1
          SCALAR1 = NUY(I1,K)
          DO I=I1+ROW_LENGTH,END_U_POINT_NO_HALO,ROW_LENGTH
            IF(NUY(I,K).LT.SCALAR1) THEN
              SCALAR1 = NUY(I,K)
            END IF
          END DO
          IF(SCALAR1.LT.0.) SCALAR1 = 0.
            DO I=I1,END_U_POINT_NO_HALO,ROW_LENGTH
            NUY(I,K) = SCALAR1
            END DO
          END DO
          END IF

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
     &                    U_MEAN_P(1,K),V_MEAN_P(1,K),
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
     &                    TIMESTEP,LATITUDE_STEP_INVERSE,
     &                    LONGITUDE_STEP_INVERSE,SEC_P_LATITUDE,
     &                    OMEGA(1,K),L_SECOND,LWHITBROM)

CL    CALL ADV_U_GD FOR V.
          CALL ADV_U_GD(V(1,KM),V(1,K),V(1,KP),
     &                    U_MEAN_P(1,K),V_MEAN_P(1,K),
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
     &                    TIMESTEP,LATITUDE_STEP_INVERSE,
     &                    LONGITUDE_STEP_INVERSE,SEC_P_LATITUDE,
     &                    OMEGA(1,K),L_SECOND,LWHITBROM)

C ---------------------------------------------------------------------
CL    SECTION 2.4    REMOVE MASS-WEIGHTING FROM INCREMENT AND ADD ONTO
CL                   FIELD TO OBTAIN INTERMEDIATE VALUE.
C ---------------------------------------------------------------------

          DO I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
            SCALAR1 = 1./(RS_U(I,K)*RS_U(I,K)
     &                    *(DELTA_AK(K)+DELTA_BK(K)*PSTAR_OLD(I)))
            U_PROV(I,K) = U(I,K)- U_FIRST_INC(I)*SCALAR1
            V_PROV(I,K) = V(I,K)-V_FIRST_INC(I)*SCALAR1
          END DO

CL    LIMITED AREA MODEL THEN FORM PROVISIONAL VALUES ON BOUNDARIES
CL    EQUAL TO FIELD VALUES AT OLD TIME LEVEL.
          DO I=1,ROW_LENGTH
            IK = U_FIELD - ROW_LENGTH + I
            U_PROV(I,K)= U(I,K)
            V_PROV(I,K)= V(I,K)
            U_PROV(IK,K)= U(IK,K)
            V_PROV(IK,K)= V(IK,K)
          END DO

      enddo
cmic$ do parallel
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
     &                  U_MEAN_P(1,K),V_MEAN_P(1,K),ETADOT_U(1,K),
     &                  ETADOT_U(1,K+1),SEC_U_LATITUDE,
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
     &                  TIMESTEP,LATITUDE_STEP_INVERSE,
     &                  LONGITUDE_STEP_INVERSE,SEC_P_LATITUDE,
     &                  OMEGA(1,K),
     &                  L_SECOND,LWHITBROM)

CL    CALL ADV_U_GD FOR V.
          CALL ADV_U_GD(V_PROV(1,KM),V_PROV(1,K),V_PROV(1,KP),
     &                  U_MEAN_P(1,K),V_MEAN_P(1,K),ETADOT_U(1,K),
     &                  ETADOT_U(1,K+1),SEC_U_LATITUDE,
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
     &                  TIMESTEP,LATITUDE_STEP_INVERSE,
     &                  LONGITUDE_STEP_INVERSE,SEC_P_LATITUDE,
     &                  OMEGA(1,K),
     &                  L_SECOND,LWHITBROM)

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
     &              LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     &              WK,U_FIELD,OMEGA(1,K),LLINTS)

C ---------------------------------------------------------------------
CL    SECTION 3.3    CALCULATE TOTAL MASS-WEIGHTED INCREMENT TO FIELD
CL                   INCLUDING CORIOLIS TERM AND ADD ONTO MASS-WEIGHTED
CL                   FIELD.
CL                   IF GLOBAL CALL POLAR_UV TO UPDATE POLAR VALUES.
CL                   IF LIMITED AREA MASS-WEIGHT BOUNDARY VALUES.
C ---------------------------------------------------------------------

! Loop over field, missing top and bottom rows
      DO I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
        SCALAR1=1.0/(RS_U(I,K)*RS_U(I,K)*
     &                        (DELTA_AK(K)+DELTA_BK(K)*PSTAR(I)))
        U_SECOND_INC(I)=U_SECOND_INC(I)*SCALAR1
        V_SECOND_INC(I)=V_SECOND_INC(I)*SCALAR1
        WK(I)=WK(I)*SCALAR1
      END DO
CL    TOTAL MASS-WEIGHTED INCREMENT IS CALCULATED INCLUDING VERTICAL
CL    CORIOIS TERM AND ADDED ONTO MASS-WEIGHTED FIELD.

! Loop over field, missing top and bottom rows
          DO I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
          U(I,K)=0.5 * (U(I,K)-U_SECOND_INC(I)+U_PROV(I,K))

          V(I,K)=0.5 * (V(I,K)-V_SECOND_INC(I)+V_PROV(I,K))
          END DO
          IF (LWHITBROM) THEN
            DO I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
              SCALAR3 = 1.0/RS_U(I,K)
              U(I,K) = U(I,K) -(F2(I) + U(I,K)*SCALAR3)*WK(I)*TIMESTEP
              V(I,K) = V(I,K) +(F1(I) - V(I,K)*SCALAR3)*WK(I)*TIMESTEP
            ENDDO
          ENDIF
CL    SET POLAR VALUES FOR OMEGA

        DO I=1,ROW_LENGTH
          OMEGA(I,K)=OMEGA(I+ROW_LENGTH,K)
          OMEGA(U_FIELD-ROW_LENGTH+I,K)=OMEGA(U_FIELD-2*ROW_LENGTH
     &      +I,K)
        END DO



CL END LOOP OVER P_LEVELS
      enddo
cmic$ end parallel

CL MASS WEIGHT THE OUTPUT FIELDS
       DO K=1,P_LEVELS
         DO I=FIRST_FLD_PT,LAST_U_FLD_PT
           U(I,K)=U(I,K)*RS_U(I,K)*RS_U(I,K)*(DELTA_AK(K)+
     &               DELTA_BK(K)*PSTAR(I))
           V(I,K)=V(I,K)*RS_U(I,K)*RS_U(I,K)*(DELTA_AK(K)+
     &               DELTA_BK(K)*PSTAR(I))
         END DO
       END DO

CL    END OF ROUTINE UV_ADV

      RETURN
      END
