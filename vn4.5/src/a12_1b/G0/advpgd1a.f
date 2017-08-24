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
CLL   SUBROUTINE ADV_P_GD -------------------------------------------
CLL
CLL   PURPOSE:   CALCULATES ADVECTION INCREMENTS TO A FIELD AT A
CLL              SINGLE MODEL LEVEL USING AN EQUATION OF THE FORM(36).
CLL              NOT SUITABLE FOR SINGLE COLUMN USE.
CLL
CLL   VERSION FOR CRAY Y-MP
CLL
CLL   WRITTEN  BY M.H MAWSON.
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL   3.4    06/08/94 Vertical advection code restructured to improve
CLL                   parallel efficiency on C90.
CLL                   Authors: A. Dickinson, D. Salmond
CLL                   Reviewer: M. Mawson
CLL
CLL   3.4   23/06/94  DEF NOWHBR replaced by LOGICAL LWHITBROM
CLL                                                  S.J.Swarbrick
!     3.5    28/03/95 MPP code: Change updateable area and
!                     remove explicit wrap around calcs.  P.Burton
!     4.1    29/04/96 Remove MPP code (new ADVPGD1C version for MPP)
!                     and add TYPFLDPT arguments
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD B. VERSION 2, DATED 18/01/90
CLL
CLL   LOGICAL COMPONENTS COVERED: P121
CLL
CLL   PROJECT TASK: P1
CLL
CLL   DOCUMENTATION:       THE EQUATION USED IS (35)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER NO. 10
CLL                        M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLLEND-------------------------------------------------------------
C
C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE ADV_P_GD
     1                   (FIELD_LOWER,FIELD,FIELD_UPPER,U,V,
     1                   ETADOT_LOWER,ETADOT_UPPER,
     2                   SEC_P_LATITUDE,FIELD_INC,NUX,NUY,P_FIELD,
     3                   U_FIELD,ROW_LENGTH,
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
     4                   ADVECTION_TIMESTEP,
     5                   LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     6                   SEC_U_LATITUDE,BRSP,
     7                   L_SECOND,LWHITBROM)

      IMPLICIT NONE

      INTEGER
     *  P_FIELD             !IN DIMENSION OF FIELDS ON PRESSSURE GRID.
     *, U_FIELD             !IN DIMENSION OF FIELDS ON VELOCITY GRID
     *, ROW_LENGTH          !IN NUMBER OF POINTS PER ROW
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

      LOGICAL
     *  L_SECOND     ! SET TO TRUE IF NU_BASIC IS ZERO.
     * ,LWHITBROM    ! SWITCH FOR WHITE & BROMLEY TERMS

      REAL
     * U(U_FIELD)           !IN ADVECTING U FIELD, MASS-WEIGHTED.
     *,V(U_FIELD)           !IN ADVECTING V FIELD, MASS-WEIGHTED.
     *,ETADOT_UPPER(P_FIELD)!IN ADVECTING VERTICAL VELOC AT K+1/2,
     *                      !   MASS-WEIGHTED.
     *,ETADOT_LOWER(P_FIELD)!IN ADVECTING VERTICAL VELOC AT K-1/2,
     *                      !   MASS-WEIGHTED.
     *,FIELD(P_FIELD)       !IN FIELD TO BE ADVECTED.
     *,FIELD_UPPER(P_FIELD) !IN FIELD TO BE ADVECTED AT LEVEL + 1 .
     *,FIELD_LOWER(P_FIELD) !IN FIELD TO BE ADVECTED AT LEVEL - 1 .
     *,NUX(P_FIELD)   !IN HOLDS PARAMETER NU FOR EAST-WEST ADVECTION.
     *,NUY(P_FIELD)   !IN HOLDS PARAMETER NU FOR NORTH-SOUTH ADVECTION.
     *,SEC_P_LATITUDE(P_FIELD) !IN HOLDS 1/COS(PHI) AT P POINTS.
     *,SEC_U_LATITUDE(U_FIELD) !IN HOLDS 1/COS(PHI) AT U POINTS.
     *,ADVECTION_TIMESTEP   !IN
     *,LATITUDE_STEP_INVERSE  !IN 1/(DELTA PHI)
     *,LONGITUDE_STEP_INVERSE !IN 1/(DELTA LAMDA)

      REAL
     * BRSP(P_FIELD)  !IN BRSP TERM AT LEVEL (SEE DOC.PAPER NO 10)

      REAL
     * FIELD_INC(P_FIELD)   !OUT HOLDS INCREMENT TO FIELD.
C
C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: 3 ARE REQUIRED

      REAL
     * WORK(P_FIELD)        ! GENERAL WORK-SPACE.
     *,U_TERM(P_FIELD)      ! HOLDS U ADVECTION TERM FROM EQUATION (35)
     *,V_TERM(P_FIELD)      ! HOLDS V ADVECTION TERM FROM EQUATION (35)
C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES

C REAL SCALARS
      REAL
     * SCALAR1,SCALAR2

C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     *  I,IJ,IK,IL,IM,J

C*L   NO EXTERNAL SUBROUTINE CALLS:------------------------------------
C*---------------------------------------------------------------------

CL    MAXIMUM VECTOR LENGTH ASSUMED IS END_P_UPDATE-START_P_UPDATE+1
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL    SECTION 1.     CALCULATE U_TERM IN EQUATION (35).
CL---------------------------------------------------------------------

C----------------------------------------------------------------------
CL    SECTION 1.1    CALCULATE TERM U D(FIELD)/D(LAMDA).
C----------------------------------------------------------------------

C CALCULATE TERM AT ALL POINTS EXCEPT LAST AND STORE IN WORK.
! Loop over field, missing top and bottom rows
      DO 110 I=START_POINT_NO_HALO,END_P_POINT_NO_HALO-1
        WORK(I) = .5*(U(I)+U(I-ROW_LENGTH))*LONGITUDE_STEP_INVERSE*
     *             (FIELD(I+1) - FIELD(I))
 110  CONTINUE


C----------------------------------------------------------------------
CL    SECTION 1.2    CALCULATE U ADVECTION TERM IN EQUATION (35).
CL                   IF L_SECOND = TRUE PERFORM SECOND ORDER ADVECTION
CL                   ONLY.
C----------------------------------------------------------------------

      IF(L_SECOND) THEN
C LIMITED AREA MODEL.

! Loop over field, missing top and bottom rows and first and
! last points.
        DO J=START_POINT_NO_HALO+1,END_P_POINT_NO_HALO-1
          U_TERM(J) = .5*(WORK(J)+WORK(J-1))
        END DO

C CORNER VALUES

        U_TERM(START_POINT_NO_HALO)=0.0
        U_TERM(END_P_POINT_NO_HALO)=0.0

      ELSE
C LIMITED AREA MODEL. VALUES NOT CALCULATED AT FIRST,SECOND,NEXT TO LAST
C AND LAST ON A ROW.

! Loop over field, missing top and bottom rows and first two points
! and last two points.
        DO 120 J=START_POINT_NO_HALO+2,END_P_POINT_NO_HALO-2
          U_TERM(J) = (1.+NUX(J))*.5*(WORK(J)+WORK(J-1))-NUX(J)*.5*
     *                  (WORK(J+1)+WORK(J-2))
 120    CONTINUE

C CALCULATE  VALUES AT SECOND AND NEXT TO LAST POINTS ON A ROW.
C THESE VALUES ARE JUST SECOND ORDER.

! Loop over first point of every row, missing top and bottom rows
        DO 124 I=START_POINT_NO_HALO,END_P_POINT_NO_HALO,ROW_LENGTH
          IK = I+LAST_ROW_PT-2  ! IK is penultimate point of row
          IL = I + 1
C SECOND POINT.
          U_TERM(IL) = .5*(WORK(IL)+WORK(I))
C NEXT TO LAST POINT.
          U_TERM(IK) = .5*(WORK(IK)+WORK(IK-1))
 124    CONTINUE
C         CORNER VALUES
C
        U_TERM(START_POINT_NO_HALO)=0.0
        U_TERM(END_P_POINT_NO_HALO)=0.0
C

      END IF

CL
CL---------------------------------------------------------------------
CL    SECTION 2.     CALCULATE V_TERM IN EQUATION (35).
CL---------------------------------------------------------------------

C----------------------------------------------------------------------
CL    SECTION 2.1    CALCULATE TERM V D(FIELD)/D(PHI).
C----------------------------------------------------------------------

C CALCULATE TERM AT ALL POINTS EXCEPT FIRST AND STORE IN WORK.
! Loop over field, missing bottom row and first point of top row
      DO 210 I=START_POINT_NO_HALO-ROW_LENGTH+1,END_P_POINT_NO_HALO
        WORK(I) = .5*(V(I)+V(I-1))*LATITUDE_STEP_INVERSE*
     *             (FIELD(I) - FIELD(I+ROW_LENGTH))
 210  CONTINUE


C----------------------------------------------------------------------
CL    SECTION 2.2    CALCULATE V ADVECTION TERM IN EQUATION (35).
CL                   IF L_SECOND = TRUE PERFORM SECOND ORDER ADVECTION
CL                   ONLY.
C----------------------------------------------------------------------

      IF(L_SECOND) THEN
C LIMITED AREA MODEL.

! Loop over field, missing top and bottom rows and first and last points
        DO I=START_POINT_NO_HALO+1,END_P_POINT_NO_HALO-1
          V_TERM(I) = .5*(WORK(I-ROW_LENGTH)+WORK(I))
        END DO

        V_TERM(START_POINT_NO_HALO)=0.0
        V_TERM(END_P_POINT_NO_HALO)=0.0

      ELSE
C LIMITED AREA MODEL.
C CALCULATE ALL VALUES EXCEPT ON ROWS NEXT TO BOUNDARIES.

! Loop over field missing top and bottom rows
        DO 220 I=START_POINT_NO_HALO+ROW_LENGTH,
     &           END_P_POINT_NO_HALO-ROW_LENGTH
          V_TERM(I) = (1.+NUY(I))*.5*(WORK(I-ROW_LENGTH)+WORK(I)) -
     *             NUY(I)*.5*(WORK(I+ROW_LENGTH)+WORK(I-2*ROW_LENGTH))
 220    CONTINUE

C CALCULATE VALUES ON SLICES NEXT TO BOUNDARIES AS SECOND ORDER.

        DO 222 I=2,ROW_LENGTH-1
          IJ = END_P_POINT_NO_HALO-ROW_LENGTH+I
          IK = START_POINT_NO_HALO+I-1
C NEXT TO NORTHERN BOUNDARY.
          V_TERM(IK) = .5*(WORK(IK-ROW_LENGTH)+WORK(IK))
C NEXT TO SOUTHERN BOUNDARY.
          V_TERM(IJ) = .5*(WORK(IJ-ROW_LENGTH)+WORK(IJ))
 222    CONTINUE
        V_TERM(START_POINT_NO_HALO) = 0.0
        V_TERM(START_POINT_NO_HALO+ROW_LENGTH-1)=0.0
C
        V_TERM(END_P_POINT_NO_HALO) = 0.0
        V_TERM(END_P_POINT_NO_HALO-ROW_LENGTH+1) = 0.0

      END IF

CL
CL---------------------------------------------------------------------
CL    SECTION 3.     CALCULATE VERTICAL FLUX AND COMBINE WITH U AND V
CL                   TERMS TO FORM INCREMENT.
CL---------------------------------------------------------------------

CL    VERTICAL FLUX ON INPUT IS .5*TIMESTEP*ETADOT*D(FIELD)/D(ETA)
CL    AT LEVEL K-1/2. AT THE END OF THIS SECTION IT IS THE SAME
CL    QUANTITY BUT AT LEVEL K+1/2.

! Loop over field missing top and bottom rows.
      DO 300 I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
        SCALAR1 = .5 * ADVECTION_TIMESTEP *
     *         ETADOT_UPPER(I) * (FIELD_UPPER(I) - FIELD(I))
        SCALAR2 = .5 * ADVECTION_TIMESTEP *
     *         ETADOT_LOWER(I) * (FIELD(I) - FIELD_LOWER(I))
        FIELD_INC(I) = ADVECTION_TIMESTEP * SEC_P_LATITUDE(I) *
     *                  (U_TERM(I)+V_TERM(I))
     &                   + SCALAR1+SCALAR2
      IF (LWHITBROM) THEN
        FIELD_INC(I) = FIELD_INC(I)
     *                  + FIELD(I)*BRSP(I)
      END IF
 300  CONTINUE


CL   LIMITED AREA MODEL SET BOUNDARY INCREMENTS
CL   TO ZERO.

! Loop over first point of each row, missing top and bottom rows.
      DO 310 I=START_POINT_NO_HALO,END_P_POINT_NO_HALO,ROW_LENGTH
        FIELD_INC(I) = 0.
        FIELD_INC(I+ROW_LENGTH-1) = 0.
 310  CONTINUE


CL    END OF ROUTINE ADV_P_GD

      RETURN
      END
