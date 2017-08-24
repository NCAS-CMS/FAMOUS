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
CLL   SUBROUTINE ADV_U_GD -------------------------------------------
CLL
CLL   PURPOSE:   CALCULATES ADVECTION INCREMENTS TO A FIELD AT A
CLL              SINGLE MODEL LEVEL USING AN EQUATION OF THE FORM(38).
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL
CLL   VERSION FOR CRAY Y-MP
CLL
CLL   WRITTEN BY M.H MAWSON.
CLL   MPP CODE ADDED BY P.BURTON
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL 4.1      29/11/95 New version of routine specifically for MPP
CLL                   Fourth order MPP code by Roar Skalin
CLL                                                      P.Burton
!LL 4.2      10/01/97 Amend calculation to prevent different answers
!LL                   with different compiler options. D. Robinson.
!LL 4.4      24/04/97 Removed QAN fix.   P.Burton
C     vn4.3    Mar. 97   T3E migration : optimisation changes
C                                       D.Salmond
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD B. VERSION 2, DATED 18/01/90
CLL
CLL   LOGICAL COMPONENTS COVERED: P122
CLL
CLL   PROJECT TASK: P1
CLL
CLL   DOCUMENTATION:       THE EQUATION USED IS (37)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER NO. 10
CLL                        M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLLEND-------------------------------------------------------------
CLL
C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE ADV_U_GD
     1                   (FIELD_LOWER,FIELD,FIELD_UPPER,U,V,
     1                    ETADOT_LOWER,ETADOT_UPPER,
     2                    SEC_U_LATITUDE,FIELD_INC,NUX,NUY,U_FIELD,
     3                    ROW_LENGTH,
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
     4                    ADVECTION_TIMESTEP,
     5                    LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     6                    SEC_P_LATITUDE,BRSP,
     7                    L_SECOND,LWHITBROM,
     &                    extended_FIELD,extended_U_FIELD,
     &                    extended_address)              

      IMPLICIT NONE

      INTEGER
     *  U_FIELD             !IN DIMENSION OF FIELDS ON VELOCITY GRID
     &, extended_U_FIELD    !IN DIMENSION of U fields with extra halo
     *, ROW_LENGTH          !IN NUMBER OF POINTS PER ROW

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
     * U(extended_U_FIELD)  !IN ADVECTING U FIELD MASS-WEIGHTED AND
     *                      ! HELD AT P POINTS. FIRST POINT OF FIELD
     *                      ! IS FIRST P POINT ON SECOND ROW OF P-GRID.
     *,V(extended_U_FIELD)  !IN ADVECTING V FIELD MASS-WEIGHTED AND
     *                      ! HELD AT P POINTS. FIRST POINT OF FIELD
     *                      ! IS FIRST P POINT ON SECOND ROW OF P-GRID.
     *,ETADOT_UPPER(U_FIELD)!IN ADVECTING VERTICAL VELOC AT K+1/2,
     *                      ! MASS-WEIGHTED.
     *,ETADOT_LOWER(U_FIELD)!IN ADVECTING VERTICAL VELOC AT K-1/2,
     *                      !   MASS-WEIGHTED.
     *,FIELD(U_FIELD)       !IN FIELD TO BE ADVECTED.
     *,FIELD_UPPER(U_FIELD) !IN FIELD TO BE ADVECTED AT LEVEL + 1 .
     *,FIELD_LOWER(U_FIELD) !IN FIELD TO BE ADVECTED AT LEVEL - 1 .
     *,NUX(U_FIELD)   !IN HOLDS PARAMETER NU FOR EAST-WEST ADVECTION.
     *,NUY(U_FIELD)   !IN HOLDS PARAMETER NU FOR NORTH-SOUTH ADVECTION.
     *,SEC_U_LATITUDE(U_FIELD) !IN HOLDS 1/COS(PHI) AT U POINTS.
     *,SEC_P_LATITUDE(U_FIELD) !IN HOLDS 1/COS(PHI) AT P POINTS.
     *,ADVECTION_TIMESTEP   !IN
     *,LATITUDE_STEP_INVERSE   !IN 1/(DELTA PHI)
     *,LONGITUDE_STEP_INVERSE  !IN 1/(DELTA LAMDA)

      REAL
     * BRSP(U_FIELD)  !IN BRSP TERM AT LEVEL+1/2 (SEE DOC.PAPER
     *                      ! NO 10)

      REAL
     * FIELD_INC(U_FIELD)   !OUT HOLDS INCREMENT TO FIELD.

      REAL
     * VERTICAL_FLUX(U_FIELD) !INOUT HOLDS VERTICAL FLUX OF FIELD
     *                        ! BETWEEN TWO LEVELS.

      REAL
     & extended_FIELD(extended_U_FIELD) ! IN field to be advected with
!                                       !    extra halos for 4th order
      INTEGER extended_address(U_FIELD)

C LOGICAL VARIABLE
      LOGICAL
     *  L_SECOND     ! SET TO TRUE IF NU_BASIC IS ZERO.
     * ,LWHITBROM    ! Switch for White & Bromley terms
C

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: 3 ARE REQUIRED

      REAL
     * WORK(U_FIELD)      ! GENERAL WORK-SPACE.
     *,U_TERM(U_FIELD)    ! HOLDS U ADVECTION TERM FROM EQUATION (37)
     *,V_TERM(U_FIELD)    ! HOLDS V ADVECTION TERM FROM EQUATION (37)
C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES

C REAL SCALARS
      REAL
     * SCALAR1,SCALAR2

C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     *  I,J

! Work space and scalars for the MPP Fourth Order Advection
       INTEGER  extended_index,  ! index for position in extended array
     &          extended_START_POINT_NO_HALO,
!                                ! start pos in extended array
     &          extended_END_U_POINT_NO_HALO,
!                                ! end pos in extended array
     &          extended_ROW_LENGTH,! row length of extended array
     &          I_start,I_end  ! loop bounds for 4th order advection

      REAL
     &          extended_WORK(extended_U_FIELD)  ! extended work space


C*L   NO EXTERNAL SUBROUTINE CALLS:------------------------------------
C*---------------------------------------------------------------------

CL    MAXIMUM VECTOR LENGTH ASSUMED IS
CL    END_U_POINT_NO_HALO-START_U_UPDATE+1
CL---------------------------------------------------------------------

      IF (.NOT. L_SECOND) THEN
! Calculate indexes in extended_arrays

        extended_ROW_LENGTH=ROW_LENGTH+2*extra_EW_Halo

        extended_START_POINT_NO_HALO=
     &    extended_address(START_POINT_NO_HALO)

        extended_END_U_POINT_NO_HALO=
     &    extended_address(END_U_POINT_NO_HALO)

      ENDIF

CL---------------------------------------------------------------------
CL    SECTION 1.     CALCULATE U_TERM IN EQUATION (37).
CL---------------------------------------------------------------------

C----------------------------------------------------------------------
CL    SECTION 1.1    CALCULATE TERM U D(FIELD)/D(LAMDA).
C----------------------------------------------------------------------

C CALCULATE TERM AT ALL POINTS EXCEPT LAST AND STORE IN WORK.

      IF (L_SECOND) THEN
! Loop over field, missing top and bottom rows, halos and last point
        DO 110 I=START_POINT_NO_HALO,END_U_POINT_NO_HALO-1
          WORK(I) = .5*(U(I+1)+U(I+1-ROW_LENGTH))*
     &              LONGITUDE_STEP_INVERSE*(FIELD(I+1) - FIELD(I))
 110    CONTINUE


      ELSE ! fourth order
! Loop over extended field, missing top and bottom rows and halos rows
        DO I=extended_START_POINT_NO_HALO-1,
     &       extended_END_U_POINT_NO_HALO+1
          extended_WORK(I) = 0.5*(U(I+1)+U(I+1-extended_ROW_LENGTH))*
     &                       LONGITUDE_STEP_INVERSE*
     &                       (extended_FIELD(I+1)-extended_FIELD(I))
        ENDDO

      ENDIF ! IF (L_SECOND)


C----------------------------------------------------------------------
CL    SECTION 1.2    CALCULATE U ADVECTION TERM IN EQUATION (37).
CL                   IF L_SECOND=TRUE ONLY DO SECOND ORDER ADVECTION.
C----------------------------------------------------------------------

      IF(L_SECOND) THEN
C LIMITED AREA MODEL. VALUES NOT CALCULATED AT FIRST,SECOND,NEXT TO LAST
C AND LAST ON A ROW.

! Loop over field, missing top and bottom rows and halos, and
! first and last points.
        DO J=START_POINT_NO_HALO+1,END_U_POINT_NO_HALO-1
          U_TERM(J) = .5*(WORK(J)+WORK(J-1))
        END DO

C CORNER VALUES
        U_TERM(START_POINT_NO_HALO)=0.0
        U_TERM(END_U_POINT_NO_HALO)=0.0
      ELSE   ! Fourth order

C LOOP OVER ALL POINTS BUT DON'T DO FIRST,SECOND AND LAST ON A ROW AS
C THEY NEED SPECIAL TREATMENT DUE TO FOURTH ORDER SCHEME.

! Loop over field, missing top and bottom rows and halos, and
! first point.
        DO 120 J=START_POINT_NO_HALO+1,END_U_POINT_NO_HALO
!        DO 120 J=START_POINT_NO_HALO+2,END_U_POINT_NO_HALO-1
          extended_index=extended_address(J)

          U_TERM(J) = (1.+NUX(J))*.5*(extended_WORK(extended_index)+
     &                                extended_WORK(extended_index-1))
     &                   -NUX(J) *.5*(extended_WORK(extended_index+1)+
     &                                extended_WORK(extended_index-2))
 120    CONTINUE

C LIMITED AREA MODEL.
C CALCULATE  VALUES AT SECOND AND NEXT TO LAST POINTS ON A ROW.
C THESE VALUES ARE JUST SECOND ORDER.

        IF (at_left_of_LPG) THEN
! Do second point along each row
          DO I=START_POINT_NO_HALO+FIRST_ROW_PT,END_U_POINT_NO_HALO,
     &         ROW_LENGTH
            extended_index=extended_address(I)

            U_TERM(I)= 0.5*(extended_WORK(extended_index)+
     &                      extended_WORK(extended_index-1))
          ENDDO
        ENDIF

! Do penultimate point along each row

        IF (at_right_of_LPG) THEN
          DO I=START_POINT_NO_HALO+LAST_ROW_PT-2,END_U_POINT_NO_HALO,
     &         ROW_LENGTH
            extended_index=extended_address(I)

            U_TERM(I)= 0.5*(extended_WORK(extended_index)+
     &                      extended_WORK(extended_index-1))
          ENDDO
        ENDIF

        U_TERM(START_POINT_NO_HALO)=0
        U_TERM(END_U_POINT_NO_HALO)=0

      END IF  ! IF (L_SECOND)

CL
CL---------------------------------------------------------------------
CL    SECTION 2.     CALCULATE V_TERM IN EQUATION (37).
CL---------------------------------------------------------------------

C----------------------------------------------------------------------
CL    SECTION 2.1    CALCULATE TERM V D(FIELD)/D(PHI).
C----------------------------------------------------------------------

C CALCULATE TERM AT ALL POINTS EXCEPT LAST AND STORE IN WORK.

      IF (L_SECOND) THEN

! Loop over field, mssing bottom row and top and bottom halos and last
! point.
        DO 210 I=START_POINT_NO_HALO-ROW_LENGTH,
     &           END_U_POINT_NO_HALO-1
          WORK(I) = .5*(V(I)+V(I+1))*LATITUDE_STEP_INVERSE*
     *               (FIELD(I) - FIELD(I+ROW_LENGTH))
 210    CONTINUE


      ELSE ! Fourth order
! Calculate WORK at the Southern halo too. This is needed for the
! computation of the Southern row

!        DO I=extended_START_POINT_NO_HALO-2*extended_ROW_LENGTH,
!     &       extended_END_U_POINT_NO_HALO+extended_ROW_LENGTH
         IF (at_top_of_LPG) THEN
           I_start=extended_address(TOP_ROW_START)
         ELSE
           I_start=extended_START_POINT_NO_HALO-2*extended_ROW_LENGTH
         ENDIF
         IF (at_base_of_LPG) THEN
           I_end=extended_END_U_POINT_NO_HALO-1
         ELSE
           I_end=extended_END_U_POINT_NO_HALO-1+extended_ROW_LENGTH
         ENDIF
         DO I=I_start,I_end
          extended_WORK(I)=0.5*(V(I)+V(I+1))*LATITUDE_STEP_INVERSE*
     &      (extended_FIELD(I)-extended_FIELD(I+extended_ROW_LENGTH))
        ENDDO
        extended_WORK(I_end+1)=extended_WORK(I_end)

      ENDIF ! L_SECOND


C----------------------------------------------------------------------
CL    SECTION 2.2    CALCULATE V ADVECTION TERM IN EQUATION (37).
CL                   IF L_SECOND=TRUE ONLY DO SECOND ORDER ADVECTION.
C----------------------------------------------------------------------

      IF(L_SECOND) THEN
C LIMITED AREA MODEL.
C CALCULATE ALL VALUES EXCEPT ON ROWS NEXT TO BOUNDARIES.

! Loop over field, missing top and bottom rows and halos, and
! first and last points.
        DO I=START_POINT_NO_HALO+1,END_U_POINT_NO_HALO-1
          V_TERM(I) = .5*(WORK(I-ROW_LENGTH)+WORK(I))
        END DO

C SET LAST POINT OF LOOP TO ZERO.
        V_TERM(END_U_POINT_NO_HALO-ROW_LENGTH) = 0.0

C CORNER VALUES

        V_TERM(START_POINT_NO_HALO)=0.0
        V_TERM(END_U_POINT_NO_HALO)=0.0

      ELSE  ! Fourth order
C LIMITED AREA MODEL.
C CALCULATE ALL VALUES EXCEPT ON ROWS NEXT TO BOUNDARIES.

! Loop over field, missing top and bottom rows and halos
        IF (at_top_of_LPG) THEN
          I_start=START_POINT_NO_HALO+ROW_LENGTH
        ELSE
          I_start=START_POINT_NO_HALO
        ENDIF
        IF (at_base_of_LPG) THEN
          I_end=END_U_POINT_NO_HALO-ROW_LENGTH
        ELSE
          I_end=END_U_POINT_NO_HALO
        ENDIF
!        DO I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
        DO I=I_start,I_end
          extended_index=extended_address(I)

          V_TERM(I) = (1.0+NUY(I))*0.5*
     &     (extended_WORK(extended_index-extended_ROW_LENGTH)
     &    + extended_WORK(extended_index))
     &                   - NUY(I) *0.5*
     &     (extended_WORK(extended_index+extended_ROW_LENGTH)
     &    + extended_WORK(extended_index-2*extended_ROW_LENGTH))
        ENDDO

C CALCULATE VALUES ON SLICES NEXT TO BOUNDARIES AS SECOND ORDER.
        IF (at_top_of_LPG) THEN
! Loop over row beneath top row
          DO I=START_POINT_NO_HALO+FIRST_ROW_PT-1,
     &         START_POINT_NO_HALO+LAST_ROW_PT-1
            extended_index=extended_address(I)
            V_TERM(I)=0.5*
     &        (extended_WORK(extended_index-extended_ROW_LENGTH)
     &       + extended_WORK(extended_index))
          ENDDO
        ENDIF

        IF (at_base_of_LPG) THEN
          DO I=END_U_POINT_NO_HALO-ROW_LENGTH+FIRST_ROW_PT,
     &         END_U_POINT_NO_HALO-ROW_LENGTH+LAST_ROW_PT
            extended_index=extended_address(I)
            V_TERM(I)=0.5*
     &        (extended_WORK(extended_index-extended_ROW_LENGTH)
     &       + extended_WORK(extended_index))
          ENDDO
        ENDIF

C         CORNER VALUES
C
        V_TERM(START_POINT_NO_HALO)=0.0
        V_TERM(START_POINT_NO_HALO+ROW_LENGTH-1)=0.0
        V_TERM(END_U_POINT_NO_HALO-ROW_LENGTH+1)=0.0
        V_TERM(END_U_POINT_NO_HALO)=0.0

      END IF

CL
CL---------------------------------------------------------------------
CL    SECTION 3.     CALCULATE VERTICAL FLUX AND COMBINE WITH U AND V
CL                   TERMS TO FORM INCREMENT.
CL---------------------------------------------------------------------

CL    VERTICAL FLUX ON INPUT IS .5*TIMESTEP*ETADOT*D(FIELD)/D(ETA)
CL    AT LEVEL K-1/2. AT THE END OF THEIS SECTION IT IS THE SAME
CL    QUANTITY BUT AT LEVEL K+1/2.

! Loop over field, missing top and bottom rows and halos
      DO 300 I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
        SCALAR1 = .5 * ADVECTION_TIMESTEP *
     *         ETADOT_UPPER(I) * (FIELD_UPPER(I) - FIELD(I))
        SCALAR2 = .5 * ADVECTION_TIMESTEP *
     *         ETADOT_LOWER(I) * (FIELD(I) - FIELD_LOWER(I))
        FIELD_INC(I) = ADVECTION_TIMESTEP * SEC_U_LATITUDE(I) *
     *                  (U_TERM(I)+V_TERM(I))
     &                   + SCALAR1+SCALAR2
      IF (LWHITBROM) THEN
        FIELD_INC(I) = FIELD_INC(I)
     *                  + FIELD(I)*BRSP(I)
      END IF
 300  CONTINUE


CL    LIMITED AREA MODEL SET BOUNDARY INCREMENTS TO ZERO.

        IF (at_left_of_LPG) THEN
          DO I=START_POINT_NO_HALO+FIRST_ROW_PT-1,
     &         END_U_POINT_NO_HALO,ROW_LENGTH
            FIELD_INC(I)=0.0
          ENDDO
        ENDIF

        IF (at_right_of_LPG) THEN
          DO I=START_POINT_NO_HALO+LAST_ROW_PT-2,
     &         END_U_POINT_NO_HALO,ROW_LENGTH
            FIELD_INC(I)=0.0
            FIELD_INC(I+1)=0.0
          ENDDO
        ENDIF


CL    END OF ROUTINE ADV_U_GD

      RETURN
      END
