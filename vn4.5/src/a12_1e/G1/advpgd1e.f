C *****************************COPYRIGHT******************************
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
CLL   SUBROUTINE ADV_P_GD -------------------------------------------
CLL
CLL   PURPOSE:   CALCULATES ADVECTION INCREMENTS TO A FIELD AT A
CLL              SINGLE MODEL LEVEL USING AN EQUATION OF THE FORM(36).
CLL              NOT SUITABLE FOR SINGLE COLUMN USE.
CLL
CLL   WAS VERSION FOR CRAY Y-MP
CLL
CLL   WRITTEN  BY M.H MAWSON.
CLL   MPP CODE ADDED BY P.BURTON
CLL
CLL  Model            Modification history:
CLL version  Date
!LL   4.4   11/08/97  New version optimised for T3E.
!LL                   Not bit-reproducible with ADVPGD1C.
CLL    4.4   04/08/97 Optimisation for T3E   D.Salmond
CLL    4.5   31/03/98 Correct uninitialised value of U_TERM which can
CLL                   cause failures for LAM with 4th order advection.
CLL                   R. Rawlins.
CLL
CLL    4.5   29/4/98   T3E Optimisation for MES D.Salmond
CLL
CLL   PROGRAMMING STANDARD:
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
     1                   (P_LEVELS,FIELD,U,V,
     1                   ETADOT,
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
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
     4                   ADVECTION_TIMESTEP,
     5                   LATITUDE_STEP_INVERSE,LONGITUDE_STEP_INVERSE,
     6                   SEC_U_LATITUDE,BRSP,
     7                   L_SECOND,LWHITBROM,
     &                   extended_FIELD,
     &                   extended_P_FIELD,extended_U_FIELD,
     &                   extended_address)

      IMPLICIT NONE

      INTEGER
     *  P_LEVELS
     *, P_FIELD             !IN DIMENSION OF FIELDS ON PRESSSURE GRID.
     *, U_FIELD             !IN DIMENSION OF FIELDS ON VELOCITY GRID
     &, extended_P_FIELD    !IN DIMESNION of P fields with extra halo
     &, extended_U_FIELD    !IN DIMESNION of U fields with extra halo
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

      LOGICAL
     *  L_SECOND     ! SET TO TRUE IF NU_BASIC IS ZERO.
     * ,LWHITBROM    ! SWITCH FOR WHITE & BROMLEY TERMS

      REAL
     * U(extended_U_FIELD,P_LEVELS)
!               !IN ADVECTING U FIELD, MASS-WEIGHTED.
     *,V(extended_U_FIELD,P_LEVELS)
!               !IN ADVECTING V FIELD, MASS-WEIGHTED.
     *,ETADOT(P_FIELD,P_LEVELS)!IN ADVECTING VERTICAL VELOC AT K+1/2,
     *                      !   MASS-WEIGHTED.
     *,FIELD(P_FIELD,P_LEVELS)       !IN FIELD TO BE ADVECTED.
     *,NUX(P_FIELD,P_LEVELS)
!               !IN HOLDS PARAMETER NU FOR EAST-WEST ADVECTION.
     *,NUY(P_FIELD,P_LEVELS)
!               !IN HOLDS PARAMETER NU FOR NORTH-SOUTH ADVECTION.
     *,SEC_P_LATITUDE(P_FIELD) !IN HOLDS 1/COS(PHI) AT P POINTS.
     *,SEC_U_LATITUDE(U_FIELD) !IN HOLDS 1/COS(PHI) AT U POINTS.
     *,ADVECTION_TIMESTEP   !IN
     *,LATITUDE_STEP_INVERSE  !IN 1/(DELTA PHI)
     *,LONGITUDE_STEP_INVERSE !IN 1/(DELTA LAMDA)

      REAL
     * BRSP(P_FIELD,P_LEVELS)
!               !IN BRSP TERM AT LEVEL (SEE DOC.PAPER NO 10)

      REAL
     * FIELD_INC(P_FIELD,P_LEVELS)   !OUT HOLDS INCREMENT TO FIELD.

      REAL
     & extended_FIELD(extended_P_FIELD,P_LEVELS)
!                     ! IN field to be advected with
!                     !    extra halos for 4th order
      INTEGER extended_address(P_FIELD)
C

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
     *  I,IJ,IK,IL,IM,J,K
      INTEGER START_ADD_base,START_ADD_top

! Work space and scalars for the MPP Fourth Order Advection
       INTEGER  info,            ! return code from comms operations
     &          extended_index,  ! index for position in extended array
     &          extended_START_POINT_NO_HALO,
!                                ! start position in extended array
     &          extended_END_P_POINT_NO_HALO,
!                                ! end position in extended array
     &          extended_ROW_LENGTH    ! row length of extended array

      REAL
     &          rot_work(ROW_LENGTH), ! work space for rotated pole rows
     &          extended_WORK(extended_P_FIELD)  ! extended work space


C*L   NO EXTERNAL SUBROUTINE CALLS:------------------------------------
C*---------------------------------------------------------------------

CL    MAXIMUM VECTOR LENGTH ASSUMED IS
CL    END_P_POINT_NO_HALO-START_POINT_NO_HALO+1
CL---------------------------------------------------------------------

      IF(L_SECOND) THEN
!  SECOND ORDER ADEVCTION

      DO K=1,P_LEVELS

CL
CL---------------------------------------------------------------------
CL    SECTION 1.     CALCULATE U_TERM IN EQUATION (35).
CL---------------------------------------------------------------------

C----------------------------------------------------------------------
CL    SECTION 1.1    CALCULATE TERM U D(FIELD)/D(LAMDA).
C----------------------------------------------------------------------

C----------------------------------------------------------------------
CL    SECTION 1.2    CALCULATE U ADVECTION TERM IN EQUATION (35).
CL                   IF L_SECOND = TRUE PERFORM SECOND ORDER ADVECTION
CL                   ONLY.
C----------------------------------------------------------------------

CL
CL---------------------------------------------------------------------
CL    SECTION 2.     CALCULATE V_TERM IN EQUATION (35).
CL---------------------------------------------------------------------

C----------------------------------------------------------------------
CL    SECTION 2.1    CALCULATE TERM V D(FIELD)/D(PHI).
C----------------------------------------------------------------------

C----------------------------------------------------------------------
CL    SECTION 2.2    CALCULATE V ADVECTION TERM IN EQUATION (35).
CL                   IF L_SECOND = TRUE PERFORM SECOND ORDER ADVECTION
CL                   ONLY.
C----------------------------------------------------------------------

CL
CL---------------------------------------------------------------------
CL    SECTION 3.     CALCULATE VERTICAL FLUX AND COMBINE WITH U AND V
CL                   TERMS TO FORM INCREMENT.
CL---------------------------------------------------------------------

CL    VERTICAL FLUX ON INPUT IS .5*TIMESTEP*ETADOT*D(FIELD)/D(ETA)
CL    AT LEVEL K-1/2. AT THE END OF THIS SECTION IT IS THE SAME
CL    QUANTITY BUT AT LEVEL K+1/2.

! Loop over field, missing top and bottom rows and halos
      if(k.ne.1.and.k.ne.P_LEVELS)then

cdir$ unroll4
      DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
        SCALAR1 = .5 * ADVECTION_TIMESTEP *
     *         ETADOT(I,K+1) * (FIELD(I,K+1) - FIELD(I,K))
        SCALAR2 = WORK(I)
        FIELD_INC(I,K) = SCALAR1 +SCALAR2
      IF (LWHITBROM) FIELD_INC(I,K) = FIELD_INC(I,K)
     *                  + FIELD(I,K)*BRSP(I,K)
        WORK(I)=SCALAR1
      ENDDO
      else if(k.eq.1) then
cdir$ unroll4
      DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
        SCALAR1 = .5 * ADVECTION_TIMESTEP *
     *         ETADOT(I,K+1) * (FIELD(I,K+1) - FIELD(I,K))
        FIELD_INC(I,K) = SCALAR1
      IF (LWHITBROM) FIELD_INC(I,K) = FIELD_INC(I,K)
     *                  + FIELD(I,K)*BRSP(I,K)
      WORK(I)=SCALAR1
      END DO
      else if(k.eq.P_LEVELS) then
cdir$ unroll4
      DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
        SCALAR2 = WORK(I)
        FIELD_INC(I,K) =  SCALAR2
      IF (LWHITBROM) FIELD_INC(I,K) = FIELD_INC(I,K)
     *                  + FIELD(I,K)*BRSP(I,K)
      END DO
      endif ! if(k.ne.1.and.k.ne.P_LEVELS)then

      DO I=START_POINT_NO_HALO+1,END_P_POINT_NO_HALO-1
        FIELD_INC(I,K) = FIELD_INC(I,K) +
     *               .25*ADVECTION_TIMESTEP * SEC_P_LATITUDE(I) *
     *                  (LONGITUDE_STEP_INVERSE*
     *                    ((U(I,K)+U(I-ROW_LENGTH,K))*
     *                           (FIELD(I+1,K)-FIELD(I,K))+
     *                     (U(I-1,K)+U(I-1-ROW_LENGTH,K))*
     *                           (FIELD(I,K)-FIELD(I-1,K)))
     *                 +
     *                   LATITUDE_STEP_INVERSE*
     *                    ((V(I-ROW_LENGTH,K)+V(I-1-ROW_LENGTH,K))*
     *               (FIELD(I-ROW_LENGTH,K) - FIELD(I,K))+
     &                     (V(I,K)+V(I-1,K))*
     *               (FIELD(I,K) - FIELD(I+ROW_LENGTH,K))))
      ENDDO




CL   LIMITED AREA MODEL SET BOUNDARY INCREMENTS
CL   TO ZERO.

       IF (at_left_of_LPG) THEN
          DO I=START_POINT_NO_HALO+FIRST_ROW_PT-1,
     &         END_P_POINT_NO_HALO,ROW_LENGTH
            FIELD_INC(I,K)=0.
          ENDDO
        ENDIF

        IF (at_right_of_LPG) THEN
          DO I=START_POINT_NO_HALO+LAST_ROW_PT-1,
     &         END_P_POINT_NO_HALO,ROW_LENGTH
            FIELD_INC(I,K)=0.
          ENDDO
      ENDIF

      ENDDO

      ELSE  !IF(L_SECOND)
!  FOURTH ORDER ADEVCTION

! Calculate indexes in extended_arrays

      extended_ROW_LENGTH=ROW_LENGTH+2*extra_EW_Halo

        extended_START_POINT_NO_HALO=
     &    extended_address(START_POINT_NO_HALO)

        extended_END_P_POINT_NO_HALO=
     &    extended_address(END_P_POINT_NO_HALO)


      DO K=1,P_LEVELS

CL
CL---------------------------------------------------------------------
CL    SECTION 1.     CALCULATE U_TERM IN EQUATION (35).
CL---------------------------------------------------------------------

C----------------------------------------------------------------------
CL    SECTION 1.1    CALCULATE TERM U D(FIELD)/D(LAMDA).
C----------------------------------------------------------------------

C CALCULATE TERM AT ALL POINTS EXCEPT LAST AND STORE IN WORK.

! Loop over extended field, missing top and bottom rows and halos rows
        DO I=extended_START_POINT_NO_HALO-1,
     &       extended_END_P_POINT_NO_HALO+1
          extended_WORK(I)=0.5*(U(I,K)+U(I-extended_ROW_LENGTH,K))*
     &                     LONGITUDE_STEP_INVERSE*
     &                     (extended_FIELD(I+1,K)-extended_FIELD(I,K))
        ENDDO



C----------------------------------------------------------------------
CL    SECTION 1.2    CALCULATE U ADVECTION TERM IN EQUATION (35).
CL                   IF L_SECOND = TRUE PERFORM SECOND ORDER ADVECTION
CL                   ONLY.
C----------------------------------------------------------------------


C LOOP OVER ALL POINTS.

! Loop over field, missing top and bottom rows and halos, and
! first point.
        DO 120 J=START_POINT_NO_HALO+1,END_P_POINT_NO_HALO
          extended_index=extended_address(J)

          U_TERM(J) = (1.+NUX(J,K))*.5*(extended_WORK(extended_index)+
     &                                extended_WORK(extended_index-1))
     &                 -  NUX(J,K) *.5*(extended_WORK(extended_index+1)+
     &                                extended_WORK(extended_index-2))
 120    CONTINUE


! Initialise first value to avoid potential flop exception failure

        U_TERM(START_POINT_NO_HALO)= 0.0

C CALCULATE  VALUES AT SECOND AND NEXT TO LAST POINTS ON A ROW.
C THESE VALUES ARE JUST SECOND ORDER.

        IF (at_left_of_LPG) THEN
! Do second point along each row
          DO I=START_POINT_NO_HALO+FIRST_ROW_PT,END_P_POINT_NO_HALO,
     &         ROW_LENGTH
            extended_index=extended_address(I)

            U_TERM(I)= 0.5*(extended_WORK(extended_index)+
     &                      extended_WORK(extended_index-1))
          ENDDO
        ENDIF

! Do penultimate point along each row

        IF (at_right_of_LPG) THEN
          DO I=START_POINT_NO_HALO+LAST_ROW_PT-2,END_P_POINT_NO_HALO,
     &         ROW_LENGTH
            extended_index=extended_address(I)

            U_TERM(I)= 0.5*(extended_WORK(extended_index)+
     &                      extended_WORK(extended_index-1))
          ENDDO
        ENDIF


CL
CL---------------------------------------------------------------------
CL    SECTION 2.     CALCULATE V_TERM IN EQUATION (35).
CL---------------------------------------------------------------------

C----------------------------------------------------------------------
CL    SECTION 2.1    CALCULATE TERM V D(FIELD)/D(PHI).
C----------------------------------------------------------------------

C CALCULATE TERM AT ALL POINTS EXCEPT FIRST AND STORE IN WORK.

! Calculate WORK at the Southern halo too. This is needed for the
! computation of the Southern row

        DO I=extended_START_POINT_NO_HALO-2*extended_ROW_LENGTH,
     &       extended_END_P_POINT_NO_HALO+extended_ROW_LENGTH
         extended_WORK(I)=0.5*(V(I,K)+V(I-1,K))*LATITUDE_STEP_INVERSE*
     &   (extended_FIELD(I,K)-extended_FIELD(I+extended_ROW_LENGTH,K))
        ENDDO


C----------------------------------------------------------------------
CL    SECTION 2.2    CALCULATE V ADVECTION TERM IN EQUATION (35).
CL                   IF L_SECOND = TRUE PERFORM SECOND ORDER ADVECTION
CL                   ONLY.
C----------------------------------------------------------------------

C LIMITED AREA MODEL.
! Calculate all values except on rows next to poles and next to the
! processor interfaces

! Loop over field, missing top and bottom rows and halos
        DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
          extended_index=extended_address(I)

          V_TERM(I) = (1.0+NUY(I,K))*0.5*
     &     (extended_WORK(extended_index-extended_ROW_LENGTH)
     &    + extended_WORK(extended_index))
     &                   - NUY(I,K) *0.5*
     &     (extended_WORK(extended_index+extended_ROW_LENGTH)
     &    + extended_WORK(extended_index-2*extended_ROW_LENGTH))
        ENDDO


C CALCULATE VALUES ON SLICES NEXT TO BOUNDARIES AS SECOND ORDER.

        IF (at_top_of_LPG) THEN
! Loop over row beneath top row, missing halos
          DO I=START_POINT_NO_HALO+FIRST_ROW_PT-1,
     &         START_POINT_NO_HALO+LAST_ROW_PT-1
            extended_index=extended_address(I)

            V_TERM(I)=0.5*
     &        (extended_WORK(extended_index-extended_ROW_LENGTH)
     &       + extended_WORK(extended_index))
          ENDDO
        ENDIF

        IF (at_base_of_LPG) THEN
! Loop over row above bottom row, missing halos
          DO I=END_P_POINT_NO_HALO-ROW_LENGTH+FIRST_ROW_PT,
     &         END_P_POINT_NO_HALO-ROW_LENGTH+LAST_ROW_PT
            extended_index=extended_address(I)
            V_TERM(I)=0.5*
     &        (extended_WORK(extended_index-extended_ROW_LENGTH)
     &       + extended_WORK(extended_index))
          ENDDO
        ENDIF


CL
CL---------------------------------------------------------------------
CL    SECTION 3.     CALCULATE VERTICAL FLUX AND COMBINE WITH U AND V
CL                   TERMS TO FORM INCREMENT.
CL---------------------------------------------------------------------

CL    VERTICAL FLUX ON INPUT IS .5*TIMESTEP*ETADOT*D(FIELD)/D(ETA)
CL    AT LEVEL K-1/2. AT THE END OF THIS SECTION IT IS THE SAME
CL    QUANTITY BUT AT LEVEL K+1/2.

! Loop over field, missing top and bottom rows and halos
      if(k.ne.1.and.k.ne.P_LEVELS)then

cdir$ unroll4
      DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
        SCALAR1 = .5 * ADVECTION_TIMESTEP *
     *         ETADOT(I,K+1) * (FIELD(I,K+1) - FIELD(I,K))
        SCALAR2 = WORK(I)
        FIELD_INC(I,K) = SCALAR1 +SCALAR2
      IF (LWHITBROM) FIELD_INC(I,K) = FIELD_INC(I,K)
     *                  + FIELD(I,K)*BRSP(I,K)
        WORK(I)=SCALAR1
      ENDDO
      else if(k.eq.1) then
cdir$ unroll4
      DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
        SCALAR1 = .5 * ADVECTION_TIMESTEP *
     *         ETADOT(I,K+1) * (FIELD(I,K+1) - FIELD(I,K))
        FIELD_INC(I,K) = SCALAR1
      IF (LWHITBROM) FIELD_INC(I,K) = FIELD_INC(I,K)
     *                  + FIELD(I,K)*BRSP(I,K)
      WORK(I)=SCALAR1
      END DO
      else if(k.eq.P_LEVELS) then
cdir$ unroll4
      DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
        SCALAR2 = WORK(I)
        FIELD_INC(I,K) =  SCALAR2
      IF (LWHITBROM) FIELD_INC(I,K) = FIELD_INC(I,K)
     *                  + FIELD(I,K)*BRSP(I,K)
      END DO
      endif ! if(k.ne.1.and.k.ne.P_LEVELS)then
cdir$ unroll4
      DO I=START_POINT_NO_HALO,END_P_POINT_NO_HALO
        FIELD_INC(I,K) = FIELD_INC(I,K) +
     *                   ADVECTION_TIMESTEP * SEC_P_LATITUDE(I) *
     *                  (U_TERM(I)+V_TERM(I))
      ENDDO



CL   LIMITED AREA MODEL SET BOUNDARY INCREMENTS
CL   TO ZERO.

       IF (at_left_of_LPG) THEN
          DO I=START_POINT_NO_HALO+FIRST_ROW_PT-1,
     &         END_P_POINT_NO_HALO,ROW_LENGTH
            FIELD_INC(I,K)=0.
          ENDDO
        ENDIF

        IF (at_right_of_LPG) THEN
          DO I=START_POINT_NO_HALO+LAST_ROW_PT-1,
     &         END_P_POINT_NO_HALO,ROW_LENGTH
            FIELD_INC(I,K)=0.
          ENDDO
        ENDIF

      ENDDO !DO K=1,P_LEVELS

      ENDIF !IF(L_SECOND)

CL    END OF ROUTINE ADV_P_GD

      RETURN
      END
