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
CLL   SUBROUTINE UV_DIF -----------------------------------------
CLL
CLL   PURPOSE:  CALCULATES DIFFUSIVE INCREMENTS TO U AND V AT
CLL             ONE MODEL LEVEL USING EQUATION (47) AND ADDS TO FIELD.
CLL             IF GLOBAL MODEL RUN THEN UPDATES POLAR VALUES.
CLL              PERFORMS FULL DEL-SQUARED CALCULATION WITH
CLL              ACROSS-POLE DIFFERENCING IN GLOBAL CASE (TCJ).
CLL              IF STEEP SLOPE THEN EFFECTIVE DIFFUSION IS ZERO. (TD)
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL   VERSION FOR CRAY Y-MP
CLL
CLL   WRITTEN BY M.H MAWSON.
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  4.0  03/02/95  RE-WRITTEN TO MAKE MORE EFFICIENT WITH TESTING FOR
CLL                 STEEP SLOPES. AUTHOR: T.DAVIES. REVIEWER: M.MAWSON
!     3.5    28/03/95 MPP code additions  P.Burton
!     4.1    07/05/96 Added MPP code and TYPFLDPT arguments  P.Burton
!     4.3    22/04/97 Added optimised vector shift 
!                     B. Carruthers
CLL
CLL   DOCUMENTATION:       THE EQUATION USED IS (47)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER
CLL                        NO. 10 M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLL                        VERSION 22 DATED 01/06/92.
CLL
CLL   SYSTEM COMPONENTS COVERED: P132
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION:       THE EQUATION USED IS (47)
CLL                        IN UNIFIED MODEL DOCUMENTATION PAPER
CLL                        NO. 10 M.J.P. CULLEN,T.DAVIES AND M.H.MAWSON
CLL                        VERSION 16 DATED 09/01/91.
CLLEND-------------------------------------------------------------

C*L   ARGUMENTS:---------------------------------------------------
      SUBROUTINE UV_DIF
     1             (FIELDU,FIELDV,RS_SQUARED_DELTAP,
     2              SEC_U_LATITUDE,START_U_UPDATE,END_U_UPDATE,
     3              ROW_LENGTH,
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
     &              P_FIELD,U_FIELD,
     4              DIFFUSION_EW,DIFFUSION_NS)

      IMPLICIT NONE

      INTEGER
     *  U_FIELD            !IN DIMENSION OF FIELDS ON VELOCITY GRID
     *, P_FIELD            !IN DIMENSION OF FIELDS ON PRESSURE GRID
     *, ROW_LENGTH         !IN NUMBER OF POINTS PER ROW
     *, START_U_UPDATE     !IN FIRST POINT TO BE UPDATED.
     *, END_U_UPDATE       !IN LAST POINT TO BE UPDATED.

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
     * FIELDU(U_FIELD)                  !INOUT. U DIFFUSION FIELD.
     *,FIELDV(U_FIELD)                  !INOUT. V DIFFUSION FIELD.

      REAL
     * RS_SQUARED_DELTAP(U_FIELD)      !IN HOLDS RS*RS*DELTA P
     *,DIFFUSION_EW(P_FIELD)  !IN EFFECTIVE EW DIFFUSION COEFFICIENT
     *,DIFFUSION_NS(P_FIELD)  !IN EFFECTIVE NS DIFFUSION COEFFICIENT
     *,SEC_U_LATITUDE(U_FIELD)         !IN 1/COS(LAT) AT U POINTS


C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: 4 ARE REQUIRED

      REAL
     * FIELD1(P_FIELD)        ! GENERAL WORKSPACE
     *,FIELD2(P_FIELD)        ! GENERAL WORKSPACE
     *,FIELD3(P_FIELD)        ! GENERAL WORKSPACE
     *,FIELD4(P_FIELD)        ! GENERAL WORKSPACE
     &,U_COPY(ROW_LENGTH)     ! copy of polar row
     &,V_COPY(ROW_LENGTH)     ! copy of polar row
C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES

C LOCAL REALS.
      REAL
     *  SCALAR
C COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     *  I,IJ,J,JI,HALF_RL
      INTEGER info

C*L   EXTERNAL SUBROUTINE CALLS: NONE---------------------------------
C*---------------------------------------------------------------------
CL    MAXIMUM VECTOR LENGTH ASSUMED IS END_U_UPDATE-START_U_UPDATE+1+
CL                                   ROW_LENGTH
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL    SECTION 1.     CALCULATE FIRST TERM IN EQUATION (47)
CL---------------------------------------------------------------------

C----------------------------------------------------------------------
CL    SECTION 1.1    CALCULATE DELTALAMBDA TERMS
C----------------------------------------------------------------------


      DO I=START_U_UPDATE+1,END_U_UPDATE
       FIELD1(I) =FIELDU(I)-FIELDU(I-1)
       FIELD2(I) =FIELDV(I)-FIELDV(I-1)
      END DO

      FIELD1(START_U_UPDATE)=FIELD1(START_U_UPDATE+1)
      FIELD2(START_U_UPDATE)=FIELD2(START_U_UPDATE+1)

C----------------------------------------------------------------------
CL    SECTION 1.3   COMPLETE DELTALAMBDA TERM
C----------------------------------------------------------------------

      DO I= START_U_UPDATE+1,END_U_UPDATE
       FIELD3(I-1)=(DIFFUSION_EW(I)*FIELD1(I)-
     &              DIFFUSION_EW(I-1)*FIELD1(I-1))*
     &              SEC_U_LATITUDE(I-1)
       FIELD4(I-1)=(DIFFUSION_EW(I)*FIELD2(I)-
     &              DIFFUSION_EW(I-1)*FIELD2(I-1))*
     &              SEC_U_LATITUDE(I-1)
      END DO

C  CORRECT END POINT
      FIELD3(END_U_UPDATE)=FIELD3(END_U_UPDATE-1)
      FIELD4(END_U_UPDATE)=FIELD4(END_U_UPDATE-1)

CL
CL---------------------------------------------------------------------
CL    SECTION 2.     CALCULATE PHI DIRECTION TERM AND ADD
CL                   ONTO FIRST TO GET TOTAL INCREMENT.
CL---------------------------------------------------------------------

! Loop over field, missing top row
      DO I=START_POINT_NO_HALO,LAST_U_VALID_PT
        FIELD1(I)=FIELDU(I-ROW_LENGTH)-FIELDU(I)
        FIELD2(I)=FIELDV(I-ROW_LENGTH)-FIELDV(I)
      END DO

C CALCULATE POLAR TERMS USING ACROSS-POLE DIFFERENCE, REMEMBERING SIGN
C CHANGE ACROSS THE POLE
C NB: EFFECTIVE COS_P_LATITUDE IS 1/4 THAT AT ADJACENT ROW

      HALF_RL = global_ROW_LENGTH/2

      IF (at_top_of_LPG) THEN
! North Pole
! Copy NP row into U/V_COPY arrays
        DO I=1,ROW_LENGTH
          J=TOP_ROW_START+I-1 ! point along North Pole row
          U_COPY(I)=FIELDU(J)
          V_COPY(I)=FIELDV(J)
        ENDDO

! and rotate these rows by half a global row length, so that item I now
! contains the value which was on the opposite side of the pole.

        CALL GCG_RVECSHIFT(ROW_LENGTH,ROW_LENGTH-2*EW_Halo,
     &                     FIRST_ROW_PT,1,HALF_RL,.TRUE.,U_COPY,
     &                     GC_ROW_GROUP,info)
        CALL GCG_RVECSHIFT(ROW_LENGTH,ROW_LENGTH-2*EW_Halo,
     &                     FIRST_ROW_PT,1,HALF_RL,.TRUE.,V_COPY,
     &                     GC_ROW_GROUP,info)

        DO I=1,ROW_LENGTH
          J=TOP_ROW_START+I-1 ! point along North Pole row
          FIELD1(J)=-(FIELDU(J)+U_COPY(I))
          FIELD2(J)=-(FIELDV(J)+V_COPY(I))
        ENDDO
      ENDIF ! (IF at_top_of_LPG)

      IF (at_base_of_LPG) THEN
! South Pole
! Copy SP row in U/V_COPY arrays
        DO I=1,ROW_LENGTH
          J=U_BOT_ROW_START+I-1 ! point along South Pole row
          U_COPY(I)=FIELDU(J)
          V_COPY(I)=FIELDV(J)
        ENDDO

! and rotate these rows by half a global row length, so that item I now
! contains the value which was on the opposite side of the pole.

        CALL GCG_RVECSHIFT(ROW_LENGTH,ROW_LENGTH-2*EW_Halo,
     &                     FIRST_ROW_PT,1,HALF_RL,.TRUE.,U_COPY,
     &                     GC_ROW_GROUP,info)
        CALL GCG_RVECSHIFT(ROW_LENGTH,ROW_LENGTH-2*EW_Halo,
     &                     FIRST_ROW_PT,1,HALF_RL,.TRUE.,V_COPY,
     &                     GC_ROW_GROUP,info)

        DO I=1,ROW_LENGTH
          J=P_BOT_ROW_START+I-1
          FIELD1(J)=U_COPY(I)+FIELDU(J-ROW_LENGTH)
          FIELD2(J)=V_COPY(I)+FIELDV(J-ROW_LENGTH)
        ENDDO
      ENDIF ! (IF at_base_of_LPG)


C----------------------------------------------------------------------
CL    SECTION 2.3    CALCULATE SECOND TERM IN EQUATION (47) AND ADD
CL                   ONTO FIRST TERM TO GET TOTAL CORRECTION.
C----------------------------------------------------------------------

      DO I= START_U_UPDATE,END_U_UPDATE
        SCALAR=SEC_U_LATITUDE(I)/RS_SQUARED_DELTAP(I)
        FIELDU(I)=(FIELD3(I)+DIFFUSION_NS(I)*FIELD1(I)
     &           -DIFFUSION_NS(I+ROW_LENGTH)*FIELD1(I+ROW_LENGTH))*
     &            SCALAR
        FIELDV(I)=(FIELD4(I)+DIFFUSION_NS(I)*FIELD2(I)
     &           -DIFFUSION_NS(I+ROW_LENGTH)*FIELD2(I+ROW_LENGTH))*
     &            SCALAR
       END DO

CL  LIMITED AREA ZERO AT LATERAL BOUNDARIES

CL    END OF ROUTINE UV_DIF

      RETURN
      END
