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
CLL   SUBROUTINE OMEGA_DIAG -----------------------------------------
CLL
CLL   PURPOSE:  CALCULATES OMEGA AT VELOCITY POINTS.
CLL   NOT SUITABLE FOR SINGLE COLUMN USE.
CLL   VERSION FOR CRAY Y-MP
CLL
CLL    Written by M.H. Mawson
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.1 02/02/93 CORRECTION OF SIGN ERROR IN DP_BY_DT TERM
CLL
CLL   4.0 12/07/95 CORRECTION OF FLOATING POINT IN ETADOT FOR
CLL                DEF GLOBAL NOT TRUE . M. Noguer
CLL   4.3 14/03/97 MPP Changes. S.D.Mullerworth
CLL   Programming Standard: Unified Model Documentation Paper No. 4,
CLL                         Standard B.
CLL
CLL   System components covered:
CLL
CLL   System task:
CLL
CLL   Documentation: The equations used are (29),(30) and (44)
CLL in Unified Model Documentation Paper No. 10 M.J.P. Cullen,T.Davies
CLL and M.H. Mawson.  Except that the variable radius of the earth RS is
CLL replaced by the constant value A.
CLL
CLLEND-------------------------------------------------------------

C*L   ARGUMENTS:---------------------------------------------------

      SUBROUTINE OMEGA_DIAG
     1                     (
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
     &                      U,V,OMEGA,SEC_P_LATITUDE,COS_U_LATITUDE,
     2                      PSTAR,PSTAR_OLD,DELTA_AK,DELTA_BK,
     3                      AK,BK,AKH,BKH,U_FIELD,P_FIELD,P_LEVELS,
     4                      ROW_LENGTH,LATITUDE_STEP_INVERSE,
     5                      LONGITUDE_STEP_INVERSE,ADVECTION_TIMESTEP)

      IMPLICIT NONE

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
     &  ROW_LENGTH            !IN NUMBER OF POINTS PER ROW
     &, P_LEVELS              !IN NUMBER OF PRESSURE LEVELS OF DATA
     &, P_FIELD               !IN NUMBER OF POINTS IN PRESSURE FIELD
     &, U_FIELD               !IN NUMBER OF POINTS IN VELOCITY FIELD

      REAL
     & U(U_FIELD,P_LEVELS)    !IN U VELOCITY
     &,V(U_FIELD,P_LEVELS)    !IN V VELOCITY
     &,PSTAR(P_FIELD)         !IN SURFACE PRESSURE AT P POINTS.
     &,PSTAR_OLD(P_FIELD)     !IN SURFACE PRESSURE AT PREVIOUS TIMESTEP
     &,SEC_P_LATITUDE(P_FIELD)!IN 1/COS(LAT) AT P POINTS
     &,COS_U_LATITUDE(U_FIELD)!IN COS(LAT) AT U POINTS
     &,LONGITUDE_STEP_INVERSE !IN 1/LONGITUDE INCREMENT
     &,LATITUDE_STEP_INVERSE  !IN 1/LATITUDE INCREMENT
     &,BK(P_LEVELS)           !IN HOLDS COEFFICIENT WHICH
     &                        !   MULTIPLIES PSTAR IN HYBRID CO-ORDS
     &                        !   AT LEVELS K
     &,AK(P_LEVELS)           !IN HOLDS FIRST COEFFICIENT
     &                        !   IN HYBRID CO-ORDS AT LEVELS K
     &,BKH(P_LEVELS+1)        !IN HOLDS COEFFICIENT WHICH
     &                        !   MULTIPLIES PSTAR IN HYBRID CO-ORDS
     &                        !   AT half LEVELS K
     &,AKH(P_LEVELS+1)        !IN HOLDS FIRST COEFFICIENT
     &                        !   IN HYBRID CO-ORDS AT  half LEVELS K
     &,DELTA_AK(P_LEVELS)     !IN AK(K+1/2)-AK(K-1/2)
     &,DELTA_BK(P_LEVELS)     !IN BK(K+1/2)-BK(K-1/2)
     &,ADVECTION_TIMESTEP     !IN


      REAL
     & OMEGA(U_FIELD,P_LEVELS) !OUT. OMEGA AT U POINTS.
C*---------------------------------------------------------------------

C*L   9 LOCAL ARRAYS NEEDED. -----------------------------------------

      REAL
     &  ETADOT(P_FIELD,P_LEVELS+1) ! HOLDS ETADOT OR PARTS THEREOF.
     &, WORK1(P_FIELD)             ! \
     &, WORK2(P_FIELD)             !  > GENERAL WORK ARRAYS.
     &, WORK3(P_FIELD)             ! /
     &, PSTAR_UV(U_FIELD)          ! PSTAR AT U POINTS.
     &, DELTA_P(U_FIELD)           ! DEPTH OF LAYER AT U POINTS.
     &, WP(U_FIELD)                ! TERM IN OMEGA EQUATION.
     &, DELTA_AKH(P_LEVELS+1)      ! AK(K) -AK(K-1)
     &, DELTA_BKH(P_LEVELS+1)      ! BK(K) - BK(K-1)
C*---------------------------------------------------------------------

C DEFINE COUNT VARIABLES FOR DO LOOPS ETC.
      INTEGER
     &  I,J,K
     &, START,END
C DEFINE LOCAL SCALARS
      REAL
     &  SCALAR


      REAL 
     &  SUM_N_ARRAY(ROW_LENGTH,P_LEVELS)
     &  ,SUM_S_ARRAY(ROW_LENGTH,P_LEVELS)
      
      INTEGER 
     & info

C*L   EXTERNAL SUBROUTINE CALLS:---------------------------------------
      EXTERNAL P_TO_UV
C*---------------------------------------------------------------------
C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
C*----------------------------------------------------------------------


CL    MAXIMUM VECTOR LENGTH ASSUMED IS P_FIELD
CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL    SECTION 0. SET UP PSTAR AT U POINTS AND DELTA_AKH,_BKH
CL---------------------------------------------------------------------
CL

CL CALL P_TO_UV TO INTERPOLATE PSTAR TO U POINTS.

      CALL P_TO_UV(PSTAR,PSTAR_UV,P_FIELD,U_FIELD,ROW_LENGTH,
     &             P_FIELD/ROW_LENGTH)
      CALL SWAPBOUNDS(PSTAR_UV,ROW_LENGTH,tot_U_ROWS,
     &  EW_Halo,NS_Halo,1)

C SET UP DELTA_AKH,_BKH
      DO K=2,P_LEVELS
        DELTA_AKH(K) = AK(K) - AK(K-1)
        DELTA_BKH(K) = BK(K) - BK(K-1)
      END DO
      DELTA_AKH(1) = 0.
      DELTA_AKH(P_LEVELS+1) = 0.
      DELTA_BKH(1) = 0.
      DELTA_BKH(P_LEVELS+1) = 0.

C LOOP OVER LEVELS
      DO 100 K=1,P_LEVELS

CL---------------------------------------------------------------------
CL    SECTION 1. CALCULATE DIVERGENCE AS IN EQUATION (30).
CL---------------------------------------------------------------------

        DO I=1,U_FIELD
          DELTA_P(I) = DELTA_AK(K)+DELTA_BK(K)*PSTAR_UV(I)
        END DO

C CALCULATE DU/D(LAMDA), STORE IN ETADOT.
        DO I=2,U_FIELD
          WORK1(I) = LONGITUDE_STEP_INVERSE*(U(I,K)*DELTA_P(I)
     &                  -U(I-1,K)*DELTA_P(I-1))
        END DO

C CALCULATE DV/D(PHI), store in work2.
        DO I=ROW_LENGTH+1,U_FIELD
          WORK2(I) = LATITUDE_STEP_INVERSE*((V(I-ROW_LENGTH,K)
     &             *COS_U_LATITUDE(I-ROW_LENGTH))*DELTA_P(I-ROW_LENGTH)
     &             -(V(I,K)*COS_U_LATITUDE(I))*DELTA_P(I))
        END DO

        WORK2(1) = 0.

C CALCULATE DIVERGENCES. STORE IN ETADOT

        DO I=START_POINT_NO_HALO+EW_Halo,END_P_POINT_NO_HALO-EW_Halo
          ETADOT(I,K)= (SEC_P_LATITUDE(I)*.5)*(WORK1(I)
     &                         + WORK1(I-ROW_LENGTH)
     &                         + WORK2(I) + WORK2(I-1))
        END DO

C ZERO DIVERGENCES ON BOUNDARIES.
        IF (at_left_of_LPG) THEN
          DO I=ROW_LENGTH+1+EW_Halo,P_FIELD,ROW_LENGTH
            ETADOT(I,K) = 0.
          ENDDO
        ENDIF 
        IF (at_right_of_LPG) THEN
          DO I=2*ROW_LENGTH-EW_Halo,P_FIELD,ROW_LENGTH
            ETADOT(I,K) = 0.
          ENDDO
        ENDIF


 100  CONTINUE


CL
CL---------------------------------------------------------------------
CL    SECTION 2. CALCULATE VERTICAL VELOCITY. EQUATION (29).
CL---------------------------------------------------------------------

      START = START_POINT_NO_HALO+EW_Halo
      END = END_P_POINT_NO_HALO-EW_Halo

C ---------------------------------------------------------------------
CL    SECTION 2.1 SUM DIVERGENCES THROUGHOUT ATMOSPHERE.
C ---------------------------------------------------------------------

C BY CODING THE SUMMATION AS FOLLOWS THE VALUES PUT INTO EACH LEVEL
C OF ETADOT ARE THE ONES NEEDED FOR THE SECOND SUMMATION TERM
C IN EQUATION 29, WHILE THE TOTAL SUM IS HELD IN ETADOT( ,1)

      DO 210 K=P_LEVELS-1,1,-1
        DO I=START,END
          ETADOT(I,K)= ETADOT(I,K)+ETADOT(I,K+1)
        END DO
 210  CONTINUE

C ---------------------------------------------------------------------
CL    SECTION 2.2 CALCULATE MASS-WEIGHTED VERTICAL VELOCITY.
C ---------------------------------------------------------------------

C DP/D(PSTAR) IS NOTHING MORE THAN THE BK COEFFICENT.

      DO 220 K= P_LEVELS,2,-1
        DO I= START, END
          ETADOT(I,K)= ETADOT(I,K) - BKH(K) * ETADOT(I,1)
        END DO
 220  CONTINUE

      DO 230 K=1,P_LEVELS
C SET NORTHERN AND SOUTHERN BOUNDARIES TO ZERO.
        IF (at_top_of_LPG) THEN
          DO I=TOP_ROW_START+EW_Halo,START-1
            ETADOT(I,K) = 0.
          ENDDO
        ENDIF
        IF (at_base_of_LPG) THEN
          DO I=END+1,LAST_P_FLD_PT-EW_Halo
            ETADOT(I,K) = 0.
          ENDDO
        ENDIF
 230  CONTINUE

CL INTERPOLATE ETADOT TO U POINTS.
CL LOOP OVER LEVELS
      CALL SWAPBOUNDS(ETADOT(1,2),ROW_LENGTH,tot_P_ROWS,EW_Halo,
     &  NS_Halo,P_LEVELS-1)
      DO K=2,P_LEVELS
        CALL P_TO_UV(ETADOT(1,K),WORK1,P_FIELD,U_FIELD,ROW_LENGTH,
     &                 P_FIELD/ROW_LENGTH)
        DO I=1,U_FIELD
          ETADOT(I,K) = WORK1(I)
        END DO
      END DO

CL SET ETADOT AT TOP AND BOTTOM TO ZERO.
      DO I=1,P_FIELD
        ETADOT(I,1) = 0.
        ETADOT(I,P_LEVELS+1) = 0.
      END DO

CL LOOP OVER LEVELS
      DO 300 K=1,P_LEVELS
CL
CL---------------------------------------------------------------------
CL    SECTION 3.     CALCULATE DP/DT
CL---------------------------------------------------------------------

        IF(BK(K).EQ.0.) THEN
C A CONSTANT PRESSURE LEVEL SO DP/DT IS ZERO.
          DO I=1,U_FIELD
            WORK3(I) = 0.
          END DO
        ELSE
C CALCULATE DP/DT.
          SCALAR = BK(k)/ADVECTION_TIMESTEP
          DO I=1,P_FIELD
            WORK2(I) = (PSTAR(I)-PSTAR_OLD(I))*SCALAR
          END DO
C INTERPOLATE DP/DT TO U POINTS.
          CALL P_TO_UV(WORK2,WORK3,P_FIELD,U_FIELD,ROW_LENGTH,
     &                 P_FIELD/ROW_LENGTH)
        END IF

CL---------------------------------------------------------------------
CL    SECTION 4.     CALCULATE U.GRAD P
CL---------------------------------------------------------------------

C----------------------------------------------------------------------
CL    SECTION 4.1    CALCULATE U DP/D(LAMDA)
C----------------------------------------------------------------------

C CALCULATE U DP/D(LAMDA) BETWEEN P POINTS
        SCALAR=(.5*LONGITUDE_STEP_INVERSE)*(BK(K)/A)
        DO I=START_POINT_NO_HALO+EW_Halo,END_P_POINT_NO_Halo-EW_Halo
          WORK1(I) = ((U(I,K)+U(I-ROW_LENGTH,K))*
     &      (PSTAR(I+1)-PSTAR(I)))*(SEC_P_LATITUDE(I)*
     &      SCALAR)
        END DO
C SET POLAR/BOUNDARY VALUES TO ZERO.
        IF (at_top_of_LPG) THEN
          DO I=TOP_ROW_START+EW_Halo,START_POINT_NO_HALO-1-EW_Halo
            WORK1(I) = 0.
          ENDDO
        ENDIF
        IF (at_base_of_LPG) THEN
          DO I=P_BOT_ROW_START+EW_Halo,LAST_P_FLD_PT-EW_Halo
            WORK1(I) = 0.
          ENDDO
        ELSE
C Calculate bottom row halo values to allow WP to be calculated
          DO I=P_BOT_ROW_START+EW_Halo,LAST_P_VALID_PT-EW_Halo
            WORK1(I) = ((U(I,K)+U(I-ROW_LENGTH,K))*
     &        (PSTAR(I+1)-PSTAR(I)))*(SEC_P_LATITUDE(I)
     &        *SCALAR)
          ENDDO
        ENDIF

C CALCULATE U DP/D(LONGITUDE) AT U POINTS
        IF (at_top_of_LPG) THEN
          DO I=TOP_ROW_START+EW_Halo,START_POINT_NO_HALO-1+EW_Halo
            WP(I) = .5*WORK1(I+ROW_LENGTH)
          ENDDO
        ENDIF
        DO I=START_POINT_NO_HALO+EW_Halo,LAST_U_FLD_PT-EW_Halo
          WP(I) = .5*(WORK1(I)+WORK1(I+ROW_LENGTH))
        END DO

C----------------------------------------------------------------------
CL    SECTION 4.2    CALCULATE V DP/D(PHI) AND HENCE U.GRAD P
C----------------------------------------------------------------------

C CALCULATE V DP/D(PHI) BETWEEN P POINTS.
        SCALAR=(.5*LATITUDE_STEP_INVERSE)*(BK(K)/A)
        DO I=TOP_ROW_START+1,LAST_U_FLD_PT
          WORK2(I) = ((V(I,K)+V(I-1,K))*
     &                (PSTAR(I)-PSTAR(I+ROW_LENGTH)))
     &             *SCALAR
        END DO

CL CALCULATE V DP/D(LAT) AT U POINTS AND ADD TO WP.
        DO I=TOP_ROW_START+1,LAST_U_FLD_PT-EW_Halo
          WORK1(I) = .5*(WORK2(I)+WORK2(I+1))
        END DO


        DO I=FIRST_FLD_PT+EW_Halo,LAST_U_FLD_PT-EW_Halo
          WP(I) = WP(I) +WORK1(I)
        END DO

CL---------------------------------------------------------------------
CL    SECTION 5.     CALCULATE OMEGA AS IN EQUATION (44).
CL---------------------------------------------------------------------

        DO I=FIRST_FLD_PT+EW_Halo,LAST_U_FLD_PT-EW_Halo
           OMEGA(I,K)= (WP(I)+WORK3(I))+((.5*(ETADOT(I,K+1)*
     &                  (DELTA_AKH(K+1)+DELTA_BKH(K+1)*PSTAR_UV(I))
     &                  +ETADOT(I,K)
     &                  *(DELTA_AKH(K)+DELTA_BKH(K)*PSTAR_UV(I))))
     &                  /(A*(DELTA_AK(K)+DELTA_BK(K)*PSTAR_UV(I))))
        END DO
C Initialise unused rows of OMEGA
        DO I=FIRST_FLD_PT+EW_Halo-1,1,-1
          OMEGA(I,K)=OMEGA(I+ROW_LENGTH,K)
        ENDDO
        DO I=LAST_U_FLD_PT-EW_Halo+1,U_FIELD
          OMEGA(I,K)=OMEGA(I-ROW_LENGTH,K)
        ENDDO

CL    LIMITED AREA MODEL SET VERTICAL VELOCITY ON BOUNDARY TO ZERO.
        IF (at_left_of_LPG) THEN
          DO I=ROW_LENGTH+1+EW_Halo,LAST_U_FLD_PT,ROW_LENGTH
            OMEGA(I,K) = 0.
          ENDDO
        ENDIF
        IF (at_right_of_LPG) THEN
          DO I=ROW_LENGTH-EW_Halo,LAST_U_FLD_PT,ROW_LENGTH
            OMEGA(I,K) = 0.
            OMEGA(I-1,K) = 0.
          ENDDO
        ENDIF

CL END LOOP OVER LEVELS.
 300  CONTINUE
      CALL SWAPBOUNDS(OMEGA(1,1),ROW_LENGTH,tot_U_ROWS,EW_Halo,
     &  NS_Halo,P_LEVELS)

CL    END OF ROUTINE OMEGA_DIAG

      RETURN
      END
