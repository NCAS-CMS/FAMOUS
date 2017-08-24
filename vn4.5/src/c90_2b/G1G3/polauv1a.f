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
CLL  SUBROUTINE POLAR_UV----------------------------------------------
CLL
CLL  Purpose:
CLL            This routine updates the polar values of u and v
CLL            stored a half-grid length from the poles by vectorially
CLL            averaging the input fields from the adjacent equatorwards
CLL            row to obtain mean cartesian winds and also calculates a
CLL            mean u and v for this row.
CLL            The polar row is set to have the same mean vorticity and
CLL            divergence as the adjacent row by COS(LATITUDE) scaling.
CLL
CLL M.Mawson    <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
!     4.1     19/06/95  Rewritten to allow multiple levels to be
!                       processed and added MPP code.   P.Burton
!LL   4.4     12/08/97  Faster non-reproducible sums added.  P.Burton
CLL
CLL  System components covered: P196
CLL
CLL  Documentation:
CLL            Section 3.6 of Unified Model Documention Paper No 10.
CLL
CLL  -----------------------------------------------------------------
C
C*L  ARGUMENTS:-------------------------------------------------------
      SUBROUTINE POLAR_UV(U,V,ROW_LENGTH,U_POINTS,LEVELS,
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
     &                    COS_LAMBDA,SIN_LAMBDA)

      IMPLICIT NONE

      INTEGER
     &  ROW_LENGTH   ! IN  Number of points per row
     &, U_POINTS     ! IN  Horizontal size of fields on U grid
     &, LEVELS       ! IN  Number of levels to process

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
     &  U(U_POINTS,LEVELS)  ! INOUT  U field to process
     &, V(U_POINTS,LEVELS)  ! INOUT  V field to process
     &, SIN_LAMBDA(ROW_LENGTH)  !  IN Sine longitude
     &, COS_LAMBDA(ROW_LENGTH)  !  IN Cosine longitude

! Local arrays
! 4 arrays containing cartesian components of polar fluxes

      REAL
     &  CART_U_NP(ROW_LENGTH,LEVELS)
     &, CART_U_SP(ROW_LENGTH,LEVELS)
     &, CART_V_NP(ROW_LENGTH,LEVELS)
     &, CART_V_SP(ROW_LENGTH,LEVELS)

! 8 arrays containing the calculated means

      REAL
     &  MEAN_U_NP(LEVELS)
     &, MEAN_U_SP(LEVELS)
     &, MEAN_V_NP(LEVELS)
     &, MEAN_V_SP(LEVELS)
     &, MEAN_CARTESIAN_U_NP(LEVELS)
     &, MEAN_CARTESIAN_U_SP(LEVELS)
     &, MEAN_CARTESIAN_V_NP(LEVELS)
     &, MEAN_CARTESIAN_V_SP(LEVELS)

! Local variables
      INTEGER
     &  K,I,I_NP,I_SP ! loop counters and indexes
     &, NP_ADJACENT_ROW_START  ! start of row below NP row
     &, SP_ADJACENT_ROW_START  ! start of row above SP row
     &, LOCAL_ROW_PTS  ! number of non-halo points in row
     &, info  ! return code from GC routines

      REAL ROW_LENGTH_RECIP,ONE_THIRD

!-----------------------------------------------------------------------


      NP_ADJACENT_ROW_START=TOP_ROW_START+ROW_LENGTH
      SP_ADJACENT_ROW_START=U_BOT_ROW_START-ROW_LENGTH
      LOCAL_ROW_PTS=LAST_ROW_PT-FIRST_ROW_PT+1

      ROW_LENGTH_RECIP=1.0/GLOBAL_ROW_LENGTH
      ONE_THIRD=1.0/3.0

! 1. Resolve u and v vectorially onto cartesian grid

! North Pole
      IF (at_top_of_LPG) THEN
        DO K=1,LEVELS
          DO I=FIRST_ROW_PT,LAST_ROW_PT
            I_NP=I+NP_ADJACENT_ROW_START-1
!           I_NP index points to points along row beneath North Pole
            CART_U_NP(I,K)=
     &        U(I_NP,K)*COS_LAMBDA(I)-V(I_NP,K)*SIN_LAMBDA(I)
            CART_V_NP(I,K)=
     &        V(I_NP,K)*COS_LAMBDA(I)+U(I_NP,K)*SIN_LAMBDA(I)
          ENDDO
        ENDDO
      ENDIF

! South Pole
      IF (at_base_of_LPG) THEN
        DO K=1,LEVELS
          DO I=FIRST_ROW_PT,LAST_ROW_PT
            I_SP=I+SP_ADJACENT_ROW_START-1
!           I_SP index points to points along row above South Pole
            CART_U_SP(I,K)=
     &        U(I_SP,K)*COS_LAMBDA(I)+V(I_SP,K)*SIN_LAMBDA(I)
            CART_V_SP(I,K)=
     &        V(I_SP,K)*COS_LAMBDA(I)-U(I_SP,K)*SIN_LAMBDA(I)
          ENDDO
        ENDDO
      ENDIF

! 2. Compute mean cartesian values at poles, and mean u and v

      DO K=1,LEVELS
        MEAN_CARTESIAN_U_NP(K)=0.0
        MEAN_CARTESIAN_U_SP(K)=0.0
        MEAN_CARTESIAN_V_NP(K)=0.0
        MEAN_CARTESIAN_V_SP(K)=0.0
        MEAN_U_NP(K)=0.0
        MEAN_U_SP(K)=0.0
        MEAN_V_NP(K)=0.0
        MEAN_V_SP(K)=0.0

      ENDDO ! K : loop over levels

      IF (at_top_of_LPG) THEN
        CALL GCG_RVECSUMF(ROW_LENGTH,LOCAL_ROW_PTS,FIRST_ROW_PT,
     &                    LEVELS,CART_U_NP,GC_ROW_GROUP,
     &                    info,MEAN_CARTESIAN_U_NP)
        CALL GCG_RVECSUMF(ROW_LENGTH,LOCAL_ROW_PTS,FIRST_ROW_PT,
     &                    LEVELS,CART_V_NP,GC_ROW_GROUP,
     &                    info,MEAN_CARTESIAN_V_NP)

        CALL GCG_RVECSUMF(U_POINTS,LOCAL_ROW_PTS,
     &                    NP_ADJACENT_ROW_START+FIRST_ROW_PT-1,
     &                    LEVELS,U,GC_ROW_GROUP,info,MEAN_U_NP)
        CALL GCG_RVECSUMF(U_POINTS,LOCAL_ROW_PTS,
     &                    NP_ADJACENT_ROW_START+FIRST_ROW_PT-1,
     &                    LEVELS,V,GC_ROW_GROUP,info,MEAN_V_NP)
      ENDIF

      IF (at_base_of_LPG) THEN
        CALL GCG_RVECSUMF(ROW_LENGTH,LOCAL_ROW_PTS,FIRST_ROW_PT,
     &                    LEVELS,CART_U_SP,GC_ROW_GROUP,
     &                    info,MEAN_CARTESIAN_U_SP)
        CALL GCG_RVECSUMF(ROW_LENGTH,LOCAL_ROW_PTS,FIRST_ROW_PT,
     &                    LEVELS,CART_V_SP,GC_ROW_GROUP,
     &                    info,MEAN_CARTESIAN_V_SP)

        CALL GCG_RVECSUMF(U_POINTS,LOCAL_ROW_PTS,
     &                    SP_ADJACENT_ROW_START+FIRST_ROW_PT-1,
     &                    LEVELS,U,GC_ROW_GROUP,info,MEAN_U_SP)
        CALL GCG_RVECSUMF(U_POINTS,LOCAL_ROW_PTS,
     &                    SP_ADJACENT_ROW_START+FIRST_ROW_PT-1,
     &                    LEVELS,V,GC_ROW_GROUP,info,MEAN_V_SP)
      ENDIF

      DO K=1,LEVELS
        MEAN_CARTESIAN_U_NP(K)=
     &    MEAN_CARTESIAN_U_NP(K)*ROW_LENGTH_RECIP
        MEAN_CARTESIAN_U_SP(K)=
     &    MEAN_CARTESIAN_U_SP(K)*ROW_LENGTH_RECIP
        MEAN_CARTESIAN_V_NP(K)=
     &    MEAN_CARTESIAN_V_NP(K)*ROW_LENGTH_RECIP
        MEAN_CARTESIAN_V_SP(K)=
     &    MEAN_CARTESIAN_V_SP(K)*ROW_LENGTH_RECIP
        MEAN_U_NP(K) = MEAN_U_NP(K)*ROW_LENGTH_RECIP
        MEAN_U_SP(K) = MEAN_U_SP(K)*ROW_LENGTH_RECIP
        MEAN_V_NP(K) = MEAN_V_NP(K)*ROW_LENGTH_RECIP
        MEAN_V_SP(K) = MEAN_V_SP(K)*ROW_LENGTH_RECIP
      ENDDO

! 3. Resolve mean values back to lat-lon grid and add in fluxes
!    Scale MEAN_U and MEAN_V by 1/3 to give uniform vorticity and
!    divergence.

! North Pole
      IF (at_top_of_LPG) THEN
        DO K=1,LEVELS
          DO I=FIRST_ROW_PT,LAST_ROW_PT

            I_NP=TOP_ROW_START+I-1
!           This points to the real North Pole row.

            U(I_NP,K) = MEAN_CARTESIAN_U_NP(K)*COS_LAMBDA(I)+
     &                  MEAN_CARTESIAN_V_NP(K)*SIN_LAMBDA(I)+
     &                  MEAN_U_NP(K)*ONE_THIRD
            V(I_NP,K) = MEAN_CARTESIAN_V_NP(K)*COS_LAMBDA(I)-
     &                  MEAN_CARTESIAN_U_NP(K)*SIN_LAMBDA(I)+
     &                  MEAN_V_NP(K)*ONE_THIRD
          ENDDO
        ENDDO
      ENDIF

! South Pole
      IF (at_base_of_LPG) THEN
        DO K=1,LEVELS
          DO I=FIRST_ROW_PT,LAST_ROW_PT

            I_SP=U_BOT_ROW_START+I-1
!           This points to the real South Pole row.

            U(I_SP,K) = MEAN_U_SP(K)*ONE_THIRD +
     &                  MEAN_CARTESIAN_U_SP(K)*COS_LAMBDA(I)-
     &                  MEAN_CARTESIAN_V_SP(K)*SIN_LAMBDA(I)
            V(I_SP,K) = MEAN_V_SP(K)*ONE_THIRD +
     &                  MEAN_CARTESIAN_V_SP(K)*COS_LAMBDA(I)+
     &                  MEAN_CARTESIAN_U_SP(K)*SIN_LAMBDA(I)
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END
