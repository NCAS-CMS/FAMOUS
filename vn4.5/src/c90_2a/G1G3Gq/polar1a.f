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
CLL  SUBROUTINE POLAR-------------------------------------------------
CLL
CLL  Purpose:  This routine updates the polar values along one level for
CLL            a primary model variable stored at p-points (ie pstar,
CLL            thetaL or qT).  This is done by adding in the
CLL            average meridional flux from the adjacent
CLL            equatorward row.
CLL
CLL M.Mawson    <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
!  4.1    28/06/95   Changed interface to allow multiple levels to
!                    be processed in single call + added MPP code.
!                                                     P.Burton
!  4.4    11/08/97   Added fast non-reproducible sums    P.Burton
!LL  4.5  27/04/98  Add Fujitsu vectorization directive.
!LL                                           RBarnes@ecmwf.int
CLL
CLL  Documentation:
CLL            Section 3.6 of Unified Model Documention Paper No 10.
CLL  -----------------------------------------------------------------
C
C*L  ARGUMENTS:-------------------------------------------------------
      SUBROUTINE POLAR(FIELD,NP_FLUX_FIELD,SP_FLUX_FIELD,
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
     &                 FIELD_SIZE,NP_FLUX_FIELD_SIZE,SP_FLUX_FIELD_SIZE,
     &                 NP_FLUX_START,SP_FLUX_START,ROW_LENGTH,
     &                 N_LEVELS)

      IMPLICIT NONE

      INTEGER
     &  FIELD_SIZE           ! IN size of single level of FIELD
     &, NP_FLUX_FIELD_SIZE   ! IN size of single level of NP_FLUX_FIELD
     &, SP_FLUX_FIELD_SIZE   ! IN size of single level of SP_FLUX_FIELD
     &, NP_FLUX_START        ! IN offset in NP_FLUX_FIELD of NP flux
     &, SP_FLUX_START        ! IN offset in SP_FLUX_FIELD of SP flux
     &, ROW_LENGTH           ! IN points per row
     &, N_LEVELS             ! IN number of levels to process

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
     &  FIELD(FIELD_SIZE,N_LEVELS)
     &                    ! INOUT primary field to be updated
     &, NP_FLUX_FIELD(NP_FLUX_FIELD_SIZE,N_LEVELS)
     &                    ! IN field containing fluxes for north pole
     &, SP_FLUX_FIELD(SP_FLUX_FIELD_SIZE,N_LEVELS)
     &                    ! IN field containing fluxes for south pole

      INTEGER info

! Local variables:
      INTEGER I,K

      REAL MEAN_NP(N_LEVELS),MEAN_SP(N_LEVELS)

      DO K=1,N_LEVELS
        MEAN_NP(K)=0.0
        MEAN_SP(K)=0.0
      ENDDO

      IF (at_top_of_LPG) THEN
        CALL GCG_RVECSUMR(NP_FLUX_FIELD_SIZE,ROW_LENGTH-2*EW_Halo,
     &                    NP_FLUX_START+EW_Halo,N_LEVELS,
     &                    NP_FLUX_FIELD,GC_ROW_GROUP,info,MEAN_NP)
        DO K=1,N_LEVELS
          MEAN_NP(K)=MEAN_NP(K)/GLOBAL_ROW_LENGTH
          DO I=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
            FIELD(I,K)=FIELD(I,K)+MEAN_NP(K)
          ENDDO
        ENDDO
      ENDIF

      IF (at_base_of_LPG) THEN
        CALL GCG_RVECSUMR(SP_FLUX_FIELD_SIZE,ROW_LENGTH-2*EW_Halo,
     &                    SP_FLUX_START+EW_Halo,N_LEVELS,
     &                    SP_FLUX_FIELD,GC_ROW_GROUP,info,MEAN_SP)
        DO K=1,N_LEVELS
          MEAN_SP(K)=MEAN_SP(K)/GLOBAL_ROW_LENGTH
          DO I=P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
            FIELD(I,K)=FIELD(I,K)+MEAN_SP(K)
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END
