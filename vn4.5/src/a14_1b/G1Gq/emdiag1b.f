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
!+ Modified version of EMDIAG1A with fewer global sums
!
! Subroutine Interface:
      SUBROUTINE ENG_MASS_DIAG (TL,U,V,AREA_P,AREA_UV,P_FIELD,
     &                          U_FIELD,ROW_LENGTH,ROWS,
     &                          DELTA_AK,DELTA_BK,AK,BK,TOT_ENERGY,
     &                          TOT_MASS_P,PART_MASS_P,P_LEVELS,PSTAR,
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
     &                          LLINTS,LWHITBROM)
      IMPLICIT NONE
!
! Description:
! Part of the energy correction suite of routines:
! To globally intergrate total energy and mass of the atmosphere
! This version reduces the total number of global sums done,
! which makes it much more efficient on MPP machines (where the
! global sum can be an expensive operation).
! Results will not bit-compare with version EMDIAG1A because of
! different rounding errors.
!
! Method:
! There are three basic sums to be calculated:
! 1) Energy
! 2) Total Mass
! 3) Partial Mass
! Sums are calculated over levels for each column (each column
! is independent) and then the three sums are calculated.
!
! Current code owner : Paul Burton
!
! History
!  Model    Date      Modification history from model version 4.1
!  version
!    4.1    7/11/95   New DECK created to make EMDIAG suitable for
!                     MPP use. P.Burton
!    4.3   18/03/97   Corrected call to CALC_RS   P.Burton
!
! Subroutine Arguments:

      LOGICAL LLINTS,
     &        LWHITBROM    ! IN use White+Bromley terms?

      INTEGER P_FIELD,    ! IN vector length of variables on P grid
     &        U_FIELD,    ! IN vector length of variables on UV grid
     &        ROW_LENGTH, ! IN number of pointer per row
     &        ROWS,       ! IN number of rows on P grid
     &        P_LEVELS    ! IN number of levels in vertical

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

      REAL    TL(P_FIELD,P_LEVELS),    ! IN temperature
     &        U(U_FIELD,P_LEVELS),     ! IN U component of wind
     &        V(U_FIELD,P_LEVELS),     ! IN V component of wind
     &        AREA_P(P_FIELD),         ! IN area of cells in P grid
     &        AREA_UV(U_FIELD),        ! IN area of cells in UV grid
     &        DELTA_AK(P_LEVELS),      ! IN  \thickness of layers
     &        DELTA_BK(P_LEVELS),      ! IN  /in eta co-ordinates
     &        AK(P_LEVELS),            ! IN  \eta co-ordinates of
     &        BK(P_LEVELS),            ! IN  /mid-layer points
     &        PSTAR(P_FIELD)           ! IN pressure at surface

      REAL    TOT_ENERGY,   ! OUT total energy of atmosphere
     &        TOT_MASS_P,   ! OUT total mass of atmosphere
     &        PART_MASS_P   ! OUT partial mass of atmosphere


! Parameters
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

C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
C*----------------------------------------------------------------------

      INTEGER N_SUMS
      PARAMETER(N_SUMS=3)  ! there are 3 sums to do
!   Define magic numbers to define the various sums
      INTEGER energy,total_mass,partial_mass
      PARAMETER( energy=1,total_mass=2,partial_mass=3)

! Local variables

      REAL PSTAR_DELBK(P_FIELD),       ! pressure at surface*DELTA_BK
     &     DELP_P(P_FIELD),            ! mass elements on P grid
     &     DELP_UV(U_FIELD),           ! mass elements on UV grid
     &     RS_P_K(P_FIELD),            ! radii on P grid
     &     RS_UV_K(U_FIELD),           ! radii on U grid
     &     WORK(P_FIELD),              ! workspace
     &     SUM_ARRAY(P_FIELD,N_SUMS),  ! the array to be summed
     &     SUM_RESULTS(N_SUMS),        ! the sums of SUM_ARRAY
     &     TS(P_FIELD)                 ! output from CALC_RS

      INTEGER START_POINT,   ! point to start sums at
     &        END_P_POINT,   ! number of points to sum on P grid
     &        END_U_POINT    ! number of points to sum on U grid

      INTEGER I,K  ! loop variables

! 1.0 Set up the range of points to sum over

! Sum over all points - missing out halos and northern/southern
! boundaries/poles
      START_POINT=FIRST_FLD_PT
      END_P_POINT=LAST_P_FLD_PT
      END_U_POINT=LAST_U_FLD_PT

! QAN fix
! Zero DELP_P and RS_P_Karray
      DO I=1,P_FIELD
        DELP_P(I)=0.0
        RS_P_K(I)=0.0
      ENDDO


      DO K=1,N_SUMS
        DO I=1,P_FIELD
          SUM_ARRAY(I,K)=0.0
        ENDDO
        SUM_RESULTS(K)=0.0
      ENDDO

! 2.0 Now the loop over levels.
!     Sum the values up over each column and store in SUM_ARRAY

      DO K=1,P_LEVELS

! 2.1 Set up arrays required for this level

        DO I=FIRST_VALID_PT,LAST_P_VALID_PT
          PSTAR_DELBK(I)=-DELTA_BK(K)*PSTAR(I)
          DELP_P(I)=-DELTA_AK(K)+PSTAR_DELBK(I)
        ENDDO

        IF (.NOT. LWHITBROM) THEN

          DO I=FIRST_VALID_PT,LAST_P_VALID_PT
           RS_P_K(I)=A
          ENDDO

        ELSE

          DO I=FIRST_VALID_PT,LAST_P_VALID_PT
            WORK(I)=RS_P_K(I) ! ie. from the last iteration or
!                             ! junk for the 1st iteration
          ENDDO

! On the first iteration (ie. first level) , the WORK array is
! ignored by CALC_RS. On subsequent iterations it will contain the
! value of RS_P_K of the level under the current level.
          CALL CALC_RS(PSTAR(FIRST_VALID_PT),AK,BK,TS(FIRST_VALID_PT),
     &                 WORK(FIRST_VALID_PT),RS_P_K(FIRST_VALID_PT),
     &                 LAST_P_VALID_PT-FIRST_VALID_PT+1,K,P_LEVELS,
     &                 LLINTS)

        ENDIF ! IF (.NOT. LWHITBROM)

        CALL P_TO_UV(DELP_P,DELP_UV,P_FIELD,U_FIELD,ROW_LENGTH,
     &               tot_P_ROWS)
        CALL P_TO_UV(RS_P_K,RS_UV_K,P_FIELD,U_FIELD,ROW_LENGTH,
     &               tot_P_ROWS)

! 2.2 Now do the sums. For each sum, first we calculate the
!     numbers to be summed, using CALC_?_SUM_ARRAY (?=ENERGY or MASS)
!     and add them into the SUM_ARRAY

! 2.2.1 Energy : CP*TL
        CALL CALC_ENERGY_SUM_ARRAY(TL(1,K),AREA_P,DELP_P,RS_P_K,CP,
     &                             P_FIELD,START_POINT_NO_HALO,
     &                             END_P_POINT_NO_HALO,
     &                             SUM_ARRAY(1,energy))

! 2.2.2 Energy : 0.5*U*U
        DO I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
          WORK(I)=U(I,K)*U(I,K)
        ENDDO

        CALL CALC_ENERGY_SUM_ARRAY(WORK,AREA_UV,DELP_UV,RS_UV_K,0.5,
     &                             U_FIELD,START_POINT_NO_HALO,
     &                             END_U_POINT_NO_HALO,
     &                             SUM_ARRAY(1,energy))

! 2.2.3 Energy : 0.5*V*V
        DO I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
          WORK(I)=V(I,K)*V(I,K)
        ENDDO

        CALL CALC_ENERGY_SUM_ARRAY(WORK,AREA_UV,DELP_UV,RS_UV_K,0.5,
     &                             U_FIELD,START_POINT_NO_HALO,
     &                             END_U_POINT_NO_HALO,
     &                             SUM_ARRAY(1,energy))


! 2.2.4 Mass : PSTAR*DELBK*(-DELAK) (total mass) :
        CALL CALC_MASS_SUM_ARRAY(DELP_P,AREA_P,RS_P_K,
     &                           P_FIELD,START_POINT_NO_HALO,
     &                           END_P_POINT_NO_HALO,
     &                           SUM_ARRAY(1,total_mass))



! 2.2.5 Mass : PSTAR*DELBK (partial mass) :
        CALL CALC_MASS_SUM_ARRAY(PSTAR_DELBK,AREA_P,RS_P_K,
     &                           P_FIELD,START_POINT_NO_HALO,
     &                           END_P_POINT_NO_HALO,
     &                           SUM_ARRAY(1,partial_mass))

      ENDDO ! K : loop over levels

! 2.3 Now finally do the global sums

! 2.3.1 Zero mass and energy of atmosphere
      TOT_MASS_P=0.0
      PART_MASS_P=0.0
      TOT_ENERGY=0.0

! 2.3.2 And do the sums:

      CALL DO_SUMS(SUM_ARRAY,P_FIELD,START_POINT_NO_HALO,
     &             END_P_POINT_NO_HALO,N_SUMS,SUM_RESULTS)

      TOT_ENERGY=SUM_RESULTS(energy)
      TOT_MASS_P=SUM_RESULTS(total_mass)
      PART_MASS_P=SUM_RESULTS(partial_mass)

      RETURN
      END
