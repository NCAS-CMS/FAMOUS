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
CLL  SUBROUTINE CAL_ENG_MASS_CORR--------------------------------------
CLL
CLL  PURPOSE : PART OF ENERGY CORRECTION SUITE OF ROUTINES
CLL            - TO CALCUATE THE NECESSARY CORRECTION TO
CLL              TEMPERATURE TO CONSERVE TOTAL ENERGY
CLL
CLL  NOT SUITABLE FOR SINGLE COLUMN MODEL USE
CLL  CODE WRITTEN FOR CRAY Y-MP BY D.GREGORY FEBRUARY 1991
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
!LL  4.1     23/04/96  MPP code : Only 1 processor to write diagnostic
!LL                    information.       P.Burton
!LL  4.3     29/04/97  Correct Write Statement. D. Robinson
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL TASK :
CLL
CLL  PROJECT TASK :
CLL
CLL  DOCUMENTATION :
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE CAL_ENG_MASS_CORR (TOT_FLUXES,TOT_ENERGY_INIT,
     1                              TOT_ENERGY_FINAL,TOT_MASS_INIT,
     2                              TOT_MASS_FINAL,PART_TOT_MASS,
     3                              P_FIELD,NPNTS,ENERGY_CORR,PSTAR,
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
     4                              DX,DY)
C
C
      IMPLICIT NONE
C
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

C*L------------------COMDECK C_PI---------------------------------------
CLL
CLL 4.0 19/09/95  New value for PI. Old value incorrect
CLL               from 12th decimal place. D. Robinson
CLL
      REAL PI,PI_OVER_180,RECIP_PI_OVER_180

      PARAMETER(
     & PI=3.14159265358979323846, ! Pi
     & PI_OVER_180 =PI/180.0,   ! Conversion factor degrees to radians
     & RECIP_PI_OVER_180 = 180.0/PI ! Conversion factor radians to
     &                              ! degrees
     & )
C*----------------------------------------------------------------------

C----------------------------------------------------------------------
C VECTOR LENGTHS
C----------------------------------------------------------------------
C
      INTEGER P_FIELD          ! IN VECTOR LENGTH OF P GRID
C
      INTEGER NPNTS            ! IN VECTOR LENGTH OF CALCULATIONS
C
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

C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
C
      REAL TOT_FLUXES           ! IN TOTAL OF FLUXES THROUGH ATMOSPHERE
C
      REAL TOT_ENERGY_INIT      ! IN TOTAL ENERGY OF ATMOSPHERE
                                ! AT START OF MODEL
C
      REAL TOT_MASS_INIT        ! IN TOTAL MASS OF ATMOSPHERE
                                ! AT START OF MODEL
C
      REAL TOT_ENERGY_FINAL     ! IN TOTAL ENERGY OF ATMOSPHERE
                                ! AT END OF ONE DAY
C
      REAL TOT_MASS_FINAL       ! IN TOTAL MASS OF ATMOSPHERE
                                ! AT END OF ONE DAY
C
      REAL PART_TOT_MASS        ! IN FACTOR TO CALCULATE PSTAR
                                ! CORRECTION
C
      REAL DX                   ! EW GRID SPACING IN DEGREES
C
      REAL DY                   ! NS GRID SPACING IN DEGREES
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE IN AND OUT
C----------------------------------------------------------------------
C
      REAL PSTAR(P_FIELD)       ! INOUT PRESSURE
C
      REAL ENERGY_CORR          ! INOUT ENERGY CORRECTION FACTOR
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C----------------------------------------------------------------------
C
      REAL CHG_ENERGY          ! CHANGE IN ENERGY
C
      REAL ERROR_ENERGY        ! ERROR IN ENERGY CALCULATION
C
      REAL PSTAR_CORR          ! CORRECTION FACTOR TO PSTAR
C
C----------------------------------------------------------------------
C INTERNAL LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER I                ! LOOP COUNTER
C
C
C----------------------------------------------------------------------
C EXTERNAL SUBROUTINE CALLS  -  NONE
C----------------------------------------------------------------------
C
C*---------------------------------------------------------------------
C
C======================================================================
C CALCULATION OF TEMPERATURE CHANGE TO BE ADDED OVER THE NEXT DAY
C TO CORRECT FOR ERROR IN ENERGY BUDGET DURING PRESENT DAY
C======================================================================
C
C
C----------------------------------------------------------------------
C CALCULATE ENERGY CHANGE DURING THE DAY
C----------------------------------------------------------------------
C
      CHG_ENERGY = TOT_ENERGY_FINAL - TOT_ENERGY_INIT
C
C
C----------------------------------------------------------------------
C CALCULATE ERROR = DIFFERENCE BETWEEN CHANGE IN TOTAL ENERGY
C DURING THE DAY AND THAT EXPECTED DUE TO THE FLUXES OF ENERGY INTO
C THE ATMOSPHERE
C-----------------------------------------------------------------------
C
      ERROR_ENERGY = TOT_FLUXES - CHG_ENERGY
C
C
C-----------------------------------------------------------------------
C CALCULATE TEMPERATURE CHANGE TO BE APPLIED OVER THE
C NEXT DAY TO CORRECT FOR THE ERROR
C-----------------------------------------------------------------------
C
      ENERGY_CORR = ERROR_ENERGY / (CP*TOT_MASS_INIT)
C
C
C----------------------------------------------------------------------
C CALCULATE RATE OF TEMPERATURE CHANGE WHICH NEEDE TO BE
C APPLIED OVER NEXT DAY TO CORRECT FOR ERROR
C----------------------------------------------------------------------
C
      ENERGY_CORR = ENERGY_CORR / 86400.0
C
C
C----------------------------------------------------------------------
C DIAGNOSTICS
C----------------------------------------------------------------------
C
      IF (MY_PROC_ID .EQ. 0) THEN
      WRITE(6,100) TOT_ENERGY_FINAL,TOT_ENERGY_INIT,
     1             CHG_ENERGY,TOT_FLUXES,ERROR_ENERGY,
     2             ENERGY_CORR*86400.0,ENERGY_CORR,
     3             (PI*ERROR_ENERGY*DX*DY)/(A*A*360.0*360.0*86400.0)
C
 100  FORMAT (1X,'FINAL TOTAL ENERGY           = ',E13.5,' J/ '/
     1        1X,'INITIAL TOTAL ENERGY         = ',E13.5,' J/ '/
     2        1X,'CHG IN TOTAL ENERGY OVER DAY = ',E13.5,' J/ '/
     3        1X,'FLUXES INTO ATM OVER DAY     = ',E13.5,' J/ '/
     4        1X,'ERROR IN ENERGY BUDGET       = ',E13.5,' J/ '/
     5        1X,'TEMP CORRECTION OVER DAY     = ',E13.5,' K    '/
     6        1X,'TEMPERATURE CORRECTION RATE  = ',E13.5,' K/S  '/
     7        1X,'FLUX CORRECTION (ATM)        = ',E13.5,' W/M2 ')
      ENDIF
C
C
C======================================================================
C CORRECTION OF PSTAR TO ENSURE MASS CONSERVATION
C======================================================================
C
C
C----------------------------------------------------------------------
C CALCULATE CORRECTION TO BE APPLIED TO PSTAR
C----------------------------------------------------------------------
C
      PSTAR_CORR = 1.0 + (TOT_MASS_INIT - TOT_MASS_FINAL) /
     1                                           PART_TOT_MASS
C
      WRITE(6,200) TOT_MASS_FINAL,TOT_MASS_INIT,
     1             PSTAR_CORR
C
 200  FORMAT (1X,'FINAL ATM MASS               = ',E13.5,' KG '/
     1        1X,'INITIAL ATM MASS             = ',E13.5,' KG '/
     2        1X,'CORRECTION FACTOR FOR PSTAR  = ',E13.5)
C
C
C----------------------------------------------------------------------
C CORRECT PSTAR TO CONSERVE GLOBAL MASS OF ATMOSPHERE
C----------------------------------------------------------------------
C
      DO I=1,NPNTS
       PSTAR(I) = PSTAR(I) * PSTAR_CORR
      END DO
C
      RETURN
      END
