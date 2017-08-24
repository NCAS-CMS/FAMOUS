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
CLL  SUBROUTINE H_GRAD--------------------------------------------------
CLL
CLL  PURPOSE:   Calculates the linear horizontal gradients of a field
CLL  Tested under compiler CFT77
CLL  Tested under OS version 5.1
CLL
CLL  Author J.T.Heming
CLL
CLL  Code version 1.0         Date 22/02/91
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL vn4.2  07/11/96:Allow this routine to be used in C93_2A (Farnon)
CLL
CLL vn4.2  18/03/97:MPP Changes - now hard-wired for U_Fields
CLL                 S.D. Mullerworth
!LL vn4.5  17/04/98 Make sure corner points initialised (FIRST_FLD_PT
!LL                 etc.) so routine is it into line with loop
!LL                 ranges in CAT S.D.Mullerworth
CLL  Logical components covered
CLL  Project TASK:
CLL
CLL  Programming standard: U M DOC  Paper NO. 4,
CLL
CLL  External documentation
CLL
CLLEND------------------------------------------------------------------
C
C*L  ARGUMENTS:---------------------------------------------------------
      SUBROUTINE H_GRAD(
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
C data and constants in
     & N,POINTS,SEC_LATITUDE,ROW_LENGTH,EW_SPACE,NS_SPACE,
C data out
     & DN_DX,DN_DY)
C*
C*L
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      INTEGER
     * POINTS        ! IN  NO OF POINTS ON REQUIRED GRID
     *,ROW_LENGTH    ! IN  NO OF COLUMNS
C-----------------------------------------------------------------------
      REAL
     * N(POINTS)           ! IN  INPUT FIELD ON REQUIRED GRID
     *,DN_DX(POINTS)       ! OUT HORIZONTAL GRADIENT IN THE X-DIRECTION
     *,DN_DY(POINTS)       ! OUT HORIZONTAL GRADIENT IN THE Y-DIRECTION
     *,SEC_LATITUDE(POINTS) ! IN  1/COS(LAT) ON REQUIRED GRID
     *,EW_SPACE             ! IN  LONGITUDE GRID SPACING
     *,NS_SPACE             ! IN  LATITUDE GRID SPACING
     *                      !     BOTH THE ABOVE IN DEGREES
C*
C*L
C-----------------------------------------------------------------------
C Local Variables
C-----------------------------------------------------------------------
      INTEGER
     * I,J                 !  LOOP COUNTERS
C-----------------------------------------------------------------------
      REAL
     * LONG_SI_OVER_TWO_A         ! LONGITUDE STEP INVERSE (IN RADIANS)
     *                            ! DIVIDED BY (2*EARTH'S RADIUS(A))
     *,LAT_SI_OVER_TWO_A          ! LATITUDE STEP INVERSE (IN RADIANS)
     *                            ! DIVIDED BY (2*EARTH'S RADIUS(A))
C-----------------------------------------------------------------------
C Constants required : A=radius of Earth,
C                      RECIP_PI_OVER_180=180/Pi
C-----------------------------------------------------------------------
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
C*----------------------------------------------------------------------
C
C     *N1  -  *N2    *N3     La  = Latitude
C          ^                 Lo  = Longitude
C         dLa                dLa = Latitude interval in radians
C          ^                 dLo = Longitude interval in radians
C     *N4  -  *N5    *N6     A   = Radius of Earth in metres
C                            N1-9= Variable at model grid-points
C     ^--dLo--^
C
C     *N7     *N8    *N9
C
C  dN        N3-N1           dN      N1-N7
C  -- = ---------------      -- = -----------
C  dX   Cos(La)*2*dLo*A      dY     2*dLa*A
C
C
C-----------------------------------------------------------------------
CL 1. Calculate 1/(2*dLo*A) and 1/(2*dLa*A)
C-----------------------------------------------------------------------
      LONG_SI_OVER_TWO_A = 0.5*RECIP_PI_OVER_180/EW_SPACE/A
      LAT_SI_OVER_TWO_A  = 0.5*RECIP_PI_OVER_180/NS_SPACE/A
C=======================================================================
C
C=======================================================================
CL 2.Calculate dN/dX and dN/dY for all points except first and last rows
C-----------------------------------------------------------------------
      DO I=START_POINT_NO_HALO,END_U_POINT_NO_HALO
C
        DN_DX(I)=(N(I+1)-N(I-1))*SEC_LATITUDE(I)*LONG_SI_OVER_TWO_A
C
        DN_DY(I)=(N(I-ROW_LENGTH)-N(I+ROW_LENGTH))*LONG_SI_OVER_TWO_A
C
      ENDDO
C-----------------------------------------------------------------------
CL 3. Calculate dN/dX and dN/dY for first row
CL    dN/dY calculated by one-sided differences (extrapolation)
C-----------------------------------------------------------------------
C
      IF (at_top_of_LPG) THEN
      DO I=FIRST_FLD_PT+1,START_POINT_NO_HALO-1
C
        DN_DX(I)=(N(I+1)-N(I-1))*SEC_LATITUDE(I)*LONG_SI_OVER_TWO_A
C
        DN_DY(I)=2.0*DN_DY(I+ROW_LENGTH)-DN_DY(I+2*ROW_LENGTH)
C
      ENDDO
! Initialise first point on row
      DN_DX(FIRST_FLD_PT)=DN_DX(FIRST_FLD_PT+1)
      DN_DY(FIRST_FLD_PT)=DN_DY(FIRST_FLD_PT+1)
      END IF
C-----------------------------------------------------------------------
CL 4. Calculate dN/dX and dN/dY for last row
CL    dN/dY calculated by one-sided differences (extrapolation)
C-----------------------------------------------------------------------
C
      IF (at_base_of_LPG) THEN
      DO I=U_BOT_ROW_START,LAST_U_FLD_PT-1
C
        DN_DX(I)=(N(I+1)-N(I-1))*SEC_LATITUDE(I)*LONG_SI_OVER_TWO_A
C
        DN_DY(I)=2.0*DN_DY(I-ROW_LENGTH)-DN_DY(I-2*ROW_LENGTH)
C
      ENDDO
! Initialise last point on row
      DN_DX(LAST_U_FLD_PT)=DN_DX(LAST_U_FLD_PT-1)
      DN_DY(LAST_U_FLD_PT)=DN_DY(LAST_U_FLD_PT-1)
      ENDIF
C=======================================================================
C     END OF SUBROUTINE H_GRAD
C=======================================================================
      RETURN
      END
C=======================================================================
