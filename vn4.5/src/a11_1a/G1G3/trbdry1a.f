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
C*LL  SUBROUTINE TRBDRY-------------------------------------------------
CLL
CLL  Purpose: Special routine to add psuedo source terms to boundary
CLL           data in limited area model.
CLL  Method:  runs along each boundary setting boundary concentration
CLL           if there is inflow. The boundary concentration is
CLL           set using a call to function BDRYV. The function BDRYV
CLL           is specific to a model configuration: the current
CLL           version (3.4) is specific to UK MES.
CLL           On outflow boundaries the concentration is set equal
CLL           to the nearest gridpoint inside the boundary.
CLL
CLL
CLL Pete Clark  <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.4:
CLL version  Date
CLL  4.2  15/08/96  Add MPP code. Remove unused variables.  RTHBarnes.
CLL
CLL  Programming standard: Unified Model Documentation Paper No 3,
CLL                        Version 7, dated 11/3/93.
CLL
CLL
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE TRBDRY(
     & AK,BK,
     & POINTS,PFIELD,UFIELD,ROW_LENGTH,
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
     & PSTAR,
     & U,V,
     & TR,TIMESTEP,ERROR
     &)
      IMPLICIT NONE
      INTEGER
     & POINTS              ! IN No. of gridpoints being processed.
     &,PFIELD              ! IN No. of points in global field (at one
C                          !    vertical level).
     &,UFIELD              ! IN No. of u points in global field (at one
C                          !    vertical level).
     &,ROW_LENGTH          ! IN Length of a row.
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
     & AK,BK             ! IN Layer  ak and bk
     &,PSTAR(PFIELD)     ! IN Surface pressure
     &,U(UFIELD),V(UFIELD) ! IN U and V component of wind
     &,TR(PFIELD)        ! INOUT Tracer field (kg per kg air).   
     &,TIMESTEP          ! IN Timestep in seconds
      INTEGER ERROR      ! OUT Error return code.
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

C
C*L  Workspace usage----------------------------------------------------
C*L  External subroutine called ----------------------------------------
C     None
C     EXTERNAL
C* Local, including SAVE'd, storage------------------------------------
C
C  (a) Scalars effectively expanded to workspace by the Cray (using
C      vector registers).
      REAL
     & DC,WDIR,WSPEED,PRESS 
      REAL
     & BDRYV ! FUNCTION giving boundary value.
C
C  (b) Others.
      INTEGER I,IU   ! Loop counters     
C-----------------------------------------------------------------------
C  Check input arguments for potential over-writing problems.
C-----------------------------------------------------------------------
      ERROR=0
      IF(POINTS.GT.PFIELD)THEN
        ERROR=1
        WRITE(6,*)'Error in TRBDRY: POINTS greater than PFIELD.'
        GOTO 9999
      ENDIF
C
      IF (at_top_of_LPG) THEN
C-----------------------------------------------------------------------
CL Loop across top row.
C-----------------------------------------------------------------------
C
      DO  I = TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
        IU=I
        IF(V(IU) .LT. 0.0)THEN
          PRESS=AK+BK*PSTAR(I)
          WDIR=ATAN2(V(IU),U(IU))*RECIP_PI_OVER_180
          WSPEED=SQRT(U(IU)*U(IU) + V(IU)*V(IU))
          DC = BDRYV(WDIR,WSPEED,PRESS)
          TR(I) = DC
        ELSE
          TR(I) = TR(I+ROW_LENGTH)
        ENDIF
      ENDDO ! Loop over points
C
      ENDIF
C
      IF (at_base_of_LPG) THEN
C-----------------------------------------------------------------------
CL Loop across bottom row.
C-----------------------------------------------------------------------
C
      DO  I = P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
        IU = I-ROW_LENGTH
        IF(V(IU) .GT. 0.0)THEN
          PRESS=AK+BK*PSTAR(I)
          WDIR=ATAN2(V(IU),U(IU))*RECIP_PI_OVER_180
          WSPEED=SQRT(U(IU)*U(IU) + V(IU)*V(IU))
          DC = BDRYV(WDIR,WSPEED,PRESS)
          TR(I) = DC
        ELSE
          TR(I) = TR(I-ROW_LENGTH)
        ENDIF
      ENDDO ! Loop over points
C
      ENDIF
C                                                                       
      IF (at_left_of_LPG) THEN
C-----------------------------------------------------------------------
CL Loop across left column
C-----------------------------------------------------------------------
C
      DO  I = TOP_ROW_START+FIRST_ROW_PT-1,
     &        P_BOT_ROW_START+FIRST_ROW_PT-1,ROW_LENGTH
        IF (I .eq. P_BOT_ROW_START+FIRST_ROW_PT-1) THEN
         IF (at_base_of_LPG) THEN
          IU = I-ROW_LENGTH
         END IF
        ELSE
          IU = I
        END IF
        IF(U(IU) .GT. 0.0)THEN
          PRESS=AK+BK*PSTAR(I)
          WDIR=ATAN2(V(IU),U(IU))*RECIP_PI_OVER_180
          WSPEED=SQRT(U(IU)*U(IU) + V(IU)*V(IU))
          DC = BDRYV(WDIR,WSPEED,PRESS)
          TR(I) = DC
        ELSE
          TR(I) = TR(I+1)
        ENDIF
      ENDDO ! Loop over points
C
      ENDIF
C                                                                       
      IF (at_right_of_LPG) THEN
C-----------------------------------------------------------------------
CL Loop across right column
C-----------------------------------------------------------------------
C
      DO  I = TOP_ROW_START+LAST_ROW_PT-1,
     &        P_BOT_ROW_START+LAST_ROW_PT-1,ROW_LENGTH
        IF (I .eq. P_BOT_ROW_START+LAST_ROW_PT-1) THEN
         IF (at_base_of_LPG) THEN
          IU = I-ROW_LENGTH
         END IF
        ELSE
          IU = I
        END IF
        IF(U(IU) .LT. 0.0)THEN
          PRESS=AK+BK*PSTAR(I)
          WDIR=ATAN2(V(IU),U(IU))*RECIP_PI_OVER_180
          WSPEED=SQRT(U(IU)*U(IU) + V(IU)*V(IU))
          DC = BDRYV(WDIR,WSPEED,PRESS)
          TR(I) = DC
        ELSE
          TR(I) = TR(I-1)
        ENDIF
      ENDDO ! Loop over points
C                                                                       
      ENDIF
C
 9999 CONTINUE ! Error exit
      RETURN
      END
C*LL  FUNCTION BDRYV----------------------------------------------------
CLL
CLL  Purpose: Special routine to add psuedo source terms to boundary
CLL           data in limited area.
CLL  PARAMETERS ARE SPECIFIC TO UK MESOSCALE MODEL
CLL  Method:  The boundary concentrations are computed using a
CLL           simple model of transport from sources outside the
CLL           model. Analysis of the source distribution outside
CLL           the UK MES shows that it can be well represented by
CLL           a line source at constant radius from the centre of
CLL           the model, with a source distribution given by the
CLL           sum of two Gaussians. Concentrations from these are
CLL           computed assuming transport using the local windspeed u
CLL           or 1 m/s, whichever is stronger, over a distance
CLL           determined from the centroid of the source distribution, x
CLL           with a linear transformation rate k from emission to
CLL           aerosol, dry deposition at a rate determined from the
CLL           dry deposition velocity vd and mean mixed layer depth h.
CLL           Thus the max concentration is given by
CLL               Q/(uh)*k/(k+vd/h)*(1-exp(-k*x/u))
CLL           The source term is assumed to decrease with level
CLL           pressure. See forthcoming documentation for details.
CLL
CLL Pete Clark  <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.4:
CLL version  Date
CLL
CLL  Programming standard: Unified Model Documentation Paper No 3,
CLL                        Version 7, dated 11/3/93.
CLL
CLL
C*L  Arguments:---------------------------------------------------------
      REAL FUNCTION BDRYV(
     & WDIR
     &,WSPEED
     &,PRESS
     &)
      REAL
     & WDIR          ! IN Wind direction : Cartesian degrees
     &,WSPEED        ! IN Wind speed m/s
     &,PRESS         ! IN Pressure
C*
C* Local, including SAVE'd, storage------------------------------------
C
C  (a) Scalars effectively expanded to workspace by the Cray (using
C      vector registers).
      REAL
     & ZANGLE1,WIDTH1 ! Centre and width of first source Gaussian.
     & ZANGLE2,WIDTH2 ! Centre and width of second source Gaussian.
     &,CMAX,CZER ! Max concentration and 'background'.
     &,WDIRN,RECIPROOT2PI
     &,MIXD,TRAVEL    ! Average mixed layer depth and travel distance.
     &,QMAX1          ! Peak height of first source Gaussian.
     &,QMAX2          ! Peak height of second source Gaussian.
     &,VD             ! Dry deposition velocity.
     &,K,K1 ! Transformation parameters.
     &,PH ! Pressure height scale.
     &,KRAT,KT
      PARAMETER(ZANGLE1=178.0)
      PARAMETER(WIDTH1=5.0)
      PARAMETER(ZANGLE2=173.0)
      PARAMETER(WIDTH2=25.0)
      PARAMETER(RECIPROOT2PI=0.3989422803)
      PARAMETER(QMAX1=4.7E5,QMAX2=9.0E4,MIXD=800.0,TRAVEL=7.7E5)
      PARAMETER(CZER=6.0)
      PARAMETER(VD=5.0E-3)
      PARAMETER(K=3.0E-6, K1=K+VD/MIXD, KRAT=K/K1,KT=-K1*TRAVEL)
      PARAMETER(PH=3.0E4)
C
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

!
!     Max concentration = Q/(uh)*k/(k+vd/h)*(1-exp(-k*x/u))
!
      CMAX=1.0/MAX(WSPEED,1.0)/MIXD
      CMAX=CMAX*KRAT*(1-EXP(KT/MAX(WSPEED,1.0)))
      WDIRN = WDIR - ZANGLE1
      IF (WDIRN .LT. -180.0) WDIRN=WDIRN+360.0
      IF (WDIRN .GT.  180.0) WDIRN=WDIRN-360.0
      WDIRN=WDIRN/WIDTH1
      BDRYV= QMAX1 * EXP(-WDIRN*WDIRN/2.0)
      WDIRN = WDIR*RECIP_PI_OVER_180 - ZANGLE2
      IF (WDIRN .LT. -180.0) WDIRN=WDIRN+360.0
      IF (WDIRN .GT.  180.0) WDIRN=WDIRN-360.0
      WDIRN=WDIRN/WIDTH2
      BDRYV= BDRYV + QMAX2 * EXP(-WDIRN*WDIRN/2.0)
!
!     Add 'background' value.
!
      BDRYV= BDRYV * CMAX + CZER
!
!     Reduce concentration with pressure altitude.
!
      BDRYV= BDRYV * EXP(-(1.E5-PRESS)/PH)
      RETURN
      END
