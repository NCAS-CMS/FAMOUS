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
CLL Subroutine calc_pv -------------------------------------------------
CLL
CLL Purpose: To compute Ertel potential vorticity on isentropic
CLL          surfaces.
CLL          Uses the Quasi-Hydrostatic equations, with complete
CLL          representation of the Coriolis terms, and no metric
CLL          terms omitted.
CLL          The shallow atmosphere approximation is not made.
CLL          Under UPDATE identifier GLOBAL, the data is
CLL          assumed periodic along the rows. Note that because
CLL          it is a diagnostic routine care needs to be taken
CLL          with missing data.
CLL          Subroutine CALCRS1 is called to compute the psuedo radius.
CLL          Values of p, rs, u, v and DTheta/Dp on theta surfaces
CLL          are then computed assuming a cubic variation,
CLL          where possible.
CLL          Two other subroutines then calculate the remaining
CLL          terms in the potential vorticity expression.
CLL          The final part of the calc_pv routine then
CLL          adds in the rotational component and
CLL          checks for missing data.
CLL
CLL Not suitable for single column use.
CLL
CLL  7/10/92 Written by Simon Anderson
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.1   14/01/93  Parameter list tidied up.
CLL   3.2   28/07/93  Change subroutine name to uppercase for
CLL                   portability.    Tracey Smith
CLL   3.4   26/05/94  Logical llints added to s/r args and passed to
CLL                                   calc_rs          S.J.Swarbrick
!LL   4.3   14/02/97  Added ARGFLDPT arguments and MPP code  P.Burton
CLL
CLL Programming Standard: UM DOC Paper3, Version 4 (05/02/92)
CLL
CLL Logical Component Covered: D415
CLL
CLL System Task: D4
CLL
CLL Documentation: U.M.D.P No 13. Derivation and Calculation of
CLL                Unified Model Potential Vorticity.
CLL                by Simon Anderson and Ian Roulstone.
CLL
CLLEND

C*L ARGUMENTS: ---------------------------------------------------------
      SUBROUTINE CALC_PV
     1                  (pstar,theta,u,v,p_field,u_field,
     2                   p_levels,row_length,
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
     3                   rmdi,ak,bk,des_theta,f3,
     & e_levels,n_levels,
     & dth_dph,
     4                   latitude_step_inverse,longitude_step_inverse,
     5                   cos_u_latitude,sec_p_latitude,
     6                   pvort,llints)


      implicit none
      logical  llints

C Input variables ------------------------------------------------------

      integer
     & p_field                 !IN    Size of field on pressure points.
     &,u_field                 !IN    Size of field on velocity points.
     &,p_levels                !IN    Number of pressure levels.
     &,row_length              !IN    Number of points in a row.
     & ,n_levels       !IN Number of levels of dtheta/dp


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
      real
     & pstar(p_field)          !IN    Surface pressure field.
     &,u(u_field,p_levels)     !IN    Primary model array for u field.
     &,v(u_field,p_levels)     !IN    Primary model array for v field.
     &,theta(p_field,p_levels) !IN       "      "     "     theta field.

      real
     & rmdi                    !IN    Real missing data indicator.
     &,ak(p_levels)            !IN    A coefficient of hybrid
     &                         !      coordinates at full levels.
     &,bk(p_levels)            !IN    B coefficient of hybrid
     &                         !      coordinates at full levels.
     & ,e_levels(n_levels)       !IN half-levels over range
     & ,dth_dph(p_field,n_levels)  !IN dtheta/dp half-levels
     &,des_theta               !IN    Value of theta we want pv on.
     &,f3(u_field)             !IN    Coriolis term.
     &,latitude_step_inverse   !IN    1/latitude increment.
     &,longitude_step_inverse  !IN    1/longitude increment.
     &,cos_u_latitude(u_field) !IN    Cosine of latitude on u field.
     &,sec_p_latitude(p_field) !IN    Secant of latitude on p field.


C Output variables -----------------------------------------------------

      real
     & pvort(p_field)          !  OUT Value of potential vorticity
     &                         !      on isentropic surface with
     &                         !      theta=des_theta.


C*----------------------------------------------------------------------
C*L Workspace Usage: 13 arrays are required.
      logical
     & mask(p_field)           ! Logical mask used to make
     &                         ! vectorising of a loop possible.
      real
     & rs(p_field,p_levels)    ! Calculated pseudo radius of Earth.
     &,eta_level(p_field)    ! Interpolated eta-value of theta level
     &,u_on_theta(u_field)     ! Interpolated value of  u on theta surf.
     &,u_p_on_theta(p_field)   ! Interpolated value of  u on theta surf.
     &,v_on_theta(u_field)     ! Interpolated value of  v on theta surf.
     &,v_p_on_theta(p_field)   ! Interpolated value of  v on theta surf.
     &,rs_on_theta(p_field)    ! Interpolated value of rs on theta surf.
     &,rs_uv_on_theta(u_field) ! Interpolated value of rs on theta surf.
     &,dtheta_dp(p_field)      ! D(theta)/D(p) on theta surface.
     &,vorticity1(p_field)     ! 1st Calculated term in the pv equation.
     &,vorticity2(p_field)     ! 2nd Calculated term in the pv equation.
     &,f3_p(p_field)           ! Interpolated f3 field on p field.

C*----------------------------------------------------------------------
C*L External subroutine calls:
      external calc_rs         ! Compute pseudo radius.
      external pv_thint        ! Interpolate variables to theta level.
      external vortic1         ! Compute the 1st vorticity term.
      external vortic2         ! Compute the 2nd vorticity term.
      external uv_to_p         ! Interpolate u grid field to p grid fld.

C*----------------------------------------------------------------------
C*L Call comdecks to get required variables:
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_OMEGA------------------------------------
C OMEGA IS MAGNITUDE OF EARTH'S ANGULAR VELOCITY
      REAL OMEGA

      PARAMETER(OMEGA=7.292116E-5)
C*----------------------------------------------------------------------


C*----------------------------------------------------------------------
C*L Define local variables.
      integer i,j,k          ! Loop variables.
      real mn                ! Mean value used in computing pole values.
      real work1(p_field)    ! General workspace.
      INTEGER info !GCG return code
      INTEGER 
     &  START ! Work variables
     &  ,END


C ----------------------------------------------------------------------
CL Section 1 Compute rs, the pseudo radius.
CL ~~~~~~~~~
C ----------------------------------------------------------------------

CL Section 1.1 Call CALC_RS to get rs for level 1.
CL ~~~~~~~~~~~ Rs is returned in rs(1,k).
CL             Ts is returned in work1. Rs at level k-1 is input as
CL             rs(1,2) as at k-1=0 the input is not used by calc_rs.

      k=1
      call calc_rs
     1            (pstar,ak,bk,work1,rs(1,2),rs(1,k),p_field,k,p_levels,
     2             llints)

CL Section 1.2 Call CALC_RS to get rs for level k.
CL             Rs is returned in rs(1,k).
CL             Ts is returned in work1. Rs at level k-1 is input as
CL             rs(1,k-1).
CL             Loop from 2 to p_levels.

      do 100 k=2,p_levels

        call calc_rs
     1          (pstar,ak,bk,work1,rs(1,k-1),rs(1,k),p_field,k,p_levels,
     2           llints)
 100  continue

C ----------------------------------------------------------------------
CL Section 2 Interpolate p, u, v and rs onto desired theta level,
CL ~~~~~~~~~ and compute Dtheta/Dp assuming cubic variation if possible.
C ----------------------------------------------------------------------

      call pv_thint
     1             (pstar,theta,rs,u,v,p_field,u_field,
     2              p_levels,row_length,
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
     3              rmdi,ak,bk,des_theta,
     4 eta_level,rs_on_theta,rs_uv_on_theta,
     & e_levels,dth_dph,n_levels,
     5              u_on_theta,v_on_theta,dtheta_dp)


C ----------------------------------------------------------------------
CL Section 3 Compute the various terms in the potential vorticity
CL ~~~~~~~~~ equation for desired theta surface.
C ----------------------------------------------------------------------

CL Section 3.1 Compute -2*omega*cos(phi)/rs D(rs)/D(phi).
CL ~~~~~~~~~~~
      call vortic1
     1            (rs_on_theta,
     2             sec_p_latitude,vorticity1,
     3             p_field,row_length,
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
     4             latitude_step_inverse,longitude_step_inverse)


CL Section 3.2 Compute 1/(rs*rs*cos(phi)) *
CL ~~~~~~~~~~~         (D(rs.v)/D(lambda)-D(rs.cos(phi).u)/D(phi)).
      call vortic2
     1            (u_on_theta,v_on_theta,rs_on_theta,
     2             rs_uv_on_theta,cos_u_latitude,
     3             sec_p_latitude,vorticity2,
     4             p_field,u_field,row_length,
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
     5             latitude_step_inverse,longitude_step_inverse)


C ----------------------------------------------------------------------
CL Section 4 Compute the potential vorticity using the vorticity terms,
CL ~~~~~~~~~ D(theta)/D(p) and the value of the coriolis term.
CL           Care needs to be taken for rmdi.
C ----------------------------------------------------------------------

CL Section 4.1 Firstly, find which points we can actually calculate
CL ~~~~~~~~~~~ potential vorticity on. To do this we test for missing
CL             data at all points in a 3x3 stencil around each point.
CL             Special care is taken at the boundaries and poles.

C Set logical mask default to .TRUE.
      do i = FIRST_FLD_PT,LAST_P_FLD_PT
        mask(i) = .true.
      end do

C Check points on non-polar rows:
      do j=FIRST_ROW,P_LAST_ROW
C First point on each row:
        i=(j-1)*row_length+FIRST_ROW_PT
        IF (at_left_of_LPG) THEN
          mask(i+EW_Halo)=.false.
        ENDIF

C Inner points on each row:
C Borders of area done separately
        START=EW_Halo+1
        END=EW_Halo
        IF (at_left_of_LPG) START=START+1
        IF (at_right_of_LPG) END=END+1
        do i=(j-1)*row_length+START,(j-1)*row_length+row_length-END
        if (   eta_level(i-row_length-1).eq.rmdi.
     &      or.eta_level(i-row_length  ).eq.rmdi.
     &      or.eta_level(i-row_length+1).eq.rmdi.
     &      or.eta_level(i-1           ).eq.rmdi.
     &      or.eta_level(i             ).eq.rmdi.
     &      or.eta_level(i+1           ).eq.rmdi.
     &      or.eta_level(i+row_length-1).eq.rmdi.
     &      or.eta_level(i+row_length  ).eq.rmdi.
     &      or.eta_level(i+row_length+1).eq.rmdi)   mask(i)=.false.
        end do
C Last point on each row:
        i=(j-1)*row_length+LAST_ROW_PT
        IF (at_right_of_LPG) THEN
          mask(i-EW_Halo)=.false.
        ENDIF
      end do

C Set mask on Northern and Southern boundaries:
      IF (at_top_of_LPG) THEN
        do i=TOP_ROW_START,TOP_ROW_START+LAST_ROW_PT-1
          mask(i)=.false.
        enddo
      ENDIF ! if at_top_of_LPG

      IF (at_base_of_LPG) THEN
        do i=P_BOT_ROW_START,P_BOT_ROW_START+LAST_ROW_PT-1
          mask(i)=.false.
        enddo
      ENDIF ! if at_base_of_LPG
! Set all the halos to false
      if (at_top_of_lpg) then
        do j=1,NS_Halo
          do i=(j-1)*row_length+1,j*row_length
            mask(i)=.false.
          enddo
        enddo
      endif ! if at_top_of_lpg

      if (at_base_of_lpg) then
        do j=tot_P_ROWS-NS_Halo+1,tot_P_ROWS
          do i=(j-1)*row_length+1,j*row_length
            mask(i)=.false.
          enddo
        enddo
      endif ! if at_base_of_lpg

      if (at_left_of_lpg) THEN
        do j=NS_Halo+1,tot_P_ROWS-NS_Halo
          do i=(j-1)*row_length+1,(j-1)*row_length+EW_Halo
            mask(i)=.false.
          enddo
        enddo
      endif ! if at_left_of_lpg

      if (at_right_of_lpg) THEN
        do j=NS_Halo+1,tot_P_ROWS-NS_Halo
          do i=j*row_length-EW_Halo+1,j*row_length
            mask(i)=.false.
          enddo
        enddo
      endif


CL Section 4.2 Interpolate fields u, v, and F3 (the Coriolis term)
CL ~~~~~~~~~~~ to the p-grid.

      call uv_to_p
     1           (u_on_theta,u_p_on_theta(row_length+1),u_field,
     2            p_field-row_length,row_length,U_LAST_ROW+1)
      call uv_to_p
     1           (v_on_theta,v_p_on_theta(row_length+1),u_field,
     2            p_field-row_length,row_length,U_LAST_ROW+1)
      call uv_to_p
     1           (f3,f3_p(row_length+1),u_field,
     2            p_field-row_length,row_length,U_LAST_ROW+1)

CL Section 4.3 Calculate the potential vorticity.
CL ~~~~~~~~~~~

      do i = START_POINT_NO_HALO,END_P_POINT_NO_HALO
        if (mask(i)) then

          pvort(i) = (vorticity1(i) + f3_p(i) + vorticity2(i)) *
     &               dtheta_dp(i) *
     &               ((2.*omega + u_p_on_theta(i)*
     &                 sec_p_latitude(i)/rs_on_theta(i))*
     &               (u_p_on_theta(i)/sec_p_latitude(i)) +
     &               (v_p_on_theta(i)*v_p_on_theta(i)/rs_on_theta(i))-g)

        else
          pvort(i) = rmdi
        end if
      enddo

CL Section 4.3.1 Compute the potential vorticity at the Northern
CL ~~~~~~~~~~~~~ and Southern boundaries.
CL               Note that under UPDATE identifier GLOBAL, as pv
CL               is a scaler this is defined to be the mean of
CL               the near pole row.
CL               If any rmdi is present in the row next to the pole,
CL               then the pole has a value rmdi.
CL               ELSE, the value returned is missing data.

C Work out the Northern boundary.
      if (at_top_of_LPG) then
        if (mask(TOP_ROW_START+LAST_ROW_PT-1)) then

          mn=0
          call GCG_RVECSUMR(
     &      row_length-2*EW_Halo,row_length-2*EW_Halo,1,1,
     &      pvort(TOP_ROW_START+FIRST_ROW_PT-1+row_length),
     &      GC_ROW_GROUP,info,mn)
          mn=mn/GLOBAL_ROW_LENGTH

          do i=TOP_ROW_START,TOP_ROW_START+row_length-1
            pvort(i)=mn
          enddo
        else
! Missing data
          do i=TOP_ROW_START,TOP_ROW_START+row_length-1
            pvort(i)=rmdi
          enddo
        endif
      endif ! if at_top_of_LPG

! Southern boundary

      if (at_base_of_LPG) then
        if (mask(P_BOT_ROW_START+FIRST_ROW_PT-1)) then

          mn=0
          call GCG_RVECSUMR(
     &      row_length-2*EW_Halo,row_length-2*EW_Halo,1,1,
     &      pvort(P_BOT_ROW_START+FIRST_ROW_PT-1-row_length),
     &      GC_ROW_GROUP,info,mn)
          mn=mn/GLOBAL_ROW_LENGTH

          do i=P_BOT_ROW_START,P_BOT_ROW_START+row_length-1
            pvort(i)=mn
          enddo
        else
! Missing data
          do i=P_BOT_ROW_START,P_BOT_ROW_START+row_length-1
            pvort(i)=rmdi
          enddo
        endif
      endif ! if at_base_of_LPG

! Fill in the N+S halo areas with rmdi
      do i=1,FIRST_FLD_PT-1
        pvort(i)=rmdi
      enddo
      do i=LAST_P_FLD_PT+1,p_field
        pvort(i)=rmdi
      enddo
      CALL SWAPBOUNDS(pvort,ROW_LENGTH,tot_P_ROWS,EW_Halo,
     &  NS_Halo,1)


      return
      end

