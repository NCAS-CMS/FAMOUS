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
CLL SUBROUTINE pv_pint ------------------------------------------------
CLL
CLL Purpose: Interpolates various fields to a desired pressure level.
CLL          Calculates the derivatives D(rs.u)/D(p), D(rs.v)/D(p)
CLL          and D(theta)/D(p).
CLL          Note that routine assumes that pressure is monotonic.
CLL
CLL Not suitable for single column use.
CLL
CLL  Model            Modification history:
CLL Version   Date
CLL   3.1   12/11/92  Written By Simon Anderson.
CLL   3.1   18/01/93  New deck at the release of version 3.1.
CLL   3.2   28/07/93  Change subroutine name to uppercase for
CLL                   portability.    Tracey Smith
!LL   4.3   17/02/97  Added ARGFLDPT arguments and MPP code  P.Burton
CLL
CLL Programming standard UM DOC Paper 3, Version 4(05/02/92) A
CLL
CLL Logical Component Covered: D415
CLL
CLL Project Task: D4
CLL
CLL Documentation: U.M.D.P. No 13. Derivation and Calculation of
CLL                Unified Model Potential Vorticity.
CLL                By Simon Anderson and Ian Roulstone.
CLL
CLLEND------------------------------------------------------------------

C*L ARGUMENTS:----------------------------------------------------------
      SUBROUTINE PV_PINT
     1                  (pstar,theta,rs,u,v,p_field,u_field,
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
     3                   rmdi,ak,bk,des_press,
     & eta_level,e_levels,dth_dph,n_levels,
     4                   theta_on_press,rs_on_press,
     5                   u_on_press,v_on_press,
     6                   drsu_dp,drsv_dp,dtheta_dp)

      implicit none

C Input variables ------------------------------------------------------

      integer
     & p_field                 !IN  Points in horizontal p field.
     &,u_field                 !IN  Points in horizontal u field.
     &,p_levels                !IN  Number of pressure levels.
     &,row_length              !IN  Number of points in a row.
     & ,n_levels           !IN Number of levels of spline

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
     & pstar(p_field)          !IN  Primary model array for surf. press.
     &,theta(p_field,p_levels) !IN  Primary model array for theta field.
     &,rs(p_field,p_levels)    !IN  Primary model array for rs.
     &,u(u_field,p_levels)     !IN  Primary model array for u field.
     &,v(u_field,p_levels)     !IN  Primary model array for v field.
     & ,e_levels(n_levels)       !IN  half-levels over range
     & ,dth_dph(p_field,n_levels)  !IN dtheta/dp half-levels

      real
     & rmdi                    !IN  Real missing data indicator.
     &,ak(p_levels)            !IN  A coefficient of hybrid coordinates
     &                         !    at full levels.
     &,bk(p_levels)            !IN  B coefficient of hybrid coordinates
     &                         !    at full levels.
     &,des_press               !IN  Desired pressure level we want
     &                         !    variables interpolated onto.


C Output variables -----------------------------------------------------

      real
     & theta_on_press(p_field) !OUT Interpolated theta field on p level.
     &,rs_on_press(p_field)    !OUT Interpolated rs field on press level
     &,u_on_press(u_field)     !OUT Interpolated u field on press level.
     &,v_on_press(u_field)     !OUT Interpolated v field on press level.
     &,drsu_dp(u_field)        !OUT Calculated derivative D(rs.u)/D(p).
     &,drsv_dp(u_field)        !OUT Calculated derivative D(rs.v)/D(p).
     &,dtheta_dp(p_field)      !OUT Calculated derivative D(theta)/D(p).
     & ,eta_level(p_field)     !OUT eta value of theta level

C*----------------------------------------------------------------------
C*L Workspace usage: 4 arrays are required.
      integer
     & base_level_eta(p_field)  ! Level pointer below desired level
     &                         ! level for each atmospheric column (set
     &                         ! to 0 if not found, and p_levels if
     &                         ! des_press is above top level).
     &                         ! Calculated at p-points.
     &,base_level_uv(u_field)  ! Base_level calculated at uv-points.

      real
     & pstar_uv(u_field)       ! Surf. Press. interpolated to uv-points.
     &,rs_uv(u_field,p_levels) ! Rs field interpolated to uv-points.

C*----------------------------------------------------------------------
C*L External Subroutine Calls:
      external p_to_uv         ! Interpolate from p-grid to uv-grid.

C*----------------------------------------------------------------------
C*L Local variables:
      integer i,j,ii,iii       ! Loop counts.

      integer imin,imax
      real zth1,zth2,zp1,zp2,zp3,zp4
      real ze1,ze2,ze3,ze4,den1,den2

C ----------------------------------------------------------------------
CL Section 1 : Find the value of the base level. This is the largest
CL ~~~~~~~~~   value of level such that  pressure(level)<des_press ,
CL             calculated here on p-points.
C ----------------------------------------------------------------------
      do 110 i = FIRST_VALID_PT,LAST_P_VALID_PT
      base_level_eta(i)=0
 110  continue
      do 120 j = 1,p_levels
        do 130 i = FIRST_VALID_PT,LAST_P_VALID_PT
          if (des_press.lt. ak(j) + bk(j)*pstar(i)) then
      base_level_eta(i)=j
          endif
 130    continue
 120  continue
C When this loop is done, base_level_p will be the value of the level
C with the highest value of p LESS than the desired pressure value.
C Base_level_p is set to 0 if no larger value is found.
C Base_level_p is set to p_levels if no smaller value is found.


C ----------------------------------------------------------------------
CL Section 2 : This section will interpolate variables held
CL ~~~~~~~~~   on the p-grid onto the desired pressure level.
CL             Used for THETA_ON_PRESS and RS_ON_PRESS at each point.
C ----------------------------------------------------------------------

C----  Given Theta as a function of p and Des_press lying between
C----  Zp1 and Zp5, calculate Theta at Des_press by fitting the
C----  Quartic Lagrange polynomial through the data and evaluating at
C----  Pressure=Des_press.
C----
C----  If level above p_levels-1 or in bottom boundary layers then
C----  we set the value to missing data.
C----
C----  If the value calculated by the quartic is out of range
C----  of the input values, then no quartic or single-valued quartic
C----  is possible between these points and linear interpolation
C----  is used between the two levels. This is usually caused by an
C----  unstable or almost unstable model profile.
C----
C----  A linear interpolation is used for rs.

      do 210 j= FIRST_VALID_PT,LAST_P_VALID_PT

      if(base_level_eta(j).lt.2.or.base_level_eta(j).gt.p_levels-1)
     & then
      eta_level(j)=rmdi
          theta_on_press(j)=rmdi
          rs_on_press(j)=rmdi
        else


C Reset iii for linear interpolation.
      iii=base_level_eta(j)
      zth1=theta(j,iii)
      zp1=ak(iii)+bk(iii)*pstar(j)
      ze1=0.00001*ak(iii)+bk(iii)
      zth2=theta(j,iii+1)
      zp2=ak(iii+1)+bk(iii+1)*pstar(j)
      ze2=0.00001*ak(iii+1)+bk(iii+1)

C Check to see if value is in range.
      eta_level(j)= (des_press-zp2)*ze1/(zp1-zp2)+
     &              (des_press-zp1)*ze2/(zp2-zp1)
      theta_on_press(j)=(des_press-zp2)*zth1/(zp1-zp2)+
     &              (des_press-zp1)*zth2/(zp2-zp1)
      rs_on_press(j)=(des_press-zp2)*rs(j,iii)/(zp1-zp2)+
     &              (des_press-zp1)*rs(j,iii+1)/(zp2-zp1)
c   reset half-level pointer to level below if desired
c   level is below full level eta
      if(eta_level(j).gt.e_levels(iii))then
        base_level_eta(j)=iii-1
      endif
        end if

 210  continue

C ----------------------------------------------------------------------
CL Section 3 : Interpolate surface pressure and Rs values from
CL ~~~~~~~~~   the p-grid to the uv-grid.
C ----------------------------------------------------------------------

      call p_to_uv
     1            (pstar,pstar_uv,p_field,u_field,
     2             row_length,P_LAST_ROW+1)
      do i=1,p_levels
        call p_to_uv
     1              (rs(1,i),rs_uv(1,i),p_field,u_field,
     2             row_length,P_LAST_ROW+1)
      end do

C ----------------------------------------------------------------------
CL Section 4 : Find the value of the base level. This is the largest
CL ~~~~~~~~~   value of level such that  pressure(level)<des_press,
CL             calculated here on uv-points.
C ----------------------------------------------------------------------
      do 410 i = FIRST_FLD_PT,LAST_U_FLD_PT
        base_level_uv(i) = 0
 410  continue
      do 420 j = 1,p_levels
        do 430 i = FIRST_FLD_PT,LAST_U_FLD_PT
          if (des_press.lt. ak(j) + bk(j)*pstar_uv(i)) then
            base_level_uv(i) = j
          endif
 430    continue
 420  continue
C When this loop is done, base_level_uv will be the value of the level
C with the highest value of p LESS than the desired pressure value.
C Base_level_uv is set to 0 if no larger value is found.
C Base_level_uv is set to p_levels if no smaller value is found.


C ----------------------------------------------------------------------
CL Section 5 : This section will interpolate variables held
CL ~~~~~~~~~   on the uv-grid onto the desired pressure level.
CL             Used for U_ON_PRESS and V_ON_PRESS at each point.
C ----------------------------------------------------------------------

C----  Given U and V as functions of P and Des_press lying
C----  between Zp3 and Zp4, calculate U and V at Des_press by
C----  by linear interpolation.
C----
C----  If level is above p_levels-1 or in bottom boundary layers then
C----  we set the value to missing data.
C----

      do 510 j=FIRST_FLD_PT,LAST_U_FLD_PT

      if(base_level_uv(j).lt.2.or.
     &   base_level_uv(j).gt.p_levels-1) then
      eta_level(j)= rmdi
          u_on_press(j)=rmdi
          v_on_press(j)=rmdi
        else

          iii=base_level_uv(j)

C Calculate pressure values at the required two levels.
          zp3 =ak(iii) + bk(iii)*pstar_uv(j)
          zp4 =ak(iii+1) + bk(iii+1)*pstar_uv(j)

C Calculate U_ON_PRESS and V_ON_PRESS.
          u_on_press(j) = (des_press-zp4)*u(j,iii)/(zp3-zp4) +
     &                    (des_press-zp3)*u(j,iii+1)/(zp4-zp3)
          v_on_press(j) = (des_press-zp4)*v(j,iii)/(zp3-zp4) +
     &                    (des_press-zp3)*v(j,iii+1)/(zp4-zp3)
        end if

 510  continue


C ----------------------------------------------------------------------
CL Section 6 : This section will calculate derivatives held
CL ~~~~~~~~~   on the p-grid on the desired pressure level.
CL             Used for DTHETA_DP at each point.
C ----------------------------------------------------------------------

C----  Given P as a function of Theta and Des_press lying between
C----  Zp1 and Zp4, calculate DTHETA_DP at Des_press by calculating
C----  dtheta_dp using centred finite-differences about each point
C----  zp1 to zp4 and then using cubic Lagrange interpolation.
C----
C----  If level above p_levels-1 or in bottom of boundary layer then
C----  we set the value to missing data.
C----

      do 610 j=FIRST_VALID_PT,LAST_P_VALID_PT

      if(base_level_eta(j).lt.2.or.base_level_eta(j).gt.p_levels-1)
     & then
       eta_level(j)= rmdi
          dtheta_dp(j)=rmdi
      elseif(base_level_eta(j).eq.p_levels-1)then
c   in penultimate half-layer set dtheta/dp to last value
       ii=base_level_eta(j)
       dtheta_dp(j)=dth_dph(j,ii)
      else
       ii=base_level_eta(j)
c  interpolate dtheta/dp to desired level by using
c  half-level values of dtheta/dp
      ze1=e_levels(ii)
      ze2=e_levels(ii+1)
      den1=ze1-ze2
      den2=ze2-ze1
C Calculate DTHETA_DP.
         dtheta_dp(j) =
     &             (eta_level(j)-ze2)*dth_dph(j,ii)/den1
     &            +(eta_level(j)-ze1)*dth_dph(j,ii+1)/den2
        end if

 610  continue


C ----------------------------------------------------------------------
CL Section 7 : This section will calculate derivatives held
CL ~~~~~~~~~   on the uv-grid on the desired pressure level.
CL             Used for DRSU_DP and DRSV_DP at each point.
C ----------------------------------------------------------------------

      do 710 j=FIRST_FLD_PT,LAST_U_FLD_PT

      if(base_level_uv(j).lt.2.or.
     &   base_level_uv(j).gt.p_levels-1) then
      eta_level(j)= rmdi
          drsu_dp(j)=rmdi
          drsv_dp(j)=rmdi
        else

          iii=base_level_uv(j)

C Calculate pressure values at the required two levels.
          zp3 =ak(iii) + bk(iii)*pstar_uv(j)
          zp4 =ak(iii+1) + bk(iii+1)*pstar_uv(j)

C Calculate DRSU_DP and DRSV_DP.
         drsu_dp(j) = (rs_uv(j,iii+1)*u(j,iii+1)-
     &                 rs_uv(j,iii)*u(j,iii))/(zp4-zp3)
         drsv_dp(j) = (rs_uv(j,iii+1)*v(j,iii+1)-
     &                 rs_uv(j,iii)*v(j,iii))/(zp4-zp3)
        end if

 710  continue
! Update halos
      CALL SWAPBOUNDS(theta_on_press,row_length,tot_P_ROWS,
     &  EW_Halo,NS_Halo,1)
      CALL SWAPBOUNDS(rs_on_press,row_length,tot_P_ROWS,
     &  EW_Halo,NS_Halo,1)
      CALL SWAPBOUNDS(u_on_press,row_length,tot_U_ROWS,
     &  EW_Halo,NS_Halo,1)
      CALL SWAPBOUNDS(v_on_press,row_length,tot_U_ROWS,
     &  EW_Halo,NS_Halo,1)
      CALL SWAPBOUNDS(drsu_dp,row_length,tot_U_ROWS,
     &  EW_Halo,NS_Halo,1)
      CALL SWAPBOUNDS(drsv_dp,row_length,tot_U_ROWS,
     &  EW_Halo,NS_Halo,1)
      CALL SWAPBOUNDS(dtheta_dp,row_length,tot_P_ROWS,
     &  EW_Halo,NS_Halo,1)
      CALL SWAPBOUNDS(eta_level,row_length,tot_P_ROWS,
     &  EW_Halo,NS_Halo,1)

      return
      end

