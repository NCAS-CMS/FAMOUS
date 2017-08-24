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
CLL Subroutine theta_pv ------------------------------------------------
CLL
CLL Purpose: To compute Potential Temperature (Theta) on
CLL          Potential Vorticity surfaces.
CLL          Subroutine can cope with an array of desired pv surfaces.
CLL          Includes a call to subroutine CALC_PV_P to calculate
CLL          pv on pressure levels, and to output the array of
CLL          values of Theta on pressure levels.
CLL          This uses the Quasi-Hydrostatic equations, with complete
CLL          representation of the Coriolis terms, and no metric
CLL          terms omitted.
CLL          The shallow atmosphere approximation is not made.
CLL          Under UPDATE identifier GLOBAL, the data is
CLL          assumed periodic along the rows. Note that because
CLL          it is a diagnostic routine, care needs to be taken
CLL          with missing data.
CLL
CLL Not suitable for single column use.
CLL
CLL  Model            Modification history:
CLL Version   Date
CLL   3.1   21/01/93  Written by Simon Anderson.
CLL   3.1   18/01/93  New deck at the release of Version 3.1.
CLL   3.2   28/07/93  Change subroutine name to uppercase and array
CLL                   switch to switch1 for portability.  Tracey Smith
CLL   3.4   27/05/94  Argument LLINTS added and passed to CALC_PV_P
CLL                                                     S.J.Swarbrick
CLL   4.3   21/03/97  MPP changes. S.D.Mullerworth
CLL
CLL
CLL Programming Standard: UM DOC Paper3, Version 4 (05/02/92)
CLL
CLL Logical Component Covered: D415
CLL
CLL System Task: D4
CLL
CLL Documentation: U.M.D.P No.13. Derivation and Calculation of
CLL                Unified Model Potential Vorticity.
CLL                By Simon Anderson and Ian Roulstone.
CLL
CLLEND
C
C*L ARGUMENTS: ---------------------------------------------------------
      SUBROUTINE THETA_PV
     1                   (pstar,theta,u,v,p_field,u_field,
     2                    p_levels,row_length,
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
     3                    rmdi,ak,bk,f3,
     & e_levels,n_levels,dth_dph,
     4                    theta_pv_levs,des_pv,theta_pv_p_levs,des_p,
     5                    latitude_step_inverse,longitude_step_inverse,
     6                    cos_u_latitude,sec_p_latitude,
     7                    theta_on_pv,llints)


      implicit none
      logical  llints

C Input variables ------------------------------------------------------

      integer
     & p_field                 !IN    Size of field on pressure points.
     &,u_field                 !IN    Size of field on velocity points.
     &,p_levels                !IN    Number of pressure levels.
     & ,n_levels           !IN Number of half levels for dtheta/dp
     &,row_length              !IN    Number of points in a row.
     &,theta_pv_levs           !IN    Number of desired pv surfaces.
     &,theta_pv_p_levs         !IN    Number of desired pressure levels.

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
     &,f3(u_field)             !IN    Coriolis term.
     & ,e_levels(n_levels)       !IN half-levels over range
     & ,dth_dph(p_field,n_levels)  !IN dtheta/dp half-levels
     &,des_pv(theta_pv_levs)   !IN    Value(s) of p.v. we want theta on.
     &,des_p(theta_pv_p_levs)  !IN    Values of pressure we want pv on.
     &,latitude_step_inverse   !IN    1/latitude increment.
     &,longitude_step_inverse  !IN    1/longitude increment.
     &,cos_u_latitude(u_field) !IN    Cosine of latitude on uv-grid.
     &,sec_p_latitude(p_field) !IN    Secant of latitude on p-grid.


C Output variables -----------------------------------------------------

      real
     & theta_on_pv(p_field,theta_pv_levs)
     &                         !  OUT Value of potential temperature
     &                         !      on p.v. surface with
     &                         !      p.v.=des_pv.

C*----------------------------------------------------------------------
C*L Workspace Usage: 4 arrays are required.
      real
     & pvort_p(p_field,theta_pv_p_levs)
     &                         ! Calculated field of p.v. on pressure
     &                         ! levels, from subroutine CALC_PV_P.
     &,theta_on_press(p_field,theta_pv_p_levs)
     &                         ! Calculated field of theta on pressure
     &                         ! levels, from subroutine CALC_PV_P.
     &,f3_p(p_field)           ! Interpolated f3 field on p grid.

      integer
     & switch1(p_field)        ! Switch to determine hemispheric
     &                         ! dependance.
C*----------------------------------------------------------------------
C*L External subroutine calls:
      external calc_pv_p       ! Compute p.v. on pressure levels.
      external th_pvint        ! Interpolate theta onto pv surfaces.
      external uv_to_p         ! Interpolate u-grid field to p-grid fld.

C*----------------------------------------------------------------------
C*L Define local variable.
      integer i,k            ! Loop variable.
      real mn                ! Mean value used in computing pole values.
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
      integer info

C ----------------------------------------------------------------------
CL Section 1 Compute p.v. on each of the desired pressure levels.
CL ~~~~~~~~~ Set switch array to check for sign of desired pv values.
C ----------------------------------------------------------------------

C Loop from 1 to theta_pv_p_levs.

      do 100 k=1,theta_pv_p_levs

C The array des_p is assumed to be ordered
C from the bottom of the atmosphere to the top.

        call calc_pv_p
     1                (pstar,theta,u,v,p_field,u_field,
     2                 p_levels,row_length,
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
     3                 rmdi,ak,bk,des_p(k),f3,
     & e_levels,n_levels,dth_dph,
     4                 latitude_step_inverse,longitude_step_inverse,
     5                 cos_u_latitude,sec_p_latitude,
     6                 pvort_p(1,k),theta_on_press(1,k),llints)

 100  continue

C Find f3_p on p-points.
C MPP: One fewer row than normal in P_field.
      call uv_to_p
     &            (f3,f3_p(row_length+1),u_field,p_field-row_length,
     &             row_length,u_field/row_length)
      CALL SWAPBOUNDS(f3_p,row_length,tot_P_ROWS,
     &  EW_Halo,NS_Halo,1) 

C Calculate f3_p at the poles.
      mn = 0.
      IF (at_top_of_LPG) THEN
        CALL GCG_RVECSUMR(ROW_LENGTH,ROW_LENGTH-2*EW_Halo,
     &    EW_Halo+TOP_ROW_START,1,f3,GC_ROW_GROUP,info,mn)
        mn=mn/GLOBAL_ROW_LENGTH
        do i=TOP_ROW_START,START_POINT_NO_HALO-1
          f3_p(i) = mn
        end do
      ENDIF

      mn = 0.
      IF (at_base_of_LPG) THEN
        CALL GCG_RVECSUMR(ROW_LENGTH,ROW_LENGTH-2*EW_Halo,
     &    EW_Halo,1,f3(U_BOT_ROW_START),GC_ROW_GROUP,info,mn)
        mn=mn/GLOBAL_ROW_LENGTH
        do i=P_BOT_ROW_START,LAST_P_FLD_PT-1
          f3_p(i) = mn
        end do
      ENDIF

C Set hemispheric switch array.
      do 101 i=1,p_field
        if (f3_p(i).ge.0.) then
          switch1(i) = 1
        else
          switch1(i) = -1
        endif
 101  continue

C Change sign of pv field in Southern hemisphere. This is required
C since user asks for values of Des_pv with Northern hemisphere sign,
C namely positive. We need to either change the sign of this Des_pv
C value depending on the hemispere we are in, or, as we do here,
C change the sign of the potential vorticity instead.

      do 102 k=1,theta_pv_p_levs
        do 103 i=1,p_field
          if (pvort_p(i,k).ne.rmdi) then
            pvort_p(i,k) = pvort_p(i,k)*switch1(i)
          endif
 103    continue
 102  continue

C ----------------------------------------------------------------------
CL Section 2 Compute theta on each pv surface using potential vorticity
CL ~~~~~~~~~ on pressure levels, and theta on pressure levels.
CL           Note that all the missing data points at limited area
CL           model boundaries, and polar rows are taken care of
CL           in the CALC_PV_P subroutine.
C ----------------------------------------------------------------------

      do 200 k=1,theta_pv_levs

        call th_pvint
     1               (theta_on_press,pvort_p,
     2                p_field,theta_pv_p_levs,rmdi,des_pv(k),
     3                theta_on_pv(1,k))

 200  continue

      return
      end

