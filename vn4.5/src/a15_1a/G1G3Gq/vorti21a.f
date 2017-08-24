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
CLL   SUBROUTINE VORTIC2 -----------------------------------------------
CLL
CLL   Purpose: To compute 'vorticity2'.
CLL            This is a term that is calculated here to be used later
CLL            in the calculation of potential vorticity.
CLL            Vorticity2 = 1/(rs*rs*cos(phi)) *
CLL                         (D(rs.v)/D(lambda)-D(rs.cos(phi).u)/D(phi)).
CLL
CLL   Not suitable for single column use.
CLL   Version for cray y-mp
CLL
CLL   7/10/92 Written by Simon Anderson
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL 3.1    14/01/93 Code inserted to get round FPP error in LAM mode.
CLL                 This code maybe removed when FPP is fixed. M. Mawson
CLL 3.2    28/07/93 Change subroutine name to uppercase for
CLL                 portability.    Tracey Smith
!LL 4.3    18/02/97 Added ARGFLDPT arguments and MPP code   P.Burton
!LL 4.4    05/09/97 One point not initialised properly. S.D.Mullerworth
CLL
CLL   Programming standard: Unified model documentation paper no. 4,
CLL                         standard b. version 2, dated 18/01/90
CLL
CLL   Logical components covered: D415
CLL
CLL   System task: D4
CLL
CLL   Documentation: U.M.D.P No 13. Derivation and Calculation of
CLL                  Unified Model Potential Vorticity.
CLL                  by Simon Anderson and Ian Roulstone.
CLL
CLLEND------------------------------------------------------------------

C*L ARGUMENTS:----------------------------------------------------------
      SUBROUTINE VORTIC2
     1                  (u_on_theta,v_on_theta,rs_on_theta,
     2                   rs_uv_on_theta,cos_u_latitude,
     3                   sec_p_latitude,vorticity2,
     4                   p_field,u_field,row_length,
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
     5                   latitude_step_inverse,longitude_step_inverse)

      implicit none

C Input variables ------------------------------------------------------

      integer
     & p_field                !IN  Number of points in pressure field.
     &,u_field                !IN  Number of points in velocity field.
     &,row_length             !IN  Number of points per row.
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
     & u_on_theta(u_field)    !IN  Mass weighted u velocity.
     &,v_on_theta(u_field)    !IN  Mass weighted v velocity*
     &                        !                      cos(latitude).
     &,rs_on_theta(p_field)   !IN  Pseudo radius at p points.
     &,rs_uv_on_theta(u_field)!IN  Pseudo radius at u points.
     &,cos_u_latitude(u_field)!IN  Cos(lat) at u points.
     &,sec_p_latitude(p_field)!IN  1/cos(lat) at p points.
     &,latitude_step_inverse  !IN  1/latitude increment.
     &,longitude_step_inverse !IN  1/longitude increment.

C Output variables -----------------------------------------------------

      real
     & vorticity2(p_field)    !OUT Term used in potential vorticity eqn.

C*----------------------------------------------------------------------
C*L Workspace usage:- 3 local arrays required.

      real
     & drsv_dlongitude(p_field)
     &,drscosphiu_dlatitude(p_field)
     &,drscosphiu_dlatitude2(u_field)
C*---------------------------------------------------------------------

C*L Define local variables:
      integer
     & i,j                    ! Loop counts.
     &,  info    ! GCOM return code
     
      REAL
     &  pole_sum(row_length)  ! array for summing around pole

      real
     & scalar                 ! Local scalar.

      real sum_n,sum_s



CL---------------------------------------------------------------------
CL    Calculate 'vorticity2'.
CL---------------------------------------------------------------------

C Calculate d(rsv)/d(lambda).
        do 110 i=TOP_ROW_START+1,LAST_U_FLD_PT
          drsv_dlongitude(i) = longitude_step_inverse*
     &                         (rs_uv_on_theta(i)*v_on_theta(i) -
     &                          rs_uv_on_theta(i-1)*v_on_theta(i-1))
 110    continue

C Calculate d(rscosphiu)/d(phi).
        do 120 i=START_POINT_NO_HALO,LAST_U_FLD_PT
          drscosphiu_dlatitude(i) = latitude_step_inverse*
     &                              (rs_uv_on_theta(i-row_length)*
     &                               cos_u_latitude(i-row_length)*
     &                               u_on_theta(i-row_length) -
     &                               rs_uv_on_theta(i)*
     &                               cos_u_latitude(i)*
     &                               u_on_theta(i))
 120    continue

C Calculate average of drscosphiu_dlatitude at p-points.
        do 130 i=START_POINT_NO_HALO+1,END_P_POINT_NO_HALO
          drscosphiu_dlatitude2(i) = drscosphiu_dlatitude(i) +
     &                               drscosphiu_dlatitude(i-1)
 130    continue

C Now do first point on each slice for
C drsv_dlongitude and drscosphiu_dlatitude2.
        i=TOP_ROW_START
! Put a sensible number in the first element (halo)
        drsv_dlongitude(i) = drsv_dlongitude(i+1)
        i=START_POINT_NO_HALO
        drscosphiu_dlatitude2(i)=drscosphiu_dlatitude2(i+1)

C Calculate vorticity2.

        do 150 j=START_POINT_NO_HALO,END_P_POINT_NO_HALO
          vorticity2(j)=sec_p_latitude(j)/
     &                  (rs_on_theta(j)*rs_on_theta(j))*
     &                       .5*(drsv_dlongitude(j)+
     &                           drsv_dlongitude(j-row_length)-
     &                           drscosphiu_dlatitude2(j))
 150    continue


C Calculate vorticity2 at poles by summing drscosphiu/d(lat) around
C Polar circle and averaging.
        scalar = latitude_step_inverse/GLOBAL_ROW_LENGTH
        if (at_top_of_LPG) then
          sum_n=0.0
          do i=1,ROW_LENGTH-2*EW_Halo
            j=TOP_ROW_START+FIRST_ROW_PT+i-2
            pole_sum(i)=-rs_uv_on_theta(j)*cos_u_latitude(j)*
     &                   u_on_theta(j)*scalar
          enddo

          CALL GCG_RVECSUMR(
     &      ROW_LENGTH-2*EW_Halo,ROW_LENGTH-2*EW_Halo,1,1,
     &      pole_sum,
     &      GC_ROW_GROUP,info,sum_n)
          do i=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
            vorticity2(i)=-sum_n*sec_p_latitude(i)/ 
     &                  (rs_on_theta(i)*rs_on_theta(i))
          enddo
        endif
        
        if (at_base_of_LPG) then
          sum_s=0.0
          do i=1,ROW_LENGTH-2*EW_Halo
            j=P_BOT_ROW_START+FIRST_ROW_PT-row_length-2+i
            pole_sum(i)=rs_uv_on_theta(j)*cos_u_latitude(j)*
     &                  u_on_theta(j)*scalar
          enddo

          CALL GCG_RVECSUMR(
     &      ROW_LENGTH-2*EW_Halo,ROW_LENGTH-2*EW_Halo,1,1,
     &      pole_sum,
     &      GC_ROW_GROUP,info,sum_s)
          do i=P_BOT_ROW_START,P_BOT_ROW_START+row_length-1
            vorticity2(i)=-sum_s*sec_p_latitude(i)/ 
     &                  (rs_on_theta(i)*rs_on_theta(i))
          enddo
        endif

! Set rest of array to sensible values
      CALL SWAPBOUNDS(vorticity2,ROW_LENGTH,tot_P_ROWS,
     &              EW_Halo,NS_Halo,1)

CL    End of routine vortic2

      return
      end
