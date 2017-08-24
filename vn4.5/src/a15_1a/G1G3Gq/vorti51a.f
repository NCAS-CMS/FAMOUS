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
CLL   SUBROUTINE VORTIC5 -----------------------------------------------
CLL
CLL   Purpose: To compute   1/(rs*cos(phi)) *
CLL                         (D(v)/D(lambda)-D(u.cos(phi))/D(phi)).
CLL            This is a term that is calculated here to be used later
CLL            in the calculation of potential vorticity
CLL            in subroutine CALC_PV_P.
CLL
CLL   Not suitable for single column use.
CLL
CLL   VERSION FOR CRAY Y-MP
CLL
CLL    Model            Modification history:
CLL   Version   Date
CLL     3.1    2/11/92  Written by Simon Anderson.
CLL     3.1   18/01/93  New deck at the release of Version 3.1.
CLL     3.2   28/07/93  Change subroutine name to uppercase for
CLL                     portability.    Tracey Smith
!LL     4.3   18/02/97  Added ARGFLDPT arguments and MPP code  P.Burton
CLL
CLL   Programming standard: Unified model documentation paper no. 4,
CLL                         standard b. version 2, dated 18/01/90
CLL
CLL   Logical components covered: D415
CLL
CLL   Project task: D4
CLL
CLL   Documentation: U.M.D.P No 13. Derivation and Calculation of
CLL                  Unified Model Potential Vorticity.
CLL                  By Simon Anderson and Ian Roulstone.
CLL
CLLEND------------------------------------------------------------------

C*L ARGUMENTS:----------------------------------------------------------
      SUBROUTINE VORTIC5
     1                  (u_on_press,v_on_press,rs_on_press,
     2                   cos_u_latitude,sec_p_latitude,
     3                   vorticity5,
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
     & u_on_press(u_field)    !IN  Mass weighted u velocity.
     &,v_on_press(u_field)    !IN  Mass weighted v velocity*
     &                        !                      cos(latitude).
     &,rs_on_press(p_field)   !IN  Pseudo radius at p-points.
     &,cos_u_latitude(u_field)!IN  Cos(lat) at uv-points.
     &,sec_p_latitude(p_field)!IN  1/cos(lat) at p-points.
     &,latitude_step_inverse  !IN  1/latitude increment.
     &,longitude_step_inverse !IN  1/longitude increment.

C Output variables -----------------------------------------------------

      real
     & vorticity5(p_field)    !OUT Term used in potential vorticity eqn.

C*----------------------------------------------------------------------
C*L Workspace usage:- 3 local arrays required.

      real
     & dv_dlongitude(p_field)
     &,dcosphiu_dlatitude(p_field)
     &,dcosphiu_dlatitude2(u_field)

C*----------------------------------------------------------------------
C*L External subroutine calls:   None.

C*---------------------------------------------------------------------
C*L Define local variables:
      integer
     & i,j                    ! Loop counts.
     &,  info    ! GC return code

      REAL
     &  pole_sum(row_length)  ! array containing polar vals to sum

      real
     & scalar                 ! Local scalar.

      real sum_n,sum_s


CL---------------------------------------------------------------------
CL    Calculate 'vorticity5'.
CL---------------------------------------------------------------------

C Calculate d(v)/d(lambda).
        do 110 i=TOP_ROW_START+1,LAST_U_FLD_PT
          dv_dlongitude(i) = longitude_step_inverse*
     &                       (v_on_press(i) - v_on_press(i-1))
 110    continue

C Calculate d(cosphiu)/d(phi).
        do 120 i=START_POINT_NO_HALO,LAST_U_FLD_PT
          dcosphiu_dlatitude(i) = latitude_step_inverse*
     &                            (cos_u_latitude(i-row_length)*
     &                             u_on_press(i-row_length) -
     &                             cos_u_latitude(i)*
     &                             u_on_press(i))
 120    continue

C Calculate average of dcosphiu_dlatitude at p-points.
        do 130 i=START_POINT_NO_HALO+1,END_P_POINT_NO_HALO
          dcosphiu_dlatitude2(i) = dcosphiu_dlatitude(i) +
     &                             dcosphiu_dlatitude(i-1)
 130    continue

C Now do first point on each slice for
C dv_dlongitude and dcosphiu_dlatitude2.
! Put a sensible number in the first element (halo)
        i=TOP_ROW_START
        dv_dlongitude(i) = dv_dlongitude(i+1)
        i=START_POINT_NO_HALO
        dcosphiu_dlatitude2(i) = dcosphiu_dlatitude2(i+1)

C Calculate vorticity5.

        do 150 j=START_POINT_NO_HALO,END_P_POINT_NO_HALO
          vorticity5(j)=sec_p_latitude(j)/rs_on_press(j)*
     &                  .5*(dv_dlongitude(j)+
     &                      dv_dlongitude(j-row_length)-
     &                      dcosphiu_dlatitude2(j))
 150    continue


C Calculate vorticity5 at poles by summing dcosphiu/dlatitude around
C Polar circle and averaging.
        scalar = latitude_step_inverse/GLOBAL_ROW_LENGTH
        if (at_top_of_LPG) then
          sum_n=0.0
          do i=1,ROW_LENGTH-2*EW_Halo
            j=TOP_ROW_START+FIRST_ROW_PT+i-2
            pole_sum(i)=-cos_u_latitude(j)*u_on_press(j)*scalar
          enddo

          CALL GCG_RVECSUMR(
     &      ROW_LENGTH-2*EW_Halo,ROW_LENGTH-2*EW_Halo,1,1,
     &      pole_sum,
     &      GC_ROW_GROUP,info,sum_n)
          do i=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
            vorticity5(i)=-sum_n*sec_p_latitude(i)/ rs_on_press(i)
          enddo
        endif

        if (at_base_of_LPG) then
          sum_s=0.0
          do i=1,ROW_LENGTH-2*EW_Halo
            j=P_BOT_ROW_START+FIRST_ROW_PT-row_length-2+i
            pole_sum(i)=cos_u_latitude(j)*u_on_press(j)*scalar
          enddo

          CALL GCG_RVECSUMR(
     &      ROW_LENGTH-2*EW_Halo,ROW_LENGTH-2*EW_Halo,1,1,
     &      pole_sum,
     &      GC_ROW_GROUP,info,sum_s)
          do i=P_BOT_ROW_START,P_BOT_ROW_START+row_length-1
            vorticity5(i)=-sum_s*sec_p_latitude(i)/ rs_on_press(i)
          enddo
        endif
! Set rest of array to sensible values
      CALL SWAPBOUNDS(vorticity5,ROW_LENGTH,tot_P_ROWS,
     &              EW_Halo,NS_Halo,1)

CL    end of routine vortic5

      return
      end

