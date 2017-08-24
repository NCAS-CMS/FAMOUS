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
CLL   SUBROUTINE VORTIC1 -----------------------------------------------
CLL
CLL   Purpose: To compute 'vorticity1'.
CLL            This is a term that is calculated here to be used later
CLL            in the calculation of potential vorticity.
CLL            Vorticity1 = -2*omega*cos(phi)/rs * D(rs)/D(phi)
CLL   Not suitable for single column use.
CLL
CLL   Version for cray y-mp
CLL
CLL   Written by Simon Anderson
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.2    28/07/93 Change subroutine name to uppercase for
CLL                   portability.    Tracey Smith
!LL   4.3    17/02/97 Added ARGFLDPT arguments and MPP code  P.Burton
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
      SUBROUTINE VORTIC1
     1                  (rs_on_theta,
     2                   sec_p_latitude,vorticity1,
     3                   p_field,row_length,
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
     4                   latitude_step_inverse,longitude_step_inverse)

      implicit none

C Input variables ------------------------------------------------------

      integer
     & p_field                !IN  Number of points in pressure field.
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
     & rs_on_theta(p_field)   !IN  Pseudo radius of earth at p points.
     &,sec_p_latitude(p_field)!IN  1/cos(lat) at p points.
     &,latitude_step_inverse  !IN  1/latitude increment.
     &,longitude_step_inverse !IN  1/longitude increment.

C Output variables -----------------------------------------------------

      real
     & vorticity1(p_field)    !OUT Term used in potential vorticity eqn.

C*----------------------------------------------------------------------
C*L Workspace usage:- 1 local array required.

      real
     & drs_dlatitude(p_field)

C*----------------------------------------------------------------------
C*L Call comdecks to get required variables:
C*L------------------COMDECK C_OMEGA------------------------------------
C OMEGA IS MAGNITUDE OF EARTH'S ANGULAR VELOCITY
      REAL OMEGA

      PARAMETER(OMEGA=7.292116E-5)
C*----------------------------------------------------------------------


C*----------------------------------------------------------------------
C*L Define local variables:
      integer
     & i,j                    ! Loop counts.
     &,  info   ! GCOM return code

      real
     & scalar                 ! Local scalar.

      real sum_n,sum_s


C ----------------------------------------------------------------------
CL    Calculate 'vorticity1'.
C ----------------------------------------------------------------------

C Calculate d(rs)/d(phi).
        do 110 i=START_POINT_NO_HALO,END_P_POINT_NO_HALO
          drs_dlatitude(i) = latitude_step_inverse*.5*
     &                       (rs_on_theta(i-row_length)-
     &                        rs_on_theta(i+row_length))
 110    continue

C Calculate vorticity1.
        do 120 i=START_POINT_NO_HALO,END_P_POINT_NO_HALO
          vorticity1(i) = -2.*omega/(sec_p_latitude(i)*
     &                               rs_on_theta(i))*drs_dlatitude(i)
 120    continue

C Calculate vorticity1 at poles by summing d(rs)/d(lat) around polar
C circle and averaging.
        scalar = .5*latitude_step_inverse/GLOBAL_ROW_LENGTH
        
        if (at_top_of_LPG) then
          sum_n=0.0
          
          CALL GCG_RVECSUMR(
     &      ROW_LENGTH-2*EW_Halo,ROW_LENGTH-2*EW_Halo,1,1,
     &      rs_on_theta(TOP_ROW_START+FIRST_ROW_PT-1+ROW_LENGTH),
     &      GC_ROW_GROUP,info,sum_n)
          sum_n=-sum_n*scalar
          
          do i=TOP_ROW_START,TOP_ROW_START+ROW_LENGTH-1
            vorticity1(i) = -sum_n*2.0*omega/(sec_p_latitude(i)*
     &                      rs_on_theta(i))
          enddo
        endif

        if (at_base_of_LPG) then
          sum_s=0.0
          
          CALL GCG_RVECSUMR(
     &      ROW_LENGTH-2*EW_Halo,ROW_LENGTH-2*EW_Halo,1,1,
     &      rs_on_theta(P_BOT_ROW_START-row_length+LAST_ROW_PT-1),
     &      GC_ROW_GROUP,info,sum_s)
          sum_s=sum_s*scalar
          
          do i=P_BOT_ROW_START,P_BOT_ROW_START+ROW_LENGTH-1
            vorticity1(i) = -sum_s*2.0*omega/(sec_p_latitude(i)*
     &                      rs_on_theta(i))
          enddo
        endif

! Set rest of array
      CALL SWAPBOUNDS(vorticity1,ROW_LENGTH,tot_P_ROWS,
     &              EW_Halo,NS_Halo,1)

CL    End of routine vortic1

      return
      end

