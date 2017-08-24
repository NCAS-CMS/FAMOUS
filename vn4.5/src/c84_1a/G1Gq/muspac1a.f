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
CLL  Routine: MULTI_SPATIAL --------------------------------------------
CLL
CLL  Purpose: Control routine for spatial processing when extracting a
CLL           timeseries within STASH.  Calls SPATIAL to extract global
CLL           mean samples from subdomains pointed to by the mother
CLL           STASHlist record using weighting and masking codes from
CLL           the STASHlist record within each subdomain.  The list of
CLL           subdomains is held as part of the STASH control file.
CLL           They may be in terms of gridpoints, or latitude/longitude
CLL           relative to the grid coordinates.  All timeseries samples
CLL           at a given step are appended to the output field.
CLL           On the first timestep it fills the entire output
CLL            vector to missing data (apart from values for the
CLL            first timestep and computes extra data for the time-
CLL            series -- tis prevents time meaning routines
CLL            failing due to uninitialised data.
CLL            however as a result of this the output vector length
CLL            will change from timestep to timestep
CLL
CLL  Author:   S.Tett/T.Johns
CLL
CLL  Tested under compiler:   cft77
CLL  Tested under OS version: UNICOS 5.1
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.1  24/02/93  Correct outstanding problems with timeseries (ST).
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL   3.3  30/03/94  Pass explicit size dimension lenout for declaring
CLL                  fieldout output array.  Author Tim Johns.
CLL   3.3  17/09/93  Correct level-dependent mass-weighting probs (TCJ).
!LL   4.3  7/3/97    Added code for MPP functionality    P.Burton
!LL   4.4  18/06/97  MPP: Fixed call to GLOBAL_TO_LOCAL_RC, where
!LL                  row and column arguments were in wrong
!LL                  order.                                 P.Burton
!LL                  MPP: All PEs must contain real data in fieldout
!LL                                                         P.Burton
!LL   4.4  6/8/97    Corrected top_left_pe calculation  P.Burton
CLL
!LL   4.4  27/11/96  New option mean timeseries. R A Stratton.          
!LL   4.5  06/01/97  Fix for MPP timeseries: Inserted sync. before
!LL                  send/receive for correct NAM/SHMEM operation
!LL                                                       P.Burton
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered: D71
CLL
CLL  Project task: D7
CLL
CLL  External documentation:
CLL    Unified Model Doc Paper C4 - Storage handling and diagnostic
CLL                                 system (STASH)
CLL
C*L  Interface and arguments: ------------------------------------------
C
      SUBROUTINE MULTI_SPATIAL(fieldin,vx,vy,vz,gr,st_grid,lcyclic,
     +     lmasswt,horiz_size,num_vert_levels,
     +     pexner,pstar,delta_ak,delta_bk,
     +     cos_p_latitude,cos_u_latitude,land,
     +     row_length,p_rows,u_rows,p_levels,
     +     fieldout,lenout,amdi,
     +     control,control_size,
     +     stash_series,series_entry_size,no_records,
     +     num_stash_levels,index_lev,level_list,start_ts,extraw,
     +     n_rows_out,n_cols_out,
     +     real_hd,len_realhd,int_hd,len_inthd,
     +     ocean,
     +     icode,errmssg)
C
      IMPLICIT NONE
C
      INTEGER vx,vy,vz          ! IN size of fieldin
      INTEGER gr                ! IN ppxref gridtype code
      INTEGER st_grid           ! IN STASH internal gridtype code
      LOGICAL lcyclic           ! IN TRUE if cyclic EW BCs
      LOGICAL lmasswt           ! IN TRUE if level-dep mass-wts OK
      INTEGER row_length,p_rows,u_rows,p_levels ! IN size parameters
      REAL fieldin(vx*vy*vz)    ! IN data field
      INTEGER horiz_size        ! OUT no. of points in horizontal slice
      INTEGER num_vert_levels   ! OUT no of horizontal slices.
      INTEGER num_stash_levels  ! IN size of vertical levels list.
      INTEGER index_lev(num_stash_levels) ! IN offsets for each horiz fi
      INTEGER level_list(num_stash_levels) ! IN level for each horiz. fi
      REAL
     &    pexner(row_length,p_rows,p_levels+1), ! IN exner pressure
     &    pstar(row_length,p_rows),             ! IN surf pressure
     &    delta_ak(p_levels),                   ! IN hybrid coords
     &    delta_bk(p_levels),                   ! IN hybrid coords
     &    cos_p_latitude(row_length,p_rows),    ! IN p-grid area fn
     &    cos_u_latitude(row_length,u_rows)     ! IN u-grid area fn
      LOGICAL land(row_length,p_rows)           ! IN land mask
      LOGICAL ocean                             ! IN true if ocean
      INTEGER lenout               ! IN max size of output field
      REAL fieldout(lenout)        ! OUT output field
      REAL amdi                    ! IN missing data indicator
      INTEGER len_realhd           ! IN size of real header
      INTEGER len_inthd            ! IN size of integer header
      REAL real_hd(len_realhd)     ! IN real header
      INTEGER int_hd(len_inthd)    ! IN integer header
      INTEGER control_size         ! IN size of control array
      INTEGER control(control_size)! IN control array (mostly not used)
      INTEGER series_entry_size    ! IN no of entries in each record.
      INTEGER no_records           ! IN no of records to process.
      INTEGER extraw               ! OUT no of words required by extra d
      INTEGER n_rows_out,n_cols_out! OUT data-set size and extent
      LOGICAL start_ts             ! IN true if first time-series timest
      INTEGER stash_series(series_entry_size,no_records) ! IN
C IN control data for calls to spatial
      INTEGER icode                       ! OUT error code
      CHARACTER*(80) errmssg              ! OUT error message
C*----------------------------------------------------------------------
C Parameters
C
CLL  Comdeck: STPARAM --------------------------------------------------
CLL
CLL  Purpose: Meaningful PARAMETER names for STASH processing routines.
CLL           Both a long name and short name have been declared, to
CLL           reduce the existence of "magic" numbers in STASH.
CLL           Format is that first the address of the item is declare in
CLL           both long and short form. example is;
CLL             integer st_item_code,s_item  !Item number (declaration)
CLL             parameter(st_item_code=3,s_item=3)
CLL
CLL  Author:   S.Tett             Date:           22 January 1991
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.5    Mar. 95  Sub-models project.
CLL                   st_model_code=28 added to STLIST addresses
CLL                                   S.J.Swarbrick
!LL   4.2    27/11/96 MPP code: Added new stlist "magic numbers" :
!LL                   st_dump_output_length, st_dump_output_addr
!LL                                                       P.Burton
!LL   4.4    23/09/97 Add st_offset_code to the STASH list
!LL                   S.D. Mullerworth
!    4.4  02/12/96 Time mean timeseries added R A Stratton.             
!    4.5  23/01/98 Added new stlist magic number
!                  st_dump_level_output_length
CLL
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered: D70
CLL
CLL  Project task: D7
CLL
CLL  External documentation:
CLL    Unified Model Doc Paper C4 - Storage handling and diagnostic
CLL                                 system (STASH)
CLLEND--------------------------------------------------------------
C
         integer st_model_code,s_modl ! Internal model number address
         parameter(st_model_code=28,s_modl=28)

         integer st_sect_no_code,s_sect ! Section Number address
     &          ,st_sect_code
         parameter(st_sect_no_code=2,s_sect=2,st_sect_code=2)

         integer st_item_code,s_item  !Item number address
         parameter(st_item_code=1,s_item=1)

         integer st_proc_no_code,s_proc ! Processing Code address
         parameter(st_proc_no_code=3,s_proc=3)

CL subsidiary codes for st_proc_no_code now

         integer st_replace_code
         parameter(st_replace_code=1)

         integer st_accum_code
         parameter(st_accum_code=2)

         integer st_time_mean_code
         parameter(st_time_mean_code=3)

         integer st_time_series_code
         parameter(st_time_series_code=4)

         integer st_max_code
         parameter(st_max_code=5)

         integer st_min_code
         parameter(st_min_code=6)

         integer st_append_traj_code
         parameter(st_append_traj_code=7)

         integer st_time_series_mean                                    
         parameter(st_time_series_mean=8)                               
                                                                        
         integer st_variance_code
         parameter(st_variance_code=9)                                  

         integer st_freq_code,s_freq ! Frequency (Input & output) addres
         parameter(st_freq_code=4,s_freq=4)

         integer st_offset_code,s_offs ! Offset for sampling
         parameter(st_offset_code=30,s_offs=30)

         integer st_start_time_code,s_times ! start timestep address
         parameter(st_start_time_code=5,s_times=5)

         integer st_end_time_code,s_timee ! end timestep address
         parameter(st_end_time_code=6,s_timee=6)

         integer st_period_code,s_period ! period in timesteps address
         parameter(st_period_code=7,s_period=7)

         integer st_infinite_time        ! infinite end/period value
         parameter(st_infinite_time=-1)

         integer st_end_of_list          ! end-of-list marker in times
         parameter(st_end_of_list=-1)

C ---------------------------- grid point stuff
         integer st_gridpoint_code,s_grid ! gridpoint info address
         parameter(st_gridpoint_code=8,s_grid=8)
CL now subsid grid point stuff
         integer stash_null_mask_code,s_nomask ! no masking done
         parameter(stash_null_mask_code=1,s_nomask=1)

         integer stash_land_mask_code,s_lndms ! land mask conds
         parameter(stash_land_mask_code=2,s_lndms=2)

         integer stash_sea_mask_code,s_seams  ! sea mask code
         parameter(stash_sea_mask_code=3,s_seams =3)

CL processing options

         integer block_size ! size of block for gridpoint code
         parameter(block_size=10)

         integer extract_top ! max code for vertical mean subroutine
         integer extract_base ! base codes for vertical mean subroutine
         parameter(extract_base=block_size*0)
         parameter(extract_top=block_size*1)

         integer vert_mean_top ! max code for vertical mean subroutine
         integer vert_mean_base ! base codes for vertical mean subroutin
         parameter(vert_mean_base=block_size*1)
         parameter(vert_mean_top=block_size*2)

         integer zonal_mean_top ! max code for zonal mean subroutine
         integer zonal_mean_base ! base codes for zonal mean subroutine
         parameter(zonal_mean_base=block_size*2)
         parameter(zonal_mean_top=block_size*3)

         integer merid_mean_top ! max code for meridional mean subroutin
         integer merid_mean_base ! base codes for meridional mean subrou
         parameter(merid_mean_base=block_size*3)
         parameter(merid_mean_top=block_size*4)

         integer field_mean_top ! max code for field mean subroutine
         integer field_mean_base ! base codes for field mean subroutine
         parameter(field_mean_base=block_size*4)
         parameter(field_mean_top=block_size*5)

         integer global_mean_top ! max code for global mean subroutine
         integer global_mean_base ! base codes for global mean subroutin
         parameter(global_mean_base=block_size*5)
         parameter(global_mean_top=block_size*6)

CL Weighting

         integer st_weight_code,s_weight ! weighting info address
         parameter(st_weight_code=9,s_weight=9)

         integer stash_weight_null_code,s_noweight ! value of null weigh
         parameter(stash_weight_null_code=0,s_noweight=0)

         integer stash_weight_area_code,s_areaweight ! value of area wei
         parameter(stash_weight_area_code=1,s_areaweight=1)

         integer stash_weight_volume_code,s_volweight
         parameter(stash_weight_volume_code=2,s_volweight=2)

         integer stash_weight_mass_code,s_massweight ! value of mass wei
         parameter(stash_weight_mass_code=3,s_massweight=3)

CL Domain definition

         integer st_north_code,s_north ! northern row address
         parameter(st_north_code=12,s_north=12)

         integer st_south_code,s_south ! southern row address
         parameter(st_south_code=13,s_south =13)

         integer st_west_code,s_west ! western column address
         parameter(st_west_code=14,s_west=14)

         integer st_east_code,s_east ! eastern row address
         parameter(st_east_code=15,s_east =15)

CL Levels

         integer st_input_bottom,s_bottom ! input bottom level address
         parameter(st_input_bottom=10,s_bottom =10)

         integer  st_special_code,s_special ! special code
         parameter(st_special_code=100,s_special=100)

         integer st_input_top,s_top          ! input top level address
         parameter(st_input_top=11,s_top=11)

         integer st_output_bottom,s_outbot   ! output bottom level addre
         parameter(st_output_bottom=21,s_outbot=21)

         integer st_output_top,s_outtop      ! output top level address
         parameter(st_output_top=22,s_outtop=22)

         integer st_model_level_code,s_model
         parameter(st_model_level_code=1,s_model=1)

         integer st_pressure_level_code,s_press ! code for pressure leve
         parameter( st_pressure_level_code=2,s_press=2)

         integer st_height_level_code,s_height ! code for height levels
         parameter(st_height_level_code=3,s_height=3)

         integer st_input_code,s_input               ! input code addres
         parameter(st_input_code=16,s_input=16)

         integer st_input_length,s_length ! input length of diagnostic
         parameter(st_input_length=17,s_length=17)             ! address

         integer st_output_code,s_output ! output code address
         parameter(st_output_code=18,s_output=18)

C Pointer to D1 addressing information
         integer st_position_in_d1,st_d1pos ! Pos of item in D1 for 
         parameter(st_position_in_d1=29,st_d1pos=29) ! relevant submodel

C Output destination options

         integer st_dump,st_secondary
         parameter(st_dump=1,st_secondary=2)

         integer st_output_length,s_outlen ! output length of diagnostic
         parameter(st_output_length=19,s_outlen=19)           ! address
         integer st_dump_output_length,s_doutlen ! output length on
         parameter(st_dump_output_length=32,s_doutlen=32)  ! dump
         integer st_dump_level_output_length,s_dlevoutlen
         parameter(st_dump_level_output_length=33,s_dlevoutlen=33)
! output length of a single level on dump

         integer st_output_addr,s_outadd ! start locn of diag after stas
         parameter(st_output_addr=20,s_outadd=20)       ! output address
         integer st_dump_output_addr,s_doutadd ! output address on
         parameter(st_dump_output_addr=31,s_doutadd=31)  ! dump

         integer st_lookup_ptr       ! ptr to dump lookup header address
         parameter(st_lookup_ptr=23)

         integer st_series_ptr ! ptr into stash_series where control dat
         parameter(st_series_ptr=24)                            ! addres

CL subsid stuff for time series
         integer series_grid_type
         parameter(series_grid_type=1)

         integer series_grid_code
         parameter(series_grid_code=0)

         integer series_long_code
         parameter(series_long_code=1)

         integer series_size
         parameter(series_size=2)

         integer series_proc_code
         parameter(series_proc_code=3)

         integer series_north
         parameter(series_north=4)

         integer series_south
         parameter(series_south=5)

         integer series_west
         parameter(series_west=6)

         integer series_east
         parameter(series_east=7)

         integer series_list_start
         parameter(series_list_start=8)

         integer series_list_end
         parameter(series_list_end=9)

         integer record_size
         parameter(record_size=9)

C Miscellaneous parameters

         integer st_macrotag   ! system/user tag field in stlist address
         parameter(st_macrotag=25)

C Pseudo-level list pointers

         integer st_pseudo_in        ! pseudo-levels input list address
         parameter(st_pseudo_in=26)

         integer st_pseudo_out       ! pseudo-levels output list address
         parameter(st_pseudo_out=27)

C Internal horizontal gridtype codes common to all diagnostics

         integer st_tp_grid,st_uv_grid, ! T-p grid, u-v grid
     &           st_cu_grid,st_cv_grid, ! C-grid (u point, v point)
     &           st_zt_grid,st_zu_grid, ! Zonal T-grid, u-grid
     &           st_mt_grid,st_mu_grid, ! Meridional T-grid, u-grid
     &           st_scalar              ! Scalar (ie. single value)
         parameter(st_tp_grid=1,
     &             st_uv_grid=2,
     &             st_cu_grid=3,
     &             st_cv_grid=4,
     &             st_zt_grid=5,
     &             st_zu_grid=6,
     &             st_mt_grid=7,
     &             st_mu_grid=8,
     &             st_scalar=9)
CLL  Comdeck: STERR ----------------------------------------------------
CLL
CLL  Purpose: PARAMETER names for STASH processing error codes;
CLL           fatal errors have positive codes, warnings negative.
CLL
CLL  Author:   S.Tett
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL   3.3  16/09/93  Add st_illegal_weight error code.
!LL                   Added st_no_data for MPP code
!LL                   (means a processor does not contain any data
!LL                    for a given subdomain)                 P.Burton
CLL
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered: D70
CLL
CLL  Project task: D7
CLL
CLL  External documentation:
CLL    Unified Model Doc Paper C4 - Storage handling and diagnostic
CLL                                 system (STASH)
C
C Warning codes
C
         integer st_upper_less_lower ! warning code for bad domain
         parameter(st_upper_less_lower=-1)

         integer st_not_supported ! warning code for unsupported routine
         parameter(st_not_supported=-2)
         integer st_no_data,st_nd ! indicates no data on a processor
         parameter(st_no_data=-3,st_nd=-3)
C
C Error codes
C
         integer st_bad_array_param ! error code for dodgy array params
         parameter(st_bad_array_param=1)

         integer st_bad_address     ! error code for address violation
         parameter(st_bad_address=2)

         integer st_unknown ! error code for unknown option
         parameter(st_unknown=3)

         integer st_bad_wraparound ! error code for illegal wraparound
         parameter(st_bad_wraparound=4)

         integer st_illegal_weight ! error code for illegal weighting
         parameter(st_illegal_weight=9)

         integer unknown_weight ! error code for an unknown weight
         parameter(unknown_weight=10)

         integer unknown_mask ! error code for an unknown mask
         parameter(unknown_mask=11)

         integer unknown_processing ! error code for unknown processing
         parameter(unknown_processing=12)

         integer nonsense ! error code for general nonsense request
         parameter(nonsense=13)

CDIR$ FIXED
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C GC - General Communication primitives package. For use on
C multiprocessor shared memory and message passing systems.
C
C
C LICENSING TERMS
C
C  GC is provided free of charge. Unless otherwise agreed with SINTEF,
C  use and redistribution in source and binary forms are permitted
C  provided that
C
C      (1) source distributions retain all comments appearing within
C          this file header, and
C
C      (2) distributions including binaries display the following
C          acknowledgement:
C
C              "This product includes software developed by SINTEF.",
C
C          in the documentation or other materials provided with the
C          distribution and in all advertising materials mentioning
C          features or use of this software.
C
C  The name of SINTEF may not be used to endorse or promote products
C  derived from this software without specific prior written
C  permission.  SINTEF disclaims any warranty that this software will
C  be fit for any specific purposes. In no event shall SINTEF be liable
C  for any loss of performance or for indirect or consequential damage
C  or direct or indirect injury of any kind. In no case shall SINTEF
C  be liable for any representation or warranty make to any third party
C  by the users of this software.
C
C
C Fortran header file. PLEASE use the parameter variables in user
C routines calling GC and NOT the numeric values. The latter are
C subject to change without further notice.
C
C---------------------------------------------- ------------------------
C $Id: gpb2f402,v 1.10 1996/11/28 20:36:24 t11pb Exp $
C (C) Jorn Amundsen, Roar Skaalin, SINTEF Industrial Mathematics.

C    4.4   30/09/97  Added code to permit the SHMEM/NAM timeout
C                    value to be set from a shell variable.
C                      Author: Bob Carruthers  Cray Research.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C     GC general options
      INTEGER GC_OK, GC_FAIL, GC_NONE, GC_ANY, GC_DONTCARE,
     $     GC_SHM_DIR, GC_SHM_GET, GC_SHM_PUT, GC_USE_GET, GC_USE_PUT
     &   , GC_NAM_TIMEOUT, GC_SHM_SAFE
      PARAMETER (GC_OK         =     0)
      PARAMETER (GC_FAIL       =    -1)
      PARAMETER (GC_NONE       =     0)
      PARAMETER (GC_ANY        =    -1)
      PARAMETER (GC_DONTCARE   =    -1)
      PARAMETER (GC_SHM_DIR    =     1)
      PARAMETER (GC_SHM_SAFE   =     2)
      PARAMETER (GC_NAM_TIMEOUT=     4)
      PARAMETER (GC_SHM_GET    = -9999)
      PARAMETER (GC_SHM_PUT    = -9998)
      PARAMETER (GC_USE_GET    = -9999)
      PARAMETER (GC_USE_PUT    = -9998)

C     GC functions
      INTEGER GC_COMLEN, GC_ISIZE, GC_RSIZE, GC_ME, GC_NPROC

C     GC groups (GCG) support
      INTEGER GC_ALLGROUP, GCG_ALL
      PARAMETER (GC_ALLGROUP = 0)
      PARAMETER (GCG_ALL = GC_ALLGROUP)

C     GC groups (GCG) functions
      INTEGER GCG_ME

C     GC reserved message tags
      INTEGER GC_MTAG_LOW, GC_MTAG_HIGH
      PARAMETER (GC_MTAG_LOW   = 999999901)
      PARAMETER (GC_MTAG_HIGH  = 999999999)

C     GCG_RALLETOALLE index parameters
      INTEGER S_DESTINATION_PE, S_BASE_ADDRESS_IN_SEND_ARRAY,
     $     S_NUMBER_OF_ELEMENTS_IN_ITEM, S_STRIDE_IN_SEND_ARRAY,
     $     S_ELEMENT_LENGTH, S_BASE_ADDRESS_IN_RECV_ARRAY,
     $     S_STRIDE_IN_RECV_ARRAY
      PARAMETER (S_DESTINATION_PE = 1)
      PARAMETER (S_BASE_ADDRESS_IN_SEND_ARRAY = 2)
      PARAMETER (S_NUMBER_OF_ELEMENTS_IN_ITEM = 3)
      PARAMETER (S_STRIDE_IN_SEND_ARRAY = 4)
      PARAMETER (S_ELEMENT_LENGTH = 5)
      PARAMETER (S_BASE_ADDRESS_IN_RECV_ARRAY = 6)
      PARAMETER (S_STRIDE_IN_RECV_ARRAY = 7)

      INTEGER R_SOURCE_PE, R_BASE_ADDRESS_IN_RECV_ARRAY,
     $     R_NUMBER_OF_ELEMENTS_IN_ITEM, R_STRIDE_IN_RECV_ARRAY,
     $     R_ELEMENT_LENGTH, R_BASE_ADDRESS_IN_SEND_ARRAY,
     $     R_STRIDE_IN_SEND_ARRAY
      PARAMETER (R_SOURCE_PE = 1)
      PARAMETER (R_BASE_ADDRESS_IN_RECV_ARRAY = 2)
      PARAMETER (R_NUMBER_OF_ELEMENTS_IN_ITEM = 3)
      PARAMETER (R_STRIDE_IN_RECV_ARRAY = 4)
      PARAMETER (R_ELEMENT_LENGTH = 5)
      PARAMETER (R_BASE_ADDRESS_IN_SEND_ARRAY = 6)
      PARAMETER (R_STRIDE_IN_SEND_ARRAY = 7)
! ------------------------ Comdeck PARVARS -------------------------
! Parameters and common blocks required by the MPP-UM
!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the MPP-UM
!
!   Two sets of parameters are set up -
!     i)  for the MPP-UM itself.
!     ii) for the interface to the Message Passing Software.
!
!   History:
!
!   Model    Date     Modification history
!  version
!   4.1      27/1/96  New comdeck based on first section of
!                     old PARVARS.   P.Burton
!   4.2      21/11/96 Add new field type parameter and
!                     magic number used in addressing to indicate
!                     if a calculation is for local data, or data
!                     on the dump on disk (ie. global data)  P.Burton
!   4.2      18/11/96 Moved MaxFieldSize to comdeck AMAXSIZE and
!                     removed Maxbuf.  P.Burton
!   4.2      18/7/96  Removed some unused variables      P.Burton
!   4.4      11/07/97 Reduced MAXPROC to 256 to save memory  P.Burton
!
! ---------------------- PARAMETERS ---------------------
!
! =======================================================
! Parameters needed for the MPP-UM
! =======================================================

      INTEGER   Ndim_max        ! maximum number of spatial dimensions
      PARAMETER (Ndim_max = 3 ) ! 3d data


      INTEGER
     &   fld_type_p           ! indicates a grid on P points
     &,  fld_type_u           ! indicates a grid on U points
     &,  fld_type_unknown     ! indicates a non-standard grid.
      PARAMETER (
     &   fld_type_p=1
     &,  fld_type_u=2
     &,  fld_type_unknown=-1)

      INTEGER
     &   local_data
     &,  global_dump_data
      PARAMETER (
     &   local_data=1        ! Used in addressing to indicate if
     &,  global_dump_data=2) ! calculation is for a local or
!                            ! global (ie. disk dump) size

! =======================================================
! Parameters needed for the Message Passing Software
! =======================================================


      INTEGER
     &   Maxproc              ! Max number of processors
      PARAMETER (
     &   MAXPROC = 256)

      INTEGER
     &   PNorth       ! North processor address in the neighbour array
     &,  PEast        ! East  processor address in the neighbour array
     &,  PSouth       ! South processor address in the neighbour array
     &,  PWest        ! West  processor address in the neighbour array
     &,  NoDomain     ! Value in neighbour array if the domain has
     &                !  no neighbor in this direction. Otherwise
     &                !  the value will be the tid of the neighbor
      PARAMETER (
     &   PNorth   = 1
     &,  PEast    = 2
     &,  PSouth   = 3
     &,  PWest    = 4
     &,  NoDomain = -1)

      INTEGER
     &   BC_STATIC            ! Static boundary conditions
     &,  BC_CYCLIC            ! Cyclic boundary conditions
      PARAMETER (
     &   BC_STATIC = 1
     &,  BC_CYCLIC = 2)

! ---------------------- End of comdeck PARPARM ---------------------
!========================== COMDECK PARCOMM ====================
!
! *** NOTE : This comdeck requires comdeck PARPARM to be *CALLed
!            first.
!
!   Description:
!
!   This COMDECK contains COMMON blocks for the MPP-UM
!
!
!   Two COMMON blocks are defined:
!     i)  UM_PARVAR holds information required by the
!         Parallel Unified Model itself
!     ii) MP_PARVAR holds information required by the interface to
!         the Message Passing Software used by the PUM
!
!   Key concepts used in the inline documentation are:
!     o GLOBAL data - the entire data domain processed by the UM
!     o LOCAL data - the fragment of the GLOBAL data which is
!       stored by this particular process
!     o PERSONAL data - the fragment of the LOCAL data which is
!       updated by this particular process
!     o HALO data - a halo around the PERSONAL data which forms
!       the LOCAL data
!
!     Acronyms used:
!     LPG - Logical Process Grid, this is the grid of logical
!           processors; each logical processor handles one of the
!           decomposed parts of the global data. It does not
!           necessarily represent a physical grid of processors.
!
!   History:
!
!   4.1      27/1/96  New comdeck based on second section of
!                     old PARVARS.   P.Burton
!   4.2     19/08/96  Removed some unused variables, and added
!                     current_decomp_type variable to allow use
!                     of flexible decompositions.
!                     Added nproc_max to indicate the max. number
!                     of processors used for MPP-UM
!                                                      P.Burton
!
! -------------------- COMMON BLOCKS --------------------
!
! =======================================================
! Common block for the Parallel Unified Model
! =======================================================

      INTEGER
     &   first_comp_pe       ! top left pe in LPG
     &,  last_comp_pe        ! bottom right pe in LPG
     &,  current_decomp_type ! current decomposition type
     &,  Offx                ! halo size in EW direction
     &,  Offy                ! halo size in NS direction
     &,  glsize(Ndim_max)    ! global data size
     &,  lasize(Ndim_max)    ! local data size
     &,  blsizep(Ndim_max)   ! personal p data area
     &,  blsizeu(Ndim_max)   ! personal u data area
     &,  datastart(Ndim_max) ! position of personal data in global data
     &                       !   (in terms of standard Fortran array
     &                       !    notation)
     &,  gridsize(Ndim_max)  ! size of the LPG in each dimension
     &,  gridpos(Ndim_max)   ! position of this process in the LPG
!                            ! 0,1,2,...,nproc_x-1 etc.

      LOGICAL
     &    atbase             ! process at the bottom of the LPG
     &,   attop              ! process at the top of the LPG
     &,   atleft             ! process at the left of the LPG
     &,   atright            ! process at the right of the LPG
! NB: None of the above logicals are mutually exclusive

      COMMON /UM_PARVAR/
     &                  first_comp_pe,last_comp_pe
     &,                 current_decomp_type,Offx, Offy
     &,                 glsize,lasize,blsizep,blsizeu
     &,                 datastart,gridsize,gridpos
     &,                 atbase,attop,atleft,atright

! =======================================================
! Common block for the Message Passing Software
! =======================================================

      INTEGER
     &  bound(Ndim_max)           ! type of boundary (cyclic or static)
     &                            !  in each direction
     &, g_lasize(Ndim_max,0:maxproc)
!                                 ! global copy of local data size
     &, g_blsizep(Ndim_max,0:maxproc)
!                                 ! global copy of personal p data area
     &, g_blsizeu(Ndim_max,0:maxproc)
!                                 ! global copy of personal u data area
     &, g_datastart(Ndim_max,0:maxproc)
!                                 ! global copy of datastart
     &, g_gridpos(Ndim_max,0:maxproc)
!                                 ! global copy of gridpos
     &, nproc                     ! number of processors in current
!                                 ! decomposition
     &, nproc_max                 ! maximum number of processors
     &, nproc_x                   ! number of processors in x-direction
     &, nproc_y                   ! number of processors in y-direction
     &, mype                      ! number of this processor
     &                            !  (starting from 0)
     &, neighbour(4)              ! array with the tids of the four
     &                            ! neighbours in the horizontal plane
     &, gc_proc_row_group         ! GID for procs along a proc row
     &, gc_proc_col_group         ! GID for procs along a proc col
     &, gc_all_proc_group         ! GID for all procs

      COMMON /MP_PARVAR/
     &                  bound
     &,                 g_lasize,g_blsizep,g_blsizeu
     &,                 g_datastart,g_gridpos
     &,                 nproc,nproc_max,nproc_x,nproc_y,mype
     &,                 neighbour,gc_proc_row_group
     &,                 gc_proc_col_group, gc_all_proc_group



! ---------------------- End of comdeck PARCOMM -----------------------
! --------------------- End of comdeck PARVARS ---------------------
C*L
C Subroutines called
C
      EXTERNAL SPATIAL
      EXTERNAL STASH_COMP_GRID,EXTRA_MAKE_VECTOR,EXTRA_TS_INFO
C*
C Local variables
C
      INTEGER fake_record(control_size) ! fake_record for SPATIAL call
      INTEGER pp_ptr ! ptr to pp_field for where output from spatial goe
      INTEGER i,j                       ! loop variable
      INTEGER data_size                 ! size of data.
      INTEGER index_size
      INTEGER base_level        !   reference model level
      INTEGER kl ! loop count for levels
      INTEGER stash_list_start ! the start address in index_levs for lev
      INTEGER stash_list_end ! the end address in index_levs for levels
      INTEGER what_process ! what kind of processing is requested
      INTEGER what_mask ! what mask is wanted.
      INTEGER extra_start ! start address in fieldout for extra data
      REAL bdx,bzx,bdy,bzy ! grid descriptors

      INTEGER
     &  proc_top_left_x ! processor co-ords of processor at top left
     &, proc_top_left_y ! corner of subarea
     &, top_left_pe     ! processor id of this processor
     &, dummy1,dummy2   ! unused return variables from
!                         GLOBAL_TO_LOCAL_RC
     &, info            ! GCOM return variable

      REAL
     &  global_mean ! global mean returned by call to SPATIAL
      COMMON /MPP_STATIC_VAR/ global_mean  !must be memory aligned
!                                           to allow fast shmem


CL----------------------------------------------------------------------
CL 1. Error checking
CL
C  Check we are in fact doing a time series. Error if not.
      IF ((control(st_proc_no_code).ne.st_time_series_code).and.        
     &   (control(st_proc_no_code).ne.st_time_series_mean)) THEN 
        icode=st_unknown
        write(errmssg,99)control(st_proc_no_code)
99      format(3x,'MULTI_SP : unexpected call to extract timeseries',i5)
        goto 999            ! jump to return
      ENDIF
CL----------------------------------------------------------------------
CL 2. Workout what kind of processing we are doing and what mask is used
CL
       what_mask=mod(control(st_gridpoint_code),block_size)
       what_process=control(st_gridpoint_code)-what_mask
       extraw=0
CL note that the first word in fieldout is assumed to be
CL the word where the first spatial domain mean for this timeseries
CL  will be stored
CL
CL Next we compute the grid discriptors -- as used by extra data
       IF (start_ts) THEN ! compute grid descrip
         extraw=6*(no_records+1)
         extra_start=control(st_output_length)-extraw+1
         CALL STASH_COMP_GRID(bzx,bzy,bdx,bdy,0,st_grid,ocean,
     &     1,1,real_hd,len_realhd,int_hd,len_inthd,extract_base+1,
     &     icode,errmssg)
       ENDIF
CL 3. Set up pseudo STASH record to be passed to SPATIAL on each call
CL    to extract a sample from the input field.
CL
      fake_record(st_input_bottom)=control(st_input_bottom)
      fake_record(st_input_top)=control(st_input_top)
      fake_record(st_weight_code)=control(st_weight_code)
C doing a field mean on this sub-domain with mask specified by control.
      pp_ptr=1
CL----------------------------------------------------------------------
CL 4. Loop over samples and extract global mean within subdomain for
CL    each, appending to output field
CL
      DO i=1,no_records ! loop over sub domains
        data_size=stash_series(series_size,i)
CL 4.1 Do preliminary verifications on stash_series
CL 4.1.1 Gridtype code
        IF (stash_series(series_grid_type,i).ne.series_grid_code) THEN
C
C NB: Latitude/longitude range conversion to gridpoint range needs to
C     be added
C
          icode=st_not_supported
          errmssg='MULTI_SP : only support grid point processing'
          goto 999
        ENDIF
CL ---------------------------------------------------------------------
CL 5. Set up the fake record domain info depending on what kind of
CL    "primary" processing is requested.
CL    As far as stash is concerned everything looks like a global mean
CL    here and it is just a question of setting up the fake record
CL    correctly.
CL
        fake_record(st_gridpoint_code)=
     &    (stash_series(series_proc_code,i)/block_size)*block_size
     &      +what_mask
        IF (what_process.eq.extract_base) THEN ! an extract
          fake_record(st_north_code)=stash_series(series_north,i)
          fake_record(st_south_code)=stash_series(series_south,i)
          fake_record(st_west_code)= stash_series(series_west,i)
          fake_record(st_east_code)= stash_series(series_east,i)
          stash_list_start=stash_series(series_list_start,i)
          stash_list_end=stash_series(series_list_end,i)
          fake_record(st_input_bottom)=stash_list_start
          fake_record(st_input_top)=stash_list_end
        ELSEIF (what_process.eq.zonal_mean_base) THEN ! a zonal_mean
          fake_record(st_north_code)=stash_series(series_north,i)
          fake_record(st_south_code)=stash_series(series_south,i)
          fake_record(st_west_code)= control(st_west_code)
          fake_record(st_east_code)= control(st_east_code)
          stash_list_start=stash_series(series_list_start,i)
          stash_list_end=stash_series(series_list_end,i)
          fake_record(st_input_bottom)=stash_list_start
          fake_record(st_input_top)=stash_list_end
        ELSEIF (what_process.eq.merid_mean_base) THEN ! a merid_mean
          fake_record(st_north_code)= control(st_north_code)
          fake_record(st_south_code)= control(st_south_code)
          fake_record(st_east_code)=stash_series(series_east,i)
          fake_record(st_west_code)=stash_series(series_west,i)
          stash_list_start=stash_series(series_list_start,i)
          stash_list_end=stash_series(series_list_end,i)
          fake_record(st_input_bottom)=stash_list_start
          fake_record(st_input_top)=stash_list_end
        ELSEIF (what_process.eq.vert_mean_base) THEN ! a vert_mean
          fake_record(st_north_code)=stash_series(series_north,i)
          fake_record(st_south_code)=stash_series(series_south,i)
          fake_record(st_east_code)=stash_series(series_east,i)
          fake_record(st_west_code)=stash_series(series_west,i)
          stash_list_start=1
          stash_list_end=num_stash_levels
          fake_record(st_input_bottom)=stash_list_start
          fake_record(st_input_top)=stash_list_end
        ELSEIF (what_process.eq.field_mean_base) THEN ! a field_mean
          fake_record(st_north_code)=control(st_north_code)
          fake_record(st_south_code)=control(st_south_code)
          fake_record(st_east_code)=control(st_east_code)
          fake_record(st_west_code)=control(st_west_code)
          stash_list_start=1
          stash_list_end=num_stash_levels
          fake_record(st_input_bottom)=stash_list_start
          fake_record(st_input_top)=stash_list_end
        ELSE ! error code...
          icode=unknown_processing
          write(errmssg,111) 'unknown processing option',what_process
          goto 999 ! jump to error return
        ENDIF
CL Check record (south > north and west < east)
        IF (fake_record(st_north_code).gt.
     +    fake_record(st_south_code))then
          write(errmssg,101)fake_record(st_north_code),
     +       fake_record(st_south_code),i
          icode=st_bad_array_param
          goto 999 ! error exit
        ENDIF
        IF (fake_record(st_west_code).gt.
     +    fake_record(st_east_code))then
          write(errmssg,102)fake_record(st_west_code),
     +       fake_record(st_east_code),i
          icode=st_bad_array_param
          goto 999 ! error exit
        ENDIF

! Figure out which processor is at the top-left of the subdomain
! that SPATIAL will process. This is the one that will send the
! global sum to PE 0 for storage

        CALL GLOBAL_TO_LOCAL_RC(gr,
     &    fake_record(st_west_code),fake_record(st_north_code),
     &    proc_top_left_x, proc_top_left_y,
     &    dummy1,dummy2)

        top_left_pe=proc_top_left_x+proc_top_left_y*nproc_x

C
C NB: At present timeseries samples are global (ie. 3D) means, so
C     there is no levels loop outside the call to SPATIAL here -
C     this may be extended at some point to allow multi-level
C     timeseries sampling inside a levels loop
C
C     n_cols_out and n_rows_out are recalculated within SPATIAL but are
C     now appropriate for an individual timeseries sample, not the whole
C     field.  They are reset for the whole field after subdomain loop.
C
        lcyclic=.false.
        base_level=control(st_input_bottom)+index_lev(1)-1
        CALL SPATIAL(fieldin,vx,vy,vz,gr,st_grid,lcyclic,lmasswt,
     +       n_cols_out,n_rows_out,base_level,
     +       level_list(stash_list_start),
     +       index_lev(stash_list_start),
     +       (stash_list_end+1-stash_list_start),
     +       pexner,pstar,delta_ak,delta_bk,
     +       cos_p_latitude,cos_u_latitude,land,
     +       row_length,p_rows,u_rows,p_levels,
     +       global_mean,1,
     +       fake_record,control_size,amdi,
     +       icode,errmssg)
        IF (icode.ne.0) goto 999 ! got some error so jump to return

! Must move the global_mean data to PE 0 which stores all timeseries
! data
! (NB. This assumes that the output from SPATIAL is just a
!      single number)

        CALL GC_SSYNC(nproc,info)


        IF (mype .EQ. top_left_pe) THEN
          CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_GET,info)
          info=GC_NONE
          CALL GC_RSEND(100,1,0,info,fieldout(pp_ptr),global_mean)
        ENDIF

        CALL GC_SSYNC(nproc,info)

        IF (mype .EQ. 0) THEN
          CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_GET,info)
          info=GC_NONE
          CALL GC_RRECV(100,1,top_left_pe,info,
     &                  fieldout(pp_ptr),global_mean)
        ELSE
          fieldout(pp_ptr)=0.0
        ENDIF


        CALL GC_SSYNC(nproc,info)
        pp_ptr=pp_ptr+(n_cols_out*n_rows_out)  ! increment the pp_ptr
C
C NB: n_cols_out and n_rows_out should both be 1 as timeseries samples
C       are currently designed to be scalar quantities only.
C
CL check on n_cols_out and n_rows_out
        IF (n_cols_out.ne.1) THEN
          errmssg='MULTI_SP : n_cols_out <> 1'
          icode=st_not_supported
          goto 999
        ENDIF
        IF (n_rows_out.ne.1) THEN
          errmssg='MULTI_SP : n_rows_out <> 1'
          icode=st_not_supported
          goto 999
        ENDIF
        IF (mype .EQ. 0) THEN
        IF (start_ts) THEN ! put the descriptive info for this record
          CALL EXTRA_MAKE_VECTOR(fake_record,control_size,
     &      i,no_records,fieldout(extra_start),extraw,bzx,bzy,bdx,bdy)
        ENDIF
        ENDIF
      ENDDO   ! end the loop over sub-domains
C
      horiz_size=pp_ptr-1
      num_vert_levels=1
CL --------------------------------------------------------------------
CL 7. If this is the first time in a time-series then
CL     put the codes describing the extra data into the extra data fld
CL     In addition set pphoriz out to the total length
CL       as well as setting the input vetor to missing
CL       where no values are set
CL----------------------------------------------------------------------

       n_cols_out=no_records
       n_rows_out=control(st_period_code)/control(st_freq_code)
       horiz_size=n_cols_out
       IF (start_ts) THEN  ! on start timestep we have entire vector
         horiz_size=n_cols_out*n_rows_out+extraw
        IF (mype .EQ. 0) THEN
         CALL EXTRA_TS_INFO(fieldout(extra_start),extraw,no_records)
         do i=no_records+1,extra_start-1
           fieldout(i)=amdi
         enddo
        ELSE
          DO i=no_records+1,lenout
            fieldout(i)=0.0
          ENDDO
        ENDIF
       ENDIF
C
999   CONTINUE ! jump here for error exit
C
111   FORMAT('MULTI_SP :  >>>FATAL ERROR <<',a40,i5,i5)
101   FORMAT('MULTI_SP : north > south',2i5,' in record ',i5)
102   FORMAT('MULTI_SP : west > east',2i5,'in record ',i5 )
      RETURN
      END
