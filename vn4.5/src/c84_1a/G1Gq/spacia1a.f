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
CLL  Routine: SPATIAL --------------------------------------------------
CLL
CLL  Purpose: Performs general spatial processing on an input field to
CLL           produce an output field or scalar within STASH.  Lower-
CLL           level routines are called to perform the various spatial
CLL           processing options.
CLL
CLL  Author:   T.Johns/S.Tett
CLL
CLL  Tested under compiler:   cft77
CLL  Tested under OS version: UNICOS 5.1
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.3  30/03/94  Explicitly declare (sub-addressed) output field
CLL                  fieldout using 'lenout' dimension.  Tim Johns.
CLL   3.3  16/09/93  Pass LOGICAL lmasswt to processing routines to
CLL                  denote that level-by-level mass-weights exist.
!LL   4.3   9/12/96  Added MPP code.
!LL                  Moved calculation of weighting and masking terms
!LL                  up from processing routines.            P.Burton
!LL   4.4   13/06/97 MPP : Where reduction spatial meaning takes place
!LL                  processors not getting results should set
!LL                  their diagnostic space to zeros.          P.Burton
!LL   4.4   22/10/97 MPP : Prevent uninitialised points when
!LL                  pstar_weight on U or C grid S.D.Mullerworth
CLL
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
CLL-----------------------------------------------------------------
C*L  Interface and arguments: ------------------------------------------
      SUBROUTINE SPATIAL(fieldin,vx,vy,vz,GR,st_grid,lcyclic,lmasswt,
     +      n_cols_out,n_rows_out,
     +      base_level,level_list,index_lev,no_of_levels,
     +      pexner,pstar,delta_ak,delta_bk,
     +      cos_p_latitude,cos_u_latitude,land,
     +      row_length,p_rows,u_rows,p_levels,
     +      fieldout,lenout,
     +      control,control_size,rmdi,
     +      icode,cmessage)
C
      IMPLICIT NONE
C
      INTEGER
     &    vx,vy,vz,                         ! IN size of fieldin
     &    lenout,                           ! IN size of fieldout
     &    GR,                               ! IN ppxref gridtype code
     &    st_grid,                          ! IN STASH gridtype code
     &    n_rows_out,                       ! OUT no. of output rows
     &    n_cols_out,                       ! OUT no. of output cols
     &    base_level,                       ! IN reference model level
     &    row_length,p_rows,u_rows,p_levels,! IN size parameters
     &    control_size,                     ! IN size of control record
     &    control(control_size),            ! IN control record
     &    icode,                            ! OUT error code 0 if ok
     &    no_of_levels,                     ! IN no of levels
     &    index_lev(no_of_levels),          ! IN index to levels
     &    level_list(no_of_levels)          ! IN model level list
      REAL
     &    fieldin(vx,vy,vz), ! IN fieldin which is to be acted on
     &    pexner(row_length,p_rows,p_levels+1), ! IN exner pressure
     &    pstar(row_length,p_rows),             ! IN surf pressure
     &    delta_ak(p_levels),                   ! IN hybrid coords
     &    delta_bk(p_levels),                   ! IN hybrid coords
     &    cos_p_latitude(row_length,p_rows),    ! IN p-grid area fn
     &    cos_u_latitude(row_length,u_rows),    ! IN u-grid area fn
     &    fieldout(lenout),                     ! OUT output field
     &    rmdi                                  ! IN  missing data indic
      LOGICAL
     &    lcyclic,                              ! IN .true. if cyclic EW
     &    lmasswt,                              ! IN  TRUE if masswts OK
     &    land(row_length,p_rows)               ! IN land mask
      CHARACTER*(*) cmessage                    ! OUT error message

C*----------------------------------------------------------------------
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
CLL  Comdeck: CPPXREF --------------------------------------------------
CLL
CLL  Purpose: Holds PARAMETERs describing structure of PP_XREF file,
CLL           and some values for valid entries.
CLL
CLL  Author    Dr T Johns
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
CLL                  1.Removes the limit on primary STASH item numbers.
CLL                  2.Removes the assumption that (section,item)
CLL                    defines the sub-model.
CLL                  3.Thus allows for user-prognostics.
CLL                  Add a PPXREF record for model number.
CLL  4.0   26/07/95  T.Johns.  Add codes for real/int/log data types.
CLL  3.5   10/3/94   Sub-Models project:
CLL                 List of PPXREF addressing codes augmented, in order
CLL                 to include all of the pre_STASH master information
CLL                 in the new PPXREF file.
CLL                 PPXREF_CODELEN increased to 38.
CLL                 PPXREF_IDLEN deleted - no longer relevant.
CLL                   S.J.Swarbrick
CLL  4.1   June 96  Wave model parameters included.
CLL                 ppx_ address parameters adjusted to allow for 
CLL                  reading option code as 4x5 digit groups.
CLL                   S.J.Swarbrick  
CLL
CLL  Logical components covered: C40
CLL
C-----------------------------------------------------------------------
C Primary file record definition
      INTEGER
     *       PPXREF_IDLEN,PPXREF_CHARLEN,PPXREF_CODELEN
     *      ,PPXREF_PACK_PROFS
      PARAMETER(
     *       PPXREF_IDLEN=2,               ! length of id in a record
! WARNING: PPXREF_CHARLEN must be an exact multiple of 4
!                         to avoid overwriting
     *       PPXREF_CHARLEN=36,            ! total length of characters
     *       PPXREF_PACK_PROFS=10,         ! number of packing profiles
     *       PPXREF_CODELEN=40)            ! total length of codes      
C Derived file record sizes
      INTEGER
     *       PPX_CHARWORD,PPX_RECORDLEN
      PARAMETER(
C            Assume that an integer is at least 4 bytes long.
C            This wastes some space and memory on 8 byte machines.
     *       PPX_CHARWORD=((PPXREF_CHARLEN+3)/4), ! i.e., ppx_charword=9
     *       PPX_RECORDLEN=
     *           PPX_CHARWORD+PPXREF_CODELEN)  ! read buffer record len
C
C-----------------------------------------------------------------------
C Addressing codes within PPXREF
      INTEGER
     &       ppx_model_number   ,
     &       ppx_section_number ,ppx_item_number    ,
     &       ppx_version_mask   ,ppx_space_code     ,
     &       ppx_timavail_code  ,ppx_grid_type      ,
     &       ppx_lv_code        ,ppx_lb_code        ,
     &       ppx_lt_code        ,ppx_lev_flag       ,
     &       ppx_opt_code       ,ppx_pt_code        ,
     &       ppx_pf_code        ,ppx_pl_code        ,
     &       ppx_ptr_code       ,ppx_lbvc_code      ,
     &       ppx_dump_packing   ,ppx_rotate_code    ,
     &       ppx_field_code     ,ppx_user_code      ,
     &       ppx_meto8_levelcode,ppx_meto8_fieldcode,
     &       ppx_cf_levelcode   ,ppx_cf_fieldcode   ,
     &       ppx_base_level     ,ppx_top_level      ,
     &       ppx_ref_LBVC_code  ,ppx_data_type      ,
     &       ppx_packing_acc    ,ppx_pack_acc
      PARAMETER(
     &       ppx_model_number   = 1,  ! Model number address
     &       ppx_section_number = 2,  ! Section number address
     &       ppx_item_number    = 3,  ! Item number address
     &       ppx_version_mask   = 4,  ! Version mask address
     &       ppx_space_code     = 5,  ! Space code address
     &       ppx_timavail_code  = 6,  ! Time availability code address
     &       ppx_grid_type      = 7,  ! Grid type code address
     &       ppx_lv_code        = 8,  ! Level type code address
     &       ppx_lb_code        = 9,  ! First level code address
     &       ppx_lt_code        =10,  ! Last level code address
     &       ppx_lev_flag       =11,  ! Level compression flag address
     &       ppx_opt_code       =12,  ! Sectional option code address
     &       ppx_pt_code        =16,  ! Pseudo dimension type address   
     &       ppx_pf_code        =17,  ! First pseudo dim code address   
     &       ppx_pl_code        =18,  ! Last pseudo dim code address    
     &       ppx_ptr_code       =19,  ! Section 0 point-back code addres
     &       ppx_dump_packing   =20,  ! Dump packing code address       
     &       ppx_lbvc_code      =21,  ! PP LBVC code address            
     &       ppx_rotate_code    =22,  ! Rotation code address           
     &       ppx_field_code     =23,  ! PP field code address           
     &       ppx_user_code      =24,  ! User code address               
     &       ppx_meto8_levelcode=25,  ! CF level code address           
     &       ppx_meto8_fieldcode=26,  ! CF field code address           
     &       ppx_cf_levelcode   =25,                                    
     &       ppx_cf_fieldcode   =26,                                    
     &       ppx_base_level     =27,  ! Base level code address         
     &       ppx_top_level      =28,  ! Top level code address          
     &       ppx_ref_lbvc_code  =29,  ! Ref level LBVC code address     
     &       ppx_data_type      =30,  ! Data type code address          
     &       ppx_packing_acc    =31,  ! Packing accuracy code address (1
     &       ppx_pack_acc       =31)                                    
C
C Valid grid type codes
      INTEGER
     &       ppx_atm_nonstd,ppx_atm_tall,ppx_atm_tland,ppx_atm_tsea,
     &       ppx_atm_uall,ppx_atm_uland,ppx_atm_usea,ppx_atm_compressed,
     &       ppx_atm_ozone,ppx_atm_tzonal,ppx_atm_uzonal,ppx_atm_rim,
     &       ppx_atm_tmerid,ppx_atm_umerid,ppx_atm_scalar,
     &       ppx_atm_cuall,ppx_atm_cvall,
     &       ppx_ocn_nonstd,ppx_ocn_tall,ppx_ocn_tcomp,ppx_ocn_tfield,
     &       ppx_ocn_uall,ppx_ocn_ucomp,ppx_ocn_ufield,
     &       ppx_ocn_tzonal,ppx_ocn_uzonal,ppx_ocn_tmerid,
     &       ppx_ocn_umerid,ppx_ocn_scalar,ppx_ocn_rim,
     &       ppx_ocn_cuall,ppx_ocn_cvall,    
     &       ppx_wam_all,ppx_wam_sea,ppx_wam_rim
C Valid rotation type codes
      INTEGER
     &       ppx_unrotated,ppx_elf_rotated
C Valid level type codes
      INTEGER
     &       ppx_full_level,ppx_half_level
C Valid data type codes
      INTEGER
     &       ppx_type_real,ppx_type_int,ppx_type_log
C Valid meto8 level type codes
      INTEGER
     &       ppx_meto8_surf
C Valid dump packing codes
      INTEGER
     &       ppx_pack_off,ppx_pack_32,ppx_pack_wgdos,ppx_pack_cfi1
C
C
C
C
      PARAMETER(
     &       ppx_atm_nonstd=0,      ! Non-standard atmos grid
     &       ppx_atm_tall=1,        ! All T points (atmos)
     &       ppx_atm_tland=2,       ! Land-only T points (atmos)
     &       ppx_atm_tsea=3,        ! Sea-only T points (atmos)
     &       ppx_atm_tzonal=4,      ! Zonal field at T points (atmos)
     &       ppx_atm_tmerid=5,      ! Merid field at T points (atmos)
     &       ppx_atm_uall=11,       ! All u points (atmos)
     &       ppx_atm_uland=12,      ! Land-only u points (atmos)
     &       ppx_atm_usea=13,       ! Sea-only u points (atmos)
     &       ppx_atm_uzonal=14,     ! Zonal field at u points (atmos)
     &       ppx_atm_umerid=15,     ! Merid field at u points (atmos)
     &       ppx_atm_scalar=17,     ! Scalar (atmos)
     &       ppx_atm_cuall=18,      ! All C-grid (u) points (atmos)
     &       ppx_atm_cvall=19,      ! All C-grid (v) points (atmos)
     &       ppx_atm_compressed=21, ! Compressed land points (atmos)
     &       ppx_atm_ozone=22,      ! Field on ozone grid (atmos)
     &       ppx_atm_rim=25,        ! Rim type field (LAM BCs atmos)
     &       ppx_ocn_nonstd=30,     ! Non-standard ocean grid
     &       ppx_ocn_tcomp=31,      ! Compressed T points (ocean)
     &       ppx_ocn_ucomp=32,      ! Compressed u points (ocean)
     &       ppx_ocn_tall=36,       ! All T points incl. cyclic (ocean)
     &       ppx_ocn_uall=37,       ! All u points incl. cyclic (ocean)
     &       ppx_ocn_cuall=38,      ! All C-grid (u) points (ocean)
     &       ppx_ocn_cvall=39,      ! All C-grid (v) points (ocean)
     &       ppx_ocn_tfield=41,     ! All non-cyclic T points (ocean)
     &       ppx_ocn_ufield=42,     ! All non-cyclic u points (ocean)
     &       ppx_ocn_tzonal=43,     ! Zonal n-c field at T points(ocean)
     &       ppx_ocn_uzonal=44,     ! Zonal n-c field at u points(ocean)
     &       ppx_ocn_tmerid=45,     ! Merid n-c field at T points(ocean)
     &       ppx_ocn_umerid=46,     ! Merid n-c field at u points(ocean)
     &       ppx_ocn_scalar=47,     ! Scalar (ocean)
     &       ppx_ocn_rim=51,        ! Rim type field (LAM BCs ocean)    
     &       ppx_wam_all=60,        ! All points (wave model)
     &       ppx_wam_sea=62,        ! Sea points only (wave model)
     &       ppx_wam_rim=65)        ! Rim type field (LAM BCs wave)
C
      PARAMETER(
     &       ppx_unrotated=0,       ! Unrotated output field
     &       ppx_elf_rotated=1)     ! Rotated ELF field
C
      PARAMETER(
     &       ppx_full_level=1,      ! Model full level
     &       ppx_half_level=2)      ! Model half level
C
      PARAMETER(
     &       ppx_type_real=1,       ! Real data type
     &       ppx_type_int=2,        ! Integer data type
     &       ppx_type_log=3)        ! Logical data type
C
      PARAMETER(
     &       ppx_meto8_surf=9999)   ! MetO8 surface type code
C
      PARAMETER(
     &       ppx_pack_off=0,        ! Field not packed (ie. 64 bit)
     &       ppx_pack_32=-1,        ! Field packed to 32 bit in dump
     &       ppx_pack_wgdos=1,      ! Field packed by WGDOS method
     &       ppx_pack_cfi1=11)      ! Field packed using CFI1 (ocean)
C
CL
CL external routines
CL
      EXTERNAL stextc ! extracts the field
      EXTERNAL stcolm ! computes the column mean
      EXTERNAL stzonm ! computes the zonal mean
      EXTERNAL stmerm ! computes the meridional mean
      EXTERNAL stglom ! computes the global mean
      EXTERNAL stfieldm ! computes the field mean
CL
CL local variables
CL
      LOGICAL lwrap                ! TRUE if output field wraparound EW
      LOGICAL lmasswt_strict       ! copy of lmasswt - but set to false
!                                  ! if mass weighting is not requested
      INTEGER xstart,ystart        ! lower left hand corner coords
      INTEGER xend,yend            ! upper right hand corner coords
      INTEGER processing_code      ! what kind of mean  will be done
      INTEGER what_level           ! what type of input level
      INTEGER what_mask            ! what mask is used
      INTEGER what_weight          ! what weighting is used

      INTEGER i,j  ! loop counters
     &,       n_rows  ! number of rows to loop over

      INTEGER
! global versions of the extracted area domain limits
     &  global_xstart,global_xend,global_ystart,global_yend

! workspace arrays containining weighting factors and masks.
      REAL
     &  area_weight(row_length,p_rows)
     &, pstar_weight(row_length,p_rows)
      LOGICAL
     &  mask(row_length,p_rows)


CL----------------------------------------------------------------------
CL 1. Set up local variables
CL
      xstart=control(st_west_code)
      xend=control(st_east_code)
      ystart=control(st_north_code)  ! NOTE: Grid is assumed to be
      yend=control(st_south_code)    !       oriented north-to-south

      global_xstart=xstart
      global_ystart=ystart
      global_xend=xend
      global_yend=yend

! and calculate what the local subdomain limits are:
        CALL GLOBAL_TO_LOCAL_SUBDOMAIN( .TRUE.,.TRUE.,
     &                                  GR,mype,
     &                                  global_ystart,global_xend,
     &                                  global_yend,global_xstart,
     &                                  ystart,xend,yend,xstart)

C Check if wraparound field
      IF (xstart.GT.xend) THEN
        IF (lcyclic) THEN
          xend=xend+row_length-2*Offx
! subtract two halos as we don't wish to include halos at the end
! and start of the row within the wrap around domain
          lwrap=.TRUE.
        ELSE
          icode=st_bad_wraparound    ! wraparound illegal unless cyclic
          GOTO 999
        ENDIF
      ELSE
        lwrap=.FALSE.
      ENDIF
      IF (global_xstart .GT. global_xend) THEN
        IF (lcyclic) THEN
          global_xend=global_xend+glsize(1)
        ELSE
          icode=st_bad_wraparound  ! wraparound illegal unless cyclic
          GOTO 999
        ENDIF
      ENDIF
      processing_code=control(st_gridpoint_code)
      what_level=control(st_input_bottom)
      what_mask=mod(processing_code,block_size)
      what_weight=control(st_weight_code)
CL
CL 1.1 Prevent masking or weighting if input field is not 2D in extent
CL     - weighting and masking is assumed to have been done outside.
CL
      IF ( (.NOT.(st_grid.EQ.st_tp_grid.OR.st_grid.EQ.st_uv_grid.OR.
     *            st_grid.EQ.st_cu_grid.OR.st_grid.EQ.st_cv_grid))
     *      .AND.(what_mask  .NE.stash_null_mask_code   .OR.
     *            what_weight.NE.stash_weight_null_code) ) THEN
       icode=st_not_supported
       cmessage='SPATIAL : Masking/weighting unsupported - non 2D field'
       GOTO 999
      ENDIF

! Check for supported weighting and masking options

      IF (.NOT. ((what_weight .EQ. stash_weight_null_code) .OR.
     &           (what_weight .EQ. stash_weight_area_code) .OR.
     &           (what_weight .EQ. stash_weight_volume_code) .OR.
     &           (what_weight .EQ. stash_weight_mass_code) ) ) THEN
        cmessage='SPATIAL : Unrecognized weighting option'
        icode=unknown_weight
        GOTO 999
      ENDIF

      IF (.NOT. ((what_mask .EQ. stash_null_mask_code) .OR.
     &           (what_mask .EQ. stash_land_mask_code) .OR.
     &           (what_mask .EQ. stash_sea_mask_code ) ) ) THEN
        cmessage='SPATIAL : Unrecognized masking option'
        icode=unknown_mask
        GOTO 999
      ENDIF

      IF (what_weight .EQ. stash_weight_volume_code) THEN
        cmessage='SPATIAL : Volume-weighting not supported'
        icode=st_illegal_weight
        GOTO 999
      ENDIF

! Set lmasswt_strict - copy of lmasswt, but set to false is mass
! weighting not requested

      lmasswt_strict=
     &  (lmasswt .AND. (what_weight .EQ. stash_weight_mass_code))

! Precalculate weighting and masking arrays
! I've used IF tests inside the loops, but since the logical
! expressions are invariant wrt i and j, the compiler will
! move them outside the DO loops. It makes the code a lot shorter!

      IF (.NOT. ((st_grid .EQ. st_tp_grid) .OR.
     &         (st_grid .EQ. st_cu_grid))) THEN
        n_rows=u_rows
      ELSE
        n_rows=p_rows
      ENDIF

! area weighting
      DO j=1,n_rows
        DO i=1,row_length
          IF (what_weight .EQ. stash_weight_null_code) THEN
            area_weight(i,j)=1.0  ! no area weighting
          ELSE ! some form of area weighting will be required
            IF (st_grid .EQ. st_tp_grid .OR.
     &          st_grid .EQ. st_cu_grid) THEN
              area_weight(i,j)=cos_p_latitude(i,j)
            ELSE
              area_weight(i,j)=cos_u_latitude(i,j)
            ENDIF
          ENDIF
        ENDDO
      ENDDO

C Ensure halos are initialised
      DO I=n_rows-1,n_rows
        DO J=1,row_length
          pstar_weight(j,i)=1.0
        ENDDO
      ENDDO
! mass weighting
      IF ((what_weight .EQ. stash_weight_null_code) .OR.
     &    (what_weight .EQ. stash_weight_area_code)) THEN
! No mass weighting is required
        DO j=1,n_rows
          DO i=1,row_length
            pstar_weight(i,j)=1.0
          ENDDO
        ENDDO
      ELSE
        IF (st_grid .EQ. st_uv_grid) THEN
          CALL P_TO_UV(pstar,pstar_weight,row_length*p_rows,
     &                 row_length*u_rows,row_length,p_rows)
        ELSEIF (st_grid .EQ. st_cu_grid) THEN
          CALL P_TO_CU(pstar,pstar_weight,row_length*p_rows,
     &                 row_length*p_rows,row_length,p_rows)
        ELSEIF (st_grid .EQ. st_cv_grid) THEN
          CALL P_TO_CV(pstar,pstar_weight,row_length*p_rows,
     &                 row_length*u_rows,row_length,p_rows)
        ELSE
          DO j=1,n_rows
            DO i=1,row_length
              pstar_weight(i,j)=pstar(i,j)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

! masking

      DO j=1,p_rows
        DO i=1,row_length
          IF (what_mask .EQ. stash_land_mask_code) THEN
            mask(i,j)=land(i,j)
          ELSEIF (what_mask .EQ. stash_sea_mask_code) THEN
            mask(i,j)=.NOT. land(i,j)
          ELSE
            mask(i,j)=.TRUE.
          ENDIF
        ENDDO
      ENDDO

CL----------------------------------------------------------------------
CL 2. Call service routine to perform required processing
CL
CL 2.1 Extract sub field (single level at a time)
CL
      IF (processing_code.lt.extract_top.and.
     +    processing_code.gt.extract_base) THEN
        n_rows_out=(yend+1)-ystart
        n_cols_out=(xend+1)-xstart
        IF (
     &   (xstart .NE. st_no_data) .AND. (xend .NE. st_no_data) .AND.
     &   (ystart .NE. st_no_data) .AND. (yend .NE. st_no_data)) THEN
        CALL STEXTC(fieldin,vx,vy,st_grid,lwrap,lmasswt_strict,
     &              xstart,ystart,xend,yend,
     &              fieldout,
     &              pstar_weight,
     &              delta_ak(base_level),delta_bk(base_level),
     &              area_weight,mask,
     &              row_length,p_rows,
     &              what_level,what_mask,what_weight,rmdi,
     &              icode,cmessage)
        ELSE  ! just set values to non NaN
          DO i=1,lenout
            fieldout(i)=0.0
          ENDDO
        ENDIF
CL
CL 2.2 Calculate column mean (over multiple levels indexed by index_lev)
CL
      ELSEIF (processing_code.lt.vert_mean_top.and.
     +        processing_code.gt.vert_mean_base) THEN
        n_rows_out=yend+1-ystart
        n_cols_out=xend+1-xstart
        IF (
     &   (xstart .NE. st_no_data) .AND. (xend .NE. st_no_data) .AND.
     &   (ystart .NE. st_no_data) .AND. (yend .NE. st_no_data)) THEN
        CALL STCOLM(fieldin,vx,vy,vz,st_grid,lwrap,lmasswt_strict,
     &              xstart,ystart,xend,yend,
     &              fieldout,index_lev,level_list,no_of_levels,
     &              pstar_weight,
     &              delta_ak,delta_bk,
     &              area_weight,mask,
     &              row_length,p_rows,
     &              what_level,what_mask,what_weight,rmdi,
     &              icode,cmessage)
        ELSE  ! just set values to non NaN
          DO i=1,lenout
            fieldout(i)=0.0
          ENDDO
        ENDIF
CL
CL 2.3 Calculate zonal mean (single level at a time)
CL
      ELSEIF (processing_code.lt.zonal_mean_top.and.
     +        processing_code.gt.zonal_mean_base) THEN
        n_rows_out=yend+1-ystart
        n_cols_out=1
        CALL STZONM(fieldin,vx,vy,st_grid,gr,lwrap,lmasswt_strict,
     &              xstart,ystart,xend,yend,
     &              global_xstart,global_ystart,global_xend,global_yend,
     &              fieldout,
     &              pstar_weight,
     &              delta_ak(base_level),delta_bk(base_level),
     &              area_weight,mask,
     &              row_length,p_rows,
     &              what_level,what_mask,what_weight,rmdi,
     &              icode,cmessage)
CL
CL 2.4 Calculate meridional mean (single level at a time)
CL
      ELSEIF (processing_code.lt.merid_mean_top.and.
     +        processing_code.gt.merid_mean_base) THEN
        n_rows_out=1
        n_cols_out=xend+1-xstart
        CALL STMERM(fieldin,vx,vy,st_grid,gr,lwrap,lmasswt_strict,
     &              xstart,ystart,xend,yend,
     &              global_xstart,global_ystart,global_xend,global_yend,
     &              fieldout,
     &              pstar_weight,
     &              delta_ak(base_level),delta_bk(base_level),
     &              area_weight,mask,
     &              row_length,p_rows,
     &              what_level,what_mask,what_weight,rmdi,
     &              icode,cmessage)
CL
CL 2.5 Calculate field mean (single level at a time)
CL
      ELSEIF (processing_code.lt.field_mean_top.and.
     +        processing_code.gt.field_mean_base) THEN
        n_rows_out=1
        n_cols_out=1
        CALL STFIELDM(fieldin,vx,vy,st_grid,gr,lwrap,lmasswt_strict,
     &              xstart,ystart,xend,yend,
     &              global_xstart,global_ystart,global_xend,global_yend,
     &              fieldout,
     &              pstar_weight,
     &              delta_ak(base_level),delta_bk(base_level),
     &              area_weight,mask,
     &              row_length,p_rows,
     &              what_level,what_mask,what_weight,rmdi,
     &              icode,cmessage)
CL
CL 2.6 Calculate global mean (over multiple levels)
CL
      ELSEIF (processing_code.lt.global_mean_top.and.
     +        processing_code.gt.global_mean_base) THEN
        n_rows_out=1
        n_cols_out=1
        CALL STGLOM(fieldin,vx,vy,vz,st_grid,gr,lwrap,lmasswt_strict,
     &              xstart,ystart,xend,yend,
     &              global_xstart,global_ystart,global_xend,global_yend,
     &              fieldout,index_lev,level_list,no_of_levels,
     &              pstar_weight,
     &              delta_ak,delta_bk,
     &              area_weight,mask,
     &              row_length,p_rows,
     &              what_level,what_mask,what_weight,rmdi,
     &              icode,cmessage)
CL
CL 2.7 Invalid processing option
CL
      ELSE
        icode=unknown_processing
        write(cmessage,111)'unknown processing option',
     +    processing_code
      ENDIF
CL
  999 CONTINUE
111   format('SPATIAL : >>FATAL ERROR <<',a40,i5)
C
      RETURN
      END
