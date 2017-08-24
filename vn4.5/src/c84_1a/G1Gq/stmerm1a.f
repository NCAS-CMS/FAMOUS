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
CLL  Routine: STMERM ---------------------------------------------------
CLL
CLL  Purpose: Calculate weighted meridional mean within a region
CLL           specified by lower left hand and upper right hand corner.
CLL           (STASH service routine).
CLL
CLL  Tested under compiler:   cft77
CLL  Tested under OS version: UNICOS 5.1
CLL
CLL  Author:   T.Johns/S.Tett
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL   3.3  16/09/93  Allow level-by-level mass-weighting if mass-weights
CLL                  are so defined, otherwise use P*.
!LL   4.3  15/01/97  Moved weighting and masking calculations up to
!LL                  SPATIAL.
!LL                  Significantly rewritten for MPP mode - meridional
!LL                  data must be gathered to a processor for
!LL                  reproducible sums to be calculated.      P.Burton
!LL   4.4  13/06/97  MPP: Set fieldout to zero for processors in
!LL                  the result subdomain area which will not
!LL                  otherwise receiveof the meridional mean. P.Burton
!LL   4.5  12/01/98  Replaced usage of shmem common block by a
!LL                  dynamic array.                   P.Burton
CLL
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered: D714
CLL
CLL  Project task: D7
CLL
CLL  External documentation:
CLL    Unified Model Doc Paper C4 - Storage handling and diagnostic
CLL                                 system (STASH)
CLL
C*L  Interface and arguments: ------------------------------------------
C
      SUBROUTINE STMERM(fieldin,vx,vy,st_grid,gr,lwrap,lmasswt,
     &                  xstart,ystart,xend,yend,
     &                  global_xstart,global_ystart,
     &                  global_xend,global_yend,
     &                  fieldout,
     &                  pstar_weight,delta_ak,delta_bk,
     &                  area_weight,mask,
     &                  row_length,p_rows,
     &                  level_code,mask_code,weight_code,rmdi,
     &                  icode,cmessage)
C
      IMPLICIT NONE
C
      INTEGER
     &    vx,vy,                                ! IN  input field size
     &    st_grid,                              ! IN  STASH grdtype code
     &    gr,                                   ! IN input fld grid
     &    xstart,ystart,                        ! IN  lower LH corner
     &    xend,yend,                            ! IN  upper RH corner
     &    global_xstart,global_ystart,           ! IN global versions of
     &    global_xend,  global_yend,             ! IN xstart etc.
     &    row_length,p_rows,                    ! IN  primary dimensions
     &    level_code,                           ! IN  input level code
     &    mask_code,                            ! IN  masking code
     &    weight_code,                          ! IN  weighting code
     &    icode                                 ! OUT error return code
      CHARACTER*(*)
     &    cmessage                              ! OUT error return msg
      LOGICAL
     &    lwrap,                                ! IN  TRUE if wraparound
     &    lmasswt,                              ! IN  TRUE if masswts OK
     &    mask(row_length,p_rows)               ! IN  mask array
      REAL
     &    fieldin(vx,vy),                       ! IN  input field
     &    fieldout(xstart:xend),                ! OUT output field
     &    pstar_weight(row_length,p_rows),      ! IN  pstar mass weight
     &    delta_ak,                             ! IN  hybrid coordinates
     &    delta_bk,                             ! IN  hybrid coordinates
     &    area_weight(row_length,p_rows),       ! IN  area weighting
! (already interpolated to the correct grid and
!  set to 1.0 where no area weighting is required)
     &    rmdi                                  ! IN  missing data indic
C*----------------------------------------------------------------------
C
C External subroutines called
C
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

C
C Local variables
C
        INTEGER i,ii,j   ! ARRAY INDICES FOR VARIABLE


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

      INTEGER

! Processor co-ordinates of processors at the corners of the
! processed subdomain
     &  proc_top_left_x,proc_top_left_y
     &, proc_bot_right_x,proc_bot_right_y

! size of the full subarea in both horizonal dimensions
     &, merid_sum_global_len_x
     &, merid_sum_global_len_y

! loop variables for loops over processors in subdomain
     &, proc_x,proc_y

! real processor x co-ordinate - when proc_x > nproc_x is just
! proc_x-nproc_x
     &, eff_proc_x

! processor id of processor (proc_x,proc_y)
     &, proc_id

! definition of the extracted subarea array on processor proc_id
     &, local_array_top_left_x,local_array_top_left_y
     &, local_array_bot_right_x,local_array_bot_right_y

! definition of the real data contained within the extracted
! subarea on processor proc_id (ie. not including halos)
     &, local_data_top_left_x,local_data_top_left_y
     &, local_data_bot_right_x,local_data_bot_right_y

! size in the x dimension of the subarea array on proc_id
     &, local_array_size_x

! length of (partial) local meridional data to be sent
     &, local_send_size_y

! offset of data to be sent from start of local row
     &, local_send_off_x

! position in the full meridional column of (partial)
! data to be sent
     &, pos_in_merid_array

! number of (partial) meridional mean columns on
! processor proc_id
     &, local_n_cols_to_send

! the first point on the column to be sent to processor proc_id
     &, local_send_off_y

! the global meridional mean number of the first local
! (partial) meridional mean column to be sent from proc_id
     &, global_merid_col_number_start

! loop counter for loop over columns to send
     &, col

! index of column in proc_id's array of local data
     &, local_col

! global meridional column number of a column
     &, global_merid_col_id

! processor which this meridional mean column will be sent to
     &, dest_proc

! column number on the destination processor
     &, work_dest_col_id

! number of items of meridional data to send and receive
     &, n_send_data , n_rec_data

! number of final meridional means to send and receive
     &, n_send_means , n_rec_means

! number of columns of (full) meridional data on this processor
     &, n_cols_full_merid_data

! size of local_sum_arrays and global_sum_arrays
     &, local_sum_array_len
     &, global_sum_array_len

! field type (P or U) of input field
     &, fld_type

! arguments for GCOM routines
     &, flag , info

! dummy variables (unused return values from subroutine calls)
     &, dummy1,dummy2


      LOGICAL

! indicates if the subarea requested for meridional meaning wraps
! over zero longitude
     &  lwrap_merid_mean

! indicates if the subdomain contains processors which hold both
! the start and end of the subdomain, which wraps over zero
! longitude
     &, lwrap_proc

! indicates that a full field is being zonal meaned
     &, fullfield

      REAL
! temporary variables used in calculation of meridional means
     &  merid_sum_top,merid_sum_bot

      INTEGER

! Send/receive maps for meridional data arrays to be summed
     &  send_data_map(7,(xend-xstart+1))
     &, rec_data_map(7,(global_xend-global_xstart+1)*nproc)

! send/receive maps for meridional means
     &, send_means_map(7,(global_xend-global_xstart+1))
     &, rec_means_map(7,(xend-xstart+1))

! Weighted version of fieldin
      REAL local_sum_array_top(xstart:xend,ystart:yend)
! Weights applied to fieldin
      REAL local_sum_array_bot(xstart:xend,ystart:yend)


      INTEGER
! Sizes of the global_sum_arrays defined below
     &  global_sum_array_sizex,global_sum_array_sizey

      REAL
! Collected versions of fieldin and the weights containing
! whole (subarea) columns of meridional data
     &  global_sum_array_top(global_ystart:global_yend,
     &                       global_xend-global_xstart+1)
     &, global_sum_array_bot(global_ystart:global_yend,
     &                       global_xend-global_xstart+1)

! Calculated meridional means on the calculating processor
     &, merid_mean_array(global_xend-global_xstart+1)



! Integer function used for obtaining field type
      INTEGER GET_FLD_TYPE


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

        REAL SUMMTOP(xstart:xend)
        REAL SUMMBOT(xstart:xend)

CL----------------------------------------------------------------------
CL 0. Initialise sums
CL
CFPP$ NOINNER R
      DO i=xstart,xend
        SUMMTOP(i)=0.0
        SUMMBOT(i)=0.0
      ENDDO
CL----------------------------------------------------------------------

! pstar_weight and area_weight arrays contain appropriate
! weighting factors, interpolated to the correct grid, for
! mass weighting and area weighting respectively. If either type
! of weighting is not required, the relevant array is set to 1.0
! The mask array contains appropriate masking

! Create arrays of weighted data suitable to be summed

! Only do the calculations if some of the subarea is contained
! within this processor
      IF ((xstart .NE. st_no_data) .AND. (xend .NE. st_no_data) .AND.
     &    (ystart .NE. st_no_data) .AND. (yend .NE. st_no_data)) THEN

        DO i=xstart,xend
          IF ( lwrap .AND. (i .GT. (lasize(1)-Offx))) THEN
            ii=i-lasize(1)+2*Offx ! miss halos on wrap around
          ELSE
            ii=i
          ENDIF
          DO j=ystart,yend
            IF (mask(ii,j)) THEN
              IF (.NOT. lmasswt) THEN
                local_sum_array_bot(i,j)=
     &            pstar_weight(ii,j)*area_weight(ii,j)
                local_sum_array_top(i,j)=
     &            fieldin(ii,j)*pstar_weight(ii,j)*area_weight(ii,j)
              ELSE
                local_sum_array_bot(i,j)=
     &            -1.0*(delta_ak+delta_bk*pstar_weight(ii,j))*
     &            area_weight(ii,j)
                local_sum_array_top(i,j)=
     &            -1.0*fieldin(ii,j)*
     &            (delta_ak+delta_bk*pstar_weight(ii,j))*
     &            area_weight(ii,j)
              ENDIF
            ELSE
              local_sum_array_bot(i,j)=0.0
              local_sum_array_top(i,j)=0.0
            ENDIF
          ENDDO
        ENDDO

      ENDIF  ! if this processor contains any of the subarea


! Initialise fieldout array - so all PE's have valid data
! (Only PEs on top of subdomain get the meridional means)
      DO i=xstart,xend
        fieldout(i)=0.0
      ENDDO

! The local_sum_arrays must be distributed so that complete
! sub-area columns exist on processors, so that a reproducible sum
! can be carried out.
! The following code calculates where the local_sum_array data
! must be sent to, and where the final answers must be sent back to

! 0.0 : Initialise variables defining the size of the arrays
!       global_sum_arrays

      global_sum_array_sizex=global_xend-global_xstart+1
      global_sum_array_sizey=global_yend-global_ystart+1

      local_sum_array_len=((xend-xstart)+1)*((yend-ystart)+1)
      global_sum_array_len=global_sum_array_sizex*
     &                     global_sum_array_sizey

! Set a logicial indicating if the area being meaned is the
! full field

      fld_type=GET_FLD_TYPE(gr)

      fullfield= ((( global_xstart .EQ. 1) .AND.
     &             ( global_xend) .EQ. glsize(1)) .AND.
     &             ( global_ystart .EQ. 1) .AND.
     &             (((fld_type .EQ. fld_type_p) .AND.
     &               (global_yend .EQ. glsize(2))) .OR.
     &              ((fld_type .EQ. fld_type_u) .AND.
     &               (global_yend .EQ. glsize(2)-1))))

! Calculate the length of the full meridional subarea

      merid_sum_global_len_x=global_xend-global_xstart+1
      merid_sum_global_len_y=global_yend-global_ystart+1

! 1.0 Find the set of processors covering the requested sub-area

      CALL GLOBAL_TO_LOCAL_RC(gr,
     &   global_xstart , global_ystart,
     &   proc_top_left_x, proc_top_left_y,
     &   dummy1,dummy2)

      CALL GLOBAL_TO_LOCAL_RC(gr,
     &   global_xend,global_yend,
     &   proc_bot_right_x, proc_bot_right_y,
     &   dummy1,dummy2)

! Set a logical to indicate if the meridional mean area required
! wraps over zero longitude

      lwrap_merid_mean=
     & ((global_xend .GT. glsize(1)) .OR.
     &   (global_xend .LT. global_xstart))

! If there is a wrap around over 0 longitude, ensure that
! proc_bot_right_x > proc_top_left_x

      IF (lwrap_merid_mean)
     &  proc_bot_right_x=proc_bot_right_x+nproc_x

! Set up a logical to indicate if a processor in the subdomain
! contains both the start and end of a meridional mean which wraps
! over zero longitude. If TRUE, some extra work is required at
! this processor as it contains data for two non-contiguous parts
! of the meridional mean

      lwrap_proc=(proc_bot_right_x .EQ. proc_top_left_x+nproc_x)

! 2.0 Loop over all the processors in the subdomain, and set
!     up the send/receive maps defining the redistribution
!     of data

      n_send_data=0            ! number of items of data to send
      n_rec_data=0             ! number of items of data to receive
      n_send_means=0           ! number of merid. means I will send
      n_rec_means=0            ! number of merid. means I will receive
      n_cols_full_merid_data=0 ! number of cols. of data I will mean

      DO proc_y=proc_top_left_y , proc_bot_right_y

        DO proc_x=proc_top_left_x , proc_bot_right_x

          eff_proc_x=MOD(proc_x,nproc_x)
          proc_id=eff_proc_x+proc_y*nproc_x

! 2.1  Find the size of the array containing the meridional
!      arrays on processor proc_id

          CALL GLOBAL_TO_LOCAL_SUBDOMAIN(.TRUE.,.TRUE.,
     &      gr,proc_id,
     &      global_ystart,global_xend,
     &      global_yend,global_xstart,
     &      local_array_top_left_y,local_array_bot_right_x,
     &      local_array_bot_right_y,local_array_top_left_x)

! 2.2 Using this information, calculate the size of this array in
!     the x dimension. If the data is wrapped round, the calculation
!     is done differently:

          IF (local_array_top_left_x .LE. local_array_bot_right_x)
     &    THEN
            local_array_size_x=
     &        local_array_bot_right_x-local_array_top_left_x+1
          ELSE
            local_array_size_x=
     &        local_array_bot_right_x-local_array_top_left_x+1+
     &        g_lasize(1,proc_id)-2*Offx
          ENDIF

! 2.3 Find out the size of the actual meridional mean data within the
!     subarea array on processor proc_id

          CALL GLOBAL_TO_LOCAL_SUBDOMAIN(.FALSE.,.FALSE.,
     &      gr,proc_id,
     &      global_ystart,global_xend,
     &      global_yend,global_xstart,
     &      local_data_top_left_y,local_data_bot_right_x,
     &      local_data_bot_right_y,local_data_top_left_x)

! 2.4 Calculate various quantities:
!     local_send_size_y  : the length of data to be sent
!     local_send_off_y   : the offset of this data from the
!                          start of column
!     pos_in_merid_array : position of this data in the full
!                          meridional array

          local_send_size_y=
     &      local_data_bot_right_y-local_data_top_left_y+1
          local_send_off_y=
     &      local_data_top_left_y-local_array_top_left_y
          pos_in_merid_array=
     &      g_datastart(2,proc_id)+local_data_top_left_y-Offy-
     &      global_ystart

! 2.5 Find the number of meridional mean columns to be sent
!     from this processor, the first column to be sent,
!     and which global meridional mean this is

          IF ((LWRAP_PROC) .AND. (proc_x .EQ. proc_top_left_x)) THEN
! Processor containing start and end of sumdomain - but here
! we're interested only in the start segment

            local_n_cols_to_send=
     &        g_lasize(1,proc_id)-local_data_top_left_x-Offx+1
            local_send_off_x=
     &        local_data_top_left_x-local_array_top_left_x
            global_merid_col_number_start=
     &        g_datastart(1,proc_id)+local_data_top_left_x-Offx-
     &        global_xstart

          ELSEIF ((LWRAP_PROC) .AND.
     &            (proc_x .EQ. proc_bot_right_x)) THEN
! Processor containing start and end of subdomain - but here
! we're interested only in the end segment

            local_n_cols_to_send=local_data_bot_right_x-Offx
            local_send_off_x=local_array_size_x-local_n_cols_to_send
            global_merid_col_number_start=
     &        merid_sum_global_len_x-local_n_cols_to_send+1

          ELSE
! all other processors

            local_n_cols_to_send=
     &        local_data_bot_right_x-local_data_top_left_x+1
            local_send_off_x=
     &        local_data_top_left_x-local_array_top_left_x
            global_merid_col_number_start=
     &        g_datastart(1,proc_id)+local_data_top_left_x-Offx-
     &        global_xstart
          ENDIF

          IF (global_merid_col_number_start .LT. 1) THEN
! This means the sub-area wraps over zero longitude - so to get
! the correct position in the array we add the global row length
            global_merid_col_number_start=
     &        global_merid_col_number_start+glsize(1)
          ENDIF


! 2.6 Loop over columns and construct send/receive maps

          DO col=1,local_n_cols_to_send

! 2.6.1 Find the local column index on proc_id, and the global
!       meridional column index of this column

            local_col=col+local_send_off_x
            global_merid_col_id=global_merid_col_number_start+col-1

! 2.6.2 and find the destination processor of this column, and
!       where on this processor it will be sent to

            dest_proc=MOD(global_merid_col_id-1,nproc)
            work_dest_col_id=((global_merid_col_id-1)/nproc)+1

! 2.6.3 If this processor is proc_id construct a send_data_map
!       entry for this column of data

            IF (mype .EQ. proc_id) THEN

              n_send_data = n_send_data+1

              send_data_map(S_DESTINATION_PE,n_send_data)=
     &          dest_proc
              send_data_map(S_BASE_ADDRESS_IN_SEND_ARRAY,
     &                      n_send_data)=
     &          local_send_off_y*local_array_size_x +
     &          local_send_off_x+col
              send_data_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,
     &                      n_send_data)=
     &          local_send_size_y
              send_data_map(S_STRIDE_IN_SEND_ARRAY,n_send_data)=
     &          local_array_size_x
              send_data_map(S_ELEMENT_LENGTH,n_send_data)=1
              send_data_map(S_BASE_ADDRESS_IN_RECV_ARRAY,
     &                      n_send_data)=
     &          (work_dest_col_id-1)*global_sum_array_sizey +
     &          pos_in_merid_array
              send_data_map(S_STRIDE_IN_RECV_ARRAY,n_send_data)=1

! 2.6.3.1 If this processor is at the top of the subarea, then it is
!         responsible for holding the final meridional mean values.
!         So we must set up a rec_means_map entry to allow the
!         meridional mean value for this column to be returned.

              IF (proc_y .EQ. proc_top_left_y) THEN

                n_rec_means = n_rec_means + 1

                rec_means_map(R_SOURCE_PE,n_rec_means)=dest_proc
                IF (fullfield) THEN ! We don't want halos
                  rec_means_map(R_BASE_ADDRESS_IN_RECV_ARRAY,
     &                          n_rec_means)=
     &              local_col-Offx
                ELSE ! halos are automatically removed
                  rec_means_map(R_BASE_ADDRESS_IN_RECV_ARRAY,
     &                          n_rec_means)=
     &              local_col
                ENDIF
                rec_means_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,
     &                        n_rec_means)= 1
                rec_means_map(R_STRIDE_IN_RECV_ARRAY,
     &                        n_rec_means)= 1
                rec_means_map(R_ELEMENT_LENGTH,n_rec_means)=1
                rec_means_map(R_BASE_ADDRESS_IN_SEND_ARRAY,
     &                        n_rec_means)=
     &            work_dest_col_id
                rec_means_map(R_STRIDE_IN_SEND_ARRAY,
     &                        n_rec_means)=1
              ENDIF

            ENDIF

! 2.6.4 If this processor is dest_proc construct a rec_data_map
!       entry for this column of data

            IF (mype .EQ. dest_proc) THEN

              IF (proc_y .EQ. proc_top_left_y)
! increment counter of full meridional columns on this processor
     &          n_cols_full_merid_data=n_cols_full_merid_data+1

              n_rec_data = n_rec_data+1

              rec_data_map(R_SOURCE_PE,n_rec_data)=
     &          proc_id
              rec_data_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_rec_data)=
     &          (work_dest_col_id-1)*global_sum_array_sizey +
     &          pos_in_merid_array
              rec_data_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_rec_data)=
     &          local_send_size_y
              rec_data_map(R_STRIDE_IN_RECV_ARRAY,n_rec_data)=1
              rec_data_map(R_ELEMENT_LENGTH,n_rec_data)=1
              rec_data_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_rec_data)=
     &          local_send_off_y*local_array_size_x +
     &          local_send_off_x+col
              rec_data_map(R_STRIDE_IN_SEND_ARRAY,n_rec_data)=
     &          local_array_size_x

! 2.6.4.1 Set up the send_means_map entry for sending the
!         resulting meridional mean of this column back to
!         the processor at the top of the subarea.
!         We only need to do this once per column (not for
!         each value of proc_y).

              IF (proc_y .EQ. proc_top_left_y) THEN

                n_send_means = n_send_means+1

                send_means_map(S_DESTINATION_PE,n_send_means)=
     &            proc_id
                send_means_map(S_BASE_ADDRESS_IN_SEND_ARRAY,
     &                         n_send_means)=work_dest_col_id
                send_means_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,
     &                         n_send_means)=1
                send_means_map(S_STRIDE_IN_SEND_ARRAY,
     &                         n_send_means)=1
                send_means_map(S_ELEMENT_LENGTH,n_send_means)=1
                IF (fullfield) THEN ! we don't want halos
                  send_means_map(S_BASE_ADDRESS_IN_RECV_ARRAY,
     &                           n_send_means)=local_col-Offx
                ELSE ! halos are automatically removed
                  send_means_map(S_BASE_ADDRESS_IN_RECV_ARRAY,
     &                           n_send_means)=local_col
                ENDIF
                send_means_map(S_STRIDE_IN_RECV_ARRAY,
     &                         n_send_means)=1
              ENDIF ! if at top of subarea

            ENDIF ! if mype .eq. dest_proc

          ENDDO ! col : loop over local columns on proc_id

        ENDDO ! proc_x : loop over processors in x dimension

      ENDDO ! proc_y : loop over processors in y dimension

! 3.0 Now the send and receive maps are set up, use
!     GCG_RALLTOALLE to redistribute the data
!     from the local_sum_arrays to the global_sum_arrays
      flag=GC_NONE ! flag argument is currently ignored by GCOM

! We do a SHM_PUT operation, because the send arrays are non-
! memory aligned, and the receiving arrays are.
      CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)
      info=GC_NONE

      CALL GCG_RALLTOALLE(
     &  local_sum_array_top , send_data_map , n_send_data ,
     &  local_sum_array_len ,
     &  global_sum_array_top , rec_data_map , n_rec_data ,
     &  global_sum_array_len ,
     &  gc_all_proc_group , flag , info)

      CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)
      info=GC_NONE

      CALL GCG_RALLTOALLE(
     &  local_sum_array_bot , send_data_map , n_send_data ,
     &  local_sum_array_len ,
     &  global_sum_array_bot , rec_data_map , n_rec_data ,
     &  global_sum_array_len ,
     &  gc_all_proc_group , flag , info)

! 4.0 Calculate mean of any meridional data on this processor

      DO i=1,n_cols_full_merid_data

        merid_sum_top=0.0
        merid_sum_bot=0.0

        DO j=global_ystart,global_yend

          merid_sum_top=merid_sum_top+
     &                     global_sum_array_top(j,i)
          merid_sum_bot=merid_sum_bot+
     &                     global_sum_array_bot(j,i)
        ENDDO

        IF (merid_sum_bot .EQ. 0.0) THEN
          merid_mean_array(i)=rmdi
        ELSE
          merid_mean_array(i)=merid_sum_top/merid_sum_bot
        ENDIF

      ENDDO

! 5.0 Send the calculated means back to the processors at the
!     top of the subarea, into the fieldout array

      flag=GC_NONE ! flag argument is currently ignored by GCOM

! We do a SHM_GET operation, because the send arrays are memory
! aligned, but the receiving arrays are not.
      CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_PUT,info)
      info=GC_NONE

      CALL GCG_RALLTOALLE(
     &  merid_mean_array , send_means_map , n_send_means ,
     &  global_sum_array_sizex,
     &  fieldout , rec_means_map , n_rec_means,
     &  (xend-xstart)+1,
     &  gc_all_proc_group , flag , info)

CL
  999 CONTINUE
      RETURN
      END
