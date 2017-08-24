C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!
!+ Scatters STASHed data from one processor to many processors
!
! Subroutine interface:
      SUBROUTINE STASH_SCATTER_FIELD (
     &  LOCAL_FIELD , GLOBAL_FIELD ,
     &  LOCAL_SIZE, GLOBAL_SIZE, LEVELS,
     &  GLOBAL_NORTH , GLOBAL_EAST_IN , GLOBAL_SOUTH , GLOBAL_WEST,
     &  GRIDTYPE_CODE ,
     &  SCATTER_PE,
     &  ICODE, CMESSAGE)

      IMPLICIT NONE


! Description:
! Takes a decomposed, STASH processed field and gathers
! it to a single processor, ready for I/O,
!
! Method:
! See in-line documentation
!
! Current code owner : P.Burton
!
! History:
!  Model    Date      Modification history from model version 4.3
!  version
!   4.3     17/09/96  New DECK created for MPP version of STASH
!                     P.Burton
!   4.3     13/06/97  More robust handling of zonal fields   P.Burton
!                     Fix for global_start_row               P.Burton
!
! Subroutine arguments:

      INTEGER
     &  LOCAL_SIZE      ! IN: size of level of LOCAL_FIELD
     &, GLOBAL_SIZE     ! IN: size of level of GLOBAL_FIELD
     &, LEVELS          ! IN: number of levels
     &, GLOBAL_NORTH    ! IN: specification of subdomain boundaries
     &, GLOBAL_EAST_IN  ! IN: ""
     &, GLOBAL_SOUTH    ! IN: ""
     &, GLOBAL_WEST     ! IN: ""
     &, GRIDTYPE_CODE   ! IN: indicates the type of grid output
     &, SCATTER_PE      ! IN: the PE to scatter global field from
     &, ICODE           ! OUT: return code, 0=OK

      REAL
     &  LOCAL_FIELD(LOCAL_SIZE,LEVELS)
!        ! OUT : local scattered data
     &, GLOBAL_FIELD(GLOBAL_SIZE,LEVELS)
!        ! IN : (PE SCATTER_PE only) - full field

      CHARACTER*80
     &  CMESSAGE        ! OUT: Error message if ICODE .NE. 0

! Parameters and common blocks
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
!====================== COMDECK AMAXSIZE ========================
! Description
!   This comdeck provides parameters giving the maximum likely sizes
!   of key UM resolution variables, useful for sizing static arrays.
!
!   History:
!   Model    Date     Modification history
!  version
!   4.2      18/11/96 New comdeck created.  P.Burton
!   4.3      24/01/97 Define MaxFieldSize to be a quarter of the
!                     SHMEM common block size.         P.Burton
!   4.4      3/7/97   Add MaxFieldSizeMes. Deborah Salmond
!   4.5     12/01/98  Added new variables, and changed sizes to
!                     correspond to global hi-res forecast - current
!                     largest configuration.                P.Burton
!                     Changed MAX_SHMEM_COMMON_SIZE to 3000000
!                     required for operational data assimilation.
!                                                           P.Burton

      INTEGER

     &  ROW_LENGTH_MAX  ! Maximum row length
     &, P_ROWS_MAX      ! Maximum number of rows
     &, HORIZ_DIM_MAX   ! MAX(ROW_LENGTH_MAX,P_ROWS_MAX)
     &, HALO_MAX        ! Maximum MPP halo width
     &, P_LEVELS_MAX    ! Maximum number of total levels
     &, Q_LEVELS_MAX    ! Maximum number of wet levels

      PARAMETER ( ROW_LENGTH_MAX = 432
     &,           P_ROWS_MAX = 325
     &,           HORIZ_DIM_MAX = 432
     &,           HALO_MAX = 2  ! fourth order double width halo
     &,           P_LEVELS_MAX = 42
     &,           Q_LEVELS_MAX = 42)

! Derived sizes

      INTEGER
     &  Max2DFieldSize
     &, Max3DFieldSize
     &, MaxHaloSize

      PARAMETER (
     &  Max2DFieldSize = ROW_LENGTH_MAX*P_ROWS_MAX
     &, Max3DFieldSize = ROW_LENGTH_MAX*P_ROWS_MAX*P_LEVELS_MAX
     &, MaxHaloSize = HORIZ_DIM_MAX*HALO_MAX
     & )

      INTEGER
     &  MAX_SHMEM_COMMON_SIZE,
     &  MaxFieldSize,
     &  MaxFieldSizeMes                                                 
      PARAMETER ( MAX_SHMEM_COMMON_SIZE = 3000000 ,
     &            MaxFieldSize   = MAX_SHMEM_COMMON_SIZE/4 ,
     &            MaxFieldSizeMes= MAX_SHMEM_COMMON_SIZE/6 )

! Local variables

      INTEGER
     &  GLOBAL_EAST   ! copy of GLOBAL_EAST_IN with wrap around s.t.
!                     ! GLOBAL_EAST > GLOBAL_ROW_LEN
     &, global_x      ! size of global data EW
     &, global_y      ! size of global data NS
     &, fld_type      ! indicates if field is on P or U grid
     &, level         ! loop index for loop over levels
     &, proc_topleft_x,proc_topleft_y   ! processors at corners of
     &, proc_botright_x,proc_botright_y ! the subarea
     &, dummy1,dummy2 ! ignored return arguments
     &, procx,procy   ! loop indexes for loops over processors
     &, eff_procx     ! real x co-ord of processor column procx
     &, procid        ! processor id of (procx,procy)
     &, local_xstart,local_xend  ! boundaries of subdomain for
     &, local_ystart,local_yend  ! processor procid
     &, local_start_row   ! first row to receive on procid
     &, local_start_col   ! first column to receive on procid
     &, sendsize_x        ! number of points on each row to send
!                         ! to procid
     &, nrows_to_send     ! number of rows to send to procid
     &, local_row_length  ! size of receiving array EW
     &, global_start_row  ! first row to send on SCATTER_PE
     &, global_start_col  ! first col. to send on SCATTER_PE
     &, global_row_length ! size of sending array EW
     &, flag,info         ! GCOM arguments

! Copies of arguments / variables used to decide if we can use the
! send/receive maps used in the last call

      INTEGER
     &  old_LOCAL_SIZE , old_GLOBAL_SIZE
     &, old_GLOBAL_NORTH , old_GLOBAL_EAST_IN
     &, old_GLOBAL_SOUTH , old_GLOBAL_WEST
     &, old_GRIDTYPE_CODE , old_SCATTER_PE
     &, old_current_decomp_type

      INTEGER
! variables defining send and receive maps to be passed to
! GC_RALL_TO_ALL, defining the data transposition
     &  send_map(7,MAXPROC*2)
     &, receive_map(7,2)
     &, n_sends,n_recvs  ! number of sends and receives


      LOGICAL
     &  wrap  ! if the subdomain wraps over 0 degree meridion
     &, wrap_special ! if there is a wrap around, which starts and
!                      ends on the same processor
     &, zonal_data   ! if this is a zonal data grid
     &, fullfield    ! if this is a full field - NOT a subarea

! Save all the variables that may be used in the next call
      SAVE
     &  old_LOCAL_SIZE , old_GLOBAL_SIZE
     &, old_GLOBAL_NORTH , old_GLOBAL_EAST_IN
     &, old_GLOBAL_SOUTH , old_GLOBAL_WEST
     &, old_GRIDTYPE_CODE , old_SCATTER_PE
     &, old_current_decomp_type
     &, send_map,receive_map,n_sends,n_recvs

! Set all the old_* variables to a number indicating they've
! not been used yet

      DATA
     &  old_LOCAL_SIZE , old_GLOBAL_SIZE
     &, old_GLOBAL_NORTH , old_GLOBAL_EAST_IN
     &, old_GLOBAL_SOUTH , old_GLOBAL_WEST
     &, old_GRIDTYPE_CODE , old_SCATTER_PE
     &, old_current_decomp_type
     &  / -1,-1,-1,-1,-1,-1,-1,-1,-1 /

! Functions

      INTEGER GET_FLD_TYPE
! ------------------------------------------------------------------

      ICODE=0

! See if there is wrap around over meridion, and if so make
! sure that GLOBAL_EAST is > glsize(1)

      GLOBAL_EAST=GLOBAL_EAST_IN
      IF (GLOBAL_EAST .GT. glsize(1)) THEN
        wrap=.TRUE.
      ELSEIF (GLOBAL_EAST .LT. GLOBAL_WEST) THEN
        wrap=.TRUE.
        GLOBAL_EAST=GLOBAL_EAST_IN+glsize(1)
      ELSE
        wrap=.FALSE.
      ENDIF

      IF ((GRIDTYPE_CODE .EQ. ppx_atm_tzonal) .OR. ! Atmos T zonal
     &   ( GRIDTYPE_CODE .EQ. ppx_atm_uzonal) .OR. ! Atmos U zonal
     &   ( GRIDTYPE_CODE .EQ. ppx_ocn_tzonal) .OR. ! Ocean T zonal
     &   ( GRIDTYPE_CODE .EQ. ppx_ocn_uzonal))     ! Atmos U zonal
     &  THEN

! This is a zonal field

        zonal_data=.TRUE.
        global_x=1

        IF ((GRIDTYPE_CODE .EQ. ppx_atm_tzonal) .OR. ! Atmos T zonal
     &      ( GRIDTYPE_CODE .EQ. ppx_ocn_tzonal))    ! Ocean T zonal
     &  THEN
          fld_type=fld_type_p
        ELSE
          fld_type=fld_type_u
        ENDIF
      ELSE

! This is a normal field

        zonal_data=.FALSE.
        global_x=glsize(1)

        fld_type=GET_FLD_TYPE(GRIDTYPE_CODE)

        IF (fld_type .EQ. fld_type_unknown) THEN
          WRITE(6,*) 'STASH_GATHER_FIELD encountered ',
     &      'field with gridtype code ',GRIDTYPE_CODE
          WRITE(6,*) 'Unable to process this field.'
          CMESSAGE='MPP : STASH_GATHER_FIELD could not process field'
          ICODE=1
          GOTO 9999
        ENDIF

      ENDIF

      IF (fld_type .EQ. fld_type_p) THEN
        global_y=glsize(2)
      ELSE
        global_y=glsize(2)-1
      ENDIF

! Set up logical indicating if this is a full field, or just
! a subdomain

      IF (zonal_data) THEN

        fullfield= ( ( GLOBAL_NORTH .EQ. 1) .AND.
     &             ( GLOBAL_SOUTH .EQ. global_y))

      ELSE

        fullfield = (( GLOBAL_WEST .EQ. 1) .AND.
     &               ( GLOBAL_EAST .EQ. global_x) .AND.
     &               ( GLOBAL_NORTH .EQ. 1) .AND.
     &               ( GLOBAL_SOUTH .EQ. global_y))

      ENDIF

! If this is a fullfield, we can simply use the standard
! SCATTER_FIELD routine

      IF (fullfield) THEN

        IF (zonal_data) THEN

          CALL SCATTER_ZONAL_FIELD( LOCAL_FIELD,GLOBAL_FIELD,
     &                              lasize(2),global_y,
     &                              LEVELS,GRIDTYPE_CODE,
     &                              SCATTER_PE)

! Don't call swapbounds for ocean zonal fields which currently
! do not have halos

          IF ((GRIDTYPE_CODE .NE. ppx_ocn_uzonal) .AND.
     &        (GRIDTYPE_CODE .NE. ppx_ocn_tzonal)) THEN
            CALL SWAPBOUNDS(LOCAL_FIELD,1,lasize(2),0,Offy,LEVELS)
          ENDIF

        ELSE

          DO level=1,LEVELS

            CALL SCATTER_FIELD( LOCAL_FIELD(1,level) ,
     &                          GLOBAL_FIELD(1,level),
     &                          lasize(1),lasize(2),
     &                          global_x,global_y,
     &                          SCATTER_PE,GC_ALL_PROC_GROUP,
     &                          info)

          ENDDO

          CALL SWAPBOUNDS(LOCAL_FIELD,lasize(1),lasize(2),
     &                    Offx,Offy,LEVELS)

         ENDIF
       ELSE
! for subdomains, life is not so easy - we must explicitly
! calculate our own send and receive maps, and use GCG_RALLTOALLE
! to shift the data around.

! If the same arguments are used as were used in the last call
! to this routine, we can just use the previously calculated
! send and receive maps, otherwise we need to calculate new maps

        IF (.NOT. (
     &    (LOCAL_SIZE .EQ. old_LOCAL_SIZE) .AND.
     &    (GLOBAL_SIZE .EQ. old_GLOBAL_SIZE) .AND.
     &    (GLOBAL_NORTH .EQ. old_GLOBAL_NORTH) .AND.
     &    (GLOBAL_EAST_IN .EQ. old_GLOBAL_EAST_IN) .AND.
     &    (GLOBAL_SOUTH .EQ. old_GLOBAL_SOUTH) .AND.
     &    (GLOBAL_WEST .EQ. old_GLOBAL_WEST) .AND.
     &    (GRIDTYPE_CODE .EQ. old_GRIDTYPE_CODE) .AND.
     &    (SCATTER_PE .EQ. old_SCATTER_PE) .AND.
     &    (current_decomp_type .EQ. old_current_decomp_type ))) THEN

          old_LOCAL_SIZE=LOCAL_SIZE
          old_GLOBAL_SIZE=GLOBAL_SIZE
          old_GLOBAL_NORTH=GLOBAL_NORTH
          old_GLOBAL_EAST_IN=GLOBAL_EAST_IN
          old_GLOBAL_SOUTH=GLOBAL_SOUTH
          old_GLOBAL_WEST=GLOBAL_WEST
          old_GRIDTYPE_CODE=GRIDTYPE_CODE
          old_SCATTER_PE=SCATTER_PE
          old_current_decomp_type=current_decomp_type

! Find out what the boundaries of the subdomain area

          CALL GLOBAL_TO_LOCAL_RC(GRIDTYPE_CODE,
     &                            GLOBAL_WEST,GLOBAL_NORTH,
     &                            proc_topleft_x,proc_topleft_y,
     &                            dummy1,dummy2)
          CALL GLOBAL_TO_LOCAL_RC(GRIDTYPE_CODE,
     &                            GLOBAL_EAST,GLOBAL_SOUTH,
     &                            proc_botright_x,proc_botright_y,
     &                            dummy1,dummy2)

! Ensure that the processor x co-ords are such that the botright_x is
! always greater than (or equal to) top_left_x.
          IF (wrap) proc_botright_x=gridsize(1)+proc_botright_x

! wrap_special is set to true if there is a wrap around which starts
! and ends on the same processor. This case requires extra work as
! the processor in question
          IF (wrap .AND. (proc_topleft_x+gridsize(1) .EQ.
     &                    proc_botright_x)) THEN
            wrap_special=.TRUE.
          ELSE
            wrap_special=.FALSE.
          ENDIF

          n_sends=0
          n_recvs=0

          DO procy=proc_topleft_y,proc_botright_y
            DO procx=proc_topleft_x,proc_botright_x

              eff_procx=MOD(procx,gridsize(1))
              procid=eff_procx+procy*gridsize(1)

              CALL GLOBAL_TO_LOCAL_SUBDOMAIN(
     &          .TRUE.,.TRUE.,
     &          GRIDTYPE_CODE,procid,
     &          GLOBAL_NORTH,GLOBAL_EAST,
     &          GLOBAL_SOUTH,GLOBAL_WEST,
     &          local_ystart,local_xend,
     &          local_yend  ,local_xstart)

! Calculate the shape of the arrays, and where to start sending/
! receiving data, and how many rows to send

              local_start_row=1
              nrows_to_send=local_yend-local_ystart+1

              global_start_row=g_datastart(2,procid)+local_ystart-Offy-
     &                         GLOBAL_NORTH
              global_row_length=GLOBAL_EAST-GLOBAL_WEST+1

! Calculate the following variables:
! local_row_length : the X dimension size of the local array
! local_send_offx  : the offset into each row to start sending from
! sendsize_x       : the number of points on each row to send
! The calculation of these numbers is different for processors
! at the start and end of a wrap_special case

              IF (wrap_special .AND. procx .EQ. proc_topleft_x) THEN
                local_row_length=g_lasize(1,procid)+local_xend-
     &                           local_xstart-2*Offx+1
                local_start_col=1
                sendsize_x=g_lasize(1,procid)-local_xstart
                global_start_col=1

              ELSEIF (wrap_special .AND. procx .EQ. proc_botright_x)
     &        THEN
                local_row_length=g_lasize(1,procid)+local_xend-
     &                           local_xstart-2*Offx+1
                local_start_col=local_row_length-local_xend+Offx+1
                sendsize_x=local_xend-Offx
                global_start_col=global_row_length-sendsize_x+1

              ELSE
                local_row_length=local_xend-local_xstart+1
                local_start_col=1
                sendsize_x=local_xend-local_xstart+1
                global_start_col=local_xstart-(Offx+1)+
     &                           g_datastart(1,procid)-GLOBAL_WEST+1
              ENDIF

              IF (global_start_col .LT. 0) THEN
! Wrapped around field, but this processor is not start or end
! processor
                global_start_col=global_start_col+glsize(1)
              ENDIF

! Now we can set up the send and receive map entries for the data on
! this processor

              IF (mype .EQ. procid) THEN  ! I need to receive some data

                  n_recvs=n_recvs+1

                receive_map(R_SOURCE_PE,n_recvs) = SCATTER_PE
                receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_recvs) =
     &            (local_start_row-1)*local_row_length +
     &            local_start_col
                receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_recvs) =
     &            nrows_to_send
                receive_map(R_STRIDE_IN_RECV_ARRAY,n_recvs) =
     &            local_row_length
                receive_map(R_ELEMENT_LENGTH,n_recvs) = sendsize_x
                receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_recvs) =
     &            (global_start_row-1)*global_row_length +
     &            global_start_col
                receive_map(R_STRIDE_IN_SEND_ARRAY,n_recvs) =
     &            global_row_length

              ENDIF ! if I'm receiving data

              IF (mype .EQ. SCATTER_PE) THEN ! I need to send data

                n_sends=n_sends+1

                send_map(S_DESTINATION_PE,n_sends) = procid
                send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,n_sends) =
     &            (global_start_row-1)*global_row_length +
     &            global_start_col
                send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,n_sends) =
     &            nrows_to_send
                send_map(S_STRIDE_IN_SEND_ARRAY,n_sends) =
     &            global_row_length
                send_map(S_ELEMENT_LENGTH,n_sends) = sendsize_x
                send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,n_sends) =
     &            (local_start_row-1)*local_row_length +
     &            local_start_col
                send_map(S_STRIDE_IN_RECV_ARRAY,n_sends) =
     &            local_row_length

              ENDIF ! if I'm sending data

            ENDDO ! procx : loop along processor row

          ENDDO ! procy : loop down processor column

        ENDIF ! if I need to recalculate my send/receive maps

! Send / receive the data using GCG_RALLTOALLE

        DO level=1,LEVELS

          flag=0  ! This is currently ignored at GCG v1.1

          CALL GC_SETOPT(GC_SHM_DIR,GC_SHM_GET,info)  ! set as scatter

          info=GC_NONE

          CALL GCG_RALLTOALLE(
     &      GLOBAL_FIELD(1,level)  ,
     &      send_map    , n_sends  ,GLOBAL_SIZE  ,
     &      LOCAL_FIELD(1,level) ,
     &      receive_map , n_recvs , LOCAL_SIZE ,
     &      GC_ALL_PROC_GROUP , flag, info)

        ENDDO

      ENDIF ! if this is a full or extracted field

 9999 CONTINUE

      RETURN
      END

