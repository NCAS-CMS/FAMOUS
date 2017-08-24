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
CLL  Subroutine DIAGDESC -----------------------------------------------
CLL
CLL  Purpose: Prints a formatted diagnostic description using the name
CLL           of a diagnostic plus it's PPXREF and STASH record.  Gives
CLL           a hardcopy record of the diagnostics included in a run.
CLL
CLL  Tested under compiler:   cft77
CLL  Tested under OS version: UNICOS 6.1
CLL
CLL  Author:   T.Johns            Date:           14 January 1992
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.1  05/02/93  Correct minor bug in printout for climate mean tag.
CLL                  Print out pseudo-level information.
CLL  3.1   3/02/93 : added comdeck CHSUNITS to define NUNITS for i/o.
CLL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
CLL                  1.Removes the limit on primary STASH item numbers.
CLL                  2.Removes the assumption that (section,item)
CLL                    defines the sub-model.
CLL                  3.Thus allows for user-prognostics.
CLL   3.4  13/01/94  Replace hardwired gridcodes by ppx_ parameters, and
CLL                  cover all options.   T. Johns
!     4.4  02/12/96 Add daily mean timeseries R. A. Stratton.    
CLL
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered: C401
CLL
CLL  Project task: C4
CLL
CLL  External documentation:
CLL    Unified Model Doc Paper C4 - Storage handling and diagnostic
CLL                                 system (STASH)
CLLEND --------------------------------------------------------------
C
C*L  Interface and arguments: ------------------------------------------
C
      SUBROUTINE DIAGDESC(seqno,name,stlist,ppxref,
     &           stash_levels,num_stash_levels,num_level_lists,
     &           stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists,
     &           sttabl,nsttims,nsttabl,
     &           stash_series,stash_series_rec_len,stash_series_len,
     &           stash_series_index,stash_ser_index_size)
C
      IMPLICIT NONE
C
      CHARACTER*36
     *    name                                  ! IN  diagnostic name
      INTEGER
     *    seqno,                                ! IN  sequence number
     *    stlist(*),                            ! IN  STASHlist record
     *    ppxref(*)                             ! IN  PPXREF record
C
C STASH levels list information
      INTEGER
     &       num_stash_levels                 ! IN Max levels in a list
     &,      num_level_lists                  ! IN Number of lists
     &,      stash_levels(num_stash_levels+1,num_level_lists)
C STASH pseudo-levels list information
      INTEGER
     &       num_stash_pseudo                 ! IN Max ps-levs in a list
     &,      num_pseudo_lists                 ! IN No of ps-lev lists
     &,      stash_pseudo_levels(num_stash_pseudo+1,num_pseudo_lists)
C STASH time list information
      INTEGER
     &       nsttims                          ! IN Max times in a list
     &,      nsttabl                          ! IN Number of lists
     &,      sttabl(nsttims,nsttabl)
C STASH timeseries information
      INTEGER
     &       stash_series_len                 ! IN Total no of records
     &,      stash_series_rec_len             ! IN Length of each record
     &,      stash_series(stash_series_rec_len,stash_series_len)
C                                             ! IN array of records
     &,      stash_ser_index_size             ! IN No of index records
     &,      stash_series_index(2,stash_ser_index_size)
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
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.5    07/04/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
!
! Declarations:
!
!  1. Internal model and submodel dump partition identifiers - fixed
!     for all experiments.
!
! Description:
!    Hold parameters defining internal model identifiers and submodel
!    data partition (ie main D1 data array and consequent dump), both
!    short and long form.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.3    26/10/93   M. Carter. Part of an extensive mod that:
!                    1.Removes the limit on primary STASH item numbers.
!                    2.Removes the assumption that (section,item)
!                      defines the sub-model.
!                    3.Thus allows for user-prognostics.
!                    Add index to submodel home dump.
! 3.5    13/03/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
!
! Declarations:
!
!   Hold parameters defining internal model identifiers and submodel
!   data partition (ie main D1 data array and consequent dump), both
!   short and long form
      INTEGER
     *   A_IM,ATMOS_IM        ! Atmosphere internal model
     *  ,O_IM,OCEAN_IM        ! Ocean      internal model
     *  ,S_IM, SLAB_IM        ! Slab       internal model
     *  ,W_IM, WAVE_IM        ! Wave       internal model
     *  ,I_IM,SEAICE_IM       ! Sea-ice    internal model
     *  ,N_IM,NATMOS_IM       ! New dynamics (Charney-Phillips grid)
!                               atmosphere internal model
!
      PARAMETER(
     *   A_IM=1,ATMOS_IM=1       ! Atmosphere internal model
     *  ,O_IM=2,OCEAN_IM=2       ! Ocean      internal model
     *  ,S_IM=3, SLAB_IM=3       ! Slab       internal model
     *  ,W_IM=4, WAVE_IM=4       ! Wave       internal model
     *  ,I_IM=5,SEAICE_IM=5      ! Sea-ice    internal model
     *  ,N_IM=6,NATMOS_IM=6      ! New dynamics (Charney-Phillips grid)
!                                  atmosphere internal model
     *)
!
      INTEGER
     *   A_SM,ATMOS_SM        ! Atmosphere submodel partition
     *  ,O_SM,OCEAN_SM        ! Ocean      submodel partition
     *  ,W_SM, WAVE_SM        ! Wave       submodel partition
     *  ,N_SM,NATMOS_SM       ! New dynamics (Charney-Phillips grid)
!                                  atmosphere internal model
!
      PARAMETER(
     *   A_SM=1,ATMOS_SM=1    ! Atmosphere submodel partition
     *  ,O_SM=2,OCEAN_SM=2    ! Ocean      submodel partition
     *  ,W_SM=4, WAVE_SM=4    ! Wave       submodel partition
     *  ,N_SM=6,NATMOS_SM=6   ! New dynamics (Charney-Phillips grid)
!                                  atmosphere internal model
     *)
!
C

!
!  2. Maximum internal model/submodel array sizes for this version.
!
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 3.5    13/07/95   Original code. D.M. Goddard
! 4.0     3/11/95   Reduce max internal model, submodel from 10 to 4
!                   to save space in model. At 4.0 the max no of 
!                   supported models is 3, 1 slot is reserved for
!                   expansion. Rick Rawlins.         
!  4.1  21/02/96  Wave model introduced as 4th sub-model.  RTHBarnes
!
! Declarations:
!
!
!  1. Maximum internal model/submodel array sizes for this version.
!
      INTEGER
     * N_INTERNAL_MODEL_MAX      ! Max no. of internal models
     *,N_SUBMODEL_PARTITION_MAX  ! Max no. of submodel dump partitions
     *,INTERNAL_ID_MAX           ! Max value of internal model id
     *,SUBMODEL_ID_MAX           ! Max value of submodel dump id

      PARAMETER(
     * N_INTERNAL_MODEL_MAX=4,                                          
     * N_SUBMODEL_PARTITION_MAX=4,
     * INTERNAL_ID_MAX=N_INTERNAL_MODEL_MAX,
     * SUBMODEL_ID_MAX=N_SUBMODEL_PARTITION_MAX)
!
!  3. Lists of internal models and their submodel dump partitions -
!     initialised by the user interface - experiment specific.
      INTEGER
     * N_INTERNAL_MODEL          ! No. of internal models
     *,N_SUBMODEL_PARTITION      ! No. of submodel partitions
     *,INTERNAL_MODEL_LIST(N_INTERNAL_MODEL_MAX) ! Internal models
     *,SUBMODEL_FOR_IM    (N_INTERNAL_MODEL_MAX) ! Submodel identifier
     *                           ! for each internal model in list
     &,SUBMODEL_FOR_SM(N_INTERNAL_MODEL_MAX) ! Submodel number for
!                                  each submodel id
!
! Namelist for information in 3.
      NAMELIST/NSUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION
     *,INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM
!
!  4. Lists calculated in model from user interface supplied arrays -
!     - experiment specific.
      INTEGER
     * N_INTERNAL_FOR_SM(SUBMODEL_ID_MAX)  ! No of internal models in
!              each submodel partition indexed by sm identifier
     *,SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION_MAX)    ! List of
!              submodel partition identifiers
     *,SUBMODEL_PARTITION_INDEX(INTERNAL_ID_MAX)  ! Submodel partition
!              identifier indexed by internal model identifier
     *,INTERNAL_MODEL_INDEX(INTERNAL_ID_MAX)      ! Sequence number of
!              internal model indexed by internal model identifier:
!              required to map from id to STASH internal model sequence
      LOGICAL
     * LAST_IM_IN_SM(INTERNAL_ID_MAX)      ! Last internal model within
!                                a submodel partition if .TRUE.,
!                                indexed by internal model id.
! Common block for information in 3. and 4.
      COMMON/SUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION,
     *     INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM,SUBMODEL_FOR_SM,
     *     N_INTERNAL_FOR_SM,SUBMODEL_PARTITION_LIST,
     *     SUBMODEL_PARTITION_INDEX,
     *     INTERNAL_MODEL_INDEX,
     *     LAST_IM_IN_SM

!
!  5. Time information specifying coupling frequencies between internal
!     models and submodels, and multipliers, indexed by sequence of
!     internal models and submodels (ie left to right along node tree).
!     {Not required at this release}.
!
! Namelists for information in 5. {Not required at this release}
!
!
!  6. Lists of coupling nodes defining coupling frequencies between
!     internal models and between submodel partitions. (Not defined
!     yet at this release).
!CALL CNODE
!
!  7. Variables dealing with general coupling switches at the control
!     level. {These will require revision at the next release when
!     coupling between internal models is dealt with more generally.
!     Logicals below are set in routine SETGRCTL.}

      LOGICAL
     * new_im   ! new internal model next group of timesteps if .true.
     *,new_sm   ! new submodel dump  next group of timesteps if .true.

      COMMON/CSUBMGRP/new_im,new_sm

      INTEGER SUBMODEL_IDENT
      COMMON/SUBMODID/SUBMODEL_IDENT                                    
C
C
C*L --------------------- Comdeck: CHSUNITS -------------------------
CLL
CLL Purpose: COMDECK defining the number of i/o units
CLL
CLL  Author : R A Stratton
CLL
CLL  Model            Modification history:
CLL version  date
CLL   3.1  03/02/93   Introduced at version 3.1
CLL   4.1  21/02/96   Increase no.of i/o units to accommodate wave 
CLL                   sub-model.  RTHBarnes.
CLL
CLL Project task:
CLL
CLL  Documentation:  Unified Model Documentation Paper
CLL                  H- History Bricks
CLL
CLLEND---------------------------------------------------------------
C
C*L Type declarations
C
      INTEGER NUNITS          ! No. of I/O units
      INTEGER NUNITS_LEN      ! length of most unit no arrays
!
!     These values must be consistent with OUTFILE_S, OUTFILE_L
!     and OUTFILE_E in comdeck VERSION.
      PARAMETER(NUNITS=149)
      PARAMETER(NUNITS_LEN=NUNITS-19)
C
C   The above parameter statements must not be altered without 
C   considering the effect on the following HISTORY COMDECKs 
C    CHISTO, CLFHIST and IHISTO.    
C  This comdeck must always preceed the above history file comdecks.
C   New file environment variable names may need to be added to     
C    CLFHIST and/or CENVIRDT (usually both) depending on manner of I/O.
CLL  Comdeck: CCONTROL -------------------------------------------------
CLL
CLL  Purpose: COMMON block for top level switches and 2nd level switches
CLL           needed by the top level (C0) and 2nd level routines, but
CLL           not held in the history COMMON block.
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.1    8/02/93 : Changed 99 to NUNITS for i/o. Note this comdeck
CLL                    must always be called after CHSUNITS so that
CLL                    NUNITS is defined.
CLL  3.1  15/02/93  Add L_Z0_OROG orographic roughness switch. R.Barnes.
CLL  3.3  09/07/93  Add L_CONVECT =F to add saved convection increments,
CLL                               =T to call conv.scheme.  R.T.H.Barnes.
CLL   3.3    13/12/93   Insert switches for half timestep
CLL                     dynamics. A.S.Lawless
CLL   3.4 23/08/94  Add switch for local -ve q correction R.A.Stratton.
CLL   3.4    16/06/94   COMMON block DEFLOGIC inserted - declares
CLL                     logical switches for control and other
CLL                     purposes - most of these have replaced *DEFs
CLL                                                     S.J.Swarbrick
CLL   3.4    1/8/94  Add control for assimilation mode S Bell
CLL  3.5  28/03/95  Sub-Models stage 1: revise History and Control file
CLL                 contents. Control expanded and subdivided by
CLL                 overall, generic and specific categories. RTHBarnes
CCL  4.1  23/02/96  Extend for new wave sub-model.  RTHBarnes.
CLL
CLL Logical components covered :
CLL
CLL External documentation: Unified Model documentation paper No
CLL                         Version
CLL
CLLEND ---------------------------------------------------------------

! ----------------------- Comdeck: CNTLALL  ----------------------------
! Description: COMDECK defining Control variables for the
!              model overall.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.0  25/10/95  Add user switch CONTROL_RESUBMIT. RTHBarnes
!  4.4  28/07/97  Add user switch LCLIMREALYR. M Gallani
!  4.4  11/10/97  Add logical switch L_AO_D1_MEMORY. D. Robinson. 
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      INTEGER
     &        MODEL_BASIS_TIME(6),! Array holding original data time
!                                 ! (prior to assimilation)
     &        MODEL_ANALYSIS_HRS, ! Model analysis time in hours since
!                                 ! Basis Time
     &        MODEL_HRS_PER_GROUP,! No. of hours in coupling period
     &        NCPU,               ! No of CPUs assigned to the program
     &        ANCIL_REFTIME(6),   ! Ref. time for updating ancillaries
     &        FT_PLOTSEL(60:69),  ! interval for plotting pp file
     &        RUN_TARGET_END(6),   ! Target end time for this run
     &        RUN_RESUBMIT_INC(6), ! Increment to be added on each
!                                  ! resubmission of the job.
     &   PP_LEN2_LOOK(20:NUNITS),  ! Number of field headers reserved
!                                  !  for non-mean PPfiles on each unit
     &   PP_PACK_CODE(20:NUNITS),  ! Internally defined PP packing code
     &   FT_STEPS(20:NUNITS),   ! Frequency of initialisation of FTunit
     &   FT_FIRSTSTEP(20:NUNITS)   ! ... starting at step number ..
!
      LOGICAL
     &       LATMOSNEXT,LOCEANNEXT,  ! Flags to select atmosphere/ocean
     &       LPP,                    ! Activate PPCTL
     &       LPP_SELECT(20:NUNITS),  ! Activate PP init'sation on unit
     &       LDUMP,                  ! Activate DUMPCTL
     &       LMEAN,                  ! Activate MEANCTL
     &       LHISTORY,               ! Update TEMP history file
     &       LPRINT,                 ! Activate PRINTCTL
     &       LINTERFACE,             ! Activate GEN_INTF
     &       LEXIT,                  ! Activate EXITCHEK
     &       LJOBRELEASE,            ! Activate JOBCTL
     &       LMEANPR(4),             ! Select printed diags from means
     &       LANCILLARY,             ! Activate UP_ANCIL
     &       LBOUNDARY,              ! Activate UP_BOUND
     &       LASSIMILATION,          ! Activate assimilation
     &       LCAL360,                ! 360-day calendar
     &       LTIMER                  ! Activate detailed TIMER routine
     &      ,L_AO_D1_MEMORY  ! T : D1 copied to memory for AO coupling
     &      ,LCLIMREALYR             ! Real-period climate means

      CHARACTER*4  EXPT_ID          ! Unique alphanumeric serial number
!                                   ! associated with model
!                                   ! (Non-Operational expts)
!                                   !
!                                   ! Operational run name
!                                   ! (Operational expts)
      CHARACTER*8  EXPT_ALIAS       ! Non unique user defined expt name
      CHARACTER*1  JOB_ID           ! Unique alphanumeric job identifier
!                                   ! used for networking
      CHARACTER*4  EXPT_ID_IN       ! Experiment ID of driving model if
!                                   ! limited-area run
      CHARACTER*4  JOB_ID_IN        ! Job ID of driving model if
!                                   ! limited-area run
      CHARACTER*14 MODEL_STATUS     ! Operational or NonOperational
      CHARACTER*14 MODEL_ASSIM_MODE ! Atmosphere,Ocean,Coupled or None
      CHARACTER*17 TIME_CONVENTION  ! Relative, Timestep, Absolute_long,
!                                    Absolute_standard or Absolute_short
      CHARACTER*1  FT_WSSEND(60:69) ! "Y" if file to be sent to HP
!
      CHARACTER*1 TYPE_LETTER_1(20:NUNITS) ! File type letter #1
      CHARACTER*1 TYPE_LETTER_2(20:NUNITS) ! File type letter #2
      CHARACTER*1 TYPE_LETTER_3(20:NUNITS) ! File type letter #3
!
      CHARACTER*1  FT_INPUT (20:NUNITS) ! "Y" if input file on unit
      CHARACTER*1  FT_OUTPUT(20:NUNITS) ! "Y" if output file on unit
      CHARACTER*1  FT_SELECT(20:NUNITS) ! "Y" if file selected for post
!                                          processing request.
      CHARACTER*1  FT_ARCHSEL(20:NUNITS) ! "Y" if file to be archived.
!
      CHARACTER*10 RUN_ASSIM_MODE      ! cf MODEL_ASSIM_MODE (Oper use)
      CHARACTER*1  CONTROL_RESUBMIT    ! User flag for auto resubmit

      NAMELIST / NLSTCALL /
     & MODEL_BASIS_TIME, MODEL_ANALYSIS_HRS,
     & MODEL_HRS_PER_GROUP,
     & NCPU, ANCIL_REFTIME, FT_PLOTSEL, RUN_TARGET_END,
     & RUN_RESUBMIT_INC, PP_LEN2_LOOK, PP_PACK_CODE,
     & FT_STEPS, FT_FIRSTSTEP,
     & LATMOSNEXT, LOCEANNEXT, LPP, LPP_SELECT, LDUMP, LMEAN,
     & LHISTORY, LPRINT, LINTERFACE, LEXIT, LJOBRELEASE,
     & LMEANPR, LANCILLARY, LBOUNDARY, LASSIMILATION,
     & LCAL360, LTIMER, L_AO_D1_MEMORY,
     & LCLIMREALYR,
     & EXPT_ID, JOB_ID, EXPT_ID_IN, JOB_ID_IN,
     & EXPT_ALIAS, MODEL_STATUS, MODEL_ASSIM_MODE,
     & TIME_CONVENTION, FT_WSSEND,
     & TYPE_LETTER_1, TYPE_LETTER_2, TYPE_LETTER_3,
     & FT_INPUT, FT_OUTPUT, FT_SELECT, FT_ARCHSEL,
     & RUN_ASSIM_MODE, CONTROL_RESUBMIT

      COMMON / CNTLCALL /
     & MODEL_BASIS_TIME, MODEL_ANALYSIS_HRS,
     & MODEL_HRS_PER_GROUP,
     & NCPU, ANCIL_REFTIME, FT_PLOTSEL, RUN_TARGET_END,
     & RUN_RESUBMIT_INC, PP_LEN2_LOOK, PP_PACK_CODE,
     & FT_STEPS, FT_FIRSTSTEP,
     & LATMOSNEXT, LOCEANNEXT, LPP, LPP_SELECT, LDUMP, LMEAN,
     & LHISTORY, LPRINT, LINTERFACE, LEXIT, LJOBRELEASE,
     & LMEANPR, LANCILLARY, LBOUNDARY, LASSIMILATION,
     & LCAL360, LTIMER, L_AO_D1_MEMORY,
     & LCLIMREALYR,
     & EXPT_ID, JOB_ID, EXPT_ID_IN, JOB_ID_IN,
     & EXPT_ALIAS, MODEL_STATUS, MODEL_ASSIM_MODE,
     & TIME_CONVENTION, FT_WSSEND,
     & TYPE_LETTER_1, TYPE_LETTER_2, TYPE_LETTER_3,
     & FT_INPUT, FT_OUTPUT, FT_SELECT, FT_ARCHSEL,
     & RUN_ASSIM_MODE, CONTROL_RESUBMIT
! ----------------------- Comdeck: CNTLGEN  ----------------------------
! Description: COMDECK defining Control variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  28/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.0   3/11/95  Move character array MEANSim to the end of the 
!                 common block to ensure that it starts correctly on a
!                 word boundary. [No problem is apparent on the Cray  
!                 if N_INTERNAL_MODEL_MAX is an even no.]
!                 Rick Rawlins                                    
!  4.1  03/04/96  Add new array DUMP_PACKim. D. Robinson
!  4.5  10/11/98  Increase number of dumps allowed at irregular 
!                 timesteps from 10 to 40: Move lengths into
!                 CNTLGEN. R Rawlins
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      INTEGER
     & DUMPTIMES_LEN1   ! Max no. of irregular times for dumps
     &,PRINTFREQ_LEN1   ! No. of areas of zonal mean prints
     &,MEANFREQ_LEN1    ! No. of time intervals for climate meaning
     &,JOBREL_LEN1      ! Max no. of irregular times for job release

      PARAMETER(
     & DUMPTIMES_LEN1 = 40 
     &,PRINTFREQ_LEN1 = 5
     &,MEANFREQ_LEN1  = 4
     &,JOBREL_LEN1    = 10
     &) 
      INTEGER
     & STEPS_PER_PERIODim(N_INTERNAL_MODEL_MAX)
     &,SECS_PER_PERIODim(N_INTERNAL_MODEL_MAX)
     &,EXITFREQim(N_INTERNAL_MODEL_MAX)  ! Number of advection
!                               timesteps between checks for model exit
     &,DUMPFREQim(N_INTERNAL_MODEL_MAX)  ! Number of steps between
!                                              atmosphere restart dumps
     &,ARCHDUMP_FREQim(N_INTERNAL_MODEL_MAX)  ! Archiving frequency
!                                                   for atmos dumps
     &,DUMPTIMESim(DUMPTIMES_LEN1,N_INTERNAL_MODEL_MAX) ! Timesteps 
!            (from start of run) at which restart dumps are written
     &,MEANFREQim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX) !Indicators 
!             for mean dump frequency
     &,MEANARCHim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX) !Switches 
!             for mean dump arch.
     &,PPSELECTim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX) !PP field 
!             selectors 
     &,ARCHPPSELim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)!Switches 
!             for pp field archive
     &,PLOTSELim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)  !Switches 
!             for chart plotting
     &,PP_LEN2_MEANim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX) !Number of
!             field headers to reserve for internal model mean PPfiles
     &,MEAN_REFTIMEim(6,N_INTERNAL_MODEL_MAX)  ! Reference time for
!                                                production of means
     &,PRINTFREQim(PRINTFREQ_LEN1,N_INTERNAL_MODEL_MAX) ! Indicators   
!             of zonal mean print frequency
     &,JOBREL_STEPim(JOBREL_LEN1,N_INTERNAL_MODEL_MAX)  ! Step numbers 
!             at which to release user-specified scripts 
     &,ARCHDUMP_OFFSETim(N_INTERNAL_MODEL_MAX)!Offset for dump archiving
     &,FT_MEANim(N_INTERNAL_MODEL_MAX)     ! Unit reserved for mean PPs
     &,DUMP_PACKim(N_INTERNAL_MODEL_MAX)  ! Packing indicator for dumps
      CHARACTER*1  MEANWSim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX) 
!                                    "Y" if mean file to be sent to HP
      LOGICAL    LLBOUTim(N_INTERNAL_MODEL_MAX)  ! Lateral b.c.'s
     &          ,LANCILim(N_INTERNAL_MODEL_MAX)  ! Ancillary files
C
      NAMELIST / NLSTCGEN /
     & STEPS_PER_PERIODim, SECS_PER_PERIODim,
     & EXITFREQim, DUMPFREQim,
     & ARCHDUMP_FREQim, DUMPTIMESim, PPSELECTim, PLOTSELim,
     & ARCHPPSELim, MEANARCHim, MEANFREQim, MEAN_REFTIMEim,
     & PRINTFREQim,  JOBREL_STEPim, ARCHDUMP_OFFSETim, PP_LEN2_MEANim,
     & FT_MEANim,
     & DUMP_PACKim,
     & MEANWSim, LLBOUTim, LANCILim
C
      COMMON / CNTLCGEN /
     & STEPS_PER_PERIODim, SECS_PER_PERIODim,
     & EXITFREQim, DUMPFREQim,
     & ARCHDUMP_FREQim, DUMPTIMESim, PPSELECTim, PLOTSELim,
     & ARCHPPSELim, MEANARCHim, MEANFREQim, MEAN_REFTIMEim,
     & PRINTFREQim,  JOBREL_STEPim, ARCHDUMP_OFFSETim, PP_LEN2_MEANim,
     & FT_MEANim,
     & DUMP_PACKim,
     &  LLBOUTim, LANCILim,
     &  MEANWSim


! ----------------------- Comdeck: CNTLATM  ----------------------------
! Description: COMDECK defining Control variables for the Atmosphere
!              internal model, and its runtime constants.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  29/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.0  17/08/95  New variables for long physics timestep,
!                 and tracer advection.  RTHBarnes.
!  4.0  7/11/95  Logical switches for convective momentum transports
!                and CAPE closure added to namelists. Pete Inness.
!  4.1  8/5/96   Logical switch for rapidly mixing boundary layer
!  4.1 28/05/96  New control switches added for Sulphur Chemistry
!                Cycle and Hydrology Schemes. D Robinson & D. Goddard.
!  4.3  18/3/97  Flag to indicate if the HadCM2 approximate treatment
!                of sulphate aerosol is being used.       William Ingram
!  4.3 14/04/97  New control switch L_OLD_PWTS for old polar geometric
!                weights, needed for HADCM2.    T Johns.
!  4.3 03/02/97  Logical L_MIXLEN for mixing in the boundary layer
!                  S Jackson
!  4.3 03/02/97  Logical switches L_XSCOMP and L_SDXS for convection
!                scheme  S Jackson
!  4.4 4/7/97    Add control for IAU  Chris Jones/Stuart Bell 
!  4.4 05/09/97  Logical LFLUX_RESET to indicate when net flux field
!                needs initialising to 0. S.D. Mullerworth
!  4.4 17/09/97  Logical switches L_CCW and L_3D_CCA added to enable
!                use of anvil package/3D conv. cloud amt. J.Gregory
!  4.4 08/09/97 Logical switches L_LSPICE, L_BL_LSPICE and L_LSPICE_BDY
!               for mixed phase precipitation.
!                                                       D.Wilson
!  4.4 10/10/97  Logical switch L_SNOW_ALBEDO.  Richard Essery   
!  4.4 10/10/97  Logical switches L_VEG_FRACS, L_TRIFFID, L_PHENOL, 
!                L_NRUN_MID_TRIF and L_TRIF_EQ for veg.  Richard Betts
!  4.5   1/07/98  Add logicals to control interactive CO2 use. C.D.Jones
!   4.5  28/04/98  Add logicals for NH3 and SOOT variables and emiss
!                                                           M Woodage
!  4.5 21/08/98  Logical switch l_ssice_albedo.  Jonathan Gregory
!  4.5 20/05/98  Logical switch L_NEG_TSTAR for negative surface 
!                temperature error check.  Richard Betts  
!  4.5 19/11/98  Add PHENOL_PERIOD and TRIFFID_PERIOD, moved from
!                NLSTCATM.  Richard Betts  
!  4.5 19/05/98  Logical switch L_PHASE_LIM for HADAM3 physics in
!                optimised convection scheme.       Julie Gregory
!  4.5 26/06/98  Logical switches L_RHCPT, L_CLD_AREA for new
!                Section 9 parametrizations.             S. Cusack
!  4.5 22/10/98  Remove redundant switch LMULTIL_HYDROL
!                Author D.M. Goddard
!  4.5 05/06/98  Add Logical switch L_VINT_TP.  D Robinson.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Parameter declarations
      INTEGER MAXSECTS            ! Max. no. of code sections
      PARAMETER (MAXSECTS=99)
!
!   Type declarations
!
      INTEGER H_SWBANDS,      ! Number of shortwave radiation bands
     &        H_LWBANDS,      ! Number of longwave radiation bands
     &        A_ADJSTEPS,         ! No. of adjustment timesteps per
!                                 ! advection step
     &        A_SW_RADSTEP,       ! Number of advection steps per
!                                 ! shortwave radiation step
     &        A_LW_RADSTEP,       ! Number of advection steps per
!                                 ! longwave radiation step
     &        A_SW_SEGMENTS,      ! No of batches used in shortwave code
     &        A_LW_SEGMENTS,      ! No of batches used in longwave code
     &        A_CONV_STEP,        ! No of advection timesteps between
!                                 ! calls to convection scheme
     &        A_CONVECT_SEGMENTS, ! No of batches in convection code
     &        A_NSET_FILTER,      ! No of advection steps after which
!                                 ! filtering wavenumber checked
     &        A_ENERGYSTEPS,      ! Number of advection steps after
!                                 ! which energy adjustment performed
     &        A_ASSIM_START_HR,   ! Time at which data assimilation
!                                 ! starts (Hours after Basis Time)
     &        A_ASSIM_END_HR      ! Time at which data assimilation
!                                 ! ends (Hours after Basis Time)
     &       ,T_IAU_START, T_IAU_END  ! IAU before and after
!
     &       ,A_SWEEPS_DYN ! No.of dynamics sweeps per physics timestep
     &       ,CALL_CHEM_FREQ ! Frequency of calls to CHEM_CTL
     &       ,PHENOL_PERIOD ! Update frequency for leaf phenology (days)
     &       ,TRIFFID_PERIOD ! Update frequency for TRIFFID (days)

      LOGICAL
     1       L_SW_RADIATE,           ! Activate SW radiation
     2       L_LW_RADIATE,           ! Activate LW radiation
     3       LADD_RADINCS,           ! Both SW and LW radiation active
     &       L_H2_SULPH,             ! HadCM2 approximate sulphate on
     4       L_SET_FILTER,           ! Recalculate filtering in dynamics
     5       LDAY,                   ! End-of-day
     6       LEXPAND_OZONE,          ! Convert zonal mean ozone to field
     6       L_COMPRESS_LAND,        ! Compress land points in physics
     *       L_NEG_THETA,            ! Test for -ve theta in dynamics
     *       L_NEG_PSTAR,            ! Test for -ve P* in dynamics
     *       L_NEG_QT,               ! Test for -ve QT over layer
     &       L_NEG_TSTAR,           ! Test for -ve surface temperature.
     *       L_FIELD_FLT,            ! Apply field filtering (atmos)
     *       L_SUPERBEE,             ! Superbee(T) or Van Leer(F)
     *                               ! limiter for tracer advection
     7       LENERGY,                ! Recalculate energy drift
     &       LFLUX_RESET,            ! Reset net flux field
     *       L_Z0_OROG,              ! T to use orog.roughness code
     *       L_RMBL,                 ! T to use rapid mixing BL code
     &       L_MIXLEN,               ! T to make mixing length above BL
     &                               ! top independent of BL depth
     *       L_CONVECT               ! T call conv.scheme, F add incrs.
     *      ,L_HALF_TIMESTEP_DYN     ! T if wind threshold exceeded
     *      ,L_HALF_TIMESTEP_DIV     ! T if diverg. threshold exceeded
     &      ,L_QT_POS_LOCAL          ! Apply -ve q correction locally.
     &      ,L_TRACER_THETAL_QT      ! T if using tracer advection
!                                    !  for thetal & qt
     &      ,L_3DVAR_BG,L_AC,L_3DVAR,L_4DVAR !Switches for assm mode
     &      ,L_OLD_PWTS              ! T if using old polar weights
     &      ,L_VINT_TP               ! T: Use V_INT_TP to output Temp
                                     ! on model levels.
C
      LOGICAL                       ! Logical switches for:
     &   LFROUDE        ,           !  Limit max grav wave amp
     &   LGWLINP        ,           !  Linear grav wave stress prof
     &   LLINTS         ,           !  Linear TS approx
     &   LWHITBROM      ,           !  White & Bromley terms
     &   LEMCORR        ,           !  Energy & mass correction
     &   LMICROPHY      ,           !  Microphysics in sw rad scheme
     &   L_MURK         , !           :Total aerosol field
     &   L_MURK_ADVECT  , !           :Aerosol advection
     &   L_MURK_SOURCE  , !Bndry      :Aerosol source & sink terms
     &   L_MURK_BDRY    , !Layer      :UK Mes bndry model
     &   L_BL_TRACER_MIX, !model      :Bndry layer tracer mixing
     &   L_MOM,                     !  convective momentum mixing
     &   L_CAPE,                    !  CAPE closure for convection      
     &   L_SDXS,                    ! Convective excess from turbulent
!                                   !  fluctuations
     &   L_XSCOMP                   ! Environmental compensation for 
!                                   !  parcel excess
     &  ,L_3D_CCA                   ! Use 3D conv cloud amount
     &  ,L_CCW                      ! Rain not inc. in conv water path
     &  ,L_PHASE_LIM                ! Select 3B physics for A05_3C
     &  ,L_CLOUD_DEEP               ! Depth criterion applied for anvils
     &  ,LSINGLE_HYDROL             ! Single level hydrology
     &  ,LMOSES                     ! MOSES hydrology only
     &  ,L_SNOW_ALBEDO              ! Prognostic snow albedo
     &  ,l_ssice_albedo             ! Sea-ice albedo affected by snow
     &  ,L_SULPC_SO2                ! Sulphur Cycle : SO2 MMR included 
     &  ,L_SULPC_DMS                ! Sulphur Cycle : DMS MMR included
     &  ,L_SULPC_OZONE              ! Sulphur Cycle : Ozone included
     &  ,L_SO2_SURFEM               ! SO2 Surface Emissions
     &  ,L_SO2_HILEM                ! SO2 High Level Emissions
     &  ,L_SO2_NATEM                ! SO2 Natural Emissions
     &  ,L_DMS_EM                   ! DMS Emissions
     &  ,L_SULPC_NH3           ! S Cycle : NH3 included
     &  ,L_NH3_EM              ! S Cycle : NH3 emiss included
     &  ,L_SOOT                ! Soot included  
     &  ,L_SOOT_SUREM          ! surface Soot emiss included
     &  ,L_SOOT_HILEM          ! elevated Soot emiss included
     &  ,L_USE_SOOT_DIRECT     ! direct radiative effects of soot
     &  ,L_USE_SULPC_DIRECT   !\Use SO4 aerosol from sulphur cycle for
     &  ,L_USE_SULPC_INDIRECT_SW !direct/indirect effect in radiation,
     &  ,L_USE_SULPC_INDIRECT_LW !the latter for both SW and LW.
     &  ,L_CLIMAT_AEROSOL           ! Switch for climatological
!                                   ! aerosols in the radiation.        
     &  ,L_RHCPT                     ! controls the use of new RHcrit
                                     ! parametrization, vn 2B of Sec 9
     &  ,L_CLD_AREA                  ! controls cloud area parametriz.
     &  ,L_IAU_DIAG                 ! controls IAU diagnostics
     &  ,L_IAU                      ! controls IAU calls
     &  ,L_IAU_RAMP                 ! controls IAU weights
     &  ,L_VEG_FRACS                ! Switch for vegetation fractions
     &  ,L_TRIFFID                  ! Switch for interactive veg model
     &  ,L_PHENOL                   ! Switch for leaf phenology
     &  ,L_NRUN_MID_TRIF            ! Switch for starting NRUN mid-way 
C                                   ! through a TRIFFID calling period
     &  ,L_TRIF_EQ                  ! Switch for running TRIFFID in
C                                   ! equilibrium mode
     &      ,L_CO2_INTERACTIVE      ! interactive 3D CO2 field for
                                    !  use with carbon cycle model
     &      ,L_CO2_EMITS            ! include surface emissions
      LOGICAL L_LSPICE              ! New cloud/precip microphysics
     &,       L_BL_LSPICE           ! Full boundary layer treatment
!                                     with new cloud/precip scheme
     &,       L_LSPICE_BDY          ! QCF present in lateral boundaries
!
      CHARACTER*5 A_ASSIM_MODE     ! Switch for BG/AC/3DVAR/4DVAR assm
!
      CHARACTER*3 H_SECT(0:MAXSECTS) ! Array of code section versions
!
      NAMELIST / NLSTCATM /
     & L_VEG_FRACS, L_TRIFFID, L_PHENOL, L_NRUN_MID_TRIF, L_TRIF_EQ,
     & PHENOL_PERIOD, TRIFFID_PERIOD,
     & H_SWBANDS, H_LWBANDS, A_ADJSTEPS,
     & A_SW_RADSTEP, A_LW_RADSTEP, A_SW_SEGMENTS, A_LW_SEGMENTS,
     & A_CONV_STEP, A_CONVECT_SEGMENTS, A_NSET_FILTER, A_ENERGYSTEPS,
     & A_ASSIM_START_HR, A_ASSIM_END_HR, A_SWEEPS_DYN, L_H2_SULPH,
     & L_SW_RADIATE, L_LW_RADIATE, LADD_RADINCS, L_SET_FILTER,
     & LDAY, LEXPAND_OZONE, L_COMPRESS_LAND,
     & L_NEG_THETA, L_NEG_PSTAR, L_NEG_QT, L_NEG_TSTAR, L_FIELD_FLT,
     & L_SUPERBEE, LENERGY, L_Z0_OROG, L_RMBL, L_MIXLEN, L_CONVECT,
     & L_HALF_TIMESTEP_DYN, L_HALF_TIMESTEP_DIV, L_QT_POS_LOCAL,
     & L_TRACER_THETAL_QT,LFLUX_RESET,
!    & L_3DVAR_BG,L_AC,L_3DVAR,L_4DVAR, in COMMON but not in NAMELIST
     & L_OLD_PWTS, L_VINT_TP,
     & LFROUDE, LGWLINP, LLINTS, LWHITBROM, LEMCORR,
     & LMICROPHY, L_MURK, L_MURK_ADVECT, L_MURK_SOURCE,
     & L_MURK_BDRY, L_BL_TRACER_MIX, L_MOM, L_CAPE, L_SDXS, L_XSCOMP,
     & L_3D_CCA, L_CCW, L_PHASE_LIM, L_CLOUD_DEEP,
     & LSINGLE_HYDROL, LMOSES,
     & L_CO2_INTERACTIVE, L_CO2_EMITS,
     & L_SNOW_ALBEDO,
     & l_ssice_albedo,
     &  L_CLIMAT_AEROSOL,
     &  L_RHCPT,L_CLD_AREA,
     & L_LSPICE,L_BL_LSPICE,L_LSPICE_BDY,
     & L_SULPC_SO2, L_SULPC_DMS, L_SULPC_OZONE,
     & L_SO2_SURFEM, L_SO2_HILEM, L_SO2_NATEM, L_DMS_EM,
     & L_SULPC_NH3,L_NH3_EM,L_SOOT,L_SOOT_SUREM,L_SOOT_HILEM,
     & L_USE_SOOT_DIRECT,
     & L_IAU,L_IAU_RAMP,L_IAU_DIAG,
     & T_IAU_START, T_IAU_END,
     & L_USE_SULPC_DIRECT, L_USE_SULPC_INDIRECT_SW,
     & L_USE_SULPC_INDIRECT_LW, CALL_CHEM_FREQ,
     & A_ASSIM_MODE, H_SECT

      COMMON / CNTLCATM /
     & L_VEG_FRACS, L_TRIFFID, L_PHENOL, L_NRUN_MID_TRIF, L_TRIF_EQ,
     & PHENOL_PERIOD, TRIFFID_PERIOD,
     & H_SWBANDS, H_LWBANDS, A_ADJSTEPS, L_H2_SULPH,
     & A_SW_RADSTEP, A_LW_RADSTEP, A_SW_SEGMENTS, A_LW_SEGMENTS,
     & A_CONV_STEP, A_CONVECT_SEGMENTS, A_NSET_FILTER, A_ENERGYSTEPS,
     & A_ASSIM_START_HR, A_ASSIM_END_HR, A_SWEEPS_DYN,
     & L_SW_RADIATE, L_LW_RADIATE, LADD_RADINCS, L_SET_FILTER,
     & LDAY, LEXPAND_OZONE, L_COMPRESS_LAND,
     & L_NEG_THETA, L_NEG_PSTAR, L_NEG_QT, L_NEG_TSTAR, L_FIELD_FLT,
     & L_SUPERBEE, LENERGY, L_Z0_OROG, L_RMBL, L_MIXLEN, L_CONVECT,
     & L_HALF_TIMESTEP_DYN, L_HALF_TIMESTEP_DIV, L_QT_POS_LOCAL,
     & L_TRACER_THETAL_QT,LFLUX_RESET,
     & L_3DVAR_BG,L_AC,L_3DVAR,L_4DVAR,
     & L_OLD_PWTS, L_VINT_TP,
     & LFROUDE, LGWLINP, LLINTS, LWHITBROM, LEMCORR,
     & LMICROPHY, L_MURK, L_MURK_ADVECT, L_MURK_SOURCE,
     & L_MURK_BDRY, L_BL_TRACER_MIX, L_MOM, L_CAPE, L_SDXS, L_XSCOMP,
     & L_3D_CCA, L_CCW, L_PHASE_LIM, L_CLOUD_DEEP,
     & LSINGLE_HYDROL, LMOSES,
     & L_CO2_INTERACTIVE, L_CO2_EMITS,
     & L_SNOW_ALBEDO,
     & l_ssice_albedo,
     &  L_CLIMAT_AEROSOL,
     &  L_RHCPT,L_CLD_AREA,
     & L_LSPICE,L_BL_LSPICE,L_LSPICE_BDY,
     & L_SULPC_SO2, L_SULPC_DMS, L_SULPC_OZONE,
     & L_SO2_SURFEM, L_SO2_HILEM, L_SO2_NATEM, L_DMS_EM,
     & L_SULPC_NH3,L_NH3_EM,L_SOOT,L_SOOT_SUREM,L_SOOT_HILEM,
     & L_USE_SOOT_DIRECT,

     & L_IAU,L_IAU_RAMP,L_IAU_DIAG,
     & T_IAU_START, T_IAU_END,
     & L_USE_SULPC_DIRECT, L_USE_SULPC_INDIRECT_SW,
     & L_USE_SULPC_INDIRECT_LW, CALL_CHEM_FREQ,
     & A_ASSIM_MODE, H_SECT
! ----------------------- Comdeck: CNTLOCN  ----------------------------
! Description: COMDECK defining Control variables for the Ocean
!              internal model.
!   This comdeck contains logical variables which are used on the
!   control of certain sections of Ocean model code
!   They replace the previous method of controlling code using *IF DEFs.
!
! Author : R.T.H.Barnes & R.Hill
!
! History:
! Version  Date      Comment.
!  3.5  29/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.1  29/05/96  include L_OZVRT     M. J. Bell 
!  4.3  8.11.96   include L_SLOPEMAX and L_COXCNVC  JMG
!  4.4  11.08.97   Remove L_OCHEQUB.    R. Hill 
!    4.4  10/09/97  Remove all references to SKIPLAND code. R.Hill
!  4.4  8.07.97   include L_FLUXD      R.Lenton
!  4.5  3.11.98    include L_OBIMOM       M. Roberts
!  4.5  3.11.98   include L_OMEDADV and L_OHUDOUT  (new outflow param)
!                 M. Roberts
!  4.5  10.11.98  New logicals: L_OISOMOM, L_OISOGM L_OISOGMSKEW
!                 L_OBIHARMGM and L_OVISHADCM4
!  4.5  3.9.98    Changes for HADCM4 sea-ice. Cresswell and Gregory
!  4.5   1/07/98  Add logical to control interactive CO2 use. C.D.Jones
CLL   4.5 G.J.Rickard include L_OFULARGE (full Large scheme),
CLL                   L_OPANDP (choice of vertical mixing),
CLL                   L_OSTATEC (density calculation choice),
CLL                   L_OUSTARWME (ustar calculation).
CLL
!  4.5  7.8.97    Removed old ocean boundary logicals L_OBGILLS,
!                 L_OBGILLN, L_OSTEVNS, L_OSTEVS and L_BOUNDSO.  
!                 Added in new logicals L_OBDY_NORTH to L_OBDY_STREAM.
!                 M.J. Bell
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      INTEGER
     &        O_CLM_START_HR,     ! Time ocean climate increments start
     &        O_CLM_END_HR,       ! Time ocean climate increments end
     &        O_INT_CLM_INC,      ! # ocean steps  } climate incs.
     &        O_INT_ANA_STP,      ! # between      } analysis steps
!
     &        O_INT_EVO_BTS,      ! # ocean steps between fwd evolution
!                                     of bathys and tesacs
     &        O_INT_VRY_BTS,      ! # ocean steps between re-calculation
!                         of future bathys and tesacs valid at this hour
     &        O_INT_WTS_ACC,      ! # ocean steps betwn accumulating wts
!
     &        O_INT_OBS_FRSH,     ! # ocean  } reading new OBS files
     &        O_INT_OBS_OUT,      ! # steps  } outputting new OBS files
     &        O_INT_OBS_STR,      ! # between} caching OBS array
     &        O_INT_FLD_STR,      ! #        } caching model fields
     &        O_ASSIM_START_HR,   ! Time at which data assimilation
!                                 ! starts (Hours after Basis Time)
     &        O_ASSIM_END_HR      ! Time at which data assimilation
!                                 ! ends (Hours after Basis Time)
!
      LOGICAL L_FLUXCORR   ! Heat & water flux correction
     &,       L_OGLOBAL    ! Global ocean
     &,       L_ICEFREEDR  ! Free Drift Sea Ice model
     &,       L_ICESIMPLE  ! Simple Advection Sea Ice model
     &,       L_HADCM4O2I  ! HADCM4 version of ocean-to-ice heat flux
     &,       L_IHANEY     ! Haney Forcing Ice
     &,       L_OADGHR2    ! Ocean assimilation diagnostics
     &,       L_OBDY_NORTH   ! Update northern lateral boundary   
     &,       L_OBDY_SOUTH   ! Update southern lateral boundary
     &,       L_OBDY_EAST    ! Update eastern lateral boundary
     &,       L_OBDY_WEST    ! Update western lateral boundary
     &,       L_OGILL_LBCS   ! Use the Gill boundary scheme
     &,       L_OFRS_LBCS    ! Use the FRS boundary scheme
     &,       L_OSTVNS_LBCS  ! Use the Stevens boundary scheme
     &,       L_OBDY_TRACER  ! Update the tracers
     &,       L_OBDY_UV      ! Update the velocities
     &,       L_OBDY_STREAM  ! Update the stream functions
     &,       L_OBDY_ICE     ! Update ice fields (snow, aice, hice)
     &,       L_OBIOLOGY   ! Effect of phytoplankton on carbon cycle
     &,       L_OCARB14    ! Calculate atmospheric C12/C14 ratio
     &,       L_OCARBON    ! Carbon cycle model
     &,       L_OCNASSM    ! Activate ocean assimilation
     &,       L_OCYCLIC    ! Cyclic boundary conditions
     &,       L_OFILTER    ! Fourier filtering for high latitudes
     &,       L_OFREESFC   ! Use free surface conditions
     &,       L_FLUXD
     &,       L_OHANEY     ! Haney Forcing heat/fresh water fluxes
     &,       L_OHMEAD     ! Mead tracer transport diagnostics
     &,       L_OICECOUP   ! Coupled model with Sea Ice
     &,       L_OIMPADDF   ! Crank-Nicholson vert. advn-diffn scheme
     &,       L_OIMPDIF    ! CN vertical diffusion scheme
     &,       L_OISLANDS   ! Include Island Routines
     &,       L_OISOPYC    ! Isopycnal diffusion scheme
     &,       L_OLATVISC   ! Latitude dependent viscosity
     &,       L_OLISS      ! Liss & Merlivat wind mixing of tracers
     &,       L_OMIXLAY    ! Wind mixing of tracers-mixed layer scheme
     &,       L_ONOCLIN    ! Barotropic solution
     &,       L_ONOPOLO    ! No sea ice at North Pole
     &,       L_OPENBC     ! Read in lateral boundary fields
     &,       L_OPSEUDIC   ! Pseudo-ice routine
     &,       L_ORICHARD   ! Evaluate & use Richardson No.
     &,       L_OROTATE    ! Coriolis force calculation
     &,       L_OSOLAR     ! Calc solar penetration for given water type
     &,       L_OSOLARAL   ! Calc sol. pen. - simplified layer structure
     &,       L_OSYMM      ! Symmetric boundary conditions
     &,       L_OVARYT     ! Varying time step with depth
     &,       L_RIVERS     ! River run-off routines
     &,       L_SEAICE     ! Include Sea Ice model
     &,       L_TRANGRID   ! Spatial interp. in coupled model
     &,       L_OCONJ     ! Whether to use conjugate gradient solver
     &,       L_UPWIND     ! Upwind differencing for tracer advection
     &,       L_OPRINT     ! Whether to print incidental ocean info
     &,       L_OSTVEW     !\
     &,       L_OPMSL      ! \
     &,       L_OTIDAL     !  \___ All for use with O. Alves free
     &,       L_OFOURW     !  /    surface modifications at V4.0
     &,       L_ODELPLUS   ! /
     &,       L_OTROPIC    !/
     &,       L_OISOMOM
     &,       L_OISOGMSKEW
     &,       L_OISOGM
     &,       L_OBIHARMGM
     &,       L_OVISHADCM4
     &,       L_OMEDOUT    ! Mediterranean outflow - 288*144 and 96*73
                           !   grids only - uses hardwired gridpoint nos
     &,       L_OCONVROUS  ! Roussenov convective adjustment
     &,       L_OEXTRAP    ! Extrapolation of vertical density gradients
     &,       L_OISOPYCGM  ! Gent and McWilliams eddy parametrisation.
     &,       L_OISOTAPER  ! Tapering of isopycnal diffusion
     &,       L_OVISBECK    ! Visbeck scheme
     &,       L_OQLARGE     ! Quadratic Large scheme
     &,       L_OFULARGE   ! FULL LARGE SCHEME
     &,       L_OPANDP     ! RI-DEPENDENT VERT MIX SCHEMES
     &,       L_OSTATEC    ! DENSITY CHOICE FOR RI-CALC
     &,       L_OUSTARWME  ! WME OR WSTRESS TO FIND USTAR

     &,       L_OZVRT      ! barotropic vorticity diagnostic switch
                           ! set by OCN_FOR_STEP (not in namelist) 
     &,       L_SLOPEMAX   ! Selects SLOPE_MAX isopycnal diffusion
     &,       L_COXCNVC    ! Selects original Cox convection scheme
     &,       L_OCOMP    ! Land pnts compressed from dump (3d fields)
     &,       L_OMEDADV
     &,       L_OHUDOUT
     &,L_REFSAL
     &,L_SALFLUXFIX
     &,L_INLANSEA
     &      ,L_CO2O_INTERACTIVE     ! interactive 3D CO2 field for
                                    !  use with carbon cycle model
     &,       L_OBIMOM  ! biharmonic momentum diffusion
! *IF DEF,OCNASSM
!    Additions to CCONTROL for ocean assimilation
!
      LOGICAL
     &       LAS_CLM_INC,    ! make increments to relax to climate
     &       LAS_ADD_INC,    ! add or subtract analysis increments
     &       LAS_ANA_STP,    ! calculate analysis increments
     &       LAS_EVO_BTS,    ! evolve bathy and tesac obs 1 step
     &       LAS_VRY_BTS,    ! estimate bathys and tesacs at this hour
     &       LAS_WTS_ACC,    ! evolve accumulated weights
     &       LAS_OBS_FRSH,   ! to refresh main OBS data set
     &       LAS_OBS_OUT,    ! output ACOBS file for incremented obs
     &       LAS_FLD_STR,    ! output model fields to cache store
     &       LAS_OBS_STR     ! output obs to cache store
! *ENDIF OCNASSM


      NAMELIST / NLSTCOCN /
     & O_CLM_START_HR, O_CLM_END_HR, O_INT_CLM_INC, O_INT_ANA_STP,
     & O_INT_EVO_BTS, O_INT_VRY_BTS, O_INT_WTS_ACC, O_INT_OBS_FRSH,
     & O_INT_OBS_OUT, O_INT_OBS_STR, O_INT_FLD_STR,
     & O_ASSIM_START_HR, O_ASSIM_END_HR, L_FLUXCORR, L_OGLOBAL,
     & L_ICEFREEDR, L_ICESIMPLE, L_IHANEY, L_HADCM4O2I, L_OADGHR2,
     & L_OBDY_NORTH,L_OBDY_SOUTH,L_OBDY_EAST,L_OBDY_WEST,
     & L_OGILL_LBCS,L_OFRS_LBCS,L_OSTVNS_LBCS,
     & L_OBDY_TRACER,L_OBDY_UV,L_OBDY_STREAM,L_OBDY_ICE,
     & L_OBIOLOGY, L_OCARB14, L_OCARBON, L_OCNASSM,
     & L_OCYCLIC, L_OFILTER, L_OFREESFC, L_FLUXD,
     & L_OHANEY, L_OHMEAD, L_OICECOUP,
     & L_OIMPADDF, L_OIMPDIF, L_OISLANDS, L_OISOPYC, L_OLATVISC,
     & L_OLISS, L_OMIXLAY, L_ONOCLIN, L_ONOPOLO, L_OPENBC, L_OPSEUDIC,
     & L_ORICHARD, L_OROTATE, L_OSOLAR, L_OSOLARAL,        
     & L_OSYMM, L_OVARYT, L_RIVERS, L_SEAICE, L_OCONJ,
     & L_TRANGRID, L_UPWIND, L_OSTVEW, L_OPMSL, L_OTIDAL, L_OPRINT,  
     & L_OFOURW, L_ODELPLUS, L_OTROPIC
     &, L_OISOMOM,L_OISOGMSKEW,L_OISOGM,L_OBIHARMGM,L_OVISHADCM4
     &, L_OMEDOUT
     &, L_OCONVROUS
     &, L_OEXTRAP,L_OISOPYCGM,L_OISOTAPER
     &, L_OVISBECK
     &,L_OBIMOM
     &, L_OQLARGE
     &,L_OCOMP
     &, L_OMEDADV,L_OHUDOUT
     &,L_REFSAL
     &,L_SALFLUXFIX
     &,L_INLANSEA 
     &, L_CO2O_INTERACTIVE
     &, L_OFULARGE,L_OPANDP,L_OSTATEC,L_OUSTARWME
! *IF DEF,OCNASSM
     &,L_SLOPEMAX,L_COXCNVC
!        additions for control of ocean assimilation
     &              ,LAS_ADD_INC,LAS_CLM_INC,LAS_ANA_STP
     &              ,LAS_EVO_BTS,LAS_VRY_BTS,LAS_WTS_ACC
     &              ,LAS_OBS_FRSH,LAS_OBS_OUT,LAS_FLD_STR,LAS_OBS_STR
! *ENDIF OCNASSM

      COMMON / CNTLCOCN /

     & O_CLM_START_HR, O_CLM_END_HR, O_INT_CLM_INC, O_INT_ANA_STP,
     & O_INT_EVO_BTS, O_INT_VRY_BTS, O_INT_WTS_ACC, O_INT_OBS_FRSH,
     & O_INT_OBS_OUT, O_INT_OBS_STR, O_INT_FLD_STR,
     & O_ASSIM_START_HR, O_ASSIM_END_HR, L_FLUXCORR, L_OGLOBAL,
     & L_ICEFREEDR, L_ICESIMPLE, L_IHANEY, L_HADCM4O2I, L_OADGHR2,
     & L_OBDY_NORTH,L_OBDY_SOUTH,L_OBDY_EAST,L_OBDY_WEST,
     & L_OGILL_LBCS,L_OFRS_LBCS,L_OSTVNS_LBCS,
     & L_OBDY_TRACER,L_OBDY_UV,L_OBDY_STREAM,L_OBDY_ICE,
     & L_OBIOLOGY, L_OCARB14, L_OCARBON, L_OCNASSM, 
     & L_OCYCLIC, L_OFILTER, L_OFREESFC, L_FLUXD,
     & L_OHANEY, L_OHMEAD, L_OICECOUP,
     & L_OIMPADDF, L_OIMPDIF, L_OISLANDS, L_OISOPYC, L_OLATVISC,
     & L_OLISS, L_OMIXLAY, L_ONOCLIN, L_ONOPOLO, L_OPENBC, L_OPSEUDIC,
     & L_ORICHARD, L_OROTATE, L_OSOLAR, L_OSOLARAL,               
     & L_OSYMM, L_OVARYT, L_RIVERS, L_SEAICE, L_OCONJ, 
     & L_TRANGRID, L_UPWIND, L_OSTVEW, L_OPMSL, L_OTIDAL, L_OPRINT,
     & L_OFOURW, L_ODELPLUS, L_OTROPIC, L_OZVRT 
     &, L_OMEDOUT,L_OISOMOM,L_OISOGMSKEW,L_OISOGM
     &, L_OBIHARMGM,L_OVISHADCM4
     &, L_OCONVROUS
     &, L_OEXTRAP,L_OISOPYCGM,L_OISOTAPER
     &, L_OVISBECK
     &,L_OBIMOM

     &, L_OQLARGE
     &,L_OCOMP
     &, L_OMEDADV,L_OHUDOUT

     &,L_REFSAL
     &,L_SALFLUXFIX
     &,L_INLANSEA 
     &, L_CO2O_INTERACTIVE
     &, L_OFULARGE,L_OPANDP,L_OSTATEC,L_OUSTARWME
! *IF DEF,OCNASSM
     &,L_SLOPEMAX,L_COXCNVC
!        additions for control of ocean assimilation
     &              ,LAS_ADD_INC,LAS_CLM_INC,LAS_ANA_STP
     &              ,LAS_EVO_BTS,LAS_VRY_BTS,LAS_WTS_ACC
     &              ,LAS_OBS_FRSH,LAS_OBS_OUT,LAS_FLD_STR,LAS_OBS_STR
! *ENDIF OCNASSM

! ----------------------- Comdeck: CNTLSLB  ----------------------------
! Description: COMDECK defining Control variables for the Slab
!              internal model.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  29/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      LOGICAL                       ! Logical switches for:
     &   L_THERM        ,  !         :Coupled model ice thermodynamics
     &   L_IDYN         ,  !Slab     :Cavitating fluid ice dynamics
     &   L_IDRIF        ,  !model    :Simple ice depth advection
     &   L_SLBADV          !

      NAMELIST / NLSTCSLB /
     & L_THERM, L_IDYN, L_IDRIF, L_SLBADV

      COMMON / CNTLCSLB /
     & L_THERM, L_IDYN, L_IDRIF, L_SLBADV
! ----------------------- Comdeck: CNTLWAV  ----------------------------
! Description: COMDECK defining Control variables for the Wave
!              internal model, and its runtime constants (if any).
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  4.1  23/02/96  New comdeck for wave sub-model.  RTHBarnes.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Parameter declarations
!??   INTEGER MAXSECTS            ! Max. no. of code sections
!??   PARAMETER (MAXSECTS=99)
!
!   Type declarations
!
      INTEGER
     & W_N_SRCE,    ! no.of source timesteps per propagation timestep
     & W_ISHALLOW,  ! 1 for deep, otherwise shallow water
     & W_IREFRACT,  ! refraction options, 0 = none,
!                   ! 1 = depth, 2 = depth & current
     & W_ICASE,     ! 1 for spherical propagation, otherwise Cartesian
     & W_IPER       ! 1 for , otherwise

      LOGICAL
     & L_WAVASSM    ! True if assimilation requested
!
      CHARACTER*5 W_ASSIM_MODE     ! Switch for BG/AC/3DVAR/4DVAR assm
!
!     CHARACTER*3 H_SECT(0:MAXSECTS) ! Array of code section versions
!
      NAMELIST / NLSTCWAV /
     & W_N_SRCE, W_ISHALLOW, W_IREFRACT, W_ICASE, W_IPER, L_WAVASSM,
     & W_ASSIM_MODE
      COMMON / CNTLCWAV /
     & W_N_SRCE, W_ISHALLOW, W_IREFRACT, W_ICASE, W_IPER, L_WAVASSM,
     & W_ASSIM_MODE

C
C Local variables
C
      CHARACTER*132 line,line1,line2 ! Encoded line of information
      CHARACTER*80  ch            ! Working character string variable
      INTEGER i1,i2,k             ! Array indices
      INTEGER j                   ! Code value
      INTEGER time_list,lev_list  ! pointers to time and levels lists
     &,       plev_list           ! pointer  to pseudo-level list
     &,       tser_list           ! pointer  to time series record list
      INTEGER ntimes              ! no of times in a time list
      INTEGER packing_profile     ! packing profile for output PPfield
C
CL----------------------------------------------------------------------
CL 0. Write header if sequence no indicates first item
CL
      IF (seqno.EQ.1) THEN
        WRITE(6,1000)
 1000   FORMAT(
     *  '   ********************************************************'/
     *  '   ********************************************************'/
     *  '   **                                                    **'/
     *  '   **    LIST OF USER-DEFINED DIAGNOSTICS IN THIS RUN    **'/
     *  '   **                                                    **'/
     *  '   ********************************************************'/
     *  '   ********************************************************'/
     *  '   **                                                    **'/
     *  '   ** NOTES:                                             **'/
     *  '   **   Time processing details are in timesteps, where  **'/
     *  '   **     ... represents "for ever".                     **'/
     *  '   **   Spatial processing domain is in gridpoints.      **'/
     *  '   **                                                    **'/
     *  '   ********************************************************'/
     *  '   ********************************************************'//
     *'=================================================================
     *==========================================================')
      ENDIF
CL----------------------------------------------------------------------
CL 1. For each diagnostic processing request in the STASHlist,
CL    print the diagnostic name followed by a summary of the processing
CL    information on 3 lines.
CL
CL 1.0 If diagnostic is not required for output, exit routine
CL
      IF (stlist(st_proc_no_code).EQ.0) GOTO 999
CL
CL 1.1 Line 1.
CL
      line=' '
C #No
      i1=2
      i2=4
      write(ch,'(i3)') seqno
      line(i1:i2)=ch(1:1+i2-i1)
C Name
      i1=i2+2
      i2=i1+36-1
      line(i1:i2)=name
C Submodel
      i1=i2+2
      i2=i1+8-1
      j=stlist(st_sect_no_code)
      IF      (ppxref(ppx_model_number).EQ.ocean_im) THEN
        ch=' OCEAN  '
      ELSE IF (ppxref(ppx_model_number).EQ. slab_im) THEN
        ch=' SLAB   '
      ELSE IF (ppxref(ppx_model_number).EQ.atmos_im) THEN
        ch=' ATMOS  '
      ELSE IF (ppxref(ppx_model_number).EQ.wave_im) THEN      
        ch=' WAVE   '                              
      ELSE
        WRITE(6,*)' Error in DIAGDES. Unknown model'
        ch=' UNKNOWN'
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
C Item
      i1=i2+2
      i2=i1+4-1
      j=stlist(st_item_code)
      write(ch,'(i4)') j
      line(i1:i2)=ch(1:1+i2-i1)
C Section
      i1=i2+2
      i2=i1+7-1
      j=stlist(st_sect_no_code)
      write(ch,'(i7)') j
      line(i1:i2)=ch(1:1+i2-i1)
C PPfcode
      i1=i2+2
      i2=i1+7-1
      j=ppxref(ppx_field_code)
      write(ch,'(i7)') j
      line(i1:i2)=ch(1:1+i2-i1)
C Datatype
      i1=i2+2
      i2=i1+8-1
      j=ppxref(ppx_data_type)
      IF (j.EQ.1.OR.j.EQ.4) THEN
        ch='  REAL  '
      ELSEIF (j.EQ.2.OR.j.EQ.5) THEN
        ch='INTEGER '
      ELSEIF (j.EQ.3) THEN
        ch='LOGICAL '
      ELSE
        ch='UNKNOWN '
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
C Gridtype
      i1=i2+2
      i2=i1+8-1
      j=ppxref(ppx_grid_type)
      IF (j.EQ.ppx_atm_nonstd.OR.j.EQ.ppx_ocn_nonstd) THEN
        ch=' NONSTD '
      ELSEIF ((j.GT.ppx_atm_nonstd.AND.j.LE.ppx_atm_tsea) .OR.
     &         j.EQ.ppx_atm_compressed.OR.j.EQ.ppx_atm_ozone) THEN
        ch=' P-GRID '
      ELSEIF (j.GE.ppx_atm_uall.AND.j.LE.ppx_atm_usea) THEN
        ch=' UV-GRID'
      ELSEIF (j.EQ.ppx_atm_cuall.OR.j.EQ.ppx_ocn_cuall) THEN
        ch=' CU-GRID'
      ELSEIF (j.EQ.ppx_atm_cvall.OR.j.EQ.ppx_ocn_cvall) THEN
        ch=' CV-GRID'
      ELSEIF (j.EQ.ppx_atm_tzonal) THEN
        ch=' PZ-GRID'
      ELSEIF (j.EQ.ppx_atm_uzonal) THEN
        ch=' UZ-GRID'
      ELSEIF (j.EQ.ppx_atm_tmerid) THEN
        ch=' PM-GRID'
      ELSEIF (j.EQ.ppx_atm_umerid) THEN
        ch=' UM-GRID'
      ELSEIF (j.EQ.ppx_atm_rim.OR.j.EQ.ppx_ocn_rim) THEN
        ch='   RIM  '
      ELSEIF (j.EQ.ppx_ocn_tcomp.OR.j.EQ.ppx_ocn_tall.OR.
     &        j.EQ.ppx_ocn_tfield) THEN
        ch=' T-GRID '
      ELSEIF (j.EQ.ppx_ocn_tzonal) THEN
        ch=' TZ-GRID'
      ELSEIF (j.EQ.ppx_ocn_uzonal) THEN
        ch=' UZ-GRID'
      ELSEIF (j.EQ.ppx_ocn_tmerid) THEN
        ch=' TM-GRID'
      ELSEIF (j.EQ.ppx_ocn_umerid) THEN
        ch=' UM-GRID'
      ELSEIF (j.EQ.ppx_ocn_ucomp.OR.j.EQ.ppx_ocn_uall.OR.
     &        j.EQ.ppx_ocn_ufield) THEN
        ch=' UV-GRID'
      ELSEIF (j.EQ.ppx_atm_scalar.OR.j.EQ.ppx_ocn_scalar) THEN
        ch=' SCALAR '
      ELSEIF (j.EQ.ppx_wam_all.OR.j.EQ.ppx_wam_sea) THEN        
        ch=' WAVE   ' 
      ELSEIF (j.EQ.ppx_wam_rim) THEN        
        ch=' RIM    ' 
      ELSE
        ch=' UNKNOWN'
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
C Leveltype
      i1=i2+2
      i2=i1+9-1
      j=ppxref(ppx_lv_code)
      IF (j.EQ.ppx_full_level) THEN
        ch='FULLLEVEL'
      ELSEIF (j.EQ.ppx_half_level) THEN
        ch='HALFLEVEL'
      ELSE
        ch='STD-LEVEL'
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
C Meto8LV
      i1=i2+2
      i2=i1+7-1
      j=ppxref(ppx_meto8_levelcode)
      write(ch,'(i7)') j
      line(i1:i2)=ch(1:1+i2-i1)
C Meto8FC
      i1=i2+2
      i2=i1+7-1
      j=ppxref(ppx_meto8_fieldcode)
      write(ch,'(i7)') j
      line(i1:i2)=ch(1:1+i2-i1)
C PackAcc
      i1=i2+2
      i2=i1+7-1
      j=stlist(st_output_code)
      IF (j.EQ.1) THEN
        IF (stlist(st_macrotag).GE.1000) THEN
          packing_profile=pp_pack_code(27)
        ELSE
          packing_profile=0
        ENDIF
      ELSEIF(j.eq.2) THEN
        packing_profile=0
      ELSEIF(j.lt.0) THEN
        packing_profile=pp_pack_code(-j)
      ELSE
        packing_profile=0
      ENDIF
      IF (packing_profile.EQ.0) THEN
        ch='       '
      ELSE
        j=ppxref(ppx_packing_acc+packing_profile-1)
        write(ch,'(i7)') j
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
C
      line1=line
CL
CL 1.2 Line 2.
CL
      line=' '
C Time-processing
      i1=2
      i2=16
      j=stlist(st_proc_no_code)
      tser_list=0
      IF (j.EQ.st_replace_code) THEN
        ch='   EXTRACT     '
      ELSEIF (j.EQ.st_accum_code) THEN
        ch=' ACCUMULATION  '
      ELSEIF (j.EQ.st_time_mean_code) THEN
        ch='  TIME MEAN    '
      ELSEIF (j.EQ.st_time_series_code) THEN
        write(ch,'(''  TIME SERIES  '')')
        tser_list=stlist(st_series_ptr)
      ELSEIF (j.EQ.st_max_code) THEN
        ch='MAX OVER PERIOD'
      ELSEIF (j.EQ.st_min_code) THEN
        ch='MIN OVER PERIOD'
      ELSEIF (j.EQ.st_append_traj_code) THEN
        ch='  TRAJECTORY   '
      ELSEIF (j.EQ.st_variance_code) THEN
        ch=' TIME VARIANCE '
      ELSEIF (j.EQ.st_time_series_mean) THEN  
        ch='MEAN TIMESERIES'                  
      ELSE
        ch='  UNKNOWN      '
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
C -From-
      i1=i2+2
      i2=i1+6-1
      IF (stlist(st_freq_code).LT.0) THEN
        ch='      '
      ELSE
        j=stlist(st_start_time_code)
        write(ch,'(i6)') j
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
C --To--
      i1=i2+2
      i2=i1+6-1
      IF (stlist(st_freq_code).LT.0) THEN
        ch='      '
      ELSE
        j=stlist(st_end_time_code)
        IF (j.EQ.st_infinite_time) THEN
          ch='  ... '
        ELSE
          write(ch,'(i6)') j
        ENDIF
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
C Frequency
      i1=i2+2
      i2=i1+9-1
      j=stlist(st_freq_code)
      IF (j.LT.0) THEN
        j=-j
        write(ch,'(''TIME LIST'')')
        time_list=j
      ELSE
        write(ch,'(i9)') j
        time_list=0
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
C Period
      i1=i2+2
      i2=i1+6-1
      IF (stlist(st_freq_code).LT.0) THEN
        ch='      '
      ELSE
        j=stlist(st_period_code)
        IF (stlist(st_proc_no_code).EQ.st_replace_code) THEN
          ch='      '
        ELSEIF (j.EQ.st_infinite_time) THEN
          ch='  ... '
        ELSE
          write(ch,'(i6)') j
        ENDIF
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
C __Source__
      i1=i2+2
      i2=i1+10-1
      j=stlist(st_input_code)
      IF (j.EQ.0) THEN
        ch='PROGNOSTIC'
      ELSEIF(j.EQ.1) THEN
        ch='  STWORK  '
      ELSEIF(j.LT.0) THEN
        j=-j
        write(ch,'(''DUMP #'',i4)') j
      ELSE
        ch=' UNKNOWN  '
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
C ___Destination___
      i1=i2+2
      i2=i1+17-1
      j=stlist(st_output_code)
      IF (j.EQ.1) THEN
        IF (stlist(st_macrotag).GE.1000) THEN
          ch='MEAN PP VIA DUMP'
        ELSEIF (stlist(st_macrotag).GT.0) THEN
          write(ch,'(''DUMP WITH TAG '',i3)') stlist(st_macrotag)
        ELSE
          ch='      DUMP       '
        ENDIF
      ELSEIF(j.eq.2) THEN
        ch='   SECONDARY     '
      ELSEIF(j.lt.0) THEN
        j=-j
        IF (j.EQ.27) THEN
          ch='MEAN PP (DIRECT) '
        ELSE
          write(ch,'(''   PP UNIT '',i2)') j
        ENDIF
      ELSE
        ch='  UNKNOWN  '
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
C
      line2=line
CL
CL 1.3 Line 3.
CL
      line=' '
C Spatial-Processing
      i1=2
      i2=19
      j=stlist(st_gridpoint_code)
      IF (j.GE.extract_base.AND.j.LT.extract_top) THEN
        ch='    FULL FIELD    '
      ELSEIF (j.GE.vert_mean_base.AND.j.LT.vert_mean_top) THEN
        ch='  VERTICAL MEAN   '
      ELSEIF (j.GE.zonal_mean_base.AND.j.LT.zonal_mean_top) THEN
        ch='   ZONAL MEAN     '
      ELSEIF (j.GE.merid_mean_base.AND.j.LT.merid_mean_top) THEN
        ch=' MERIDIONAL MEAN  '
      ELSEIF (j.GE.field_mean_base.AND.j.LT.field_mean_top) THEN
        ch=' FIELD MEAN - 2D  '
      ELSEIF (j.GE.global_mean_base.AND.j.LT.global_mean_top) THEN
        ch=' GLOBAL MEAN - 3D '
      ELSE
        ch='  ** UNKNOWN **   '
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
C Levels-domain
      i1=i2+2
      i2=i1+13-1
      j=stlist(st_output_bottom)
      lev_list=0
      IF (j.EQ.st_special_code) THEN
        ch='STANDARD LEV '
      ELSEIF (j.gt.0) THEN
        write(ch,'(''LEVELS '',i2,''-'',i2)') j,stlist(st_output_top)
      ELSEIF (j.lt.0) THEN
        j=-j
        write(ch,'('' LEVELS LIST '')')
        lev_list=j
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
C Pseudo-levels
      i1=i2+2
      i2=i1+15-1
      j=stlist(st_pseudo_out)
      plev_list=0
      IF (j.GT.0) THEN
        ch='PSEUDO-LEV LIST'
        plev_list=j
      ELSE
        ch='     NONE      '
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
C Horizontal-domain.....
      i1=i2+2
      i2=i1+23-1
      write(ch,'(''ROW:'',i3,''-'',i3,'' COL:'',i3,''-'',i3)')
     *  stlist(st_north_code),stlist(st_south_code),
     *  stlist(st_west_code),stlist(st_east_code)
      line(i1:i2)=ch(1:1+i2-i1)
C Weighting
      i1=i2+2
      i2=i1+9-1
      j=stlist(st_weight_code)
      IF (j.EQ.stash_weight_null_code) THEN
        ch='  NONE   '
      ELSEIF (j.EQ.stash_weight_area_code) THEN
        ch='  AREA   '
      ELSEIF (j.EQ.stash_weight_volume_code) THEN
        ch=' VOLUME  '
      ELSEIF (j.EQ.stash_weight_mass_code) THEN
        ch='  MASS   '
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
C Masking
      i1=i2+2
      i2=i1+7-1
      j=mod(stlist(st_gridpoint_code),block_size)
      IF (j.EQ.stash_null_mask_code) THEN
        ch=' NONE  '
      ELSEIF (j.EQ.stash_land_mask_code) THEN
        ch=' LAND  '
      ELSEIF (j.EQ.stash_sea_mask_code) THEN
        ch='  SEA  '
      ELSE
        ch='UNKNOWN'
      ENDIF
      line(i1:i2)=ch(1:1+i2-i1)
CL
CL 1.4 Print the main part of the summary
CL
      WRITE(6,1010) line1,line2,line
 1010 FORMAT(' #No ',
     *      'Diagnostic Description-------------- Submodel Item Section
     *PPfcode Datatype Gridtype Leveltype MetO8lv MetO8fc Packacc'/
     *  a124/
     *' Time-processing -From- --To-- Frequency Period --Source-- ---Des
     *tination---                                               '/
     *  a124/
     *' Spatial-processing Levels-domain -Pseudo-levels- ---Horizontal-d
     *omain--- Weighting Masking                                '/
     *  a124)
CL
CL 1.5 Print associated time and levels lists if appropriate
CL
CL 1.5.1 Time list
CL
      IF (time_list.NE.0) THEN
        DO j=1,nsttims
          IF (sttabl(j,time_list).EQ.st_end_of_list) THEN
            ntimes=j-1
            GOTO 210
          ENDIF
        ENDDO
  210   CONTINUE
        WRITE(6,'('' ***** TIME LIST ***** '',i3,
     &            '' times are as follows:-'')') ntimes
        i1=1
        i2=8
        DO j=1,ntimes
          IF (i1.EQ.1) line=' '
          WRITE(ch,'(1x,i7)') sttabl(j,time_list)
          line(i1:i2)=ch(1:8)
          i1=i1+8
          i2=i2+8
          IF (i2.GT.80) THEN
            i1=1
            i2=8
            WRITE(6,'(a80)') line
          ENDIF
        ENDDO
        IF (i2.LE.80) WRITE(6,'(a80)') line
      ENDIF
CL
CL 1.5.2 Levels list
CL
      IF (lev_list.NE.0) THEN
        write(6,'('' ***** LEVELS LIST ***** '',i3,
     &          '' levels are as follows:-'')') stash_levels(1,lev_list)
        i1=1
        i2=8
        DO j=2,1+stash_levels(1,lev_list)
          IF (i1.EQ.1) line=' '
          write(ch,'(1x,i7)') stash_levels(j,lev_list)
          line(i1:i2)=ch(1:8)
          i1=i1+8
          i2=i2+8
          IF (i2.GT.80) THEN
            i1=1
            i2=8
            write(6,'(a80)') line
          ENDIF
        ENDDO
        IF (i2.LE.80) write(6,'(a80)') line
      ENDIF
CL
CL 1.5.3 Pseudo-levels list
CL
      IF (plev_list.NE.0) THEN
        write(6,'('' ***** PSEUDO-LEVELS LIST ***** '',i3,
     &          '' pseudo-levels are as follows:-'')')
     &          stash_pseudo_levels(1,plev_list)
        i1=1
        i2=8
        DO j=2,1+stash_pseudo_levels(1,plev_list)
          IF (i1.EQ.1) line=' '
          write(ch,'(1x,i7)') stash_pseudo_levels(j,plev_list)
          line(i1:i2)=ch(1:8)
          i1=i1+8
          i2=i2+8
          IF (i2.GT.80) THEN
            i1=1
            i2=8
            write(6,'(a80)') line
          ENDIF
        ENDDO
        IF (i2.LE.80) write(6,'(a80)') line
      ENDIF
CL
CL 1.5.4 Time series subdomain record list
CL
      IF (tser_list.NE.0) THEN
        i1=stash_series_index(1,tser_list)
        i2=stash_series_index(2,tser_list)
        WRITE(6,'('' ***** TIME SERIES ***** '',i3,
     & '' subdomain records are as follows:-''/
     & '' Record      North/South       West/ East     Bottom/  Top'')')
     &    i2
        DO j=1,i2
          WRITE(6,'(3x,i4,1x,3(5x,i5,1x,i5,1x))')
     &        j,(stash_series(3+k,i1+j-1),k=1,6)
        ENDDO
      ENDIF
CL
CL 1.5.5 Print final ruler line
CL
      WRITE(6,1020)
 1020 FORMAT(
     *'=================================================================
     *==========================================================')
C
 999  CONTINUE
      RETURN
      END
