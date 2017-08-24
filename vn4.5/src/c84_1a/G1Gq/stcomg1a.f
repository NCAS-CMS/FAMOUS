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
CLL  SUBROUTINE STASH_COMP_GRID----------------------
CLL
CLL  Compute grid descriptors
CLL  out_bdx,out_bdy,out_bzx and out_bzy from inputs
CLL  1)  samples
CLL  2)  grid type code
CLL  3)  ocean (true if an ocean grid)
CLL  4)  real and integer headers
CLL  5)  w_mostcol and n_mostrow
CLL  6) processing code
CLL
CLL  Tested under compiler CFT77
CLL  Tested under OS version 6.1
CLL
CLL  Author Simon Tett (Based on  code in pp_head by Tim Johns)
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL  Logical components covered: D40
CLL
CLL  Project TASK: C4
CLL
CLL  Programming standard: U M DOC  Paper NO. 4,
CLL
CLL  External documentation  C4
CLL
CLLEND-------------------------------------------------------------

C
C*L  INTERFACE and ARGUMENTS:------------------------------------------
      SUBROUTINE STASH_COMP_GRID(
     1  out_bzx,out_bzy,out_bdx,out_bdy,
     2  samples,st_grid,ocean,
     3  w_mostcol,n_mostrow,
     4  realhd,len_realhd,inthd,len_inthd,gr,
     5  ICODE,CMESSAGE)
C
      IMPLICIT NONE


C
      LOGICAL  OCEAN          !IN TRUE if processing an ocean diagnostic
C
      INTEGER samples                ! IN no of samples in period (times
C
      INTEGER
     *  ICODE             !OUT   Return code from the routine
     *, LEN_REALHD        !IN    Length of the Real Constants
     *, LEN_INTHD         !IN    Length of integer constants
     *, INTHD(LEN_INTHD)  !IN    Integer constants
C
      INTEGER
     *  st_grid           !IN    STASH horizontal grid type
     *, N_MOSTROW         !IN    The most nrthrly row.
     *, W_MOSTCOL         !IN    The most westerly column
     *, gr                !IN    The type of processing done
C
      REAL
     *  REALHD(LEN_REALHD)  !IN  Real header
      REAL
     *  out_bzx,out_bdx,out_bzy,out_bdy  ! OUT grid descriptors
      CHARACTER*(*)
     *  cmessage           ! OUT error messages
C*
C
CL Local Variables
C
      INTEGER mean_code
C
C*---------------------------------------------------------------------
C*L------------------ COMDECK LOOKADD ----------------------------------
CLL
CLL Purpose : Contains information about the format
CLL           of the PP header
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   4.0  12/09/95   Change NPERIODS to LBUSER3, BRSVD1 to BULEV,
CLL                   BRSVD2 to BHULEV and definitions for BRLEV and
CLL                   BHRLEV. Corresponding changes made to STWORK1A
CLL                   and PPHEAD1A. (Andrew Brady)   
CLL  4.0  12/10/95  Change item 45 from lbuser7 to model_code. RTHBarnes
CLL
CLL Programming standard :
CLL
CLL Logical components covered : F092
CLL
CLL Project task :
CLL
CLL External documentation:
CLL
CLLEND -----------------------------------------------------------------
C
      INTEGER
C Validity time
     &       LBYR,       ! Year
     &       LBMON,      ! Month
     &       LBDAT,      ! Day of month
     &       LBHR,       ! Hour
     &       LBMIN,      ! Minute
     &       LBDAY       ! Day number

C Data time

      INTEGER
     &       LBYRD,      ! Year
     &       LBMOND,     ! Month
     &       LBDATD,     ! Day of month
     &       LBHRD,      ! Hour
     &       LBMIND,     ! Minute
     &       LBDAYD      ! Day number

      INTEGER
     &       LBTIM,      ! Time indicator
     &       LBFT,       ! Forcast period (hours)
     &       LBLREC,     ! Length of data record
     &       LBCODE,     ! Grid type code
     &       LBHEM,      ! Hemisphere indicator
     &       LBROW,      ! Number of rows in grid
     &       LBNPT,      ! Number of points per row
     &       LBEXT,      ! Length of extra data
     &       LBPACK,     ! Packing method indicator
     &       LBREL       ! Header release number

      INTEGER
     &       LBFC,       ! Field code
     &       LBCFC,      ! Second field code
     &       LBPROC,     ! Processing code
     &       LBVC,       ! Vertical coordinate type
     &       LBRVC,      ! Coordinate type for reference level
     &       LBEXP,      ! Experiment number
     &       LBEGIN,     ! Start record
     &       LBNREC,     ! No of records-Direct access only
     &       LBPROJ,     ! Met-O-8 projection number
     &       LBTYP,      ! Met-O-8 field type
     &       LBLEV,      ! Met-O-8 level code
     &       LBRSVD1,    ! Reserved for future PP-package use
     &       LBRSVD2,    ! Reserved for future PP-package use
     &       LBRSVD3,    ! Reserved for future PP-package use
     &       LBRSVD4,    ! Reserved for future PP-package use
     &       LBSRCE      ! =1111 to indicate following apply to UM
      INTEGER
     &       DATA_TYPE,  ! Indicator for real/int or timeseries
     &       NADDR,      ! Start address in DATA_REAL or DATA_INT
     &       LBUSER3,    ! Free for user-defined function   
     &       ITEM_CODE,  ! Stash code
     &       LBPLEV,     ! Pseudo-level indicator (if defined)
     &       LBUSER6,    ! Free for user-defined function
     &       MODEL_CODE ! internal model identifier
      INTEGER
     &       BULEV,      ! Upper level boundary (Bk for ATMOS)
     &       BHULEV,     ! Upper level boundary (Ak for ATMOS)   
     &       BRSVD3,     ! Reserved for future PP-package use
     &       BRSVD4,     ! Reserved for future PP-package use
     &       BDATUM,     ! Datum value
     &       BACC,       ! (Packed fields) Packing accuracy
     &       BLEV,       ! Level
     &       BRLEV,      ! Lower level boundary (Bk for ATMOS)   
     &       BHLEV,      ! (Hybrid levels) A-level of value
     &       BHRLEV,     ! Lower level boundary (Ak for ATMOS)   
     &       BPLAT,      ! Real latitude of 'pseudo' N Pole
     &       BPLON,      ! Real longitude of 'pseudo' N Pole
     &       BGOR,       ! Grid orientation
     &       BZY,        ! Zeroth latitude
     &       BDY,        ! Latitude interval
     &       BZX,        ! Zeroth longitude
     &       BDX,        ! Longitude interval
     &       BMDI,       ! Missing data indicator
     &       BMKS        ! M,K,S scaling factor

C Mapping of MPP_LOOKUP; analogous to mapping in PP header

      INTEGER
     &       P_NADDR,    ! Address on local PE
     &       P_LBLREC    ! Local length of record

      PARAMETER (
     &       P_NADDR=1,
     &       P_LBLREC=2)
C*----------------------------------------------------------------------
C NADDR IS LOCATION IN PP-HEADER (LOOKUP) FOR START POSN OF VARIABLE
C ITEM_CODE is the location in PP header for a code defined as
C           (section number)*1000+item number
C DATA_TYPE is the location in the PP header defining data as REAL or
C           INTEGER.
C LBNPT is the location defining the number of points per row
C
      PARAMETER(
C Validity time
     &       LBYR=1,
     &       LBMON=2,
     &       LBDAT=3,
     &       LBHR=4,
     &       LBMIN=5,
     &       LBDAY=6,
C Data time
     &       LBYRD=7,
     &       LBMOND=8,
     &       LBDATD=9,
     &       LBHRD=10,
     &       LBMIND=11,
     &       LBDAYD=12)

      PARAMETER (
     &       LBTIM=13,
     &       LBFT=14,
     &       LBLREC=15,
     &       LBCODE=16,
     &       LBHEM=17,
     &       LBROW=18,
     &       LBNPT=19,
     &       LBEXT=20,
     &       LBPACK=21,
     &       LBREL=22,
     &       LBFC=23,
     &       LBCFC=24,
     &       LBPROC=25,
     &       LBVC=26,
     &       LBRVC=27)

      PARAMETER (
     &       LBEXP=28,
     &       LBEGIN=29,
     &       LBNREC=30,
     &       LBPROJ=31,
     &       LBTYP=32,
     &       LBLEV=33,
     &       LBRSVD1=34,
     &       LBRSVD2=35,
     &       LBRSVD3=36,
     &       LBRSVD4=37,
     &       LBSRCE=38,
     &       DATA_TYPE=39,
     &       NADDR=40,
     &       LBUSER3=41,    
     &       ITEM_CODE=42,
     &       LBPLEV=43,
     &       LBUSER6=44,
     &       MODEL_CODE=45)

      PARAMETER (
     &       BULEV=46,
     &       BHULEV=47, 
     &       BRSVD3=48,
     &       BRSVD4=49,
     &       BDATUM=50,
     &       BACC=51,
     &       BLEV=52,
     &       BRLEV=53,
     &       BHLEV=54,
     &       BHRLEV=55,
     &       BPLAT=56,
     &       BPLON=57,
     &       BGOR=58,
     &       BZY=59,
     &       BDY=60,
     &       BZX=61,
     &       BDX=62,
     &       BMDI=63,
     &       BMKS=64)

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
C*L  WORKSPACE USAGE:-------------------------------------------------
C   DEFINE LOCAL WORKSPACE ARRAYS: None
C
C*---------------------------------------------------------------------
C
CLL   Construct PP header
      IF (samples.GT.0) THEN   ! Indicates a timeseries/trajectory
        OUT_BZX=0.0
        OUT_BDX=0.0
        OUT_BZY=0.0
        OUT_BDY=0.0
      ELSE
        IF (OCEAN) THEN       !   set OUT_BZY,OUT_BZX,OUT_BDY,OUT_BDX fo
          IF (st_grid.EQ.st_uv_grid) THEN
            OUT_BZY=REALHD(3)-REALHD(2)/2.0
            OUT_BZX=REALHD(4)-REALHD(1)/2.0
          ELSEIF (st_grid.EQ.st_tp_grid) THEN
            OUT_BZY=REALHD(3)-REALHD(2)
            OUT_BZX=REALHD(4)-REALHD(1)
          ELSEIF (st_grid.EQ.st_cu_grid) THEN
            OUT_BZY=REALHD(3)-REALHD(2)
            OUT_BZX=REALHD(4)-REALHD(1)/2.0
          ELSEIF (st_grid.EQ.st_cv_grid) THEN
            OUT_BZY=REALHD(3)-REALHD(2)/2.0
            OUT_BZX=REALHD(4)-REALHD(1)
          ENDIF
          IF (REALHD(32).GT.REALHD(29)) THEN !   greater than RMDI
            OUT_BDY=0.0
            OUT_BDX=REALHD(32)
          ELSE
            OUT_BDY=REALHD(2)
            OUT_BDX=REALHD(1)
          ENDIF
        ELSE                 !   set OUT_BZY,OUT_BZX,OUT_BDY,OUT_BDX for
          IF(st_grid.EQ.st_uv_grid.OR.st_grid.EQ.st_cv_grid) THEN
            OUT_BZY=REALHD(3)+REALHD(2)/2.0 ! UV pts
          ELSE
            OUT_BZY=REALHD(3)+REALHD(2) ! Zeroth Lat OUT_BZY
          ENDIF
C
          IF(st_grid.EQ.st_uv_grid.OR.st_grid.EQ.st_cu_grid) THEN
            OUT_BZX=REALHD(4)-REALHD(1)/2.0 !UV points
          ELSE
            OUT_BZX=REALHD(4)-REALHD(1) ! Zeroth Long OUT_BZX
          ENDIF
          OUT_BDX=REALHD(1) ! Long intvl OUT_BDX
          OUT_BDY=-REALHD(2) ! Lat intvl OUT_BDY
        ENDIF
C
C Add on offset for fields not starting from the origin
C
        OUT_BZY=OUT_BZY
     &       +(N_MOSTROW-1)*OUT_BDY
        OUT_BZX=OUT_BZX
     &     +(W_MOSTCOL-1)*OUT_BDX
        IF(OUT_BZX.GE.360.0)
     *     OUT_BZX=OUT_BZX-360.0
C
C If horizontal averaging has been applied to the output field,
C set OUT_BDX and/or OUT_BDY to the full domain extent
C
        mean_code=(GR/block_size)*block_size
        IF (mean_code.EQ.zonal_mean_base .OR.
     &      mean_code.EQ.field_mean_base .OR.
     &      mean_code.EQ.global_mean_base) THEN
          OUT_BDX=REAL(INTHD(6))*REALHD(1)
        ENDIF
        IF (mean_code.EQ.merid_mean_base .OR.
     &      mean_code.EQ.field_mean_base .OR.
     &      mean_code.EQ.global_mean_base) THEN
          OUT_BDY=REAL(INTHD(7))*REALHD(2)
        ENDIF
      ENDIF
C
  999 CONTINUE
      RETURN
      END
