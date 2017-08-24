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
CLL  SUBROUTINE PPHEAD------------------------------------------
CLL
CLL  Creates a 64 word PP header from the the following:-
CLL  1)  PP_XREF (PP cross-reference array record for this sect/item)
CLL  2)  FIXED length header
CLL  3)  INTEGER constants array
CLL  4)  REAL constants array
CLL  5)  Some input arguments
CLL
CLL  Tested under compiler CFT77
CLL  Tested under OS version 5.1
CLL
CLL T.Johns     <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL   3.2  27/05/93  Code for new real missing data indicator. (TCJ)
CLL   3.5  05/06/95  Remove PP_XREF from argument list and call
CLL                  EXPPXI instead, for submodels work. K Rogers
CLL   4.0  12/09/95  LBUSER(3) [PP_INT_HEAD(LBUSER3)] set to 0.
CLL                  LBCODE set to 31300 + 20 (for Gregorian calendar)
CLL                  or 31300 + 23 (for any other calendar type), if
CLL                  the field is a timeseries.
CLL                  BRLEV, BHRLEV, BULEV[BRSVD1] and BHULEV[BRSVD2]
CLL                  contain lower level boundary and upper level bndry
CLL                  information. Above changes agreed by the WGDUM in
CLL                  first half of 1994. Code for new LBEXP experiment
CLL                  name encoding. Also removed RUN_INDIC_OP from arg
CLL                  list as it is called from CHISTORY  (Andy Brady)
CLL  4.0  12/10/95  Set Lookup(model_code) to internal model ident. RTHB
CLL  4.1  18/04/96  RUN_ID now declared in CHISTORY.  RTHBarnes.
CLL  4.1    Apr. 96  Rationalise *CALLs  S.J.Swarbrick
!LL  4.3    14/02/97 Correct bug where ocean models can try to access
!LL                  uninitialised BKH array               P.Burton
!LL  4.5    14/05/98 Put the correct data type into PP header
!LL                                                  P.Burton
!LL  4.5  14/10/97   Set correct packing type for platform
!LL                  Author D.M. Goddard
!LL  4.5    02/09/98 Set Projection No for High Res Global. D. Robinson.
CLL
CLL  Programming standard: U M DOC  Paper NO. 4,
CLL
CLL  Logical components covered: D40
CLL
CLL  Project TASK: C4
CLL
CLL  External documentation  C4
CLL
CLLEND-------------------------------------------------------------

C
C*L  INTERFACE and ARGUMENTS:------------------------------------------
      SUBROUTINE PP_HEAD(
! COMDECK ARGPPX     
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
!
     & PPXI,PPXC,ppxRecs,
     & PPXPTR,    
! End of comdeck
     *    im_ident,FIXHD,INTHD,REALHD,
     1    LEN_FIXHD,LEN_INTHD,LEN_REALHD,IE,IS,GR,
     2    lfullfield,LEVEL,pseudo_level,
     3    samples,start,start_or_verif_time,end_or_data_time,pp_len,
     4    extraw,PP_INT_HEAD,PP_REAL_HEAD,N_COLS_OUT,NUM_WORDS,
     5    LEN_BUF_WORDS,N_ROWS_OUT,NROW_IN,SROW_IN,WCOL_IN,ECOL_IN,
     5    lbproc_comp,
     6    sample_prd,FCST_PRD,COMP_ACCRCY,PACKING_TYPE,
     7    st_grid,IWA,AK,BK,AKH,BKH,T_levels,LevIndex,ROTATE,ELF,
     8    OCEAN,OCN_DZ,OCN_KM,
     9    ICODE,CMESSAGE)
C*----------------------------------------------------------------
      IMPLICIT NONE


      CHARACTER*(80) CMESSAGE !OUT OUT MESSAGE FROM ROUTINE
C
      LOGICAL
     *  start         ! IN flag to control update for verif/start time
     *, OCEAN          !IN TRUE if processing an ocean diagnostic
     *, lfullfield     !IN TRUE if output field on full horiz domain
C
      INTEGER
     *  start_or_verif_time(7) ! IN verif time/start time for means etc
     *, end_or_data_time(7)    ! IN data time/end time for means etc
     *, samples                ! IN no of samples in period (timeseries)
C
      INTEGER
     *  ICODE             !IN    Return code from the routine
     *, im_ident          !IN    Internal model identifier
     *, PP_LEN            !IN    Length of the lookup table
     *, LEN_FIXHD         !IN    Length of the Fixed Length constants
     *, LEN_INTHD         !IN    Length of the Integer Constants
     *, LEN_REALHD        !IN    Length of the Real Constants
     *, FIXHD(LEN_FIXHD)  !IN    Array of Fixed Constants
     *, INTHD(LEN_INTHD)  !IN    Array of Integer Constants
     *, OCN_KM            !IN    number of ocean model levels
C
      INTEGER
     *  st_grid           !IN    STASH horizontal grid type
     *, T_levels          !IN    No of model Press/Temp levels
     *, LevIndex          !IN    level index  
     *, N_ROWS_OUT        !IN    PPHORIZ_OUT=N_ROWS_OUT*N_COLS_OUT+extra
     *, N_COLS_OUT        !IN    PPHORIZ_OUT=N_COLS_OUT*N_ROWS_OUT+extra
     *, NROW_IN,SROW_IN   !IN    The most nrthrly/southerly row.
     *, WCOL_IN,ECOL_IN   !IN    The most westerly/easterly column
     *, pseudo_level      !IN    Output PP pseudo-level
     *, NUM_OUT           !IN    Number of compressed (32 BIT) words
     *, COMP_ACCRCY       !IN    PACKING ACCURACY IN POWER OF 2
     *, PACKING_TYPE      !IN   0 = No packing, 1 = WGDOS, 3 = GRIB
      INTEGER
     *  U_ROWS            !IN    NO OF U,V, ROWS
     *, P_ROWS            !IN    PRESS/TEMP ROWS
     *, NUM_WORDS         !IN    Number of 64 Bit words to hold DATA
     &, extraw            !IN    Number of extra-data words
     *, LEN_BUF_WORDS     !IN    Number of 64 Bit words (rounded to 512)
     *, IWA               !IN    Start word address.
     *, IE                !IN    Item Number
     *, IS                !IN    Section Number
     *, GR                !IN    Grid point code
     *, FCST_PRD          !IN    Forecast period
     *, LBPROC_COMP(14)   !IN    Subcomponents(0/1) to make up LBPROC
     *, PP_INT_HEAD(PP_LEN)          !OUT  Integer Lookup table
C
      REAL
     *  PP_REAL_HEAD(PP_LEN)!OUT Real Lookup table
     *, REALHD(LEN_REALHD)  !IN  Real header
     *, LEVEL               !IN  Output PP level(REAL)
     *, sample_prd          !IN  Sampling period in hours for time mean
     *, AK(T_levels)        !IN  Hybrid coord Ak at full level
     *, BK(T_levels)        !IN  Hybrid coord Bk at full level
     *, AKH(T_levels+1)     !IN  Hybrid coord Ak at half level
     *, BKH(T_levels+1)     !IN  Hybrid coord Bk at half level   
     *, OCN_DZ(OCN_KM)      !IN  ocean depths at KM levels
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
! COMDECK PPXLOOK
! Description:
!
!   Declares ppxref look-up arrays used by the UM and associated
!    arrays and parameters.
!   Comdecks CSUBMODL,CPPXREF must be *CALLed before this         
!    comdeck
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       May. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.0       Dec. 95   Replace dynamic dim ppxRecs with 
!                     NUM_DIAG_MAX in PPXC   N. Farnon
! 4.1       July 96   *CALL VERSION introduced - NUM_DIAG_MAX made
!                      equal to NDIAGP.
!                     NUM_USR_DIAG_MAX increased from 200 to 300
!                      (just in case).    
! 4.4       03/11/97  Removed MKPPXRF *DEF references. K Rogers
! 4.4       04/11/97  Changed -RECON def line to allow for other small
!                     execs which had used the RECON def. K Rogers
!
! Declarations:

! Global parameters:
! COMDECK VERSION
! Description:                                                          
!   STASH parameter definitions
!                                                                       
! Current code owner: S.J.Swarbrick                                     
!                                                                       
! History:                                                              
! Version   Date      Comment                                           
! -------   ----      -------                                           
! 3.5       Mar. 95   Original code.  S.J.Swarbrick                     
! 4.0                                 S.J.Swarbrick
! 4.1       Apr. 96   Rationalise MDI  S.J.Swarbrick
!  4.1  29/02/96  Increase OUTFILE_E.  RTHBarnes.
!  4.2  27/11/96  MPP code : Increase NELEMP   P.Burton
!  4.3  04/06/97  Increase NELEMP for D1 addressing S.D.Mullerworth
!  4.4  04/06/97  Increase NELEMP for sampling offset. S.D.Mullerworth 
!  4.5  28/01/98  Increade NELEMP for MPP code.   P.Burton
!  4.5  18/09/98  Modify name of VERSION common block to stop potential
!                 clashes with Fortran variable names          P.Burton
!  4.5  30/09/98  Increase NRECDP from 600 to 800. D. Robinson.
!                                                                       
! Declarations:                                                         
      INTEGER   NSECTP           !Max. no. of STASH sections            
      PARAMETER(NSECTP=99)       ! per internal model (44 in practise)  
      INTEGER   NITEMP           !Max. no. of STASH items per           
      PARAMETER(NITEMP=512)      ! section                              
      INTEGER   NRECDP           !Max. no. of STASH list records        
      PARAMETER(NRECDP=800)      ! (prognostic + diagnostic)
      INTEGER   NTIMEP           !Max. no. of output times tables       
      PARAMETER(NTIMEP=100)      ! in STASHC                            
      INTEGER   NPROFTP          !Max. no. of time profiles 
      PARAMETER(NPROFTP=100)     ! in STASHC                            
      INTEGER   NPROFDP          !Max. no. of domain profiles/levels    
      PARAMETER(NPROFDP=100)     ! lists in STASHC (used for both)      
      INTEGER   NTimSerP         !Max. total no. of time series 
      PARAMETER(NTimSerP=1500)    ! in STASHC
      INTEGER   tsdp             !Max. no. time series per
      PARAMETER(tsdp=250)        ! domain profile 
      INTEGER   NPROFUP          !Max. no. of useage profiles           
      PARAMETER(NPROFUP=40)      ! in STASHC                            
      INTEGER   NLEVP            !Max. no. of levels in a               
      PARAMETER(NLEVP=50)        ! levels list                          
      INTEGER   NPSLEVP          !Max. no. of pseudo levels in a        
      PARAMETER(NPSLEVP=40)      ! pseudo levels list                   
      INTEGER   NPSLISTP         !Max. no. of pseudo levels lists       
      PARAMETER(NPSLISTP=40)     ! in STASHC                            
      INTEGER   NDIAGP           !Max. no. non-blank records in         
      PARAMETER(NDIAGP=1800)     ! PPXREF file  
      INTEGER   NDIAGPM          !Same as NRECDP                        
      PARAMETER(NDIAGPM=NRECDP)  ! (will be tidied)                     
      INTEGER   NELEMP           !No. of elements in a ppxref record    
      PARAMETER(NELEMP=33)
      INTEGER   NLEVP_S
      PARAMETER(NLEVP_S=NLEVP*6+1)
      INTEGER   NLEVLSTSP
      PARAMETER(NLEVLSTSP=NPROFDP)                                      
      INTEGER   NMEANP           !No. of meaning periods 
      PARAMETER(NMEANP=4)
! OUTFILE_S, OUTFILE_L and OUTFILE_E must be consistent with
! NUNITS and NUNITS_LEN in comdeck CHSUNITS.
      INTEGER   OUTFILE_S        !Range of                              
      PARAMETER(OUTFILE_S=20)    ! output file                          
      INTEGER   OUTFILE_E        ! numbers                              
      PARAMETER(OUTFILE_E=149)   ! 
      INTEGER   OUTFILE_L
      PARAMETER(OUTFILE_L=OUTFILE_E-OUTFILE_S+1)
!Global scalar:    
      CHARACTER*55 STASH_SET     !Names of stasets files                
!Common block:                                                          
      COMMON/common_VERSION/ STASH_SET
C-----------------------------------------------------------------

! No. of STASH items per section
      INTEGER      PPXREF_ITEMS
        PARAMETER (PPXREF_ITEMS    =NITEMP) 
! No. of STASH sections per internal model
      INTEGER      PPXREF_SECTIONS
        PARAMETER (PPXREF_SECTIONS =NSECTP-55)    
! Max. number of non-null records in ppxref file (>1200) 
      INTEGER      NUM_DIAG_MAX
        PARAMETER (NUM_DIAG_MAX    =NDIAGP)
! Max. number of user-defined ppxref records allowed
      INTEGER      NUM_USR_DIAG_MAX
        PARAMETER (NUM_USR_DIAG_MAX=300)

! No. of ppxref records read into PPXI,PPXC (for dyn. allocation)
      INTEGER      ppxRecs

! Global arrays:
! ppxref look-up arrays
      INTEGER   PPXI(ppxRecs,PPXREF_CODELEN)
      CHARACTER PPXC(NUM_DIAG_MAX,PPXREF_CHARLEN)
! Arrays for temporary storage of user-ppxref records -
!   used to transfer these records from STASH_PROC into U_MODEL
      INTEGER   PPXI_U(NUM_USR_DIAG_MAX,PPXREF_CODELEN)
      CHARACTER PPXC_U(NUM_USR_DIAG_MAX,PPXREF_CHARLEN)
! Array of flags to indicate origin of ppxref record
! 'P' for ppxref file; 'U' for user-stash master file
      CHARACTER OriginFlag(NUM_DIAG_MAX)
! Array of indices to identify which ppxref record corresponds to
!   any given row of PPXI, PPXC
      INTEGER   RowIndex(NUM_DIAG_MAX)
! Pointer array for PPXI, PPXC arrays  
      INTEGER PPXPTR
     & (N_INTERNAL_MODEL    ,0:PPXREF_SECTIONS ,PPXREF_ITEMS)

! Common block:
      COMMON/PPX_INT/ RowIndex,PPXI_U
      COMMON/PPX_CHA/ OriginFlag,PPXC_U
! - End --------------------------------------------------------------
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------
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
C*L --------------------- Comdeck: CHISTORY ----------------------------
CLL
CLL  Purpose: COMMON block for history data needed by top level (C0)
CLL           routines, and passed from run to run.  Mostly set by
CLL           the User Interface.
CLL
CLL           Note that CHISTORY *CALLs ALL individual history comdecks
CLL
CLL  Author : A. Sangster
CLL
CLL  Model            Modification history
CLL version  Date
CLL  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
CLL                 contents.  RTHBarnes.
CLL
CLL  Documentation:  Unified Model Documentation Paper
CLL                  H- History Bricks
CLLEND----------------------------------------------------------------
C*
CCC   *CALL CHSUNITS
! ----------------------- Comdeck: IHISTO   ----------------------------
! Description: COMDECK defining Integer History variables for the
!              model overall.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      INTEGER
     & MODEL_DATA_TIME(6),     ! Array containing model data time
!                              ! (Same as MODEL_BASIS_TIME/MODEL
!                              !  ANALYSIS_HRS depending whether
!                              !  before/after assimilation)
     & RUN_MEANCTL_RESTART,    ! Indicator for next mean period
!                              ! to be processed
     & RUN_INDIC_OP            ! Indicator of operational run type
C
      INTEGER
     & RUN_RESUBMIT_TARGET(6), ! Final target date for the run
C
     & FT_LASTFIELD(20:NUNITS) ! Last field written/read per FT unit
C
C
C History Common Block for overall model integers variables.
C
      COMMON /IHISTO/
     & MODEL_DATA_TIME,
     & RUN_MEANCTL_RESTART, RUN_INDIC_OP,
     & RUN_RESUBMIT_TARGET, FT_LASTFIELD
C
      NAMELIST /NLIHISTO/
     & MODEL_DATA_TIME,
     & RUN_MEANCTL_RESTART, RUN_INDIC_OP,
     & RUN_RESUBMIT_TARGET, FT_LASTFIELD
C
! ----------------------- Comdeck: CHISTO   ----------------------------
! Description: COMDECK defining Character History variables for the
!              model overall.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.1  18/04/96  Add RUN_IN for qxhistreport.  RTHBarnes.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      CHARACTER*10 RUN_HIST_TYPE       ! Type of history file
      CHARACTER*8  RUN_TYPE            ! Type of run
      CHARACTER*14 RUN_COMPCODE        ! Run completion code
      CHARACTER*14 RUN_LAST_MEAN       ! Last mean dump created by run
C APPREARS UNUSED                      ! for pp fields
      CHARACTER*1  RUN_MEANS_TO_DO     ! Flag indicating the run stopped
C                                      ! before creating next mean dump
      CHARACTER*1  RUN_OCEAN_FIRST     ! Flag set to true if ocean to be
C                                      ! run first
      CHARACTER*8  RUN_JOB_NAME        ! Jobname this run
      CHARACTER*5  RUN_ID              ! Expt./Job id for this run
      CHARACTER*1  RUN_RESUBMIT        ! Flag controlling auto resubmit
      CHARACTER*12 RUN_RESUBMIT_Q      ! Job queue to which resubmit run
      CHARACTER*20 RUN_RESUBMIT_TIME   ! Time at which run resubmits
      CHARACTER*6  RUN_RESUBMIT_CPU    ! Time limit for resubmitted job
      CHARACTER*6  RUN_RESUBMIT_MEMORY ! Resubmitted job's memory limit
      CHARACTER*2  RUN_RESUBMIT_PRTY   ! Resubmitted job intra q prty
      CHARACTER*8  RUN_RESUBMIT_JOBNAME! Resubmitted jobname
      CHARACTER*1  FT_ACTIVE(20:NUNITS) ! "Y" if file partly written
C
C
C History Common Block for overall model character variables.
C
      COMMON /CHISTO/
     & RUN_HIST_TYPE, RUN_TYPE, RUN_COMPCODE, RUN_LAST_MEAN,
     & RUN_MEANS_TO_DO, RUN_OCEAN_FIRST, RUN_JOB_NAME, RUN_ID, 
     & RUN_RESUBMIT, RUN_RESUBMIT_Q, RUN_RESUBMIT_TIME,
     & RUN_RESUBMIT_CPU, RUN_RESUBMIT_MEMORY, RUN_RESUBMIT_PRTY,
     & RUN_RESUBMIT_JOBNAME, FT_ACTIVE
C
      NAMELIST /NLCHISTO/
     & RUN_HIST_TYPE, RUN_TYPE, RUN_COMPCODE, RUN_LAST_MEAN,
     & RUN_MEANS_TO_DO, RUN_OCEAN_FIRST, RUN_JOB_NAME, RUN_ID, 
     & RUN_RESUBMIT, RUN_RESUBMIT_Q, RUN_RESUBMIT_TIME,
     & RUN_RESUBMIT_CPU, RUN_RESUBMIT_MEMORY, RUN_RESUBMIT_PRTY,
     & RUN_RESUBMIT_JOBNAME, FT_ACTIVE

! ----------------------- Comdeck: IHISTG   ----------------------------
! Description: COMDECK defining Integer History variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      INTEGER
     & LENGTH(N_INTERNAL_MODEL_MAX) ! No. of tsteps completed this run
     &,ACTUAL_ENDT(6,N_INTERNAL_MODEL_MAX) ! Model end time this run
!     These 2 appears to be purely diagnostic, and not really used.

     &,H_STEPim(N_INTERNAL_MODEL_MAX) ! History block copy of A/O_STEP
!                                   ! held in COMDECK CTIME
     &,H_GROUPim(N_INTERNAL_MODEL_MAX) ! No of steps in coupling period
     &,MEAN_OFFSETim(N_INTERNAL_MODEL_MAX) ! No of means activated
     &,OFFSET_DUMPSim(N_INTERNAL_MODEL_MAX) ! Offset between
!                 MEAN_REFTIME and model basis time  (in model dumps)
     &,MEAN_NUMBERim(N_INTERNAL_MODEL_MAX)  ! No of mean periods chosen
     &,RUN_MEANCTL_INDICim(4,N_INTERNAL_MODEL_MAX) ! Indicators used to
!    correct logical units are used for atmos/ocean partial sum dump I/O
!
C
C History Common Block for generic model integer variables.
C
      COMMON /IHISTG/
     & H_STEPim, H_GROUPim, MEAN_OFFSETim, OFFSET_DUMPSim,
     & MEAN_NUMBERim, RUN_MEANCTL_INDICim
C
      NAMELIST /NLIHISTG/
     & H_STEPim, H_GROUPim, MEAN_OFFSETim, OFFSET_DUMPSim,
     & MEAN_NUMBERim, RUN_MEANCTL_INDICim
C
! ----------------------- Comdeck: CHISTG   ----------------------------
! Description: COMDECK defining Character History variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.4  30/05/97  Added vars LASTATMim, CURRATMim, LASTDMPim.  K Rogers
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      CHARACTER*14 END_DUMPim(N_INTERNAL_MODEL_MAX)!most recent dumpname
      CHARACTER*80 RESTARTim(N_INTERNAL_MODEL_MAX) !current restart dump
      CHARACTER*14 SAFEDMPim(N_INTERNAL_MODEL_MAX) 
! Name of old safe restart dump
      CHARACTER*14 NEWSAFEim(N_INTERNAL_MODEL_MAX)  
! Name of new safe restart dump
      CHARACTER*14 LASTATMim(N_INTERNAL_MODEL_MAX) ! Keep name of last
!                                                  ! atmos restart dump
!                                                  ! until ocean dump
      CHARACTER*14 CURRATMim(N_INTERNAL_MODEL_MAX) ! Keep name of 
!                                                  ! current atmos 
!                                                  ! restart dump
      CHARACTER*14 LASTDMPim(N_INTERNAL_MODEL_MAX) ! Keep name of last
!                                                  ! atmos/ocean dumps
!                                                  ! until meaning done

!
!
! History Common Block for generic model characters variables.
!
      COMMON /CHISTG/
     & END_DUMPim, RESTARTim
     & , SAFEDMPim, NEWSAFEim, LASTATMim, CURRATMim, LASTDMPim
C
      NAMELIST /NLCHISTG/
     & END_DUMPim, RESTARTim
     & , SAFEDMPim, NEWSAFEim, LASTATMim, CURRATMim, LASTDMPim
C
C*L --------------------- Comdeck: CLFHIST  ----------------------------
CLL
CLL  Purpose: COMDECK defining unit numbers relevant to history file
CLL           and variables used to hold the logical to physical
CLL           file associations made within the model
CLL
CLL  Author : A. Sangster
CLL
CLL  Documentation:  Unified Model Documentation Paper
CLL                  H- History Bricks
CLL                  Version 5  18/6/90
CLL
CLL  Model             Modification history from model version 3.0
CLL version  Date
CLL
CLL  3.4  30/09/94  Add files MURKFILE,OUSRANCL,OUSRMULT at 109,113,114
CLL  3.4  05/09/94  Add files USRANCIL,USRMULTI at unit nos. 111,112.
CLL
CLL  3.3  22/11/93  Add file SOURCES at unit number 110. R.T.H.Barnes.
CLL 3.2     28/05/93  Add file BAS_IND at unit number 58. M.Carter.
CLL  Vn3.0  12/02/93 - Variables PERTURB and TRANSP equivalenced to unit
CLL                    numbers 37, and 97 respectively. C.S. Douglas
CLL  3.4  1/8/94     Revised Obs file specification: Stuart Bell
CLL  3.5  01/05/95  Sub-models stage 1: History/control files. RTHBarnes
!    4.0  22/09/95  Added units for Spectral data for Radiation scheme.
!                                        (J. M. Edwards)
CLL  4.1  11/03/96  Introduce Wave sub-model.  RTHBarnes.
!    4.1  26/02/96  Associate new env. variables SO2NATEM and CHEMOXID
!                   with unit nos. 115 & 116. Rename SOURCES to 
!                   SULPEMIS. D. Robinson.
!  4.3   18/3/97  Add aerosol forcings of climate change.  Will Ingram
!  4.4   4/7/97   Add ANLINCR  Chris Jones/Stuart Bell
CLL  4.4   12/9/97  Associate ancillary file EVs for initial surface 
CLL                 type fracs, initial vegetation state and vegetation 
CLL                 disturbance with unit no.s 135-137 R. Betts
CLL  4.4  17/10/97  Associate env var. CACHED with Unit 138. D Robinson 
CLL  4.5  22/04/98  Add new ancillary file for soot emissions:  
CLL                 SOOTEMIS - in I/O unit 139. R.Rawlins
CLL  4.5  29/07/98  Add new variables ALABCOU5/6/7/8. D. Robinson.
CLL  4.5  17/08/98  Add new variables OLABCOU1/2/3/4. Remove
CLL                 OLABCOUT. D. Robinson.
CLL
CLL  Type declarations
CLL
CLL
CLL  Logical Filenames used in the model
CLL
      CHARACTER*80 HKFILE,PPXREF,CONFIG,STASHCTL,NAMELIST,OUTPUT,
     *             OUTPUT2,MCTL,ICTL,PHIST,IHIST,THIST,FTXX,
     *             CACHE1,CACHE2,ASWAP,OSWAP,AOTRANS,
     2             AINITIAL,ASTART,ARESTART,AOPSUM1,AOPSUM2,AOPSUM3,
     *             AOPSUM4,AOMEAN,SSU,
     3             OZONE,SMCSNOWD,DSOILTMP,SOILTYPE,VEGTYPE,SSTIN,
     *             SICEIN,PERTURB,MASK,
     4             OINITIAL,OSTART,ORESTART,AOPSTMP1,AOPSTMP2,AOPSTMP3,
     *             AOPSTMP4,
     5             WFIN,HFLUXIN,PMEIN,ICEFIN,AIRTMP,
     &             SWSPECTD,
     6             PP0,PP1,PP2,PP3,PP4,PP5,PP6,PP7,PP8,PP9,
     &             OBS01,OBS02,OBS03,OBS04,OBS05,
     &             OBS06,OBS07,OBS08,OBS09,OBS10,
     8             LWSPECTD,WAVEOUT,SURGEOUT,MESOUT,STRATOUT,WFOUT,     
     &          HFLUXOUT,FLXCROUT,PMEOUT,ICEFOUT,MOSOUT,SSTOUT,SICEOUT,
     *             CURNTOUT,ALABCIN,OROG,OLABCIN,OCNDEPTH,CURNTIN,
     *             FLUXCORR,SLABHCON,ATMANL,OCNANL,BAS_IND
     &             ,TRANSP,ATRACER,OTRACER,SULPEMIS,USRANCIL,USRMULTI,
     *             OUSRANCL,OUSRMULT,MURKFILE,
     *             ALABCOU1,ALABCOU2,ALABCOU3,ALABCOU4
     &            ,ALABCOU5,ALABCOU6,ALABCOU7,ALABCOU8
     &            ,OLABCOU1,OLABCOU2,OLABCOU3,OLABCOU4
     &            ,ANLINCR
     &            ,WINITIAL,WSTART,WRESTART,WAVANL,WAVANCIN  
     &            ,SO2NATEM,CHEMOXID,AEROFCG,FRACINIT,VEGINIT,DISTURB
     &            ,CACHED,SOOTEMIS
     &            ,CO2EMITS
C
      CHARACTER*80 MODEL_FT_UNIT ! Array holding FORTRAN unit file
C                                ! associations details for each unit
C
      INTEGER
     *        MCTL_UNIT,         ! Master control namelist file unit
     *        ICTL_UNIT,         ! Interim control namelist file unit
     *        PHIST_UNIT,        ! Permanent history file unit
     *        IHIST_UNIT,        ! Interim history file unit
     *        THIST_UNIT,        ! Temporary history file unit
     *        FTXX_UNIT,         ! Logical/physical file associations
     *        HKFILE_UNIT        ! Operational houskeeping file unit
C*
C  Parameters specifying unit numbers relevant to control/history tasks
C
      PARAMETER(HKFILE_UNIT= 1)
      PARAMETER(MCTL_UNIT  = 8)
      PARAMETER(ICTL_UNIT  = 9)
      PARAMETER(PHIST_UNIT =10)
      PARAMETER(IHIST_UNIT =11)
      PARAMETER(THIST_UNIT =12)
      PARAMETER(FTXX_UNIT  =13)
!
! Namelist of all permissible logical files.
!
      NAMELIST / NLCFILES /
     &             HKFILE,PPXREF,CONFIG,STASHCTL,NAMELIST,OUTPUT,
     &             OUTPUT2,MCTL,ICTL,PHIST,IHIST,THIST,FTXX,
     &             CACHE1,CACHE2,ASWAP,OSWAP,AOTRANS,
     &             AINITIAL,ASTART,ARESTART,AOPSUM1,AOPSUM2,AOPSUM3,
     &             AOPSUM4,AOMEAN,SSU,
     &             OZONE,SMCSNOWD,DSOILTMP,SOILTYPE,VEGTYPE,SSTIN,
     &             SICEIN,PERTURB,MASK,
     &             OINITIAL,OSTART,ORESTART,AOPSTMP1,AOPSTMP2,AOPSTMP3,
     &             AOPSTMP4,
     &             WFIN,HFLUXIN,PMEIN,ICEFIN,AIRTMP,
     &             SWSPECTD,
     &             PP0,PP1,PP2,PP3,PP4,PP5,PP6,PP7,PP8,PP9,
     &             OBS01,OBS02,OBS03,OBS04,OBS05,
     &             OBS06,OBS07,OBS08,OBS09,OBS10,
     &             LWSPECTD,WAVEOUT,SURGEOUT,MESOUT,STRATOUT,WFOUT,     
     &          HFLUXOUT,FLXCROUT,PMEOUT,ICEFOUT,MOSOUT,SSTOUT,SICEOUT,
     &             CURNTOUT,ALABCIN,OROG,OLABCIN,OCNDEPTH,CURNTIN,
     &             FLUXCORR,SLABHCON,ATMANL,OCNANL,BAS_IND
     &             ,TRANSP,ATRACER,OTRACER,SULPEMIS,USRANCIL,USRMULTI,
     &             OUSRANCL,OUSRMULT,MURKFILE,
     &             ALABCOU1,ALABCOU2,ALABCOU3,ALABCOU4
     &            ,ALABCOU5,ALABCOU6,ALABCOU7,ALABCOU8
     &            ,OLABCOU1,OLABCOU2,OLABCOU3,OLABCOU4
     &            ,ANLINCR
     &            ,WINITIAL,WSTART,WRESTART,WAVANL,WAVANCIN 
     &            ,SO2NATEM,CHEMOXID,AEROFCG,FRACINIT,VEGINIT,DISTURB
     &            ,CACHED,SOOTEMIS
     &            ,CO2EMITS
C
C Common block definition
C
      COMMON/CLFHIST/MODEL_FT_UNIT(NUNITS)
C
C  Equivalence logical filenames within array MODEL_FT_UNIT
C
      EQUIVALENCE
     *(HKFILE    ,MODEL_FT_UNIT(1)  ),(PPXREF     ,MODEL_FT_UNIT(2)  ),
     *(CONFIG    ,MODEL_FT_UNIT(3)  ),(STASHCTL   ,MODEL_FT_UNIT(4)  ),
     *(NAMELIST  ,MODEL_FT_UNIT(5)  ),(OUTPUT     ,MODEL_FT_UNIT(6)  ),
     *(OUTPUT2   ,MODEL_FT_UNIT(7)  ),(MCTL       ,MODEL_FT_UNIT(8)  ),
     *(ICTL      ,MODEL_FT_UNIT(9)  ),(PHIST      ,MODEL_FT_UNIT(10) ),
     *(IHIST     ,MODEL_FT_UNIT(11) ),(THIST      ,MODEL_FT_UNIT(12) ),
     *(FTXX      ,MODEL_FT_UNIT(13) ),
     *(CACHE1    ,MODEL_FT_UNIT(15) ),(CACHE2     ,MODEL_FT_UNIT(16) ),
     *(AOTRANS   ,MODEL_FT_UNIT(17) ),(ASWAP      ,MODEL_FT_UNIT(18) ),
     *(OSWAP     ,MODEL_FT_UNIT(19) ),(AINITIAL   ,MODEL_FT_UNIT(20) ),
     *(ASTART    ,MODEL_FT_UNIT(21) ),(ARESTART   ,MODEL_FT_UNIT(22) ),
     *(AOPSUM1   ,MODEL_FT_UNIT(23) ),(AOPSUM2    ,MODEL_FT_UNIT(24) ),
     *(AOPSUM3   ,MODEL_FT_UNIT(25) )
C
      EQUIVALENCE
     *(AOPSUM4   ,MODEL_FT_UNIT(26) ),(AOMEAN     ,MODEL_FT_UNIT(27) ),
     *(ATMANL    ,MODEL_FT_UNIT(28) ),(SSU        ,MODEL_FT_UNIT(29) ),
     *(OZONE     ,MODEL_FT_UNIT(30) ),(SMCSNOWD   ,MODEL_FT_UNIT(31) ),
     *(DSOILTMP  ,MODEL_FT_UNIT(32) ),(SOILTYPE   ,MODEL_FT_UNIT(33) ),
     *(VEGTYPE   ,MODEL_FT_UNIT(34) ),(SSTIN      ,MODEL_FT_UNIT(35) ),
     *(SICEIN    ,MODEL_FT_UNIT(36) ),(PERTURB    ,MODEL_FT_UNIT(37) ),
     *(CURNTIN   ,MODEL_FT_UNIT(38) ),(MASK       ,MODEL_FT_UNIT(39) ),
     *(OINITIAL  ,MODEL_FT_UNIT(40) ),(OSTART     ,MODEL_FT_UNIT(41) ),
     *(ORESTART  ,MODEL_FT_UNIT(42) ),(AOPSTMP1   ,MODEL_FT_UNIT(43) ),
     *(AOPSTMP2  ,MODEL_FT_UNIT(44) ),(AOPSTMP3   ,MODEL_FT_UNIT(45) ),
     *(AOPSTMP4  ,MODEL_FT_UNIT(46) ),(OCNANL     ,MODEL_FT_UNIT(47) ),
     *(ATRACER   ,MODEL_FT_UNIT(48) ),(OTRACER    ,MODEL_FT_UNIT(49) ),
     *(WFIN      ,MODEL_FT_UNIT(50) )
C
      EQUIVALENCE
     *(HFLUXIN   ,MODEL_FT_UNIT(51) ),(PMEIN      ,MODEL_FT_UNIT(52) ),
     *(ICEFIN    ,MODEL_FT_UNIT(53) ),(AIRTMP     ,MODEL_FT_UNIT(54) ),
     *                                (FLUXCORR   ,MODEL_FT_UNIT(56) ),
     *(SWSPECTD  ,MODEL_FT_UNIT(57) ),(BAS_IND    ,MODEL_FT_UNIT(58) ), 
     *(SLABHCON  ,MODEL_FT_UNIT(59) ),(PP0        ,MODEL_FT_UNIT(60) ),
     *(PP1       ,MODEL_FT_UNIT(61) ),(PP2        ,MODEL_FT_UNIT(62) ),
     *(PP3       ,MODEL_FT_UNIT(63) ),(PP4        ,MODEL_FT_UNIT(64) ),
     *(PP5       ,MODEL_FT_UNIT(65) ),(PP6        ,MODEL_FT_UNIT(66) ),
     *(PP7       ,MODEL_FT_UNIT(67) ),(PP8        ,MODEL_FT_UNIT(68) ),
     &(PP9       ,MODEL_FT_UNIT(69) ),(OBS01      ,MODEL_FT_UNIT(70) ),
     &(OBS02     ,MODEL_FT_UNIT(71) ),(OBS03      ,MODEL_FT_UNIT(72) ),
     &(OBS04     ,MODEL_FT_UNIT(73) ),(OBS05      ,MODEL_FT_UNIT(74) )
C
      EQUIVALENCE
     &(OBS06     ,MODEL_FT_UNIT(75) ),(OBS07      ,MODEL_FT_UNIT(76) ),
     &(OBS08     ,MODEL_FT_UNIT(77) ),(OBS09      ,MODEL_FT_UNIT(78) ),
     &(OBS10     ,MODEL_FT_UNIT(79) ),(LWSPECTD   ,MODEL_FT_UNIT(80) ), 
     *(WAVEOUT   ,MODEL_FT_UNIT(81) ),(SURGEOUT   ,MODEL_FT_UNIT(82) ),
     *(MESOUT    ,MODEL_FT_UNIT(83) ),(STRATOUT   ,MODEL_FT_UNIT(84) ),
     *(WFOUT     ,MODEL_FT_UNIT(85) ),(HFLUXOUT   ,MODEL_FT_UNIT(86) ),
     *(PMEOUT    ,MODEL_FT_UNIT(87) ),(ICEFOUT    ,MODEL_FT_UNIT(88) ),
     &(MOSOUT    ,MODEL_FT_UNIT(89) ),
     *(SSTOUT    ,MODEL_FT_UNIT(91) ),(SICEOUT    ,MODEL_FT_UNIT(92) ),
     *(CURNTOUT  ,MODEL_FT_UNIT(93) ),(FLXCROUT   ,MODEL_FT_UNIT(94) ),
     *(ALABCIN   ,MODEL_FT_UNIT(95) ),(OROG       ,MODEL_FT_UNIT(96) ),
     *(TRANSP    ,MODEL_FT_UNIT(97) ),(OLABCIN    ,MODEL_FT_UNIT(98) ),
     *(OCNDEPTH  ,MODEL_FT_UNIT(99) ),
     &(OLABCOU1  ,MODEL_FT_UNIT(100)),(OLABCOU2   ,MODEL_FT_UNIT(101)),
     &(OLABCOU3  ,MODEL_FT_UNIT(102)),(OLABCOU4   ,MODEL_FT_UNIT(103)),
     &(ANLINCR   ,MODEL_FT_UNIT(108)),(MURKFILE   ,MODEL_FT_UNIT(109)),
     &(SULPEMIS  ,MODEL_FT_UNIT(110)),(USRANCIL   ,MODEL_FT_UNIT(111)),
     *(USRMULTI  ,MODEL_FT_UNIT(112)),(OUSRANCL   ,MODEL_FT_UNIT(113)),
     *(OUSRMULT  ,MODEL_FT_UNIT(114)),(SO2NATEM   ,MODEL_FT_UNIT(115)),
     &(CHEMOXID  ,MODEL_FT_UNIT(116)),(AEROFCG    ,MODEL_FT_UNIT(117)),
     *(CO2EMITS  ,MODEL_FT_UNIT(118)),
     *(WINITIAL  ,MODEL_FT_UNIT(130)),(WSTART     ,MODEL_FT_UNIT(131)), 
     *(WRESTART  ,MODEL_FT_UNIT(132)),(WAVANL     ,MODEL_FT_UNIT(133)), 
     *(WAVANCIN  ,MODEL_FT_UNIT(134)),(FRACINIT   ,MODEL_FT_UNIT(135)),
     *(VEGINIT   ,MODEL_FT_UNIT(136)),(DISTURB    ,MODEL_FT_UNIT(137)),
     &(CACHED    ,MODEL_FT_UNIT(138)),(SOOTEMIS   ,MODEL_FT_UNIT(139)),
     &(ALABCOU1  ,MODEL_FT_UNIT(140)),(ALABCOU2   ,MODEL_FT_UNIT(141)),
     &(ALABCOU3  ,MODEL_FT_UNIT(142)),(ALABCOU4   ,MODEL_FT_UNIT(143)),
     &(ALABCOU5  ,MODEL_FT_UNIT(144)),(ALABCOU6   ,MODEL_FT_UNIT(145)),
     &(ALABCOU7  ,MODEL_FT_UNIT(146)),(ALABCOU8   ,MODEL_FT_UNIT(147)) 
C

        EXTERNAL EXPPXI
        EXTERNAL EXPT_ENC  

C*L  WORKSPACE USAGE:-------------------------------------------------
C   DEFINE LOCAL WORKSPACE ARRAYS: None
C
C*---------------------------------------------------------------------
C    DEFINE LOCAL VARIABLES
      REAL
     *  ocn_depth     !     depth of ocean at level
     *, ocn_depth_h   !     depth of ocean at half level  
      INTEGER
     *  PP_LBFC       !     M08 Level code
     *, PP_LBTYP      !     M08 Field type code
     *, PP_LBLEV      !     M08 Field level code
     *, PP_IPROJ      !     M08 Projection number
     *, PP_LBVC       !     Vertical coord type
     *, II            !     Local Counter
     *, int_level     !     integer value of level
     *, K             !     local counter
     *, IA,IB,IC      !     Component codes to make up LBTIM
     *, mean_code     !     spatial averaging code derived from GR
     *, lvcode        !     lv code
     *, EXPPXI        !     Function to extract ppxref info
     *, EXPTCODE      !     integer coded experiment name  

      LOGICAL
     *  ELF,
     *  ROTATE
C
CLL   Construct PP header
C
C  Timestamps ----------------------------------------------------------
C
CL
CL Set up time info dependent on start flag.
CL For all but time series start will be TRUE so all time information
CL will be set up from FIXHD in effect, but for time series start
CL will be set up by TEMPORAL and passed in, so that dump headers are
CL set correctly for such fields.
CL Note: end_or_data_time will be updated from current model time in
CL       FIXHD(28-34) for time means/accumulations etc.
CL
      IF (start) THEN    ! start timestep so update start time
        PP_INT_HEAD(LBYR)=start_or_verif_time(1)
        PP_INT_HEAD(LBMON)=start_or_verif_time(2)
        PP_INT_HEAD(LBDAT)=start_or_verif_time(3)
        PP_INT_HEAD(LBHR)=start_or_verif_time(4)
        PP_INT_HEAD(LBMIN)=start_or_verif_time(5)
        PP_INT_HEAD(LBDAY)=start_or_verif_time(7)
      ENDIF
      PP_INT_HEAD(LBYRD)=end_or_data_time(1)
      PP_INT_HEAD(LBMOND)=end_or_data_time(2)
      PP_INT_HEAD(LBDATD)=end_or_data_time(3)
      PP_INT_HEAD(LBHRD)=end_or_data_time(4)
      PP_INT_HEAD(LBMIND)=end_or_data_time(5)
      PP_INT_HEAD(LBDAYD)=end_or_data_time(7)
C
C  Secondary time information ------------------------------------------
C
C LBTIM is 100*IA+10*IB+IC - this encodes the time processing type
C
      IA=INT(sample_prd)           ! Sampling period in whole hours
      IF(sample_prd.eq.0.0) THEN   ! NB: may be a fraction of an hour
        IB=1                       ! Forecast field
      ELSE
        IF (IA.EQ.0) THEN
          IA=1                     ! 0 < sample_prd < 1 counts as 1 hour
        ENDIF
        IB=2                       ! Time mean or accumulation
      ENDIF
      IC=FIXHD(8)                  ! Calendar (1: Gregorian, 2: 360 day)
C
      PP_INT_HEAD(LBTIM)=100*IA+10*IB+IC
      PP_INT_HEAD(LBFT)=FCST_PRD
C
C  Data length ---------------------------------------------------------
C
      PP_INT_HEAD(LBLREC)=NUM_WORDS
C
C  Grid code (determined from dump fixed-length header) ----------------
C
      IF (samples.EQ.0) THEN
C       Field is not a timeseries
        IF(FIXHD(4).LT.100) THEN
          PP_INT_HEAD(LBCODE)=1   ! Regular lat/long grid
        ELSE
          PP_INT_HEAD(LBCODE)=101 ! lat/long grid non-std polar axis
        ENDIF
      ELSE
C       Field is a timeseries
        PP_INT_HEAD(LBCODE)=31300
        IF (FIXHD(8).EQ.1) THEN
C         Calendar --  1: Gregorian
          PP_INT_HEAD(LBCODE)=PP_INT_HEAD(LBCODE)+20
        ELSEIF (FIXHD(8).EQ.2) THEN
C         Calendar -- 360 day (Model Calendar)
          PP_INT_HEAD(LBCODE)=PP_INT_HEAD(LBCODE)+23
        ELSE
C         Unknown calendar. Fail.
          ICODE=2
      CMESSAGE='PPHEAD: unknown calender type in fixhd(8)'             
        ENDIF
      ENDIF
C
C  Hemispheric subregion indicator -------------------------------------
C
      IF (samples.GT.0 .OR. .NOT.lfullfield) THEN
C  Field is a timeseries/trajectory or subdomain of the full model area
        PP_INT_HEAD(LBHEM)=3
      ELSEIF (FIXHD(4).LT.100) THEN
C  Otherwise, use the value for the full model area encoded in the dump
        PP_INT_HEAD(LBHEM)=FIXHD(4)
      ELSE
        PP_INT_HEAD(LBHEM)=FIXHD(4)-100
      ENDIF
C
C  Field dimensions (rows x cols) --------------------------------------
C
      PP_INT_HEAD(LBROW)=N_ROWS_OUT
      PP_INT_HEAD(LBNPT)=N_COLS_OUT
C
C  'Extra data' length (now accomodates timeseries sampling data) ------
C
      PP_INT_HEAD(LBEXT)=extraw
C
C  Packing method indicator (new definition introduced at vn2.8)--------
       IF(PACKING_TYPE.EQ.1)THEN    ! WGDOS packing
         PP_INT_HEAD(LBPACK)=00001
       ELSEIF(PACKING_TYPE.EQ.3)THEN ! GRIB packing
         PP_INT_HEAD(LBPACK)=00003
       ELSEIF(PACKING_TYPE.EQ.0)THEN ! No packing
         PP_INT_HEAD(LBPACK)=00000
       ELSE
         ICODE=1
         CMESSAGE='PPHEAD  Packing type undefined'
         PP_INT_HEAD(LBPACK)=00000
      ENDIF
C
C  PP header release no ------------------------------------------------
C
      PP_INT_HEAD(LBREL)=2
C
C  Primary fieldcode (some hardwiring for ELF winds) -------------------
C  Secondary fieldcode not used currently
C
      PP_LBFC=EXPPXI(im_ident, is, ie, ppx_field_code,
! COMDECK ARGPPX     
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
!
     & PPXI,PPXC,ppxRecs,
     & PPXPTR,    
! End of comdeck
     &               icode, cmessage)
      IF(ELF.AND..NOT.ROTATE) THEN  ! ELF winds are in x,y direction
        IF(PP_LBFC.EQ.56) PP_LBFC=48
        IF(PP_LBFC.EQ.57) PP_LBFC=49
      ENDIF
      PP_INT_HEAD(LBFC)=PP_LBFC
      PP_INT_HEAD(LBCFC)=0
C
C  Processing code (encodes several things in one field) ---------------
C
      PP_INT_HEAD(LBPROC)=0
      DO II=14,1,-1
        PP_INT_HEAD(LBPROC)=PP_INT_HEAD(LBPROC)*2+LBPROC_COMP(II)
      ENDDO
C
C  Vertical coordinate type --------------------------------------------
C  Vertical coordinate type for reference level not coded
C
      PP_LBVC=EXPPXI(im_ident, is, ie, ppx_lbvc_code,
! COMDECK ARGPPX     
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
!
     & PPXI,PPXC,ppxRecs,
     & PPXPTR,    
! End of comdeck
     &               icode, cmessage)
      PP_INT_HEAD(LBVC)=PP_LBVC
      PP_INT_HEAD(LBRVC)=0
C
C  Experiment number coded from EXPT_ID and JOB_ID for non
C  operational set to RUN_INDIC_OP for operational use.
C
      IF (MODEL_STATUS.NE.'Operational') THEN
        RUN_ID(1:4)=EXPT_ID
        RUN_ID(5:5)=JOB_ID
C  Function EXPT_ENC will encode the run_id into a unique integer
        CALL EXPT_ENC(RUN_ID,EXPTCODE,ICODE,CMESSAGE)
C  We do not return here. We wait until the end of the subroutine.
        PP_INT_HEAD(LBEXP)=EXPTCODE          ! LBEXP
      ELSE
        PP_INT_HEAD(LBEXP)=RUN_INDIC_OP      ! LBEXP (ITAB)
      ENDIF 
C
C  Direct access dataset start address and no of records ---------------
C
      PP_INT_HEAD(LBEGIN)=IWA
      PP_INT_HEAD(LBNREC)=LEN_BUF_WORDS
C
C  Operational fieldsfile projection no, fieldtype + level codes -------
C  These are hardwired according to model resolution
C
      IF(INTHD(6).EQ.192) THEN
        PP_IPROJ=802
      ELSE IF(INTHD(6).EQ.288) THEN
        PP_IPROJ=800
      ELSE IF(INTHD(6).EQ.96) THEN
        PP_IPROJ=870
       ELSE IF(INTHD(6).EQ.432) THEN
         PP_IPROJ=800
      ELSE
        PP_IPROJ=900
      ENDIF
      PP_LBTYP=EXPPXI(im_ident, is, ie, ppx_meto8_fieldcode,
! COMDECK ARGPPX     
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
!
     & PPXI,PPXC,ppxRecs,
     & PPXPTR,    
! End of comdeck
     &               icode, cmessage)
      lvcode=EXPPXI(im_ident, is, ie, ppx_lv_code,
! COMDECK ARGPPX     
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
!
     & PPXI,PPXC,ppxRecs,
     & PPXPTR,    
! End of comdeck
     &             icode, cmessage)
      IF(LEVEL.EQ.-1.0) THEN
        PP_LBLEV=EXPPXI(im_ident, is, ie, ppx_meto8_levelcode,
! COMDECK ARGPPX     
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
!
     & PPXI,PPXC,ppxRecs,
     & PPXPTR,    
! End of comdeck
     &               icode, cmessage)   ! levelcode 9999 or 8888
      ELSE
        IF (im_ident .eq. atmos_im) THEN
          IF (lvcode.eq.ppx_half_level.and.BKH(LevIndex).eq.1.0) THEN
! NB: If BK indicates surface hybrid level, reset LBLEV to correspond
            PP_LBLEV=ppx_meto8_surf
          ELSE
            PP_LBLEV=LEVEL+0.00001
          ENDIF
        ELSE
          PP_LBLEV=LEVEL+0.00001
        ENDIF
      ENDIF
      PP_INT_HEAD(LBPROJ)=PP_IPROJ
      PP_INT_HEAD(LBTYP)=PP_LBTYP
      PP_INT_HEAD(LBLEV)=PP_LBLEV
C
C  Reserved slots for future expansion ---------------------------------
C
      PP_INT_HEAD(LBRSVD1)=0
      PP_INT_HEAD(LBRSVD2)=0
      PP_INT_HEAD(LBRSVD3)=0
      PP_INT_HEAD(LBRSVD4)=0
C
C  Spare for user's use ------------------------------------------------
C
      PP_INT_HEAD(LBSRCE)=1111
C
! Data type - extract from PPXREF
      PP_INT_HEAD(DATA_TYPE)=EXPPXI(im_ident, is, ie, ppx_data_type,
! COMDECK ARGPPX     
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
!
     & PPXI,PPXC,ppxRecs,
     & PPXPTR,    
! End of comdeck
     &                              icode, cmessage)
C
C  Address within dump or PP file --------------------------------------
C
      PP_INT_HEAD(NADDR)=IWA
C
C  LBUSER3 is not currently used (ie set to 0).   
C
      PP_INT_HEAD(LBUSER3)=0  
C
C  STASH section/item code ---------------------------------------------
C
      PP_INT_HEAD(ITEM_CODE)=IS*1000+IE
C
C  STASH pseudo-level (for fields which have pseudo-levels defined) ----
C
      PP_INT_HEAD(LBPLEV)=pseudo_level
C
C  Spare for user's use ------------------------------------------------
C
      PP_INT_HEAD(LBUSER6)=0
      PP_INT_HEAD(MODEL_CODE) = im_ident
C
C  Reserved for future PP package use ----------------------------------
C
      PP_REAL_HEAD(BRSVD3)=0.0
      PP_REAL_HEAD(BRSVD4)=0.0
      PP_REAL_HEAD(BDATUM)=0.0
      PP_REAL_HEAD(BACC)=COMP_ACCRCY ! packing accuracy stored as real
C
C  Vertical grid description -------------------------------------------
C  Level and reference level
C
      IF(PP_LBVC.GE.126.AND.PP_LBVC.LE.139) THEN ! Special codes
C                                                  (surf botttom,
C                                                   top all zero)
        PP_REAL_HEAD(BLEV)=0.0
        PP_REAL_HEAD(BHLEV)=0.0
        PP_REAL_HEAD(BRLEV)=0.0
        PP_REAL_HEAD(BHRLEV)=0.0
        PP_REAL_HEAD(BULEV)=0.0
        PP_REAL_HEAD(BHULEV)=0.0
      ELSEIF(PP_LBVC.EQ.9) THEN      ! Hybrid/ETA levels  
        lvcode=EXPPXI(im_ident, is, ie, ppx_lv_code,
! COMDECK ARGPPX     
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
!
     & PPXI,PPXC,ppxRecs,
     & PPXPTR,    
! End of comdeck
     &    icode, cmessage)
        IF (lvcode.EQ.ppx_half_level) THEN ! model levels
          PP_REAL_HEAD(BLEV)=BKH(LevIndex)
          PP_REAL_HEAD(BHLEV)=AKH(LevIndex)
          IF(LevIndex.eq.1) THEN
C           This case for surface eta diags. Halflevel below
C           surface does not exist.
            PP_REAL_HEAD(BRLEV)=BKH(LevIndex)
            PP_REAL_HEAD(BHRLEV)=AKH(LevIndex)
          ELSE
            PP_REAL_HEAD(BRLEV)=BK(LevIndex-1)
            PP_REAL_HEAD(BHRLEV)=AK(LevIndex-1)
          ENDIF                               
          IF(LevIndex.eq.T_levels+1) THEN 
C           This case for eta diags at top of atmosphere.
C           Half level above toa does not exist.
            PP_REAL_HEAD(BULEV)=BKH(LevIndex)
            PP_REAL_HEAD(BHULEV)=AKH(LevIndex)
          ELSE
            PP_REAL_HEAD(BULEV)=BK(LevIndex)
            PP_REAL_HEAD(BHULEV)=AK(LevIndex)
          ENDIF 
        ELSE                                             ! half levels
          PP_REAL_HEAD(BLEV)=BK(LevIndex)
          PP_REAL_HEAD(BHLEV)=AK(LevIndex)
          PP_REAL_HEAD(BRLEV)=BKH(LevIndex)
          PP_REAL_HEAD(BHRLEV)=AKH(LevIndex)
          PP_REAL_HEAD(BULEV)=BKH(LevIndex+1)
          PP_REAL_HEAD(BHULEV)=AKH(LevIndex+1)
        ENDIF
      ELSEIF (PP_LBVC.EQ.2.AND.OCEAN) THEN ! Depth levels
        PP_REAL_HEAD(BHRLEV)=0.0
        PP_REAL_HEAD(BHULEV)=0.0
        int_level=level
        lvcode=EXPPXI(im_ident, is, ie, ppx_lv_code,
! COMDECK ARGPPX     
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
!
     & PPXI,PPXC,ppxRecs,
     & PPXPTR,    
! End of comdeck
     &               icode, cmessage)
C ocn_depth defined for ocean full levels (e.g. temperature), and
C ocn_depth_h for ocean half-levels (e.g. vertical velocity)
        ocn_depth=0.5*OCN_DZ(1)
        ocn_depth_h=0.
        IF (int_level.GT.1) THEN
          DO K=2,int_level
C           Loop over levels calculating half levels as we go.
            ocn_depth=ocn_depth+0.5*(OCN_DZ(K-1)+OCN_DZ(K))

            ocn_depth_h=ocn_depth_h+OCN_DZ(K-1)
          END DO
        ENDIF
        IF (lvcode.EQ.ppx_half_level) THEN
          PP_REAL_HEAD(BLEV)=ocn_depth_h
          PP_REAL_HEAD(BRLEV)=ocn_depth
          IF (int_level.EQ.1) THEN
            PP_REAL_HEAD(BULEV)=0.0 ! This level would be a
C                                     half level above the ocean.
C                                     Set to zero.
          ELSE
            PP_REAL_HEAD(BULEV)=ocn_depth_h-0.5*OCN_DZ(int_level-1)
          ENDIF
        ELSE
          PP_REAL_HEAD(BLEV)=ocn_depth
          PP_REAL_HEAD(BRLEV)=ocn_depth_h+OCN_DZ(int_level)
          PP_REAL_HEAD(BULEV)=ocn_depth_h
        ENDIF

        PP_REAL_HEAD(BHLEV)=0.0
      ELSE
        PP_REAL_HEAD(BLEV)=LEVEL
        PP_REAL_HEAD(BHLEV)=0.0
        PP_REAL_HEAD(BRLEV)=0.0  ! The boundary levels
        PP_REAL_HEAD(BHRLEV)=0.0 ! are not known
        PP_REAL_HEAD(BULEV)=0.0  ! for pressure
        PP_REAL_HEAD(BHULEV)=0.0 ! levels.
      ENDIF
C
C  Horizontal grid description -----------------------------------------
C  Position of pole (from dump fixed-length header)
C  Grid orientation (hardwired 0.0)
C  Origin and spacing of grid (depends on output grid type)
C
      PP_REAL_HEAD(BPLAT)=REALHD(5)
      PP_REAL_HEAD(BPLON)=REALHD(6)
      PP_REAL_HEAD(BGOR)=0.0
      IF (samples.GT.0) THEN   ! Indicates a timeseries/trajectory
        PP_REAL_HEAD(BZX)=0.0
        PP_REAL_HEAD(BDX)=0.0
        PP_REAL_HEAD(BZY)=0.0
        PP_REAL_HEAD(BDY)=0.0
      ELSE
        IF (OCEAN) THEN       !   set BZY,BZX,BDY,BDX for ocean
          IF (st_grid.EQ.st_uv_grid .OR. st_grid.EQ.st_zu_grid
     &        .OR. st_grid.EQ.st_mu_grid) THEN
            PP_REAL_HEAD(BZY)=REALHD(3)-REALHD(2)/2.0
            PP_REAL_HEAD(BZX)=REALHD(4)-REALHD(1)/2.0
          ELSEIF (st_grid.EQ.st_tp_grid .OR. st_grid.EQ.st_zt_grid
     &       .OR. st_grid.EQ.st_mt_grid .OR. st_grid.EQ.st_scalar) THEN
            PP_REAL_HEAD(BZY)=REALHD(3)-REALHD(2)
            PP_REAL_HEAD(BZX)=REALHD(4)-REALHD(1)
          ELSEIF (st_grid.EQ.st_cu_grid) THEN
            PP_REAL_HEAD(BZY)=REALHD(3)-REALHD(2)
            PP_REAL_HEAD(BZX)=REALHD(4)-REALHD(1)/2.0
          ELSEIF (st_grid.EQ.st_cv_grid) THEN
            PP_REAL_HEAD(BZY)=REALHD(3)-REALHD(2)/2.0
            PP_REAL_HEAD(BZX)=REALHD(4)-REALHD(1)
          ENDIF
          IF (REALHD(32).GT.REALHD(29)) THEN !   greater than RMDI
            PP_REAL_HEAD(BDY)=0.0
            PP_REAL_HEAD(BDX)=REALHD(32)
          ELSE
            PP_REAL_HEAD(BDY)=REALHD(2)
            PP_REAL_HEAD(BDX)=REALHD(1)
          ENDIF
        ELSE                 !   set BZY,BZX,BDY,BDX for atmos
          IF(st_grid.EQ.st_uv_grid.OR.st_grid.EQ.st_cv_grid .OR.
     &       st_grid.EQ.st_zu_grid.OR.st_grid.EQ.st_mu_grid) THEN
            PP_REAL_HEAD(BZY)=REALHD(3)+REALHD(2)/2.0 ! UV pts
          ELSE
            PP_REAL_HEAD(BZY)=REALHD(3)+REALHD(2) ! Zeroth Lat BZY
          ENDIF
C
          IF(st_grid.EQ.st_uv_grid.OR.st_grid.EQ.st_cu_grid .OR.
     &       st_grid.EQ.st_zu_grid.OR.st_grid.EQ.st_mu_grid) THEN
            PP_REAL_HEAD(BZX)=REALHD(4)-REALHD(1)/2.0 !UV points
          ELSE
            PP_REAL_HEAD(BZX)=REALHD(4)-REALHD(1) ! Zeroth Long BZX
          ENDIF
          PP_REAL_HEAD(BDX)=REALHD(1) ! Long intvl BDX
          PP_REAL_HEAD(BDY)=-REALHD(2) ! Lat intvl BDY
        ENDIF
C
C Add on offset for fields not starting from the origin
C
        PP_REAL_HEAD(BZY)=PP_REAL_HEAD(BZY)
     &                    +(NROW_IN-1)*PP_REAL_HEAD(BDY)
        PP_REAL_HEAD(BZX)=PP_REAL_HEAD(BZX)
     &                    +(WCOL_IN-1)*PP_REAL_HEAD(BDX)
        IF(PP_REAL_HEAD(BZX).GE.360.0)
     *     PP_REAL_HEAD(BZX)=PP_REAL_HEAD(BZX)-360.0
C
C If horizontal averaging has been applied to the output field,
C set BDX and/or BDY to the full (sub)domain extent which was processed.
C If the input field was intrinsically non-2D (eg. zonal), assume that
C the collapsed dimension(s) covered the full model domain.
C
        mean_code=(GR/block_size)*block_size
        IF (st_grid.EQ.st_zt_grid .OR. st_grid.EQ.st_zu_grid
     &      .OR. st_grid.EQ.st_scalar) THEN
          PP_REAL_HEAD(BDX)=REAL(INTHD(6))*PP_REAL_HEAD(BDX)
        ELSEIF (mean_code.EQ.zonal_mean_base .OR.
     &      mean_code.EQ.field_mean_base .OR.
     &      mean_code.EQ.global_mean_base) THEN
          PP_REAL_HEAD(BDX)=ABS(REAL(ECOL_IN-WCOL_IN))*PP_REAL_HEAD(BDX)
        ENDIF
C
        IF (st_grid.EQ.st_mt_grid .OR. st_grid.EQ.st_mu_grid
     &      .OR. st_grid.EQ.st_scalar) THEN
          PP_REAL_HEAD(BDY)=REAL(INTHD(7))*PP_REAL_HEAD(BDY)
        ELSEIF (mean_code.EQ.merid_mean_base .OR.
     &      mean_code.EQ.field_mean_base .OR.
     &      mean_code.EQ.global_mean_base) THEN
          PP_REAL_HEAD(BDY)=ABS(REAL(NROW_IN-SROW_IN))*PP_REAL_HEAD(BDY)

        ENDIF
      ENDIF
C
C Missing data indicator (from PARAMETER) ------------------------------
C MKS scaling factor (unity as model uses SI units throughout)
C
      PP_REAL_HEAD(BMDI)=RMDI
      PP_REAL_HEAD(BMKS)=1.0
C
  999 CONTINUE
      RETURN
      END
