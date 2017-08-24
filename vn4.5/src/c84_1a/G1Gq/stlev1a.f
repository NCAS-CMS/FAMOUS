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
CLL  Subroutine STLEVELS -----------------------------------------------
CLL
CLL  Purpose: Generate a level index from STASHrecord and level_lists
CLL           and number of levels tailored to a particular diagnostic.
CLL           Also set levels and pseudo-levels information for encoding
CLL           PPheader details.  (This subroutine based on a merger
CLL           between GEN_INDEX and PP_COMPUTE_LEVEL).
CLL                  New subroutine STLEVELS is based on GEN_INDEX and
CLL                  PP_COMPUTE_LEVEL with merged functionality.
CLL           A general note as levels list is an integer
CLL           real values are multiplied by a 1000.0.
CLL           When computing the real value of the level for the
CLL           pp header it is necessary to divide by a 1000.0.
CLL           Levels that are affected by this are theta, pressure and
CLL           height. S. Anderson.
CLL
CLL  Author:   T.Johns
CLL
CLL  Tested under compiler:   cft77
CLL  Tested under OS version: UNICOS 5.1
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.1  14/01/93  Set pseudo_level to 0 not IMDI if no pseudo-level.
CLL        29/01/93  Correct FORMAT statements.
CLL  3.1     14/01/93 Include PV levels as levels divided by 1000.
CLL   3.2  19/04/93  Correct roundoff error for LEVEL type conversion.
CLL                  1.0E-10 is added after REAL divide by 1000 (TCJ).
CLL  4.0  14/12/95  Correct long-standing error in input levels range 
CLL                 to output levels list conversion.  RTHBarnes.
!    4.4  02/12/96 Time mean timeseries added R A Stratton.             
CLL
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered : C4?
CLL
CLL  Project task: C4
CLL
CLL  External documentation : UMDP no C4
CLL
C*L  Interface and arguments: ------------------------------------------
C
      SUBROUTINE STLEVELS(stash_control,stash_control_size,
     +     stash_levels,num_stash_levels,num_level_lists,
     +     stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists,
     +     max_stash_levs,num_levs_in,num_levs_out,index_size,
     +     index_lev,level_list,
     +     lbvcl,ak,bk,level,pseudo_level,ak_lev,bk_lev,
     +     icode,cmessage)
C
      IMPLICIT NONE
C
      INTEGER
     &       stash_control_size, ! IN size of stash control record
     &       stash_control(stash_control_size),! IN  stash control
     &       num_stash_levels,   ! IN max. no of hts for a levels list
     &       num_level_lists,    ! IN max. no of level lists
     &       stash_levels(num_stash_levels+1,num_level_lists), ! IN
C                                !    lookup table for level lists
     &       num_stash_pseudo,num_pseudo_lists,! IN dims of pseudo_levs
     &       stash_pseudo_levels(num_stash_pseudo+1,num_pseudo_lists),
C                                ! IN lookup table for pseudo-lev lists
     &       max_stash_levs,     ! IN max. no of output levels
     &       num_levs_in,        ! OUT no of levels in input data
     &       num_levs_out,       ! OUT no of levels in output data
     &       index_size,         ! OUT no of levels in levels index
     &       index_lev(max_stash_levs), ! OUT index of output level
C                                               relative to input level
     &       level_list(max_stash_levs), ! OUT value of model level
     &       pseudo_level(max_stash_levs), ! OUT Value of pseudo levels
     &       lbvcl,              ! IN  vertical coordinate PP code
     &       icode               ! OUT error code
      REAL
     &       ak(*),                 ! IN  Hybrid Ak value on model levs
     &       bk(*),                 ! IN  Hybrid Bk value on model levs
     &       level(max_stash_levs), ! OUT Value of output levels (real)
     &       ak_lev(max_stash_levs),! OUT Hybrid Ak value on output levs
     &       bk_lev(max_stash_levs) ! OUT Hybrid Bk value on output levs
      CHARACTER*(*)
     &       cmessage            ! OUT error message
C*----------------------------------------------------------------------
C Parameters
C
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
C
C Local variables
C
      INTEGER
     &       index_pseudo_lev(max_stash_levs), ! Pseudo-level 1D index
     &       num_pseudo_in,num_pseudo_out,     ! Number of pseudo levs
     &       k2,ml,kl,                 ! loop counts
     &       NI,NO,                    ! Number In/Out
     &       indx1,                    ! index count
     &       ilev,                     ! Integer level/pseudo-level
     &       what_mean,what_proc       ! Meaning and processing code
C
C First compute the index for physical levels
C
      IF(STASH_CONTROL(st_input_bottom).LT.0) THEN ! Input LEVELS list
        NI=-STASH_CONTROL(st_input_bottom)
        NUM_LEVS_IN=STASH_LEVELS(1,NI)
        IF(STASH_CONTROL(st_output_bottom).LT.0) THEN ! LEVELS LIST out
          NO=-STASH_CONTROL(st_output_bottom)
          NUM_LEVS_OUT=STASH_LEVELS(1,NO)
          INDX1=0
          DO ML=1,NUM_LEVS_OUT
            ilev=STASH_LEVELS(ML+1,NO)    !  Level required
            DO KL=1,NUM_LEVS_IN
              IF(STASH_LEVELS(KL+1,NI).EQ.ilev) THEN
                INDX1=INDX1+1
                INDEX_LEV(INDX1)=KL   ! Relative position of Input to Ou
                level_list(indx1)=ilev
                GOTO 400
              ENDIF
            ENDDO
            ICODE=nonsense
            WRITE(CMESSAGE,101) 'Output level ',ilev,
     &                          ' not found in input levels list'
            GOTO 999
 400        CONTINUE
            ENDDO
        ELSE           !  Output as a Level range
          NUM_LEVS_OUT=STASH_CONTROL(st_output_top)-
     &                 STASH_CONTROL(st_output_bottom)+1
          ilev=STASH_CONTROL(st_output_bottom) !1st output model level
          DO KL=1,NUM_LEVS_IN
            IF(STASH_LEVELS(KL+1,NI).EQ.ilev) THEN
              INDEX_LEV(1)=KL ! Relative posn of Input to the 1st level
              level_list(1)=ilev
              GOTO 401
            ENDIF
          ENDDO
          ICODE=nonsense
          WRITE(CMESSAGE,101) 'Output bottom model level ',ilev,
     &                        ' not found in input levels list'
          GOTO 999
 401      CONTINUE
          DO KL=2,NUM_LEVS_OUT
            INDEX_LEV(KL)=INDEX_LEV(KL-1)+1
            level_list(kl)=level_list(kl-1)+1
          ENDDO
        ENDIF
      ELSEIF(STASH_CONTROL(st_input_bottom).EQ.100) THEN !Special level
          NUM_LEVS_IN=1
          NUM_LEVS_OUT=1
          INDEX_LEV(1)=1
          level_list(1)=1 ! could be worth setting to some nonsense no.
      ELSE     !  Input is Model level range
        NUM_LEVS_IN=STASH_CONTROL(st_input_top)-
     &              STASH_CONTROL(st_input_bottom)+1
        IF(STASH_CONTROL(st_output_bottom).LT.0) THEN ! LEVELS LIST out
          NO=-STASH_CONTROL(st_output_bottom)
          NUM_LEVS_OUT=STASH_LEVELS(1,NO)
          INDX1=0
          DO ML=1,NUM_LEVS_OUT
            ilev=STASH_LEVELS(ML+1,NO)    ! Output level reqd
            DO KL=1,NUM_LEVS_IN
              IF((STASH_CONTROL(st_input_bottom)+KL-1).EQ.ilev) THEN
                INDX1=INDX1+1
                INDEX_LEV(INDX1)=KL   ! Relative posn of output to inpt
                level_list(INDX1)=ilev
                GOTO 402
              ENDIF
            ENDDO
            ICODE=nonsense
            WRITE(CMESSAGE,101) 'Output model level ',ilev,
     &                          ' not in input model level range'
            GOTO 999
 402        CONTINUE
          ENDDO
        ELSE     !   Output as model level range
C Do some consistency checks here to ensure valid processing request
C output bottom should be greater or equal to input bottom
          IF (stash_control(st_output_bottom).lt.
     +       stash_control(st_input_bottom)) THEN
            icode=nonsense
            write(cmessage,103)'bad level spec, bot input>output',
     +       stash_control(st_input_bottom),
     +       stash_control(st_output_bottom)
            goto 999 ! jump to error
          ELSEIF (stash_control(st_output_top).gt.
     +         stash_control(st_input_top)) THEN
            icode=nonsense
            write(cmessage,103)'bad level spec, top input<output',
     +        stash_control(st_input_top),
     +        stash_control(st_output_top)
              goto 999 ! jump to error
          ENDIF
          NUM_LEVS_OUT=STASH_CONTROL(st_output_top)-
     &                 STASH_CONTROL(st_output_bottom)+1
          INDEX_LEV(1)=STASH_CONTROL(st_output_bottom)-
     &                 STASH_CONTROL(st_input_bottom)+1
          level_list(1)=stash_control(st_output_bottom)
          DO kl=2,NUM_LEVS_OUT
            INDEX_LEV(kl)=INDEX_LEV(kl-1)+1
            level_list(kl)=level_list(kl-1)+1
          ENDDO
        ENDIF
      ENDIF
      index_size=num_levs_out
      IF (num_levs_out.gt.num_levs_in) THEN   ! things very badly wrong
        icode=nonsense
        write(cmessage,103)'asking for num_levs_out>num_levs_in',
     +   num_levs_out,num_levs_in
        goto 999 ! jump to return
      ENDIF
C
C Next, compute actual (physical) levels for encoding PPheaders
C
      IF (STASH_CONTROL(st_output_bottom).LT.0) THEN ! Levels List ?
        NO=-STASH_CONTROL(st_output_bottom)     ! Index of Levels list
C Pressure or Height or Theta levels or PV levels?
        IF(LBVCL.EQ.8.OR.LBVCL.EQ.1.OR.LBVCL.EQ.19
     &                             .OR.LBVCL.EQ.82) THEN

          DO ML=1,NUM_LEVS_OUT
            LEVEL(ML)=REAL(STASH_LEVELS(ML+1,NO))*0.001+1.0E-10
          ENDDO
        ELSE
          DO ML=1,NUM_LEVS_OUT
            LEVEL(ML)=REAL(STASH_LEVELS(ML+1,NO))
          ENDDO
        ENDIF
      ELSEIF (STASH_CONTROL(st_output_bottom).EQ.st_special_code) THEN
C Special level
        DO ML=1,NUM_LEVS_OUT
          LEVEL(ML)=-1.0
        ENDDO
      ELSE
        DO ML=1,NUM_LEVS_OUT
          LEVEL(ML)=REAL(STASH_CONTROL(st_output_bottom)+ML-1)
        ENDDO
      ENDIF
C
      IF (lbvcl.eq.9) THEN
        DO ML=1,NUM_LEVS_OUT
          ilev=INT(LEVEL(ML))
          ak_lev(ML)=ak(ilev)
          bk_lev(ML)=bk(ilev)
        ENDDO
      ELSE
        DO ML=1,NUM_LEVS_OUT
          ak_lev(ML)=0.0
          bk_lev(ML)=0.0
        ENDDO
      ENDIF
C
C Now reset the number of output levels to 1 if vertical compression is
C to be done in SPATIAL.  NB: index_lev and level_list need to be filled
C with values corresponding to the full range of levels processed.
C
      what_proc=STASH_CONTROL(st_proc_no_code)
      what_mean=(STASH_CONTROL(st_gridpoint_code)/block_size)*block_size
      IF(what_mean.EQ.vert_mean_base .OR. what_mean.EQ.global_mean_base
     &   .OR. what_proc.EQ.st_time_series_code
     &   .OR. what_proc.EQ.8                                            
     &   .OR. what_proc.EQ.st_append_traj_code) num_levs_out=1
C
C Next compute the index for pseudo levels, if there are any
C
      IF(STASH_CONTROL(st_pseudo_in).GT.0) THEN ! Input PSEUDO_LEVELS
        NI=STASH_CONTROL(st_pseudo_in)
        num_pseudo_in=STASH_PSEUDO_LEVELS(1,NI)
        IF(STASH_CONTROL(st_pseudo_out).GT.0) THEN ! Output PSEUDO_LEVS
          NO=STASH_CONTROL(st_pseudo_out)
          num_pseudo_out=STASH_PSEUDO_LEVELS(1,NO)
          INDX1=0
          DO ML=1,NUM_PSEUDO_OUT
            ilev=STASH_PSEUDO_LEVELS(ML+1,NO)   !  Level required
            DO KL=1,NUM_PSEUDO_IN
              IF(STASH_PSEUDO_LEVELS(KL+1,NI).EQ.ilev) THEN
                INDX1=INDX1+1
                INDEX_PSEUDO_LEV(INDX1)=KL
                pseudo_level(indx1)=ilev
                GOTO 500
              ENDIF
            ENDDO
            ICODE=nonsense
            WRITE(CMESSAGE,101) 'Output pseudo level ',ilev,
     &                          ' not found in input levels list'
            GOTO 999
 500        CONTINUE
          ENDDO
        ELSE  ! Illegal combination
          ICODE=nonsense
          WRITE(CMESSAGE,101) 'Input pseudo level list ',NI,
     &         ' has illegal output pseudo levels list'
          GOTO 999
        ENDIF
      ELSE  ! Only levels lists are supported for pseudo levels
        num_pseudo_out=0
      ENDIF
C
C Next expand the separate indexes and physical levels arrays into
C combined arrays if necessary, taking care not to overwrite earlier
C parts of the arrays.  If no pseudo-levels, set pseudo-level to 0.
C
      IF (num_pseudo_out.GT.0) THEN
        DO K2=num_pseudo_out,1,-1
          DO ML=1,num_levs_out
            INDEX_LEV(ML+(K2-1)*num_levs_out)=
     &        (INDEX_PSEUDO_LEV(K2)-1)*num_levs_in+INDEX_LEV(ML)
            level(ML+(K2-1)*num_levs_out)=level(ML)
            ak_lev(ML+(K2-1)*num_levs_out)=ak_lev(ML)
            bk_lev(ML+(K2-1)*num_levs_out)=bk_lev(ML)
          ENDDO
          DO ML=num_levs_out,1,-1
            pseudo_level(ML+(K2-1)*num_levs_out)=pseudo_level(K2)
          ENDDO
        ENDDO
        num_levs_out=num_levs_out*num_pseudo_out
      ELSE
        DO ML=1,num_levs_out
          pseudo_level(ML)=0
        ENDDO
      ENDIF
C
999   CONTINUE ! jump here for error return
 101  FORMAT('STLEVELS : ',a,i6,a)
 103  FORMAT('STLEVELS : >> FATAL ERROR <<',a,2i5)
      RETURN
      END

