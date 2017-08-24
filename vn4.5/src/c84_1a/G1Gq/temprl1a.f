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
CLL  Routine: TEMPORAL -------------------------------------------------
CLL
CLL  Purpose: Control routine to handle temporal processing options
CLL           within STASH.  Its input and output arguments look like
CLL           1D arrays (ie. all the data should be in contiguous areas
CLL           of memory).  Lower level service routines are called to
CLL           perform the individual processing options.
CLL
CLL  Tested under compiler:   cft77
CLL  Tested under OS version: UNICOS 5.1
CLL
CLL  Author:   S.Tett
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL   3.1  24/02/93  Change name of variable 'end' to 'last_ts' (ST).
!     4.4  25/11/96  Add processing code option 8 - daily mean          
!                    timeseries. R A Stratton.                          
CLL
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered: D72
CLL
CLL  Project task: D7
CLL
CLL  External documentation:
CLL    Unified Model Doc Paper C4 - Storage handling and diagnostic
CLL                                 system (STASH)
CLL
C*L  Interface and arguments: ------------------------------------------
C
      SUBROUTINE TEMPORAL(variable,result,size,extra_size,
     &  control,control_size,ocean,
     +  timestep,error,errmssg,start,amdi)
C
      IMPLICIT NONE
C
      INTEGER size                  ! IN  size of arrays
      REAL variable(size)           ! IN  data array
      REAL result(size)             ! OUT output array
      INTEGER extra_size            ! IN size of extra data
      INTEGER control_size          ! IN  size of control
      INTEGER control(control_size) ! IN  control
      INTEGER timestep              ! IN  present value of timestep
      INTEGER error                 ! OUT error code
      CHARACTER*(*) errmssg         ! OUT error message
      REAL amdi                     ! IN  missing data indicator
      LOGICAL ocean                 ! IN  true if ocean diagnostic
      LOGICAL start                 ! OUT true if start timestep
C*----------------------------------------------------------------------
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
C Subroutines called
C
      EXTERNAL staccum,stmax,stmin
C
C Local variables
C
      LOGICAL masking        ! indicator for masking (ie. missing data)
      INTEGER proc_code      ! value of processing code
      INTEGER mask_code      ! value of masking code
      REAL divisor           ! divisor for the time mean (1/period)
      INTEGER mod_period     ! timesteps since start modulo period.
      INTEGER start_time     ! value of start time
      INTEGER i              ! loop counter
      LOGICAL last_ts        ! true if end timestep
      INTEGER proc_size      ! size of data to be processed
CL---------------------------------------------------------------------
CL 1. Set processing option code and select appropriate service routine
CL
      proc_size=size-extra_size
      proc_code=control(st_proc_no_code)
C
C  Replace (null processing)
C
      IF (proc_code.eq.st_replace_code) THEN
        DO i=1,size
          result(i)=variable(i)
        ENDDO
        start=(control(st_start_time_code).eq.timestep)
C
C  Mean/accumulation
C
      ELSEIF (proc_code.eq.st_accum_code.or.
     +        proc_code.eq.st_time_mean_code) THEN
        start_time=control(st_start_time_code)
        IF (control(st_period_code).EQ.st_infinite_time) THEN
          start=(timestep.eq.start_time)
          last_ts=.FALSE.
        ELSE
          mod_period=mod(timestep-start_time,control(st_period_code))
          start=(mod_period.eq.0)
          last_ts=(mod_period.eq.(control(st_period_code)-
     &                        control(st_freq_code)))
        ENDIF
        mask_code=control(st_gridpoint_code)
        mask_code=mod(mask_code,block_size)
        masking=(mask_code.ne.stash_null_mask_code).or.ocean
        IF (start) THEN      ! first timestep.
          DO i=1,size
            result(i)=variable(i)
          ENDDO
        ELSE
          CALL STACCUM(variable,result,proc_size,masking,amdi)
          DO i=proc_size+1,size
            result(i)=variable(i) ! copy over the extra data (if any)
          ENDDO
        ENDIF
C  Normalise at end of mean period
        IF (last_ts.and.proc_code.eq.st_time_mean_code) THEN
          divisor=(float(control(st_freq_code))/
     &             float(control(st_period_code)))
C If field is masked test for MDI, otherwise don't
          IF (masking) THEN
            DO i=1,proc_size
              IF (result(i).ne.amdi) THEN
                result(i)=result(i)*divisor
              ENDIF
            ENDDO
          ELSE
            DO i=1,proc_size
              result(i)=result(i)*divisor
            ENDDO
          ENDIF
        ENDIF
C
C  Maximum
C
      ELSEIF (proc_code.eq.st_max_code) THEN
        start_time=control(st_start_time_code)
        mod_period=mod(timestep-start_time,control(st_period_code))
        start=(mod_period.eq.0)
        IF (start) THEN
          DO i=1,size
            result(i)=variable(i)
          ENDDO
        ELSE
          mask_code=control(st_gridpoint_code)
          mask_code=mod(mask_code,block_size)
          masking=(mask_code.ne.stash_null_mask_code).or.ocean
          CALL STMAX(variable,result,proc_size,masking,amdi)
          DO i=proc_size+1,size
            result(i)=variable(i) ! copy over the extra data (if any)
          ENDDO
        ENDIF
C
C  Minimum
C
      ELSEIF (proc_code.eq.st_min_code) THEN
        start_time=control(st_start_time_code)
        mod_period=mod(timestep-start_time,control(st_period_code))
        start=(mod_period.eq.0)
        IF (start) THEN
          DO i=1,size
            result(i)=variable(i)
          ENDDO
        ELSE
          mask_code=control(st_gridpoint_code)
          mask_code=mod(mask_code,block_size)
          masking=(mask_code.ne.stash_null_mask_code).or.ocean
          CALL STMIN(variable,result,proc_size,masking,amdi)
          DO i=proc_size+1,size
            result(i)=variable(i) ! copy over the extra data (if any)
          ENDDO
        ENDIF
C
C  Timeseries (append)
C
      ELSEIF (proc_code.eq.st_time_series_code) THEN
        DO i=1,size
C Note that on start timestep this will include the extra data
          result(i)=variable(i)
        ENDDO
        start_time=control(st_start_time_code)
        mod_period=mod(timestep-start_time,control(st_period_code))
        start=(mod_period.eq.0)
        last_ts=(mod_period.eq.(control(st_period_code)-
     &                      control(st_freq_code)))
C
C  Append trajectories
C
      ELSEIF (proc_code.eq.st_append_traj_code) THEN
        start_time=control(st_start_time_code)
        mod_period=mod(timestep-start_time,control(st_period_code))
        start=(mod_period.eq.0)
        last_ts=(mod_period.eq.(control(st_period_code)-
     &                      control(st_freq_code)))
        error=st_not_supported
        write(errmssg,100)' do not support append trajects'
        goto 999
!                                                                       
!  Timeseries (append) - option 8 daily mean
!                                                                       
      ELSEIF (proc_code.eq.st_time_series_mean) THEN      
                                                                        
        DO i=1,size                                                     
C Note that on start timestep this will include the extra data          
          result(i)=variable(i)                                         
        ENDDO                                                           
        start_time=control(st_start_time_code)                          
        mod_period=mod(timestep-start_time,control(st_period_code))     
        start=(mod_period.eq.0)                                         
        last_ts=(mod_period.eq.(control(st_period_code)-                
     &                      control(st_freq_code)))                     
C
C  Error condition
C
      ELSE
        error=unknown_processing
        write(errmssg,101)' unknown processing code',proc_code
        goto 999
      ENDIF
C
999   CONTINUE   ! jump for errors
C
100   FORMAT('TEMPORAL : >>> FATAL ERROR <<<',a30)
101   FORMAT('TEMPORAL : >>> FATAL ERROR <<<',a30,i5)
C
      RETURN
      END
