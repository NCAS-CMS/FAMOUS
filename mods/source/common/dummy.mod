*ID DUMMY
*/
*/  Dummy routines to remove unresolved externals linking errors
*/
*DECK DUMMY1

*IF DEF,A18_1A,OR,DEF,A18_2A,OR,DEF,O35_1A
C
*ELSE
      subroutine buffin_shmem()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'buffin_shmem'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'buffin_shmem'
      call abort()

      stop
      end

      subroutine buffin_acobs()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'buffin_acobs'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'buffin_acobs'
      call abort()

      stop
      end

      subroutine swapbounds_shmem()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'swapbounds_shmem'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'swapbounds_shmem'
      call abort()

      stop
      end

      subroutine swapbounds_sum()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'swapbounds_sum'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'swapbounds_sum'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A18_1A,OR,DEF,A18_2A,OR,DEF,RECON
C
*ELSE
      subroutine stratq()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'stratq'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'stratq'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A18_1A,OR,DEF,A18_2A
C
*ELSE
      subroutine ac()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'ac'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'ac'
      call abort()

      stop
      end

      subroutine var_umprocessing()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'var_umprocessing'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'var_umprocessing'
      call abort()

      stop
      end

      subroutine var_umsetup()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'var_umsetup'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'var_umsetup'
      call abort()

      stop
      end

      subroutine ac_init()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'ac_init'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'ac_init'
      call abort()

      stop
      end

      subroutine iau_ctl()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'iau_ctl'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'iau_ctl'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A11_1A
C
*ELSE
      subroutine set_trac()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'set_trac'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'set_trac'
      call abort()

      stop
      end

      subroutine trac_vert_adv()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'trac_vert_adv'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'trac_vert_adv'
      call abort()

      stop
      end

      subroutine trac_adv()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'trac_adv'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'trac_adv'
      call abort()

      stop
      end

      subroutine trbdry()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'trbdry'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'trbdry'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A11_1A,OR,DEF,A03_7A,OR,DEF,A03_6A
C
*ELSE
      subroutine trsrce()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'trsrce'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'trsrce'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A01_1A,OR,DEF,A01_1B,OR,DEF,A01_2A,OR,DEF,A01_2B
C
*ELSE
      subroutine swlkin()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'swlkin'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'swlkin'
      call abort()

      stop
      end

      subroutine swdkdi()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'swdkdi'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'swdkdi'
      call abort()

      stop
      end

      subroutine swrad()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'swrad'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'swrad'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A02_1A,OR,DEF,A02_1B,OR,DEF,A02_1C
C
*ELSE
      subroutine lwlkin()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'lwlkin'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'lwlkin'
      call abort()

      stop
      end

      subroutine lwrad()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'lwrad'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'lwrad'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A70_1A,OR,DEF,A70_1B
C
*IF DEF,A01_3A,OR,DEF,A02_3A
C
*ELSE
      subroutine r2_sw_specin()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'r2_sw_specin'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'r2_sw_specin'
      call abort()

      stop
      end

      subroutine r2_lw_specin()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'r2_lw_specin'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'r2_lw_specin'
      call abort()

      stop
      end

      subroutine tropin()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'tropin'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'tropin'
      call abort()

      stop
      end

      subroutine r2_global_cloud_top()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'r2_global_cloud_top'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'r2_global_cloud_top'
      call abort()

      stop
      end
*ENDIF
*ELSE
      subroutine r2_sw_specin()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'r2_sw_specin'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'r2_sw_specin'
      call abort()

      stop
      end

      subroutine r2_lw_specin()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'r2_lw_specin'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'r2_lw_specin'
      call abort()

      stop
      end

      subroutine tropin()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'tropin'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'tropin'
      call abort()

      stop
      end

      subroutine r2_global_cloud_top()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'r2_global_cloud_top'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'r2_global_cloud_top'
      call abort()

      stop
      end

      subroutine gas_calc()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'gas_calc'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'gas_calc'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A01_3A
C
*ELSE
      subroutine r2_swrad()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'r2_swrad'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'r2_swrad'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A02_3A
C
*ELSE
      subroutine r2_lwrad()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'r2_lwrad'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'r2_lwrad'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A09_2A,OR,DEF,A09_2B
C
*ELSE
      subroutine area_cld()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'area_cld'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'area_cld'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A07_1A,OR,DEF,A07_1B
C
*ELSE
      subroutine vdif_ctl()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'vdif_ctl'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'vdif_ctl'
      call abort()

      stop
      end

      subroutine vert_dif()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'vert_dif'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'vert_dif'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A09_2B
C
*ELSE
      subroutine rhcrit_calc()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'rhcrit_calc'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'rhcrit_calc'
      call abort()

      stop
      end
*ENDIF

*IF DEF,CONTROL,AND,DEF,WAVE
C
*ELSE
      subroutine stwvgt()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'stwvgt'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'stwvgt'
      call abort()

      stop
      end
*ENDIF

*IF DEF,CONTROL,AND,DEF,OCEAN
C
*ELSE
      subroutine stocgt()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'stocgt'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'stocgt'
      call abort()

      stop
      end
*ENDIF

*IF DEF,CONTROL,AND,DEF,OCNASSM
C
*ELSE
      subroutine oc_ac_ctl()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'oc_ac_ctl'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'oc_ac_ctl'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A03_7A
C
*ELSE
      subroutine rad_moses()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'rad_moses'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'rad_moses'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A17_1A
C
*ELSE
      subroutine new2old()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'new2old'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'new2old'
      call abort()

      stop
      end

      subroutine sootscav()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'sootscav'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'sootscav'
      call abort()

      stop
      end

      subroutine sulphur()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'sulphur'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'sulphur'
      call abort()

      stop
      end

      subroutine gravsett()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'gravsett'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'gravsett'
      call abort()

      stop
      end
*ENDIF

*IF DEF,O35_1A
C
*ELSE
      subroutine oa_zero()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'oa_zero'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'oa_zero'
      call abort()

      stop
      end

      subroutine oa_int_lev()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'oa_int_lev'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'oa_int_lev'
      call abort()

      stop
      end

      subroutine oa_int_1d()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'oa_int_1d'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'oa_int_1d'
      call abort()

      stop
      end
*ENDIF

*IF DEF,C90_1A,OR,DEF,C90_2A,OR,DEF,C90_2B
C
*ELSE
      subroutine p_to_cv()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'p_to_cv'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'p_to_cv'
      call abort()

      stop
      end

      subroutine p_to_cu()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'p_to_cu'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'p_to_cu'
      call abort()

      stop
      end

      subroutine p_to_uv()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'p_to_uv'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'p_to_uv'
      call abort()

      stop
      end
*ENDIF

*IF DEF,SEAICE
C
*ELSE
      subroutine cnvstop()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'cnvstop'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'cnvstop'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A87_1A
C
*ELSE
      subroutine zonm_atm()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'zonm_atm'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'zonm_atm'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A14_1A,OR,DEF,A14_1B
C
*ELSE
      subroutine init_emcorr()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'init_emcorr'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'init_emcorr'
      call abort()

      stop
      end

      subroutine add_eng_corr()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'add_eng_corr'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'add_eng_corr'
      call abort()

      stop
      end

      subroutine eng_mass_diag()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'eng_mass_diag'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'eng_mass_diag'
      call abort()

      stop
      end

      subroutine cal_eng_mass_corr()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'cal_eng_mass_corr'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'cal_eng_mass_corr'
      call abort()

      stop
      end

      subroutine flux_diag()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'flux_diag'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'flux_diag'
      call abort()

      stop
      end
*ENDIF

*IF DEF,A06_1A,OR,DEF,A06_2A,OR,DEF,A06_3A,OR,DEF,A06_3B
C
*ELSE
      subroutine gwav_intctl()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'gwav_intctl'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'gwav_intctl'
      call abort()

      stop
      end

*ENDIF

*IF DEF,T3E,OR,DEF,CRAY
C
*ELSE
      subroutine ibm2cri()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'ibm2cri'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'ibm2cri'
      call abort()

      stop
      end

      subroutine cri2ibm()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'cri2ibm'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'cri2ibm'
      call abort()

      stop
      end

      subroutine strmov()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'strmov'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'strmov'
      call abort()

      stop
      end

      subroutine movbit()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'movbit'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'movbit'
      call abort()

      stop
      end
*ENDIF

      subroutine buffout_shmem()

      write(0,*) 'Error you have called an undefined subroutine ',
     :           'buffout_shmem'
      write(6,*) 'Error you have called an undefined subroutine ',
     :           'buffout_shmem'
      call abort()

      stop
      end
