*/=============================================
*/Mod to enable a perturbed physics ensemble
*/of ocean carbon cycle parameters to
*/be performed. It is based on mods from 
*/Ben Booth but has been adapted to take 
*/account of the necessary different coding 
*/standard needed for FAMOUS UM version 4.5.3  
*/versus HadCM3 Um version 4.7. Comments
*/within this mod refer to the line below the
*/comment.
*/
*/Jonny Williams, May 2011
*/
*/=============================================
*ID ID_OCEAN_CC
*/
*DECLARE OBIOCONST
*D OBIOCONST.25,91
      COMMON /OCEAN_CC/ grow_sat,mort_sat,Q10H,ref_t
     &    ,c2chl,c2n_p,c2n_z,c2n_d
     &    ,alpha,alphmx,psmax,resp_rate,phyto_min,pmort_max
     &    ,z_mort_1,z_mort_2,z_remin,z_detrit
     &    ,remin_rate_shallow,remin_rate_deep
     &    ,graze_max,beta_p,beta_dt,sink_rate_dt
     &    ,par,c_mol_wt,chl2pig
     &    ,holling_coef,graze_sat,graze_threshold
     &    ,n2be_p,n2be_z,n2be_d,mw_nitrogen,mw_carbon
     &    ,c2n_redfield
     &    ,rain_ratio
     &    ,lysocline
*/
*DECLARE READNLST
*I GRB1F305.511
*CALL OBIOCONST
*I NT071293.25
       NAMELIST /OCEAN_CC/ grow_sat,mort_sat,Q10H,ref_t
     &    ,c2chl,c2n_p,c2n_z,c2n_d
     &    ,alpha,alphmx,psmax,resp_rate,phyto_min,pmort_max
     &    ,z_mort_1,z_mort_2,z_remin,z_detrit
     &    ,remin_rate_shallow,remin_rate_deep
     &    ,graze_max,beta_p,beta_dt,sink_rate_dt
     &    ,par,c_mol_wt,chl2pig
     &    ,holling_coef,graze_sat,graze_threshold
     &    ,n2be_p,n2be_z,n2be_d,mw_nitrogen,mw_carbon
     &    ,c2n_redfield
     &    ,rain_ratio
     &    ,lysocline
*I ORH1F305.227
      grow_sat = 0.1
      mort_sat = 0.1
      Q10H = 1.0
      ref_t = 10.0
      c2chl = 40.0
      c2n_p = 6.625
      c2n_z = 5.625
      c2n_d = 7.500
      c2n_redfield = 6.625
      alpha=0.02
      alphmx=alpha*2.602
      psmax = 0.6
      resp_rate = 0.02
      phyto_min = 0.01
      pmort_max = 0.05
      z_mort_1 = 0.02
      z_mort_2 = 0.3
      z_remin = 0.667
      z_detrit=1.0 - z_remin
      remin_rate_shallow = 0.1
      remin_rate_deep = 0.02
      graze_sat = 0.75
      graze_threshold = 0.1
      holling_coef = 2.0
      graze_max = 1.0
      mw_nitrogen = 14.01
      mw_carbon = 12.01
      n2be_p = (mw_nitrogen+mw_carbon*c2n_p)
     &                   /(mw_nitrogen+mw_carbon*c2n_redfield)
      n2be_z = (mw_nitrogen+mw_carbon*c2n_z)
     &                   /(mw_nitrogen+mw_carbon*c2n_redfield)
      n2be_d = (mw_nitrogen+mw_carbon*c2n_d)
     &                   /(mw_nitrogen+mw_carbon*c2n_redfield)
      beta_p = 0.7
      beta_dt = 0.5
      sink_rate_dt=10.0
      par = 0.41
      c_mol_wt = 12.0
      chl2pig = 0.8
      rain_ratio = 0.0070
      lysocline = 15
       
      REWIND(5) 
      READ(5,OCEAN_CC)
*/
*DECLARE BIOLOGY
*D BIOLOGY.190,201
      const1_p = 1.0 - beta_p
      const1_dt = 1.0 - beta_dt
      const2_p = c2n_p - beta_p*c2n_z
      const2_dt = c2n_d - beta_dt*c2n_z
      const3_p = 1.0 - beta_p*c2n_p/c2n_z
      const3_dt = 1.0 - beta_dt*c2n_d/c2n_z
      const4_p = c2n_p * const1_p
      const4_dt = c2n_d * const1_dt
      const5 = c2n_p - c2n_d
      const6 = c2n_p / c2n_d
      const7 =  c2n_z/c2n_d
      const8 = 1.0 - const7
