*DECLARE FTSA1A
*D ARE2F404.308
     &     S, S_SEA, RGRAIN, SOOT,
*I AWI2F400.5
     &     ,S_SEA(L1)                         ! Snow amount on sea-ice
*D AJG1F405.62
               if (s_sea(j).gt.0.0) then
*D AJG1F405.69
     &           +(snow_albedo-alphab)*(1.0-exp(-maskd*s_sea(j)))
*DECLARE RAD_CTL1
*I RAD_CTL1.321 
     &      D1(JSNODEP_SEA+JS),
*I AWI1F403.293
     &      D1(JSNODEP_SEA+JS),
*DECLARE INITA2O1
*D INITA2O1.506,INITA2O1.507
CL 2.8 Snow depth on atmos grid (specifically over sea ice)
      JA_SNOWDEPTH=JSNODEP_SEA
*DECLARE STATMPT1
*I GDR4F305.312
      JSNODEP_SEA    = SI( 95,Sect_No,im_index) ! Snow depth on sea ice
*DECLARE TYPPTRA
*I GRB0F304.269  
     &       JSNODEP_SEA,                ! snow depth on sea-ice 
*B ABX1F404.46
     &  JSNODEP_SEA,
