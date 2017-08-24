*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  cap_sflintsea.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*DECLARE SFLINT7A
*I SFLINT7A.370
c IF SQRTCD_K too small, CDTEMP1 explodes and CDR10M goes NaN
            IF ( CDR10M(I).ne.CDR10M(I) ) CDR10M(I)=0.
            IF ( abs(CDR10M(I)).gt.1e6 ) CDR10M(I)=0.
*I SFLINT7A.410
c IF SQRTCD_K too small, CDTEMP1 explodes and CHR1P5M goes -Inf
            IF ( abs(CHR1P5M(I)).gt.1e6 ) CHR1P5M(I)=0.
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  ccouple_moses2.2
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*IDENT ccouple
*/
*/ Coupling mod for FAMOUS with MOSES 2.2 & coastal tiling. 
*/
*/ Based on the mod ccouple (t20rt) with updated STASH codes etc.
*/ 
*/ Decks modified: DOARAV1,TRANA2O1,TRANO2A1,SWAPA2O2,SWAPO2A2,INITA2O1
*/
*/ Annette Osprey 3 September 2007
*/
*/
*DC DOARAV1
*D OJG1F403.3
     &,DATA_TARG,ADJUST_TARG,ICODE,CMESSAGE)                            
*I OJG1F403.28    
     &,ADJUST_TARG(LROW_TARG,*)
*D DOARAV1.82,DOARAV1.83   
          IF (MASK_TARG(IX2,IY2).EQV.WANT) THEN                         
          IF (COUNT_TARG(IX2,IY2).NE.0) THEN                            
*I DOARAV1.91    
            ELSE
            IF(ADJUST.EQ.5)THEN
            TEMP_TARG=0.0
            ELSEIF(ADJUST.EQ.6)THEN
            TEMP_TARG=1.0
            ELSE
            TEMP_TARG=DATA_TARG(IX2,IY2)
            ENDIF
            ENDIF
*I OJG1F403.54    
            ELSEIF (ADJUST.EQ.3) THEN
              ADJUST_TARG(IX2,IY2)=DATA_TARG(IX2,IY2)-TEMP_TARG
            ELSEIF (ADJUST.EQ.4) THEN
              IF (TEMP_TARG.EQ.0) THEN
                ADJUST_TARG(IX2,IY2)=1.0
              ELSE
                ADJUST_TARG(IX2,IY2)=DATA_TARG(IX2,IY2)/TEMP_TARG
              ENDIF
            ELSEIF (ADJUST.EQ.5) THEN
              DATA_TARG(IX2,IY2)=DATA_TARG(IX2,IY2)+TEMP_TARG
            ELSEIF (ADJUST.EQ.6) THEN
              DATA_TARG(IX2,IY2)=DATA_TARG(IX2,IY2)*TEMP_TARG
*DC TRANA2O1
*I TRANA2O1.40    
     + FRAC,
*I TRANA2O1.120   
     + FRAC(ICOLS,JROWS),    ! IN FRACTIONAL LAND MASK
*I TRANA2O1.169   
     +,ADJUSTMENT(ICOLS,JROWS) ! RATIO/DELTA stored for 2nd do_areaver 
     +                       ! call. Assumes atm size:ICOLS & JROWS
     +,DUMMY(IMT,JMT)
*I OJG1F403.58    
     &,index_back(maxl)
*I OJG1F403.59    
     &,backweight(maxl)
*D OJG1F403.60
*D OJG1F403.64
*I OJG1F403.67    
*IF DEF,CSRV_TWO_WAY
C
C    New arrays required to invert the coupling fields during 
C    conservative interpolation
C
      REAL
     & wmixinv(imt,jmt)      ! wind mixing power 
     &,blueinv(imt,jmt)      ! penetrative solar
     &,heatinv(imt,jmt)      ! non penetrative heat fluxes
     &,pmininv(imt,jmt)      ! precipitation - evaporation
     &,riverinv(imt,jmt)     ! river outflow
     &,snowinv(imt,jmt)      ! snowfall
     &,sublminv(imt,jmt)     ! sublimation
     &,btmltinv(imt,jmt)     ! sea ice diffusive heat flux
     &,tpmltinv(imt,jmt)     ! sea ice top melt heat flux
C-----------------------------------------------------------------------
C
C    The coupling is done in 3 stages:
C
C    1  Atmospheric fields are interpolated onto the ocean grid
C    2  The ocean fields are back interpolated onto the atmospheric 
C       grid, and the ratio/difference between the original atmosphere. 
C       and the back-interpolated atmosphere is calculated.
C    3  The ratio/difference from 2 is then applied to the ocean 
C       field to give the correct oceanic values for conservation to 
C       be maintained.
C
C    The logic is as follows (assuming the atmosphere and ocean grids 
C    are orientated oppositely in latitude).
C
C    Pre_areaver calculates weights for A==>O and O==>A. The weights 
C    are valid for atmospheric mask amasktp, and inverted ocean, with 
C    inverted mask omaskd.
C
C    Stage 2: use Do_areaver from O==>A. Invert =true, results go onto 
C             atmospheric mask.
C
C    Stage 3: use Do_areaver from A==>O. Invert =false. The weights will
C             assume the ocean is inverted, therefore it must be 
C             flipped before do_areaver is used, and omaskd must 
C             be specified. The resulting ocean field will be upside dow
C             therefore it must be inverted again after do_areaver has 
C             been called.
C
C
C   Summary
C
C   *    pre_areaver called twice, A==>O(inverted) and O(inverted)==>A
C   *    field X is interpolated from A==>O
C   *    do_areaver O(inv)==>A produces the differences/ratios on the 
C        atmospheric mask
C   *    ocean field is inverted
C   *    do_areaver A==>O(inv) calculates new ocean field
C   *    ocean field is inverted
C=======================================================================
C
*I OOM3F403.4     
       REAL VALMAX,VALMIN             ! Field A==>O test variables
       REAL WRKMAX,WRKMIN
       PARAMETER(WRKMAX=-1.0e10,WRKMIN=1.0e10)
*I TRANA2O1.308   
      VALMAX=WRKMAX
      VALMIN=WRKMIN
      DO J=1,JMTM1
      DO I=1,IMT
      IF(USTRSOUT(I,J).GT.VALMAX)VALMAX=USTRSOUT(I,J)
      IF(USTRSOUT(I,J).LT.VALMIN)VALMIN=USTRSOUT(I,J)
      ENDDO
      ENDDO
*     WRITE(6,*)'USTRSOUT A==>O  MAX/MIN ',VALMAX,VALMIN
      VALMAX=WRKMAX
      VALMIN=WRKMIN
      DO J=1,JMTM1
      DO I=1,IMT
      IF(VSTRSOUT(I,J).GT.VALMAX)VALMAX=VSTRSOUT(I,J)
      IF(VSTRSOUT(I,J).LT.VALMIN)VALMIN=VSTRSOUT(I,J)
      ENDDO
      ENDDO
*     WRITE(6,*)'VSTRSOUT A==>O  MAX/MIN ',VALMAX,VALMIN
*I OJG1F403.73    
*IF DEF,CSRV_TWO_WAY
      call pre_areaver(irt,xuo,jmt,ocyd,global,imt,.false.
     &,omaskd,icols,xua,jrows,yua,.true.,.true.
     &,lenl,count_a,base_a,index_arav,weight,icode,cmessage)
      lenl=maxl
      call pre_areaver(icols,xua,jrows,yua,.true.,icols,.false.
     &,amasktp,irt,xuo,jmt,ocyd,global,.true.
     &,lenl,count_o,base_o,index_back,backweight,icode,cmessage)
*ELSE
*I TRANA2O1.374   
*ENDIF                                                                  
*I OJG1F403.80    
*IF DEF,CSRV_TWO_WAY
      call do_areaver(irt,jmt,imt,invert,wmixout,icols,jrows
     &,count_a,base_a,icols,.false.,amasktp,index_arav,weight,4
     &,wmixin,adjustment,icode,cmessage)
      if (invert.eqv..true.) then
        do j=1,jmt
          do i=1,imt
            wmixinv(i,j)=wmixout(i,jmt+1-j)
          enddo
        enddo
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &   count_o,base_o,imt,.false.,omaskd,index_back,backweight,6,
     &   wmixinv,dummy,icode,cmessage)
        do j=1,jmt
          do i=1,imt
            wmixout(i,j)=wmixinv(i,jmt+1-j)
          enddo
        enddo
      else 
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,6,
     &    wmixout,dummy,icode,cmessage)
      endif
*ELSE
*D OJG1F403.83
     &,wmixin,adjustment,icode,cmessage)                                
*ENDIF
*D CJG6F401.95
     &,work,adjustment,icode,cmessage)                                  
*I TRANA2O1.387   
      VALMAX=WRKMAX
      VALMIN=WRKMIN
      DO J=1,JMT
      DO I=1,IMT
      IF(OMASKD(I,J).EQV..FALSE.)THEN
      IF(WMIXOUT(I,J).GT.VALMAX)VALMAX=WMIXOUT(I,J)
      IF(WMIXOUT(I,J).LT.VALMIN)VALMIN=WMIXOUT(I,J)
      ENDIF
      ENDDO
      ENDDO
*     WRITE(6,*)'WMIXOUT A==>O  MAX/MIN ',VALMAX,VALMIN
*I OJG1F403.87    
*IF DEF,CSRV_TWO_WAY
      call do_areaver(irt,jmt,imt,invert,blueout,icols,jrows
     &,count_a,base_a,icols,.false.,amasktp,index_arav,weight,4
     &,bluein,adjustment,icode,cmessage)
      if (invert.eqv..true.) then
        do j=1,jmt
          do i=1,imt
            blueinv(i,j)=blueout(i,jmt+1-j)
          enddo
        enddo
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,6,
     &    blueinv,dummy,icode,cmessage)
        do j=1,jmt
          do i=1,imt
            blueout(i,j)=blueinv(i,jmt+1-j)
          enddo
        enddo
      else
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,6,
     &    blueout,dummy,icode,cmessage)
      endif
*ELSE
*D OJG1F403.90
     &,bluein,adjustment,icode,cmessage)                                
*ENDIF
*D CJG6F401.102
     &,work,adjustment,icode,cmessage)                                  
*I TRANA2O1.412   
      VALMAX=WRKMAX
      VALMIN=WRKMIN
      DO J=1,JMT
      DO I=1,IMT
      IF(OMASKD(I,J).EQV..FALSE.)THEN
      IF(BLUEOUT(I,J).GT.VALMAX)VALMAX=BLUEOUT(I,J)
      IF(BLUEOUT(I,J).LT.VALMIN)VALMIN=BLUEOUT(I,J)
      ENDIF
      ENDDO
      ENDDO
*     WRITE(6,*)'BLUEOUT A==>O  MAX/MIN ',VALMAX,VALMIN
*I OJG1F403.94    
*IF DEF,CSRV_TWO_WAY
      call do_areaver(irt,jmt,imt,invert,heatflux,icols,jrows
     &,count_a,base_a,icols,.false.,amasktp,index_arav,weight,3
     &,worka,adjustment,icode,cmessage)
      if (invert.eqv..true.) then
        do j=1,jmt
          do i=1,imt
            heatinv(i,j)=heatflux(i,jmt+1-j)
          enddo
        enddo
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,5,
     &    heatinv,dummy,icode,cmessage)
        do j=1,jmt
          do i=1,imt
            heatflux(i,j)=heatinv(i,jmt+1-j)
          enddo
        enddo
      else
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,5,
     &    heatflux,dummy,icode,cmessage)
      endif
*ELSE
*D OJG1F403.97
     &,worka,adjustment,icode,cmessage)                                 
*ENDIF
*D CJG6F401.109
     &,work,adjustment,icode,cmessage)                                  
*I TRANA2O1.438   
      VALMAX=WRKMAX
      VALMIN=WRKMIN
      DO J=1,JMT
      DO I=1,IMT
      IF(OMASKD(I,J).EQV..FALSE.)THEN
      IF(HEATFLUX(I,J).GT.VALMAX)VALMAX=HEATFLUX(I,J)
      IF(HEATFLUX(I,J).LT.VALMIN)VALMIN=HEATFLUX(I,J)
      ENDIF
      ENDDO
      ENDDO
*     WRITE(6,*)'HEATFLUX A==>O  MAX/MIN ',VALMAX,VALMIN
*I OJG1F403.101   
*IF DEF,CSRV_TWO_WAY
      call do_areaver(irt,jmt,imt,invert,pminuse,icols,jrows
     &,count_a,base_a,icols,.false.,amasktp,index_arav,weight,3
     &,worka,adjustment,icode,cmessage)
      if (invert.eqv..true.) then
        do j=1,jmt
          do i=1,imt
            pmininv(i,j)=pminuse(i,jmt+1-j)
          enddo
        enddo
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,5,
     &    pmininv,dummy,icode,cmessage)
        do j=1,jmt
          do i=1,imt
            pminuse(i,j)=pmininv(i,jmt+1-j)
          enddo
        enddo
      else
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,5,
     &    pminuse,dummy,icode,cmessage)
      endif
*ELSE
*D OJG1F403.104
     &,worka,adjustment,icode,cmessage)                                 
*ENDIF
*D CJG6F401.116
     &,work,adjustment,icode,cmessage)                                  
*I TRANA2O1.469   
      VALMAX=WRKMAX
      VALMIN=WRKMIN
      DO J=1,JMT
      DO I=1,IMT
      IF(OMASKD(I,J).EQV..FALSE.)THEN
      IF(PMINUSE(I,J).GT.VALMAX)VALMAX=PMINUSE(I,J)
      IF(PMINUSE(I,J).LT.VALMIN)VALMIN=PMINUSE(I,J)
      ENDIF
      ENDDO
      ENDDO
*     WRITE(6,*)'PMINUSE A==>O  MAX/MIN ',VALMAX,VALMIN
*D TRANA2O1.488
c  Note River run off if there is some land
          IF (FRAC(I,J).GT.0.) THEN
*I TRANA2O1.492   
     &        FRAC(I,J)/(1.-FRAC(K,L)) *
*I OJG1F403.108   
*IF DEF,CSRV_TWO_WAY
      call do_areaver(irt,jmt,imt,invert,riverout,icols,jrows
     &,count_a,base_a,icols,.false.,amasktp,index_arav,weight,4
     &,worka,adjustment,icode,cmessage)
      if (invert.eqv..true.) then
        do j=1,jmt
          do i=1,imt
            riverinv(i,j)=riverout(i,jmt+1-j)
          enddo
        enddo
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,6,
     &    riverinv,dummy,icode,cmessage)
        do j=1,jmt
          do i=1,imt
            riverout(i,j)=riverinv(i,jmt+1-j)
          enddo
        enddo
      else
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,6,
     &    riverout,dummy,icode,cmessage)
      endif
*ELSE
*D OJG1F403.111
     &,worka,adjustment,icode,cmessage)                                 
*ENDIF
*D CJG6F401.123
     &,work,adjustment,icode,cmessage)                                  
*I TRANA2O1.506   
      VALMAX=WRKMAX
      VALMIN=WRKMIN
      DO J=1,JMT
      DO I=1,IMT
      IF(OMASKD(I,J).EQV..FALSE.)THEN
      IF(RIVEROUT(I,J).GT.VALMAX)VALMAX=RIVEROUT(I,J)
      IF(RIVEROUT(I,J).LT.VALMIN)VALMIN=RIVEROUT(I,J)
      ENDIF
      ENDDO
      ENDDO
*     WRITE(6,*)'RIVEROUT A==>O  MAX/MIN ',VALMAX,VALMIN
*I OJG1F403.115   
*IF DEF,CSRV_TWO_WAY
      call do_areaver(irt,jmt,imt,invert,snowout,icols,jrows
     &,count_a,base_a,icols,.false.,amasktp,index_arav,weight,4
     &,worka,adjustment,icode,cmessage)
      if (invert.eqv..true.) then
        do j=1,jmt
          do i=1,imt
            snowinv(i,j)=snowout(i,jmt+1-j)
          enddo
        enddo
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,6,
     &    snowinv,dummy,icode,cmessage)
        do j=1,jmt
          do i=1,imt
            snowout(i,j)=snowinv(i,jmt+1-j)
          enddo
        enddo
      else
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,6,
     &    snowout,dummy,icode,cmessage)
      endif
*ELSE
*D OJG1F403.118
     &,worka,adjustment,icode,cmessage)                                 
*ENDIF
*D CJG6F401.131
     &,work,adjustment,icode,cmessage)                                  
*I TRANA2O1.531   
      VALMAX=WRKMAX
      VALMIN=WRKMIN
      DO J=1,JMT
      DO I=1,IMT
      IF(OMASKD(I,J).EQV..FALSE.)THEN
      IF(SNOWOUT(I,J).GT.VALMAX)VALMAX=SNOWOUT(I,J)
      IF(SNOWOUT(I,J).LT.VALMIN)VALMIN=SNOWOUT(I,J)
      ENDIF
      ENDDO
      ENDDO
*     WRITE(6,*)'SNOWOUT A==>O  MAX/MIN ',VALMAX,VALMIN
*I OJG1F403.122   
*IF DEF,CSRV_TWO_WAY
      call do_areaver(irt,jmt,imt,invert,sublmout,icols,jrows
     &,count_a,base_a,icols,.false.,amasktp,index_arav,weight,3
     &,sublmin,adjustment,icode,cmessage)
      if (invert.eqv..true.) then
        do j=1,jmt
          do i=1,imt
            sublminv(i,j)=sublmout(i,jmt+1-j)
          enddo
        enddo
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,5,
     &    sublminv,dummy,icode,cmessage)
        do j=1,jmt
          do i=1,imt
            sublmout(i,j)=sublminv(i,jmt+1-j)
          enddo
        enddo
      else
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,5,
     &    sublmout,dummy,icode,cmessage)
      endif
*ELSE
*D OJG1F403.125
     &,sublmin,adjustment,icode,cmessage)                               
*ENDIF
*D CJG6F401.138
     &,work,adjustment,icode,cmessage)                                  
*I TRANA2O1.544   
      VALMAX=WRKMAX
      VALMIN=WRKMIN
      DO J=1,JMT
      DO I=1,IMT
      IF(OMASKD(I,J).EQV..FALSE.)THEN
      IF(SUBLMOUT(I,J).GT.VALMAX)VALMAX=SUBLMOUT(I,J)
      IF(SUBLMOUT(I,J).LT.VALMIN)VALMIN=SUBLMOUT(I,J)
      ENDIF
      ENDDO
      ENDDO
*     WRITE(6,*)'SUBLMOUT A==>O  MAX/MIN ',VALMAX,VALMIN
*I OJG1F403.129   
*IF DEF,CSRV_TWO_WAY
      call do_areaver(irt,jmt,imt,invert,btmltout,icols,jrows
     &,count_a,base_a,icols,.false.,amasktp,index_arav,weight,3
     &,btmltin,adjustment,icode,cmessage)
      if (invert.eqv..true.) then
        do j=1,jmt
          do i=1,imt
            btmltinv(i,j)=btmltout(i,jmt+1-j)
          enddo
        enddo
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,5,
     &    btmltinv,dummy,icode,cmessage)
        do j=1,jmt
          do i=1,imt
            btmltout(i,j)=btmltinv(i,jmt+1-j)
          enddo
        enddo
      else
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,5,
     &    btmltout,dummy,icode,cmessage)
      endif
*ELSE
*D OJG1F403.132
     &,btmltin,adjustment,icode,cmessage)                               
*ENDIF
*D CJG6F401.145
     &,work,adjustment,icode,cmessage)                                  
*I TRANA2O1.557   
      VALMAX=WRKMAX
      VALMIN=WRKMIN
      DO J=1,JMT
      DO I=1,IMT
      IF(OMASKD(I,J).EQV..FALSE.)THEN
      IF(BTMLTOUT(I,J).GT.VALMAX)VALMAX=BTMLTOUT(I,J)
      IF(BTMLTOUT(I,J).LT.VALMIN)VALMIN=BTMLTOUT(I,J)
      ENDIF
      ENDDO
      ENDDO
*     WRITE(6,*)'BTMLTOUT A==>O  MAX/MIN ',VALMAX,VALMIN
*I OJG1F403.136   
*IF DEF,CSRV_TWO_WAY
      call do_areaver(irt,jmt,imt,invert,tpmltout,icols,jrows
     &,count_a,base_a,icols,.false.,amasktp,index_arav,weight,4
     &,tpmltin,adjustment,icode,cmessage)
      if (invert.eqv..true.) then
        do j=1,jmt
          do i=1,imt
            tpmltinv(i,j)=tpmltout(i,jmt+1-j)
          enddo
        enddo
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,6,
     &    tpmltinv,dummy,icode,cmessage)
        do j=1,jmt
          do i=1,imt
            tpmltout(i,j)=tpmltinv(i,jmt+1-j)
          enddo
        enddo
      else
        call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &    count_o,base_o,imt,.false.,omaskd,index_back,backweight,6,
     &    tpmltout,dummy,icode,cmessage)
      endif
*ELSE
*D OJG1F403.139
     &,tpmltin,adjustment,icode,cmessage)                               
*ENDIF
*D CJG6F401.152
     &,work,adjustment,icode,cmessage)                                  
*I TRANA2O1.570   
      VALMAX=WRKMAX
      VALMIN=WRKMIN
      DO J=1,JMT
      DO I=1,IMT
      IF(OMASKD(I,J).EQV..FALSE.)THEN
      IF(TPMLTOUT(I,J).GT.VALMAX)VALMAX=TPMLTOUT(I,J)
      IF(TPMLTOUT(I,J).LT.VALMIN)VALMIN=TPMLTOUT(I,J)
      ENDIF
      ENDDO
      ENDDO
*     WRITE(6,*)'TPMLTOUT A==>O  MAX/MIN ',VALMAX,VALMIN
*DC TRANO2A1
*D TRANO2A1.140
     +,I,J,K                 ! Working indices                          
*I TRANO2A1.202   
     &,DUMMY(IMT,JMT)          ! Dummy argument (for DO_AREAVER)
*I TRANO2A1.210   
      REAL VALMAX,VALMIN      ! Max & min values passed O--> A
      REAL WRKMAX,WRKMIN      ! Max & min work start points
      PARAMETER(WRKMAX=-1.0e10,WRKMIN=1.0e10)
*D OJG1F403.146
     &,0,UOUT,DUMMY,ICODE,CMESSAGE)                                     
*D OJG1F403.150
     &,0,VOUT,DUMMY,ICODE,CMESSAGE)  
*D OJG1F403.154
     &,0,AICEOUT,DUMMY,ICODE,CMESSAGE)                                  
*D TRANO2A1.536
          IF(INVERT)THEN
          K=JMT+1-J
          ELSE
          K=J
          ENDIF
          IF (OMASK(I,K).EQV..FALSE.) THEN                              
*D OJG1F403.158
     &,0,WORK_A,DUMMY,ICODE,CMESSAGE)                                   
*D TRANO2A1.570
        IF (AMASKTP(I,J).EQV..FALSE.) THEN                              
*D TRANO2A1.584
        IF(INVERT)THEN
        K=JMT+1-J
        ELSE
        K=J
        ENDIF
        IF (OMASK(I,K).EQV..FALSE.)                                     
*D OJG1F403.162
     &,0,HSNOWOUT,DUMMY,ICODE,CMESSAGE)                                 
*D OJG1F403.166
     &,0,WORK_A,DUMMY,ICODE,CMESSAGE)                                   
*D TRANO2A1.630
C*IF DEF,SEAICE                                                         
*D TRANO2A1.640,TRANO2A1.654  
C     DO 420 J = 1,JROWS                                                
C       DO 415 I = 1,ICOLS                                              
C         IF (WORK_A(I,J) .NE. RMDI) THEN                               
C           IF (AICEOUT(I,J) .EQ. 0.0) THEN                             
C             TSTAROUT(I,J) = WORK_A(I,J) + ZERODEGC                    
C           ELSEIF (AICEREF(I,J) .GE. AICEMIN) THEN                     
C             TSTAROUT(I,J) = TFS + ( AICEOUT(I,J)/AICEREF(I,J) )       
C    +                             *( TSTAROUT(I,J) - TFS )             
C           ELSE                                                        
C             TSTAROUT(I,J) = TFS                                       
C           ENDIF                                                       
C         ENDIF                                                         
C415     CONTINUE                                                       
C420   CONTINUE                                                         
C*ELSE                                                                  
*D TRANO2A1.657
          IF (AMASKTP(I,J) .EQV..FALSE.) THEN                           
*D TRANO2A1.660,TRANO2A1.662  
 425     CONTINUE                                                       
 430   CONTINUE                                                         
C*ENDIF
*DC SWAPA2O2
*B SWAPA2O2.159
      REAL FLAND_GLO(G_P_FIELD),FLAND_LOC(P_FIELD)
      LOGICAL AMASKTP(G_P_FIELD)
*I SWAPA2O2.345   
*IF DEF,TRANGRID,OR,DEF,RIVERS
      CALL FROM_LAND_POINTS(FLAND_LOC,D1(JFRAC_LAND),
     & atmos_landmask_local,
     & lasize(1)*lasize(2),number_of_landpts_out)

      CALL GATHER_FIELD(FLAND_LOC,FLAND_GLO,
     & lasize(1),lasize(2),glsize(1),glsize(2),
     & gather_pe,GC_ALL_PROC_GROUP,info)
      IF(mype.EQ.gather_pe) THEN
        DO I=1,G_P_FIELD
        AMASKTP(I)=(FLAND_GLO(I).EQ.1.0)
        ENDDO
      ENDIF
*ENDIF
*D SWAPA2O2.481   
     + AMASKTP,FLAND_GLO,
*DC SWAPO2A2
*I SWAPO2A2.162    
*IF DEF,TRANGRID
      LOGICAL
     &       AMASK(G_P_FIELD)
      REAL FLAND_GLO(G_P_FIELD),FLAND_LOC(P_FIELD)
     & ,TSTAR_SSI(P_FIELD)
      INTEGER NLANDPT
*ENDIF
*I SWAPO2A2.319   
*IF DEF,TRANGRID
C
C Unpack the fractional land mask, reassign land fraction of 0.0 to the 
C points where there is no land (previously missing data) and 
C define a mask which is .true. over all points that are land only.
C Transo2a will take effect only where this mask is .false. (some 
C sea).
C
      CALL FROM_LAND_POINTS(FLAND_LOC,D1(JFRAC_LAND),
     & atmos_landmask_local,
     & lasize(1)*lasize(2),NLANDPT)

*     write(6,*) 'lasize(1,2),nlandpt: ',lasize(1),lasize(2),nlandpt
*     write(6,*) 'atmos_landmask_local: ',atmos_landmask_local
*     write(6,*) 'FLAND_LOC: ',FLAND_LOC

      do i=1,p_field
        if (fland_loc(i).eq.rmdi) fland_loc(i)=0
      enddo
*     write(6,*) 'FLAND_LOC zeroed: ',FLAND_LOC


      CALL GATHER_FIELD(FLAND_LOC,FLAND_GLO,
     & lasize(1),lasize(2),glsize(1),glsize(2),
     & gather_pe,GC_ALL_PROC_GROUP,info)
      IF(info.NE.0) THEN      ! Check return code
         CMESSAGE='SWAPO2A : ERROR in gather of FLAND_GLO'
         ICODE=20
         GO TO 999
      ENDIF

      if(mype.eq.gather_pe) then
      DO I=1,G_P_FIELD
        AMASK(I)=(FLAND_GLO(I).EQ.1.0)
*       WRITE(6,*)'FLAND_GLO',AMASK(I),FLAND_GLO(I)
      ENDDO
      ENDIF
*ENDIF
*D SWAPO2A2.332
     & XUO,XTO,YUO,YTO,XTA,XUA,YTA,YUA,AMASK,                           
*I SWAPO2A2.407
*/ need to add this to the MPP version, since omitted initially
*/ back on the standard atmosphere decomposition
C
C     Diagnose the new grid box means for TSTAR based upon the
C     open sea surface temperature (field 507/TSTAR_SEA)
C     which is passed by TRANSO2A
C
      DO I=0,P_FIELD-1
        TSTAR_SSI(I+1)=0.0
c        IF(atmos_landmask_local(i+1).EQV..FALSE.)THEN
        IF(fland_loc(i+1).NE.1.)THEN
c       IF(AMASK(I+1).EQV..FALSE.)THEN
c          IF(LTLEADS.EQV..FALSE.) D1(JTSTAR_SEA+I)=TFS
          IF((LTLEADS.EQV..FALSE.).AND.(D1(JICE_FRACTION+I).GT.0))
     &          D1(JTSTAR_SEA+I)=TFS
          D1(JTSTAR_SICE+I)=AMIN1(D1(JTSTAR_SICE+I),TFS)
          TSTAR_SSI(I+1) =D1(JICE_FRACTION+I) *D1(JTSTAR_SICE+I)+
     &                 (1.0-D1(JICE_FRACTION+I))*D1(JTSTAR_SEA+I)
        ENDIF
        D1(JTSTAR+I)=FLAND_LOC(I+1) *D1(JTSTAR_LAND+I)+
     &          (1.0-FLAND_LOC(I+1))*TSTAR_SSI(I+1)
      ENDDO

      CALL SWAPBOUNDS(D1(JTSTAR),lasize(1),lasize(2),offx,offy,
     &                   swap_levels)
      CALL SET_SIDES(D1(JTSTAR),lasize(1)*lasize(2),lasize(1),
     &                   swap_levels,fld_type_p)

      CALL SWAPBOUNDS(D1(JTSTAR_SICE),lasize(1),lasize(2),offx,offy,
     &                   swap_levels)
      CALL SET_SIDES(D1(JTSTAR_SICE),lasize(1)*lasize(2),lasize(1),
     &                   swap_levels,fld_type_p)

      CALL SWAPBOUNDS(D1(JTSTAR_SEA),lasize(1),lasize(2),offx,offy,
     &                   swap_levels)
      CALL SET_SIDES(D1(JTSTAR_SEA),lasize(1)*lasize(2),lasize(1),
     &                   swap_levels,fld_type_p)


*DC INITA2O1
*D GSS2F305.117
      CALL FINDPTR(ATMOS_IM, 3,392,                                     
*D INITA2O1.99
        ICODE=3392                                                      
*D GSS2F305.118
      CALL FINDPTR(ATMOS_IM, 3,394,                                     
*D INITA2O1.121
        ICODE=3394                                                      
*D GSS2F305.121
      CALL FINDPTR(ATMOS_IM, 1,260,                                     
*D INITA2O1.179
        ICODE=1260                                                      
*D GSS2F305.129
      CALL FINDPTR(ATMOS_IM, 3,353,                                     
*D INITA2O1.329
        ICODE=3353                                                      
*D INITA2O1.449
      JA_TSTAR=JTSTAR_SEA                                               
*D INITA2O1.451
        ICODE=507                                                       
*D INITA2O1.522
        ICODE=507                                                       
*/
*/       - End of mod ccouple
*/
*/ ------------------------------------------------------------------
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  co2_coupling_landfix_060912.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID CO2_COUPLE
*/ Enable coupling of ocean-atmosphere CO2 for FAMOUS. 4.5.1 code as it was
*/ only worked for HadCM3LC, where ocean and atmosphere grids are identical
*/  - coupling needs to follow other FAMOUS fields in terms of interpolation,
*/ area averaging and coastal tiling. In addition, the GATHER_FIELD in 
*/ SWAPO2A just hung on HPCx/Ruby for nprocs>1
*/
*/ 16-12-11: multiply land co2 by fractional land coverage to get coastal
*/ fluxes correct - bug common to all UM until 2011
*/
*/ 06-09-12: move triffid counter to correct place (ABX5F406), include
*/ mass carbon ->CO2 conversion when updating CO2 concentrations from
*/ land carbon
*DECLARE ATMPHY1
*I GDR5F305.16
C-----------------------------------------------------------------------
C Increment counter for number of atmosphere timesteps since last
C call to TRIFFID vegetation model
C-----------------------------------------------------------------------
      A_INTHD(23) = A_INTHD(23) + 1
*D ABX1F404.288,291
*DECLARE BL_CTL1
*I ACN1F405.6 
     &   FLAND_GLO(P_FIELDDA),
*D ACN1F405.30,ACN1F405.31
            LAND_CO2_L(I) =   (RESP_S(I) ! soil respiration
     &                      - WORK8(I))  ! NPP
c    !!!!THIS SHOULD BE *(M_CO2/M_CARBON), by comparison with 
c    !!!!UPSTREAM CODE
c    !&                      * (M_CO2 / M_CARBON)
     &                      * (28.966*1.5194)/12.0
*I ACN1F405.35
          CALL FROM_LAND_POINTS(FLAND_GLO,D1(JFRAC_LAND),
     &                         D1(JLAND),P_FIELD,LAND_FIELD)
*D ACN1F405.43
          CO2_FLUX(I) = CO2_FLUX(I) + D1(J_CO2_EMITS+I-1)*FLAND_GLO(I)
*D ACN1F405.54
          CO2_FLUX(I) = CO2_FLUX(I) + LAND_CO2(I)*FLAND_GLO(I)
*/
*DECLARE INITA2O1
*D CCN1F405.77,CCN1F405.80
*DECLARE SWAPO2A2
*D CCN1F405.251,253
c
c GATHER_FIELD doesn't appear to work for input lengths that aren't
c glsize(1),glsize(2),lasize(1),lasize(2) - this call has been changed 
c to a GENERAL_GATHER. We need the D1 index number to do this
c 
        do i=1,N_OBJ_D1_MAX
          if (d1_addr(d1_address,i,ocean_sm).eq.JO_co2flux) exit
        end do

        CALL GENERAL_GATHER_FIELD(D1(JO_co2flux),CO2FLUX,
     &  (lasize(1)-2)*lasize(2),(glsize(1)-2)*glsize(2),
     &   D1_ADDR(d1_object_type,i,ocean_sm),gather_pe,
     &   ICODE,CMESSAGE)
*D CCN1F405.262
c the resultant GATHERed field appears to be shifted one j line
c up relative to other imt*jmt fields - so shift things back down
          do j=1,g_jmt-1
*D CCN1F405.265 
              O_CO2FLUX(i+(j-1)*g_imt) = CO2FLUX(i+(j)*(g_imt-2))
*I CCN1F405.267
          do i=1,g_imt-2
            O_CO2FLUX(i+(g_jmt-1)*g_imt) = RMDI
          enddo  ! i
*I CCN1F405.271
c nasty - at the beginning of a model run, the pure-diagnostic CO2 flux
c field may not have any data. It's initialised as RMDI, but these aren't 
c detected over the sea and they destroy the atmospheric pCO2.
          do i=1,g_imt*g_jmt
            if (stepim(o_im).eq.0
     &     .and.o_flddepc(i).gt.0
     &     .and.o_co2flux(i).eq.RMDI
     &         ) then
              o_co2flux(i)=0.
              write(6,*)"RESETTING OCEAN CO2 FLUX TO 0"
            end if
          end do
*I CCN1F405.302 
c need to multiply by the sea-fraction to end up conserving area total
c co2 passed between O and A
        if (mype.eq.gather_pe) then
        DO I=1,G_P_FIELD
          A_CO2FLUX(I)=A_CO2FLUX(I)*(1-FLAND_GLO(I))
        END DO
        end if
*DECLARE TRANA2O1
*I CCN1F405.185
     +,ATMCO2_INV(IMT,JMT)   ! atm. CO2
*D CCN1F405.187
*D CCN1F405.191,CCN1F405.193
*B CCN1F405.199
*IF -DEF,TRANGRID 
*I CCN1F405.200
*ELSE
c follow method used for other fields to get a conserved atmos CO2 
c on the finer ocean grid.
c
        CALL H_INT_BL(JROWS,ICOLS,SEAPOINTS,INDEXBL,INDEXBR,atmco2
     &  ,             WEIGHTBL,WEIGHTBR,WEIGHTTL,WEIGHTTR,WORK)
        CALL POST_H_INT(NCOASTAL,INDEXA,INDEXO,ATPOINTS,atmco2,AMINT
     &  ,SEAPOINTS,WORK,OCPOINT,RMDI,OCPOINTS,atmco2_OUT)
        call do_areaver(irt,jmt,imt,invert,atmco2_OUT,icols,jrows
     &  ,count_a,base_a,icols,.false.,amasktp,index_arav,weight,4
     &  ,atmco2,adjustment,icode,cmessage)
        if (invert.eqv..true.) then
          do j=1,jmt
          do i=1,imt
            atmco2_inv(i,j)=atmco2_out(i,jmt+1-j)
          enddo
          enddo
          call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &      count_o,base_o,imt,.false.,omaskd,index_back,backweight,6,
     &      atmco2_inv,dummy,icode,cmessage)
          do j=1,jmt
          do i=1,imt
            atmco2_out(i,j)=atmco2_inv(i,jmt+1-j)
          enddo
          enddo
        else
          call do_areaver(icols,jrows,icols,.false.,adjustment,irt,jmt,
     &      count_o,base_o,imt,.false.,omaskd,index_back,backweight,6,
     &      atmco2_out,dummy,icode,cmessage)
        endif
*ENDIF
*I CCN1F405.207
            else
              ATMCO2_OUT(i,j) = RMDI
*D CCN1F405.212 

*DECLARE TRANO2A1
*D CCN1F405.326
*B CCN1F405.332
*IF -DEF,TRANGRID
*I CCN1F405.333
*ELSE
c do simple area-averaging onto atmos grid 
c
        CALL COPYO2A(IMT,JMT,co2fluxin,.FALSE.,OMASK
     &,.FALSE.,INVERT,IMT,WORK_O)

        CALL DO_AREAVER(IRT,JMT,
     &                  IMT,.FALSE.,work_o,
     &                  ICOLS,JROWS,
     &                  COUNT_A,BASE_A,
     &                  ICOLS,.FALSE.,AMASKTP,
     &                  POINT_O,WEIGHT,0,
     &                  co2fluxout,DUMMY,
     &                  ICODE,CMESSAGE)

*ENDIF
*D CCN1F405.343
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  ice_drift_fix
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID ice_drift_fix
*DECLARE ICEFREED
*I ICEFREED.126
     & ,utest(imt,jmtm1)    !Velocity vestor to carry values across
     & ,vtest(imt,jmtm1)    !cnvstop, to avoid o-i stress problems
*I ICEFREED.207
       utest(i,j)=uice(i,j)
       vtest(i,j)=vice(i,j)
*D ICEFREED.220,ICEFREED.221
               u1 = utest(i,j) - ucurrent(i,j)
               v1 = vtest(i,j) - vcurrent(i,j)
*D ICEFREED.228,ICEFREED.229
               mfx(i,j) = mf*vtest(i,j)
               mfy(i,j) = -mf*utest(i,j)
*COMPILE CNVSTOP
*DECLARE BLOKINIT
*D ODC1F405.21
     & ISX(IMT_idr,JMTM1_idr)        ! IN Stress under sea ice
*D OLA0F404.20,23
        CALL SWAPBOUNDS(WSX_LEADS,IMT,JMTM1_idr,O_EW_HALO,O_NS_HALO,1)
        CALL SWAPBOUNDS(WSY_LEADS,IMT,JMTM1_idr,O_EW_HALO,O_NS_HALO,1)
        CALL SWAPBOUNDS(ISX,IMT,JMTM1_idr,O_EW_HALO,O_NS_HALO,1)
        CALL SWAPBOUNDS(ISY,IMT,JMTM1_idr,O_EW_HALO,O_NS_HALO,1)
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  icedrift_090912.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*DECLARE ICEDRIFT
*D ICEDRIFT.286,ICEDRIFT.289
         aice(1,j)    = aice(imtm1,j)
         aice(imt,j)  = aice(2,j)
         hsnow(1,j)   = hsnow(imtm1,j)
         hsnow(imt,j) = hsnow(2,j)
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  inlansea_biogeo_safer.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*DECLARE ROWCALC
*I OOM1F405.90
     &,D1(joc_vflux_mask),RMDI
*DECLARE TRACER
*I OOM1F405.151
     &,vflux_mask,RMDI
*B TRACER.201
      REAL vflux_mask(IMT,JMT),RMDI
*B TRACER.1267
c do for actual inland seas, as specified by the salmask
      if (L_OCARBON) then
         DO I=1,IMT
           IF (vflux_mask(i,j).lt.0 .AND. FM(I,1).gt.0) THEN
               do n=3,nt
               do k=1,km
                 TA(I,k,n)=RMDI
               end do
               end do
               CO2_FLUX(I)=0.
               VTCO2_FLUX(I)=0.
               VALK_FLUX(I)=0.
               INVADE(I)=0.
               EVADE(I)=0.
           ENDIF
         END DO
      end if
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  moses2.2_ctile_060912.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*IDENT MOSES_CTILE
*/ 
*/ MOSES 2.2 + Coastal Tiling for FAMOUS
*/ 
*/ Based on the MOSES 2.2 mods for vn4.5:
*/   abx0f406, abx1f406, abx2f406, abx3f406, abx4f406, amv1f406,
*/   apa1f406, are1f406, are2f406, are3f406, newdecks, scatter_fix,
*/   update_m22.mdk 
*/
*/ With coastal tiling for MOSES 2.2 derived from:
*/   ang0f503.mf77 (Nic Gedney, Apr 01)
*/  
*/ Also includes:
*/   adb1f406: radiation mod (J.M. Edwards, May 99)
*/   Tibet_snow.mf77
*/   runoff_MII.mf77
*/  
*/ Should be used with: 
*/   ccouple_moses2.2: FAMOUS coupling mod
*/   LANDTOGLOBAL, PTOUV_LAND, PTOUV_SEA: coastal tiling mods
*/   amv2f406: *optional* MOSES 2.2 mod for radiative canopy
*/   arerf406_ctile: reconfiguration mod
*/ 
*/ Decks modified: ADDRES1, ATMPHY1, ATMSTEP1, BLEND_H, BL_CTL1, 
*/   BL_IC7A, BOUYTQ5B, BOUYTQ6A, CLDCTL1, CNTLATM, COMPT2A, CRADINCS, 
*/   C_ROUGH, DARCY5A, DECAY2A, DESCENT, DFPLN3A, DIAG3A, EXFXUV5B, 
*/   FCDCH6A, FCDCH7A, FILL3A, FTSA1A, FXCA3A,G ROWT2A, HYDCON5A, 
*/   HYDROL7A, HYDR_CT1, HYD_IC7A, IMCLTQ7A, INITIAL1, INITVEG1, 
*/   LEAF7A, LOTKA2A, LWRAD3A, NAMSIZE, NVEGPARM, PHIMH6A, PHYSIO7A, 
*/   PPARM1A, PRELIM1, PSLIMS1, RAD_CTL1, READFL1A, ROOTFR7A, RPANCA1A, 
*/   SCREEN7A, SEED, SETMODL1, SFEVAP7A, SFEXCH7A, SFFLUX7A, SFLINT6A, 
*/   SFLINT7A, SFMELT7A, SFREST7A, SFRIB7A, SFSNOW7A, SICEHT5B, 
*/   SMCEXT7A, SOILHT7A, SOILHY5A, SOILMC7A, SPARM1A, STATMPT1, 
*/   STDEV7A, SURFHY7A, SWRAD3A, TRIF, TRIFD2A, TYPPTRA, TYPSIZE, 
*/   U_MODEL1, UPANCIL1, VEG1A, VEG2A, VEG_CTL1, VEG_IC1A, VEG_IC2A
*/ 
*/ Decks added: ALBPFT, ALBSNOW, BDY_EXPL7A, BDY_EXPL8A, BDY_IMPL8A,
*/   CANCAP7A, GAUSS7A, IMBLPT18A, IMBLPT28A, IMSFPT8A, RESTILE7A,
*/   SFEXPL8A, SFIMPL8A, SICEHT7A, SOILEV7A, TILEALB
*/
*/ Annette Osprey 25 January 2008
*/
*/
*DECLARE ADDRES1
*D ABX2F404.79,ABX2F404.80   
! All surface tiles
          ILOUT = NTILES
*DECLARE ATMPHY1
*D ARE1F404.16
     &     TILE_FIELDDA,NTILESDA,                                       
*I AJS1F401.177   
     &        NTILESDA,             ! dynamic allocation: NTILES        
*D ARE2F404.3,ARE2F404.8    
*D ARE1F404.20,ARE1F404.26   
     &     ECAN_TILE(TILE_FIELDDA,NTILESDA)! Canopy evaporation from    
!                                          ! land tiles                 
     &    ,MELT_TILE(TILE_FIELDDA,NTILESDA)! Snowmelt on land tiles     
     &    ,LW_DOWN(P_FIELDDA)              ! Downward LW radiation      
     &    ,DOLR(P_FIELDDA)                 ! TOA - surface upward LW rad
     &    ,SW_TILE(TILE_FIELDDA,NTILESDA)  ! SW radiation on land tiles 
     &    ,TILE_FRAC(TILE_FIELDDA,NTILESDA)! Land-surface type fractions
*D ARE1F404.37
*D ARE1F404.39,ARE1F404.52   
*I ATMPHY1.142   
!  Similar for land surface albedos - different values are needed for   
!    direct & diffuse sunlight if the HadCM2 approximate treatment of   
!    sulphate aerosol is being used:                                    
      IF ( L_H2_SULPH ) THEN                                            
         NLALBS = 2                                                     
       ELSE                                                             
         NLALBS = 1                                                     
      ENDIF                                                             
                                                 
*D ARE2F404.10
     &     ,NTILES,TILE_FIELDDA,DOLR,LW_DOWN,SW_TILE
     &     ,D1(JLAND_ALB+1),D1(JSICE_ALB+1)                    
*D ARN1F404.92
     &     ,WORKF,BL_LEVELS,L_RADHEAT,RADHEAT_DIM1,NLALBS
*D AWI1F403.120,AWI1F403.127  
*D ARE2F404.12,ARE2F404.17   
      SAL_DIM = 1                                                       
*IF DEF,A03_7A                                                          
! Extra workspace required if MOSES II is used                          
      SAL_DIM = P_FIELD                                                 
*ENDIF                                                                  
*D ARE2F404.18
     &  NTILES,TILE_FIELDDA,DOLR,LW_DOWN,SW_TILE,
     &  D1(JLAND_ALB+1),D1(JSICE_ALB+1),                       
*D ARE1F404.53,ARE1F404.55   
     &            NTILES,TILE_FIELDDA,TILE_PTS,TILE_INDEX,              
     &            DOLR,LW_DOWN,SW_TILE,ECAN_TILE,MELT_TILE,TILE_FRAC,   
*D ARE1F404.56,ARE1F404.57   
     &              NTILES,TILE_FIELDDA,TILE_PTS,TILE_INDEX,            
     &              ECAN_TILE,MELT_TILE,TILE_FRAC,                      
*DECLARE ATMSTEP1
*I ACB2F405.41    
      ELSEIF ( (H_SECT(3).EQ."08A") )THEN
        L_RADHEAT = .TRUE.
        RADHEAT_DIM1 = P_FIELDDA
        TILE_FIELD = LAND_FIELD
*D ARE1F404.13
     &              TILE_FIELD,NTILES,
*DECLARE BLEND_H
*D BLEND_H.5
      PARAMETER (LB = 20.0)                                            
*DECLARE BL_CTL1
*D AJS1F401.219
     &           SOIL_EVAPORATION,SURF_HT_FLUX_LAND,SURF_RADFLUX,
*D ARE1F404.59,ARE1F404.61   
     &           NTILESDA,TILE_FIELDDA,TILE_PTS,TILE_INDEX,             
     &           OLR,LW_DOWN,SW_TILE,ECAN_TILE,MELT_TILE,TILE_FRAC,     
*D AJS1F401.228
     &       SURF_HT_FLUX_LAND(P_FIELDDA),         !                    
*I AJS1F401.229   
     &       SURF_HT_FLUX(P_FIELDDA),              
     &       SURF_HT_FLUX_SICE(P_FIELDDA),           
*I ARE1F404.65    
     &       NTILESDA,                             ! IN                 
*D ARE1F404.69,ARE1F404.74   
     &       OLR(P_FIELDDA),                       ! IN                 
     &       LW_DOWN(P_FIELDDA),                   ! IN                 
     &       SW_TILE(TILE_FIELDDA,NTILESDA),       ! IN                 
     &       ECAN_TILE(TILE_FIELDDA,NTILESDA),     ! OUT                
     &       MELT_TILE(TILE_FIELDDA,NTILESDA),     ! OUT                
     &       TILE_FRAC(TILE_FIELDDA,NTILESDA)      ! OUT                
*D ARE1F404.77,ARE1F404.78   
     &       EI_TILE(TILE_FIELDDA,NTILESDA),                            
     &       ESOIL_TILE(TILE_FIELDDA,NTILESDA),                         
     &       FTL_TILE(TILE_FIELDDA,NTILESDA),                           
     &       GS_TILE(TILE_FIELDDA,NTILESDA),                            
     &       LE_TILE(TILE_FIELDDA,NTILESDA),                            
*I ARE1F404.80    
     &       Q1P5M_TILE(TILE_FIELDDA,NTILESDA),                         
     &       T1P5M_TILE(TILE_FIELDDA,NTILESDA),                         
     &       RAD_TILE(TILE_FIELDDA,NTILESDA),                           
*D ARE1F404.83
     &       RIB_TILE(TILE_FIELDDA,NTILESDA)                            
*D ARE1F404.85
     &   RHO_ARESIST_TILE(TILE_FIELDDA,NTILESDA),                       
*D ARE1F404.87
     &   ARESIST_TILE(TILE_FIELDDA,NTILESDA),                           
*D ARE1F404.89
     &   RESIST_B_TILE(TILE_FIELDDA,NTILESDA),                          
*D ABX1F405.200
*D ABX1F405.205
     & ,     PLLTILE(NTILESDA)  ! pseudolevel list for surface tiles
     & ,     L_CTILE            ! Coastal tiling switch for  
!                               ! BL_INTCT    
*I APC5F400.3     
     &    .OR.SF(328,3)     ! Needed for T at 1.5 m over land tiles     
*I APC5F400.4     
     &    .OR.SF(329,3)     ! Needed for Q at 1.5 m over land tiles     
*D AJS1F401.256
        SURF_HT_FLUX_LAND(I) = 0.0                                      
*I ACN1F405.19    

! Set coastal tiling flag
      IF ( H_SECT(3).EQ.'07A' .OR. H_SECT(3).EQ.'08A') THEN
        L_CTILE = .TRUE.
      ENDIF                                                 
 
*D AYY1F404.66
     & L_BL_LSPICE,L_MOM,L_MIXLEN,L_CTILE,                              
*D ARN0F405.26
*I ANG1F405.2     
! OUT Coastal tiling diagnostics :
                                   
     & STASHWORK(SI(339,3,im_index)),
     & STASHWORK(SI(391,3,im_index)),STASHWORK(SI(392,3,im_index)),
     & STASHWORK(SI(393,3,im_index)),STASHWORK(SI(394,3,im_index)),
     & STASHWORK(SI(389,3,im_index)),STASHWORK(SI(390,3,im_index)),
     & STASHWORK(SI(347,3,im_index)),STASHWORK(SI(353,3,im_index)),
     & STASHWORK(SI(343,3,im_index)),STASHWORK(SI(381,3,im_index)),
     & STASHWORK(SI(395,3,im_index)),
                                                                       
*I AJS1F401.307   
     & SURF_HT_FLUX_LAND,SURF_HT_FLUX_SICE,                      
*D ABX1F405.218,ARE1F404.95   
     & A_INTHD(23),                                                     
     & L_PHENOL,L_TRIFFID,L_NEG_TSTAR,NTILES,                           
     & D1(JCANHT_PFT),D1(JCAN_WATER_TYP),D1(JCATCH_TYP),                
*D ARE1F404.97
     & LW_DOWN,SW_TILE,D1(JZ0_TYP),                                     
*I ACN1F405.20    
     & D1(JFRAC_LAND),                           
*D ARE1F404.99
     & OLR,D1(JSNODEP_TYP),D1(JTSTAR_TYP),                              
*I ARE1F404.101   
     & D1(JTSTAR_LAND),D1(JTSTAR_SEA),D1(JTSTAR_SICE),                  
*D ARE1F404.103
     & ECAN_TILE,EI_TILE,ESOIL_TILE,FTL_TILE,                           
     & GS_TILE,LE_TILE,MELT_TILE,RAD_TILE,                              
*D ARE1F404.106,ABX1F405.220  
     & RIB_TILE,Q1P5M_TILE,T1P5M_TILE,TILE_INDEX,TILE_PTS,TILE_FRAC,    
*I AJS1F401.403   
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I AJS1F401.404   
      CALL RES_TILE(                                                    
     &        P_FIELD,LAND_FIELD,LAND_LIST,NTILES,TILE_INDEX            
     &       ,TILE_PTS,.TRUE.,ARESIST,ARESIST_TILE,D1(JCAN_WATER_TYP)   
     &       ,D1(JCATCH_TYP),GS_TILE,RESIST_B_TILE,D1(JSNODEP_TYP)      
     &       ,TILE_FRAC,RESB_SO2,RESS_SO2,RES_FACTOR)                   
*ENDIF                                                                  
!                                                                       
*I AWO3F405.66    
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I AWO3F405.67    
      CALL RES_TILE(                                                    
     &        P_FIELD,LAND_FIELD,LAND_LIST,NTILES,TILE_INDEX            
     &       ,TILE_PTS,.TRUE.,ARESIST,ARESIST_TILE,D1(JCAN_WATER_TYP)   
     &       ,D1(JCATCH_TYP),GS_TILE,RESIST_B_TILE,D1(JSNODEP_TYP)      
     &       ,TILE_FRAC,RESB_NH3,RESS_NH3,RES_FACTOR)                   
*ENDIF                                                                  
!                                                                       
*I AJS1F401.469   
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I AJS1F401.470   
      CALL RES_TILE(                                                    
     &        P_FIELD,LAND_FIELD,LAND_LIST,NTILES,TILE_INDEX            
     &       ,TILE_PTS,.FALSE.,ARESIST,ARESIST_TILE,D1(JCAN_WATER_TYP)  
     &       ,D1(JCATCH_TYP),GS_TILE,RESIST_B_TILE,D1(JSNODEP_TYP)      
     &       ,TILE_FRAC,RESB_SO4_AIT,RESS_SO4_AIT,RES_FACTOR)           
*ENDIF                                                                  
!                                                                       
*I AJS1F401.519   
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I AJS1F401.520   
      CALL RES_TILE(                                                    
     &        P_FIELD,LAND_FIELD,LAND_LIST,NTILES,TILE_INDEX            
     &       ,TILE_PTS,.FALSE.,ARESIST,ARESIST_TILE,D1(JCAN_WATER_TYP)  
     &       ,D1(JCATCH_TYP),GS_TILE,RESIST_B_TILE,D1(JSNODEP_TYP)      
     &       ,TILE_FRAC,RESB_SO4_ACC,RESS_SO4_ACC,RES_FACTOR)           
*ENDIF                                                                  
!                                                                       
*I AJS1F401.569   
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I AJS1F401.570   
      CALL RES_TILE(                                                    
     &        P_FIELD,LAND_FIELD,LAND_LIST,NTILES,TILE_INDEX            
     &       ,TILE_PTS,.FALSE.,ARESIST,ARESIST_TILE,D1(JCAN_WATER_TYP)  
     &       ,D1(JCATCH_TYP),GS_TILE,RESIST_B_TILE,D1(JSNODEP_TYP)      
     &       ,TILE_FRAC,RESB_SO4_DIS,RESS_SO4_DIS,RES_FACTOR)           
*ENDIF                                                                  
!                                                                       
*I AWO3F405.175   
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I AWO3F405.176   
      CALL RES_TILE(                                                    
     &        P_FIELD,LAND_FIELD,LAND_LIST,NTILES,TILE_INDEX            
     &       ,TILE_PTS,.FALSE.,ARESIST,ARESIST_TILE,D1(JCAN_WATER_TYP)  
     &       ,D1(JCATCH_TYP),GS_TILE,RESIST_B_TILE,D1(JSNODEP_TYP)      
     &       ,TILE_FRAC,RESB_FreshSoot,RESS_Soot,RES_FACTOR)            
*ENDIF                                                                  
!                                                                       
*I AWO3F405.235   
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I AWO3F405.236   
      CALL RES_TILE(                                                    
     &        P_FIELD,LAND_FIELD,LAND_LIST,NTILES,TILE_INDEX            
     &       ,TILE_PTS,.FALSE.,ARESIST,ARESIST_TILE,D1(JCAN_WATER_TYP)  
     &       ,D1(JCATCH_TYP),GS_TILE,RESIST_B_TILE,D1(JSNODEP_TYP)      
     &       ,TILE_FRAC,RESB_AgedSoot,RESS_Soot,RES_FACTOR)             
*ENDIF                                                                  
!                                                                       
*I AWO3F405.282   
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I AWO3F405.283   
      CALL RES_TILE(                                                    
     &        P_FIELD,LAND_FIELD,LAND_LIST,NTILES,TILE_INDEX            
     &       ,TILE_PTS,.FALSE.,ARESIST,ARESIST_TILE,D1(JCAN_WATER_TYP)  
     &       ,D1(JCATCH_TYP),GS_TILE,RESIST_B_TILE,D1(JSNODEP_TYP)      
     &       ,TILE_FRAC,RESB_SootInCloud,RESS_Soot,RES_FACTOR)          
*ENDIF                                                                  
!                                                                       
*I APBGF401.39    
          D1(JTSTAR_LAND+TOP_ROW_START+I-2)=0.0                         
          D1(JTSTAR_SEA+TOP_ROW_START+I-2)=0.0                          
          D1(JTSTAR_SICE+TOP_ROW_START+I-2)=0.0                         
*I APBGF401.47    
          D1(JTSTAR_LAND+P_BOT_ROW_START+I-2)=0.0                       
          D1(JTSTAR_SEA+P_BOT_ROW_START+I-2)=0.0                        
          D1(JTSTAR_SICE+P_BOT_ROW_START+I-2)=0.0                       
*I APB2F401.96    
C     
       CALL POLAR(D1(JTSTAR_LAND),D1(JTSTAR_LAND),D1(JTSTAR_LAND),      
*CALL ARGFLDPT                                                          
     &           P_FIELD,P_FIELD,P_FIELD,                               
     &           TOP_ROW_START+ROW_LENGTH,                              
     &           P_BOT_ROW_START-ROW_LENGTH,                            
     &           ROW_LENGTH,1)
C     
       CALL POLAR(D1(JTSTAR_SEA),D1(JTSTAR_SEA),D1(JTSTAR_SEA),      
*CALL ARGFLDPT                                                          
     &           P_FIELD,P_FIELD,P_FIELD,                               
     &           TOP_ROW_START+ROW_LENGTH,                              
     &           P_BOT_ROW_START-ROW_LENGTH,                            
     &           ROW_LENGTH,1)
C     
      CALL POLAR(D1(JTSTAR_SICE),D1(JTSTAR_SICE),D1(JTSTAR_SICE),       
*CALL ARGFLDPT                                                          
     &           P_FIELD,P_FIELD,P_FIELD,                               
     &           TOP_ROW_START+ROW_LENGTH,                              
     &           P_BOT_ROW_START-ROW_LENGTH,                            
     &           ROW_LENGTH,1)
*I APBGF401.57    
          D1(JTSTAR_LAND+TOP_ROW_START+I-2)=                            
     &      D1(JTSTAR_LAND+TOP_ROW_START+ROW_LENGTH+I-2)                
          D1(JTSTAR_SEA+TOP_ROW_START+I-2)=                             
     &      D1(JTSTAR_SEA+TOP_ROW_START+ROW_LENGTH+I-2)                 
          D1(JTSTAR_SICE+TOP_ROW_START+I-2)=                            
     &      D1(JTSTAR_SICE+TOP_ROW_START+ROW_LENGTH+I-2)                
*I APBGF401.66    
          D1(JTSTAR_LAND+P_BOT_ROW_START+I-2)=                          
     &      D1(JTSTAR_LAND+P_BOT_ROW_START-ROW_LENGTH+I-2)              
          D1(JTSTAR_SEA+P_BOT_ROW_START+I-2)=                           
     &      D1(JTSTAR_SEA+P_BOT_ROW_START-ROW_LENGTH+I-2)               
          D1(JTSTAR_SICE+P_BOT_ROW_START+I-2)=                          
     &      D1(JTSTAR_SICE+P_BOT_ROW_START-ROW_LENGTH+I-2)              
*I BL_CTL1.298   
! Coastal tiling diagnostics

C Item 337:
C Land heat flux from surface to level 1 (land mean) (W/m2)

      IF (SF(337,3)) THEN                                               
                                                                        
        CALL COPYDIAG(STASHWORK(SI(337,3,im_index)),SURF_HT_FLUX_LAND,  
     &       FIRST_POINT,LAST_POINT,P_FIELD,ROW_LENGTH,                 
     &       im_ident,3,337,                                            
*CALL ARGPPX                                                            
     &       ICODE,CMESSAGE)                                            
                                                                        
        IF (ICODE .GT. 0) GOTO 9999                                     
                                                                        
      END IF                                                            

C Item 338:
C Net surface sea-ice heat flux (sea mean) (W/m2) 

      IF (SF(338,3)) THEN                                               
                                                                        
        CALL COPYDIAG(STASHWORK(SI(338,3,im_index)),SURF_HT_FLUX_SICE,  
     &       FIRST_POINT,LAST_POINT,P_FIELD,ROW_LENGTH,                 
     &       im_ident,3,338,                                            
*CALL ARGPPX                                                            
     &       ICODE,CMESSAGE)                                            
                                                                        
        IF (ICODE .GT. 0) GOTO 9999                                     
                                                                        
      END IF                                                            

! End of coastal tiling diagsnostics
                                                                        
*D ARE1F405.18,AJS1F401.704  
*D ABX1F405.221
CL ITEM 287: CANOPY EVAPORATION ON TILES                                
*D ABX1F405.224
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
*D ABX1F405.226
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.232,ABX1F405.233  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*D ABX1F405.244
CL ITEM 288: TRANSPIRATION + SOIL EVAPORATION ON TILES                  
*D ABX1F405.247
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
*D ABX1F405.249
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.255,ABX1F405.256  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*D ABX1F405.293
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
*D ABX1F405.295
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.301,ABX1F405.302  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*D ABX1F405.368
CL ITEM 294: BULK RICHARDSON NUMBER ON TILES                            
*D ABX1F405.371
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
*D ABX1F405.373
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.379,ABX1F405.380  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*I ABX1F405.389   
CL ITEM 314: SURFACE NET RADIATION ON TILES                             
*D ABX1F405.391,ABX1F405.405  
      IF (SF(314,3)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,314,3,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.411,ABX1F405.412  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*D ABX1F405.415,ABX1F405.416  
     &          STASHWORK(SI(314,3,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),RAD_TILE(1,PSLEVEL),                         
*D ABX1F405.422
*D ABX1F405.426
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
*D ABX1F405.428
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.434,ABX1F405.435  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*D ABX1F405.524
CL ITEM 321: CANOPY WATER CONTENT ON TILES                              
*D ABX1F405.527
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
*D ABX1F405.529
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.535,ABX1F405.536  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*D ABX1F405.540
     &           *P_FIELD),D1(JCAN_WATER_TYP+((PSLEVEL-1)*LAND_FIELD)), 
*D ABX1F405.547
CL ITEM 322: CANOPY CAPACITY ON TILES                                   
*D ABX1F405.550
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
*D ABX1F405.552
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.558,ABX1F405.559  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*D ABX1F405.563
     &           *P_FIELD),D1(JCATCH_TYP+((PSLEVEL-1)*LAND_FIELD)),     
*D ABX1F405.570,ABX1F405.578  
*D ABX1F405.582
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
*D ABX1F405.584
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
*D ABX1F405.590,ABX1F405.591  
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
*I ABX1F405.623   
                                                                        
CL ITEM 328: TEMPERATURE AT 1.5M OVER TILES                             
                                                                        
      IF (SF(328,3)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,328,3,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)                                            
        IF (ICODE.GT.0) THEN                                            
          RETURN                                                        
        END IF                                                          
        PSLEVEL_OUT=0                                                   
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1                                   
            CALL FROM_LAND_POINTS (                                     
     &          STASHWORK(SI(328,3,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),T1P5M_TILE(1,PSLEVEL),                       
     &           D1(JLAND),P_FIELD,LAND_FIELD)                          
          END IF                                                        
        END DO                                                          
      END IF                                                            
                                                                        
                                                                        
CL ITEM 329: SPECIFIC HUMIDITY AT 1.5M OVER TILES                       
                                                                        
      IF (SF(329,3)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,329,3,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)                                            
        IF (ICODE.GT.0) THEN                                            
          RETURN                                                        
        END IF                                                          
        PSLEVEL_OUT=0                                                   
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1                                   
            CALL FROM_LAND_POINTS (                                     
     &          STASHWORK(SI(329,3,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),Q1P5M_TILE(1,PSLEVEL),                       
     &           D1(JLAND),P_FIELD,LAND_FIELD)                          
          END IF                                                        
        END DO                                                          
      END IF                                                            
                                                                        
                                                                        
CL ITEM 330: SURFACE LATENT HEAT FLUX ON TILES                          
                                                                        
      IF (SF(330,3)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,330,3,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)                                            
        IF (ICODE.GT.0) THEN                                            
          RETURN                                                        
        END IF                                                          
        PSLEVEL_OUT=0                                                   
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1                                   
            CALL FROM_LAND_POINTS (                                     
     &          STASHWORK(SI(330,3,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),LE_TILE(1,PSLEVEL),                          
     &           D1(JLAND),P_FIELD,LAND_FIELD)                          
          END IF                                                        
        END DO                                                          
      END IF                                                            
                                                                        
                                                                        
CL ITEM 331: SUBLIMATION ON TILES                                       
                                                                        
      IF (SF(331,3)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,331,3,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)                                            
        IF (ICODE.GT.0) THEN                                            
          RETURN                                                        
        END IF                                                          
        PSLEVEL_OUT=0                                                   
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1                                   
            CALL FROM_LAND_POINTS (                                     
     &          STASHWORK(SI(331,3,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),EI_TILE(1,PSLEVEL),                          
     &           D1(JLAND),P_FIELD,LAND_FIELD)                          
          END IF                                                        
        END DO                                                          
      END IF                                                            
                                                                        
                                                                        
CL ITEM 332: TOA outward LW radiation after boundary layer              
        IF (SF(332,3)) THEN                                             
        CALL COPYDIAG(STASHWORK(SI(332,3,im_index)),OLR,                
     &       FIRST_POINT,LAST_POINT,P_FIELD,ROW_LENGTH,                 
     &       im_ident,3,332,                                            
*CALL ARGPPX                                                            
     &       ICODE,CMESSAGE)                                            
                                                                        
        IF (ICODE .GT. 0) GOTO 9999                                     
        END IF                                                          
*DECLARE BL_IC7A
*D BL_IC7A.2
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*D BL_IC7A.76,BL_IC7A.77   
     & CO2_MMR,PHOTOSYNTH_ACT_RAD,PSTAR,RAD_SICE,                       
     & TIMESTEP,L_RMBL,L_BL_LSPICE,L_MOM,L_MIXLEN,L_CTILE,              
*I BL_IC7A.86    
! OUT Coastal tiling diagnostics :              
     & RIB_SSI,TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,VSHR_LAND,VSHR_SSI,
     & E_SSI,EI_SICE,FTL_SSI,RADNET_SICE,FLANDG,                        
                                                                        
*D BL_IC7A.102
     & SURF_HT_FLUX,SURF_HT_FLUX_LAND,SURF_HT_FLUX_SICE,
     & ZH,T1_SD,Q1_SD,ERROR,                               
*D ABX1F405.754
     & ASTEPS_SINCE_TRIFFID,                                            
     & L_PHENOL,L_TRIFFID,L_NEG_TSTAR,NTILES,                           
*D BL_IC7A.109
     & FRAC,LW_DOWN,SW_TILE,Z0V_TILE,                                   
*I ACN1F405.114   
     & FLAND,                          
*D BL_IC7A.111
     & OLR,SNOW_TILE,TSTAR_TILE,                                        
*I BL_IC7A.112   
     & TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,                  
*D BL_IC7A.114,ABX1F405.755  
     & ECAN_TILE,EI_TILE,ESOIL_TILE,FTL_TILE,GS_TILE,LE_TILE,MELT_TILE, 
     & RADNET_TILE,G_LEAF,GPP_FT,NPP_FT,RESP_P_FT,RESP_S,RESP_W_FT,     
*D BL_IC7A.117,ABX1F405.756  
     & RIB_TILE,Q1P5M_TILE,T1P5M_TILE,TILE_INDEX,TILE_PTS,TILE_FRAC,
*I BL_IC7A.186   
     &,ASTEPS_SINCE_TRIFFID      ! IN Number of atmospheric             
!                                !    timesteps since last call         
!                                !    to TRIFFID.                       
*I ACN1F405.116   
     &,NTILES                    ! IN No. of land-surface tiles         
*D BL_IC7A.197,BL_IC7A.198  
     &,CANOPY_TILE(LAND_FIELD,NTILES)                                   
!                                ! IN Surface/canopy water for          
*D BL_IC7A.200
     &,CATCH_TILE(LAND_FIELD,NTILES)                                    
*D BL_IC7A.202
!                                !    land tiles (kg per sq m)          
*I BL_IC7A.203   
     &,FRAC(LAND_FIELD,NTYPE)    ! IN Fractions of surface types.       
     &,FLAND(LAND_FIELD)         ! IN Land fraction on land tiles.    
     &,FLANDG(P_FIELD)           ! IN Land fraction on all tiles.     
*I BL_IC7A.210   
     &,LW_DOWN(P_FIELD)          ! IN Surface downward LW radiation     
!                                !    (W/m2).                           
*D BL_IC7A.227,BL_IC7A.229  
     &,SW_TILE(LAND_FIELD,NTILES)! IN Surface net SW radiation on land  
!                                !    tiles (W/m2).                     
*D BL_IC7A.233
     &,Z0V_TILE(LAND_FIELD,NTILES)!IN Tile roughness lengths (m).       
*D BL_IC7A.266,BL_IC7A.271  
     &,RAD_SICE(P_FIELD)         ! IN Surface net SW and downward LW    
!                                !    radiation on sea-ice fraction     
!                                !    (W/sq m, positive downwards).     
*I ABX1F405.767   
     &,L_CTILE                  ! IN true if model to be run with 
!                               !    coastal tiling       
*I BL_IC7A.308   
     &,OLR(P_FIELD)              ! IN    TOA - surface upward LW on     
!                                !       last radiation timestep        
!                                ! OUT   Corrected TOA outward LW       
     &,SNOW_TILE(LAND_FIELD,NTILES)                                     
!                                ! INOUT Snow on tiles (kg/m2).         
*I BL_IC7A.312   
     &,TSTAR_LAND(P_FIELD)       ! OUT   Land mean sfc temperature (K)
     &,TSTAR_SEA(P_FIELD)        ! IN    Open sea sfc temperature (K).
     &,TSTAR_SICE(P_FIELD)       ! INOUT Sea-ice sfc temperature (K). 
*D BL_IC7A.314
     &,TSTAR_TILE(LAND_FIELD,NTILES)                                    
*D BL_IC7A.357,BL_IC7A.358  
     &,EI_TILE(LAND_FIELD,NTILES)!OUT EI for land tiles                 
     &,ESOIL_TILE(LAND_FIELD,NTILES)                                    
                                ! OUT ES for land tiles                 
*D BL_IC7A.367
     &,FTL_TILE(LAND_FIELD,NTILES)                                      
*I BL_IC7A.368   
     &,E_SSI(P_FIELD)           ! OUT   Surface FQW for mean sea.       
     &,EI_SICE(P_FIELD)         ! OUT   Sea-ice sumblimation            
!                               !       (sea mean).                     
     &,FTL_SSI(P_FIELD)         ! OUT sea mean surface heat flux        
*I BL_IC7A.371   
     &,GS_TILE(LAND_FIELD,NTILES)!OUT Surface conductance for           
!                               !     land tiles.                       
*I BL_IC7A.373   
     &,LE_TILE(LAND_FIELD,NTILES)!OUT Surface latent heat flux for      
!                               !     land tiles (W/m2).                
     &,MELT_TILE(LAND_FIELD,NTILES)                                     
!                               ! OUT Snowmelt on tiles (kg/m2/s).      
*I BL_IC7A.377   
     &,RADNET_TILE(LAND_FIELD,NTILES)                                   
!                               ! OUT Tile surface net radiation.       
*D BL_IC7A.389
     &,RIB_TILE(LAND_FIELD,NTILES)                                      
!                               ! OUT RIB for land tiles.               
     &,RIB_LAND(P_FIELD)        !     Land mean bulk Richardson no.  
!                                        for lowest layer.              
     &,RIB_SSI(P_FIELD)         ! OUT Sea mean bulk Richardson no.   
!                                        for lowest layer.              
*I BL_IC7A.395   
     &,SURF_HT_FLUX_LAND(P_FIELD)                               
!                               ! OUT Net downward heat flux at         
!                               !     surface over land                 
!                               !     fraction of gridbox (W/m2).       
     &,SURF_HT_FLUX_SICE(P_FIELD)                               
!                               ! OUT Net downward heat flux at         
!                               !     surface over sea-ice              
!                               !     fraction of gridbox (W/m2).       
*I BL_IC7A.400   
     &,TAUX_LAND(U_FIELD)       ! OUT W'ly compt of land sfc wind    
!                               !     stress (N/sq m). (On UV-grid    
!                               !     with first and last rows       
!                               !     undefined or, at present,      
!                               !     set to missing data            
     &,TAUX_SSI(U_FIELD)        ! OUT W'ly compt of mean sea sfc wind
!                               !     stress (N/sq m). (On UV-grid    
!                               !     with first and last rows       
!                               !     undefined or, at present,      
!                               !     set to missing data            
*D ABX1F405.784,ABX1F405.787  
     &,TAUY_LAND(U_FIELD)       ! OUT S'ly compt of land sfc wind    
!                               !     stress (N/sq m).  On UV-grid;   
!                               !     comments as per TAUX.          
     &,TAUY_SSI(U_FIELD)        ! OUT S'ly compt of mean sea sfc wind
!                               !     stress (N/sq m).  On UV-grid;   
!                               !     comments as per TAUX.          
     &,TILE_FRAC(LAND_FIELD,NTILES)                                     
!                               ! OUT Tile fractions. Equal to surface  
!                               !     type fractions if NTILES=NTYPE,   
!                               !     or 1 if NTILES=1.                 
*I BL_IC7A.405   
     &,VSHR_LAND(P_FIELD)       ! OUT Magnitude of land sfc-to-lowest
!                                     atm level wind shear (m per s).   
     &,VSHR_SSI(P_FIELD)        ! OUT Mag. of mean sea sfc-to-lowest 
!                                     atm level wind shear (m per s).   
*D BL_IC7A.413
     &,RHO_ARESIST_TILE(LAND_FIELD,NTILES)                              
*D BL_IC7A.415
     &,ARESIST_TILE(LAND_FIELD,NTILES)                                  
*D BL_IC7A.417
     &,RESIST_B_TILE(LAND_FIELD,NTILES)                                 
*I BL_IC7A.418   
     &,RADNET_SICE(P_FIELD)     ! OUT Sea-ice surface net radiation.  
*I BL_IC7A.436   
     &,Q1P5M_TILE(LAND_FIELD,NTILES)                                    
!                               ! OUT Q1P5M over land tiles.            
*I BL_IC7A.437   
     &,T1P5M_TILE(LAND_FIELD,NTILES)                                    
!                               ! OUT T1P5M over land tiles.            
*D BL_IC7A.444,BL_IC7A.445  
     & ECAN_TILE(LAND_FIELD,NTILES)                                     
                      ! OUT ECAN for land tiles                         
*D BL_IC7A.455,BL_IC7A.460  
*D BL_IC7A.480,BL_IC7A.482  
     & CANHC_TILE(LAND_FIELD,NTILES)                                    
!                               ! LOCAL Areal heat capacity of canopy   
!                               !       for land tiles (J/K/m2).        
     &,VFRAC_TILE(LAND_FIELD,NTILES)                                    
!                               ! LOCAL Fractional canopy coverage for  
!                               !       land tiles.                     
     &,WT_EXT_TILE(LAND_FIELD,SM_LEVELS,NTILES)                         
*D BL_IC7A.485
!                               !       soil layer by each tile.        
*D BL_IC7A.491,BL_IC7A.494  
                                                                        
      REAL                                                              
     & ALPHA1(LAND_FIELD,NTILES)! LOCAL Mean gradient of saturated      
!                               !       specific humidity with respect  
!                               !       to temperature between the      
!                               !       bottom model layer and tile     
!                               !       surfaces                        
     &,ALPHA1_SICE(P_FIELD)     ! LOCAL ALPHA1 for sea-ice.             
     &,ASHTF(P_FIELD)           ! LOCAL Coefficient to calculate        
!                               !       surface heat flux into soil or  
!                               !       sea-ice.                        
     &,ASHTF_TILE(LAND_FIELD,NTILES)                                    
!                               ! LOCAL Coefficient to calculate        
!                               !       surface heat flux into land     
!                               !       tiles.                          
     &,DTRDZ(P_FIELD,BL_LEVELS) ! LOCAL -g.dt/dp for model layers.      
     &,FLAKE(LAND_FIELD,NTILES) ! LOCAL Lake fraction.                  
     &,FQW_TILE(LAND_FIELD,NTILES)                                      
!                               ! LOCAL Surface FQW for land tiles      
     &,FQW_ICE(P_FIELD)         ! LOCAL Surface FQW for sea-ice         
     &,FTL_ICE(P_FIELD)         ! LOCAL Surface FTL for sea-ice         
     &,QW(P_FIELD,BL_LEVELS)    ! LOCAL Total water content             
     &,TL(P_FIELD,BL_LEVELS)    ! LOCAL Ice/liquid water temperature    
     &,TSTAR_TILE_OLD(LAND_FIELD,NTILES)                                
!                               ! LOCAL Tile surface temperatures at    
!                               !       beginning of timestep.          
     &,FRACA(LAND_FIELD,NTILES) ! LOCAL Fraction of surface moisture    
!                               !       flux with only aerodynamic      
!                               !       resistance for snow-free land   
!                               !       tiles.                          
     &,RESFS(LAND_FIELD,NTILES) ! LOCAL Combined soil, stomatal         
!                               !       and aerodynamic resistance      
!                               !       factor for fraction (1-FRACA)   
!                               !       of snow-free land tiles.        
     &,RESFT(LAND_FIELD,NTILES) ! LOCAL Total resistance factor.        
!                               !       FRACA+(1-FRACA)*RESFS for       
!                               !       snow-free land, 1 for snow.     
     &,RHOKH_TILE(LAND_FIELD,NTILES)                                    
!                               ! LOCAL Surface exchange coefficients   
!                               !       for land tiles                  
     &,RHOKH_SICE(P_FIELD)      ! LOCAL Surface exchange coefficients   
!                               !       for sea and sea-ice             
     &,RHOKPM(LAND_FIELD,NTILES)! LOCAL Land surface exchange coeff.    
     &,RHOKPM_SICE(P_FIELD)     ! LOCAL Sea-ice surface exchange coeff. 
     &,H_BLEND_OROG(P_FIELD)    ! LOCAL Blending height used as part of 
!                               !       effective roughness scheme      
     &,Z0H(P_FIELD)             ! LOCAL Roughness length for heat and   
!                               !       moisture (m).                   
     &,Z0H_TILE(LAND_FIELD,NTILES)                                      
!                               ! LOCAL Tile roughness lengths for heat 
!                               !       and moisture (m).               
     &,Z0M(P_FIELD)             ! LOCAL Roughness length for            
!                               !       momentum (m).                   
     &,Z0M_TILE(LAND_FIELD,NTILES)                                      
!                               ! LOCAL Tile roughness lengths for      
!                               !       momentum.                       
     &,Z0M_EFF(P_FIELD)         ! LOCAL Effective grid-box roughness    
!                               !       length for momentum             
     &,CDR10M_UV(U_FIELD)       ! LOCAL Ratio of CD's reqd for          
!                               !       calculation of 10 m wind. On    
!                               !       UV-grid; comments as per RHOKM. 
     &,CHR1P5M(LAND_FIELD,NTILES)!LOCAL Ratio of coefffs for            
!                               !       calculation of 1.5m temp for    
!                               !       land tiles.                     
     &,CHR1P5M_SICE(P_FIELD)    ! LOCAL CHR1P5M for sea and sea-ice     
!                               !       (leads ignored).                
     &,CT_CTQ(P_FIELD,BL_LEVELS)! LOCAL Coefficient in T and q          
!                                       tri-diagonal implicit matrix    
     &,CQ_CM(U_FIELD,BL_LEVELS) ! LOCAL Coefficient in U and V          
!                                       tri-diagonal implicit matrix    
     &,DQW(P_FIELD,BL_LEVELS)   ! LOCAL BL increment to q field         
     &,DTL(P_FIELD,BL_LEVELS)   ! LOCAL BL increment to T field         
     &,DU(U_FIELD,BL_LEVELS)    ! LOCAL BL increment to u wind field    
     &,DV(U_FIELD,BL_LEVELS)    ! LOCAL BL increment to v wind field    
     &,RDZ(P_FIELD,BL_LEVELS)   ! LOCAL RDZ(,1) is the reciprocal of    
!                               !       the height of level 1, i.e. of  
!                               !       the middle of layer 1.  For     
!                               !       K > 1, RDZ(,K) is the           
!                               !       reciprocal of the vertical      
!                               !       distance from level K-1 to      
!                               !       level K.                        
     &,RDZUV(U_FIELD,BL_LEVELS) ! LOCAL RDZ (K > 1) on UV-grid.         
!                               !       Comments as per RHOKM (RDZUV).  
     &,FB_SURF(P_FIELD)         ! Surface flux buoyancy over density    
!                               ! (m^2/s^3)                             
!                                                                       
     &,U_S(P_FIELD)             ! Surface friction velocity (m/s)       
     &,TV1_SD(P_FIELD)          ! Standard deviation of turbulent       
!                               ! fluctuations of surface layer         
!                               ! virtual temperature (K).              
     &,E_LAND(P_FIELD)          ! LOCAL FQW over mean land      
     &,EI_LAND(P_FIELD)         ! LOCAL EI over mean land       
     &,FTL_LAND(P_FIELD)        ! LOCAL FTL over mean land      
     &,FLANDG_UV(U_FIELD)       ! Land frac (on UV-grid, with 1st 
!                               ! and last rows undefined or, at 
!                               ! present, set to "missing data")
     &,TSTAR_SSI(P_FIELD)       ! LOCAL Sea mean sfc temperature (K).
                                                                        
*I BL_IC7A.497   
     &,J                        ! LOCAL Tile point index                
*D BL_IC7A.598
*I BL_IC7A.632   
! ----------------------------------------------------------------------
! Set Coastal tiling dependent prognostics:                             
! ----------------------------------------------------------------------
      IF(L_CTILE)THEN                                                   
                                                                        
        DO I=P1,P1+P_POINTS-1
          FLANDG(I)=0.0                                          
          IF(ICE_FRACT(I).LE.0.0)THEN                               
            TSTAR_SSI(I)=TSTAR_SEA(I)                             
          ELSE                                                        
            TSTAR_SSI(I)=ICE_FRACT(I)*TSTAR_SICE(I)             
     &        +(1.0-ICE_FRACT(I))*TSTAR_SEA(I)                    
          ENDIF                                                       
        ENDDO           
        DO L=1,LAND1+LAND_PTS-1                                         
          I = LAND_INDEX(L)                                             
          FLANDG(I)=FLAND(L)                                            
        ENDDO                                                           
                                                                        
      ENDIF                                                             
                                                                        
*D BL_IC7A.634
! Call TILEPTS to calculate TILE_PTS and TILE_INDEX for surface types   
*D BL_IC7A.650
     & N_ROWS,FIRST_ROW,ROW_LENGTH,FLANDG,                              
*D BL_IC7A.653
     & VSHR,VSHR_LAND,VSHR_SSI,Z1                                       
*D BL_IC7A.662
     & P_FIELD,SM_LEVELS,NTILES,TILE_PTS,TILE_INDEX,                    
*I BL_IC7A.665   
     & CANHC_TILE,VFRAC_TILE,FLAKE,                                     
*D ABX1F405.802
     & RESP_P,RESP_P_FT,RESP_S,RESP_W_FT,SMC,WT_EXT_TILE                
*I BL_IC7A.668   
                                                                        
!---------------------------------------------------------------------- 
! If TRIFFID is being used apply any correction to the land-atmosphere  
! fluxes on the first timestep after the last TRIFFID call. Such a      
! correction will typically be associated with a total depletion of     
! carbon or with maintanence of the seed fraction. The corrections      
! are stored in the accumulation variables after the call to TRIFFID.   
! The correction is added to the instantaneous land-atmosphere fluxes   
! (so that the atmospheric carbon budget is corrected) but is not       
! included in the accumulation variables which drive TRIFFID, since     
! this has already been dealt with during the last TRIFFID call.        
!---------------------------------------------------------------------- 
      IF (L_TRIFFID.AND.(ASTEPS_SINCE_TRIFFID.EQ.1)) THEN               
        DO N=1,NPFT                                                     
          DO L=LAND1,LAND1+LAND_PTS-1                                   
            NPP_FT(L,N)=NPP_FT(L,N)+NPP_FT_ACC(L,N)/TIMESTEP            

            if (n.eq.1) NPP(L)=0.
            NPP(L)=NPP(L)+FRAC(L,N)*NPP_FT(L,N)

            RESP_P_FT(L,N)=RESP_P_FT(L,N)-NPP_FT_ACC(L,N)/TIMESTEP      
            NPP_FT_ACC(L,N)=-NPP_FT_ACC(L,N)                            
          ENDDO                                                         
        ENDDO                                                           
        DO L=LAND1,LAND1+LAND_PTS-1                                     
          RESP_S(L)=RESP_S(L)+RESP_S_ACC(L)/TIMESTEP                    
          RESP_S_ACC(L)=-RESP_S_ACC(L)                                  
        ENDDO                                                           
      ENDIF                                                             
*D BL_IC7A.685
*D BL_IC7A.687
! Reset TILE_PTS and TILE_INDEX and set tile fractions to 1 if aggregate
! tiles are used (NTILES=1).                                            
! Otherwise, set tile fractions to surface type fractions.              
*D BL_IC7A.689
      DO N=1,NTILES                                                     
*D BL_IC7A.691
          TILE_FRAC(L,N) = 0.                                           
*D BL_IC7A.695,BL_IC7A.703  
      IF (NTILES.EQ.1) THEN                                             
        TILE_PTS(1) = LAND_PTS                                          
        DO L=LAND1,LAND1+LAND_PTS-1                                     
          TILE_FRAC(L,1) = 1.                                           
          TILE_INDEX(L+1-LAND1,1) = L                                   
        ENDDO                                                           
      ELSE                                                              
        DO N=1,NTYPE                                                    
          DO J=1,TILE_PTS(N)                                            
            L = TILE_INDEX(J,N)                                         
            TILE_FRAC(L,N) = FRAC(L,N)                                  
          ENDDO                                                         
        ENDDO                                                           
      ENDIF                                                             
*D BL_IC7A.706,BL_IC7A.713  
! Call boundary layer routine                                           
*D BL_IC7A.716
      CALL SF_EXPL (                                                    
*D BL_IC7A.719,BL_IC7A.720  
     & P_FIELD,U_FIELD,LAND_FIELD,ROW_LENGTH,                           
*D BL_IC7A.724
     & AK(1),BK(1),AKH(1),BKH(1),DELTA_AK(1),DELTA_BK(1),               
*D BL_IC7A.732,BL_IC7A.735  
     & NTILES,TILE_INDEX,TILE_PTS,SM_LEVELS,                            
     & CANHC_TILE,CANOPY_TILE,CATCH_TILE,FLAKE,GS_TILE,HCON,            
     & HO2R2_OROG,FLAND,FLANDG,
     & SNOW_TILE,SIL_OROG_LAND,SMVCST,STHF,STHU,             
     & TILE_FRAC,VFRAC_TILE,Z0V_TILE,                                   
*D BL_IC7A.738
     & ICE_FRACT,U_0,V_0,                                               
*D BL_IC7A.741
     & CF(1,1),QCF(1,1),QCL(1,1),                                       
*D BL_IC7A.744,ABX1F405.835  
     & PSTAR,LW_DOWN,RAD_SICE,SW_TILE,TIMESTEP,
     & VSHR,VSHR_LAND,VSHR_SSI,ZH,                 
     & Q(1,1),T(1,1),T_DEEP_SOIL,TI,                                    
     & TSTAR,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,
     & TSTAR_TILE,U(1,1),V(1,1),                                  
     & L_BL_LSPICE,                                                     
*D BL_IC7A.748
     & SFME,SQ1P5,ST1P5,SU10,SV10,                                      
*D BL_IC7A.751,BL_IC7A.752  
     & Z0MSEA,                                                          
*D BL_IC7A.755,BL_IC7A.757  
     & CD,CH,E_SEA,QW(1,1),TL(1,1),FQW(1,1),                            
     & FTL(1,1),FTL_TILE,LE_TILE,H_SEA,RADNET_SICE,RADNET_TILE,         
     & RHOKM(1,1),RIB,RIB_TILE,TAUX(1,1),TAUY(1,1),                     
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           
*D BL_IC7A.760,BL_IC7A.761  
     & FME,                                                             
*D BL_IC7A.768,BL_IC7A.769  
! OUT data required for 4D-VAR :                                        
     & RHO_CD_MODV1,                                                    
*D BL_IC7A.772,BL_IC7A.774  
     & FB_SURF,U_S,T1_SD,Q1_SD,TV1_SD,                                  
                                                                        
! OUT data required elsewhere in boundary layer or surface code         
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,DTRDZ(1,1),FQW_TILE,         
     & FQW_ICE,FTL_ICE,TSTAR_TILE_OLD,FRACA,RESFS,RESFT,                
     & RHOKH(1,1),RHOKH_TILE,RHOKH_SICE,RHOKPM,RHOKPM_SICE,             
     & Z1,H_BLEND_OROG,Z0H,Z0H_TILE,Z0M,Z0M_TILE,Z0M_EFF,               
     & CDR10M_UV,CHR1P5M,CHR1P5M_SICE,                                  
     & FLANDG_UV,                                               
*I BL_IC7A.778   
                                                                        
      CALL BDY_EXPL (                                                   
                                                                        
! IN values defining field dimensions and subset to be processed :      
     & P_FIELD,U_FIELD,ROW_LENGTH,                                      
     & N_P_ROWS,N_U_ROWS,P_POINTS,P1,U_POINTS,U1,                       
                                                                        
! IN values defining vertical grid of model atmosphere :                
     & BL_LEVELS,P_LEVELS,AK,BK,AKH,BKH,DELTA_AK,DELTA_BK,              
     & EXNER,                                                           
                                                                        
! IN sea/sea-ice data :                                                 
     & U_0,V_0,                                                         
                                                                        
! IN cloud data :                                                       
     & CF,QCF,QCL,CCA,CCB,CCT,                                          
                                                                        
! IN everything not covered so far :                                    
     & PSTAR,RAD_HR,RADHR_DIM1,                                         
     & FB_SURF,U_S,T1_SD,Q1_SD,TV1_SD,                                  
     & H_BLEND_OROG,Z0M_EFF,                                            
     & TIMESTEP,L_BL_LSPICE,L_MOM,                                      
                                                                        
! INOUT data :                                                          
     & Q,T,U,V,ZH,                                                      
                                                                        
! OUT Diagnostic not requiring STASH flags :                            
     & QW,TL,FQW,FTL,                                                   
     & RHOKH,RHOKM,                                                     
     & TAUX,TAUY,ZHT,                                                   
     & BL_TYPE_1,BL_TYPE_2,BL_TYPE_3,BL_TYPE_4,BL_TYPE_5,BL_TYPE_6,     
                                                                        
! OUT data required for tracer mixing :                                 
     & NRML,                                                            
                                                                        
! OUT data required for 4D-VAR :                                        
     & RHO_KM,                                                          
                                                                        
! OUT data required elsewhere in UM system :                            
     & DTRDZ,RDZ,RDZUV,                                                 
     & DU,DV,CT_CTQ,DQW,DTL,CQ_CM,                                      
                                                                        
! LOGICAL LTIMER                                                        
     & LTIMER                                                           
     & )                                                                
                                                                        
                                                                        
      CALL SF_IMPL (                                                    
                                                                        
! IN values defining field dimensions and subset to be processed :      
     & P_FIELD,U_FIELD,LAND_FIELD,ROW_LENGTH,                           
     & P_POINTS,P1,LAND1,LAND_PTS,U_POINTS,U1,                          
                                                                        
! IN soil/vegetation/land surface data :                                
     & LAND_INDEX,LAND_MASK,                                            
     & NTILES,TILE_INDEX,TILE_PTS,SM_LEVELS,                            
     & CANHC_TILE,CANOPY_TILE,FLAKE,SMC,                                
     & TILE_FRAC,WT_EXT_TILE,                                           
     & FLAND,FLANDG,                                                    
                                                                        
! IN sea/sea-ice data :                                                 
     & DI,ICE_FRACT,U_0,V_0,                                            
                                                                        
! IN everything not covered so far :                                    
     & PSTAR,LW_DOWN,RAD_SICE,SW_TILE,TIMESTEP,                         
     & T_DEEP_SOIL,QW(1,1),TL(1,1),U(1,1),V(1,1),RHOKM(1,1),            
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,                             
     & DTRDZ(1,1),DU(1,1),DV(1,1),FQW_TILE,FQW_ICE,FTL_ICE,             
     & TSTAR_TILE_OLD,                                                  
     & FRACA,RESFS,RESFT,RHOKH_TILE,RHOKH_SICE,RHOKPM,RHOKPM_SICE,      
     & Z1,Z0H,Z0H_TILE,Z0M,Z0M_TILE,CDR10M_UV,CHR1P5M,CHR1P5M_SICE,     
     & CT_CTQ(1,1),DQW(1,1),DTL(1,1),CQ_CM(1,1),                        
     & L_NEG_TSTAR,                                                     
     & FLANDG_UV,                                               
                                                                        
! IN STASH flags :-                                                     
     & SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,                            
                                                                        
! INOUT data :                                                          
     & TI,TSTAR,
     & TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,                       
     & TSTAR_TILE,SNOW_TILE,                                   
     & LE_TILE,RADNET_SICE,RADNET_TILE,                                 
     & E_SEA,FQW(1,1),FTL(1,1),FTL_TILE,H_SEA,OLR,TAUX(1,1),TAUY(1,1),  
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           
                                                                        
! OUT Diagnostic not requiring STASH flags :                            
     & ECAN,EI_TILE,ESOIL_TILE,                                         
     & SEA_ICE_HTF,SURF_HT_FLUX,SURF_HT_FLUX_LAND,SURF_HT_FLUX_SICE,    
                                                                        
! OUT diagnostic requiring STASH flags :                                
     & SICE_MLT_HTF,SNOMLT_SURF_HTF,LATENT_HEAT,                        
     & Q1P5M,Q1P5M_TILE,T1P5M,T1P5M_TILE,U10M,V10M,                     
                                                                        
! OUT data required elsewhere in UM system :                            
     & ECAN_TILE,EI,ES,EXT,SNOWMELT,MELT_TILE,                          
     & ERROR,                                                           
                                                                        
! LOGICAL LTIMER                                                        
     & LTIMER                                                           
     & )                                                                
                                                                        
                                                                        
      CALL BDY_IMPL (                                                   
                                                                        
! IN values defining field dimensions and subset to be processed :      
     & P_FIELD,P1,U_FIELD,U1,P_POINTS,U_POINTS,ROW_LENGTH,              
                                                                        
! IN values defining vertical grid of model atmosphere :                
     & BL_LEVELS,                                                       
                                                                        
! IN data :                                                             
     & RHOKH,RHOKM,RDZ,RDZUV,                                           
                                                                        
! INOUT data :                                                          
     & Q,T,U,V,QW,TL,FQW,FTL,TAUX,TAUY,                                 
     & DU,DV,CT_CTQ,DQW,DTL,CQ_CM,                                      
                                                                        
! LOGICAL LTIMER                                                        
     & LTIMER                                                           
     & )                                                                

                                                                        
      DO I=P1,P1+P_POINTS-1                                             
        RIB_LAND(I)=0.0                                             
        RIB_SSI(I)=0.0                                              
        FTL_LAND(I)=0.0                                             
        FTL_SSI(I)=0.0                                              
        E_LAND(I)=0.0                                               
        E_SSI(I)=0.0                                                
        EI_LAND(I)=0.0                                              
        EI_SICE(I)=0.0                                              
      ENDDO                                                             
                                                                        
      DO N=1,NTILES                                                     
        DO J=1,TILE_PTS(N)                                              
          L = TILE_INDEX(J,N)                                           
          I = LAND_INDEX(L)                          
          RIB_LAND(I)=RIB_LAND(I) +                                 
     &      RIB_TILE(L,N)*TILE_FRAC(L,N)
          FTL_LAND(I)=FTL_LAND(I) +                                 
     &      FTL_TILE(L,N)*TILE_FRAC(L,N)                                
          E_LAND(I)=E_LAND(I) +                                     
     &      FQW_TILE(L,N)*TILE_FRAC(L,N)                                
          EI_LAND(I)=EI_LAND(I) +                                   
     &      EI_TILE(L,N)*TILE_FRAC(L,N)                                 
        ENDDO
      ENDDO                                                             
                                                                        
      DO I=P1,P1+P_POINTS-1                                             
        IF(FLANDG(I).LT.1.0)THEN                                        
            RIB_SSI(I)=(RIB(I)-RIB_LAND(I)*FLANDG(I))        
     &        /(1.0-FLANDG(I))                                          
            FTL_SSI(I)=(FTL(I,1)-FTL_LAND(I)*FLANDG(I))         
     &        /(1.0-FLANDG(I))                                        
            E_SSI(I)=(FQW(I,1)-E_LAND(I)*FLANDG(I))             
     &        /(1.0-FLANDG(I))                                        
            EI_SICE(I)=(EI(I)-EI_LAND(I)*FLANDG(I))             
     &        /(1.0-FLANDG(I))
        ENDIF                                        
      ENDDO
                                                                        
*DECLARE BOUYTQ5B
*D BOUYTQ5B.45,ADM3F404.359  
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,P,CF,T,TL,Q,QCF,QCL
     &,BT,BQ,BF,BT_CLD,BQ_CLD,A_QS,A_DQSDT,DQSDT
*I ADM3F404.366   
                                                                        
      REAL
     & BT_CLD(P_FIELD,BL_LEVELS) 
!                             ! DUMMY Used in 6A boundary layer scheme
     &,BQ_CLD(P_FIELD,BL_LEVELS)
!                             ! DUMMY Used in 6A boundary layer scheme
     &,A_QS(P_FIELD,BL_LEVELS)
!                             ! DUMMY Used in 6A boundary layer scheme
     &,A_DQSDT(P_FIELD,BL_LEVELS)
!                             ! DUMMY Used in 6A boundary layer scheme
     &,DQSDT(P_FIELD,BL_LEVELS)
!                             ! DUMMY Used in 6A boundary layer scheme
*DECLARE BOUYTQ6A
*D BOUYTQ6A.2
*IF DEF,A03_8A                                                          
*D BOUYTQ6A.44,BOUYTQ6A.46   
     &,P,CF,T,TL,Q,QCF,QCL
     &,BT,BQ,BF,BT_CLD,BQ_CLD,A_QS,A_DQSDT,DQSDT
     &,L_BL_LSPICE,LTIMER
*I BOUYTQ6A.92    
                                                                        
      REAL
     & CF(P_FIELD,BL_LEVELS)  ! DUMMY Used in 7A boundary layer scheme
     &,TL(P_FIELD,BL_LEVELS)  ! DUMMY Used in 7A boundary layer scheme
     &,BF(P_FIELD,BL_LEVELS)  ! DUMMY Used in 7A boundary layer scheme

      LOGICAL
     & L_BL_LSPICE            ! DUMMY Used in 7A boundary layer scheme
*DECLARE CLDCTL1
*D ARE2F404.22
     &           NTILESDA,TILE_FIELDDA,DOLR,LW_DOWN,SW_TILE,            
     &           LAND_ALB,SICE_ALB,
*D ARN1F404.102
     &           L_RADHEAT,RADHEAT_DIM1,NLALBS,P_FIELDDA,               
*I @DYALLOC.747   
     &       NTILESDA,     !  COPY OF NTILES                            
     &       TILE_FIELDDA, ! Set to LAND_FIELD for MOSES II, 1 otherwise
*I ARN1F404.105   
     &       NLALBS,     !  Number of fields of land surface albedo
*D ARE2F404.23,ARE2F404.28   
     &       LW_DOWN(P_FIELDDA),            ! Surface downward LW       
C                                           ! (W/sq m)                  
     &       DOLR(P_FIELDDA),               ! TOA - surf upward LW rad  
     &       SW_TILE(TILE_FIELDDA,NTILESDA),! Surface net SW radiation  
C                                           ! on land tiles (W/sq m)    
     &       LAND_ALB(P_FIELDDA,NLALBS),    ! IN: Mean land albedo      
     &       SICE_ALB(P_FIELDDA,NLALBS),    ! IN: Mean sea-ice albedo   
*I @DYALLOC.755   

C
C MOSES 2.2 coastal tiling variables
      REAL
     &       FRACSOLID(P_FIELDDA),     ! Solid surface fraction in gridb
     &       ICE_FRACT(P_FIELDDA),     ! Ice fraction
     &       FLANDG(P_FIELDDA),        ! Land fraction 
     &       SW_NET_LAND(P_FIELDDA),   ! SW net local flux over land
     &       SW_NET_SICE(P_FIELDDA)    ! SW net local flux over sea-ice 

      REAL
     &       ALBSOLID                  ! Mean solid surface albedo

      LOGICAL 
     &       L_CTILE                   ! Coastal tiling switch
                                       ! Set to .TRUE.
                    
*D ARE2F404.29
     &     ,RADINCS((P_FIELDDA*(P_LEVELSDA+2+NTILESDA)+511)/512*512*2)  
*D ARE2F404.33
     &       N,J,                                                       
*D ARE2F404.37,CLDCTL1.112  
                                                                        
*D ARE2F404.38
        LEN=(P_FIELDDA*(P_LEVELS+2+NTILES)+511)/512*512*2               
C       no. words for LW/SW                                             
*D ARE2F404.41
        LEN=(P_FIELDDA*(P_LEVELS+2+NTILES)+511)/512*512                 
C       offset to 2nd RADINCS                                           
*D ARE2F404.43
C zenith angle adjustment and surface net SW on tiles)                  
*I CLDCTL1.387   


        IF ( H_SECT(3) .EQ. '07A' .OR.                                  
     &       H_SECT(3) .EQ. '08A' ) THEN                                
          L_CTILE=.TRUE. 
        ENDIF          

C  Set up GLOBAL fractional land field:                                 
        CALL LAND_TO_GLOBAL                                             
     &    (D1(JLAND),D1(JFRAC_LAND),FLANDG,LAND_PTS,P_FIELDDA)          

C Set up ice fraction field and solid fraction field                    
        DO I=1,P_FIELD
          ICE_FRACT(I) = D1(JICE_FRACTION+I-1)   
          FRACSOLID(I) = FLANDG(I) + (1.0-FLANDG(I))*ICE_FRACT(I)  
        ENDDO     

*I CLDCTL1.394   
        IF (L_CTILE) THEN
          
          DO I=1,P_FIELD
            SW_NET_LAND(I) = 0.0
            SW_NET_SICE(I) = 0.0
            SURF_RADFLUX(I) = 0.0       
          ENDDO

          DO I=FIRST_POINT,LAST_POINT                                   
                                                                        
! land_alb is only set on points where there is incoming sw radiation   
! at the previous timestep, therefore it will be zero over some         
! land points                                                           
            IF (FRACSOLID(I).GT.0.0) THEN
          
              IF (FLANDG(I).GT.0.0.AND.LAND_ALB(I,1).LE.0.0) THEN
                SW_NET_LAND(I) = 0.0
                SW_NET_SICE(I) = 0.0                      
              ELSE                                         
                ALBSOLID = ( FLANDG(I) * LAND_ALB(I,1) +
     &            (1.0-FLANDG(I)) * ICE_FRACT(I) * SICE_ALB(I,1) )  
     &              /FRACSOLID(I)                                
                                                                        
                IF (FLANDG(I).GT.0.0) THEN
                  SW_NET_LAND(I) = RADINCS(I)                
     &              * COS_ZENITH_ANGLE(I) / FRACSOLID(I)   
     &              * (1.0-LAND_ALB(I,1))/(1.0-ALBSOLID)   
                ENDIF
            
                IF (ICE_FRACT(I).GT.0.0) THEN
                  SW_NET_SICE(I) = RADINCS(I)      
     &              * COS_ZENITH_ANGLE(I) / FRACSOLID(I)  
     &              * (1.0-SICE_ALB(I,1))/(1.0-ALBSOLID)  
                
                SURF_RADFLUX(I) = SW_NET_SICE(I) * ICE_FRACT(I)   
                ENDIF                                                   
              ENDIF 
                                                                      
            ENDIF 
          ENDDO                                          
        ELSE ! If not L_CTILE          
          DO I=FIRST_POINT,LAST_POINT                                   
            SURF_RADFLUX(I) =                                           
     &          RADINCS(I) * COS_ZENITH_ANGLE(I) + RADINCS(I+LEN) 
          ENDDO
        ENDIF

*D CLDCTL1.396,CLDCTL1.397  
*I CLDCTL1.398   

C                                                                       
C Output SW land and sea-ice:                                   

        IF(SF(257,1)) THEN                                              
          DO I=FIRST_POINT, LAST_POINT                                  
            STASHWORK(SI(257,1,im_index)+I-1) = SW_NET_LAND(I)          
          ENDDO                                                         
          CALL EXTDIAG (STASHWORK,si(1,1,im_index),SF(1,1),257,257,INT9,
     &       ROW_LENGTH, STLIST, LEN_STLIST, STINDEX(1,1,1,im_index), 2,
     &       STASH_LEVELS, NUM_STASH_LEVELS+1,                          
     &       STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,                     
     &       im_ident,1,                                                
*CALL ARGPPX                                                            
     &     ICODE, CMESSAGE)                                             
        ENDIF                                                           
C                                                                       
C                                                                       
        IF(SF(258,1)) THEN                                              
          DO I=FIRST_POINT, LAST_POINT                                  
            STASHWORK(SI(258,1,im_index)+I-1) = SW_NET_SICE(I)          
          ENDDO                                                         
          CALL EXTDIAG (STASHWORK,si(1,1,im_index),SF(1,1),258,258,INT9,
     &       ROW_LENGTH, STLIST, LEN_STLIST, STINDEX(1,1,1,im_index), 2,
     &       STASH_LEVELS, NUM_STASH_LEVELS+1,                          
     &       STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,                     
     &       im_ident,1,                                                
*CALL ARGPPX                                                            
     &       ICODE, CMESSAGE)                                           
        ENDIF                                                           
*D ARE2F404.46
        IF ( H_SECT(3) .EQ. '07A' .OR.                                  
     &       H_SECT(3) .EQ. '08A' ) THEN                                
*D ABX1F405.97,ABX1F405.103  
C Set the net SW flux on land tiles                                     
          DO N=1,NTILES                                                 
            DO L=LAND1,LAND1+LAND_PTS-1                                 
              I = LAND_LIST(L)                                          
              J = I + (P_LEVELS + 1 + N)*P_FIELD                        
              SW_TILE(L,N) = RADINCS(J)*COS_ZENITH_ANGLE(I)             
            ENDDO                                                       
          ENDDO                                                         
*D ARE2F404.58,ARE1F405.6    
C Set the surface downward and TOA outward LW fluxes                    
           DO I=FIRST_POINT,LAST_POINT                                  
             J = I + (P_LEVELS + 2)*P_FIELD                             
             LW_DOWN(I) = RADINCS(J+LEN)                                
             DOLR(I) = RADINCS(J+LEN+P_FIELD)                           
*I ARE1F405.8     
                                                                        
C Overwrite SURF_RADFLUX for sea-ice with net SW + downward LW          
           DO I=FIRST_POINT,LAST_POINT                                  
             SURF_RADFLUX(I) = SURF_RADFLUX(I) +
     &                           ICE_FRACT(I)*LW_DOWN(I)                
           ENDDO                                                        
                                                                        
*IF DEF,MPP                                                             
          CALL SWAPB_LAND(SW_TILE,LAND_FIELD,P_FIELD,                   
     &                    ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,            
     &                    NTILES,LAND_LIST)                             
          CALL SWAPBOUNDS(LW_DOWN,ROW_LENGTH,P_ROWS,                    
     &                    EW_Halo,NS_Halo,1)                            
          CALL SWAPBOUNDS(DOLR,ROW_LENGTH,P_ROWS,                       
     &                    EW_Halo,NS_Halo,1)                            
*ENDIF                                                                  
*D ABX1F405.110
     &       //'encountered in CLDCTL.'                                 
*DECLARE CNTLATM
*I CNTLATM.112   
     &   LTLEADS        ,           !  Let leads temp vary if .TRUE.    
*D CNTLATM.136
     & LFROUDE, LGWLINP, LLINTS, LWHITBROM, LEMCORR, LTLEADS,           
*D CNTLATM.166
     & LFROUDE, LGWLINP, LLINTS, LWHITBROM, LEMCORR, LTLEADS,           
*DECLARE COMPT2A
*D COMPT2A.62,COMPT2A.66   
     &,FRACN,FRACM                ! WORK Fractions used in the spreading
C                                 !      calculation.
*D COMPT2A.96,COMPT2A.97   

        FRACN=FRAC(L,N)
        FRACN=MAX(FRACN,FRAC_SEED)

        FRACM=FRAC(L,M)
        FRACM=MAX(FRACM,FRAC_SEED)

        P1 = GAMMA/FRACN-FORW*DB_DFRAC(L,N,N)                           
        P2 = GAMMA/FRACM-FORW*DB_DFRAC(L,M,M)                           
*D COMPT2A.109
*D COMPT2A.126
*D COMPT2A.148,COMPT2A.149  

        FRACN=FRAC(L,N)
        FRACN=MAX(FRACN,FRAC_SEED)

        DENOM = GAMMA/FRACN-FORW*DB_DFRAC(L,N,N)           
*D COMPT2A.178,COMPT2A.179  

        FRACN=FRAC(L,N)
        FRACN=MAX(FRACN,FRAC_SEED)

        FRACM=FRAC(L,M)
        FRACM=MAX(FRACM,FRAC_SEED)

        P1 = GAMMA/FRACN-FORW*DB_DFRAC(L,N,N)                           
        P2 = GAMMA/FRACM-FORW*DB_DFRAC(L,M,M)                           
*D COMPT2A.191
*D COMPT2A.208
*DECLARE CRADINCS
*D ARE2F404.69
     &        RADINCS ( (P_FIELD*(P_LEVELS+2+9)+511)/512*512*2 )   
*DECLARE C_ROUGH
*I ARN1F404.63    
*IF DEF,A03_8A
!!----------------------------------------------------------------------
!!!-----------COMDECK C_ROUGH FOR SUBROUTINE SF_EXCH----------
! Z0FSEA = roughness length for free convective heat and moisture
!          transport over the sea (m).
!          DUMMY VARIABLE - Only used in 7A boundary layer scheme
! Z0HSEA = roughness length for heat and moisture transport
!          over the sea (m).
! Z0MIZ  = roughness length for heat, moisture and momentum over
!          the Marginal Ice Zone (m).
! Z0SICE = roughness length for heat, moisture and momentum over
!          sea-ice (m).
      REAL Z0FSEA,Z0HSEA,Z0MIZ,Z0SICE

      PARAMETER(Z0FSEA = 1.3E-3,
     &          Z0HSEA = 4.0E-5,
     &          Z0MIZ  = 1.0E-1,
     &          Z0SICE = 3.0E-3)
!!----------------------------------------------------------------------
*ENDIF
*DECLARE DARCY5A
*I DARCY5A.24    
     &,                 DWFLUX_DSTHU1,DWFLUX_DSTHU2
*I DARCY5A.43    
!  4.6      2/99     Extended to provide derivatives. Peter Cox
*I DARCY5A.84    
     &,DWFLUX_DSTHU1(NPNTS) ! OUT The rate of change of the explicit
!                           !     flux with STHU1 (kg/m2/s).
     &,DWFLUX_DSTHU2(NPNTS) ! OUT The rate of change of the explicit
!                           !     flux with STHU2 (kg/m2/s).
*I DARCY5A.88    

      REAL
     & DTHK_DTH1,DTHK_DTH2  ! WORK DTHETAK/DTHETA(1:2).
     &,DK_DTH1,DK_DTH2      ! WORK DK/DTHETA(1:2) (kg/m2/s).
     &,PD                   ! WORK Hydraulic potential difference (m).
*I DARCY5A.97    
     &,DK_DTHK(NPNTS)       ! WORK The rate of change of K with THETAK
!                           !      (kg/m2/s).
*I DARCY5A.99    
     &,DPSI_DTH(NPNTS,2)    ! WORK The rate of change of PSI with
!                           !      THETA(1:2) (m).       
*D DARCY5A.121
            DPSI_DTH(I,N)=0.0   
          ELSE
*D DARCY5A.123,DARCY5A.124  
            DPSI_DTH(I,N)=-B(I)*PSI(I,N)/THETA(I,N)
*I DARCY5A.136   
      DTHK_DTH1=DZ2/(DZ1+DZ2)
      DTHK_DTH2=DZ1/(DZ2+DZ1)

*D DARCY5A.140
      CALL HYD_CON (NPNTS,SOIL_PTS,SOIL_INDEX,B,KS,THETAK,K,DK_DTHK
     &,             LTIMER)          
*D DARCY5A.147
        PD=(2.0*(PSI(I,2)-PSI(I,1))/(DZ2+DZ1)+1)                   
        WFLUX(I)=K(I)*PD 

!-----------------------------------------------------------------------
! Calculate the rate of change of WFLUX with respect to the STHU1 and
! STHU2.
!-----------------------------------------------------------------------
        DK_DTH1=DK_DTHK(I)*DTHK_DTH1
        DK_DTH2=DK_DTHK(I)*DTHK_DTH2
        DWFLUX_DSTHU1(I)=DK_DTH1*PD-2*K(I)*DPSI_DTH(I,1)/(DZ1+DZ2)     
        DWFLUX_DSTHU2(I)=DK_DTH2*PD+2*K(I)*DPSI_DTH(I,2)/(DZ1+DZ2)
*DECLARE DECAY2A
*D DECAY2A.52,DECAY2A.56   
*D DECAY2A.69
*DECLARE DESCENT
*D DESCENT.4
     + DENOM_MIN                  ! Minimum value for the denominator
C                                 ! of the update equation. Ensures 
C                                 ! that gradient descent does not  
C                                 ! lead to an unstable solution.  
     +,GAMMA_EQ                   ! Inverse timestep for gradient       
*D ABX1F405.1722
      PARAMETER(DENOM_MIN=1.0E-6, GAMMA_EQ = 1.0E-1, ITER_EQ = 10)      
*DECLARE DFPLN3A
*I DFPLN3A.25    
     &   , T_SOLID, T_SEA, L_CTILE                                      
*D ADB1F401.54,ADB1F401.55   
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION             
     &   , FLANDG, PLANCK_FREEZE_SEA, PLANCK_LEADS_SEA                  
*I DFPLN3A.54    
     &   , L_CTILE                                                      
!             Coastal tiling switch                                     
*I DFPLN3A.63    
     &   , T_SOLID(NPD_PROFILE)                                         
!             TEMPERATURES AT SOLID SURFACE                             
     &   , T_SEA(NPD_PROFILE)                                           
!             SURFACE TEMPERATURE OVER OPEN SEA                         
*I DFPLN3A.74    
     &   , FLANDG(NPD_PROFILE)                                          
!             GATHERED LAND FRACTION                                    
*D ADB1F401.61,ADB1F401.64   
     &     N_FRAC_SOL_POINT                                             
!             NUMBER OF POINTS WITH FRACTIONAL ICE COVER/LAND COVER     
     &   , I_FRAC_SOL_POINT(NPD_PROFILE)                                
!             INDICES OF POINTS WITH FRACTIONAL ICE COVER/LAND COVER    
*D ADB1F401.67
!             FRACTION OF SEA-ICE IN SEA PORTION OF GRID BOX            
*I ADB1F401.71    
     &   , PLANCK_LEADS_SEA(NPD_PROFILE)                                
!             PLANCK FUNCTION OVER SEA LEADS                            
*D ADB1F401.83,ADB1F401.85   
!             PLANCKIAN OF SOLID SURFACE (GATHERED OVER SOLID POINTS)   
     &   , SEAFRAC(NPD_PROFILE)                                         
!             FRACTION OF OPEN SEA IN GRID-BOX                          
*I DFPLN3A.144   
      IF(L_CTILE)THEN                                                   
!     SOURCE AT THE OPEN SEA SURFACE.                                   
!     CALCULATE OVER ALL POINTS EVEN THOUGH IT IS                       
!     OVERWRITTEN WHERE THERE IS LAND OR SEA-ICE.                       
      DO L=1, N_PROFILE                                                 
         T_RATIO(L)=T_SEA(L)/T_REF_PLANCK                               
         IF(FLANDG(L).GE.1.0.OR.ICE_FRACTION(L).GE.1.0)                 
     &      T_RATIO(L)=T_SOLID(L)/T_REF_PLANCK                          
                                                                        
         PLANCK_GROUND(L)=THERMAL_COEFFICIENT(N_DEG_FIT)                
      ENDDO                                                             
      DO J=N_DEG_FIT-1, 0, -1                                           
         DO L=1, N_PROFILE                                              
            PLANCK_GROUND(L)=PLANCK_GROUND(L)*T_RATIO(L)                
     &         +THERMAL_COEFFICIENT(J)                                  
         ENDDO                                                          
      ENDDO                                                             
!                                                                       
C                                                                       
! Initialise to zero:                                                   
      DO LG=1,NPD_PROFILE                                               
        PLANCK_LEADS_SEA(LG)=0.0                                        
      ENDDO                                                             
!     SET THE SOURCE FUNCTION OVER OPEN SEA LEADS.                      
!     DETERMINE THE TEMPERATURE OF THE NON-SEA FRACTION.                
!     CALCULATE THE SOURCE FUNCTION AT POINTS WITH SOLID SURFACE        
      DO L=1, N_FRAC_SOL_POINT                                          
         LG=I_FRAC_SOL_POINT(L)                                         
         PLANCK_LEADS_SEA(LG)=PLANCK_GROUND(LG)                         
         SEAFRAC(LG)=(1.-FLANDG(LG))*(1.0E+00-ICE_FRACTION(LG))         
         T_RATIO(L)=T_SOLID(LG)/T_REF_PLANCK                            
         PLANCK_GROUND_G(L)=THERMAL_COEFFICIENT(N_DEG_FIT)              
      ENDDO                                                             
      DO J=N_DEG_FIT-1, 0, -1                                           
         DO L=1, N_FRAC_SOL_POINT                                       
            PLANCK_GROUND_G(L)=PLANCK_GROUND_G(L)*T_RATIO(L)            
     &         +THERMAL_COEFFICIENT(J)                                  
         ENDDO                                                          
      ENDDO                                                             
!     DETERMINE THE OVERALL PLANCKIAN FUNCTION OF THE SURFACE.          
      DO L=1, N_FRAC_SOL_POINT                                          
         LG=I_FRAC_SOL_POINT(L)                                         
         PLANCK_GROUND(LG)=(1.-SEAFRAC(LG))*PLANCK_GROUND_G(L)          
     &      +PLANCK_LEADS_SEA(LG)*SEAFRAC(LG)                           
      ENDDO                                                             
                                                                        
      ELSE                      !End of L_CTILE                         
                                                                        
*D ADB1F401.102,ADB1F401.103  
      DO L=1, N_FRAC_SOL_POINT                                          
         LG=I_FRAC_SOL_POINT(L)                                         
*D ADB1F401.109
      DO L=1, N_FRAC_SOL_POINT                                          
*D ADB1F401.114
         DO L=1, N_FRAC_SOL_POINT                                       
*D ADB1F401.121,ADB1F401.122  
      DO L=1, N_FRAC_SOL_POINT                                          
         LG=I_FRAC_SOL_POINT(L)                                         
*I ADB1F401.126   
      ENDIF                                                             
                                                                        
*DECLARE DIAG3A
*D DIAG3A.29
!     Dummy arguments                                                   
*D DIAG3A.31,DIAG3A.32   
     &    N                                                             
!           Length of array                                             
*D DIAG3A.34,DIAG3A.35   
     &    X(N)                                                          
!           Array to be zeroed                                          
*D DIAG3A.37
!     Local variables                                                   
*D DIAG3A.39,DIAG3A.40   
     &    I                                                             
!           loop variable                                               
*D DIAG3A.45
        X(I)=0.0E+00                                                    
*D DIAG3A.72,DIAG3A.78   
     &  , SEA_FLUX                                                      
     &  , L_SURFACE_DOWN_FLUX, SURFACE_DOWN_FLUX                        
     &  , L_SURF_DOWN_CLR, SURF_DOWN_CLR                                
     &  , L_SURF_UP_CLR, SURF_UP_CLR                                    
     &  , L_FLUX_BELOW_690NM_SURF                                       
     &  , FLUX_BELOW_690NM_SURF, FL_SEA_BELOW_690NM_SURF                
     &  , L_MOSES_II                                                    
     &  , SURF_VIS_DIR, SURF_VIS_DIF, SURF_NIR_DIR, SURF_NIR_DIF        
     &  , NPD_PROFILE                                                   
     & )                                                               
*D DIAG3A.85
!     Dummy arguments                                                   
*D DIAG3A.87
!     Dimensions of arrays                                              
*D DIAG3A.89,DIAG3A.90   
     &    NPD_PROFILE                                                   
!           Maximum number of atmospheric profiles                      
*D DIAG3A.93,DIAG3A.94   
     &    N_PROFILE                                                     
!           Number of atmospheric profiles                              
*D DIAG3A.96
!     Switches for diagnostics:                                         
*D DIAG3A.98,DIAG3A.105  
     &    L_FLUX_BELOW_690NM_SURF                                       
!           Flux below 690nm at surface to be calculated                
     &  , L_MOSES_II                                                   
!           Surface SW fluxes required for MOSES II 
     &  , L_SURFACE_DOWN_FLUX                                           
!           Downward surface flux required                              
     &  , L_SURF_DOWN_CLR                                               
!           Calculate downward clear flux                               
     &  , L_SURF_UP_CLR                                                 
!           Calculate upward clear flux                                 
*D DIAG3A.107
!     Surface fluxes for coupling or diagnostic use                     
*D DIAG3A.109,DIAG3A.118  
     &    SEA_FLUX(NPD_PROFILE)                                         
!           Net downward flux into sea                                  
     &  , SURFACE_DOWN_FLUX(NPD_PROFILE)                                
!           Downward flux at surface                                    
     &  , SURF_DOWN_CLR(NPD_PROFILE)                                    
!           Clear-sky downward flux at surface                          
     &  , SURF_UP_CLR(NPD_PROFILE)                                      
!           Clear-sky upward flux at surface                            
     &  , FLUX_BELOW_690NM_SURF(NPD_PROFILE)                            
!          GRID BOX MEAN NET SURFACE FLUX BELOW 690NM                   
     &  , FL_SEA_BELOW_690NM_SURF(NPD_PROFILE)                          
!          OPEN SEA NET SURFACE FLUX BELOW 690NM                        
     &  , SURF_VIS_DIR(NPD_PROFILE)                                     
!           Downward surface direct beam visible flux                   
     &  , SURF_VIS_DIF(NPD_PROFILE)                                     
!           Downward surface diffuse visible flux                       
     &  , SURF_NIR_DIR(NPD_PROFILE)                                     
!           Downward surface direct beam near-infrared flux             
     &  , SURF_NIR_DIF(NPD_PROFILE)                                     
!           Downward surface diffuse near-infrared flux                 
*D DIAG3A.125
        CALL R2_ZERO_1D(N_PROFILE, SURFACE_DOWN_FLUX)                   
*D DIAG3A.129
        CALL R2_ZERO_1D(N_PROFILE, SURF_DOWN_CLR)                       
*D DIAG3A.133
        CALL R2_ZERO_1D(N_PROFILE, SURF_UP_CLR)                         
*D DIAG3A.137
        CALL R2_ZERO_1D(N_PROFILE, FLUX_BELOW_690NM_SURF)               
        CALL R2_ZERO_1D(N_PROFILE, FL_SEA_BELOW_690NM_SURF)             
*I DIAG3A.139   
      IF (L_MOSES_II) THEN                                              
        CALL R2_ZERO_1D(N_PROFILE, SURF_VIS_DIR)                        
        CALL R2_ZERO_1D(N_PROFILE, SURF_VIS_DIF)                        
        CALL R2_ZERO_1D(N_PROFILE, SURF_NIR_DIR)                        
        CALL R2_ZERO_1D(N_PROFILE, SURF_NIR_DIF)                        
      ENDIF                                                             
*I ADB1F401.131   
!       4.6             07-05-99                Calculation of fluxes   
!                                               into the sea changed    
!                                               to use albedos for      
!                                               sea. Code for obsolete  
!                                               net solvers removed.    
!                                               (J. M. Edwards)         
!       5.3             22-02-02    Add two new "blue band" fluxes.     
!                                               (N. Gedney)             
*D DIAG3A.163,DIAG3A.176  
      SUBROUTINE R2_COUPLE_DIAG(N_PROFILE, ISOLIR                       
     &  , ALBEDO_FIELD_DIFF, ALBEDO_FIELD_DIR                           
     &  , ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR                               
     &  , FLANDG, ICE_FRACTION                                          
     &  , PLANCK_FREEZE_SEA, PLANCK_LEADS_SEA                           
     &  , PLANCK_AIR_SURFACE, THERMAL_SOURCE_GROUND                     
     &  , FLUX_DOWN, FLUX_UP, FLUX_DIRECT                               
     &  , FLUX_DOWN_CLEAR, FLUX_UP_CLEAR, FLUX_DIRECT_CLEAR             
     &  , WEIGHT_690NM                                                  
     &  , SEA_FLUX                                                      
     &  , L_SURFACE_DOWN_FLUX, SURFACE_DOWN_FLUX                        
     &  , L_SURF_DOWN_CLR, SURF_DOWN_CLR                                
     &  , L_SURF_UP_CLR, SURF_UP_CLR                                    
     &  , L_FLUX_BELOW_690NM_SURF                                       
     &  , FLUX_BELOW_690NM_SURF, FL_SEA_BELOW_690NM_SURF                
     &  , L_MOSES_II, L_CTILE                                           
     &  , SURF_VIS_DIR, SURF_VIS_DIF, SURF_NIR_DIR, SURF_NIR_DIF        
     &  , NPD_PROFILE                                                   
     &  )                                                               
*D DIAG3A.183,DIAG3A.184  
!     Comdecks included                                                 
!     Spectral regions                                                  
*D DIAG3A.187
!     Dummy Arguments                                                   
*D DIAG3A.189
!     Dimensions of arrays                                              
*D DIAG3A.191,DIAG3A.192  
     &    NPD_PROFILE                                                   
!           Maximum number of atmospheric profiles                      
*D DIAG3A.195,DIAG3A.198  
     &    N_PROFILE                                                     
!           Number of atmospheric profiles                              
     &  , ISOLIR                                                        
!           Spectral region                                             
*D DIAG3A.200
!     Logical switches for the code                                     
*D DIAG3A.202,DIAG3A.203  
     &    L_NET                                                         
!           Flag for net fluxes                                         
*D DIAG3A.205
!     Switches for diagnostics:                                         
*D DIAG3A.207,DIAG3A.214  
     &    L_FLUX_BELOW_690NM_SURF                                       
!           Flux below 690nm at surface to be calculated                
     &  , L_MOSES_II                                                    
!           Surface SW required for MOSES II                            
     &  , L_CTILE                                                       
!           Switch for coastal tiling                                   
     &  , L_SURFACE_DOWN_FLUX                                           
!           Downward surface flux required                              
     &  , L_SURF_DOWN_CLR                                               
!           Calculate downward clear flux                               
     &  , L_SURF_UP_CLR                                                 
!           Calculate upward clear flux                                 
*D DIAG3A.216
!     Albedos                                                           
*D DIAG3A.218,DIAG3A.225  
     &    ALBEDO_FIELD_DIFF(NPD_PROFILE)                                
!           Diffuse albedo meaned over grid box                         
     &  , ALBEDO_FIELD_DIR(NPD_PROFILE)                                 
!           Direct albedo meaned over grid box                          
     &  , ALBEDO_SEA_DIFF(NPD_PROFILE)                                  
!           Diffuse albedo of open sea                                  
     &  , ALBEDO_SEA_DIR(NPD_PROFILE)                                   
!           Direct albedo meaned of open sea                            
*D DIAG3A.228,ADB1F401.136  
     &    THERMAL_SOURCE_GROUND(NPD_PROFILE)                            
!           Thermal source at ground                                    
     &  , PLANCK_AIR_SURFACE(NPD_PROFILE)                               
!           Planck function at near-surface air temperature in band     
*D ADB1F401.137,ADB1F401.142  
!     Arguments relating to sea ice.                                    
*D ADB1F401.144,ADB1F401.148  
     &    PLANCK_FREEZE_SEA                                             
!           Planck function over freezing sea                           
     &   , PLANCK_LEADS_SEA(NPD_PROFILE)                                
!           Planck function over sea leads                              
     &   , FLANDG(NPD_PROFILE)                                          
!            Land fraction                                              
     &   , ICE_FRACTION(NPD_PROFILE)                                    
!            FRACTION OF SEA-ICE IN SEA PORTION OF GRID BOX!            
*D DIAG3A.232,DIAG3A.233  
     &    WEIGHT_690NM                                                  
!           Weighting applied to band for region below 690 nm           
*D DIAG3A.235
!     Calculated fluxes                                                 
*D DIAG3A.237,DIAG3A.248  
     &    FLUX_DOWN(NPD_PROFILE)                                        
!           Total downward or net flux at surface                       
     &  , FLUX_DIRECT(NPD_PROFILE)                                      
!           Direct solar flux at surface                                
     &  , FLUX_UP(NPD_PROFILE)                                          
!           Upward flux at surface                                      
     &  , FLUX_DOWN_CLEAR(NPD_PROFILE)                                  
!           Total clear-sky downward or net flux at surface             
     &  , FLUX_UP_CLEAR(NPD_PROFILE)                                    
!           Clear-sky upward flux at surface                            
     &  , FLUX_DIRECT_CLEAR(NPD_PROFILE)                                
!           Clear-sky direct solar flux at surface                      
*D DIAG3A.251
!     Surface fluxes for coupling or diagnostic use                     
*D DIAG3A.253,DIAG3A.262  
     &    SEA_FLUX(NPD_PROFILE)                                         
!           Net downward flux into sea                                  
     &  , SURFACE_DOWN_FLUX(NPD_PROFILE)                                
!           Downward flux at surface                                    
     &  , SURF_DOWN_CLR(NPD_PROFILE)                                    
!           Clear-sky downward flux at surface                          
     &  , SURF_UP_CLR(NPD_PROFILE)                                      
!           Clear-sky upward flux at surface                            
     &  , FLUX_BELOW_690NM_SURF(NPD_PROFILE)                            
!          GRID BOX MEAN NET SURFACE FLUX BELOW 690NM                   
     &  , FL_SEA_BELOW_690NM_SURF(NPD_PROFILE)                          
!          OPEN SEA NET SURFACE FLUX BELOW 690NM                        
     &  , SURF_VIS_DIR(NPD_PROFILE)                                     
!           Downward surface direct beam visible flux                   
     &  , SURF_VIS_DIF(NPD_PROFILE)                                     
!           Downward surface diffuse visible flux                       
     &  , SURF_NIR_DIR(NPD_PROFILE)                                     
!           Downward surface direct beam near-infrared flux             
     &  , SURF_NIR_DIF(NPD_PROFILE)                                     
!           Downward surface diffuse near-infrared flux                 
*D DIAG3A.265
!     Local variables                                                   
*D DIAG3A.267,DIAG3A.268  
     &    L                                                             
!           Loop variable                                               
*D DIAG3A.272,DIAG3A.313  
!     This is the flux into the sea over the ice-free parts of the      
!     grid-box. The model is that the downward fluxes are uniform       
!     across the grid-box, but the upward fluxes are not. At this       
!     stage no weighting by the actual ice-free area is carried out:    
!     that will be done later.                                          
      IF (ISOLIR.EQ.IP_SOLAR) THEN                                      
        DO L=1, N_PROFILE                                               
          IF (FLANDG(L).LT.1.0.AND.ICE_FRACTION(L).LT.1.0) THEN         
            SEA_FLUX(L)=SEA_FLUX(L)                                     
     &        +FLUX_DOWN(L)*(1.0E+00-ALBEDO_SEA_DIFF(L))                
     &        +FLUX_DIRECT(L)*(ALBEDO_SEA_DIFF(L)-ALBEDO_SEA_DIR(L))    
          ENDIF                                                         
        ENDDO                                                           
      ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN                             
        IF(L_CTILE)THEN                                                 
        DO L=1, N_PROFILE                                               
          IF (FLANDG(L).LT.1.0.AND.ICE_FRACTION(L).LT.1.0) THEN         
            SEA_FLUX(L)=SEA_FLUX(L)                                     
     &        +(1.0E+00-ALBEDO_SEA_DIFF(L))                             
     &        *(FLUX_DOWN(L)+PLANCK_AIR_SURFACE(L)                      
     &        -PLANCK_LEADS_SEA(L))                                     
          ENDIF                                                         
        ENDDO                                                           
        ELSE                                                            
        DO L=1, N_PROFILE                                               
          IF (FLANDG(L).LT.1.0) THEN                                    
            SEA_FLUX(L)=SEA_FLUX(L)                                     
     &        +(1.0E+00-ALBEDO_SEA_DIFF(L))                             
     &        *(FLUX_DOWN(L)+PLANCK_AIR_SURFACE(L)                      
     &        -PLANCK_FREEZE_SEA)                                       
          ENDIF                                                         
        ENDDO                                                           
        ENDIF                                                           
*D DIAG3A.317,DIAG3A.338  
        IF (ISOLIR.EQ.IP_SOLAR) THEN                                    
          DO L=1, N_PROFILE                                             
            SURFACE_DOWN_FLUX(L)=SURFACE_DOWN_FLUX(L)                   
     &        +FLUX_DOWN(L)                                             
          ENDDO                                                         
        ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN                           
          DO L=1, N_PROFILE                                             
            SURFACE_DOWN_FLUX(L)=SURFACE_DOWN_FLUX(L)                   
     &        +FLUX_DOWN(L)+PLANCK_AIR_SURFACE(L)                       
          ENDDO                                                         
        ENDIF                                                           
*D DIAG3A.342,DIAG3A.365  
        IF (ISOLIR.EQ.IP_SOLAR) THEN                                    
          DO L=1, N_PROFILE                                             
            SURF_DOWN_CLR(L)=SURF_DOWN_CLR(L)                           
     &        +FLUX_DOWN_CLEAR(L)                                       
          ENDDO                                                         
        ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN                           
          DO L=1, N_PROFILE                                             
            SURF_DOWN_CLR(L)=SURF_DOWN_CLR(L)                           
     &        +FLUX_DOWN_CLEAR(L)+PLANCK_AIR_SURFACE(L)                 
          ENDDO                                                         
        ENDIF                                                           
*D DIAG3A.369,DIAG3A.393  
        IF (ISOLIR.EQ.IP_SOLAR) THEN                                    
          DO L=1, N_PROFILE                                             
            SURF_UP_CLR(L)=SURF_UP_CLR(L)                               
     &        +FLUX_UP_CLEAR(L)                                         
          ENDDO                                                         
        ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN                           
          DO L=1, N_PROFILE                                             
            SURF_UP_CLR(L)=SURF_UP_CLR(L)                               
     &        +FLUX_UP_CLEAR(L)+PLANCK_AIR_SURFACE(L)                   
          ENDDO                                                         
        ENDIF                                                           
*D DIAG3A.396
!     This diagnostic is available only in the solar region. Over       
!     sea-points it refers only to the flux over the open sea           
!     (see the comments about sea_flux above). Over land, both the      
!     upward and downward fluxes are taken as uniform.                  
*D DIAG3A.398,DIAG3A.406  
        IF (ISOLIR.EQ.IP_SOLAR) THEN                                    
         IF (L_CTILE) THEN                                              
          DO L=1, N_PROFILE                                             
            FLUX_BELOW_690NM_SURF(L)=FLUX_BELOW_690NM_SURF(L)           
     &          +WEIGHT_690NM*(FLUX_DOWN(L)-FLUX_UP(L))                 
            IF(FLANDG(L).LT.1.0.AND.ICE_FRACTION(L).LT.1.0) THEN        
              FL_SEA_BELOW_690NM_SURF(L)                                
     &          =FL_SEA_BELOW_690NM_SURF(L)                             
     &          +WEIGHT_690NM                                           
     &          *(FLUX_DOWN(L)*(1.0E+00-ALBEDO_SEA_DIFF(L))             
     &          +FLUX_DIRECT(L)*(ALBEDO_SEA_DIFF(L)                     
     &          -ALBEDO_SEA_DIR(L)))                                    
            ENDIF                                                       
          ENDDO                                                         
         ELSE                                                           
          DO L=1, N_PROFILE                                             
            IF (FLANDG(L).GT.0.0) THEN                                  
              FLUX_BELOW_690NM_SURF(L)=FLUX_BELOW_690NM_SURF(L)         
     &          +WEIGHT_690NM*(FLUX_DOWN(L)-FLUX_UP(L))                 
*D DIAG3A.408,DIAG3A.411  
              FLUX_BELOW_690NM_SURF(L)=FLUX_BELOW_690NM_SURF(L)         
     &          +WEIGHT_690NM                                           
     &          *(FLUX_DOWN(L)*(1.0E+00-ALBEDO_SEA_DIFF(L))             
     &          +FLUX_DIRECT(L)*(ALBEDO_SEA_DIFF(L)-ALBEDO_SEA_DIR(L))) 
*I DIAG3A.412   
          ENDDO                                                         
*I DIAG3A.413   
        ENDIF                                                           
*I DIAG3A.415   
!     SURFACE SHORTWAVE DIAGNOSTICS REQUIRED FOR MOSES II  
      IF (L_MOSES_II) THEN                                              
        IF (ISOLIR.EQ.IP_SOLAR) THEN                   
          DO L=1, N_PROFILE                                             
            SURF_VIS_DIR(L) = SURF_VIS_DIR(L) +                         
     &                        WEIGHT_690NM*FLUX_DIRECT(L)               
            SURF_NIR_DIR(L) = SURF_NIR_DIR(L) +                         
     &                        (1. - WEIGHT_690NM)*FLUX_DIRECT(L)        
            SURF_VIS_DIF(L) = SURF_VIS_DIF(L) +                         
     &                        WEIGHT_690NM*(FLUX_DOWN(L)-FLUX_DIRECT(L))
            SURF_NIR_DIF(L) = SURF_NIR_DIF(L) +                         
     &                 (1. - WEIGHT_690NM)*(FLUX_DOWN(L)-FLUX_DIRECT(L))
          ENDDO                                                         
        ENDIF                                                           
      ENDIF                                                             
*DECLARE EXFXUV5B
*D ACB1F405.4
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*D EXFXUV5B.138,EXFXUV5B.149  
*DECLARE FCDCH6A
*D FCDCH6A.2
*IF DEF,A03_8A
*D FCDCH6A.21
!!!   SUBROUTINES FCDCH_SEA AND FCDCH_LAND-----------------------------
*I FCDCH6A.38    

!     SUBROUTINE FCDCH_SEA--------------------------------------------- 
!                                                                       
!     Transfer coefficients for sea, sea-ice and leads                  
!                                                                       
!     ----------------------------------------------------------------- 
*D FCDCH6A.40,FCDCH6A.43   
      SUBROUTINE FCDCH_SEA(
     & P_POINTS,P_FIELD,P1,LAND_MASK,
     & RIB,DB,VSHR,Z0M,Z0H,Z0F,ZH,Z1_UV,Z1_TQ,
     & CDV,CHV,V_S,RECIP_L_MO,LTIMER
*D FCDCH6A.54
*D FCDCH6A.70,FCDCH6A.72   
*D FCDCH6A.79,FCDCH6A.81   
*D FCDCH6A.84,FCDCH6A.86   
*I FCDCH6A.89    

      REAL
     & RIB(P_FIELD)  ! DUMMY Used in 7A boundary layer scheme
     &,Z0F(P_FIELD)  ! DUMMY Used in 7A boundary layer scheme
*D FCDCH6A.104
      EXTERNAL TIMER,PHI_M_H_SEA
*D FCDCH6A.124
      INTEGER I     ! Loop counter; horizontal field index.
*D FCDCH6A.130
*D FCDCH6A.142
*D FCDCH6A.144,FCDCH6A.152  
      DO I=P1,P1+P_POINTS-1                                             
        IF ( .NOT. LAND_MASK(I) ) THEN
*D FCDCH6A.164
        ENDIF  ! LAND_MASK
*D FCDCH6A.166

      CALL PHI_M_H_SEA (P_POINTS,P_FIELD,P1,LAND_MASK,
*D FCDCH6A.171,FCDCH6A.181  
      DO I=P1,P1+P_POINTS-1                                             
        IF ( .NOT. LAND_MASK(I) ) THEN
*D FCDCH6A.186
            V_S(I) = BETA *
*D FCDCH6A.188
*D FCDCH6A.194
*D FCDCH6A.196
          CHV(I) = ( VKMAN / PHI_H(I) ) * V_S(I)
*D FCDCH6A.198,FCDCH6A.200  
        ENDIF  ! LAND_MASK
*D FCDCH6A.207,FCDCH6A.217  
        DO I=P1,P1+P_POINTS-1                                           
          IF ( .NOT. LAND_MASK(I) ) THEN
*D FCDCH6A.220
*D FCDCH6A.224
*D FCDCH6A.227
*D FCDCH6A.231
          ENDIF  ! LAND_MASK
*D FCDCH6A.233
        CALL PHI_M_H_SEA (P_POINTS,P_FIELD,P1,LAND_MASK,
*D FCDCH6A.238
*D FCDCH6A.240,FCDCH6A.249  
        DO I=P1,P1+P_POINTS-1                                           
          IF ( .NOT. LAND_MASK(I) ) THEN
            CHV(I) = ( VKMAN / PHI_H(I) ) * V_S(I)
*D FCDCH6A.251,FCDCH6A.253  
          ENDIF  ! LAND_MASK
*I FCDCH6A.255   

!-----------------------------------------------------------------------
!! Set CD's and CH's to be dimensionless paremters
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1
        IF ( .NOT. LAND_MASK(I) ) THEN
          CDV(I) = CDV(I) / VSHR(I)
          CHV(I) = CHV(I) / VSHR(I)
        ENDIF  ! LAND_MASK
      ENDDO


      IF (LTIMER) THEN                                                  
        CALL TIMER('FCDCH   ',4)                                        
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               


!!  Arguments:--------------------------------------------------------- 
      SUBROUTINE FCDCH_LAND(
     & P_FIELD,LAND_FIELD,TILE_PTS,TILE_INDEX,LAND_INDEX,
     & RIB,DB,VSHR,Z0M,Z0H,Z0F,ZH,Z1_UV,Z1_TQ,WIND_PROFILE_FACTOR,      
     & CDV,CHV,CDV_STD,V_S,V_S_STD,RECIP_L_MO,LTIMER                    
     &)                                                                 
      IMPLICIT NONE                                                     
                                                                        
      INTEGER                                                           
     & P_FIELD            ! IN Size of field on p-grid.
     &,LAND_FIELD         ! IN Number of land points.
     &,TILE_PTS           ! IN Number of tile points.
     &,TILE_INDEX(LAND_FIELD)
!                         ! IN Index of tile points.
     &,LAND_INDEX(P_FIELD)! IN Index of land points.
                                                                        
      LOGICAL                                                           
     & LTIMER                                                           
                                                                        
                                                                        
      REAL                                                              
     & DB(LAND_FIELD)! IN Buoyancy difference between surface and lowest
!                    !    temperature and humidity level in the         
!                    !    atmosphere (m/s^2).                           
     &,VSHR(P_FIELD) ! IN Wind speed difference between the surface and 
!                    !    the lowest wind level in the atmosphere (m/s).
     +,Z0M(LAND_FIELD)! IN Roughness length for momentum transport (m).
     +,Z0H(LAND_FIELD)! IN Roughness length for heat and moisture (m).
     +,ZH(P_FIELD)   ! IN Depth of boundary layer (m).                  
     +,Z1_UV(P_FIELD)! IN Height of lowest wind level (m).              
     +,Z1_TQ(P_FIELD)! IN Height of lowest temperature and humidity     
!                    !    level (m).                                    
     &,WIND_PROFILE_FACTOR(LAND_FIELD)
!                    ! IN for adjusting the surface transfer            
!                    !    coefficients to remove form drag effects.     
                                                                        
      REAL                                                              
     & CDV(LAND_FIELD)! OUT Surface transfer coefficient for momentum
!                    !     including orographic form drag (m/s).        
     +,CHV(LAND_FIELD)! OUT Surface transfer coefficient for
!                    !     heat, moisture & other scalars (m/s).        
     &,CDV_STD(LAND_FIELD)
!                    ! OUT Surface transfer coefficient for momentum    
!                    !     excluding orographic form drag (m/s).        
     &,V_S(LAND_FIELD)! OUT Surface layer scaling velocity
!                    !     including orographic form drag (m/s).        
     &,V_S_STD(LAND_FIELD)
!                    ! OUT Surface layer scaling velocity               
!                    !     excluding orographic form drag (m/s).        
     &,RECIP_L_MO(LAND_FIELD)
!                    ! OUT Reciprocal of the Monin-Obukhov length
!                    !     (m^-1).

      REAL
     & RIB(LAND_FIELD)  ! DUMMY Used in 7A boundary layer scheme
     &,Z0F(LAND_FIELD)  ! DUMMY Used in 7A boundary layer scheme

                                                                        
!*L  Workspace usage----------------------------------------------------
!                                                                       
!     Local work arrays.                                                
!                                                                       
      REAL                                                              
     & PHI_M(LAND_FIELD)! Monin-Obukhov stability function for momentum
!                     ! integrated to the model's lowest wind level.    
     &,PHI_H(LAND_FIELD)! Monin-Obukhov stability function for scalars
!                     ! integrated to the model's lowest temperature    
!                     ! and humidity level.                             
!                                                                       
!*----------------------------------------------------------------------
                                                                        
      EXTERNAL TIMER,PHI_M_H_LAND
                                                                        
!*----------------------------------------------------------------------
!  Common and local constants.                                          
*CALL C_VKMAN                                                           
      REAL BETA,THIRD                                                   
      PARAMETER (                                                       
     & BETA=0.08,   ! Tunable parameter in the surface layer scaling    
!                   ! velocity formula (multiplying the turbulent       
!                   ! convective scaling velocity).                     
     + THIRD=1./3.  ! One third.                                        
     +)                                                                 
      INTEGER N_ITS ! Number of iterations for Monin-Obukhov length     
!                   ! and stability functions.                          
      PARAMETER (                                                       
     & N_ITS=5                                                          
     &)                                                                 
!                                                                       
!  Define local variables                                               
!                                                                       
      INTEGER I,J,L ! Loop counter; horizontal field index.             
      INTEGER IT    ! Iteration loop counter.                           
                                                                        
      REAL                                                              
     & B_FLUX       ! Surface bouyancy flux over air density.           
     &,U_S          ! Surface friction velocity (effective value).      
     &,U_S_STD      ! Surface friction velocity (standard value).       
     &,W_S          ! Surface turbulent convective scaling velocity.    
                                                                        
      IF (LTIMER) THEN                                                  
        CALL TIMER('FCDCH   ',3)                                        
      ENDIF                                                             
                                                                        
!                                                                       
!-----------------------------------------------------------------------
!! 1. Set initial values for the iteration.                             
!-----------------------------------------------------------------------
!                                                                       
      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        I = LAND_INDEX(L)
        IF (DB(L) .LT. 0.0 .AND. VSHR(I) .LT. 2.0) THEN
!-----------------------------------------------------------------------
!         Start the iteration from the convective limit.
!-----------------------------------------------------------------------
          RECIP_L_MO(L) = -VKMAN/(BETA*BETA*BETA*ZH(I))
        ELSE
!-----------------------------------------------------------------------
!         Start the iteration from neutral values.
!-----------------------------------------------------------------------
          RECIP_L_MO(L) = 0.0
        ENDIF
        ENDDO
      CALL PHI_M_H_LAND (P_FIELD,LAND_FIELD,TILE_PTS,
     &                   TILE_INDEX,LAND_INDEX,
     &                   RECIP_L_MO,Z1_UV,Z1_TQ,Z0M,Z0H,
     &                   PHI_M,PHI_H,
     &                   LTIMER)
!                                                                       
      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        I = LAND_INDEX(L)
        IF (DB(L) .LT. 0.0 .AND. VSHR(I) .LT. 2.0) THEN
!-----------------------------------------------------------------------
!         Start the iteration from the convective limit.
!-----------------------------------------------------------------------
          V_S_STD(L) = BETA *
     &        SQRT( BETA * ( VKMAN / PHI_H(L) ) * ZH(I) * (-DB(L)) )
          V_S(L) = V_S_STD(L)
        ELSE
!-----------------------------------------------------------------------
!         Start the iteration from neutral values.
!-----------------------------------------------------------------------
          V_S(L) = ( VKMAN / PHI_M(L) ) * VSHR(I)
          V_S_STD(L) = V_S(L) * WIND_PROFILE_FACTOR(L)
        ENDIF
        CHV(L) = ( VKMAN / PHI_H(L) ) * V_S_STD(L)
        CDV(L) = ( VKMAN / PHI_M(L) ) * V_S(L)
        CDV_STD(L) = CDV(L) * ( V_S_STD(L) / V_S(L) ) *
     &                        WIND_PROFILE_FACTOR(L)
      ENDDO                                                             
!-----------------------------------------------------------------------
!! 2. Iterate to obtain sucessively better approximations for CD & CH.  
!-----------------------------------------------------------------------
      DO IT = 1,N_ITS                                                   
!                                                                       
        DO J=1,TILE_PTS
          L = TILE_INDEX(J)
          I = LAND_INDEX(L)
          B_FLUX = -CHV(L) * DB(L)
          U_S = SQRT( CDV(L) * VSHR(I) )
          U_S_STD = SQRT( CDV_STD(L) * VSHR(I) )
          IF (DB(L) .LT. 0.0) THEN
            W_S = (ZH(I) * B_FLUX)**THIRD
            V_S(L) = SQRT(U_S*U_S + BETA*BETA*W_S*W_S)
            V_S_STD(L) = SQRT(U_S_STD*U_S_STD + BETA*BETA*W_S*W_S)
          ELSE
            V_S(L) = U_S
            V_S_STD(L) = U_S_STD
          ENDIF
          RECIP_L_MO(L) = -VKMAN * B_FLUX /
     &                     (V_S(L)*V_S(L)*V_S(L))
        ENDDO                                                           
        CALL PHI_M_H_LAND (P_FIELD,LAND_FIELD,TILE_PTS,
     &                     TILE_INDEX,LAND_INDEX,
     &                     RECIP_L_MO,Z1_UV,Z1_TQ,Z0M,Z0H,
     &                     PHI_M,PHI_H,
     &                     LTIMER)
!                                                                       
                                                                        
        DO J=1,TILE_PTS
          L = TILE_INDEX(J)
          I = LAND_INDEX(L)
          CHV(L) = ( VKMAN / PHI_H(L) ) * V_S_STD(L)
          CDV(L) = ( VKMAN / PHI_M(L) ) * V_S(L)
          CDV_STD(L) = CDV(L) * ( V_S_STD(L) / V_S(L) ) *
     &                          WIND_PROFILE_FACTOR(L)
        ENDDO                                                           
      ENDDO ! Iteration loop
                                                                        
!-----------------------------------------------------------------------
!! Set CD's and CH's to be dimensionless paremters
!-----------------------------------------------------------------------
      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        I = LAND_INDEX(L)
        CDV(L) = CDV(L) / VSHR(I)
        CDV_STD(L) = CDV_STD(L) / VSHR(I)
        CHV(L) = CHV(L) / VSHR(I)
      ENDDO

*DECLARE FCDCH7A
*D FCDCH7A.46,FCDCH7A.48   
      SUBROUTINE FCDCH_SEA (P_POINTS,P_FIELD,P1,FLANDG,              
     &                      RIB,DB,VSHR,Z0M,Z0H,Z0F,                    
     &                      ZH,Z1_UV,Z1_TQ,CD,CH,                       
     &                      V_S,RECIP_L_MO,LTIMER)                      
*D FCDCH7A.59
*D FCDCH7A.62
     & FLANDG(P_FIELD)                                          
!                         ! IN Land fraction                         
     &,RIB(P_FIELD)       ! IN Bulk Richardson number.                  
*I FCDCH7A.75    
      REAL                                                              
     & DB(P_FIELD)   ! DUMMY Used in 6A boundary layer scheme           
     &,VSHR(P_FIELD) ! DUMMY Used in 6A boundary layer scheme           
     &,ZH(P_FIELD)   ! DUMMY Used in 6A boundary layer scheme           
     &,V_S(P_FIELD)  ! DUMMY Used in 6A boundary layer scheme           
     &,RECIP_L_MO(P_FIELD)                                              
!                    ! DUMMY Used in 6A boundary layer scheme           
                                                                        
                                                                        
*D FCDCH7A.121
        IF ( FLANDG(I).LT.1.0 ) THEN                                  
*D FCDCH7A.176
        ENDIF ! SEA_MASK                                              
*D FCDCH7A.194,FCDCH7A.195  
     & RIB,DB,VSHR,Z0M,Z0H,Z0F,ZH,Z1_UV,Z1_TQ,WIND_PROFILE_FACTOR,      
     & CD,CH,CD_STD,V_S,V_S_STD,RECIP_L_MO,LTIMER                       
*I FCDCH7A.229   
                                                                        
      REAL                                                              
     & DB(LAND_FIELD)       ! DUMMY Used in 6A boundary layer scheme    
     &,VSHR(P_FIELD)        ! DUMMY Used in 6A boundary layer scheme    
     &,ZH(P_FIELD)          ! DUMMY Used in 6A boundary layer scheme    
     &,V_S(LAND_FIELD)      ! DUMMY Used in 6A boundary layer scheme    
     &,V_S_STD(LAND_FIELD)  ! DUMMY Used in 6A boundary layer scheme    
     &,RECIP_L_MO(LAND_FIELD)                                           
!                           ! DUMMY Used in 6A boundary layer scheme    
                                                                        
*DECLARE FILL3A
*I ADB1F402.149   
!       5.3             25-04-01   Gather land, sea and                 
!                                  sea-ice temperatures and             
!                                  land fraction. Replace TFS           
!                                  with general sea temperature.        
!                                       (N. Gedney)                     
*D ADB1F401.205,ADB1F401.206  
     &   , PSTAR, TSTAR, TSTAR_SOLID, TSTAR_SEA
     &   , AB, BB, AC, BC, PEXNER, TAC                    
     &   , P, T, T_BDY, T_SURFACE, T_SOLID, T_SEA, D_MASS               
*I FILL3A.466   
     &   , TSTAR_SOLID(NPD_FIELD)                                       
!             SOLID SURFACE TEMPERATURES                                
     &   , TSTAR_SEA(NPD_FIELD)                                         
!             OPEN SEA SURFACE TEMPERATURES                             
*I ADB1F401.208   
     &   , T_SOLID(NPD_PROFILE)                                         
!             GATHERED TEMPERATURE OF SOLID SURFACE                     
     &   , T_SEA(NPD_PROFILE)                                           
!             GATHERED OPEN SEA TEMPERATURE                             
*I ADB1F401.213   
            T_SOLID(L)=TSTAR_SOLID(LG)                                  
            T_SEA(L)=TSTAR_SEA(LG)                                      
*DECLARE FTSA1A
*D ARE2F404.294,ARE2F404.304  
*D ARE2F404.307,AJG1F405.36   
     &     LAND, FLANDG, AICE,                                          
     &     TSTAR, TSTAR_SICE, SFA, MDSA, COSZ, S,                       
     &     ALPHAM,ALPHAC,ALPHAB,DTICE,L_MOSES_II,L_SSICE_ALBEDO,        
*D ARE2F404.309
     &     L1, L2, SA_LAND, SA_SICE, SAOS)                              
*D ARE2F404.310,ARE2F404.312  
*D ARE2F404.313,ARE2F404.314  
     &    ,L_MOSES_II                     ! .TRUE. if MOSES II land     
!                                         ! surface is selected.        
*I ADB1F400.16    
     &     FLANDG(L1),                    ! Land fraction               
*D ARE2F404.315
*D ARE2F404.316
     &     TSTAR_SICE(L1),                ! Seaice surface temperature  
*D ARE2F404.317,ARE2F404.318  
*D ARE2F404.319,FTSA1A.50   
     &     SA_LAND(L1),                   ! Surface Albedos for Land.   
!                                         ! (Not output for MOSESII).   
     &     SA_SICE(L1),                   ! Surface Albedos for seaice  
*D FTSA1A.58,FTSA1A.60   
      REAL DSA                           ! Deep-snow albedo (alphasubD) 
                                                                        
*D ARE2F404.321,ARE2F404.328  
     &     SNOW_ALBEDO,                   ! Snow albedo                 
*D ARE2F404.329
      PARAMETER ( MASKD = 0.2 )                                         
*I AWA1F304.1407  
      DO J=1, L2                                                        
        SA_LAND(J)=0.0                                                  
        SA_SICE(J)=0.0                                                  
      ENDDO                                                             
*D ARE2F404.330,ARE2F404.334  
! Land surface albedos are calculated by routine TILEALB if MOSES II    
! is selected                                                           
      IF (.NOT. L_MOSES_II) THEN                                        
*D ARE2F404.336
*D ARE2F404.337,ARE2F404.353  
          IF ( TSTAR(J) .LT. TCLAND ) THEN                              
             DSA = MDSA(J)                                              
           ELSE                                                         
            DSA=MDSA(J) + KLAND * (SFA(J)-MDSA(J)) * (TSTAR(J)-TCLAND)  
*D ARE2F404.355,ARE2F404.384  
          SA_LAND(J) = SFA(J) + (DSA-SFA(J)) * ( 1. - EXP(-MASKD*S(J)) )
          SAOS(J,1) = SA_LAND(J)                                        
          SAOS(J,2) = SA_LAND(J)                                        
*D ARE2F404.387,ARE2F404.407  
      ENDIF                                                             
*D ARE2F404.410
        IF (FLANDG(J).LT.1.0) THEN                                      
*D FTSA1A.103
             SA_SICE(J) = 0.                                            
*D FTSA1A.105,FTSA1A.106  
*D AJG1F405.63,AJG1F405.64   
                 if (tstar_sice(j).gt.tcice) then                       
                   snow_albedo=ice2+ice1*tstar_sice(j)                  
*D AJG1F405.68
                 sa_sice(j)=alphab                                      
*D AJG1F405.71
                 sa_sice(j)=alphab                                      
*D FTSA1A.108,FTSA1A.109  
             IF ( TSTAR_SICE(J) .LT. TCICE ) THEN                       
                SA_SICE(J) = ALPHAC                                     
*D FTSA1A.111
                SA_SICE(J) = ICE1 * TSTAR_SICE(J) + ICE2                
*I FTSA1A.114   
                                                                        
*I FTSA1A.115   
                                                                        
*DECLARE FXCA3A
*I ADB2F404.557   
!       4.6             10-05-98                Land flag added to      
!                                               the argument list for   
!                                               diagnostics.            
!                                               (J. M. Edwards)         
*I FXCA3A.17    
!                                                                       
*D FXCA3A.38
     &   , P, T, T_GROUND, T_SOLID, T_SEA, T_LEVEL, D_MASS              
*D ADB1F401.433,FXCA3A.80   
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION             
     &   , ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR, FLANDG                      
*D ADB1F401.435
     &   , FLUX_BELOW_690NM_SURF, FL_SEA_BELOW_690NM_SURF               
     &   , L_MOSES_II, L_CTILE                                          
     &   , SURF_VIS_DIR, SURF_VIS_DIF, SURF_NIR_DIR, SURF_NIR_DIF       
*I FXCA3A.232   
     &   , T_SOLID(NPD_PROFILE)                                         
!             TEMPERATURE OF SOLID SURFACE                              
     &   , T_SEA(NPD_PROFILE)                                           
!             SURFACE TEMPERATURE OVER OPEN SEA                         
*I ADB1F401.438   
     &   , L_MOSES_II                                                   
!             SURFACE SW REQUIRED FOR MOSES II                          
     &   , L_CTILE                                                      
!             COASTAL TILING SWITCH                                     
*I FXCA3A.458   
                                                                        
      REAL     !, INTENT(IN)                                            
     &     FLANDG(NPD_PROFILE)                                          
!            Land fraction                                              
*D ADB1F401.440,ADB1F401.443  
     &     N_FRAC_SOL_POINT                                             
!             NUMBER OF POINTS WITH FRACTIONAL ICE/LAND COVER           
     &   , I_FRAC_SOL_POINT(NPD_PROFILE)                                
!             INDICES OF POINTS WITH FRACTIONAL ICE/LAND COVER          
*D ADB1F401.450
!          GRID BOX MEAN NET SURFACE FLUX BELOW 690NM                   
     &   , FL_SEA_BELOW_690NM_SURF(NPD_PROFILE)                         
!          OPEN SEA NET SURFACE FLUX BELOW 690NM                        
     &   , SURF_VIS_DIR(NPD_PROFILE)                                    
!             DOWNWARD SURFACE DIRECT BEAM VISIBLE FLUX                 
     &   , SURF_VIS_DIF(NPD_PROFILE)                                    
!             DOWNWARD SURFACE DIFFUSE VISIBLE FLUX                     
     &   , SURF_NIR_DIR(NPD_PROFILE)                                    
!             DOWNWARD SURFACE DIRECT BEAM NEAR-INFRARED FLUX           
     &   , SURF_NIR_DIF(NPD_PROFILE)                                    
!             DOWNWARD SURFACE DIFFUSE NEAR-INFRARED FLUX               
*I ADB1F401.454   
     &   , PLANCK_LEADS_SEA(NPD_PROFILE)                                
!             PLANCK FUNCTION OVER SEA LEADS                            
*D ADB1F401.462
     &   , L_FLUX_BELOW_690NM_SURF                                      
     &   , FLUX_BELOW_690NM_SURF, FL_SEA_BELOW_690NM_SURF               
     &   , L_MOSES_II                                                   
     &   , SURF_VIS_DIR, SURF_VIS_DIF, SURF_NIR_DIR, SURF_NIR_DIF       
*I FXCA3A.1090  
     &         , T_SOLID, T_SEA, L_CTILE                                
*D ADB1F401.464,ADB1F401.465  
     &         , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION       
     &         , FLANDG, PLANCK_FREEZE_SEA, PLANCK_LEADS_SEA            
*D FXCA3A.1515
         CALL R2_COUPLE_DIAG(N_PROFILE, ISOLIR                          
*D ADB1F401.467,ADB1F401.468  
     &      , FLANDG, ICE_FRACTION                                      
     &      , PLANCK_FREEZE_SEA, PLANCK_LEADS_SEA                       
*D FXCA3A.1519,ADB1F401.470  
     &      , FLUX_TOTAL_BAND(1, 2*N_LAYER+2)                           
     &      , FLUX_TOTAL_BAND(1, 2*N_LAYER+1)                           
*D FXCA3A.1522,ADB1F401.471  
     &      , FLUX_TOTAL_CLEAR_BAND(1, 2*N_LAYER+2)                     
     &      , FLUX_TOTAL_CLEAR_BAND(1, 2*N_LAYER+1)                     
*D ADB1F401.473
     &      , L_FLUX_BELOW_690NM_SURF                                   
     &      , FLUX_BELOW_690NM_SURF, FL_SEA_BELOW_690NM_SURF            
     &      , L_MOSES_II, L_CTILE                                       
     &      , SURF_VIS_DIR, SURF_VIS_DIF, SURF_NIR_DIR, SURF_NIR_DIF    
*D FXCA3A.1549
*DECLARE GROWT2A
*D ABX1F405.1629
*D GROWT2A.58,GROWT2A.62   
*D GROWT2A.92
*DECLARE HYDCON5A
*I HYDCON5A.23    
     &,                   DK_DTHK
*I HYDCON5A.69    
     &,DK_DTHK(NPNTS)   ! OUT The rate of change of K with THETAK
!                       !     (kg/m2/s).
*I HYDCON5A.81    
        DK_DTHK(I)=0.0
*I HYDCON5A.83    
          DK_DTHK(I)=(2*B(I)+3)*KS(I)*(THETAK(I)**(2*B(I)+2))     
*DECLARE HYDROL7A
*D HYDROL7A.24,HYDROL7A.25   
     &                   LICE_PTS,LICE_INDEX,SOIL_PTS,SOIL_INDEX,
     &                   NTILES,TILE_PTS,TILE_INDEX,
*D HYDROL7A.27,HYDROL7A.35   
     &                   E_CANOPY,EXT,HCAP,HCON,INFIL_TILE,LS_RAIN,     
     &                   LS_SNOW,MELT_TILE,SATCON,SATHH,
     &                   SURF_HT_FLUX,TSTAR_TILE,FRAC, 
     &                   TIMESTEP,V_SAT,V_WILT,L_SNOW_ALBEDO,
     &                   STF_HF_SNOW_MELT,STF_SUB_SURF_ROFF,            
     &                   CAN_WCNT,RGRAIN,SMCL,SNOW_TILE,STHF,STHU,TSOIL,
     &                   CAN_WCNT_GB,HF_SNOW_MELT,INFIL,LYING_SNOW,SMC,
     &                   SNOW_MELT,SNOMLT_SUB_HTF,
*D HYDROL7A.48,HYDROL7A.51   

*D HYDROL7A.54,HYDROL7A.56   
*D HYDROL7A.101
     &,NTILES              ! IN Number of tiles.                        
*D HYDROL7A.118,HYDROL7A.119  
     &,TILE_PTS(NTILES)    ! IN Number of tile points.                  
     &,TILE_INDEX(NPNTS,NTILES)                                         
*D HYDROL7A.124,HYDROL7A.126  
     &,CAN_CPY(NPNTS,NTILES)!IN Canopy/surface capacity of              
!                          !    land tiles (kg/m2).           
*D HYDROL7A.129,HYDROL7A.130  
     &,E_CANOPY(NPNTS,NTILES)                                          
!                          ! IN Canopy evaporation from
*I HYDROL7A.135   
     &,INFIL_TILE(NPNTS,NTILES)                                         
!                          ! IN Maximum surface infiltration            
*I HYDROL7A.137   
     &,MELT_TILE(NPNTS,NTILES)
!                          ! IN Snowmelt on tiles (kg/m2/s).
*D HYDROL7A.141,HYDROL7A.148  
     &,SURF_HT_FLUX(NPNTS) ! IN Net downward surface heat flux (W/m2) 
     &,TSTAR_TILE(NPNTS,NTILES)
!                          ! IN Surface temperature (K).                
     &,FRAC(NPNTS,NTILES)  ! IN Tile fractions.
*D HYDROL7A.159,HYDROL7A.160  
     & CAN_WCNT(NPNTS,NTILES)                                          
!                          ! INOUT Canopy water content for
*D HYDROL7A.162,HYDROL7A.163  
     &,RGRAIN(NPNTS,NTILES)! INOUT Snow grain size (microns).           
*D HYDROL7A.166
     &,SNOW_TILE(NPNTS,NTILES)
!                          ! INOUT Snowmass on tiles (kg/m2). 
*D HYDROL7A.174
*I HYDROL7A.179   
     &,HF_SNOW_MELT(NPNTS)  ! OUT Gridbox snowmelt heat flux (W/m2).    
*I HYDROL7A.181   
     &,LYING_SNOW(NPNTS)    ! OUT Gridbox snowmass (kg/m2). 
*D HYDROL7A.202,HYDROL7A.208  
*D HYDROL7A.212,HYDROL7A.214  
*D HYDROL7A.217
     & SFSNOW,SURF_HYD,SOIL_HYD,SOIL_HTC,ICE_HTC,SOILMC
*D HYDROL7A.226
! Add snowfall to snow mass and update snow grain size            
*D HYDROL7A.228
      CALL SFSNOW(NPNTS,NTILES,TILE_PTS,TILE_INDEX,
     &            CON_SNOW,LS_SNOW,FRAC,TSTAR_TILE,TIMESTEP,
     &            RGRAIN,SNOW_TILE,L_SNOW_ALBEDO,LTIMER)
*D HYDROL7A.231
! Calculate the gridbox-mean snow mass and melt rate 
*D HYDROL7A.233,HYDROL7A.237  
      DO I=1,NPNTS  
        LYING_SNOW(I) = 0.
        SNOW_MELT(I) = 0. 
      ENDDO               
      DO N=1,NTILES                                                     
        DO J=1,TILE_PTS(N)                                              
          I = TILE_INDEX(J,N)
          LYING_SNOW(I) = LYING_SNOW(I) + FRAC(I,N)*SNOW_TILE(I,N)
          SNOW_MELT(I) = SNOW_MELT(I) + FRAC(I,N)*MELT_TILE(I,N)
        ENDDO                                                           
      ENDDO  
*I HYDROL7A.244   
          SNOMLT_SUB_HTF(I) = 0. 
*D HYDROL7A.247,HYDROL7A.256  
*D HYDROL7A.262,HYDROL7A.263  
      CALL SURF_HYD (NPNTS,NTILES,TILE_PTS,TILE_INDEX, 
*D HYDROL7A.265
     &               MELT_TILE,SNOW_MELT,TIMESTEP,
*I HYDROL7A.271   
      IF (SOIL_PTS.NE.0) THEN 
*D HYDROL7A.274
     &               SUB_SURF_ROFF,SMCL,STHU,SURF_ROFF,W_FLUX,
*I HYDROL7A.275   

      ELSE 

!-----------------------------------------------------------------------
! If required by STASH flag and there are no soil points, 
! set sub-surface runoff to zero.
!-----------------------------------------------------------------------
        IF(STF_SUB_SURF_ROFF) THEN   
          DO I=1,NPNTS  
            SUB_SURF_ROFF(I)=0.0
          ENDDO 
        ENDIF   

      ENDIF    
*D HYDROL7A.281,HYDROL7A.283  
        CALL SOIL_HTC (NPNTS,NSHYD,NTILES,SOIL_PTS,SOIL_INDEX,
     &                 TILE_PTS,TILE_INDEX,B,DZSOIL,FRAC,HCAP,HCON,
     &                 SATHH,SNOW_TILE,SURF_HT_FLUX,TIMESTEP,V_SAT,   
*D HYDROL7A.292
     &                SURF_HT_FLUX,TIMESTEP,                           
*D HYDROL7A.294,HYDROL7A.298  
*DECLARE HYDR_CT1
*D ARE1F404.113,ARE1F404.114  
     &           NTILESDA,TILE_FIELDDA,TILE_PTS,TILE_INDEX,
     &           CAN_EVAP_TILE,MELT_TILE,TILE_FRAC,  
*I ARE1F404.118   
     &       NTILESDA,                                                  
*D ARE1F404.122,ARE1F404.125  
     &       CAN_EVAP_TILE(TILE_FIELDDA,NTILESDA),                      
     &       MELT_TILE(TILE_FIELDDA,NTILESDA),
     &       TILE_FRAC(TILE_FIELDDA,NTILESDA)
*D ARE2F404.413
*I ACB1F304.20    
     &     ,PSLEVEL     !  loop counter for pseudolevels             
     &     ,PSLEVEL_OUT !  index for pseudolevels sent to STASH
*I ABX1F405.937   
     &,     PLLTILE(NTILESDA) ! pseudolevel list for surface tiles
*D ARE2F404.414,ARE2F404.415  
*D ARE2F404.416
     &   D1(JCANOPY_WATER),D1(JRGRAIN_TYP),L_SNOW_ALBEDO,SNODEP_LAND,
*D ARE1F404.142,ARE1F404.146  
     &   NTILES,TILE_PTS,TILE_INDEX,                                    
     &   D1(JCATCH_TYP),CAN_EVAP_TILE,                                  
     &   TILE_FRAC,D1(JINFIL_TYP),MELT_TILE,           
     &   D1(JTSTAR_TYP),D1(JCAN_WATER_TYP),D1(JSNODEP_TYP),
*D ARE2F404.417,ARE2F404.419  
*I AJS1F401.951   

CL ITEM 236 SNOW AMOUNT ON TILES                             
                                                                        
      IF (SF(236,8)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,236,8,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)                                            
        IF (ICODE.GT.0) THEN                                            
          RETURN                                                        
        END IF                                                          
        PSLEVEL_OUT=0                                                   
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1                                   
            CALL FROM_LAND_POINTS (                                     
     &          STASHWORK(SI(236,8,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),D1(JSNODEP_TYP+((PSLEVEL-1)*LAND_FIELD)),    
     &           D1(JLAND),P_FIELD,LAND_FIELD)
          END IF                                                        
        END DO                                                          
      END IF    

CL ITEM 237 SNOW MELT RATE ON TILES                             
                                                                        
      IF (SF(237,8)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,237,8,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)                                            
        IF (ICODE.GT.0) THEN                                            
          RETURN                                                        
        END IF                                                          
        PSLEVEL_OUT=0                                                   
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1                                   
            CALL FROM_LAND_POINTS (                                     
     &          STASHWORK(SI(237,8,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),MELT_TILE(1,PSLEVEL),        
     &           D1(JLAND),P_FIELD,LAND_FIELD)                          
          END IF                                                        
        END DO                                                          
      END IF

CL ITEM 238 SNOW GRAIN SIZE ON TILES                             
                                                                        
      IF (SF(238,8)) THEN                                               
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         
     &       STLIST(1,STINDEX(1,238,8,im_index)),                       
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              
     &       ICODE,CMESSAGE)                                            
        IF (ICODE.GT.0) THEN                                            
          RETURN                                                        
        END IF                                                          
        PSLEVEL_OUT=0                                                   
        DO PSLEVEL=1,NTILES                                             
          IF (PLLTILE(PSLEVEL)) THEN                                    
            PSLEVEL_OUT=PSLEVEL_OUT+1                                   
            CALL FROM_LAND_POINTS (                                     
     &          STASHWORK(SI(238,8,im_index)+(PSLEVEL_OUT-1)            
     &           *P_FIELD),D1(JRGRAIN_TYP+((PSLEVEL-1)*LAND_FIELD)),    
     &           D1(JLAND),P_FIELD,LAND_FIELD)                          
          END IF                                                        
        END DO                                                          
      END IF
                         
*DECLARE HYD_IC7A
*D HYD_IC7A.71
     &     NTILES,TILE_PTS,TILE_INDEX,                                  
*D HYD_IC7A.73,HYD_IC7A.74   
     &     FRAC,INFIL_TILE,MELT_TILE,TSTAR_TILE,       
     &     CAN_WCNT_TILE,SNOW_TILE,
*D HYD_IC7A.86,HYD_IC7A.87   
*D HYD_IC7A.102,HYD_IC7A.103  
     &   NTILES,                  ! IN Number of tiles. 
     &   TILE_PTS(NTILES),        ! IN Number of tile points.           
     &   TILE_INDEX(LAND,NTILES)  ! IN Index of tile points.            
*D HYD_IC7A.109,HYD_IC7A.113  
     &  CAN_CPY_TILE(LAND,NTILES),! IN Canopy/surface capacity of       
!                                 !    land tiles (kg/m2)     
     &  CAN_EVAP_TILE(LAND,NTILES),!IN Canopy evaporation from
*D HYD_IC7A.119
     &  FRAC(LAND,NTILES),        ! IN Tile fractions                   
*I HYD_IC7A.121   
     &  INFIL_TILE(LAND,NTILES),  ! IN Maximum surface infiltration     
!                                 !    rate for tiles (kg/m2/s)     
*I HYD_IC7A.123   
     &  MELT_TILE(LAND,NTILES),   ! IN Surface snowmelt on tiles
!                                 !    (kg/m2/s)       
*D HYD_IC7A.127,HYD_IC7A.133  
     &  SURF_HT_FLUX(LAND),       ! IN Net downward surface heat flux
!                                 !    (W/m2)
     &  TSTAR_TILE(LAND,NTILES),  ! IN Tile surface temperatures (K)    
*D HYD_IC7A.152,HYD_IC7A.156  
     &  CAN_WCNT_TILE(LAND,NTILES),!INOUT Canopy water content for      
!                                 !       land tiles (Kg/m2)  
     &  RGRAIN(LAND,NTILES),      ! INOUT Snow grain size (microns) 
     &  SNOW_TILE(LAND,NTILES),   ! INOUT Snowmass on tiles (kg/m2).
*D HYD_IC7A.163,HYD_IC7A.164  
     &  T_DEEP_SOIL(LAND,ST_LEVELS)!INOUT Deep soil temps. (K) 
*D HYD_IC7A.176
     &  SNODEP(LAND),             ! OUT Snow depth (Kg of water)        
*I HYD_IC7A.207   
     &  SNOW_SUBLIMATION(LAND),
*D HYD_IC7A.210
*D HYD_IC7A.229
     &     LICE_PTS,LICE_INDEX,SOIL_PTS,SOIL_INDEX,NTILES,              
*D HYD_IC7A.233,HYD_IC7A.241  
     &     EXT,HCAP,HCON,INFIL_TILE,LS_RAIN,LS_SNOW,MELT_TILE,
     &     SATCON,SATHH,
     &     SURF_HT_FLUX,TSTAR_TILE,                      
     &     FRAC,TIMESTEP,VSAT,VWILT,L_SNOW_ALBEDO,
     &     STF_HF_SNOW_MELT,STF_SUB_SURF_ROFF,                        
     &     CAN_WCNT_TILE,RGRAIN,SMCL,                  
     &     SNOW_TILE,STHF,STHU,T_DEEP_SOIL,  
     &     CAN_WCNT_GB,HF_SNOW_MELT,INFIL,SNODEP,SMC,SNOW_MELT, 
*D HYD_IC7A.243
*DECLARE IMCLTQ7A
*D IMCLTQ7A.47
     &                      LAND_INDEX,NTILES,TILE_INDEX,TILE_PTS,      
*D IMCLTQ7A.49
     &                      ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,        
*D IMCLTQ7A.53
     &                      SNOW_TILE,TILE_FRAC,                        
*D IMCLTQ7A.70,IMCLTQ7A.72   
     &,NTILES                      ! IN Number of land surface tiles.   
     &,TILE_PTS(NTILES)            ! IN Number of tiles.                
     &,TILE_INDEX(LAND_FIELD,NTILES)!IN Index of tile points.           
*D IMCLTQ7A.75
     & ALPHA1(LAND_FIELD,NTILES)   ! IN Gradient of saturated specific  
*D IMCLTQ7A.81,IMCLTQ7A.84   
!                                  !    heat flux into sea-ice (W/m2/K).
     &,ASHTF_TILE(LAND_FIELD,NTILES)!IN Coefficient to calculate heat
!                                  !    flux into land tiles (W/m2/K).  
*D IMCLTQ7A.91
     &,RESFT(LAND_FIELD,NTILES)    ! IN Total resistance factor         
*D IMCLTQ7A.94
     &,RHOKH_1(LAND_FIELD,NTILES)  ! IN Land surface exchange coeff.    
*D IMCLTQ7A.97
     &,RHOKPM(LAND_FIELD,NTILES)   ! IN Land surface exchange coeff.    
*D IMCLTQ7A.99,IMCLTQ7A.100  
     &,SNOW_TILE(LAND_FIELD,NTILES)! IN Lying snow on land tiles (kg/m2)
     &,TILE_FRAC(LAND_FIELD,NTILES)! IN Tile fractions.                 
*D IMCLTQ7A.111
     &,FQW_TILE(LAND_FIELD,NTILES) ! INOUT Surface flux of QW for land  
*D IMCLTQ7A.117
     &,FTL_TILE(LAND_FIELD,NTILES) ! INOUT H/Cp for land tiles.         
*I IMCLTQ7A.157   
     &,LAT_HT   ! Latent heat of evaporation for snow-free land 
!               ! or sublimation for snow-covered land and ice.
*D IMCLTQ7A.224,IMCLTQ7A.225  
! Land tiles                                                  
      DO N=1,NTILES                                                    
*D IMCLTQ7A.232
     &                      RESFT(L,N)*RHOKPM(L,N) *
     &                      (CP*RHOKH_1(L,N) + ASHTF_TILE(L,N))
*D IMCLTQ7A.234,IMCLTQ7A.244  
*D IMCLTQ7A.274,IMCLTQ7A.275  
! Land tiles                                                  
      DO N=1,NTILES                                                    
*D IMCLTQ7A.279
          LAT_HT = LC
          IF (SNOW_TILE(L,N).GT.0.) LAT_HT = LS   
          CT1(I) = CT1(I) - LAT_HT*GAMMA(1)*DTRDZ(I,1)*TILE_FRAC(L,N) * 
*D IMCLTQ7A.282,IMCLTQ7A.283  
     &          RHOKPM(L,N)*( LAT_HT*ALPHA1(L,N)*RESFT(L,N)*RHOKH_1(L,N)
     &                                               + ASHTF_TILE(L,N) )
*D IMCLTQ7A.285,IMCLTQ7A.296  
*D IMCLTQ7A.400,IMCLTQ7A.401  
! Land tiles                                                  
      DO N=1,NTILES                                                    
*I IMCLTQ7A.404   
          LAT_HT = LC
          IF (SNOW_TILE(L,N).GT.0.) LAT_HT = LS
*D IMCLTQ7A.406,IMCLTQ7A.407  
     &                             ( LAT_HT*RESFT(L,N)*RHOKH_1(L,N) *
     &                               (DQW(I,1) - ALPHA1(L,N)*DTL(I,1)) 
     &                                      - ASHTF_TILE(L,N)*DTL(I,1) )
*D IMCLTQ7A.409,IMCLTQ7A.410  
     &                    RHOKPM(L,N)*( CP*RHOKH_1(L,N) * 
     &                               (DQW(I,1) - ALPHA1(L,N)*DTL(I,1)) 
     &                                      + ASHTF_TILE(L,N)*DQW(I,1) )
*D IMCLTQ7A.414,IMCLTQ7A.428  
*DECLARE INITIAL1
*D ABX1F404.224
      IF (L_VEG_FRACS .AND. STEPim(a_im) == 0 ) THEN
*DECLARE INITVEG1
*I INITVEG1.47    
!   4.6   16/02/99   Correct initialisation of accumulated leaf turnover
!                    rate when both TRIFFID and phenology in use.
!                    Richard Betts
*D INITVEG1.88,INITVEG1.90   
*D INITVEG1.125,INITVEG1.136  
      CALL SPARM (LAND_FIELD,LAND1,LAND_PTS,NTILES,TILE_PTS,TILE_INDEX, 
     &            D1(JFRAC_TYP),D1(JCANHT_PFT),           
     &            D1(JLAI_PFT),D1(JSAT_SOIL_COND),
     &            D1(JCATCH_TYP),D1(JINFIL_TYP),D1(JZ0_TYP))
*D INITVEG1.187
      ENDIF
*I INITVEG1.188   
      IF (L_PHENOL) THEN
*D ABX1F405.72,ABX1F405.75   
        DO L = 1,LAND_FIELD
          D1(JG_LF_PFT_ACC+L-1) = 0.0
        ENDDO
*DECLARE LEAF7A
*D LEAF7A.174
        RD(L) = FDC3 * VCM(L)   
*D LEAF7A.189
        RD(L) = FDC3 * VCM(L)  
*D LEAF7A.465
        RD(L) = FDC4 * VCM(L)  
*D LEAF7A.481
        RD(L) = FDC4 * VCM(L)   
*DECLARE LOTKA2A
*D ABX1F405.1644
     &,                 C_VEG,FORW,FRAC_VS,FRAC_AGRIC,GAMMA,LAI,PC_S    
*I ABX1F405.1655  
     &,FRAC_AGRIC(LAND_FIELD)     ! IN Fraction of agriculture.    
*D LOTKA2A.53,ABX1F405.1651 
*I LOTKA2A.153   
C----------------------------------------------------------------------
C Exclude non-crop types from agricultural regions
C----------------------------------------------------------------------
*D LOTKA2A.157
          SPACE(L,N)=1.0-NOSOIL(L)-FRAC_AGRIC(L)*(1-CROP(N))
     &                            -FRAC_MIN*(NPFT-K)  
*D LOTKA2A.176,LOTKA2A.177  
          B(L,N) = PC_S(L,N)*SPACE(L,N)/C_VEG(L,N)-G_AREA(N)   
*DECLARE LWRAD3A
*I ADB2F404.631   
!       4.6             10-05-98                Land flag passed to     
!                                               FLUX_CALC.              
!                                               (J. M. Edwards)         
*I ASK1F405.289   
!                                                                       
*D LWRAD3A.27
     &   , TAC, PEXNER, TSTAR, TSTAR_SOLID, TSTAR_SEA, L_CTILE 
     &   , PSTAR, AB, BB, AC, BC                    
*D LWRAD3A.33
     &   , LAND, FLANDG, ICE_FRACTION                                   
*I LWRAD3A.219   
     &   , L_CTILE                                                      
!             COASTAL TILING SWITCH                                     
      REAL      !, INTENT(IN)                                           
     &     FLANDG(NPD_PROFILE)                                          
!            Land fraction                                              
*I LWRAD3A.224   
     &   , TSTAR_SOLID(NPD_FIELD)                                       
!             SOLID SURFACE TEMPERATURE                                 
     &   , TSTAR_SEA(NPD_FIELD)                                         
!             OPEN SEA SURFACE TEMPERATURES                             
*D LWRAD3A.226
!             SEA ICE FRACTION OF SEA PORTION OF GRID BOX               
*I LWRAD3A.321   
     &   , T_SOLID(NPD_PROFILE)                                         
!             GATHERED TEMPERATURE OF SOLID SURFACE                     
     &   , T_SEA(NPD_PROFILE)                                           
!             GATHERED OPEN SEA TEMPERATURE                             
*D ADB1F401.537,ADB1F401.540  
     &     N_FRAC_SOL_POINT                                             
!             NUMBER OF POINTS WITH FRACTIONAL ICE/LAND COVER           
     &   , I_FRAC_SOL_POINT(NPD_PROFILE)                                
!             INDICES OF POINTS WITH FRACTIONAL ICE/LAND COVER          
*D ADB1F401.556,ADB1F401.557  
     &   , PSTAR, TSTAR, TSTAR_SOLID, TSTAR_SEA 
     &   , AB, BB, AC, BC, PEXNER, TAC                    
     &   , P, T, T_BDY, T_SURFACE, T_SOLID, T_SEA, D_MASS               
*D ADB1F401.568
     &   , FLANDG                                                     
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION           
*D ADB1F401.574
     &   , P, T, T_SURFACE, T_SOLID, T_SEA, T_BDY, D_MASS               
*D ADB1F401.575,LWRAD3A.629  
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION           
     &   , ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR, FLANDG, LWSEA               
*D LWRAD3A.631
     &   , L_DUMMY, DUMMY, DUMMY, DUMMY                           
     &   , L_DUMMY, L_CTILE, DUMMY, DUMMY, DUMMY, DUMMY         
*D LWRAD3A.706
         IF (FLANDG(L).EQ.1.0.OR.ICE_FRACTION(L).EQ.1.0) THEN           
*D LWRAD3A.708
         ELSE IF (FLANDG(L).LT.TOL_TEST.AND.
     &     ICE_FRACTION(L).LT.TOL_TEST) THEN                     
*D ADB1F401.576,ADB1F401.577  
!           LWSEA MUST BE SCALED BY THE FRACTION OF OPEN SEA TO         
!           TOTAL SEA FOR CONSISTENCY WITH UPPER LEVELS IN THE MODEL.   
*D LWRAD3A.714
            LWOUT(L, 1)=LWOUT(L, 1)-(1.0E+00-FLANDG(L))*LWSEA(L)        
*D ADB1F401.579
     &   , FLANDG                                                       
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION             
*I ADB1F401.584   
!              FRACTION OF SEA-ICE IN SEA PART OF GRID BOX.             
     &   , FLANDG(NPD_PROFILE)                                          
!                  LAND FRACTION                                        
                                                                        
*D ADB1F401.587,ADB1F401.590  
     &     N_FRAC_SOL_POINT                                             
!             NUMBER OF POINTS WITH FRACTIONAL ICE/LAND COVER           
     &   , I_FRAC_SOL_POINT(NPD_PROFILE)                                
!             INDICES OF POINTS WITH FRACTIONAL ICE/LAND COVER          
*D ADB1F401.598,ADB1F401.599  
!     SET THE FRACTIONAL OPEN SEA COVERAGE. POINTS ARE REQUIRED WHERE   
!     THIS IS NEITHER 0 NOR 1.                                          
*D ADB1F401.601
       SEARCH_ARRAY(L)=(1.0E+00-FLANDG(L))*(1.0E+00-ICE_FRACTION(L))    
       SEARCH_ARRAY(L)=SEARCH_ARRAY(L)*(1.-SEARCH_ARRAY(L))             
*D GSS2F402.246
      N_FRAC_SOL_POINT=0                                                
*D GSS2F402.249,GSS2F402.250  
          N_FRAC_SOL_POINT                  =N_FRAC_SOL_POINT+1         
          I_FRAC_SOL_POINT(N_FRAC_SOL_POINT)=L                          
*DECLARE NAMSIZE
*I RB300993.113   
     & ,NTILES
*DECLARE NVEGPARM
*D NVEGPARM.16,NVEGPARM.18   
      DATA CATCH_NVG   /  0.50,  0.00,  0.00,  0.00 /    
      DATA GS_NVG      /  0.00,  0.00,  1E-2,   1E6 /
      DATA INFIL_NVG   /  0.10,  0.00,  0.50,  0.00 /
*DECLARE PHIMH6A
*D PHIMH6A.2
*IF DEF,A03_8A                                                          
*D PHIMH6A.21
!!!   SUBROUTINES PHI_M_H_SEA AND PHI_M_H_LAND ------------------------
*D PHIMH6A.39,PHIMH6A.40   
      SUBROUTINE PHI_M_H_SEA(
     & P_POINTS,P_FIELD,P1,LAND_MASK,
*D PHIMH6A.52
*D PHIMH6A.84
      INTEGER I       ! Loop counter; horizontal field index.           
*D PHIMH6A.104
      DO I=P1,P1+P_POINTS-1                                             
*D PHIMH6A.106,PHIMH6A.111  
        IF ( .NOT. LAND_MASK(I) ) THEN
*D PHIMH6A.113,PHIMH6A.115  
*D PHIMH6A.162
        ENDIF  ! LAND_MASK
                                                                        
      ENDDO                                                             
                                                                        
      IF (LTIMER) THEN                                                  
        CALL TIMER('PHI_M_H ',4)                                        
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               

!!!                                                                     
!*L  Arguments:---------------------------------------------------------
      SUBROUTINE PHI_M_H_LAND(
     & P_FIELD,LAND_FIELD,TILE_PTS,TILE_INDEX,LAND_INDEX,
     & RECIP_L_MO,Z_UV,Z_TQ,Z0M,Z0H,PHI_M,PHI_H,LTIMER                  
     &)                                                                 
      IMPLICIT NONE                                                     
                                                                        
      INTEGER                                                           
     & P_FIELD            ! IN Size of field on p-grid.
     &,LAND_FIELD         ! IN Number of land points.
     &,TILE_PTS           ! IN Number of tile points.
     &,TILE_INDEX(LAND_FIELD)
!                         ! IN Index of tile points.
     &,LAND_INDEX(P_FIELD)! IN Index of land points.
                                                                        
      LOGICAL                                                           
     & LTIMER                                                           
                                                                        
                                                                        
      REAL                                                              
     & RECIP_L_MO(LAND_FIELD)
!                       ! IN Reciprocal of the Monin-Obukhov length 
!                       !     (m^-1).
     &,Z_UV(P_FIELD)    ! IN Height of wind level above roughness
!                       !    height (m)
     &,Z_TQ(P_FIELD)    ! IN Height of temperature, moisture and scalar
!                       !    lev above the roughness height (m).        
     &,Z0M(LAND_FIELD)  ! IN Roughness length for momentum (m).
     &,Z0H(LAND_FIELD)  ! IN Roughness length for heat/moisture/scalars
!                       !    (m)
!                                                                       
      REAL                                                              
     & PHI_M(LAND_FIELD)! OUT Stability function for momentum.
     &,PHI_H(LAND_FIELD)! OUT Stability function for
!                       !     heat/moisture/scalars.
!                                                                       
!*L  Workspace usage----------------------------------------------------
!    No work areas are required.                                        
!                                                                       
!*----------------------------------------------------------------------
!*L  External subprograms called:                                       
                                                                        
      EXTERNAL TIMER                                                    
                                                                        
!*----------------------------------------------------------------------
!  Common and local physical constants.                                 
!                                                                       
!  None.                                                                
!                                                                       
!  Define local variables.                                              
!                                                                       
      INTEGER I,J,L     ! Loop counter; horizontal field index.
!                                                                       
      REAL                                                              
     & PHI_MN         ! Neutral value of stability function for momentum
     &,PHI_HN         ! Neutral value of stability function for scalars.
     &,ZETA_UV        ! Temporary in calculation of PHI_M.              
     &,ZETA_0M        ! Temporary in calculation of PHI_M.              
     &,ZETA_TQ        ! Temporary in calculation of PHI_H.              
     &,ZETA_0H        ! Temporary in calculation of PHI_H.              
     &,X_UV_SQ        ! Temporary in calculation of PHI_M.              
     &,X_0M_SQ        ! Temporary in calculation of PHI_M.              
     &,X_UV           ! Temporary in calculation of PHI_M.              
     &,X_0M           ! Temporary in calculation of PHI_M.              
     &,Y_TQ           ! Temporary in calculation of PHI_H.              
     &,Y_0H           ! Temporary in calculation of PHI_H.              
                                                                        
      IF (LTIMER) THEN                                                  
        CALL TIMER('PHI_M_H ',3)                                        
      ENDIF                                                             
                                                                        
      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        I = LAND_INDEX(L)
!                                                                       
!-----------------------------------------------------------------------
!! 1. Calculate neutral values of PHI_M and PHI_H.                      
!-----------------------------------------------------------------------
!                                                                       
        PHI_MN = LOG( (Z_UV(I) + Z0M(L)) / Z0M(L) )
        PHI_HN = LOG( (Z_TQ(I) + Z0M(L)) / Z0H(L) )
!                                                                       
!-----------------------------------------------------------------------
!! 2. Calculate stability parameters.                                   
!-----------------------------------------------------------------------
!                                                                       
        ZETA_UV = (Z_UV(I) + Z0M(L)) * RECIP_L_MO(L)
        ZETA_TQ = (Z_TQ(I) + Z0M(L)) * RECIP_L_MO(L)
        ZETA_0M = Z0M(L) * RECIP_L_MO(L)
        ZETA_0H = Z0H(L) * RECIP_L_MO(L)
!                                                                       
!-----------------------------------------------------------------------
!! 3. Calculate PHI_M and PHI_H for neutral and stable conditions.      
!-----------------------------------------------------------------------
!                                                                       
        IF (RECIP_L_MO(L) .GE. 0.0) THEN
          PHI_M(L) = PHI_MN + 4.0 * (ZETA_UV - ZETA_0M)
          PHI_H(L) = PHI_HN +
     &               (1.0 + 2.0*ZETA_TQ) * (1.0 + 2.0*ZETA_TQ) -
     &               (1.0 + 2.0*ZETA_0H) * (1.0 + 2.0*ZETA_0H)
!                                                                       
!-----------------------------------------------------------------------
!! 4. Calculate PHI_M and PHI_H for unstable conditions.                
!-----------------------------------------------------------------------
!                                                                       
        ELSE
                                                                        
          X_UV_SQ = SQRT(1.0 - 16.0*ZETA_UV)
          X_0M_SQ = SQRT(1.0 - 16.0*ZETA_0M)
          X_UV = SQRT(X_UV_SQ)
          X_0M = SQRT(X_0M_SQ)
          PHI_M(L) = PHI_MN - 2.0*LOG( (1.0+X_UV) / (1.0+X_0M) )
     &                    - LOG( (1.0+X_UV_SQ) / (1.0+X_0M_SQ) )
     &                    + 2.0*( ATAN(X_UV) - ATAN(X_0M) )
                                                                        
          Y_TQ = SQRT(1.0 - 16.0*ZETA_TQ)
          Y_0H = SQRT(1.0 - 16.0*ZETA_0H)
          PHI_H(L) = PHI_HN - 2.0*LOG( (1.0+Y_TQ) / (1.0+Y_0H) )
                                                                        
        ENDIF
*DECLARE PHYSIO7A
*D PHYSIO7A.33
     &,                   P_FIELD,NSHYD,NTILES,TILE_PTS,TILE_INDEX      
*D PHYSIO7A.36,ABX1F405.839  
     &,                   V_CRIT,V_SAT,V_WILT,WIND,Z0_TILE,Z1
     &,                   CANHC_TILE,VFRAC_TILE,FLAKE,G_LEAF,GS,GS_TILE
     &,                   GPP,GPP_FT,NPP,NPP_FT,RESP_P,RESP_P_FT       
     &,                   RESP_S,RESP_W_FT,SMCT,WT_EXT_TILE)          
*I PHYSIO7A.56    
     &,NTILES                     ! IN Number of surface tiles.
*D PHYSIO7A.67
     &,FRAC(LAND_FIELD,NTYPE)     ! IN Surface type fractions.          
*D PHYSIO7A.78
     &,TSTAR_TILE(LAND_FIELD,NTILES)                                    
*D PHYSIO7A.91
     &,Z0_TILE(LAND_FIELD,NTILES) ! IN Tile roughness lengths (m).
*D ABX1F405.840,PHYSIO7A.98   
     & CANHC_TILE(LAND_FIELD,NTILES)
!                                 ! OUT Areal heat capacity of canopy
!                                 !     for land tiles (J/K/m2).
     &,FLAKE(LAND_FIELD,NTILES)   ! OUT Lake fraction.
     &,G_LEAF(LAND_FIELD,NPFT)    ! OUT Leaf turnover rate (/360days).  
     &,GS_TILE(LAND_FIELD,NTILES) ! OUT Surface conductance for
*D PHYSIO7A.113
     &,VFRAC_TILE(LAND_FIELD,NTILES)
!                                 ! OUT Fractional canopy coverage for
!                                 !     land tiles.
     &,WT_EXT_TILE(LAND_FIELD,NSHYD,NTILES)
!                                 ! OUT Fraction of evapotranspiration  
*D PHYSIO7A.115
!                                 !     soil layer by each tile.
*D PHYSIO7A.118
     & CANHC(LAND_FIELD)          ! WORK Canopy heat capacity (J/K/m2).
     &,CH_TYPE(LAND_FIELD,NTYPE)  ! WORK CANHC for surface types. 
     &,F_ROOT(NSHYD)              ! WORK Fraction of roots in each soil 
*I PHYSIO7A.120   
     &,GSOIL(LAND_FIELD)          ! WORK Bare soil conductance.
     &,GS_TYPE(LAND_FIELD,NTYPE)  ! WORK Conductance for surface types.
*I PHYSIO7A.126   
     &,TSTAR(LAND_FIELD)          ! WORK Surface temperature (K).
     &,VFRAC(LAND_FIELD)          ! WORK Fractional canopy coverage. 
     &,VF_TYPE(LAND_FIELD,NTYPE)  ! WORK VFRAC for surface types. 
     &,WT_EXT(LAND_FIELD,NSHYD)   ! WORK Gridbox-mean WT_EXT.  
     &,WT_EXT_TYPE(LAND_FIELD,NSHYD,NTYPE) 
!                                 ! WORK WT_EXT for surface types. 
     &,Z0(LAND_FIELD)             ! WORK Roughness length (m).
*I PHYSIO7A.149   
          DO N=1,NTYPE 
            WT_EXT_TYPE(L,K,N)=0.0
          ENDDO                        
*D PHYSIO7A.174
          GS_TYPE(L,N)=GS(L)
*I PHYSIO7A.185   
        CANHC(L)=0.0
        VFRAC(L)=0.0
      ENDDO

      DO L=LAND1,LAND1+LAND_PTS-1 
        GSOIL(L) = 0.
        IF (V_CRIT(L).GT.0.)                     
     &    GSOIL(L) = GS_NVG(SOIL-NPFT)*(STHU(L,1)*V_SAT(L)/V_CRIT(L))**2
*I PHYSIO7A.194   
        IF (NTILES.EQ.1) THEN                                           
          DO L=1,LAND_FIELD
            TSTAR(L) = TSTAR_TILE(L,1)
            Z0(L) = Z0_TILE(L,1)
          ENDDO 
        ELSE
          DO L=1,LAND_FIELD
            TSTAR(L) = TSTAR_TILE(L,N)
            Z0(L) = Z0_TILE(L,N)
          ENDDO
        ENDIF                                                      

*D PHYSIO7A.198,PHYSIO7A.199  
     &,               F_ROOT,STHU,V_CRIT,V_SAT,V_WILT 
     &,               WT_EXT_TYPE(1,1,N),FSMC)
*D PHYSIO7A.203
     &,             RIB,WIND,Z0,Z0,Z1,RA)                   
*D PHYSIO7A.208
     &,               Q1,RA,TSTAR
*D PHYSIO7A.210
     &,               RESP_W_FT(1,N),GS_TYPE(1,N))

        CALL SOIL_EVAP (LAND_FIELD,NSHYD,TILE_PTS(N),TILE_INDEX(1,N)
     &,                 GSOIL,LAI(1,N),GS_TYPE(1,N),WT_EXT_TYPE(1,1,N))
*D PHYSIO7A.213
     &,                N,FSMC,TSTAR,G_LEAF(1,N)) 

        CALL CANCAP (LAND_FIELD,TILE_PTS(N),TILE_INDEX(1,N),N
     &,              HT(1,N),LAI(1,N),CH_TYPE(1,N),VF_TYPE(1,N))
*D PHYSIO7A.218,PHYSIO7A.219  
! Non-vegetated surface types 
*D PHYSIO7A.221,PHYSIO7A.228  
      DO N=NPFT+1,NTYPE                                             
*D PHYSIO7A.231
          GS_TYPE(L,N) = GS_NVG(N-NPFT)
          DO K=1,NSHYD 
            WT_EXT_TYPE(L,K,N) = 0.                   
          ENDDO                                        
*I PHYSIO7A.232   
      ENDDO 
*I PHYSIO7A.233   
! Copy soil conductance and add bare soil fraction to extraction from
! surface layer
      N = SOIL
      DO J=1,TILE_PTS(N)                                              
        L=TILE_INDEX(J,N)                                             
        GS_TYPE(L,N) = GSOIL(L)
        WT_EXT_TYPE(L,1,N) = 1.                          
      ENDDO  

!----------------------------------------------------------------------
! Canopy heat capacity and coverage set to 0 for non-vegetated surfaces 
!---------------------------------------------------------------------- 
      DO N=NPFT+1,NTYPE                                               
        DO J=1,TILE_PTS(N)                                              
          L=TILE_INDEX(J,N)
          CH_TYPE(L,N) = 0.
          VF_TYPE(L,N) = 0.
        ENDDO
*D PHYSIO7A.245

      IF (NTILES.EQ.1) THEN
        DO N=1,NTYPE
          DO J=1,TILE_PTS(N)                                            
            L=TILE_INDEX(J,N)                                           
            CANHC(L) = CANHC(L) + FRAC(L,N)*CH_TYPE(L,N) 
            GS(L) = GS(L) + FRAC(L,N)*GS_TYPE(L,N)                     
            VFRAC(L) = VFRAC(L) + FRAC(L,N)*VF_TYPE(L,N)
            DO K=1,NSHYD               
              WT_EXT(L,K) = WT_EXT(L,K) + FRAC(L,N)*WT_EXT_TYPE(L,K,N)
            ENDDO                                  
          ENDDO                                                         
        ENDDO
        DO L=LAND1,LAND1+LAND_PTS-1
          FLAKE(L,1) = FRAC(L,7)
          GS_TILE(L,1) = 0.
          IF (FLAKE(L,1).LT.1.)
     &      GS_TILE(L,1) = GS(L) / (1. - FLAKE(L,1))
          CANHC_TILE(L,1) = CANHC(L)                                    
          VFRAC_TILE(L,1) = VFRAC(L)
          DO K=1,NSHYD
            WT_EXT_TILE(L,K,1) = WT_EXT(L,K)
          ENDDO                     
        ENDDO
      ELSE
        DO N=1,NTYPE                                                    
          DO J=1,TILE_PTS(N)                                            
            L=TILE_INDEX(J,N)
            FLAKE(L,N) = 0.
            GS_TILE(L,N) = GS_TYPE(L,N) 
            CANHC_TILE(L,N) = CH_TYPE(L,N)                              
            VFRAC_TILE(L,N) = VF_TYPE(L,N)
            DO K=1,NSHYD
              WT_EXT_TILE(L,K,N) = WT_EXT_TYPE(L,K,N)
            ENDDO                        
          ENDDO                                                         
        ENDDO 
        N = 7    ! Lake tile
*D PHYSIO7A.248
          FLAKE(L,N) = 1.
*D PHYSIO7A.250
      ENDIF
*DECLARE PPARM1A
*D PPARM1A.26,PPARM1A.28   
     &,                      HT,LAI,SATCON,CATCH_T,INFIL_T,Z0_T)
*D PPARM1A.43,PPARM1A.44   
     & HT(LAND_FIELD)             ! IN Vegetation height (m).           
*D PPARM1A.46,PPARM1A.47   
     &,SATCON(LAND_FIELD)         ! IN Saturated hydraulic conductivity
!                                 !    (kg/m2/s). 
*I PPARM1A.48    
     &,INFIL_T(LAND_FIELD)        ! OUT Maximum surface infiltration 
!                                 !     rate (kg/m2/s).         
*D PPARM1A.50
*D PPARM1A.64,PPARM1A.66   
*I PPARM1A.68    
        INFIL_T(L) = INFIL_F(N) * SATCON(L)                
*DECLARE PRELIM1
*I ABX1F405.2     
!   4.6    19/01/99    Use vegetation sampling frequencies in
!                      atmosphere timesteps, and set start times to
!                      first phenology or TRIFFID timestep.  Requires
!                      call to CNTLGEN to provide information on
!                      timestep length.  Richard Betts
*I GSS3F401.804   
*CALL CNTLGEN
*I PRELIM1.95    
! Local vectors:
      INTEGER SECS_PER_ASTEP ! Number of seconds in atmos timestep
      INTEGER A_PHENOL_STEP  ! Leaf phenology period in atmos timesteps
      INTEGER A_TRIFFID_STEP ! TRIFFID period in atmos timesteps

*I PRELIM1.102   
! 0.1  Convert PHENOL_PERIOD and TRIFFID_PERIOD to atmosphere timesteps.

      SECS_PER_ASTEP = FLOAT(SECS_PER_PERIODim(atmos_sm))/
     &                 FLOAT(STEPS_PER_PERIODim(atmos_sm))

      A_PHENOL_STEP = PHENOL_PERIOD*(86400.0/SECS_PER_ASTEP)
      A_TRIFFID_STEP = TRIFFID_PERIOD*(86400.0/SECS_PER_ASTEP)

*D ABX1F405.5,ABX1F405.6    
          ELSE IF((ITIMA.EQ.14).AND.(A_PHENOL_STEP.NE.1)) THEN
            IMD=MOD(LIST_S(st_start_time_code,NRECS),A_PHENOL_STEP)
*D ABX1F405.8
     &      LIST_S(st_start_time_code,NRECS)+A_PHENOL_STEP-IMD
*D ABX1F405.10,ABX1F405.11   
          ELSE IF((ITIMA.EQ.15).AND.(A_TRIFFID_STEP.NE.1)) THEN
            IMD=MOD(LIST_S(st_start_time_code,NRECS),A_TRIFFID_STEP)
*D ABX1F405.13
     &      LIST_S(st_start_time_code,NRECS)+A_TRIFFID_STEP-IMD
*D ABX1F405.19
              LIST_S(st_freq_code,NRECS)=A_PHENOL_STEP
*D ABX1F405.21
     &      (MOD(LIST_S(st_freq_code,NRECS),A_PHENOL_STEP).NE.0) THEN
*D ABX1F405.23
     &       'PRELIM: INCORRECT SAMPLING FOR A_PHENOL_STEP . FREQ=',
*D ABX1F405.32
              LIST_S(st_freq_code,NRECS)=A_TRIFFID_STEP
*D ABX1F405.34
     &      (MOD(LIST_S(st_freq_code,NRECS),A_TRIFFID_STEP).NE.0) THEN
*D ABX1F405.36
     &       'PRELIM: INCORRECT SAMPLING FOR A_TRIFFID_STEP . FREQ=',
*D ABX1F405.46
              LIST_S(st_freq_code,NRECS)=A_PHENOL_STEP
*D ABX1F405.48
              LIST_S(st_freq_code,NRECS)=A_TRIFFID_STEP
*DECLARE PSLIMS1
*D ABX1F404.16,ABX1F404.17   
!Direct vegetation parametrization: all surface tiles
        ILAST=NTILES
*DECLARE RAD_CTL1
*D ARE2F404.72
     &             NTILESDA,TILE_FIELDDA,DOLR,LW_DOWN,SW_TILE,
     &             LAND_ALB,SICE_ALB,           
*I ARN1F404.123   
     &       NTILESDA,      ! and NTILES                                
*D ARE2F404.74,ARE2F404.75   
     &       SAL_DIM,     ! IN Set to P_FIELD for MOSES II,             
C                         !    1 otherwise                              
     &       TILE_FIELDDA,! IN Set to LAND_FIELD for MOSES II,          
C                         !    1 otherwise                              
*I GSS1F304.768   
     &      DOLR(P_FIELDDA),                                            
     &      LW_DOWN(P_FIELDDA),                                         
     &      SW_TILE(TILE_FIELDDA,NTILESDA),
     &      LAND_ALB(P_FIELDDA,NLALBS),   ! Mean land albedo            
     &      SICE_ALB(P_FIELDDA,NLALBS),   ! Mean sea-ice albedo         
*D ARE2F404.76,ARE2F404.81   
*I ARE2F404.83    
*CALL C_0_DG_C                                                          
*I ADB2F404.909   
     &     ,TILE_ALBEDO                                                 
*D ARE2F404.85,ARE2F404.86   
C zenith angle adjustment, net surface SW on tiles and downward LW      
     &      RADINCS((P_FIELDDA*(P_LEVELSDA+2+NTILESDA)+511)/512*512),   
*D ARE2F404.87,ARE2F404.88   
     &      LAND_ALBEDO(SAL_DIM,4),                                     
*I AWI1F403.143   
     &      FLANDG(P_FIELDDA),            ! Land fraction               
*D ARE2F404.89,ARE2F404.92   
     &     ,ALB_TILE(TILE_FIELDDA,NTILESDA,4)                           
     &     ,ICE_FRACT(P_FIELDDA)                                     
     &     ,TILE_FRAC(TILE_FIELDDA,NTYPE)                               
     &     ,TSTAR_TILE(TILE_FIELDDA,NTILESDA)                           
     &     ,SURF_DOWN_SW(P_FIELDDA,4)                                   
     &     ,T_SOL_RAD(P_FIELDDA)     ! Effective surface radiative temp
     &     ,TSTAR_SICE(P_FIELDDA)    ! Sea-ice sfc temperature (K)
     &     ,FRACSOLID(P_FIELDDA)     ! Solid surface fraction in gridbox
     &     ,SW_NET_LAND(P_FIELDDA)   ! SW net local flux over land
     &     ,SW_NET_SICE(P_FIELDDA)   ! SW net local flux over sea-ice  
     &     ,SW_NET_RTS(P_FIELDDA)    ! net SW on tiles
     &     ,SW_DOWN_RTS(P_FIELDDA)   ! down SW on tiles
*I RAD_CTL1.123   
                                                                        
! Index arrays for MOSES II                                             
      INTEGER                                                           
     & TILE_PTS(NTYPE)               ! Number of land points which      
                                     ! include the nth surface type     
     &,TILE_INDEX(TILE_FIELDDA,NTYPE)! Indices of land points which     
                                     ! include the nth surface type     
*I ADB1F401.772   
     &  ,L_SURF_DOWN_FLUX     !Logical to calculate surface downward    
!                             !fluxes                                   
*I ACN2F405.47    
     &  ,L_CTILE              !coastal tiling switch      
*I ADB1F400.70    
      REAL                                                              
     &       ALBSOLID    !Mean solid surface albedo                     
!                                                                       
*I ADB2F404.914   
      LOGICAL 
     &      LAND0P5(P_FIELDDA) ! LOCAL Land mask 
!                              !   (TRUE if land fraction >0.5)  

*D ARE2F404.96,ARE2F404.99   
!  Set index arrays and flags for MOSES II, and copy tile fractions     
!  and surface temperatures from D1                                     
      IF ( H_SECT(3) .EQ. "07A" .OR.                                    
     &     H_SECT(3) .EQ. "08A" ) THEN                                  
        CALL TILEPTS(P_FIELD,LAND_FIELD,LAND1,LAND_PTS,D1(JFRAC_TYP),   
     &               TILE_PTS,TILE_INDEX)                               
        DO N=1,NTYPE                                                    
          DO POINT=1,TILE_PTS(N)                                        
            L = TILE_INDEX(POINT,N)                                     
            J = (N-1)*LAND_FIELD + L - 1                                
            TILE_FRAC(L,N) = D1(JFRAC_TYP+J)                            
          ENDDO                                                         
        ENDDO                                                           
        DO N=1,NTILES                                                   
          DO L=LAND1,LAND1+LAND_PTS-1                                   
            J = (N-1)*LAND_FIELD + L - 1                                
            TSTAR_TILE(L,N) = D1(JTSTAR_TYP+J)                          
          ENDDO                                                         
        ENDDO                                                           
        L_MOSES_II=.TRUE.
        L_CTILE=.TRUE.               
        L_SURF_DOWN_FLUX=.TRUE.                                         
*D ARE2F404.101
        L_MOSES_II=.FALSE.                                              
        L_SURF_DOWN_FLUX=SF(235, 1)                                     
*D RAD_CTL1.173
*D ARE2F404.104,ARE2F404.105  
*I AJS1F401.961   
        SW_NET_LAND(I) = 0.0
        SW_NET_SICE(I) = 0.0    
        LAND_ALB(I,1) = 0.0                                         
*I RAD_CTL1.180   
!  Set up GLOBAL fractional land field:                                 
      CALL LAND_TO_GLOBAL                                               
     & (D1(JLAND),D1(JFRAC_LAND),FLANDG,LAND_PTS,P_FIELDDA)   
                                                          
!  Set up FRACSOLID (solid surface fraction in grid-box),
!  ICE_FRACT (ice fraction in grid-box) and
!  LAND0P5 (set to TRUE where land fraction greater than or equal
!  to 0.5) 
      DO I=1,P_FIELD
        ICE_FRACT(I) = D1(JICE_FRACTION+I-1)   
        FRACSOLID(I) = FLANDG(I) + (1.0-FLANDG(I))*ICE_FRACT(I)
        LAND0P5(I) = .FALSE.
      ENDDO     

      DO I=1,P_FIELD                                              
        IF (FLANDG(I).GE.0.5) THEN
          LAND0P5(I) = .TRUE.
        ENDIF                        
      ENDDO  
     
*D ARE2F404.108
        LEN=(P_FIELDDA*(P_LEVELS+2+NTILES)+511)/512*512                 
C                                                  !no words for SW incs
*D ARE2F404.109,ARE2F404.110  
*D RAD_CTL1.311,ABX1F405.141  
*D @DYALLOC.3042,ARE2F404.132  
     &      D1(JLAND+JS),FLANDG(FIRST_POINT),
     &      D1(JICE_FRACTION+JS),D1(JTSTAR+JS),D1(JTSTAR_SICE+JS),      
*D ARE2F404.133
*D AJG1F405.31
     &      ALPHAM,ALPHAC,ALPHAB,DTICE,L_MOSES_II,L_SSICE_ALBEDO,       
*D ARE2F404.134
*D ARE2F404.135,GHM5F405.5    
     &      LAND_ALB(FIRST_POINT,1),SICE_ALB(FIRST_POINT,1),            
*I GHM5F405.9     
            IF ( L_MOSES_II ) THEN                                      
!-----------------------------------------------------------------------
! Calculate MOSES II tile albedos, then reset TILE_PTS and TILE_INDEX   
! and set tile fractions to 1 if aggregate tiles are used (NTILES=1).   
!-----------------------------------------------------------------------
              CALL TILE_ALBEDO (                                        
     &          P_FIELD,LAND_FIELD,LAND1,LAND_PTS,LAND_LIST,NTILES,     
     &          TILE_PTS,TILE_INDEX,L_SNOW_ALBEDO,D1(JSOIL_ALB),        
     &          COS_ZENITH_ANGLE,TILE_FRAC,D1(JLAI_PFT),D1(JRGRAIN_TYP),
     &          D1(JSNODEP_TYP),TSTAR_TILE,D1(JZ0_TYP),                 
     &          ALB_TILE,LAND_ALBEDO )                                  
                                                                        
              IF (NTILES.EQ.1) THEN                                     
                TILE_PTS(1) = LAND_PTS                                  
                DO L=LAND1,LAND1+LAND_PTS-1                             
                  TILE_FRAC(L,1) = 1.                                   
                  TILE_INDEX(L+1-LAND1,1) = L                           
                ENDDO                                                   
              ENDIF                                                     
                                                                        
            ENDIF                                                       
*D ADB1F404.2
*D ADB1F404.4,ARE2F404.137  
     &       (H_SECT(3).EQ."07A").OR.                                   
     &       (H_SECT(3).EQ."08A") ) THEN                                
*D ARE2F404.139
*D ARE2F404.140
         DO LEVEL=0,P_LEVELS+1+NTILES                                   
*D ARE2F404.142
            IF( L_MOSES_II ) THEN                                       
*D GHM5F405.76,ARE2F404.150  
     &        LAND_ALBEDO(FIRST_POINT_SAL,1), L_CTILE,                  
     &        LAND_ALB(FIRST_POINT,1), SICE_ALB(FIRST_POINT,1),
     &        FLANDG(FIRST_POINT), OPEN_SEA_ALBEDO(FIRST_POINT,1),    
     &        D1(JICE_FRACTION+JS),D1(JLAND+JS),LAND0P5(FIRST_POINT),
     &        D1(JSNODEP+JS),         
!                       MOSES II flag and array dimension               
     &        L_MOSES_II, SAL_DIM,                                      
*D ADB1F401.823,ADB1F400.198  
     &        STASHWORK(JS+SI(204,1,im_index)),
     &        STASHWORK(JS+SI(259,1,im_index)), 
     &        STASHWORK(JS+SI(260,1,im_index)), L_FLUX_BELOW_690NM_SURF,
     &        STASHWORK(JS+SI(235,1,im_index)), L_SURF_DOWN_FLUX,       
*I ADB2F404.1008  
     &        SURF_DOWN_SW(FIRST_POINT,1),                              
*I RAD_CTL1.472   
          DO BAND=1,4                                                   
            DO POINT = FIRST_POINT, LAST_POINT                          
              SURF_DOWN_SW(POINT,BAND) = 0.                             
            ENDDO                                                       
          ENDDO                                                         
*D ARE2F404.151
CL Set up downward surface SW components if required for MOSES II       
!   Calculate SW_NET_RTS, SW_DOWN_RTS and LAND_ALBEDO
*D ARE2F404.153,ARE2F404.155  
      IF ( L_MOSES_II ) THEN
*I ARE2F404.156   
        DO I=FIRST_POINT, LAST_POINT
          LAND_ALB(I,1) = 0.0
        ENDDO
        DO L=LAND1,LAND1+LAND_PTS-1
          I = LAND_LIST(L)
          SW_DOWN_RTS(I) = 0.0
          SW_NET_RTS(I) = 0.0
        END DO                                        
        DO N=1,NTILES                                                   
          DO I=FIRST_POINT,LAST_POINT                                   
            J = I + (P_LEVELS + 1 + N)*P_FIELD                          
            RADINCS(J) = 0.                                             
          ENDDO                                                         
        ENDDO
               
        DO N=1,NTILES                                                 
          DO POINT=1,TILE_PTS(N)                                      
            L = TILE_INDEX(POINT,N)                                   
            I = LAND_LIST(L)                                          
            J = I + (P_LEVELS + 1 + N)*P_FIELD        
            DO BAND=1,4                                                 
              RADINCS(J) = RADINCS(J) + (1. - ALB_TILE(L,N,BAND)) *     
     &                                              SURF_DOWN_SW(I,BAND)
            ENDDO
            SW_NET_RTS(I)=SW_NET_RTS(I)+RADINCS(J)*TILE_FRAC(L,N)
          ENDDO                                                         
        ENDDO

        IF (L_CTILE) THEN
          DO L=LAND1,LAND1+LAND_PTS-1
            I = LAND_LIST(L)
            DO BAND=1,4
              SW_DOWN_RTS(I) = SW_DOWN_RTS(I) + SURF_DOWN_SW(I,BAND)
            ENDDO
            IF (SW_DOWN_RTS(I).GT.0.0) THEN
              LAND_ALB(I,1) = 1.0-SW_NET_RTS(I)/SW_DOWN_RTS(I)
            ENDIF 
          ENDDO
        ENDIF
                  
      ENDIF                                                        
*D ARE2F404.157,ARE2F404.158  
CL and net surface on tiles                                             
        LEN=(P_FIELDDA*(P_LEVELS+2+NTILES)+511)/512*512                 
C                                                                       
*D RAD_CTL1.622,RAD_CTL1.623  
      IF (L_CTILE) THEN
        DO I=FIRST_POINT,LAST_POINT                                     
          NET_ATM_FLUX(I) = NETSW(I)                                    
*D RAD_CTL1.626
!  Mutliply sea field by sea fraction. 
     &                   *(1.-FLANDG(I))       
        END DO 
      ELSE
        DO I=FIRST_POINT,LAST_POINT                                     
          NET_ATM_FLUX(I) = NETSW(I)                                    
     &                   - MEAN_COSZ(I) * RADINCS(I)                    
     &                   - STASHWORK(SI(203,1,im_index)+I-1)  
        ENDDO
      ENDIF                                                   
*D RAD_CTL1.695,RAD_CTL1.697  
! Only set on sea-ice points if MOSES2:                                 
      IF (L_CTILE) THEN
        DO I=FIRST_POINT,LAST_POINT                                     
*D ARE2F404.159,ARE2F404.160  
! land_alb is only set on points where there is incoming sw radiation   
! at the previous timestep, therefore it will be zero over some         
! land points                                                           
          IF (FRACSOLID(I).GT.0.0) THEN
*D ABX1F405.142
            IF (FLANDG(I).GT.0.0.AND.LAND_ALB(I,1).LE.0.0) THEN
              SW_NET_LAND(I) = 0.0
              SW_NET_SICE(I) = 0.0                                      
            ELSE                                         
              ALBSOLID = ( FLANDG(I) * LAND_ALB(I,1) +
     &          (1.0-FLANDG(I)) * ICE_FRACT(I) * SICE_ALB(I,1) )        
     &            /FRACSOLID(I)                                         
*D ABX1F405.144,ABX1F405.145  
              IF (FLANDG(I).GT.0.0) THEN
                SW_NET_LAND(I) = RADINCS(I)                             
     &              * COS_ZENITH_ANGLE(I) / FRACSOLID(I)                
     &              * (1.0-LAND_ALB(I,1))/(1.0-ALBSOLID)                
              ENDIF
            
              IF (ICE_FRACT(I).GT.0.0) THEN
                SW_NET_SICE(I) = RADINCS(I)                             
     &              * COS_ZENITH_ANGLE(I) / FRACSOLID(I)                
     &              * (1.0-SICE_ALB(I,1))/(1.0-ALBSOLID)                
                                                                        
                SURF_RADFLUX(I) = SW_NET_SICE(I) * ICE_FRACT(I)
              ENDIF                                                     
            ENDIF 
                                                                      
          ENDIF 
        ENDDO                                                           
      ELSE                                                              
*D ABX1F405.147,ABX1F405.160  
          SURF_RADFLUX(I) = RADINCS(I) * COS_ZENITH_ANGLE(I)            
        END DO                                                          
*I RAD_CTL1.698   
      IF ( L_MOSES_II ) THEN                                            
        DO N=1,NTILES                                                   
          DO POINT=1,TILE_PTS(N)                                        
            L = TILE_INDEX(POINT,N)                                     
            I = LAND_LIST(L)                                            
            J = I + (P_LEVELS + 1 + N)*P_FIELD                          
            SW_TILE(L,N) = RADINCS(J) * COS_ZENITH_ANGLE(I)             
          END DO                                                        
        END DO                                                          
      ENDIF                                                             
                                                                        
*D RAD_CTL1.705
       IF(SF(202,1)) THEN                                              
*D RAD_CTL1.712,RAD_CTL1.715  
*I RAD_CTL1.716   
          IF (L_CTILE) THEN                                      
            DO I = FIRST_POINT,LAST_POINT                               
              STASHWORK(SI(201,1,im_index)+I-1)=RADINCS(I)*MEAN_COSZ(I)+
!  Multiply sea field by sea fraction  
     &        STASHWORK(SI(203,1,im_index)+I-1)*(1.-FLANDG(I))          
            END DO                                                      
          ELSE
            DO I = FIRST_POINT,LAST_POINT                               
              STASHWORK(SI(201,1,im_index)+I-1)=RADINCS(I)*MEAN_COSZ(I)+
     &        STASHWORK(SI(203,1,im_index)+I-1)
            END DO
          ENDIF

        END IF 

C Calculate new diagnostics for down sw land and sea-ice portions: 
    
        IF(SF(257,1)) THEN                                              
          DO I=FIRST_POINT, LAST_POINT                                  
            STASHWORK(SI(257,1,im_index)+I-1) = SW_NET_LAND(I)          
          ENDDO                                                         
        END IF                                                          
                                                                        
        IF(SF(258,1)) THEN                                              
          DO I=FIRST_POINT, LAST_POINT                                  
            STASHWORK(SI(258,1,im_index)+I-1) = SW_NET_SICE(I)          
          ENDDO                                                         
*D ARE2F404.174
        OFFSET=(P_FIELDDA*(P_LEVELS+2+NTILES)+511)/512*512              
*D ARE2F404.175
        LEN=(P_FIELDDA*(P_LEVELS+2+NTILES)+511)/512*512  !no words for L
*D ARE2F404.178
          T_SOL_RAD(I) = 0.0 
          TSTAR_SICE(I) = D1(JTSTAR_SICE+I-1)
*D ARE2F404.180
        IF ( L_CTILE ) THEN
          DO I=1,P_FIELD                                            
            T_SOL_RAD(I) = (1.0-FLANDG(I))*                             
     &            ICE_FRACT(I)*TSTAR_SICE(I)**4                         
          ENDDO
          DO N=1,NTILES                                                 
            DO POINT=1,TILE_PTS(N)    
              L = TILE_INDEX(POINT,N)                                   
              I = LAND_LIST(L)                                          
              T_SOL_RAD(I) = T_SOL_RAD(I) + 
     &            FLANDG(I)*TILE_FRAC(L,N)*TSTAR_TILE(L,N)**4
            ENDDO
          ENDDO                     
          DO  I=1,P_FIELD
            IF (FRACSOLID(I).GT.0.0) THEN
              T_SOL_RAD(I)=(T_SOL_RAD(I)/FRACSOLID(I))**0.25            
            ENDIF                                                       
          ENDDO
                                                                        
        ELSEIF ( L_MOSES_II ) THEN                                      
*D ARE2F404.183,ARE2F404.184  
            T_SOL_RAD(I) = 0.                                           
*D ARE2F404.186,ABX1F405.163  
          DO N=1,NTILES                                                 
            DO POINT=1,TILE_PTS(N)                                      
              L = TILE_INDEX(POINT,N)                                   
*D ARE2F404.189,ARE2F404.191  
              T_SOL_RAD(I) = T_SOL_RAD(I) + TILE_FRAC(L,N) *            
     &                                      TSTAR_TILE(L,N)**4          
*D ARE2F404.196
            T_SOL_RAD(I) = T_SOL_RAD(I)**0.25                           
*I ARE2F404.197   
        ENDIF                                                           
                                                                        
        IF ( L_MOSES_II ) THEN                                          
           L_SURF_DOWN_FLUX=.TRUE.                                      
        ELSE                                                            
           L_SURF_DOWN_FLUX=SF(207, 2)                                  
*D ARE2F404.200
         DO LEVEL=0,P_LEVELS+1+NTILES                                   
*I APBBF401.69    
           LW_DOWN(I)=0.0                                               
*D ARE2F404.201
     &      D1(JP_EXNER(1)+JS_LOCAL(I)),T_SOL_RAD(FP_LOCAL(I)),         
*D ARE2F404.202
     &        D1(JP_EXNER(1)+JS_LOCAL(I)),D1(JTSTAR+JS_LOCAL(I)),
     &        T_SOL_RAD(FP_LOCAL(I)),D1(JTSTAR_SEA+JS_LOCAL(I)),
     &        L_CTILE,        
*D ADB1F400.322
C Only want the 0.5 threshold LAND mask and fractional land:            
     &        LAND0P5(FP_LOCAL(I)),FLANDG(JS_LOCAL(I)+1),           
*D ADB1F400.334
     &        LW_DOWN(FP_LOCAL(I)),L_SURF_DOWN_FLUX,                    
*D ARE2F404.203
C Store downward surface LW radiation flux if required for MOSES II     
C DOLR is TOA outward LW - surface upward LW for land and sea-ice       
*D ARE2F404.205,ARE2F404.207  
      IF ( L_MOSES_II ) THEN                                            
                                                                        
        DO I=FIRST_POINT,LAST_POINT
          DOLR(I) = OLR(I) 
C         If no land in grid box and some sea-ice...                    
          IF (FLANDG(I).EQ.0.0 .AND. ICE_FRACT(I).GT.0.0) THEN 
            DOLR(I) = DOLR(I) - ICE_FRACT(I)*SBCON*TSTAR_SICE(I)**4     
          ENDIF                                                         
        ENDDO                                                           
                                                                        
        DO L=LAND1,LAND1+LAND_PTS-1                                     
          I = LAND_LIST(L)                                              
          DOLR(I) = DOLR(I) - FRACSOLID(I)*SBCON*T_SOL_RAD(I)**4 
        ENDDO                                                           
                                                                        
        DO I=FIRST_POINT,LAST_POINT                                     
          J = I + (P_LEVELS + 2)*P_FIELD + OFFSET                       
          RADINCS(J) = LW_DOWN(I)                                       
          RADINCS(J+P_FIELD) = DOLR(I)                                  
        ENDDO                                                           
                                                                        
      ENDIF                                                             
*D RAD_CTL1.906,RAD_CTL1.907  
      IF (L_CTILE) THEN
        DO I=FIRST_POINT,LAST_POINT                                     
          NET_ATM_FLUX(I) = - OLR(I)                                    
*D RAD_CTL1.910
!  Mutliply sea field by sea fraction
     &                   *(1.-FLANDG(I))     
        END DO
      ELSE
          NET_ATM_FLUX(I) = - OLR(I)                                    
     &                   - RADINCS(I+OFFSET)                            
     &                   - STASHWORK(SI(203,2,im_index)+I-1)
      ENDIF                                                            
*I RAD_CTL1.940   
      IF (L_MOSES_II) THEN                                              
        DO I=FIRST_POINT,LAST_POINT                                     
          SURF_RADFLUX(I) = SURF_RADFLUX(I) +                           
     &                      D1(JICE_FRACTION+I-1)*LW_DOWN(I)            
        ENDDO                                                           
      ELSE                                                              
*I RAD_CTL1.943   
      ENDIF                                                             
*D ARE2F404.211,ARE2F404.235  
*D RAD_CTL1.964,RAD_CTL1.967  
          IF  (L_CTILE) THEN
            DO I = FIRST_POINT,LAST_POINT                               
                STASHWORK(SI(201,2,im_index)+I-1) = RADINCS(I+OFFSET)+
!  Multiply sea field by sea fraction    
     &          STASHWORK(SI(203,2,im_index)+I-1)*(1.-FLANDG(I)) 
            END DO 
          ELSE
            DO I = FIRST_POINT,LAST_POINT                               
                STASHWORK(SI(201,2,im_index)+I-1) = RADINCS(I+OFFSET)+  
     &          STASHWORK(SI(203,2,im_index)+I-1)
            END DO
          ENDIF                                                       
*I RAD_CTL1.1005  
        IF(SF(207,2)) THEN                                              
                                                                        
CL  Downward Surface radiative flux                                     
                                                                        
          CALL COPYDIAG (STASHWORK(SI(207,2,im_index)),LW_DOWN,         
     &        FIRST_POINT,LAST_POINT,P_FIELD,ROW_LENGTH,                
     &        im_ident,2,207,                                           
*CALL ARGPPX                                                            
     &        ICODE,CMESSAGE)                                           
                                                                        
          IF (ICODE .GT. 0) RETURN                                      
                                                                        
        END IF                                                          
                                                                        
*D ABX1F405.185
        CALL SWAPB_LAND(SW_TILE,LAND_FIELD,P_FIELD,                     
     &                  ROW_LENGTH,P_ROWS,EW_Halo,NS_Halo,              
     &                  NTILES,LAND_LIST)                               
        CALL SWAPBOUNDS(LW_DOWN,ROW_LENGTH,P_ROWS,                      
*D ABX1F405.187
        CALL SWAPBOUNDS(DOLR,ROW_LENGTH,P_ROWS,                         
*I ABX1F405.188   
      ENDIF
      IF (L_CTILE) THEN
        CALL SWAPBOUNDS(LAND_ALB(1,1),ROW_LENGTH,P_ROWS, 
     &                  EW_Halo,NS_Halo,NLALBS)                         
        CALL SWAPBOUNDS(SICE_ALB(1,1),ROW_LENGTH,P_ROWS,
     &                  EW_Halo,NS_Halo,NLALBS)        
*DECLARE READFL1A
*I GSI1F405.366   
      ELSE IF (LOOKUP(LBPACK,K).EQ.0.AND.
     &         ppxref_grid_type.EQ.ppx_atm_compressed) THEN
! Grid code suggests land-only field but field is unpacked,
! so use grid code for land data on full field
        fake_D1_ADDR(d1_grid_type)=ppx_atm_tland


*DECLARE ROOTFR7A
*D ROOTFR7A.72
      PARAMETER (P=1.0)
*DECLARE RPANCA1A
*I GDG0F401.1361  
     &                     FLAND_CTILE,                                 
     &                     TSTAR_LAND_CTILE,TSTAR_SEA_CTILE,            
     &                     TSTAR_SICE_CTILE,                            
*D RPANCA1A.111,RPANCA1A.113  
     &       ICE_FRACTION(P_FIELD), 
!                               !INOUT  Ice frac of sea part of grid
!                               !       box, updated if requested   
*IF -DEF,RECON                                                     
     &       FLAND_CTILE(LAND_FIELD),                                   
!                                !IN  Fractional land on land pts.   
*ENDIF                                                                  
     &       FLAND_G(P_FIELD),   !WORK Frac land over all points.       
     &       TSTAR(P_FIELD),     !INOUT  TSTAR:updated if requested     
     &       TSTAR_LAND_CTILE(P_FIELD),                                 
!                                !INOUT  as above, but for land.        
     &       TSTAR_SEA_CTILE(P_FIELD),                                  
!                                !INOUT  as above, but for open sea.    
     &       TSTAR_SICE_CTILE(P_FIELD),                                 
!                                !INOUT  as above, but for sea-ice.     
*D RPANCA1A.120
     &       LAND(P_FIELD),      ! WORK LAND mask                       
     &       SEA(P_FIELD),       ! WORK SEA mask                        
     &       LTSTAR_SICE         ! IN TRUE if TSTAR_SICE has been read i
!                                ! from input dump.                     
!                                ! If FALSE set to TSTAR_SEA.           
*I RPANCA1A.132   
*CALL CNTLATM
*I RPANCA1A.165   
     &,      TSTAR_LAND(P_FIELD)!Temporary store for land surface temp.
     &,      TSTAR_SEA(P_FIELD) !as above, but for open sea.           
     &,      TSTAR_SICE(P_FIELD)!as above, but for sea-ice.            
     &,      TSTAR_SSI(P_FIELD) !as above, but for sea mean.           
*I UDG4F402.255   
     &       L,                 ! Land index                        
*I GRS2F404.12    
     &       ,L_CTILE           ! Coastal tiling switch
C                               ! always true for 7A/8A BL
*I RPANCA1A.266   
!     Set coastal tiling flag
      IF ( H_SECT(3).EQ.'07A' .OR. H_SECT(3).EQ.'08A') THEN
        L_CTILE = .TRUE.
      ENDIF                                                 
                                                                        
!     Set up surface temperatures:                                      
                                                                        
       IF(L_CTILE)THEN                                                  
         DO I=1,P_FIELD                                                 
            TSTAR_LAND(I)=TSTAR_LAND_CTILE(I)                           
            TSTAR_SEA(I)=TSTAR_SEA_CTILE(I)                             
            TSTAR_SICE(I)=TSTAR_SICE_CTILE(I)                           
            IF(ICE_FRACTION(I).LE.0.0)THEN                              
              TSTAR_SSI(I)=TSTAR_SEA(I)                                 
            ELSE                                                        
              TSTAR_SSI(I)=ICE_FRACTION(I)*TSTAR_SICE(I)                
     &          +(1.0-ICE_FRACTION(I))*TSTAR_SEA(I)                     
            ENDIF                                                       
         ENDDO                                                          
       ELSE                                                             
         DO I=1,P_FIELD                                                 
            TSTAR_LAND(I)=TSTAR(I)                                      
            TSTAR_SSI(I)=TSTAR(I)                                       
         ENDDO                                                          
       ENDIF                                                            
                                                                        
                                                                        
*I RPANCA1A.332   
! Read in fractional land field                                         
*IF DEF,RECON                                                      
       IF(L_CTILE)THEN                                                  
        FILE=FILEANCIL(48)                                             
        NFTIN=FTNANCIL(FILE)                                            
        WRITE(6,*)'READING IN LAND FRACTION'                            
        CALL READFLDS(NFTIN,1,NLOOKUP(48),LOOKUP(1,LOOKUP_START(FILE)),
     &                LEN1_LOOKUP,ANCIL1,P_FIELD,FIXHD(1,FILE),         
*CALL ARGPPX                                              
     &                      ICODE,CMESSAGE)                             
        WRITE(6,*)'READ IN LAND FRACTION'                               
                                                                        
      DO I=1,P_FIELD                                                    
        FLAND_G(I)=0.0                                                  
        IF(LAND(I))FLAND_G(I)=ANCIL1(I)                                 
! If land or sea fraction is less than machine tolerance print warning  
          IF(LAND(I).AND.FLAND_G(I).LE.1.42E-14)THEN                
           WRITE(6,*)'*****************WARNING********************'     
           WRITE(6,*)'LAND FRACTION IS LESS THAN MACHINE TOLERANCE'     
          ENDIF                                                         
          IF(.NOT.LAND(I).AND.1.0-FLAND_G(I).LE.1.42E-14)THEN       
           WRITE(6,*)'*****************WARNING********************'     
           WRITE(6,*)'SEA FRACTION IS LESS THAN MACHINE TOLERANCE'      
          ENDIF                                                         
!                                                                       
          IF(FLAND_G(I).LE.0.0.AND.LAND(I))THEN                         
           WRITE(6,*)'*ERROR* a) LAND FRAC & LAND MASK ARE INCONSISTENT'
           ICODE = 800                                                  
           CMESSAGE='REPLANCA:ERROR:LAND FRAC & MASK ARE INCONSISTENT'  
          ENDIF                                                         
          IF(FLAND_G(I).GT.0.0.AND..NOT.LAND(I))THEN                    
           WRITE(6,*)'*ERROR* b) LAND FRAC & LAND MASK ARE INCONSISTENT'
           ICODE = 801                                                  
           CMESSAGE='REPLANCA:ERROR:LAND FRAC & MASK ARE INCONSISTENT'  
          ENDIF                                                         
         ENDDO                                                          
                                                                        
        ELSE                     ! Not coastal tiling:                  
         DO I=1,P_FIELD                                                 
          IF(LAND(I))THEN                                               
            FLAND_G(I)=1.0                                              
          ELSE                                                          
            FLAND_G(I)=0.0                                              
          ENDIF                                                         
         ENDDO                                                          
        ENDIF                                                           
*ENDIF                                                                 
*IF -DEF,RECON                                                     
! Set up global fractional land field                                   
         IF(L_CTILE)THEN                                                
           L=0                                                          
           DO I=1,P_FIELD                                               
             FLAND_G(I)=0.0                                             
             IF(LAND(I))THEN                                            
               L=L+1                                                    
               FLAND_G(I)=FLAND_CTILE(L)                                
            ENDIF                                                       
           ENDDO                                                        
         ELSE                                                           
           DO I=1,P_FIELD                                               
             IF(LAND(I))THEN                                            
               FLAND_G(I)=1.0                                           
             ELSE                                                       
               FLAND_G(I)=0.0                                           
            ENDIF                                                       
           ENDDO                                                        
         ENDIF                                                          
!                                                                       
*ENDIF                                                                  
                                                                        
        DO I=1,P_FIELD                                                  
          SEA(I)=.FALSE.                                                
          IF(FLAND_G(I).LT.1.0)SEA(I)=.TRUE.                            
        ENDDO                                                           
*D GRS2F404.52
            IF(SEA(I)) THEN                                             
*D RPANCA1A.874
c      write(6,*)'cmt rpanca1a.873 ',field,lookup_step(field),i2,i1,    
c     +lookup_start(file),nlookup(field),level,nlookups,                
c     +lookup(item_code,i1+LOOKUP_START(FILE)-1),                       
c     +lookup(item_code,i2+LOOKUP_START(FILE)-1)                        
cc                                                                      
cc      do i=1,45                                                       
cc       write(6,*)'cmtloop',i,lookup(i,i1+lookup_start(file)-1),       
cc     +lookup(i,i2+lookup_start(file)-1)                               
cc      enddo                                                           
cmt                                                                     
cmt  Current code checks that data doesn't go beyond end of file        
cmt  but doesn't check that data is the same type. Additional check     
cmt  made which if failed, forces code to go back to start of file.     
cmt  This avoids problems at year end.                                  
cmt                                                                     
        IF (I1.LE.FIXHD(152,FILE)                                       
     +  .AND. (LOOKUP(ITEM_CODE,I1+LOOKUP_START(FILE)-1) .EQ.           
     +         STASHANCIL(FIELD)) ) THEN                                
cmt                                                                     
*I GRS2F404.186   
            IF(.NOT.LTLEADS)THEN                                        
*I RPANCA1A.968   
            ELSE                                                        
            CALL T_INT (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,           
     &              TIME,P_FIELD)                                       
            ENDIF                                                       
*D RPANCA1A.1025
            IF(SEA(I).AND.FIELD.EQ.26) THEN                             
*D RPANCA1A.1040
                IF(TSTAR_LAND(I).GT.TM) TSTAR_LAND(I)=TM                
*D RPANCA1A.1053,RPANCA1A.1054 
                IF(TSTAR_LAND(I).GT.TM.AND.ANCIL_DATA(I).GT.0.0) THEN   
                  TSTAR_LAND(I)=TM                                      
*D RPANCA1A.1074
            IF (SEA(I)) THEN                                            
*D RPANCA1A.1081,RPANCA1A.1085 
          IF(.NOT.LTLEADS)THEN                                          
            DO I=1,P_FIELD                                              
              IF(ICE_FRACTION(I).GT.0.0) THEN                           
                TSTAR_SSI(I)=AMIN1(TSTAR_SSI(I),TFS)                    
              ENDIF                                                     
            END DO                                                      
          ENDIF                                                         
*D RPANCA1A.1101,RPANCA1A.1102 
            IF(L_CTILE.OR.ICE_FRACTION(I).EQ.0.0)THEN                   
              IF (SEA(I)) THEN                                          
                IF (L_SSTANOM) THEN                                     
*D RPANCA1A.1104
                TSTAR_SEA(I)=ANCIL_DATA(I)+TSTAR_ANOM(I)                
*D RPANCA1A.1106
                TSTAR_ANOM(I)=TSTAR_SEA(I)-ANCIL_DATA(I)                
*D RPANCA1A.1108,RPANCA1A.1109 
                ELSE                                                    
                  TSTAR_SEA(I)=ANCIL_DATA(I)                            
                END IF                                                  
                IF(ICE_FRACTION(I).EQ.0.0)TSTAR_SSI(I)=TSTAR_SEA(I)     
*D GJT1F304.121
            IF (SEA(I)) THEN                                            
*D TJ240293.40
              D1(D1_ANCILADD(FIELD)+I-1)=TSTAR_LAND(I)                  
                                                                        
*D RPANCA1A.1121
            IF(SEA(I)) THEN                                             
*D RPANCA1A.1141
            IF(SEA(I)) THEN                                             
*I RPANCA1A.1161  
      IF(L_CTILE)THEN                                                   
        DO I=1,P_FIELD                                                  
          IF(SEA(I).AND.ICE_FRACTION(I).GT.0.0)THEN                     
            IF(LTLEADS.OR.LAMIPII)THEN                                  
                                                                        
              TSTAR_SSI(I)=ICE_FRACTION(I)*TSTAR_SICE(I)                
     &          +(1.-ICE_FRACTION(I))*TSTAR_SEA(I)                      
                                                                        
            ELSE                                                        
                                                                        
              TSTAR_SEA(I)=TFS                                          
              TSTAR_SICE(I)=(TSTAR_SSI(I)                               
     &          -(1.-ICE_FRACTION(I))*TSTAR_SEA(I))/ICE_FRACTION(I)     
                                                                        
            ENDIF                                                       
          ENDIF                                                         
C                                                                       
          TSTAR(I)=FLAND_G(I)*TSTAR_LAND(I)                             
     &      +(1.-FLAND_G(I))*TSTAR_SSI(I)                               
        ENDDO                                                           
      ELSE                                                              
        DO I=1,P_FIELD                                                  
          IF(LAND(I))THEN                                               
            TSTAR(I)=TSTAR_LAND(I)                                      
          ELSE                                                          
            TSTAR(I)=TSTAR_SSI(I)                                       
          ENDIF                                                         
        ENDDO                                                           
      ENDIF                                                             
                                                                        
!     Set up surface temperatures:                                      
      IF(L_CTILE)THEN                                                   
        DO I=1,P_FIELD                                                  
          TSTAR_LAND_CTILE(I)=TSTAR_LAND(I)                             
          TSTAR_SEA_CTILE(I)=TSTAR_SEA(I)                               
          TSTAR_SICE_CTILE(I)=TSTAR_SICE(I)                             
        ENDDO                                                           
      ENDIF                                                             
!                                                                       
*DECLARE SCREEN7A
*D SCREEN7A.24,SCREEN7A.26   
*D SCREEN7A.36,SCREEN7A.37   
     & P_POINTS,P_FIELD,P1,LAND1,LAND_PTS,LAND_FIELD,NTILES,            
     & LAND_INDEX,TILE_INDEX,TILE_PTS,FLANDG,                        
*D SCREEN7A.39
     & TILE_FRAC,TL_1,TSTAR_SSI,TSTAR_TILE,                             
*D SCREEN7A.41
     & Q1P5M,Q1P5M_TILE,T1P5M,T1P5M_TILE                                
*D SCREEN7A.54
     &,NTILES               ! IN Number of tiles per land point.        
*D SCREEN7A.56
     &,TILE_INDEX(LAND_FIELD,NTILES)                                    
*D SCREEN7A.58
     &,TILE_PTS(NTILES)     ! IN Number of tile points.                 
*D SCREEN7A.61,SCREEN7A.62   
     & SQ1P5                ! IN STASH flag for 1.5-metre sp humidity.  
*D SCREEN7A.66
     & FLANDG(P_FIELD)      ! IN Fraction of gridbox which is land.     
     &,CHR1P5M(LAND_FIELD,NTILES)                                       
*D SCREEN7A.74
     &,RESFT(LAND_FIELD,NTILES)                                         
*D SCREEN7A.76
     &,TILE_FRAC(LAND_FIELD,NTILES)                                     
*D SCREEN7A.80,SCREEN7A.81   
     &,TSTAR_SSI(P_FIELD)   ! IN Sea/sea-ice mean sfc temperature (K).  
     &,TSTAR_TILE(LAND_FIELD,NTILES)                                    
*D SCREEN7A.85
     &,Z0H_TILE(LAND_FIELD,NTILES)                                      
*D SCREEN7A.89
     &,Z0M_TILE(LAND_FIELD,NTILES)                                      
*I SCREEN7A.95    
     &,Q1P5M_TILE(LAND_FIELD,NTILES)                                    
!                           ! OUT Q1P5M over land tiles.                
*I SCREEN7A.97    
     &,T1P5M_TILE(LAND_FIELD,NTILES)                                    
!                           ! OUT T1P5M over land tiles.                
*D SCREEN7A.105,SCREEN7A.106  
*D SCREEN7A.129,SCREEN7A.131  
          IF (FLANDG(I).LT.1.0 ) THEN                                 
            T1P5M(I) = (1.-FLANDG(I))*                              
     &        (TSTAR_SSI(I) - GRCP*Z1P5M +                            
     &        CHR1P5M_SICE(I) *  (TL_1(I) - TSTAR_SSI(I) +        
     &          GRCP*(Z1(I)+Z0M(I)-Z0H(I))))                      
*D SCREEN7A.135
        DO N=1,NTILES                                                   
          DO L=1,LAND_FIELD                                             
            T1P5M_TILE(L,N) = 0.                                        
          ENDDO                                                         
*D SCREEN7A.139,SCREEN7A.140  
            T1P5M_TILE(L,N) = TSTAR_TILE(L,N) - GRCP*Z1P5M +            
     &                        CHR1P5M(L,N)*( TL_1(I) - TSTAR_TILE(L,N) +
*D SCREEN7A.142
            T1P5M(I) = T1P5M(I)                                     
     &        + FLANDG(I)*TILE_FRAC(L,N)*T1P5M_TILE(L,N)              
*D SCREEN7A.153
        CALL QSAT(QS(P1),TSTAR_SSI(P1),PSTAR(P1),P_POINTS)              
*D SCREEN7A.156
          IF (FLANDG(I).LT.1.0 ) THEN                                 
*D SCREEN7A.158
            Q1P5M(I) = (1.-FLANDG(I))*                              
     &        (QW_1(I) + CER1P5M*( QW_1(I) - QS(I) ))             
*D SCREEN7A.167
        DO N=1,NTILES                                                   
          DO L=1,LAND_FIELD                                             
            Q1P5M_TILE(L,N) = 0.                                        
          ENDDO                                                         
*D SCREEN7A.174,SCREEN7A.175  
            Q1P5M(I) = Q1P5M(I)                                     
     &        + FLANDG(I)*TILE_FRAC(L,N)*Q1P5M_TILE(L,N)              
*DECLARE SEED
*D SEED.4,SEED.5    
     + FRAC_MIN                   ! Minimum areal fraction for PFTs.    
     +,FRAC_SEED                  ! "Seed" fraction for PFTs.
      PARAMETER(FRAC_MIN = 1.0E-6, FRAC_SEED = 0.01)                 
*DECLARE SETMODL1
*I GAV0F405.6     
!   4.6     22/01/99   Remove lines which overwrite H_VERS to 0 for
!                      section 19 for all internal models.
!                      Richard Betts
*D GSS3F401.1081,GSS3F401.1083 
*DECLARE SFEVAP7A
*D SFEVAP7A.2
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*D SFEVAP7A.50,SFEVAP7A.55   
     & P_POINTS,P_FIELD,P1,LAND1,LAND_PTS,LAND_FIELD,NTILES,            
     & LAND_INDEX,TILE_INDEX,TILE_PTS,NSHYD,LTIMER,FLAND,               
     & ASHTF_TILE,CANOPY,DTRDZ_1,FLAKE,FRACA,SNOW_TILE,RESFS,           
     & RESFT,RHOKH_1,TILE_FRAC,SMC,WT_EXT_TILE,TIMESTEP,                
     & FQW_1,FQW_TILE,FTL_1,FTL_TILE,TSTAR_TILE,                        
     & ECAN,ECAN_TILE,ELAKE_TILE,ESOIL,ESOIL_TILE,EI_TILE,EXT           
*D SFEVAP7A.68
     &,NTILES                ! IN Number of tiles per land point.       
*D SFEVAP7A.70
     &,TILE_INDEX(LAND_FIELD,NTILES)                                    
*D SFEVAP7A.72
     &,TILE_PTS(NTILES)      ! IN Number of tile points.                
*D SFEVAP7A.79,SFEVAP7A.84   
     & FLAND(LAND_FIELD)     ! IN Fraction of gridbox which is land.    
     &,ASHTF_TILE(LAND_FIELD,NTILES)                                    
!                            ! IN Coefficient to calculate surface      
!                            !    heat flux into soil.                  
     &,CANOPY(LAND_FIELD,NTILES)                                        
!                            ! IN Surface/canopy water on land          
!                            !    tiles (kg/m2).                        
*D SFEVAP7A.86
     &,FLAKE(LAND_FIELD,NTILES)                                         
!                            ! IN Lake fraction.                        
     &,FRACA(LAND_FIELD,NTILES)                                         
*D SFEVAP7A.89,SFEVAP7A.91   
!                            !    for land tiles.                       
     &,SNOW_TILE(LAND_FIELD,NTILES)                                     
!                            ! IN Lying snow amount on tiles (kg/m2).   
     &,RESFS(LAND_FIELD,NTILES)                                         
*D SFEVAP7A.94,SFEVAP7A.95   
!                            !    of land tiles.                        
     &,RESFT(LAND_FIELD,NTILES)                                         
*D SFEVAP7A.98
     &,RHOKH_1(LAND_FIELD,NTILES)                                       
*D SFEVAP7A.100
     &,TILE_FRAC(LAND_FIELD,NTILES)                                     
*D SFEVAP7A.103
     &,WT_EXT_TILE(LAND_FIELD,NSHYD,NTILES)                             
*D SFEVAP7A.105
!                            !    extracted from each soil layer        
!                            !    by each tile.                         
*D SFEVAP7A.110
     &,FQW_TILE(LAND_FIELD,NTILES)                                      
*D SFEVAP7A.113
     &,FTL_TILE(LAND_FIELD,NTILES)                                      
*D SFEVAP7A.115,SFEVAP7A.119  
     &,TSTAR_TILE(LAND_FIELD,NTILES)                                    
*D SFEVAP7A.126,SFEVAP7A.127  
     &,ECAN_TILE(LAND_FIELD,NTILES)                                     
!                            ! OUT ECAN for land tiles.                 
     &,ELAKE_TILE(LAND_FIELD,NTILES)                                    
!                            ! OUT Lake evaporation.                    
*D SFEVAP7A.131,SFEVAP7A.132  
     &,ESOIL_TILE(LAND_FIELD,NTILES)                                    
!                            ! OUT ESOIL for land tiles.                
     &,EI_TILE(LAND_FIELD,NTILES)                                       
!                            ! OUT Sublimation from snow or land-ice    
!                            !     (kg per sq m per s).                 
*I SFEVAP7A.136   
*CALL C_0_DG_C                                                          
*D SFEVAP7A.143,SFEVAP7A.145  
     &,E_TILE_OLD(LAND_FIELD,NTILES)                                    
!                            ! Surface moisture flux before adjustment. 
     &,LE_TILE_OLD(LAND_FIELD,NTILES)                                   
!                            ! Surf latent heat flux before adjustment. 
*D SFEVAP7A.152,SFEVAP7A.153  
                                                                        
*D SFEVAP7A.165
      DO N=1,NTILES                                                     
*D SFEVAP7A.168
          E_TILE_OLD(L,N) = FQW_TILE(L,N)                               
          IF (SNOW_TILE(L,N) .GT. 0.) THEN                              
            LE_TILE_OLD(L,N) = (LC + LF)*FQW_TILE(L,N)                  
          ELSE                                                          
            LE_TILE_OLD(L,N) = LC*FQW_TILE(L,N)                         
          ENDIF                                                         
*D SFEVAP7A.172
      DO N=1,NTILES                                                     
*I SFEVAP7A.175   
          ELAKE_TILE(L,N) = 0.                                          
          EI_TILE(L,N) = 0.                                             
*D SFEVAP7A.179,SFEVAP7A.182  
*D SFEVAP7A.184
! Sublimation from snow-covered land tiles                              
*D SFEVAP7A.186,SFEVAP7A.197  
      DO N=1,NTILES                                                     
        DO J=1,TILE_PTS(N)                                              
          L = TILE_INDEX(J,N)                                           
          IF (SNOW_TILE(L,N) .GT. 0.) THEN                              
            EI_TILE(L,N) =  FQW_TILE(L,N)                               
            EDT = EI_TILE(L,N)*TIMESTEP                                 
            IF ( EDT .GT. SNOW_TILE(L,N) )                              
     &        EI_TILE(L,N) = SNOW_TILE(L,N) / TIMESTEP                  
            FQW_TILE(L,N) = FQW_TILE(L,N) -  EI_TILE(L,N)               
          ENDIF                                                         
        ENDDO                                                           
*D SFEVAP7A.202
*D SFEVAP7A.209
      DO N=1,NTILES                                                     
*D SFEVAP7A.214,SFEVAP7A.215  
            ECAN_TILE(L,N) = (1. - FLAKE(L,N)) *                        
     &                       FRACA(L,N) * FQW_TILE(L,N) / RESFT(L,N)    
            ESOIL_TILE(L,N) = (1. - FLAKE(L,N)) *                       
     &                        (1. - FRACA(L,N))*RESFS(L,N)*FQW_TILE(L,N)
*I SFEVAP7A.216   
            ELAKE_TILE(L,N) = FLAKE(L,N)*FQW_TILE(L,N) / RESFT(L,N)     
*D SFEVAP7A.219
              ESOIL_TILE(L,N) =  (1. - FLAKE(L,N)) *                    
     &                           (1. - FRACA(L,N)*CANOPY(L,N)/EDT) *    
*D SFEVAP7A.223,SFEVAP7A.225  
          ELSEIF (SNOW_TILE(L,N).LE.0.) THEN                            
            IF (TSTAR_TILE(L,N).GE.TM) THEN                             
              ECAN_TILE(L,N) = (1. - FLAKE(L,N))*FQW_TILE(L,N)          
              ELAKE_TILE(L,N) = FLAKE(L,N)*FQW_TILE(L,N)                
            ELSE                                                        
              EI_TILE(L,N) =  FQW_TILE(L,N)                             
            ENDIF                                                       
*D SFEVAP7A.239
          DO N=1,NTILES                                                 
*D SFEVAP7A.249
          EXT(L,K) = 0.                                                 
        ENDDO                                                           
      ENDDO                                                             
                                                                        
      DO K=1,NSHYD                                                      
        DO N=1,NTILES                                                   
          DO J=1,TILE_PTS(N)                                            
            L = TILE_INDEX(J,N)                                         
            EXT(L,K) = EXT(L,K) + TILE_FRAC(L,N)*WT_EXT_TILE(L,K,N)     
     &                                          *ESOIL_TILE(L,N)        
          ENDDO                                                         
*D SFEVAP7A.262,SFEVAP7A.263  
      DO N=1,NTILES                                                     
*D ABX1F405.912,SFEVAP7A.266  
          DIFF_LAT_HTF = (LC + LF)*EI_TILE(L,N) + LC*ECAN_TILE(L,N)     
     &                    + LC*ESOIL_TILE(L,N) + LC*ELAKE_TILE(L,N)     
     &                    - LE_TILE_OLD(L,N)                            
*D SFEVAP7A.268
     &                        ( 1. + ASHTF_TILE(L,N)/(CP*RHOKH_1(L,N)) )
*D SFEVAP7A.270
          DTSTAR = - (DIFF_LAT_HTF + DIFF_SENS_HTF) / ASHTF_TILE(L,N)   
*D SFEVAP7A.273
          DFQW(L) = DFQW(L) + TILE_FRAC(L,N)*( ECAN_TILE(L,N) +         
     &                  ESOIL_TILE(L,N) + EI_TILE(L,N) + ELAKE_TILE(L,N)
     &                  - E_TILE_OLD(L,N) )                             
*D SFEVAP7A.275,SFEVAP7A.289  
*D SFEVAP7A.298,SFEVAP7A.301  
        FTL_1(I) = FTL_1(I) + FLAND(L)*DFTL(L)                      
        FQW_1(I) = FQW_1(I) + FLAND(L)*DFQW(L)                      
*DECLARE SFEXCH7A
*D SFEXCH7A.2
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*D SFEXCH7A.49,SFEXCH7A.54   
     & P_POINTS,P_FIELD,P1,LAND1,LAND_PTS,LAND_FIELD,NTILES,LAND_INDEX, 
     & TILE_INDEX,TILE_PTS,FLAND,FLANDG,                                
     & BQ_1,BT_1,CANHC_TILE,CANOPY,CATCH,DZSOIL,FLAKE,GC,HCONS,         
     & HO2R2_OROG,ICE_FRACT,SNOW_TILE,PSTAR,QW_1,RADNET,RADNET_TILE,    
     & SIL_OROG,SMVCST,TILE_FRAC,TIMESTEP,T_1,Q_1,QCF_1,QCL_1,          
     & TL_1,TI,TS1,
     & TSTAR_TILE,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,            
     & VFRAC_TILE,VSHR_LAND,VSHR_SSI,ZH,Z0_TILE,Z1_UV,Z1_TQ,LAND_MASK,  
*D SFEXCH7A.56
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,CD,CH,CDR10M,                
*D SFEXCH7A.59
     & RESFS,RESFT,RIB,RIB_TILE,                                        
     & FB_SURF,U_S,Q1_SD,T1_SD,TV1_SD,Z0M_EFF,                          
*D SFEXCH7A.62
     & RHO_CD_MODV1,RHOKH_1,RHOKH_1_SICE,RHOKM_1,RHOKM_LAND,RHOKM_SSI,
     & RHOKPM,RHOKPM_SICE,    
*D SFEXCH7A.76
     &,NTILES                ! IN Number of land tiles per land point.  
*D SFEXCH7A.78
     &,TILE_INDEX(LAND_FIELD,NTILES)                                    
*D SFEXCH7A.80
     &,TILE_PTS(NTILES)      ! IN Number of tile points.                
*D SFEXCH7A.87
     &,CANHC_TILE(LAND_FIELD,NTILES)                                    
!                            ! IN Areal heat capacity of canopy for     
!                            !    land tiles (J/K/m2).                  
     &,CANOPY(LAND_FIELD,NTILES)                                        
*D SFEXCH7A.90
     &,CATCH(LAND_FIELD,NTILES)                                         
*D SFEXCH7A.92
!                            !    of land tiles (kg/m2).                
*D SFEXCH7A.95
     &,FLAKE(LAND_FIELD,NTILES)                                         
!                            ! IN Lake fraction.                        
     &,GC(LAND_FIELD,NTILES) ! IN "Stomatal" conductance to evaporation 
*D SFEXCH7A.102
     &,FLAND(LAND_FIELD)     ! IN Land fraction on land tiles.          
     &,FLANDG(P_FIELD)       ! IN Land fraction on all tiles.           
     &,SNOW_TILE(LAND_FIELD,NTILES)                                     
!                            ! IN Lying snow on land tiles (kg/m2).     
*D SFEXCH7A.104
*D SFEXCH7A.107,SFEXCH7A.110  
     &,RADNET(P_FIELD)       ! IN Sea-ice net surface radiation (W/m2)  
     &,RADNET_TILE(LAND_FIELD,NTILES)                                   
!                            ! IN Land tile net surface radiation (W/m2)
*D SFEXCH7A.115
     &,TILE_FRAC(LAND_FIELD,NTILES)                                     
*I SFEXCH7A.117   
     &,T_1(P_FIELD)          ! IN Atmospheric temperature (K).          
     &,Q_1(P_FIELD)          ! IN Specific humidity ( kg/kg air).       
     &,QCF_1(P_FIELD)        ! IN Cloud ice (kg per kg air)             
     &,QCL_1(P_FIELD)        ! IN Cloud liquid water (kg                
!                            !    per kg air).                          
*D SFEXCH7A.124,SFEXCH7A.126  
     &,TSTAR_TILE(LAND_FIELD,NTILES)                                    
*D SFEXCH7A.128,SFEXCH7A.129  
     &,TSTAR_LAND(P_FIELD)   ! IN Land mean surface temperature (K).    
     &,TSTAR_SEA(P_FIELD)    ! IN Open sea surface temperature (K).     
     &,TSTAR_SICE(P_FIELD)   ! IN Sea-ice surface temperature (K).      
     &,TSTAR_SSI(P_FIELD)    ! IN Mean sea surface temperature (K).     
     &,VFRAC_TILE(LAND_FIELD,NTILES)                                    
!                            ! IN Fractional canopy coverage for        
!                            !    land tiles.                           
     &,VSHR_LAND(P_FIELD)    ! IN Magnitude of land sfc-to-lowest-level 
*D SFEXCH7A.131
     &,VSHR_SSI(P_FIELD)     ! IN Mag. of mean sea sfc-to-lowest-level  
!                            !    wind shear                            
     &,ZH(P_FIELD)           ! IN Height above surface of top of        
!                            !    boundary layer (metres).              
     &,Z0_TILE(LAND_FIELD,NTILES)                                       
*D SFEXCH7A.133
*D SFEXCH7A.159
     & ALPHA1(LAND_FIELD,NTILES)                                        
*D SFEXCH7A.166,SFEXCH7A.168  
!                            !     flux into sea-ice (W/m2/K)           
     &,ASHTF_TILE(LAND_FIELD,NTILES)                                    
!                            ! OUT Coefficient to calculate surface heat
!                            !     flux into land tiles (W/m2/K)        
*I SFEXCH7A.170   
     &,CD_SSI(P_FIELD)       ! OUT Bulk transfer coefficient for        
!                            !      momentum over sea mean.             
*D SFEXCH7A.172
     &,CH_SSI(P_FIELD)       ! OUT Bulk transfer coefficient for heat   
!                            !    and/or moisture over sea mean.        
*D SFEXCH7A.179
     &,CHR1P5M(LAND_FIELD,NTILES)                                       
*D SFEXCH7A.190
     &,FQW_TILE(LAND_FIELD,NTILES)                                      
*D SFEXCH7A.195
     &,FTL_TILE(LAND_FIELD,NTILES)                                      
*D SFEXCH7A.198
     &,FRACA(LAND_FIELD,NTILES)                                         
*D SFEXCH7A.201
!                            !     for land tiles.                      
*I SFEXCH7A.206   
     &,RESFS(LAND_FIELD,NTILES)                                         
!                            ! OUT Combined soil, stomatal and          
!                            !     aerodynamic resistance factor for    
!                            !     fraction 1-FRACA of land tiles       
     &,RESFT(LAND_FIELD,NTILES)                                         
!                            ! OUT Total resistance factor              
!                            !     FRACA+(1-FRACA)*RESFS for snow-free  
!                            !     tiles, 1 for snow and land-ice.      
     &,RIB(P_FIELD)          ! OUT Mean bulk Richardson number for      
!                            !     lowest layer                         
     &,RIB_TILE(LAND_FIELD,NTILES)                                      
!                            ! OUT RIB for land tiles.                  
     &,FB_SURF(P_FIELD)      ! OUT Surface flux buoyancy over           
!                            !     density (m^2/s^3)                    
     &,U_S(P_FIELD)          ! OUT Surface friction velocity (m/s)      
*D SFEXCH7A.210,SFEXCH7A.221  
*I SFEXCH7A.224   
     &,TV1_SD(P_FIELD)       ! OUT Standard deviation of turbulent      
!                            !     fluctuations of surface layer        
!                            !     virtual temperature (K).             
*D SFEXCH7A.229
     &,Z0H_TILE(LAND_FIELD,NTILES)                                      
*D SFEXCH7A.233
     &,Z0M_TILE(LAND_FIELD,NTILES)                                      
*D SFEXCH7A.238
     &,RHO_ARESIST_TILE(LAND_FIELD,NTILES)                              
*D SFEXCH7A.240
     &,ARESIST_TILE(LAND_FIELD,NTILES)                                  
*D SFEXCH7A.242
     &,RESIST_B_TILE(LAND_FIELD,NTILES)                                 
*D SFEXCH7A.249
     &,RHOKH_1(LAND_FIELD,NTILES)                                       
*D SFEXCH7A.257
     &,RHOKM_LAND(P_FIELD)   ! OUT For land momentum. NB: This is output
!                            !     on UV-grid, but with the first and   
!                            !      last rows set to "missing data".    
     &,RHOKM_SSI(P_FIELD)    ! OUT For mean sea mom. NB: This is output 
!                            !     on UV-grid, but with the first and   
!                            !     last rows set to "missing data".     
     &,RHOKPM(LAND_FIELD,NTILES)                                        
*I SFEXCH7A.271   
*CALL C_EPSLON                                                          
*D SFEXCH7A.284
*D SFEXCH7A.287
*D SFEXCH7A.296
      EXTERNAL SF_OROG,SF_OROG_GB,QSAT,SF_RESIST,TIMER,                 
     &         SFL_INT_LAND,SFL_INT_SEA,                                
*D SFEXCH7A.316
*I SFEXCH7A.320   
     &,Z0M_SEA(P_FIELD)            ! Open sea roughness length for      
!                                  ! momentum transport.                
     &,DB_SEA(P_FIELD)             ! Buoyancy difference for sea points 
     &,V_S_SEA(P_FIELD)            ! Surface layer scaling velocity     
!                                  ! for sea points (m/s).              
     &,RECIP_L_MO_SEA(P_FIELD)     ! Reciprocal of the Monin-Obukhov    
!                                  ! length for sea points (m^-1).      
*I SFEXCH7A.324   
     &,CD_LAND(P_FIELD)            ! Bulk transfer coefficient for      
!                                  !      momentum over land.           
*D SFEXCH7A.333
*I SFEXCH7A.335   
     &,DB_ICE(P_FIELD)             ! Buoyancy difference for sea ice    
     &,V_S_ICE(P_FIELD)            ! Surface layer scaling velocity     
!                                  ! for sea ice (m/s).                 
     &,V_S_MIZ(P_FIELD)            ! Surface layer scaling velocity     
!                                  ! for marginal sea ice (m/s).        
     &,RECIP_L_MO_ICE(P_FIELD)     ! Reciprocal of the Monin-Obukhov    
!                                  ! length for sea ice (m^-1).         
     &,RECIP_L_MO_MIZ(P_FIELD)     ! Reciprocal of the Monin-Obukhov    
!                                  ! length for marginal sea ice (m^-1).
     &,RHO_ARESIST_LAND(P_FIELD)   ! Land mean of rho_aresist_tile      
!
*D SFEXCH7A.342
     & CD_STD(LAND_FIELD,NTILES)   ! Local drag coefficient for calc    
*D SFEXCH7A.344,SFEXCH7A.345  
     &,CD_TILE(LAND_FIELD,NTILES)  ! Drag coefficient                   
     &,CH_TILE(LAND_FIELD,NTILES)  ! Transfer coefficient for heat and  
*I SFEXCH7A.350   
     &,FZ0(LAND_FIELD)             ! Aggregation function for Z0.       
*D SFEXCH7A.352,SFEXCH7A.353  
     &,QSTAR_TILE(LAND_FIELD,NTILES)!Surface saturated sp humidity.     
     &,RHOKM_1_TILE(LAND_FIELD,NTILES)                                  
*D SFEXCH7A.355
     &,WIND_PROFILE_FACTOR(LAND_FIELD,NTILES)                           
*D SFEXCH7A.359,SFEXCH7A.360  
     &,Z0_GB(LAND_FIELD)           ! GBM roughness length               
     &,Z0M_EFF_TILE(LAND_FIELD,NTILES)                                  
*D SFEXCH7A.362
     &,Z0F_TILE(LAND_FIELD,NTILES) !Roughness length for free convective
*I SFEXCH7A.363   
     &,DB_TILE(LAND_FIELD,NTILES)  ! Buoyancy difference for surface    
!                                  ! tile                               
     &,V_S_TILE(LAND_FIELD,NTILES) ! Surface layer scaling velocity     
!                                  ! for tiles (m/s).                   
     &,V_S_STD(LAND_FIELD,NTILES)  ! Surface layer scaling velocity     
!                                  ! for tiles excluding orographic     
!                                  ! form drag (m/s).                   
     &,RECIP_L_MO_TILE(LAND_FIELD,NTILES)                               
!                                  ! Reciprocal of the Monin-Obukhov    
!                                  ! length for tiles (m^-1).           
*D ABX1F405.902
      DO N=1,NTILES                                                     
*D SFEXCH7A.389
        IF ( ICE_FRACT(I).GT.0.0 .AND. FLANDG(I).LT.1.0 ) THEN         
*D SFEXCH7A.400,SFEXCH7A.402  
        RHOSTAR(I) = PSTAR(I) / 
     &    ( R*(FLANDG(I)*TSTAR_LAND(I) +                            
     &    (1.-FLANDG(I))*TSTAR_SSI(I)) )                            
*D SFEXCH7A.404,SFEXCH7A.408  
*D SFEXCH7A.412
      CALL QSAT(QSTAR_ICE(P1),TSTAR_SICE(P1),PSTAR(P1),P_POINTS)        
*D SFEXCH7A.417
      DO N=1,NTILES                                                     
*I SFEXCH7A.437   
        DB_SEA(I) = 0.                                                  
        DB_ICE(I) = 0.                                                  
*D SFEXCH7A.443
      DO N=1,NTILES                                                     
*D SFEXCH7A.446,SFEXCH7A.450  
          IF ( SNOW_TILE(L,N).GT.0. .AND. SMVCST(L).NE.0. ) THEN        
            Z0 = Z0_TILE(L,N) - 4.0E-4*SNOW_TILE(L,N)                   
            ZETA1 = MIN( 5.0E-4 , Z0_TILE(L,N)  )                       
*I SFEXCH7A.451   
          ELSE                                                          
            Z0M_TILE(L,N) = Z0_TILE(L,N)                                
*I SFEXCH7A.455   
          DB_TILE(L,N) = 0.                                             
*D SFEXCH7A.459
      DO N=1,NTILES                                                     
*D SFEXCH7A.470
! of Richardson number. RESFT=1 for snow and land-ice.                  
*D SFEXCH7A.473,SFEXCH7A.474  
      DO N=1,NTILES                                                     
*D SFEXCH7A.486,SFEXCH7A.487  
     &   CANOPY(1,N),CATCH(1,N),CHN,DQ,EPDT,FLAKE(1,N),GC(1,N),         
     &   SNOW_TILE(1,N),VSHR_LAND,FRACA(1,N),RESFS(1,N),RESFT(1,N),     
     &   LTIMER                                                         
*D SFEXCH7A.489,SFEXCH7A.494  
*D SFEXCH7A.503,SFEXCH7A.506  
     & P_POINTS,P_FIELD,P1,FLANDG,NSICE,SICE_INDEX,                  
     & BQ_1,BT_1,ICE_FRACT,QSTAR_ICE,QSTAR_SEA,QW_1,TL_1,TSTAR_SICE,    
     & TSTAR_SEA,VSHR_SSI,Z0_ICE,Z0H_SEA,Z0_ICE,Z0MSEA,Z1_TQ,Z1_UV,     
     & RIB_SEA,RIB_ICE,DB_SEA,DB_ICE,LTIMER                             
*D SFEXCH7A.510
      DO N=1,NTILES                                                     
*D SFEXCH7A.514,SFEXCH7A.515  
     &   TSTAR_TILE(1,N),VSHR_LAND,Z0H_TILE(1,N),Z0M_TILE(1,N),         
     &   Z1_TQ,Z1_UV,                                                   
     &   RIB_TILE(1,N),DB_TILE(1,N),LTIMER                              
*D SFEXCH7A.524
      DO N=1,NTILES                                                     
*D SFEXCH7A.538,SFEXCH7A.540  
      CALL FCDCH_SEA(P_POINTS,P_FIELD,P1,FLANDG,                     
     &               RIB_ICE,DB_ICE,VSHR_SSI,Z0_ICE,Z0_ICE,Z0_ICE,      
     &               ZH,Z1_UV,Z1_TQ,                                    
     &               CD_ICE,CH_ICE,V_S_ICE,                             
     &               RECIP_L_MO_ICE,LTIMER)                             
*D SFEXCH7A.543,SFEXCH7A.545  
      CALL FCDCH_SEA(P_POINTS,P_FIELD,P1,FLANDG,                     
     &               RIB_ICE,DB_ICE,VSHR_SSI,Z0_MIZ,Z0_MIZ,Z0_MIZ,      
     &               ZH,Z1_UV,Z1_TQ,                                    
     &               CD_MIZ,CH_MIZ,V_S_MIZ,                             
     &               RECIP_L_MO_MIZ,LTIMER)                             
*D SFEXCH7A.548,SFEXCH7A.550  
      CALL FCDCH_SEA(P_POINTS,P_FIELD,P1,FLANDG,                     
     &               RIB_SEA,DB_SEA,VSHR_SSI,Z0MSEA,Z0H_SEA,Z0F_SEA,    
     &               ZH,Z1_UV,Z1_TQ,                                    
     &               CD_SEA,CH_SEA,V_S_SEA,                             
     &               RECIP_L_MO_SEA,LTIMER)                             
*D SFEXCH7A.553
      DO N=1,NTILES                                                     
*D SFEXCH7A.556,SFEXCH7A.558  
     &   RIB_TILE(1,N),DB_TILE(1,N),VSHR_LAND,                          
     &   Z0M_EFF_TILE(1,N),Z0H_TILE(1,N),Z0F_TILE(1,N),ZH,              
     &   Z1_UV,Z1_TQ,WIND_PROFILE_FACTOR(1,N),                          
     &   CD_TILE(1,N),CH_TILE(1,N),CD_STD(1,N),                         
     &   V_S_TILE(1,N),V_S_STD(1,N),RECIP_L_MO_TILE(1,N),LTIMER         
*D SFEXCH7A.563,SFEXCH7A.564  
!!  4.1 Recalculate RESFT using "true" CH and EPDT for land tiles       
*D SFEXCH7A.567
      DO N=1,NTILES                                                     
*D SFEXCH7A.572
          EPDT(L) = - RHOSTAR(I)*CH_TILE(L,N)*VSHR_LAND(I)
     &      *DQ(L)*TIMESTEP    
*D SFEXCH7A.576,SFEXCH7A.578  
     &   CANOPY(1,N),CATCH(1,N),CH_TILE(1,N),DQ,EPDT,FLAKE(1,N),        
     &   GC(1,N),SNOW_TILE(1,N),VSHR_LAND,FRACA(1,N),                   
     &   RESFS(1,N),RESFT(1,N),                                         
     &   LTIMER)                                                        
*I SFEXCH7A.585   
        CD_LAND(I) = 0.                                               
        CD_SSI(I) = 0.                                                
        CH_SSI(I) = 0.                                                
*D SFEXCH7A.590
        IF ( FLANDG(I).LT.1.0 ) THEN                                   
*D SFEXCH7A.592,SFEXCH7A.595  
            CD_SSI(I) = ( ICE_FRACT(I)*CD_MIZ(I) +                
     &              (0.7-ICE_FRACT(I))*CD_SEA(I) ) / 0.7  ! P2430.5 
            CH_SSI(I) = ( ICE_FRACT(I)*CH_MIZ(I) +                
     &              (0.7-ICE_FRACT(I))*CH_SEA(I) ) / 0.7  ! P2430.4 
*D SFEXCH7A.597,SFEXCH7A.600  
            CD_SSI(I) = ( (1.0-ICE_FRACT(I))*CD_MIZ(I) +          
     &              (ICE_FRACT(I)-0.7)*CD_ICE(I) ) / 0.3  ! P2430.7 
            CH_SSI(I) = ( (1.0-ICE_FRACT(I))*CH_MIZ(I) +          
     &              (ICE_FRACT(I)-0.7)*CH_ICE(I) ) / 0.3  ! P2430.7 
*I SFEXCH7A.601   
          CD(I)=(1.-FLANDG(I))*CD_SSI(I)                          
          CH(I)=(1.-FLANDG(I))*CH_SSI(I)                          
*D SFEXCH7A.607,SFEXCH7A.608  
      DO N=1,NTILES                                                     
       DO J=1,TILE_PTS(N)                                              
*D SFEXCH7A.611,SFEXCH7A.612  
          CD_LAND(I) = CD_LAND(I) + TILE_FRAC(L,N)*CD_TILE(L,N)     
          CD(I) = CD(I) + FLANDG(I)* TILE_FRAC(L,N)*CD_TILE(L,N)  
          CH(I) = CH(I) + FLANDG(I)*TILE_FRAC(L,N)*CH_TILE(L,N) 
*I SFEXCH7A.625   
        RHO_ARESIST_LAND(I)=0.0                                      
*D SFEXCH7A.628
        RHOKM_LAND(I) = 0.                                            
        RHOKM_SSI(I) = 0.                                             
*D SFEXCH7A.631,SFEXCH7A.636  
        IF ( FLANDG(I).LT.1.0 ) THEN                                  
          RHOKM_SSI(I) = RHOSTAR(I)*CD_SSI(I)*VSHR_SSI(I)   ! P243.124  
          RHOKH_1_SICE(I) = RHOSTAR(I)*CH_SSI(I)*VSHR_SSI(I)    
!                                                           ! P243.125  
          RHO_ARESIST(I) = RHOSTAR(I)*CD_SSI(I)*VSHR_SSI(I)     
          ARESIST(I) =  1. / (CD_SSI(I) * VSHR_SSI(I))            
          RESIST_B(I)= (CD_SSI(I)/CH_SSI(I) - 1.0) * ARESIST(I) 
*D SFEXCH7A.642
      DO N=1,NTILES                                                     
*D SFEXCH7A.651,SFEXCH7A.655  
          RHOKM_1_TILE(L,N) = RHOSTAR(I)*CD_TILE(L,N)*VSHR_LAND(I)  
!                                                         ! P243.124    
          RHOKM_LAND(I) = RHOKM_LAND(I) +                           
     &           TILE_FRAC(L,N)*RHOKM_1_TILE(L,N)                       
          RHOKH_1(L,N) = RHOSTAR(I)*CH_TILE(L,N)*VSHR_LAND(I)       
!                                                         ! P243.125    
          RHO_ARESIST_TILE(L,N) = RHOSTAR(I) * CD_STD(L,N)            
     &                * VSHR_LAND(I)                                  
          ARESIST_TILE(L,N) = 1. / ( CD_STD(L,N) * VSHR_LAND(I) )     
*I SFEXCH7A.657   
          IF (RESIST_B_TILE(L,N) .LT. 0.) RESIST_B_TILE(L,N) = 0.       
          RHO_ARESIST_LAND(I) = RHO_ARESIST_LAND(I) +               
     &                     TILE_FRAC(L,N)*RHO_ARESIST_TILE(L,N)         
*I SFEXCH7A.660   
      DO L=LAND1,LAND1+LAND_PTS-1                                       
        I = LAND_INDEX(L)                                               
        RHO_ARESIST(I) = FLANDG(I)*RHO_ARESIST_LAND(I) +          
     &                      (1.0-FLANDG(I))*RHO_ARESIST(I)          
        ARESIST(I) = RHOSTAR(I) / RHO_ARESIST(I)                        
      ENDDO                                                             
                                                                        
*I SFEXCH7A.661   
        RHOKM_1(I)= FLANDG(I) * RHOKM_LAND(I) +                   
     &                     (1.0-FLANDG(I)) * RHOKM_SSI(I)           
*D SFEXCH7A.667,SFEXCH7A.668  
!!  moisture.                                                           
*D SFEXCH7A.674
*D SFEXCH7A.677
! Sea and sea-ice                                                       
      CALL SF_FLUX_SEA (                                                
     & P_POINTS,P_FIELD,P1,NSICE,SICE_INDEX,FLANDG,                  
     & ICE_FRACT,QS1,QSTAR_ICE,QSTAR_SEA,QW_1,RADNET,RHOKH_1_SICE,      
     & TI,TL_1,TSTAR_SICE,TSTAR_SEA,Z0_ICE,Z0_ICE,Z0H_SEA,Z0MSEA,Z1_TQ, 
     & ALPHA1_SICE,ASHTF,E_SEA,FQW_ICE,FQW_1,FTL_ICE,FTL_1,H_SEA,       
     & RHOKPM_SICE,LTIMER                                               
     & )                                                                
                                                                        
! Land tiles                                                            
      DO N=1,NTILES                                                     
*D SFEXCH7A.684,SFEXCH7A.694  
      DO N=1,NTILES                                                     
        CALL SF_FLUX_LAND (                                             
     &   P_FIELD,LAND_FIELD,TILE_PTS(N),FLAND,
     &   LAND_INDEX,TILE_INDEX(1,N),     
     &   CANHC_TILE(1,N),DZSOIL,HCONS,QS1,QSTAR_TILE(1,N),QW_1,         
     &   RADNET_TILE(1,N),RESFT(1,N),RHOKH_1(1,N),SMVCST,SNOW_TILE(1,N),
     &   TILE_FRAC(1,N),TIMESTEP,TL_1,TS1,TSTAR_TILE(1,N),              
     &   VFRAC_TILE(1,N),Z0H_TILE(1,N),Z0M_EFF_TILE(1,N),Z1_TQ,         
     &   FQW_1,FTL_1,                                                   
     &   ALPHA1(1,N),ASHTF_TILE(1,N),FQW_TILE(1,N),FTL_TILE(1,N),       
     &   RHOKPM(1,N),LTIMER                                             
     & )                                                                
*D SFEXCH7A.696,SFEXCH7A.727  
*D SFEXCH7A.742,SFEXCH7A.743  
     & P_POINTS,P_FIELD,P1,FLANDG,                                   
     & BQ_1,BT_1,FQW_1,FTL_1,ICE_FRACT,RHOKM_SSI,RHOSTAR,VSHR_SSI,      
*D SFEXCH7A.749
      DO N=1,NTILES                                                     
*D SFEXCH7A.751
     &   P_FIELD,LAND_FIELD,TILE_PTS(N),
     &   LAND_INDEX,TILE_INDEX(1,N),FLAND,     
*D SFEXCH7A.753
     &   RHOSTAR,VSHR_LAND,Z0M_TILE(1,N),Z1_TQ,TILE_FRAC(1,N),          
*I SFEXCH7A.756   
                                                                        
!-----------------------------------------------------------------------
!! Calculate scaling parameters required for new boundary layer scheme  
!-----------------------------------------------------------------------
                                                                        
      DO I=P1,P1+P_POINTS-1                                             
        U_S(I) = SQRT(FLANDG(I)*CD_LAND(I)*VSHR_LAND(I)         
     &    +(1.-FLANDG(I))*CD_SSI(I)*VSHR_SSI(I))                  
        FB_SURF(I) = G * ( BT_1(I)*FTL_1(I) +                           
     &                     BQ_1(I)*FQW_1(I) ) / RHOSTAR(I)              
        IF (FB_SURF(I) .GT. 0.0) THEN                                   
          TV1_SD(I) = T_1(I) *                                          
     &                ( 1.0 + C_VIRTUAL*Q_1(I) - QCL_1(I) - QCF_1(I) ) *
     &                ( BT_1(I)*T1_SD(I) + BQ_1(I)*Q1_SD(I) )           
          IF (TV1_SD(I) .LT. 0.0) THEN                                  
            TV1_SD(I)=0.0                                               
          ENDIF                                                         
        ELSE                                                            
            TV1_SD(I)=0.0                                               
        ENDIF                                                           
      ENDDO                                                             
                                                                        
*D SFEXCH7A.776,SFEXCH7A.777  
        IF (FLANDG(I).LT.1.0) THEN                                    
          TAU = RHOKM_SSI(I) * VSHR_SSI(I)             ! P243.130   
*D SFEXCH7A.779
     &      TAU = RHOSTAR(I) * CD_SEA(I) * VSHR_SSI(I) * VSHR_SSI(I)  
*D SFEXCH7A.801,SFEXCH7A.802  
        IF ( .NOT.LAND_MASK(I) )H_BLEND_OROG(I) = H_BLEND_MIN       
        IF (FLANDG(I).LT.1.0) THEN                                    
*D SFEXCH7A.817
      DO L = LAND1,LAND1+LAND_PTS-1                                     
        Z0_GB(L) = 0.                                                   
        FZ0(L) = 0.
      ENDDO
                                                                        
! Copy sea/sea-ice roughness lengths and Richardson numbers onto the    
! the land point tiles.                                                 
                                                                        
! Weight so that gridbox mean can be calculated:                        
      DO L=LAND1,LAND1+LAND_PTS-1          
        I = LAND_INDEX(L)                         
        IF(FLAND(L).LT.1.0)THEN                                         
          RIB(I) = (1.0-FLAND(L))*RIB(I)                            
          FZ0(L) = (1.0-FLAND(L)) / (LOG(LB/Z0M(I))**2)               
        ENDIF                                                           
      ENDDO                                                             
                                                                        
      DO N=1,NTILES                                                     
*D SFEXCH7A.821
          RIB(I) = RIB(I) + FLAND(L)*TILE_FRAC(L,N)*RIB_TILE(L,N)  
          FZ0(L) = FZ0(L) 
     &      + FLAND(L)*TILE_FRAC(L,N) / (LOG(LB/Z0M_TILE(L,N))**2)
*D SFEXCH7A.824,SFEXCH7A.825  
*D SFEXCH7A.827,SFEXCH7A.833  
        Z0_GB(L) = LB * EXP( - SQRT(1./FZ0(L)) )                        
*D SFEXCH7A.852
          IF (FLANDG(I).LT.1.0 .AND. ICE_FRACT(I).LE.0. )           
     &      Z0M_SEA(I) = Z0MSEA(I)  
          IF (FLANDG(I).LT.1.0 .AND. ICE_FRACT(I).GT.0. ) THEN
*I SFEXCH7A.856   
            Z0M_SEA(I) = Z0_ICE(I)                                  
            V_S_SEA(I) = V_S_ICE(I)                                     
            RECIP_L_MO_SEA(I) = RECIP_L_MO_ICE(I)                       
*D SFEXCH7A.861,SFEXCH7A.864  
     &   P_POINTS,P_FIELD,P1,FLANDG,                                 
     &   VSHR_SSI,CD_SEA,CH_SEA,RIB,Z0M_SEA,Z0H_SEA,Z0F_SEA,Z1_UV,      
     &   RECIP_L_MO_SEA,V_S_SEA,                                        
     &   SU10,SV10,ST1P5,SQ1P5,                                         
     &   CDR10M,CHR1P5M_SICE,LTIMER                                     
*D SFEXCH7A.868
        DO N=1,NTILES                                                   
*D SFEXCH7A.870,SFEXCH7A.875  
     &     P_FIELD,LAND_FIELD,TILE_PTS(N),
     &     TILE_INDEX(1,N),LAND_INDEX,FLANDG,   
     &     VSHR_LAND,CD_STD(1,N),CD_TILE(1,N),CH_TILE(1,N),
     &     RIB_TILE(1,N),    
     &     TILE_FRAC(1,N),WIND_PROFILE_FACTOR(1,N),                     
     &     Z0M_EFF_TILE(1,N),Z0M_TILE(1,N),Z0H_TILE(1,N),               
     &     Z0F_TILE(1,N),Z1_UV,RECIP_L_MO_TILE(1,N),                    
     &     V_S_TILE(1,N),V_S_STD(1,N),                                  
     &     SU10,SV10,ST1P5,SQ1P5,                                       
     &     CDR10M,CHR1P5M(1,N),LTIMER                                   
*DECLARE SFFLUX7A
*D SFFLUX7A.35,SFFLUX7A.37   
     & P_FIELD,LAND_FIELD,TILE_PTS,FLAND,LAND_INDEX,TILE_INDEX,         
     & CANHC,DZSOIL,HCONS,QS1,QSTAR,QW_1,RADNET,RESFT,RHOKH_1,SMVCST,   
     & SNOW,TILE_FRAC,TIMESTEP,TL_1,TS1,TSTAR,VFRAC,Z0H,Z0M_EFF,Z1_TQ,  
*D SFFLUX7A.39
     & ALPHA1,ASHTF,FQW_1,FTL_1,RHOKPM,LTIMER                           
*D SFFLUX7A.55,SFFLUX7A.57   
     & FLAND(LAND_FIELD)                                                
     &,CANHC(LAND_FIELD)   ! IN Areal heat capacity of canopy (J/K/m2). 
     &,DZSOIL              ! IN Soil or land-ice surface layer          
!                          !    thickness (m).                          
     &,HCONS(LAND_FIELD)   ! IN Soil thermal conductivity (W/m/K).      
*D SFFLUX7A.62
     &,RADNET(LAND_FIELD)  ! IN Net surface radiation (W/m2) positive   
*I SFFLUX7A.65    
     &,SMVCST(LAND_FIELD)  ! IN Volumetric saturation point             
!                          !    - zero at land-ice points.              
     &,SNOW(LAND_FIELD)    ! IN Lying snow amount (kg/m2).              
*I SFFLUX7A.67    
     &,TIMESTEP            ! IN Timestep (s).                           
*I SFFLUX7A.71    
     &,VFRAC(LAND_FIELD)   ! IN Fractional canopy coverage.             
*D SFFLUX7A.81
     & ASHTF(LAND_FIELD)   ! OUT Coefficient to calculate surface       
!                          !     heat flux into soil (W/m2/K).          
     &,ALPHA1(LAND_FIELD)  ! OUT Gradient of saturated specific humidity
*I SFFLUX7A.92    
*CALL CSIGMA                                                            
*CALL C_SOILH                                                           
*D SFFLUX7A.94
*I SFFLUX7A.109   
     &,DS_RATIO            ! 2 * snowdepth / depth of top soil layer.   
*I SFFLUX7A.110   
     &,LH                  ! Latent heat (J/K/kg).                      
*I SFFLUX7A.139   
        ASHTF(L) = 2.0 * HCONS(L) / DZSOIL                              
        IF (SNOW(L).GT.0.0 .AND. SMVCST(L).NE.0.) THEN                  
          DS_RATIO = 2.0 * SNOW(L) / (RHO_SNOW * DZSOIL)                
          IF (DS_RATIO.LE.1.0) THEN                                     
            ASHTF(L) =  ASHTF(L) /                                      
     &                         (1. + DS_RATIO*(HCONS(L)/SNOW_HCON - 1.))
          ELSE                                                          
            ASHTF(L) =  ASHTF(L)*SNOW_HCON / HCONS(L)                   
          ENDIF                                                         
        ENDIF                                                           
        ASHTF(L) = (1. - VFRAC(L))*ASHTF(L) + CANHC(L)/TIMESTEP +       
     &             4*(1. + VFRAC(L))*SBCON*TS1(L)**3                    
      ENDDO                                                             
                                                                        
      DO J=1,TILE_PTS                                                   
        L = TILE_INDEX(J)                                               
*D SFFLUX7A.142
        LH = LC                                                         
        IF (SNOW(L) .GT. 0.) LH = LS                                    
        RHOKPM(L) = RHOKH_1(L) / ( ASHTF(L)  +                          
*D SFFLUX7A.144
        RAD_REDUC = RADNET(L) - ASHTF(L) * ( TL_1(I) - TS1(L)           
*I SFFLUX7A.145   
     &                        + CANHC(L)*(TSTAR(L) - TS1(L)) / TIMESTEP 
*D SFFLUX7A.149
     &                                + (CP*RHOKH_1(L) + ASHTF(L))*DQ1 )
*D SFFLUX7A.152,SFFLUX7A.153  
        FTL_1_GB(I) = FTL_1_GB(I) + FLAND(L)*TILE_FRAC(L)*FTL_1(L)      
        FQW_1_GB(I) = FQW_1_GB(I) + FLAND(L)*TILE_FRAC(L)*FQW_1(L)      
*D SFFLUX7A.171,SFFLUX7A.175  
     & P_POINTS,P_FIELD,P1,NSICE,SICE_INDEX,FLANDG,                  
     & ICE_FRACT,QS1,QSTAR_ICE,QSTAR_SEA,QW_1,RADNET,RHOKH_1,TI,        
     & TL_1,TSTAR_SICE,TSTAR_SEA,Z0H_ICE,Z0M_ICE,Z0H_SEA,Z0M_SEA,Z1_TQ, 
     & ALPHA1,ASHTF,E_SEA,FQW_ICE,FQW_1,FTL_ICE,FTL_1,H_SEA,RHOKPM,     
     & LTIMER)                                                          
*D SFFLUX7A.188
*D SFFLUX7A.191,SFFLUX7A.192  
     & FLANDG(P_FIELD)                                          
*D SFFLUX7A.205
     &,TSTAR_SICE(P_FIELD)  ! IN Sea-ice surface temperature (K).       
*I SFFLUX7A.218   
     &,ASHTF(P_FIELD)      ! OUT Coefficient to calculate surface       
!                          !     heat flux into sea-ice (W/m2/K).       
*I SFFLUX7A.231   
*CALL C_KAPPAI                                                          
*I SFFLUX7A.233   
*CALL CSIGMA                                                            
*D SFFLUX7A.273
        D_T = TSTAR_SICE(I) - TL_1(I)                                   
*I SFFLUX7A.282   
        ASHTF(I) = 2 * KAPPAI / DE + 4*SBCON*TI(I)**3                   
*D SFFLUX7A.286
        IF ( FLANDG(I).LT.1.0 ) THEN                                  
*D SFFLUX7A.310,SFFLUX7A.311  
          FTL_1(I) = (1.-FLANDG(I))*(FTL_ICE(I) + H_SEA(I) / CP)        
          FQW_1(I) = (1.-FLANDG(I))*(FQW_ICE(I) + E_SEA(I))             
*DECLARE SFLINT6A
*D SFLINT6A.2
*IF DEF,A03_8A                                                          
*D SFLINT6A.21
!!!  SUBROUTINES SFL_INT_SEA AND SFL_INT_LAND--------------------------
*D SFLINT6A.45,SFLINT6A.50   
      SUBROUTINE SFL_INT_SEA (
     & P_POINTS,P_FIELD,P1,LAND_MASK
     &,VSHR,CD,CH,RIB,Z0M,Z0H,Z0F,Z1
     &,RECIP_L_MO,V_S
     +,SU10,SV10,ST1P5,SQ1P5
     &,CDR10M,CHR1P5M,LTIMER
*D ARN0F405.1829,ARN0F405.1830 
     &,VSHR(P_FIELD)     ! IN Wind speed difference between the
!                        !    surface and the lowest wind level in
!                        !    the atmosphere (m/s).
*D ARN0F405.1831,SFLINT6A.67   
*D SFLINT6A.72
                                                                        
      LOGICAL                                                           
     & LAND_MASK(P_FIELD)        ! IN T for land points, F otherwise.
     &,SU10                      ! IN 10m U-wind diagnostic flag        
     &,SV10                      ! IN 10m V-wind diagnostic flag        
     &,ST1P5                     ! IN screen temp diagnostic flag       
     &,SQ1P5                     ! IN screen specific humidity          
!                                !    diagnostic flag                   
     +,LTIMER                    ! IN TIMER diagnostics flag            
! Output variables                                                      
!                                                                       
      REAL                                                              
     + CDR10M(P_FIELD)   ! OUT interpolation coefficicent for 10m wind  
     +,CHR1P5M(P_FIELD)  ! OUT Interpolation coefficient for 1.5m       
!                        !     temperature                              

      REAL
     & RIB(P_FIELD)    ! DUMMY Used in 7A boundary layer scheme
     &,Z0F(P_FIELD)    ! DUMMY Used in 7A boundary layer scheme
     &,Z1(P_FIELD)     ! DUMMY Used in 7A boundary layer scheme

!*                                                                      
!*L---------------------------------------------------------------------
      EXTERNAL TIMER , PHI_M_H                                          
!*                                                                      
!*L---------------------------------------------------------------------
!    Local and other symbolic constants :-                              
*CALL C_VKMAN                                                           
      REAL Z_OBS_TQ,Z_OBS_WIND                                          
      PARAMETER (                                                       
     + Z_OBS_TQ = 1.5    ! Height of screen observations of temperature 
!                        ! and humidity.                                
     +,Z_OBS_WIND = 10.0 ! Height of surface wind observations.         
     +)                                                                 
!                                                                       
!  Define local storage.                                                
!                                                                       
!  (a) Local work arrays.                                               
!                                                                       
      REAL                                                              
     & Z_WIND(P_FIELD)     ! Height of wind observations.               
     &,Z_TEMP(P_FIELD)     ! Height of temperature and humidity         
!                          ! observations.                              
     &,PHI_M_OBS(P_FIELD)  ! Monin-Obukhov stability function for       
!                          ! momentum integrated to the wind observation
!                          ! height.                                    
     &,PHI_H_OBS(P_FIELD)  ! Monin-Obukhov stability function for       
!                          ! scalars integrated to their observation    
!                          ! height.                                    
!                                                                       
!  (b) Scalars.                                                         
!                                                                       
      INTEGER                                                           
     + I       ! Loop counter (horizontal field index).                 
!*                                                                      
      IF (LTIMER) THEN                                                  
        CALL TIMER('SFL_INT   ',3)                                      
      ENDIF                                                             
!                                                                       
!-----------------------------------------------------------------------
!! 1. If diagnostics required calculate M-O stability functions at      
!!    observation heights.                                              
!-----------------------------------------------------------------------
                                                                        
      IF (SU10 .OR. SV10 .OR. ST1P5 .OR. SQ1P5) THEN                    
        DO I=P1,P1+P_POINTS-1                                           
          Z_WIND(I) = Z_OBS_WIND                                        
          Z_TEMP(I) = Z_OBS_TQ + Z0H(I) - Z0M(I)
        ENDDO                                                           
        CALL PHI_M_H_SEA (P_POINTS,P_FIELD,P1,LAND_MASK,
     &                    RECIP_L_MO,Z_WIND,Z_TEMP,Z0M,Z0H,
     &                    PHI_M_OBS,PHI_H_OBS,LTIMER)
      ENDIF                                                             
                                                                        
!-----------------------------------------------------------------------
!! 2. If diagnostics required calculate interpolation coefficient       
!!    for 1.5m screen temperature and specific humidity.                
!-----------------------------------------------------------------------
!                                                                       
      IF (ST1P5 .OR. SQ1P5) THEN                                        
        DO I=P1,P1+P_POINTS-1                                           
          IF ( .NOT. LAND_MASK(I) ) THEN
            CHR1P5M(I) = CH(I) * VSHR(I) * PHI_H_OBS(I)/(VKMAN*V_S(I))
          ENDIF
        ENDDO                                                           
      ENDIF                                                             
!                                                                       
!-----------------------------------------------------------------------
!! 3. If diagnostics required calculate interpolation coefficient       
!!    for 10m winds.                                                    
!-----------------------------------------------------------------------
!                                                                       
      IF ( SU10 .OR. SV10 ) THEN
        DO I=P1,P1+P_POINTS-1                                           
          IF ( .NOT. LAND_MASK(I) ) THEN
            CDR10M(I) = CD(I) * VSHR(I) * PHI_M_OBS(I)/(VKMAN*V_S(I))
          ENDIF
        ENDDO                                                           
      ENDIF                                                             
!                                                                       
      IF (LTIMER) THEN                                                  
        CALL TIMER('SFL_INT ',4)                                        
      ENDIF                                                             
      RETURN                                                            
      END                                                               

!!!                                                                     
!!!---------------------------------------------------------------------
!*L  Arguments :-                                                       
      SUBROUTINE SFL_INT_LAND (
     & P_FIELD,LAND_FIELD,TILE_PTS,TILE_INDEX,LAND_INDEX
     &,VSHR,CD_STD,CD,CH,RIB,TILE_FRAC,WIND_PROFILE_FACTOR
     &,Z0M,Z0M_STD,Z0H,Z0F,Z1
     &,RECIP_L_MO,V_S,V_S_STD
     &,SU10,SV10,ST1P5,SQ1P5
     &,CDR10M,CHR1P5M,LTIMER
     +)                                                                 
      IMPLICIT NONE                                                     
                                                                        
      INTEGER                                                           
     & P_FIELD            ! IN Size of field on p-grid.
     &,LAND_FIELD         ! IN Number of land points.
     &,TILE_PTS           ! IN Number of tile points.
     &,TILE_INDEX(LAND_FIELD)
!                         ! IN Index of tile points.
     &,LAND_INDEX(P_FIELD)! IN Index of land points.
                                                                        
      REAL                                                              
     + Z0M(LAND_FIELD)      ! IN Roughness length for momentum (m).     
     +,Z0H(LAND_FIELD)      ! IN Roughness length for heat and          
!                        !    moisture (m).                             
     &,Z0M_STD(LAND_FIELD)  ! IN Roughness length for momentum without  
!                        !    orographic component (m).                 
     &,VSHR(P_FIELD)        ! IN Wind speed difference between the
!                           !    surface and the lowest wind level in
!                           !    the atmosphere (m/s).
     &,CD(LAND_FIELD)       ! IN Surface drag coefficient.              
     &,CH(LAND_FIELD)       ! IN Surface transfer coefficient for heat a
!                        !    moisture.                                 
     &,CD_STD(LAND_FIELD)   ! IN Surface drag coefficient excluding     
!                        !    orographic from drag.                     
     &,TILE_FRAC(LAND_FIELD)                                            
!                        ! IN Tile fraction.                            
     &,RECIP_L_MO(LAND_FIELD)                                           
!                        ! IN Reciprocal of the Monin-Obukhov length (m)
     &,V_S(LAND_FIELD)      ! IN Surface layer scaling velocity includin
!                        !    orographic form drag (m/s).
     &,V_S_STD(LAND_FIELD)  ! IN Surface layer scaling velocity excludin
*D SFLINT6A.86
     +,CHR1P5M(LAND_FIELD)  ! OUT Interpolation coefficient for 1.5m    
*D SFLINT6A.88,SFLINT6A.89   

      REAL
     & RIB(LAND_FIELD)   ! DUMMY Used in 7A boundary layer scheme
     &,WIND_PROFILE_FACTOR(LAND_FIELD)
!                        ! DUMMY Used in 7A boundary layer scheme
     &,Z0F(LAND_FIELD)   ! DUMMY Used in 7A boundary layer scheme
     &,Z1(P_FIELD)       ! DUMMY Used in 7A boundary layer scheme

*D SFLINT6A.112
     &,PHI_M_OBS(LAND_FIELD)  ! Monin-Obukhov stability function for    
*D SFLINT6A.115
     &,PHI_H_OBS(LAND_FIELD)  ! Monin-Obukhov stability function for    
*D SFLINT6A.118,SFLINT6A.119  
*D SFLINT6A.124,SFLINT6A.126  
     + I,J,L       ! Loop counter (horizontal field index).
*D SFLINT6A.138,SFLINT6A.139  
        DO J=1,TILE_PTS
          L = TILE_INDEX(J)
          I = LAND_INDEX(L)
*D SFLINT6A.141,SFLINT6A.142  
          Z_TEMP(I) = Z_OBS_TQ + Z0H(L) - Z0M(L)
*D SFLINT6A.144
        CALL PHI_M_H_LAND (P_FIELD,LAND_FIELD,TILE_PTS,
     &                     TILE_INDEX,LAND_INDEX,
*D SFLINT6A.155,ARN0F405.1838 
        DO J=1,TILE_PTS
          L = TILE_INDEX(J)
          I = LAND_INDEX(L)
          CHR1P5M(L) = CH(L) * VSHR(I) * PHI_H_OBS(L)/(VKMAN*V_S_STD(L))
*D SFLINT6A.166,ARN0F405.1841 
        DO J=1,TILE_PTS
          L = TILE_INDEX(J)
          I = LAND_INDEX(L)
          CDR10M(I) = CDR10M(I) + TILE_FRAC(L) *
     &                CD(L) * VSHR(I) * PHI_M_OBS(L)/(VKMAN*V_S(L))
*D ARN0F405.1844
        DO J=1,TILE_PTS
          L = TILE_INDEX(J)
          I = LAND_INDEX(L)
          Z_TEMP(I) = Z_OBS_TQ + Z0H(L) - Z0M_STD(L)
        ENDDO
        CALL PHI_M_H_LAND (P_FIELD,LAND_FIELD,TILE_PTS,
     &                     TILE_INDEX,LAND_INDEX,
*D ARN0F405.1847,ARN0F405.1848 
        DO J=1,TILE_PTS
          L = TILE_INDEX(J)
          I = LAND_INDEX(L)
          CDR10M(I) = CDR10M(I) + TILE_FRAC(L) *
     &                CD_STD(L) * VSHR(I) * PHI_M_OBS(L)/
     &                     (VKMAN*V_S_STD(L))
*DECLARE SFLINT7A
*D SFLINT7A.60,SFLINT7A.64   
     & P_FIELD,LAND_FIELD,TILE_PTS,TILE_INDEX,LAND_INDEX,FLANDG,        
     & VSHR,CD_STD,CD,CH,RIB,TILE_FRAC,WIND_PROFILE_FACTOR,             
     & Z0M_EFF,Z0M,Z0H,Z0F,Z1,                                          
     & RECIP_L_MO,V_S,V_S_STD,                                          
     & SU10,SV10,ST1P5,SQ1P5,                                           
     & CDR10M,CHR1P5M,LTIMER                                            
*D SFLINT7A.79
     & FLANDG(P_FIELD)   ! IN Land fraction                             
     &,CD_STD(LAND_FIELD)! IN Surface drag coefficient for shear stress 
*I SFLINT7A.115   
      REAL                                                              
     & VSHR(P_FIELD)          ! DUMMY Used in 6A boundary layer scheme  
     &,RECIP_L_MO(LAND_FIELD) ! DUMMY Used in 6A boundary layer scheme  
     &,V_S(LAND_FIELD)        ! DUMMY Used in 6A boundary layer scheme  
     &,V_S_STD(LAND_FIELD)    ! DUMMY Used in 6A boundary layer scheme  
                                                                        
*D SFLINT7A.198
     &                FLANDG(I)*TILE_FRAC(L)*CD10*WIND_PROFILE_FACTOR(L)
*D SFLINT7A.254,SFLINT7A.256  
     & P_POINTS,P_FIELD,P1,FLANDG,                                   
     & VSHR,CD,CH,RIB,Z0M,Z0H,Z0F,Z1,                                   
     & RECIP_L_MO,V_S,                                                  
     & SU10,SV10,ST1P5,SQ1P5,                                           
     & CDR10M,CHR1P5M,LTIMER                                            
*D SFLINT7A.267
     & FLANDG(P_FIELD)    ! IN Land fraction.                         
     &,CD(P_FIELD)        ! IN Effective surface drag coefficient,      
*D SFLINT7A.281,SFLINT7A.282  
     & SU10               ! IN 10m U-wind diagnostic flag               
*I SFLINT7A.293   
      REAL                                                              
     & VSHR(P_FIELD)       ! DUMMY Used in 6A boundary layer scheme     
     &,RECIP_L_MO(P_FIELD) ! DUMMY Used in 6A boundary layer scheme     
     &,V_S(P_FIELD)        ! DUMMY Used in 6A boundary layer scheme     
                                                                        
*D SFLINT7A.347
          IF ( FLANDG(I).LT.1.0 ) THEN                                
*I SFLINT7A.371   
            CDR10M(I) = CDR10M(I)*(1.0-FLANDG(I))
                                                                       
*D SFLINT7A.383
          IF ( FLANDG(I).LT.1.0 ) THEN                                
*DECLARE SFMELT7A
*D SFMELT7A.2
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*I SFMELT7A.25    
C                                                                       
C           Modified for MOSES II. RE 4/4/00                            
C                                                                       
*D SFMELT7A.28,SFMELT7A.34   
     & POINTS,P_FIELD,P1,LAND_FIELD,NTILES,LAND_INDEX                   
     &,TILE_INDEX,TILE_PTS,LAND_MASK,LTIMER,SIMLT,SMLT,FLANDG           
     &,ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,DTRDZ_1,ICE_FRACT            
     &,RHOKH_1,RHOKH_1_SICE,TILE_FRAC,TIMESTEP                          
     &,EI_TILE,FQW_1,FQW_ICE,FTL_1,FTL_ICE,FTL_TILE                     
     &,TSTAR_SEA,TSTAR_SSI,TSTAR_TILE,SNOW_TILE                         
     &,EI_LAND,EI_SICE                                                  
     &,SICE_MLT_HTF,SNOMLT_SURF_HTF,SNOWMELT,MELT_TILE                  
*D SFMELT7A.45,SFMELT7A.47   
     &,LAND_INDEX(P_FIELD)  ! IN Index of land points.                  
     &,NTILES               ! IN Number of tiles per land point.        
     &,TILE_INDEX(LAND_FIELD,NTILES)                                    
!                           ! IN Index of tile points.                  
     &,TILE_PTS(NTILES)     ! IN Number of tile points.                 
*D SFMELT7A.56
     & FLANDG(P_FIELD)      ! IN Fraction of gridbox which is land.     
     &,ALPHA1(LAND_FIELD,NTILES)                                        
!                           ! IN Gradients of saturated specific        
*D SFMELT7A.59
!                           !    and land tile surfaces.                
*D SFMELT7A.63,SFMELT7A.64   
     &,ASHTF_TILE(LAND_FIELD,NTILES)                                    
!                           ! IN Coefficient to calculate surface       
!                           !    heat flux into soil.                   
*D SFMELT7A.68,SFMELT7A.69   
     &,RHOKH_1(LAND_FIELD,NTILES)                                       
!                           ! IN Surface exchange coeffs for land tiles.
*D SFMELT7A.72,SFMELT7A.73   
     &,TILE_FRAC(LAND_FIELD,NTILES)                                     
!                           ! IN Tile fractions.                        
*D SFMELT7A.77
     & EI_TILE(LAND_FIELD,NTILES)                                       
!                           ! INOUT Sublimation for land tiles (kg/m2/s)
     &,FQW_1(P_FIELD)       ! INOUT GBM surface moisture flux (kg/m2/s).
*D SFMELT7A.79
*D SFMELT7A.81,SFMELT7A.88   
     &,FTL_ICE(P_FIELD)     ! INOUT FTL for sea-ice.                    
     &,FTL_TILE(LAND_FIELD,NTILES)                                      
!                           ! INOUT FTL for land tiles.                 
     &,TSTAR_SEA(P_FIELD)   ! IN Open sea surface temperature (K).      
     &,TSTAR_SSI(P_FIELD)   ! INOUT Sea mean surface temperature (K).   
     &,TSTAR_TILE(LAND_FIELD,NTILES)                                    
!                           ! INOUT Land tile surface temperatures (K). 
     &,SNOW_TILE(LAND_FIELD,NTILES)                                     
!                           ! INOUT Lying snow on land tiles (kg/m2).   
*D SFMELT7A.91,SFMELT7A.92   
     & EI_LAND(P_FIELD)     ! OUT Sublimation from lying snow (kg/m2/s).
     &,EI_SICE(P_FIELD)     ! OUT Sublimation from sea-ice (kg/m2/s).   
     &,MELT_TILE(LAND_FIELD,NTILES)                                     
!                           ! OUT Surface snowmelt on tiles (kg/m2/s).  
*D SFMELT7A.98
     &,SNOWMELT(P_FIELD)    ! OUT GBM surface snowmelt (kg/m2/s).       
*D SFMELT7A.101
*I SFMELT7A.120   
     &,N                    ! Loop counter - tile index.                
cccccc Tibet Snow mod cccccc
      REAL MASKD
      PARAMETER( MASKD = 0.1 )
cccccccccccccccccccccccccccc
*D SFMELT7A.130
        EI_LAND(I) = 0.0                                              
        EI_SICE(I) = 0.0                                              
*D SFMELT7A.133,SFMELT7A.135  
      DO N=1,NTILES                                                     
        DO L=1,LAND_FIELD                                               
          MELT_TILE(L,N) = 0.                                           
        ENDDO                                                           
      ENDDO                                                             
                                                                        
*D SFMELT7A.137
!  Melt snow on land tiles if TSTAR_TILE is greater than TM.            
*D SFMELT7A.139,SFMELT7A.153  
      DO N=1,NTILES                                                     
        DO J=1,TILE_PTS(N)                                              
          L = TILE_INDEX(J,N)                                           
          I = LAND_INDEX(L)                                             
          SNOW_MAX = MAX( 0.0, SNOW_TILE(L,N) - EI_TILE(L,N)*TIMESTEP ) 
          IF ( SNOW_MAX.GT.0.0 .AND. TSTAR_TILE(L,N).GT.TM ) THEN       
            LCMELT = (CP + LC*ALPHA1(L,N))*RHOKH_1(L,N)                 
     &               + ASHTF_TILE(L,N)                                  
            LSMELT = LCMELT + LF*ALPHA1(L,N)*RHOKH_1(L,N)               
cccccc Tibet Snow mod cccccc
c           DTSTAR = - MIN( TSTAR_TILE(L,N) - TM ,                      
c    &                      LF*SNOW_MAX / (LCMELT*TIMESTEP) )           
            DTSTAR = - MIN( (TSTAR_TILE(L,N) - TM) *
     &                      (1.0 - EXP(-MASKD*SNOW_MAX)) ,
     &                      LF*SNOW_MAX / (LCMELT*TIMESTEP) )           
cccccccccccccccccccccccccccc
            MELT_TILE(L,N) = - LSMELT*DTSTAR / LF                       
            DFTL = CP*RHOKH_1(L,N)*DTSTAR                               
            DFQW = ALPHA1(L,N)*RHOKH_1(L,N)*DTSTAR                      
            FTL_TILE(L,N) = FTL_TILE(L,N) + DFTL                        
            EI_TILE(L,N) = EI_TILE(L,N) + DFQW                          
            TSTAR_TILE(L,N) = TSTAR_TILE(L,N) + DTSTAR                  
*D SFMELT7A.157,SFMELT7A.166  
            DFTL = TILE_FRAC(L,N)*DFTL                                  
            DFQW = TILE_FRAC(L,N)*DFQW                                  
            FTL_1(I) = FTL_1(I) + FLANDG(I)*DFTL                  
            FQW_1(I) = FQW_1(I) + FLANDG(I)*DFQW                  
          ENDIF                                                         
          EI_LAND(I) = EI_LAND(I) + TILE_FRAC(L,N)*EI_TILE(L,N)     
        ENDDO                                                           
*I SFMELT7A.168   
!-----------------------------------------------------------------------
!  Increment snow by sublimation and melt                               
!-----------------------------------------------------------------------
      DO N=1,NTILES                                                     
        DO J=1,TILE_PTS(N)                                              
          L = TILE_INDEX(J,N)                                           
          I = LAND_INDEX(L)                                             
          SNOW_TILE(L,N) = SNOW_TILE(L,N) -                             
     &                     (EI_TILE(L,N) + MELT_TILE(L,N))*TIMESTEP     
          SNOWMELT(I) = SNOWMELT(I) + TILE_FRAC(L,N)*MELT_TILE(L,N)     
        ENDDO                                                           
      ENDDO                                                             
      IF (SMLT) THEN                                                    
        DO I=P1,P1+POINTS-1                                             
          SNOMLT_SURF_HTF(I) = LF*SNOWMELT(I)                           
        ENDDO                                                           
      ENDIF                                                             
                                                                        
*D SFMELT7A.170
        IF ( FLANDG(I).LT.1.0 .AND. ICE_FRACT(I).GT.0.0 ) THEN      
*D SFMELT7A.172
!   Melt sea-ice if TSTAR > TSTARMAX                                    
*D SFMELT7A.174,SFMELT7A.176  
          EI_SICE(I) = FQW_ICE(I)                                   
          TSTARMAX = ICE_FRACT(I)*TM                                  
     &         + (1.0 - ICE_FRACT(I))*TSTAR_SEA(I)                  
          IF ( TSTAR_SSI(I) .GT. TSTARMAX ) THEN                      
*D SFMELT7A.178,SFMELT7A.179  
     &                       + ICE_FRACT(I)*GAMMA(1)*DTRDZ_1(I) )   
            DTSTAR = TSTARMAX - TSTAR_SSI(I)                          
*D SFMELT7A.181
     &                                                 + ASHTF(I)     
*D SFMELT7A.184,SFMELT7A.185  
            TSTAR_SSI(I) = TSTARMAX                                   
*D SFMELT7A.187,SFMELT7A.196  
            FTL_1(I) = FTL_1(I) + (1.0-FLANDG(I))*DFTL            
            FQW_1(I) = FQW_1(I) + (1.0-FLANDG(I))*DFQW            
            EI_SICE(I) = EI_SICE(I) + DFQW                          
            FTL_ICE(I) = FTL_ICE(I) + DFTL                          
            FQW_ICE(I) = FQW_ICE(I) + DFQW                          
*DECLARE SFREST7A
*D SFREST7A.30
     & CANOPY,CATCH,CH,DQ,EPDT,FLAKE,GC,SNOW,VSHR,                      
*I SFREST7A.56    
     &,FLAKE(LAND_FIELD)   ! IN Lake fraction.                          
*I SFREST7A.58    
     &,SNOW(LAND_FIELD)    ! IN Lying snow amount (kg per sq metre).    
*I SFREST7A.75    
cccccc Tibet Snow mod cccccccc
      REAL MASKD
      PARAMETER( MASKD = 0.1 )
cccccccccccccccccccccccccccccc
                                                                        
*D SFREST7A.93,SFREST7A.94   
! Set to 1 for negative moisture flux or snow-covered land              
! (no surface/stomatal resistance to condensation).                     
*D SFREST7A.97,SFREST7A.98   
cccccc Tibet Snow mod cccccccc
        IF(SNOW(L) .GT. 0.) FRACA(L) = 1.0 - EXP(-MASKD*SNOW(L))
cccccccccccccccccccccccccccccc
        IF (DQ(L).LT.0. .AND. SNOW(L).LE.0.) FRACA(L) = 0.0             
        IF (DQ(L).LT.0. .AND. SNOW(L).LE.0. .AND. CATCH(L).GT.0.)       
*D SFREST7A.106
        RESFT(L) = FLAKE(L) + (1. - FLAKE(L)) *                         
     &                        ( FRACA(L) + (1. - FRACA(L))*RESFS(L) )   
*DECLARE SFRIB7A
*D SFRIB7A.2
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*D SFRIB7A.37
     & RIB,DB,LTIMER                                                    
*I SFRIB7A.72    
     &,DB(LAND_FIELD)      ! OUT Buoyancy difference between surface    
!                          !     and lowest atmospheric level.          
*D SFRIB7A.113,SFRIB7A.114  
        DB(L) = G*(BT_1(I)*DTEMP(L) + BQ_1(I)*RESFT(L)*DQ(L))           
        RIB(L) = Z1_UV(I)*DB(L) / ( VSHR(I)*VSHR(I) )                   
*D SFRIB7A.130,SFRIB7A.131  
     & P_POINTS,P_FIELD,P1,FLANDG,NSICE,SICE_INDEX,                  
     & BQ_1,BT_1,ICE_FRACT,QSTAR_ICE,QSTAR_SEA,QW_1,TL_1,TSTAR_SICE,    
*D SFRIB7A.133
     & RIB_SEA,RIB_ICE,DB_SEA,DB_ICE,LTIMER                             
*D SFRIB7A.147
*D SFRIB7A.150
     & FLANDG(P_FIELD)     ! IN Land fraction on all pts.               
     &,BQ_1(P_FIELD)       ! IN A buoyancy parameter for lowest atm     
*D SFRIB7A.163
     &,TSTAR_SICE(P_FIELD)  ! IN Surface temperature of sea-ice (K).    
*I SFRIB7A.183   
     &,DB_SEA(P_FIELD)     ! OUT Buoyancy difference between surface    
!                          !     and lowest atmospheric level over      
!                          !     sea or sea-ice leads.                  
     &,DB_ICE(P_FIELD)     ! OUT Buoyancy difference between surface    
!                          !     and lowest atmospheric level over      
!                          !     sea-ice.                               
*D SFRIB7A.207
        IF ( FLANDG(I).LT.1.0 ) THEN                                   
*D SFRIB7A.212,SFRIB7A.213  
          DB_SEA(I) = G*( BT_1(I)*DTEMP + BQ_1(I)*DQ )                  
          RIB_SEA(I) = Z1_UV(I)*DB_SEA(I) / ( VSHR(I)*VSHR(I) )         
*D SFRIB7A.220
        DTEMP = TL_1(I) - TSTAR_SICE(I)                                 
*D SFRIB7A.223,SFRIB7A.224  
        DB_ICE(I) = G*( BT_1(I)*DTEMP + BQ_1(I)*DQ )                    
        RIB_ICE(I) = Z1_UV(I)*DB_ICE(I) / ( VSHR(I) * VSHR(I) )         
*DECLARE SFSNOW7A
*D SFSNOW7A.23,SFSNOW7A.28   
CLL  Purpose:  adds the large-scale and convective snowfall to the
CLL            snowdepth;  
*D SFSNOW7A.45,SFSNOW7A.49   
     & NPNTS,NTILES,TILE_PTS,TILE_INDEX,
     & CONV_SNOW,LS_SNOW,TILE_FRAC,TSTAR_TILE,TIMESTEP,      
     & RGRAIN,SNOW_TILE,L_SNOW_ALBEDO,
     & LTIMER)            
*D SFSNOW7A.55,SFSNOW7A.56   
     &,NTILES               ! IN Number of tiles.
     &,TILE_PTS(NTILES)     ! IN Number of tile points.
     &,TILE_INDEX(NPNTS,NTILES)
!                           ! IN Index of tile points.  
*D SFSNOW7A.61,SFSNOW7A.67   
     &,TILE_FRAC(NPNTS,NTILES)
                            ! IN Tile fractions.                        
     &,TSTAR_TILE(NPNTS,NTILES)
!                           ! IN Tile surface temperatures (K).
*D SFSNOW7A.69
*D SFSNOW7A.72,SFSNOW7A.80   
     & RGRAIN(NPNTS,NTILES) ! INOUT Snow grain size (microns).
     &,SNOW_TILE(NPNTS,NTILES)
!                           ! INOUT Snow on the ground (kg/m2). 
*D SFSNOW7A.83,SFSNOW7A.84   
     & L_SNOW_ALBEDO        ! IN Flag for prognostic snow albedo.       
*D SFSNOW7A.87,SFSNOW7A.89   
*CALL C_0_DG_C
*D SFSNOW7A.93,SFSNOW7A.95   
     & R0                  ! Grain size for fresh snow (microns).       
*D SFSNOW7A.98,SFSNOW7A.100  
     &,SNOWFALL(NPNTS)     ! Snowfall in timestep (kg/m2).
*D SFSNOW7A.102
      INTEGER I,J,N        ! Loop counters.                             
*D SFSNOW7A.109,SFSNOW7A.154  
*D SFSNOW7A.161,SFSNOW7A.162  
        SNOWFALL(I) = TIMESTEP*(LS_SNOW(I) + CONV_SNOW(I))
      ENDDO

      DO N=1,NTILES
        DO J=1,TILE_PTS(N)
          I = TILE_INDEX(J,N) 
          SNOW_TILE(I,N) = SNOW_TILE(I,N) + SNOWFALL(I)
        ENDDO
*D SFSNOW7A.169,SFSNOW7A.171  
      DO N=1,NTILES
        DO J=1,TILE_PTS(N)
          I = TILE_INDEX(J,N) 
          IF ( SNOW_TILE(I,N) .GT. 0.) THEN                             
*D SFSNOW7A.173,SFSNOW7A.174  
            IF (TSTAR_TILE(I,N) .LT. TM) THEN                           
              IF (RGRAIN(I,N) .LT. 150.) THEN                           
*D SFSNOW7A.177
                RATE = 0.23E6*EXP(-3.7E4/(8.13451*TSTAR_TILE(I,N)))     
*D SFSNOW7A.180,SFSNOW7A.183  
            RGRAIN(I,N) = SQRT( RGRAIN(I,N)**2 
     &                          + (RATE/3.14159)*TIMESTEP )  
     &                              - (RGRAIN(I,N) - R0)*SNOWFALL(I)/2.5
            RGRAIN(I,N) = MIN( RMAX, RGRAIN(I,N) )                      
            RGRAIN(I,N) = MAX( R0, RGRAIN(I,N) )                        
*D SFSNOW7A.185
            RGRAIN(I,N) = R0                                            
*I SFSNOW7A.187   
      ENDDO                                                           
*DECLARE SICEHT5B
*I SICEHT5B.44    
!!!
!!! Note: 
!!!   Routine renamed to avoid clash with SICEHT7A which replaces 
!!!   this deck in MOSES2.2 7A boundary layer. 
!!!
*D SICEHT5B.48
      SUBROUTINE SICE_HTF_5B (
*DECLARE SMCEXT7A
*D SMCEXT7A.24
     &,                   F_ROOT,STHU,V_CRIT,V_SAT,V_WILT
*D SMCEXT7A.67
*D SMCEXT7A.137,SMCEXT7A.138  
     &      WT_EXT(I,N) = F_ROOT(N)*FSMC_L(I,N)/FSMC(I)
  
*DECLARE SOILHT7A
*D SOILHT7A.24,SOILHT7A.26   
     & NPNTS,NSHYD,NTILES,SOIL_PTS,SOIL_INDEX,TILE_PTS,TILE_INDEX       
     &,BEXP,DZ,FRAC,HCAP,HCON,SATHH,SNOW_TILE,SURF_HT_FLUX,TIMESTEP
     &,V_SAT,W_FLUX,SMCL,STHU,STHF,TSOIL                        
*I SOILHT7A.37    
!     GAUSS - to solve matrix equation for temperature increments       
*I GPB8F405.40    
!  4.6     10/99     Implicit updating. P. Cox, R. Essery  
*I SOILHT7A.63    
     &,NTILES               ! IN Number of tiles.       
*I SOILHT7A.72    
     &,TILE_PTS(NTILES)     ! IN Number of tile points.
     &,TILE_INDEX(NPNTS,NTILES)
!                           ! IN Index of tile points.                  
*D SOILHT7A.75
     & BEXP(NPNTS)          ! IN Clapp-Hornberger exponent.             
*I SOILHT7A.76    
     &,FRAC(NPNTS,NTILES)   ! IN Tile fractions.    
*D SOILHT7A.83,SOILHT7A.84   
     &,SNOW_TILE(NPNTS,NTILES)
!                           ! IN Lying snow on tiles (kg/m2). 
     &,SURF_HT_FLUX(NPNTS)  ! IN Net downward surface heat flux (W/m2). 
*D SOILHT7A.108
     &,ITER_PTS             ! WORK Number of soil points which require  
*I SOILHT7A.116   
     &,SI_TILE              ! WORK Tile snow insulation factor.
     &,SNOW_DEPTH           ! WORK Depth of lying snow (m). 
*D SOILHT7A.123
     & ITER_INDEX(NPNTS)    ! WORK Array of soil points which require   
*D SOILHT7A.133,SOILHT7A.134  
*D SOILHT7A.142
     &,DTSLMAX(NPNTS,NSHYD) ! WORK Maximum value of DTSL (K/timestep).  
     &,DTSLMIN(NPNTS,NSHYD) ! WORK Minimum value of DTSL (K/timestep).
     &,HCAPT(NPNTS,NSHYD)  ! WORK The total volumetric heat capacity   
*I SOILHT7A.149   
     &,HADV(NPNTS,NSHYD)    ! WORK Heat flux due to moisture advection
!                           !      (W/m2).
     &,SIFACT(NPNTS)        ! WORK Snow insulation factor.
*I SOILHT7A.168   
      LOGICAL
     & ITER(NPNTS)          ! WORK .T. on points requiring iterations.

!-----------------------------------------------------------------------
! Variables required for the implicit calculation.
!-----------------------------------------------------------------------
      REAL
     & DHFLUX_DTSL1(NPNTS,0:NSHYD),DHFLUX_DTSL2(NPNTS,0:NSHYD)
     &,DHADV_DTSL0(NPNTS,NSHYD),DHADV_DTSL1(NPNTS,NSHYD)
     &,DHADV_DTSL2(NPNTS,NSHYD) ! WORK Rate of change of the explicit
!                           ! fluxes with the layer temperatures
!                           ! (W/m2/K).  
     &,A(NPNTS,NSHYD),B(NPNTS,NSHYD),C(NPNTS,NSHYD),D(NPNTS,NSHYD)
!                           ! WORK Matrix elements.
     &,GAMCON               ! WORK Forward timestep weighting constant.
! Local parameters:
      REAL
     & GAMMA                ! Forward timestep weighting.
      PARAMETER (GAMMA=1.0)    
*D SOILHT7A.172
     & HEAT_CON,GAUSS                                                   
*I SOILHT7A.216   
!-----------------------------------------------------------------------
! Initialise the array of points for which calculations are required.   
!-----------------------------------------------------------------------
      DO I=1,NPNTS
        ITER(I)=.FALSE.
      ENDDO

      DO J=1,SOIL_PTS                                                   
        I=SOIL_INDEX(J)
        ITER(I)=.TRUE.                                                  
      ENDDO                  
                                                                        
*I SOILHT7A.238   
! Calculate the snow insulation factor                               
!--------------------------------------------------------------------
      DO I=1,NPNTS                                               
        SIFACT(I) = 0.                                              
      ENDDO         
      DO N=1,NTILES
        DO J=1,TILE_PTS(N)
          I = TILE_INDEX(J,N)
          SI_TILE = 1.
          IF (SNOW_TILE(I,N).GT.0. .AND. V_SAT(I).NE.0.) THEN
            SNOW_DEPTH = SNOW_TILE(I,N) / RHO_SNOW
            IF (SNOW_DEPTH .LE. (0.5*DZ(1))) THEN
              SI_TILE = 1. / ( 1. + 2*SNOW_DEPTH/(DZ(1) + DZ(2)) )
            ELSE 
              SI_TILE =(DZ(1) + DZ(2)) /
     &                 ( HC(I,1)*(2*SNOW_DEPTH - DZ(1))/SNOW_HCON
     &                   + 2*DZ(1) + DZ(2) )
            ENDIF
          ENDIF
          SIFACT(I) = SIFACT(I) + FRAC(I,N)*SI_TILE
        ENDDO
      ENDDO
                                                                        
!--------------------------------------------------------------------   
*I SOILHT7A.245   
          DHFLUX_DTSL1(I,N)=HC(I,N)*2.0/(DZ(N+1)+DZ(N))
          DHFLUX_DTSL2(I,N)=-HC(I,N)*2.0/(DZ(N+1)+DZ(N))                
*D ARE1F405.44
        H_FLUX(I,0) = SURF_HT_FLUX(I)
        H_FLUX(I,1) = SIFACT(I)*H_FLUX(I,1)
        DHFLUX_DTSL1(I,NSHYD)=0.0
        DHFLUX_DTSL2(I,NSHYD)=0.0
        DHFLUX_DTSL1(I,0)=0.0
        DHFLUX_DTSL2(I,0)=0.0
        DHFLUX_DTSL1(I,1)=SIFACT(I)*DHFLUX_DTSL1(I,1)
        DHFLUX_DTSL2(I,1)=SIFACT(I)*DHFLUX_DTSL2(I,1)               
*D SOILHT7A.262
          HADV(I,N)=HCAPW*DZ(N)* 
*I SOILHT7A.264   
          DHADV_DTSL0(I,N)=HCAPW*DZ(N)*W_FLUX(I,N-1)/(DZ(N)+DZ(N-1))    
          DHADV_DTSL1(I,N)=HCAPW*DZ(N)*                                 
     &                     (-W_FLUX(I,N-1)/(DZ(N)+DZ(N-1))
     &                      +W_FLUX(I,N)/(DZ(N)+DZ(N+1)))            
          DHADV_DTSL2(I,N)=-HCAPW*DZ(N)*W_FLUX(I,N)/(DZ(N)+DZ(N+1))     
*D SOILHT7A.271
        HADV(I,1)=HCAPW*DZ(1)*
*D SOILHT7A.274
        DHADV_DTSL0(I,1)=0.0                
        DHADV_DTSL1(I,1)=HCAPW*DZ(1)*                                   
     &    (-W_FLUX(I,0)/DZ(1)+W_FLUX(I,1)/(DZ(1)+DZ(2)))            
        DHADV_DTSL2(I,1)=-HCAPW*DZ(1)*W_FLUX(I,1)/(DZ(1)+DZ(2))
        HADV(I,NSHYD)=HCAPW*DZ(NSHYD)*
*I SOILHT7A.276   
        DHADV_DTSL0(I,NSHYD)=HCAPW*DZ(NSHYD)*W_FLUX(I,NSHYD-1)
     &                       /(DZ(NSHYD)+DZ(NSHYD-1))                
        DHADV_DTSL1(I,NSHYD)=-HCAPW*DZ(NSHYD)*W_FLUX(I,NSHYD-1)
     &                       /(DZ(NSHYD)+DZ(NSHYD-1))
        DHADV_DTSL2(I,NSHYD)=0.0                           
*D ACB2F405.15,ACB2F405.17   
! unfrozen. Check that (SMCLSAT/SMCL)**BEXP will not overflow when SMCL 
! is very small. The function EPSILON  gives the number of type (real)  
! of SMCL that is negligeable compared to 1.                            
*D SOILHT7A.289
     &             *(SMCLSAT(I,N)/SMCL(I,N))**(BEXP(I))                 
*D SOILHT7A.295,SOILHT7A.297  
          DHSL0(I,N)=TIMESTEP*(H_FLUX(I,N-1)-H_FLUX(I,N)+HADV(I,N))     
*I SOILHT7A.300   
        ENDDO
      ENDDO                                          
                                                                        
*D SOILHT7A.302
! Iteration loop                                                     
*I SOILHT7A.303   
      DO M=1,MMAX                                                       
      

!-----------------------------------------------------------------------
! Define the array of points which fail to meet the flux criterion.     
!-----------------------------------------------------------------------
      ITER_PTS=0                                                        
      DO J=1,SOIL_PTS                                                   
        I=SOIL_INDEX(J)                                                 
                                                                        
        IF (ITER(I)) THEN                                   
          ITER_PTS=ITER_PTS+1                                 
          ITER_INDEX(ITER_PTS)=I                                        
        ENDIF                                                           
        ITER(I)=.FALSE.
                                                                        
      ENDDO                                                             

!-----------------------------------------------------------------------
! Update calculations at these points.  
!-----------------------------------------------------------------------
      DO N=1,NSHYD

        DO J=1,ITER_PTS
          I=ITER_INDEX(J)

*I SOILHT7A.307   
          DTSLMAX(I,N)=1.0E4-TSL(I,N)
          DTSLMIN(I,N)=-ZERODEGC-TSL(I,N)

*I SOILHT7A.309   
            DTSLMIN(I,N)=TMAX(I,N)-TSL(I,N)  
*D SOILHT7A.315,SOILHT7A.316  
     &               /(BEXP(I)*SATHH(I)*RHO_WATER*DZ(N))              
     &            *(-DPSIDT*TSL(I,N)/SATHH(I))**(-1.0/BEXP(I)-1.0)     
            DTSLMAX(I,N)=TMAX(I,N)-TSL(I,N)
*D SOILHT7A.319
          HCAPT(I,N)=HCAP(I)+(HCAPW-HCAPI)*SMCLU(I,N)/DZ(N)         
*D SOILHT7A.323,SOILHT7A.324  
*D SOILHT7A.327,SOILHT7A.328  
! Calculate the matrix elements required for the implicit update.
*D SOILHT7A.330,SOILHT7A.372  
          GAMCON=GAMMA*TIMESTEP/(HCAPT(I,N)*DZ(N))
          A(I,N)=-GAMCON*(DHFLUX_DTSL1(I,N-1)+DHADV_DTSL0(I,N))
          B(I,N)=1.0-GAMCON*(DHFLUX_DTSL2(I,N-1)-DHFLUX_DTSL1(I,N)
     &                                          +DHADV_DTSL1(I,N))
          C(I,N)=GAMCON*(DHFLUX_DTSL2(I,N)+DHADV_DTSL2(I,N))
          D(I,N)=1.0/(HCAPT(I,N)*DZ(N))*DHSL(I,N)                       
*D SOILHT7A.375,SOILHT7A.389  
      ENDDO
*D SOILHT7A.393
! Solve the triadiagonal matrix equation.
*D SOILHT7A.395,SOILHT7A.398  
      CALL GAUSS(NSHYD,NPNTS,ITER_PTS,ITER_INDEX,A,B,C,D
     &,          DTSLMIN,DTSLMAX,DTSL)
*D SOILHT7A.400,SOILHT7A.402  
!-----------------------------------------------------------------------
! Diagnose the implicit DHSL
!-----------------------------------------------------------------------
      DO N=2,NSHYD-1
        DO J=1,ITER_PTS
          I=ITER_INDEX(J)
          DHSL(I,N)=DHSL(I,N)-DZ(N)*HCAPT(I,N)*(A(I,N)*DTSL(I,N-1)
     &                    +(B(I,N)-1)*DTSL(I,N)+C(I,N)*DTSL(I,N+1)) 
        ENDDO
      ENDDO
*D SOILHT7A.404,SOILHT7A.413  
      DO J=1,ITER_PTS
        I=ITER_INDEX(J)
        DHSL(I,1)=DHSL(I,1)-DZ(1)*HCAPT(I,1)*(
     &                  +(B(I,1)-1)*DTSL(I,1)+C(I,1)*DTSL(I,2)) 
        DHSL(I,NSHYD)=DHSL(I,NSHYD)-DZ(NSHYD)*HCAPT(I,NSHYD)*
     &  (A(I,NSHYD)*DTSL(I,NSHYD-1)+(B(I,NSHYD)-1)*DTSL(I,NSHYD)) 
      ENDDO
*D SOILHT7A.415,SOILHT7A.417  
!-----------------------------------------------------------------------
! Update the layer temperatures
!-----------------------------------------------------------------------
      DO N=1,NSHYD
        DO J=1,ITER_PTS
          I=ITER_INDEX(J)
*D SOILHT7A.419
*D SOILHT7A.436,SOILHT7A.445  
*D SOILHT7A.453
     &                *(-DPSIDT*TSL(I,N)/SATHH(I))**(-1.0/BEXP(I))     
*D SOILHT7A.465
!-----------------------------------------------------------------------
! Calculate the error in flux conservation                              
!-----------------------------------------------------------------------
          CEACUR(I)=ABS(DHSL(I,N))/TIMESTEP                          

          IF (CEACUR(I) .GT. FACUR) THEN                                
            ITER(I)=.TRUE.
          ENDIF

*I SOILHT7A.466   
      ENDDO
                                                                        
!-----------------------------------------------------------------------
! End of iteration loop  
!-----------------------------------------------------------------------
*DECLARE SOILHY5A
*D SOILHY5A.23,SOILHY5A.25   
      SUBROUTINE SOIL_HYD (NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX              
     &,                    BEXP,DZ,EXT,FW,KS,SATHH,TIMESTEP,V_SAT       
     &,                    SLOW_RUNOFF,SMCL,STHU,SURF_ROFF,W_FLUX       
     &,                    STF_SLOW_RUNOFF,LTIMER                       
*I SOILHY5A.48    
!  4.6      2/99     Modified for implicit updating.  Peter Cox         
*D SOILHY5A.82
     & BEXP(NPNTS)          ! IN Clapp-Hornberger exponent.             
*I SOILHY5A.111   
     &,SURF_ROFF(NPNTS)     ! INOUT Surface runoff (kg/m2/s).           
*I SOILHY5A.117   
      REAL                                                              
     & GAMCON               ! WORK Constant (s/mm).                     
                                                                        
*D SOILHY5A.120,SOILHY5A.124  
     & A(NPNTS,NSHYD)       ! WORK Matrix elements corresponding        
!                           !      to the coefficients of DSTHU(n-1).   
     &,B(NPNTS,NSHYD)       ! WORK Matrix elements corresponding        
!                           !      to the coefficients of DSTHU(n).     
     &,C(NPNTS,NSHYD)       ! WORK Matrix elements corresponding        
!                           !      to the coefficients of DSTHU(n+1).   
     &,D(NPNTS,NSHYD)       ! WORK Matrix elements corresponding        
!                           !      to the RHS of the equation.          
     &,DSMCL(NPNTS,NSHYD)   ! WORK Soil moisture increment              
!                           !      (kg/m2/timestep).                    
     &,DSTHU(NPNTS,NSHYD)   ! WORK Increment to STHU (/timestep).       
     &,DSTHUMIN(NPNTS,NSHYD)! WORK Minimum value of DSTHU.              
     &,DSTHUMAX(NPNTS,NSHYD)! WORK Maximum value of DSTHU.              
     &,DWFLUX_DSTHU1(NPNTS,NSHYD) ! WORK The rate of change of the expli
!                           !      flux with STHU1 (kg/m2/s).           
     &,DWFLUX_DSTHU2(NPNTS,NSHYD) ! WORK The rate of change of the expli
!                           !      flux with STHU2 (kg/m2/s).           
*D SOILHY5A.127,SOILHY5A.128  
                                                                        
! Local parameters:                                                     
      REAL                                                              
     & GAMMA                ! Forward timestep weighting.               
      PARAMETER (GAMMA=1.0)                                             
*D SOILHY5A.141,GRB0F405.538  
*I SOILHY5A.146   
          DSTHUMIN(I,N)=-STHU(I,N)                                      
          DSTHUMAX(I,N)=1.0-SMCL(I,N)/SMCLSAT(I,N)                      
          DWFLUX_DSTHU1(I,N)=0.0                                        
          DWFLUX_DSTHU2(I,N)=0.0                                        
*D SOILHY5A.159
! Calculate the Darcian fluxes and their dependencies on the soil       
! moisture contents.                                                    
*I SOILHY5A.160   
      CALL HYD_CON (NPNTS,SOIL_PTS,SOIL_INDEX,BEXP,KS,STHU(1,NSHYD)     
     &,             W_FLUX(1,NSHYD),DWFLUX_DSTHU1(1,NSHYD)              
     &,             LTIMER)                                             
                                                                        
*I SOILHY5A.161   
        CALL DARCY (NPNTS,SOIL_PTS,SOIL_INDEX,BEXP,KS,SATHH             
     &,             STHU(1,N-1),DZ(N-1),STHU(1,N),DZ(N),W_FLUX(1,N-1)   
     &,             DWFLUX_DSTHU1(1,N-1),DWFLUX_DSTHU2(1,N-1)           
     &,             LTIMER)                                             
       ENDDO                                                            
                                                                        
!-----------------------------------------------------------------------
! Calculate the explicit increments.                                    
!-----------------------------------------------------------------------
ccccc runoff_MII mod cccc
c     DO N=NSHYD,1,-1
      DO N=1,NSHYD,1
ccccccccccccccccccccccccc
*D SOILHY5A.164,SOILHY5A.166  
          DSMCL(I,N)=(W_FLUX(I,N-1)-W_FLUX(I,N)-EXT(I,N))*TIMESTEP      
*D SOILHY5A.169,SOILHY5A.170  
! Limit the explicit fluxes to prevent supersaturation.                 
*D SOILHY5A.172,SOILHY5A.275  
          IF (DSMCL(I,N).GT.(SMCLSAT(I,N)-SMCL(I,N))) THEN              
            DSMCL(I,N)=SMCLSAT(I,N)-SMCL(I,N)                           
ccccc runoff_MII mod cccc
c           W_FLUX(I,N-1)=DSMCL(I,N)/TIMESTEP+W_FLUX(I,N)+EXT(I,N)
            W_FLUX(I,N)=W_FLUX(I,N-1)-DSMCL(I,N)/TIMESTEP-EXT(I,N)
ccccccccccccccccccccccccc
*D SOILHY5A.282,SOILHY5A.283  
! Calculate the matrix elements required for the implicit update.       
*D SOILHY5A.285,GRB0F405.548  
*D SOILHY5A.288,SOILHY5A.289  
        GAMCON=GAMMA*TIMESTEP/SMCLSAT(I,1)                              
        A(I,1)=0.0                                                      
        B(I,1)=1.0+GAMCON*DWFLUX_DSTHU1(I,1)                            
        C(I,1)=GAMCON*DWFLUX_DSTHU2(I,1)                                
        D(I,1)=DSMCL(I,1)/SMCLSAT(I,1)                                  
      ENDDO                                                             
*D SOILHY5A.291
      DO N=2,NSHYD                                                      
        DO J=1,SOIL_PTS                                                 
          I=SOIL_INDEX(J)                                               
          GAMCON=GAMMA*TIMESTEP/SMCLSAT(I,N)                            
          A(I,N)=-GAMCON*DWFLUX_DSTHU1(I,N-1)                           
          B(I,N)=1.0-GAMCON*(DWFLUX_DSTHU2(I,N-1)-DWFLUX_DSTHU1(I,N))   
          C(I,N)=GAMCON*DWFLUX_DSTHU2(I,N)                              
          D(I,N)=DSMCL(I,N)/SMCLSAT(I,N)                                
        ENDDO                                                           
      ENDDO                                                             
*D SOILHY5A.293,SOILHY5A.295  
!-----------------------------------------------------------------------
! Solve the triadiagonal matrix equation.                               
!-----------------------------------------------------------------------
      CALL GAUSS(NSHYD,NPNTS,SOIL_PTS,SOIL_INDEX,A,B,C,D                
     &,          DSTHUMIN,DSTHUMAX,DSTHU)                               
*D SOILHY5A.297
!-----------------------------------------------------------------------
! Diagnose the implicit fluxes.                                         
!-----------------------------------------------------------------------
      DO N=1,NSHYD                                                      
        DO J=1,SOIL_PTS                                                 
          I=SOIL_INDEX(J)                                               
          DSMCL(I,N)=DSTHU(I,N)*SMCLSAT(I,N)                            
          W_FLUX(I,N)=W_FLUX(I,N-1)-EXT(I,N)-DSMCL(I,N)/TIMESTEP        
        ENDDO                                                           
      ENDDO                                                             
                                                                        
!-----------------------------------------------------------------------
! Update the prognostic variables.                                      
!-----------------------------------------------------------------------
      DO N=1,NSHYD                                                      
        DO J=1,SOIL_PTS                                                 
          I=SOIL_INDEX(J)                                               
          SMCLU(I,N)=SMCLU(I,N)+DSMCL(I,N)                              
          SMCL(I,N)=SMCL(I,N)+DSMCL(I,N)                                
          STHU(I,N)=SMCLU(I,N)/SMCLSAT(I,N)                             
        ENDDO                                                           
*I SOILHY5A.311   
                                                                        
!-----------------------------------------------------------------------
! Update surface runoff diagnostic.                                     
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS                                                   
        I=SOIL_INDEX(J)                                                 
        SURF_ROFF(I)=SURF_ROFF(I)+(FW(I)-W_FLUX(I,0))                   
      ENDDO                                                             
*DECLARE SOILMC7A
*D SOILMC7A.77
            SMC(I) = SMC(I) + RHO_WATER * ( ZSMC - Z1 ) * 
*DECLARE SPARM1A
*D SPARM1A.27,SPARM1A.29   
      SUBROUTINE SPARM (LAND_FIELD,LAND1,LAND_PTS,NTILES
     &,                 TILE_PTS,TILE_INDEX,FRAC,HT,LAI,SATCON
     &,                 CATCH_TILE,INFIL_TILE,Z0_TILE)           
*I SPARM1A.38    
     &,NTILES                ! IN Number of surface tiles.
*D SPARM1A.45,SPARM1A.46   
     & FRAC(LAND_FIELD,NTYPE)     ! IN Fractional cover of each       
*I SPARM1A.49    
     &,SATCON(LAND_FIELD)         ! IN Saturated hydraulic
!                                 !    conductivity (kg/m2/s).        
*D SPARM1A.52,SPARM1A.58   
     & CATCH_TILE(LAND_FIELD,NTILES)! OUT Canopy capacity for each tile 
!                                   !     (kg/m2).
     &,INFIL_TILE(LAND_FIELD,NTILES)! OUT Maximum surface infiltration 
!                                   !     rate for each tile (kg/m2/s).
     &,Z0_TILE(LAND_FIELD,NTILES)   ! OUT Roughness length for each     
!                                   !     tile (m).                     
*D SPARM1A.60,SPARM1A.63   
     & CATCH(LAND_FIELD)          ! WORK GBM canopy capacity (kg/m2). 
     &,CATCH_T(LAND_FIELD,NTYPE)  ! WORK Capacities for types.
*I SPARM1A.64    
     &,INFIL(LAND_FIELD)          ! WORK GBM infiltration rate(kg/m2/s).
     &,INFIL_T(LAND_FIELD,NTYPE)  ! WORK Infiltration rates for types.
     &,Z0(LAND_FIELD)             ! WORK GBM roughness length (m).
     &,Z0_T(LAND_FIELD,NTYPE)     ! WORK Roughness lengths for types.   
*D SPARM1A.68,SPARM1A.74   
*D SPARM1A.85,SPARM1A.87   
     &,                 HT(1,N),LAI(1,N),SATCON
     &,                 CATCH_T(1,N),INFIL_T(1,N),Z0_T(1,N))
*D SPARM1A.96,SPARM1A.98   
          CATCH_T(L,N) = CATCH_NVG(N-NPFT)
          INFIL_T(L,N) = INFIL_NVG(N-NPFT)*SATCON(L)
*D SPARM1A.103,SPARM1A.106  
      IF (NTILES .EQ. 1) THEN
!----------------------------------------------------------------------
! Form means and copy to tile arrays if required for aggregate tiles
!----------------------------------------------------------------------
        DO L=1,LAND_FIELD
          CATCH(L) = 0.0
          FZ0(L) = 0.0
          INFIL(L) = 0.0
          Z0(L) = 0.0
*D SPARM1A.108
*I SPARM1A.109   
        DO N=1,NTYPE
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            FZ0(L) = FZ0(L) + FRAC(L,N) / (LOG(LB / Z0_T(L,N)))**2
          ENDDO
        ENDDO
        DO L=LAND1,LAND1+LAND_PTS-1
          Z0(L) = LB * EXP(-SQRT(1. / FZ0(L)))
        ENDDO

        DO N=1,NTYPE
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            CATCH(L) = CATCH(L) + FRAC(L,N) * CATCH_T(L,N)
            INFIL(L) = INFIL(L) + FRAC(L,N) * INFIL_T(L,N)
          ENDDO
        ENDDO

        DO L=1,LAND_FIELD
!         Canopy capacity is average over non-lake surface types
          CATCH_TILE(L,1) = 0.
          IF (FRAC(L,7).LT.1.)
     &      CATCH_TILE(L,1) = CATCH(L) / (1. - FRAC(L,7))
          INFIL_TILE(L,1) = INFIL(L)
          Z0_TILE(L,1) = Z0(L)
        ENDDO  

      ELSE
*D SPARM1A.111
! Copy surface-type arrays to tiles if separate tiles used
*D SPARM1A.113,SPARM1A.118  
        DO N=1,NTYPE
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            CATCH_TILE(L,N) = CATCH_T(L,N)
            INFIL_TILE(L,N) = INFIL_T(L,N)
            Z0_TILE(L,N) = Z0_T(L,N)
          ENDDO
        ENDDO
*D SPARM1A.120,SPARM1A.141  
      ENDIF
*DECLARE STATMPT1
*I GSM3F404.7     

! Coastal tiling fields
      JFRAC_LAND     = SI(505,Sect_No,im_index)                         
      JTSTAR_LAND    = SI(506,Sect_No,im_index)                         
      JTSTAR_SEA     = SI(507,Sect_No,im_index)                         
      JTSTAR_SICE    = SI(508,Sect_No,im_index)                         
      JSICE_ALB      = SI(509,Sect_No,im_index)                         
      JLAND_ALB      = SI(510,Sect_No,im_index)                         
                             
*D ABX1F404.74,ABX1F404.78   
      JCAN_WATER_TYP= SI(229,Sect_No,im_index) ! Canopy water content   
C                                              ! on tiles               
      JCATCH_TYP    = SI(230,Sect_No,im_index) ! Canopy capacity on     
C                                              ! tiles                  
      JRGRAIN_TYP   = SI(231,Sect_No,im_index) ! Snow grain size on     
C                                              ! tiles                  
*I ABX1F404.81    
      JSNODEP_TYP   = SI(235,Sect_No,im_index) ! Tiled snow depth       
      JINFIL_TYP    = SI(236,Sect_No,im_index) ! Max tile infilt rate   
*DECLARE STDEV7A
*D STDEV7A.2
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
*D STDEV7A.34
     & P_POINTS,P_FIELD,P1,FLANDG,                                   
*D STDEV7A.49
*D STDEV7A.52
     & FLANDG(P_FIELD)       ! IN Land fraction.                   
     &,BQ_1(P_FIELD)         ! IN Buoyancy parameter.                   
*D STDEV7A.91
        IF ( FLANDG(I).LT.1.0 ) THEN                                   
*D STDEV7A.100,STDEV7A.101  
            T1_SD(I) = MAX ( 0.0 , 
     &          (1.-FLANDG(I))*1.93*FTL_1(I) / (RHOSTAR(I)*WS1) )   
            Q1_SD(I) = MAX ( 0.0 , 
     &          (1.-FLANDG(I))*1.93*FQW_1(I) / (RHOSTAR(I)*WS1) )   
*D STDEV7A.118,STDEV7A.119  
     & P_FIELD,LAND_FIELD,TILE_PTS,LAND_INDEX,TILE_INDEX,FLAND,         
     & BQ_1,BT_1,FQW_1,FTL_1,RHOKM_1,RHOSTAR,VSHR,Z0M,Z1_TQ,TILE_FRAC,  
*D STDEV7A.136
     & FLAND(LAND_FIELD)                                                
     &,BQ_1(P_FIELD)         ! IN Buoyancy parameter.                   
*I STDEV7A.145   
     &,TILE_FRAC(LAND_FIELD) ! IN Tile fraction.                        
*D STDEV7A.182,STDEV7A.188  
          T1_SD(I) = T1_SD(I) + MAX ( 0.0 ,                             
     &      FLAND(L)*TILE_FRAC(L)*1.93*FTL_1(L) / (RHOSTAR(I)*WS1) )    
          Q1_SD(I) = Q1_SD(I) + MAX ( 0.0 ,                             
     &      FLAND(L)*TILE_FRAC(L)*1.93*FQW_1(L) / (RHOSTAR(I)*WS1) )    
*DECLARE SURFHY7A
*D SURFHY7A.55,SURFHY7A.56   
      SUBROUTINE SURF_HYD (NPNTS,NTILES,TILE_PTS,TILE_INDEX,            
*D SURFHY7A.58
     &                     LS_RAIN,MELT_TILE,SNOW_MELT,TIMESTEP,        
*D SURFHY7A.66,SURFHY7A.68   
     &,NTILES               ! IN Number of tiles.                       
     &,TILE_PTS(NTILES)     ! IN Number of tile points.                 
     &,TILE_INDEX(NPNTS,NTILES)                                         
*D SURFHY7A.70,SURFHY7A.71   
*D SURFHY7A.74,SURFHY7A.79   
     & CAN_CPY(NPNTS,NTILES)! IN Canopy capacity for land tiles (kg/m2).
     &,E_CANOPY(NPNTS,NTILES)!IN Canopy evaporation (kg/m2/s).
     &,FRAC(NPNTS,NTILES)   ! IN Tile fractions.                        
     &,INFIL(NPNTS,NTILES)  ! IN Infiltration rate (kg/m2/s).           
*D SURFHY7A.82,SURFHY7A.83   
     &,MELT_TILE(NPNTS,NTILES)
!                           ! IN Snow melt on tiles (kg/m2/s).        
     &,SNOW_MELT(NPNTS)     ! IN GBM snow melt (kg/m2/s).               
*D SURFHY7A.87,SURFHY7A.89   
     & CAN_WCNT(NPNTS,NTILES)!INOUT Tile canopy water contents (kg/m2).
*D SURFHY7A.118,SURFHY7A.135  
*D SURFHY7A.137
      DO N=1,NTILES                                                   
*D SURFHY7A.146
      DO N=1,NTILES

! Surface runoff of snowmelt, assumed to cover 100% of tile
        CALL FRUNOFF (NPNTS,TILE_PTS(N),TILE_INDEX(1,N),1.,             
     &                CAN_CPY(1,N),CAN_CPY(1,N),INFIL(1,N),   
     &                MELT_TILE(1,N),FRAC(1,N),TIMESTEP,                
     &                SURF_ROFF) 
*D SURFHY7A.188
        DO J=1,TILE_PTS(N)                                              
          I = TILE_INDEX(J,N)                                       
*DECLARE SWRAD3A
*I ADB2F404.1504  
!       4.6             10-05-98                Land flag passed to     
!                                               FLUX_CALC.              
!                                               (J. M. Edwards)         
*D ARE2F404.239,ARE2F404.241  
     &   , LAND_ALBEDO, L_CTILE                                         
     &   , LAND_ALB, SICE_ALB, FLANDG                                   
     &   , OPEN_SEA_ALBEDO, ICE_FRACTION, LAND, LAND0P5, LYING_SNOW
!                       MOSES II flag and array dimension               
     &   , L_MOSES_II, SAL_DIM                                          
*D SWRAD3A.46
     &   , FLUX_BELOW_690NM_SURF,  FL_SOLID_BELOW_690NM_SURF            
     &   , FL_SEA_BELOW_690NM_SURF, L_FLUX_BELOW_690NM_SURF             
*I SWRAD3A.74    
     &   , SURF_DOWN_SW                                                 
*D ARE2F404.242,ARE2F404.243  
     &   , LAND0P5(NPD_FIELD)                                           
!             LAND MASK (TRUE if land fraction >0.5)                    
c     &   , SEA(NPD_FIELD)                                              
c!             SEA MASK                                           
*D ARE2F404.246
!             DIMENSION FOR LAND_ALBEDO (MOSES II)                      
*D SWRAD3A.255,ADB1F401.1045 
!             FRACTION OF SEA ICE IN SEA PORTION OF GRID BOX            
     &   , LAND_ALBEDO(SAL_DIM,4)                                       
!             MOSES II LAND SURFACE ALBEDO FIELDS                       
     &   , LAND_ALB(NPD_FIELD)                                          
!             SURFACE ALBEDO OF LAND                                    
     &   , SICE_ALB(NPD_FIELD)                                          
!             SURFACE ALBEDO OF SEA-ICE                                 
     &   , FLANDG(NPD_FIELD)                                            
!             Fractional land                                           
     &   , FLANDG_G(NPD_PROFILE)                                        
!             Gathered Fractional land                                  
     &   , ICE_FRACTION_G(NPD_PROFILE)                                  
!             Gathered Fractional sea-ice                               
*I SWRAD3A.276   
!             WEIGHTED BY (OPEN SEA)/(TOTAL SEA) FRACTION               
*I SWRAD3A.278   
     &   , SURF_DOWN_SW(NPD_FIELD, 4)                                   
!             SURFACE DOWNWARD SHORTWAVE RADIATION COMPONENTS           
!             (*,1) - DIRECT BEAM VISIBLE                               
!             (*,2) - DIFFUSE VISIBLE                                   
!             (*,3) - DIRECT BEAM NEAR-IR                               
!             (*,4) - DIFFUSE NEAR-IR                                   
*I ADB1F401.1048  
     &   , L_MOSES_II                                                   
!             SURFACE SW FLUXES REQUIRED FOR MOSES II                   
     &   , L_CTILE                                                      
!             SWITCH FOR COASTAL TILING                                 
*I ADB1F401.1050  
!             NB: ONLY USED FOR NON MOSESII RUNS.                       
     &   , FL_SOLID_BELOW_690NM_SURF(NPD_FIELD)                         
!             SOLID SFC NET SURFACE FLUX BELOW 690NM                    
     &   , FL_SEA_BELOW_690NM_SURF(NPD_FIELD)                           
!             OPEN SEA NET SURFACE FLUX BELOW 690NM                     
!             (AT POINTS WHERE THERE THIS IS SEA-ICE THIS IS            
!             WEIGHTED BY THE FRACTION OF OPEN SEA.)                    
*I SWRAD3A.432   
     &   , LAND0P5_G(NPD_PROFILE)                                       
!             GATHERED LAND MASK (TRUE if land fraction >0.5)           
*D SWRAD3A.485
!             GATHERED GRID BOX MEAN DOWNWARD SURFACE FLUX              
!             BELOW 690NM                                               
     &   , FL_SEA_BELOW_690NM_SURF_G(NPD_PROFILE)                       
!             GATHERED OPEN SEA NET SURFACE FLUX BELOW 690NM            
*I SWRAD3A.486   
!     TEMPORARY FIELDS ASSOCIATED WITH 690NM FLUX OVER                  
!     MEAN SOLID SURF.                                                  
      REAL                                                              
     &     ALBEDOSOLID                                                  
     &    ,FRACSOLID                                                    
     &    ,FLUXSOLID                                                    
!                                                                       
     &   , SURF_VIS_DIR_G(NPD_PROFILE)                                  
!             GATHERED DOWNWARD SURFACE DIRECT BEAM VISIBLE FLUX        
     &   , SURF_VIS_DIF_G(NPD_PROFILE)                                  
!             GATHERED DOWNWARD SURFACE DIFFUSE VISIBLE FLUX            
     &   , SURF_NIR_DIR_G(NPD_PROFILE)                                  
!             GATHERED DOWNWARD SURFACE DIRECT BEAM NEAR-INFRARED FLUX  
     &   , SURF_NIR_DIF_G(NPD_PROFILE)                                  
!             GATHERED DOWNWARD SURFACE DIFFUSE NEAR-INFRARED FLUX      
!                                                                       
*D ADB1F401.1069,ADB1F401.1070 
     &     N_FRAC_SOL_POINT                                             
     &   , I_FRAC_SOL_POINT(NPD_PROFILE)                                
*D SWRAD3A.530,SWRAD3A.531  
*D ARE2F404.252
      IF ( L_FLUX_BELOW_690NM_SURF .OR. L_MOSES_II ) THEN               
*D ARE2F404.264,ARE2F404.265  
     &   , L_MICROPHYSICS, L_MOSES_II, SAL_DIM, L_CTILE                 
     &   , LAND, LAND0P5, OPEN_SEA_ALBEDO                               
     &   , LAND_ALB, SICE_ALB                                           
     &   , FLANDG, ICE_FRACTION                                         
     &   , LAND_ALBEDO, WEIGHT_690NM                                    
*D SWRAD3A.545
     &   , LAND_G, LAND0P5_G, FLANDG_G, ICE_FRACTION_G                  
     &   , ALBEDO_SEA_DIFF_G, ALBEDO_SEA_DIR_G                          
*D ADB2F404.1571
     &   , PSTAR, DUMMY, DUMMY, DUMMY
     &   , AB, BB, AC, BC                                 
*D ADB2F404.1573
     &   , P, T, DUMMY, DUMMY, DUMMY, DUMMY, D_MASS                     
*D ADB1F402.725
     &      , LAND0P5, LYING_SNOW, PSTAR, AB, BB, TRINDX                
*D AYY1F404.372
     &   , L_CLOUD_WATER_PARTITION,  LAND0P5_G                          
*D ADB2F404.1604
     &   , P, T, DUMMY, DUMMY, DUMMY, DUMMY, D_MASS                     
*D ADB1F401.1105,SWRAD3A.760  
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION_G           
     &   , ALBEDO_SEA_DIFF_G, ALBEDO_SEA_DIR_G, FLANDG_G                
*D SWRAD3A.764
     &   , FLUX_BELOW_690NM_SURF_G, FL_SEA_BELOW_690NM_SURF_G           
     &   , L_MOSES_II, L_CTILE                                          
     &   , SURF_VIS_DIR_G,SURF_VIS_DIF_G,SURF_NIR_DIR_G,SURF_NIR_DIF_G  
*I SWRAD3A.795   
      IF (L_MOSES_II) THEN                                              
         DO I=1, 4                                                      
            CALL R2_ZERO_1D(N_PROFILE, SURF_DOWN_SW(1, I))              
         ENDDO                                                          
      ENDIF                                                             
*I SWRAD3A.891   
        IF(L_CTILE)THEN                                                 
           CALL R2_ZERO_1D(N_PROFILE, FL_SOLID_BELOW_690NM_SURF)        
           CALL R2_ZERO_1D(N_PROFILE, FL_SEA_BELOW_690NM_SURF)          
          DO L=1, NLIT                                                  
            IF (FLANDG(LIST(L)).LT.1.0) THEN                            
              FL_SEA_BELOW_690NM_SURF(LIST(L))                          
     &          =FL_SEA_BELOW_690NM_SURF_G(L)                           
     &         *(1.0E+00-ICE_FRACTION(LIST(L)))                         
            ENDIF                                                       
            FRACSOLID=FLANDG(LIST(L))                                   
     &        +(1.-FLANDG(LIST(L)))*ICE_FRACTION(LIST(L))               
            IF(FRACSOLID.GT.0.0)THEN                                    
              FL_SOLID_BELOW_690NM_SURF(LIST(L))=                       
     &         (FLUX_BELOW_690NM_SURF_G(L)-                             
     &         (1.-FLANDG(LIST(L)))*FL_SEA_BELOW_690NM_SURF(LIST(L)))   
     &         /FRACSOLID                                               
            ENDIF                                                       
          ENDDO                                                         
        ELSE                                                            
*I ADB1F401.1111  
                                                                        
*I SWRAD3A.895   
        ENDIF                                                           
                                                                        
      ENDIF                                                             
!                                                                       
!                                                                       
!     COMPONENTS OF DOWNWARD FLUX AT THE SURFACE FOR MOSES II           
!                                                                       
      IF (L_MOSES_II) THEN                                              
        DO L=1, NLIT                                                    
           SURF_DOWN_SW(LIST(L),1) = SURF_VIS_DIR_G(L)                  
           SURF_DOWN_SW(LIST(L),2) = SURF_VIS_DIF_G(L)                  
           SURF_DOWN_SW(LIST(L),3) = SURF_NIR_DIR_G(L)                  
           SURF_DOWN_SW(LIST(L),4) = SURF_NIR_DIF_G(L)                  
        ENDDO                                                           
*I SWRAD3A.957   
! Weight open sea flux with open sea fraction over total sea            
      IF(L_CTILE)THEN                                                   
        DO L=1, NLIT                                                    
          IF (FLANDG(LIST(L)).LT.1.0) THEN                              
            SWSEA(LIST(L))=(1.0E+00-ICE_FRACTION(LIST(L)))              
     &         *SEA_FLUX_G(L)                                           
            SWOUT(LIST(L), 1)=SWOUT(LIST(L), 1)                         
     &        -(1.0E+00-FLANDG(LIST(L)))*SWSEA(LIST(L))                 
           ELSE                                                         
            SWSEA(LIST(L))=0.0                                          
           ENDIF                                                        
         ENDDO                                                          
       ENDIF                                                            
       IF(.NOT.L_CTILE)THEN                                             
*I SWRAD3A.965   
      ENDIF                                                             
*D AJS1F401.1424
!     TOTAL DOWNWARD FLUX OF PHOTOSYTHETICALLY ACTIVE RADIATION         
*I AJS1F401.1428  
        IF (L_MOSES_II) THEN                                            
          DO L=1, NLIT                                                  
            SWOUT(LIST(L),NLEVS+2) = SURF_VIS_DIR_G(L) +                
     &                               SURF_VIS_DIF_G(L)                  
          ENDDO                                                         
        ELSE       
*D AJS1F401.1430,AJS1F401.1432 
             IF(LAND(L))THEN                                            
               SWOUT(L, NLEVS+2)=FLUX_BELOW_690NM_SURF(L) /             
     &           (1 - LAND_ALB(L))                                      
              ELSE                                                      
               SWOUT(L, NLEVS+2)=FLUX_BELOW_690NM_SURF(L) /             
     &           (1 - SICE_ALB(L))                                      
              ENDIF                                                     
                                                                        
          ENDDO                                                         
        ENDIF                                                           
*I SWRAD3A.976   
!                                                                       
!     DIVIDE SURFACE DOWNWARD SW COMPONENTS BY COSINE OF SOLAR ZENITH   
!     ANGLE FOR MOSES II.                                               
      IF (L_MOSES_II) THEN                                              
        DO I=1, 4                                                       
          DO L=1, N_PROFILE                                             
            SURF_DOWN_SW(L, I) = SURF_DOWN_SW(L, I) /                   
     &                          (COSZIN(L)*LIT(L)+TOL_MACHINE)          
          ENDDO                                                         
        ENDDO                                                           
      ENDIF                                                             
*D ARE2F404.273,ARE2F404.274  
     &   , L_MICROPHYSICS, L_MOSES_II, SAL_DIM, L_CTILE                 
     &   , LAND, LAND0P5, OPEN_SEA_ALBEDO                               
     &   , LAND_ALB, SICE_ALB                                           
     &   , FLANDG, ICE_FRACTION                                         
     &   , LAND_ALBEDO, WEIGHT_690NM                                    
*D SWRAD3A.1008
     &   , LAND_G, LAND0P5_G, FLANDG_G, ICE_FRACTION_G                  
     &   , ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR                              
*D ARE2F404.276
!             DIMENSION OF LAND_ALBEDO                                  
*I SWRAD3A.1057  
     &   , LAND0P5(NPD_FIELD)                                           
!             LAND MASK (TRUE if land fraction >0.5)                    
*D SWRAD3A.1061,ARE2F404.280  
     &   , FLANDG(NPD_FIELD)                                            
     &   , FLANDG_G(NPD_PROFILE)                                        
!             GATHERED LAND FRACTION                                    
     &   , ICE_FRACTION_G(NPD_PROFILE)                                  
!             GATHERED SEA-ICE FRACTION IN SEA PORTION OF GRID BOX      
     &   , LAND_ALB(NPD_FIELD)                                          
     &   , SICE_ALB(NPD_FIELD)    
     &   , LAND_ALBEDO(SAL_DIM,4)          
!             MOSES II LAND SURFACE ALBEDO FIELDS                       
*D SWRAD3A.1064
!             FRACTION OF SEA ICE IN SEA PORTION OF GRID BOX            
*D ARE2F404.283,ARE2F404.284  
     &   , L_MOSES_II                                                   
     &   , L_CTILE                                                      
!             FLAG FOR COASTAL TILING                                   
*I SWRAD3A.1084  
     &   , LAND0P5_G(NPD_PROFILE)                                       
!             GATHERED LAND MASK (TRUE if land fraction >0.5)           
*I SWRAD3A.1113  
            LAND0P5_G(L)=LAND0P5(LIST(L))                               
            FLANDG_G(L)=FLANDG(LIST(L))                                 
            ICE_FRACTION_G(L)=ICE_FRACTION(LIST(L))                     
*D SWRAD3A.1119,SWRAD3A.1121 
!     THERE IS A COMBINATION OF OPEN SEA, SEA-ICE AND LAND. SEPARATE    
!     ALBEDOS ARE PROVIDED FOR FOR OPEN SEA. BAND-DEPENDENT COPIES      
!     OF THE ALBEDOS MUST BE MADE FOR CALCULATING COUPLING FLUXES.      
*D SWRAD3A.1128
            ALBEDO_SEA_DIR(L, I)=0.0E+00                                
            ALBEDO_SEA_DIFF(L, I)=0.0E+00                               
            ALBEDO_FIELD_DIFF(L, I)=0.0E+00                             
            ALBEDO_FIELD_DIR(L, I)=0.0E+00                              
                                                                        
            IF (FLANDG(LIST(L)).LT.1.0) THEN                            
*D SWRAD3A.1130
     &            =SICE_ALB(LIST(L))*ICE_FRACTION(LIST(L))              
*D SWRAD3A.1134
     &            =SICE_ALB(LIST(L))*ICE_FRACTION(LIST(L))              
*D SWRAD3A.1139,ARE2F404.285  
            ENDIF                                                       
            IF (FLANDG(LIST(L)).GT.0.0) THEN                            
               IF ( L_MOSES_II ) THEN                                   
              IF ( L_CTILE ) THEN                                       
                ALBEDO_FIELD_DIFF(L,I) =                                
     &            (1.-FLANDG_G(L))*ALBEDO_FIELD_DIFF(L, I) +            
     &            FLANDG_G(L)*(WEIGHT_690NM(I)*LAND_ALBEDO(LIST(L),2)   
     &               + (1. - WEIGHT_690NM(I))*LAND_ALBEDO(LIST(L),4))   
                ALBEDO_FIELD_DIR(L,I) =                                 
     &            (1.-FLANDG_G(L))*ALBEDO_FIELD_DIR(L, I) +             
     &            FLANDG_G(L)*(WEIGHT_690NM(I)*LAND_ALBEDO(LIST(L),1)   
     &               + (1. - WEIGHT_690NM(I))*LAND_ALBEDO(LIST(L),3))   
               ELSE                                                     
*D ARE2F404.287,ARE2F404.288  
     &                            WEIGHT_690NM(I)*LAND_ALBEDO(LIST(L),2)
     &                   + (1. - WEIGHT_690NM(I))*LAND_ALBEDO(LIST(L),4)
*D ARE2F404.290,ARE2F404.291  
     &                            WEIGHT_690NM(I)*LAND_ALBEDO(LIST(L),1)
     &                   + (1. - WEIGHT_690NM(I))*LAND_ALBEDO(LIST(L),3)
               ENDIF                                                    
*D SWRAD3A.1140,SWRAD3A.1141 
! For non MOSES_II cannot have coastal tiling, therefore                
! must be completely land:                                              
               ALBEDO_FIELD_DIFF(L, I)=LAND_ALB(LIST(L))                
               ALBEDO_FIELD_DIR(L, I)=LAND_ALB(LIST(L))                 
*D SWRAD3A.1142,SWRAD3A.1143 
                                                                        
*DECLARE TRIF
*I TRIF.7     
     +,CROP(NPFT)                 ! 1 for crop type, 0 for non-crop.
     +,ORIENT(NPFT)               ! 1 for horizontal, 0 for spherical.
*I TRIF.11    
     +,ALNIR(NPFT)                ! Leaf reflection coefficient for
C                                 ! near infra-red.                  
     +,ALPAR(NPFT)                ! Leaf reflection coefficient for
C                                 ! PAR.                  
*I TRIF.53    
     +,OMNIR(NPFT)                ! Leaf scattering coefficient for
C                                 ! near infra-red.
*I TRIF.67    
      DATA CROP    /      0,     0,     1,     1,     0 /
      DATA ORIENT  /      0,     0,     0,     0,     0 /  
*I TRIF.68    
      DATA ALNIR   /   0.45,  0.35,  0.58,  0.58,  0.58 /               
      DATA ALPAR   /   0.10,  0.07,  0.10,  0.10,  0.10 /
*D ABX1F405.1732,TRIF.73   
      DATA DGL_DM  /    0.0,   0.0,   0.0,   0.0,   0.0 /               
      DATA DGL_DT  /    9.0,   9.0,   0.0,   0.0,   9.0 /               
*D ABX1F405.1733
      DATA FSMC_OF /   0.00,  0.00,  0.00,  0.00,  0.00 /               
*D TRIF.79
      DATA G_AREA  /  0.005, 0.004,  0.25,  0.25,  0.05 /               
*D ABX1F405.1734,TRIF.82   
      DATA G_LEAF_0/   0.25,  0.25,  0.25,  0.25,  0.25 /               
      DATA G_ROOT  /   0.25,  0.25,  0.25,  0.25,  0.25 /               
*D TRIF.85,TRIF.87   
      DATA LAI_MAX /   9.00,  9.00,  4.00,  4.00,  4.00 /               
      DATA LAI_MIN /   3.00,  3.00,  1.00,  1.00,  1.00 /               
      DATA NL0     /  0.040, 0.030, 0.060, 0.030, 0.030 /  
*I TRIF.90    
      DATA OMNIR   /   0.70,  0.45,  0.83,  0.83,  0.83 /             
*D TRIF.92,ABX1F405.1737 
      DATA SIGL    / 0.0375,0.1000,0.0250,0.0500,0.0500 /               
      DATA TLEAF_OF/ 273.15,243.15,258.15,258.15,243.15 /               
*DECLARE TRIFD2A
*D ABX1F405.1591
     &,                   FRAC_VS,FRAC_AGRIC,G_LEAF,NPP,RESP_S,RESP_W   
*D TRIFD2A.61,ABX1F405.1600 
     &,FRAC_AGRIC(LAND_FIELD)     ! IN Fraction of agriculture.    
*I TRIFD2A.89    
     &,FRAC_FLUX                  ! WORK PFT fraction to be used
C                                 !      in the calculation of
C                                 !      the gridbox mean fluxes.
*D ABX1F405.1618
     &,           C_VEG,FORW,FRAC_VS,FRAC_AGRIC,GAMMA,LAI_BAL,PC_S     
*D TRIFD2A.165
C type                                              
*D TRIFD2A.172,TRIFD2A.174  
          FRAC_FLUX=FRAC(L,N)-(1.0-FORW)*DFRAC(L,N)
          LIT_C(L,N) = NPP(L,N)-GAMMA/FRAC_FLUX*(C_VEG(L,N)*FRAC(L,N)
     &               -(C_VEG(L,N)-DCVEG(L,N))*(FRAC(L,N)-DFRAC(L,N)))
          LIT_C_T(L) = LIT_C_T(L)+FRAC_FLUX*LIT_C(L,N)                  
*DECLARE TYPPTRA
*D ABX1F404.24
*I TYPPTRA.36    
     &       JFRAC_LAND,             ! Land fraction in grid box
     &       JTSTAR_LAND,            ! Land surface temperature
     &       JTSTAR_SEA,             ! Sea surface temperature
     &       JTSTAR_SICE,            ! Sea-ice surface temperature
     &       JSICE_ALB,              ! Sea-ice albedo
     &       JLAND_ALB,              ! Mean land albedo
*D ABX1F404.41,ABX1F404.43   
     &       JCAN_WATER_TYP,         ! Canopy water content on tiles    
     &       JCATCH_TYP,             ! Canopy capacity on tiles         
     &       JINFIL_TYP,             ! Max infiltration rate on tiles   
     &       JRGRAIN_TYP,            ! Snow grain size on tiles         
     &       JSNODEP_TYP,            ! Snow depth on tiles              
*D ABX1F404.46,AJS1F401.24   
     &  JSNSOOT, JTSTAR_ANOM, 
     &  JFRAC_LAND, JTSTAR_LAND, JTSTAR_SEA, JTSTAR_SICE,
     &  JSICE_ALB, JLAND_ALB,
     &  JZH, JZ0, JLAND, JICE_FRACTION,                    
*D ABX1F404.50,ABX1F404.51   
     &  JRSP_S_ACC, JTSNOW, JCAN_WATER_TYP, JCATCH_TYP, JINFIL_TYP,     
     &  JRGRAIN_TYP, JSNODEP_TYP, JTSTAR_TYP, JZ0_TYP                   
*DECLARE TYPSIZE
*I TYPSIZE.25    
     &      ,NTILES               ! IN: No of land surface tiles
*I TYPSIZE.127   
     & NTILES,   
*DECLARE U_MODEL1
*D ARE2F404.530
     &    RADINCS ( (P_FIELDDA*(P_LEVELSDA+2+9)+511)/512*512*2 )
*DECLARE UPANCIL1
*I GDG0F401.1492  
     &                D1(JFRAC_LAND),                                   
     &                D1(JTSTAR_LAND),D1(JTSTAR_SEA),                   
     &                D1(JTSTAR_SICE),                                  
*DECLARE VEG1A
*D ABX3F405.29
     &,              LAND1,LAND_PTS,LAND_INDEX,NTILES,P_ROWS,ROW_LENGTH 
*D VEG1A.26
     &,              ATIMESTEP,SATCON                                
*D ABX1F405.1337
     &,              CATCH_T,INFIL_T,Z0_T                    
*I VEG1A.63    
     &,NTILES                ! IN Number of land-surface tiles.
*D VEG1A.81,VEG1A.82   
     & ATIMESTEP                    ! IN Atmospheric timestep (s).
     &,SATCON(LAND_FIELD)           ! IN Saturated hydraulic
!                                   !    conductivity (kg/m2/s).        
*D VEG1A.90,VEG1A.93   
     &,CATCH_T(LAND_FIELD,NTILES)   ! OUT Canopy capacity for tiles 
!                                   !     (kg/m2).  
     &,INFIL_T(LAND_FIELD,NTILES)   ! OUT Maximum surface infiltration 
!                                   !     rate for tiles (kg/m2/s).     
*D VEG1A.94,VEG1A.97   
     &,Z0_T(LAND_FIELD,NTILES)      ! OUT Roughness length for tiles (m)
*D VEG1A.114,VEG1A.115  
*D ABX1F405.1348,VEG1A.132  
      DO N=1,NTILES
*I VEG1A.134   
          INFIL_T(L,N)=0.0  
          Z0_T(L,N)=0.0 
*D VEG1A.193,VEG1A.203  
      CALL SPARM (LAND_FIELD,LAND1,LAND_PTS,NTILES,TILE_PTS,TILE_INDEX  
     &,           FRAC,HT,LAI,SATCON,CATCH_T,INFIL_T,Z0_T)
*DECLARE VEG2A
*D ABX3F405.74
     &,              LAND1,LAND_PTS,LAND_INDEX,NTILES,P_ROWS,ROW_LENGTH
*D VEG2A.29
     &,              ATIMESTEP,FRACA,SATCON                    
*D VEG2A.33
     &,              CATCH_T,INFIL_T,Z0_T                    
*I VEG2A.69    
     &,NTILES                ! IN Number of land-surface tiles.         
*D VEG2A.94,VEG2A.97   
     & ATIMESTEP                    ! IN Atmospheric timestep (s). 
     &,FRACA(LAND_FIELD)            ! IN Fraction of agriculture. 
     &,SATCON(LAND_FIELD)           ! IN Saturated hydraulic
C                                   !    conductivity (kg/m2/s).        
*D VEG2A.114,VEG2A.121  
     &,CATCH_T(LAND_FIELD,NTILES)   ! OUT Canopy capacity for tiles 
C                                   !     (kg/m2). 
     &,INFIL_T(LAND_FIELD,NTILES)   ! OUT Maximum surface infiltration 
!                                   !     rate for tiles (kg/m2/s).
     &,Z0_T(LAND_FIELD,NTILES)      ! OUT Roughness length for tiles (m)
*D VEG2A.150,ABX1F405.1445 
     &,RATIO                        ! WORK Ratio of fractional 
!                                   !      coverage before to that 
!                                   !      after TRIFFID. 
     &,DCS(LAND_FIELD)              ! WORK Change in soil carbon
!                                   !      (kg C/m2).
     &,FRAC_AGRIC(LAND_FIELD)       ! WORK Fraction of agriculture as
!                                   !      seen by TRIFFID.          
     &,FRAC_OLD(LAND_FIELD,NTYPE)   ! WORK Fractions of surface types 
!                                   !      before the call to TRIFFID.
     &,FRAC_VS(LAND_FIELD)          ! WORK Total fraction of gridbox    
!                                   !      covered by veg or soil.      
*D ABX1F405.1454,VEG2A.165  
*D VEG2A.170,VEG2A.173  

      LOGICAL
     & AGRIC                        ! .T. for TRIFFID to see agriculture
!                                   ! .F. for natural vegetation.
      PARAMETER (AGRIC=.TRUE.)
*D VEG2A.190
      DO N=1,NTILES                                                     
*I ABX1F405.1459  
          CATCH_T(L,N)=0.0
          INFIL_T(L,N)=0.0
*D VEG2A.196
*D VEG2A.198,VEG2A.207  
        RESP_S_DR(L)=0.0
*I ABX1F405.1462  
        FRAC_AGRIC(L) = 0.0
*D ABX3F405.96,ABX3F405.99   
*D ABX3F405.113
      CALL SWAPB_LAND(FRACA,LAND_FIELD,P_FIELD,                        
*D VEG2A.313
C Define the agricultural regions.
*I VEG2A.314   
        IF (AGRIC) THEN
*D VEG2A.316
            FRAC_AGRIC(L)=FRACA(L)
*I VEG2A.317   
        ENDIF
*D ABX1F405.1490,ABX1F405.1492 
C-----------------------------------------------------------------------
C Take copies of TRIFFID input variables for output as diagnostics.     
C-----------------------------------------------------------------------
*I ABX1F405.1497  
            FRAC_OLD(L,N)=FRAC(L,N)
*I ABX1F405.1501  
          DCS(L)=CS(L)
*D VEG2A.327
          FORW=0.0        
*D ABX1F405.1511
     &,                 FRAC_VS,FRAC_AGRIC,G_LEAF_DR,NPP_DR,RESP_S_DR   
*I ABX1F405.1522  

*D VEG2A.337
C Reset the accumulation fluxes to zero in equilibrium mode.
*I VEG2A.338   
        IF (L_TRIF_EQ) THEN

*D VEG2A.346
              RESP_W_AC(L,N)=0.0      
*I VEG2A.349   
C-----------------------------------------------------------------------
C Reset the accumulation fluxes to the TRIFFID-diagnosed corrections
C in dynamic mode. Such corrections will be typically associated with 
C the total depletion of a carbon reservoir or with the maintenance
C of the plant seed fraction.
C-----------------------------------------------------------------------
        ELSE

          DO L=LAND1,LAND1+LAND_PTS-1                                   
            DCS(L)=CS(L)-DCS(L)
            RESP_S_DR(L)=LIT_C_MN(L)-GAMMA*DCS(L)
            RESP_S_AC(L)=(RESP_S_DR(L)-RESP_S_DR_OUT(L))/GAM_TRIF   
          ENDDO                                                         
                                                                        
          DO N=1,NPFT                                                   
            DO J=1,TILE_PTS(N)                                          
              L=TILE_INDEX(J,N)                                         
              RATIO=FRAC_OLD(L,N)/FRAC(L,N)
              NPP_AC(L,N)=RATIO*(NPP_DR(L,N)-NPP_DR_OUT(L,N))/GAM_TRIF  
              RESP_W_AC(L,N)=RATIO*(RESP_W_DR(L,N)-RESP_W_DR_OUT(L,N))
     &                             /GAM_TRIF   
            ENDDO                                                       
          ENDDO                                                         

        ENDIF
                                                                        
C-----------------------------------------------------------------------
C Reset the accumulated leaf turnover rates to zero.
C-----------------------------------------------------------------------
*D VEG2A.372,VEG2A.382  
      CALL SPARM (LAND_FIELD,LAND1,LAND_PTS,NTILES,TILE_PTS,TILE_INDEX  
     &,           FRAC,HT,LAI,SATCON,CATCH_T,INFIL_T,Z0_T)           
*DECLARE VEG_CTL1
*I VEG_CTL1.93    
     &,I                   ! Loop counter for STASHWORK
*I VEG_CTL1.124   
! Initialise STASH workspace to zero to prevent STASH failure
      DO I=1,INT3
         STASHWORK(I)=0.0
       ENDDO

*D ABX3F405.8
     &            LAND_PTS,LAND_LIST,NTILES,P_ROWS,ROW_LENGTH,          
*D VEG_CTL1.130
     &            SECS_PER_STEPim(atmos_im),D1(JDISTURB),
     &            D1(JSAT_SOIL_COND), 
*D VEG_CTL1.134,VEG_CTL1.135  
     &            D1(JCANHT_PFT),D1(JCATCH_TYP),D1(JINFIL_TYP),     
     &            D1(JZ0_TYP),                                  
*D ABX1F405.1302
       CALL STASH(a_sm,a_im,19,STASHWORK,
*DECLARE VEG_IC1A
*D ABX3F405.13
     &,                 LAND1,LAND_PTS,LAND_INDEX,NTILES,P_ROWS
     &,                 ROW_LENGTH     
*D VEG_IC1A.28
     &,                 ATIMESTEP,FRAC_DISTURB,SATCON                 
*D VEG_IC1A.32
     &,                 CATCH_T,INFIL_T,Z0_T                 
*I VEG_IC1A.67    
     &,NTILES                ! IN Number of land-surface tiles.
*D VEG_IC1A.88,VEG_IC1A.89   
     & ATIMESTEP                    ! IN Atmospheric timestep (s).      
*I VEG_IC1A.91    
     &,SATCON(LAND_FIELD)           ! IN Saturated hydraulic
!                                   !    conductivity (kg/m2/s).        
*D VEG_IC1A.108,VEG_IC1A.111  
     &,CATCH_T(LAND_FIELD,NTILES)   ! OUT Canopy capacity for tiles 
!                                   !     (kg/m2).
     &,INFIL_T(LAND_FIELD,NTILES)   ! OUT Maximum surface infiltration 
!                                   !     rate for tiles (kg/m2/s).     
*D VEG_IC1A.112,VEG_IC1A.115  
     &,Z0_T(LAND_FIELD,NTILES)      ! OUT Roughness length for tiles (m)
*D ABX3F405.25
     &,        LAND1,LAND_PTS,LAND_INDEX,NTILES,P_ROWS,ROW_LENGTH       
*D VEG_IC1A.127
     &,        ATIMESTEP,SATCON                                       
*D ABX1F405.1334
     &,        CATCH_T,INFIL_T,Z0_T                          
*DECLARE VEG_IC2A
*D ABX3F405.58
     &,                 LAND1,LAND_PTS,LAND_INDEX,NTILES,P_ROWS
     &,                 ROW_LENGTH          
*D VEG_IC2A.29
     &,                 ATIMESTEP,FRAC_DISTURB,SATCON                 
*D VEG_IC2A.33
     &,                 CATCH_T,INFIL_T,Z0_T                 
*I VEG_IC2A.68    
     &,NTILES                       ! IN Number of land-surface tiles.  
*D VEG_IC2A.89,VEG_IC2A.90   
     & ATIMESTEP                    ! IN Atmospheric timestep (s).      
*I VEG_IC2A.92    
     &,SATCON(LAND_FIELD)           ! IN Saturated hydraulic
!                                   !    conductivity (kg/m2/s).       
*D VEG_IC2A.111,VEG_IC2A.112  
     &,CATCH_T(LAND_FIELD,NTILES)   ! OUT Canopy capacity for tiles 
C                                   !     (kg/m2). 
     &,INFIL_T(LAND_FIELD,NTILES)   ! OUT Maximum surface infiltration 
!                                   !     rate for tiles (kg/m2/s).     
*D VEG_IC2A.113,VEG_IC2A.116  
     &,Z0_T(LAND_FIELD,NTILES)      ! OUT Roughness length for tiles (m)
*D ABX3F405.70
     &,        LAND1,LAND_PTS,LAND_INDEX,NTILES,P_ROWS,ROW_LENGTH       
*D VEG_IC2A.131
     &,        ATIMESTEP,FRAC_DISTURB,SATCON                          
*D VEG_IC2A.135
     &,        CATCH_T,INFIL_T,Z0_T                          
*DECLARE VSHRZ7A
*D VSHRZ7A.33
     & N_ROWS,FIRST_ROW,ROW_LENGTH,FLANDG, 
*D VSHRZ7A.36
     & VSHR,VSHR_LAND,VSHR_SSI,Z1         
*D ABX1F405.836
     & FLANDG(P_FIELD)             ! IN Land fraction on all tiles
     &,AKH(2)                      ! IN Hybrid 'A' for layer 1.         
*I VSHRZ7A.72    
     &,VSHR_LAND(P_FIELD)          ! OUT Magnitude of surface-to-lowest 
!                                  !     atm level wind shear (m per s).
     &,VSHR_SSI(P_FIELD)           ! OUT Magnitude of surface-to-lowest 
!                                  !     atm level wind shear (m per s).
*I VSHRZ7A.149   
      IF(FLANDG(I).LT.1.0)THEN                                      
*D VSHRZ7A.153
        VSHR_SSI(I) = SQRT(VSHR2)                                   
      ELSE                                                            
        VSHR_SSI(I) = 0.0                                           
      ENDIF                                                           
C                                                                       
      IF(FLANDG(I).GT.0.0)THEN                                      
        VSHR2 = MAX (1.0E-6 , U_1_P(I)*U_1_P(I)                     
     &    + V_1_P(I)*V_1_P(I))                                      
        VSHR_LAND(I) = SQRT(VSHR2)                                    
      ELSE                                                            
        VSHR_LAND(I) = 0.0                                          
      ENDIF                                                           
C                                                                       
      VSHR(I)= FLANDG(I)*VSHR_LAND(I)                           
     &  + (1.0 - FLANDG(I))*VSHR_SSI(I)                           
*DECK CANCAP7A
*IF DEF,A03_7A,OR,DEF,A03_8A
C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1999, METEOROLOGICAL OFFICE, All Rights Reserved.
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
! Routine to calculate the heat capacity of a given PFT from its LAI
!
! Written by Richard Essery (September 1999)
C**********************************************************************
      SUBROUTINE CANCAP (LAND_FIELD,VEG_PTS,VEG_INDEX,FT
     &,                  HT,LAI,CANHC,VFRAC)

      IMPLICIT NONE

       INTEGER
     & LAND_FIELD                 ! IN Total number of land points.
     &,VEG_PTS                    ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_FIELD)      ! IN Index of vegetated points.

      INTEGER
     & FT                         ! IN Plant functional type.

      REAL
     & HT(LAND_FIELD)             ! IN Vegetation height (m).
     &,LAI(LAND_FIELD)            ! IN Leaf area index.
     &,CANHC(LAND_FIELD)          ! OUT Areal heat capacity of
!                                 !     vegetation canopy (J/K/m2).
     &,VFRAC(LAND_FIELD)          ! OUT Fractional canopy coverage.

      REAL
     & LAI_BAL(LAND_FIELD)        ! WORK Leaf area index in balanced
!                                 !      growth state.
     &,LEAF(LAND_FIELD)           ! WORK Leaf biomass (kg C/m2).
     &,WOOD(LAND_FIELD)           ! WORK Woody biomass (kg C/m2).

      INTEGER
     & J,L                        ! WORK Loop counters.

!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
      REAL
     & HLEAF                      ! Specific heat capacity of leaves
!                                 ! (J / K / kg Carbon).
     &,HWOOD                      ! Specific heat capacity of wood
!                                 ! (J / K / kg Carbon).
      PARAMETER ( HLEAF=5.7E4, HWOOD=1.1E4 )

*CALL NSTYPES
*CALL PFTPARM
*CALL TRIF
*CALL MOSES_OPT

      DO J=1,VEG_PTS
        L = VEG_INDEX(J)
        CANHC(L) = 0.
        VFRAC(L) = 0.
      ENDDO

      IF (CAN_MODEL .EQ. 2) THEN
!     Radiative canopy without heat capacity
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)
        CANHC(L) = 0.
        VFRAC(L) = 1. - EXP(-KEXT(FT)*LAI(L))
      ENDDO

      ELSEIF (CAN_MODEL .EQ. 3) THEN
!     Radiative canopy with heat capacity
        DO J=1,VEG_PTS
          L = VEG_INDEX(J)
          LAI_BAL(L) = ( A_WS(FT)*ETA_SL(FT)*HT(L) /
     &                   A_WL(FT) )**(1.0/(B_WL(FT)-1))
          LEAF(L) = SIGL(FT)*LAI_BAL(L)
          WOOD(L) = A_WL(FT)*(LAI_BAL(L)**B_WL(FT))
          CANHC(L) = HLEAF*LEAF(L) + HWOOD*WOOD(L)
          VFRAC(L) = 1. - EXP(-KEXT(FT)*LAI(L))
        ENDDO

      ENDIF

      RETURN
      END
*ENDIF
*DECK GAUSS7A
*IF DEF,A08_5A,OR,DEF,A08_7A
C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1999, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!    SUBROUTINE GAUSS--------------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE GAUSS (NLEVS,NPNTS,SOIL_PTS,SOIL_INDEX,A,B,C,D
     &,                 XMIN,XMAX,X)

      IMPLICIT NONE
!
! Description:
!     Solves a tridiagnonal matrix equation of the form:
!
!             A(n) X(n-1) + B(n) X(n) + C(n) X(n+1) = D(n)
!
!     by Gausian elimination, assuming boundary conditions:
!
!             A(1) = 0.0    at the top
!             C(N) = 0.0    at the bottom.
!                                                          (Cox, 2/99)
!
!
! Documentation :
!
! Current Code Owner : Peter Cox
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.6      2/99     New deck.  Peter Cox
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered:
! System Task:
!

! Global variables:

! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER
     & NLEVS                ! IN Number of levels.
     &,NPNTS                ! IN Number of gridpoints.
     &,SOIL_PTS             ! IN Number of soil points.


!   Array arguments with intent(IN) :
      INTEGER
     & SOIL_INDEX(NPNTS)    ! IN Array of soil points.

      REAL
     & A(NPNTS,NLEVS)       ! IN Matrix elements corresponding
!                           !    to the coefficients of X(n-1).
     &,B(NPNTS,NLEVS)       ! IN Matrix elements corresponding
!                           !    to the coefficients of X(n).
     &,C(NPNTS,NLEVS)       ! IN Matrix elements corresponding
!                           !    to the coefficients of X(n+1).
     &,D(NPNTS,NLEVS)       ! IN Matrix elements corresponding
!                           !    to the RHS of the equation.
     &,XMIN(NPNTS,NLEVS)    ! IN Minimum permitted value of X.
     &,XMAX(NPNTS,NLEVS)    ! IN Maximum permitted value of X.

!   Array arguments with intent(OUT) :
      REAL
     & X(NPNTS,NLEVS)       ! OUT Solution.

! Local scalars:
      INTEGER
     & I,J,N                ! WORK Loop counters.

! Local arrays:
      REAL
     & ADASH(NPNTS,NLEVS),BDASH(NPNTS,NLEVS),CDASH(NPNTS,NLEVS)
     &,DDASH(NPNTS,NLEVS)   ! WORK Transformed matrix elements

!-----------------------------------------------------------------------
! By default set the implicit increment to the explicit increment
! (for when denominators vanish).
!-----------------------------------------------------------------------
      DO N=1,NLEVS
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          X(I,N)=D(I,N)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Upward Sweep: eliminate "C" elements by replacing nth equation with:
!                  B'(n+1)*Eq(n)-C(n)*Eq'(n+1)
! where "'" denotes a previously tranformed equation. The resulting
! equations take the form:
!                A'(n) X(n-1) + B'(n) X(n) = D'(n)
! (NB. The bottom boundary condition implies that the NLEV equation does
!  not need transforming.)
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        ADASH(I,NLEVS)=A(I,NLEVS)
        BDASH(I,NLEVS)=B(I,NLEVS)
        DDASH(I,NLEVS)=D(I,NLEVS)
      ENDDO

      DO N=NLEVS-1,1,-1
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          ADASH(I,N)=BDASH(I,N+1)*A(I,N)
          BDASH(I,N)=BDASH(I,N+1)*B(I,N)-C(I,N)*ADASH(I,N+1)
          DDASH(I,N)=BDASH(I,N+1)*D(I,N)-C(I,N)*DDASH(I,N+1)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Top boundary condition: A(1) = 0.0 , allows X(1) to be diagnosed
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        IF (BDASH(I,1).NE.0.0) THEN
          X(I,1)=DDASH(I,1)/BDASH(I,1)
        ELSE
          WRITE(6,*)'Hydrology routine GAUSS : Layer 1 WARNING!'
        ENDIF
        X(I,1)=MAX(X(I,1),XMIN(I,1))
        X(I,1)=MIN(X(I,1),XMAX(I,1))
      ENDDO

!-----------------------------------------------------------------------
! Downward Sweep: calculate X(n) from X(n-1):
!                X(n) = (D'(n) - A'(n) X(n-1)) / B'(n)
!-----------------------------------------------------------------------
      DO N=2,NLEVS
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          IF (BDASH(I,N).NE.0.0) THEN
            X(I,N)=(DDASH(I,N)-ADASH(I,N)*X(I,N-1))/BDASH(I,N)
          ELSE
            WRITE(6,*)'Hydrology routine GAUSS : Layer ',N,' WARNING!'
          ENDIF
          X(I,N)=MAX(X(I,N),XMIN(I,N))
          X(I,N)=MIN(X(I,N),XMAX(I,N))
        ENDDO
      ENDDO

      RETURN
      END
*ENDIF
*DECK RESTILE7A
*IF DEF,A03_7A,OR,DEF,A03_8A
C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1999, METEOROLOGICAL OFFICE, All Rights Reserved.
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
C*LL  SUBROUTINE RES_TILE ----------------------------------------------
CLL
CLL  Purpose: Calculate gridbox-mean resistance factor used by TR_MIX
CLL           for 7A tiled land surface.
CLL
CLL  Model           Modification history:
CLL version  Date
CLL   4.6    11/99   New deck. R. Essery
CLL
C*----------------------------------------------------------------------
      SUBROUTINE RES_TILE (P_FIELD,LAND_FIELD,LAND_INDEX,
     &                     NTILES,TILE_INDEX,TILE_PTS,SOLUBLE,
     &                     ARESIST,ARESIST_TILE,CANOPY,CATCH,GS_TILE,
     &                     RESIST_B_TILE,SNOW_TILE,TILE_FRAC,RESB,RESS,
     &                     RES_FACTOR)

      IMPLICIT NONE

      INTEGER
     & P_FIELD               ! IN Total number of P-grid points.
     &,LAND_FIELD            ! IN Total number of land points.
     &,LAND_INDEX(P_FIELD)   ! IN Index of land points.
     &,NTILES                ! IN Number of land tiles.
     &,TILE_INDEX(LAND_FIELD,NTILES)
!                            ! IN Index of tile points.
     &,TILE_PTS(NTILES)      ! IN Number of tile points.

      LOGICAL
     & SOLUBLE               ! IN .TRUE. for soluble aerosols.

      REAL
     & ARESIST(P_FIELD)      ! IN GBM aerodynamic resistance (s/m).
     &,ARESIST_TILE(LAND_FIELD,NTILES)
!                            ! IN 1/(CD_STD*VSHR) on land tiles.
     &,CANOPY(LAND_FIELD,NTILES)
!                            ! IN Surface water on land tiles (kg/m2).
     &,CATCH(LAND_FIELD,NTILES)
!                            ! IN Surface capacity (max. surface water)
!                            !    of land tiles (kg/m2).
     &,GS_TILE(LAND_FIELD,NTILES)
!                            ! IN Surface conductance for land tiles.
     &,RESIST_B_TILE(LAND_FIELD,NTILES)
!                            ! IN (1/CH-1/CD_STD)/VSHR on land tiles.
     &,SNOW_TILE(LAND_FIELD,NTILES)
!                            ! IN Snow mass on land tiles (kg/m2).
     &,TILE_FRAC(LAND_FIELD,NTILES)
                             ! IN Tile fractions.
     &,RESB                  ! IN Rb(aerosol) / Rb(H2O).
     &,RESS                  ! IN Rs(aerosol) / Rs(H2O).
     &,RES_FACTOR(P_FIELD)   ! OUT Ra/(Ra+Rb+Rs) for dry deposition.

      REAL
     & DAMP_FACTOR(LAND_FIELD,NTILES)
!                            ! Canopy moistening factor
     &,RS_TILE(LAND_FIELD,NTILES)
!                            ! Surface reistance for land tiles.
     &,STR_RESIST_B          ! Rb for aerosol.
     &,STR_RESIST_S          ! Rs for aerosol.
     &,ASNOW                 ! Parameter for snow fraction calculation.
     &,COND_LIM              ! Low limit for canopy conductance.
     &,R_SNOW                ! Resistance to dry deposition over snow.
     &,SNOW_F                ! Snow cover fraction.
      PARAMETER (ASNOW=0.2, COND_LIM=1.0E-3, R_SNOW=1.0E3)

      INTEGER
     & I           ! Loop counter (horizontal field index).
     &,J           ! Loop counter (tile field index).
     &,L           ! Loop counter (land point field index).
     &,N           ! Loop counter (tile index).


      DO N=1,NTILES
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          DAMP_FACTOR(L,N) = 1.0
          IF (GS_TILE(L,N) .GT. COND_LIM) THEN
            RS_TILE(L,N) = 1. / GS_TILE(L,N)
          ELSE
            RS_TILE(L,N) = 1. / COND_LIM
          ENDIF
        ENDDO
      ENDDO

      IF (SOLUBLE) THEN
        DO N=1,NTILES
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            IF( (CATCH(L,N) .GT. 0.01) .AND.
     &          (CANOPY(L,N) .GT. 0.0) ) THEN
              IF( CANOPY(L,N) .LE. CATCH(L,N) ) THEN
                DAMP_FACTOR(L,N) = 1. - 0.66667*CANOPY(L,N)/CATCH(L,N)
              ELSE
                DAMP_FACTOR(L,N) = 0.33333
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      DO L=1,LAND_FIELD
        I = LAND_INDEX(L)
        RES_FACTOR(I) = 0.
      ENDDO

      DO N=1,NTILES
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          I = LAND_INDEX(L)
          STR_RESIST_B = RESB*RESIST_B_TILE(L,N)
          STR_RESIST_S = RESS*RS_TILE(L,N)*DAMP_FACTOR(L,N)
          IF (SNOW_TILE(L,N).GT.0.) THEN
            SNOW_F = 1. - EXP(-ASNOW*SNOW_TILE(L,N))
            STR_RESIST_S = 1. /
     &      (SNOW_F/R_SNOW + (1.-SNOW_F)/STR_RESIST_S)
          ENDIF
          RES_FACTOR(I) = RES_FACTOR(I) + ARESIST(I)*TILE_FRAC(L,N) /
     &                   (ARESIST_TILE(L,N)+STR_RESIST_B+STR_RESIST_S)
        ENDDO
      ENDDO

      RETURN
      END
*ENDIF
*DECK TILEALB                                                           
C *****************************COPYRIGHT******************************  
C (c) CROWN COPYRIGHT 2000, METEOROLOGICAL OFFICE, All Rights Reserved. 
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
! Routine to calculate albedos of land-surface tiles and gridbox-mean   
! albedo for MOSES II.                                                  
!                                                                       
! Richard Essery (March 2000)                                           
!                                                                       
C ******************************COPYRIGHT****************************** 
      SUBROUTINE TILE_ALBEDO (                                          
     & P_FIELD,LAND_FIELD,LAND1,LAND_PTS,LAND_INDEX,NTILES,TILE_PTS,    
     & TILE_INDEX,L_SNOW_ALBEDO,                                        
     & ALBSOIL,COSZ,FRAC,LAI_IN,RGRAIN,SNOW_TILE,TSTAR_TILE,Z0_TILE,    
     & ALB_TILE,LAND_ALBEDO                                             
     & )                                                                
                                                                        
      IMPLICIT NONE                                                     
                                                                        
*CALL NSTYPES                                                           
                                                                        
      INTEGER                                                           
     & P_FIELD                     ! IN Total number of grid points.    
     &,LAND_FIELD                  ! IN No. of land points.             
     &,LAND1                       ! IN First land point to be processed
     &,LAND_PTS                    ! IN No. of land pts to be processed.
     &,NTILES                      ! IN Number of surface tiles.        
     &,LAND_INDEX(P_FIELD)         ! IN Index of land points.           
     &,TILE_PTS(NTYPE)             ! IN Number of land points which     
!                                  !    include the nth surface type.   
     &,TILE_INDEX(LAND_FIELD,NTYPE)! IN Indices of land points which    
!                                  !    include the nth surface type.   
                                                                        
      LOGICAL                                                           
     & L_SNOW_ALBEDO               ! IN .TRUE. for spectral albedos     
!                                  !    and prognostic snow albedo.     
                                                                        
      REAL                                                              
     & ALBSOIL(LAND_FIELD)         ! IN Soil albedo.                    
     &,COSZ(P_FIELD)               ! IN Cosine of the zenith angle.     
     &,FRAC(LAND_FIELD,NTYPE)      ! IN Fractional cover of each        
!                                  !    surface type.                   
     &,LAI_IN(LAND_FIELD,NPFT)     ! IN Leaf area index.                
     &,RGRAIN(LAND_FIELD,NTILES)   ! IN Snow grain size on tiles        
!                                  !    (microns).                      
     &,SNOW_TILE(LAND_FIELD,NTILES)! IN Lying snow on tiles (kg/m2).    
     &,TSTAR_TILE(LAND_FIELD,NTILES)!IN Tile surface temperatures (K).  
     &,Z0_TILE(LAND_FIELD,NTILES)  ! IN Surface roughness on tiles (m). 
                                                                        
      REAL                                                              
     & ALB_TILE(LAND_FIELD,NTILES,4)!OUT Albedos for surface tiles.     
!                                  !     (*,*,1) - Direct beam visible  
!                                  !     (*,*,2) - Diffuse visible      
!                                  !     (*,*,3) - Direct beam near-IR  
!                                  !     (*,*,4) - Diffuse near-IR      
     &,LAND_ALBEDO(P_FIELD,4)      ! OUT GBM albedos.                   
                                                                        
      REAL                                                              
     & ALBSNC(LAND_FIELD,NTYPE)    ! Snow-covered albedo of surf types. 
     &,ALBSNF(LAND_FIELD,NTYPE)    ! Snow-free albedo of surf types.    
     &,ALB_TYPE(LAND_FIELD,NTYPE,4)! Albedos of surface types.          
     &,ALB_SNOW(LAND_FIELD,NTYPE,4)! Snow albedos.                      
     &,LAI(LAND_FIELD,NPFT)        ! Adjusted leaf area index.          
     &,SNOW(LAND_FIELD)            ! Copy of SNOW_TILE.                 
     &,TSTAR(LAND_FIELD)           ! Copy of TSTAR_TILE.                
     &,Z0(LAND_FIELD)              ! Copy of Z0_TILE.                   
     &,DSA                         ! Deep-snow albedo.                  
     &,FLIT                        ! Weighting factor for albedo.       
     &,FSNOW                       ! Weighting factor for albedo.       
                                                                        
      INTEGER                                                           
     & BAND,I,J,L,N                ! Loop counters                      
                                                                        
*CALL C_SOILH                                                           
*CALL C_0_DG_C                                                          
      REAL                                                              
     & DTLAND,KLAND,TCLAND,MASKD                                        
      PARAMETER( DTLAND = 2., KLAND = 0.3/DTLAND, TCLAND = TM-DTLAND,   
ccccc Tibet Snow mod ccccc
     &           MASKD = 0.1 )
cccccccccccccccccccccccccc
                                                                        
*CALL PFTPARM                                                           
*CALL NVEGPARM                                                          
                                                                        
      DO N=1,NTILES                                                     
        DO BAND=1,4                                                     
          DO L=1,LAND_FIELD                                             
            ALB_TILE(L,N,BAND) = 0.                                     
          ENDDO                                                         
        ENDDO                                                           
      ENDDO                                                             
      DO N=1,NTYPE                                                      
        DO BAND=1,4                                                     
          DO L=1,LAND_FIELD                                             
            ALB_TYPE(L,N,BAND) = 0.                                     
            ALB_SNOW(L,N,BAND) = 0.                                     
          ENDDO                                                         
        ENDDO                                                           
      ENDDO                                                             
                                                                        
! Impose minimum LAI for bare vegetation                                
      DO N=1,NPFT                                                       
        DO J=1,TILE_PTS(N)                                              
          L = TILE_INDEX(J,N)                                           
          LAI(L,N) = MAX( LAI_IN(L,N), 0.5 )                            
        ENDDO                                                           
      ENDDO                                                             
                                                                        
      IF (L_SNOW_ALBEDO) THEN                                           
!---------------------------------------------------------------------- 
! Spectral albedo scheme with prognostic snow albedo                    
!---------------------------------------------------------------------- 
                                                                        
! Set albedos of vegetated surface types                                
        CALL ALBPFT(P_FIELD,LAND_FIELD,LAND_INDEX,TILE_INDEX,TILE_PTS,  
     &              ALBSOIL,COSZ,LAI,ALB_TYPE)                          
                                                                        
! Set albedos of non-vegetated surface types                            
        DO BAND=1,4                                                     
          DO N=NPFT+1,NTYPE                                             
            DO J=1,TILE_PTS(N)                                          
              L = TILE_INDEX(J,N)                                       
              ALB_TYPE(L,N,BAND) = ALBSNF_NVG(N-NPFT)                   
              IF ( ALBSNF_NVG(N-NPFT).LT.0. )        ! Soil tile        
     &          ALB_TYPE(L,N,BAND) = ALBSOIL(L)                         
            ENDDO                                                       
          ENDDO                                                         
        ENDDO                                                           
                                                                        
! Calculate snow albedos                                                
        CALL ALBSNOW(P_FIELD,LAND_FIELD,LAND1,LAND_PTS,LAND_INDEX,      
     &               NTILES,TILE_INDEX,TILE_PTS,                        
     &               COSZ,RGRAIN,SNOW_TILE,ALB_SNOW)                    
                                                                        
! Adjust surface type albedos for snow cover                            
        DO L=LAND1,LAND1+LAND_PTS-1                                     
          SNOW(L) = SNOW_TILE(L,1)                                      
          Z0(L) = Z0_TILE(L,1)                                          
        ENDDO                                                           
        DO N=1,NTYPE                                                    
          IF (NTILES.NE.1) THEN                                         
            DO J=1,TILE_PTS(N)                                          
              L = TILE_INDEX(J,N)                                       
              SNOW(L) = SNOW_TILE(L,N)                                  
              Z0(L) = Z0_TILE(L,N)                                      
            ENDDO                                                       
          ENDIF                                                         
          DO J=1,TILE_PTS(N)                                            
            L = TILE_INDEX(J,N)                                         
            IF ( SNOW(L) .GT. 0.) THEN                                  
              FSNOW = SNOW(L) / (SNOW(L) + 10.*RHO_SNOW*Z0(L))          
              DO BAND=1,4                                               
                ALB_TYPE(L,N,BAND) = FSNOW*ALB_SNOW(L,N,BAND) +         
     &                               (1. - FSNOW)*ALB_TYPE(L,N,BAND)    
              ENDDO                                                     
            ENDIF                                                       
          ENDDO                                                         
        ENDDO                                                           
                                                                        
      ELSE                                                              
!---------------------------------------------------------------------- 
! Non-spectral albedo scheme with diagnosed snow albedo                 
!---------------------------------------------------------------------- 
                                                                        
! Set albedos of vegetated surface types                                
        DO N=1,NPFT                                                     
          DO J=1,TILE_PTS(N)                                            
            L = TILE_INDEX(J,N)                                         
            FLIT = 1.0 - EXP(-KEXT(N)*LAI(L,N))                         
            ALBSNC(L,N) = ALBSNC_MIN(N)*(1 - FLIT) + ALBSNC_MAX(N)*FLIT 
            ALBSNF(L,N) = ALBSOIL(L)*(1 - FLIT) + ALBSNF_MAX(N)*FLIT    
          ENDDO                                                         
        ENDDO                                                           
                                                                        
! Set albedos of non-vegetated surface types                            
        DO N=NPFT+1,NTYPE                                               
          DO J=1,TILE_PTS(N)                                            
            L = TILE_INDEX(J,N)                                         
            ALBSNC(L,N) = ALBSNC_NVG(N-NPFT)                            
            ALBSNF(L,N) = ALBSNF_NVG(N-NPFT)                            
            IF ( ALBSNF_NVG(N-NPFT).LT.0. ) ALBSNF(L,N) = ALBSOIL(L)    
          ENDDO                                                         
        ENDDO                                                           
                                                                        
! Adjust surface type albedos for snow cover                            
        DO L=LAND1,LAND1+LAND_PTS-1                                     
          TSTAR(L) = TSTAR_TILE(L,1)                                    
          SNOW(L) = SNOW_TILE(L,1)                                      
        ENDDO                                                           
        DO N=1,NTYPE                                                    
          IF (NTILES.NE.1) THEN                                         
            DO J=1,TILE_PTS(N)                                          
              L = TILE_INDEX(J,N)                                       
              TSTAR(L) = TSTAR_TILE(L,N)                                
              SNOW(L) = SNOW_TILE(L,N)                                  
            ENDDO                                                       
          ENDIF                                                         
          DO J=1,TILE_PTS(N)                                            
            L = TILE_INDEX(J,N)                                         
            IF ( TSTAR(L) .LT. TCLAND ) THEN                            
              DSA = ALBSNC(L,N)                                         
            ELSEIF ( TSTAR(L) .LT. TM ) THEN                            
              DSA = ALBSNC(L,N) + KLAND*(ALBSNF(L,N) - ALBSNC(L,N))     
     &                                 *(TSTAR(L) - TCLAND)             
            ELSE                                                        
              DSA = ALBSNC(L,N) + KLAND*(ALBSNF(L,N) - ALBSNC(L,N))     
     &                                 *(TM - TCLAND)                   
            ENDIF                                                       
            ALB_TYPE(L,N,1) = ALBSNF(L,N) + (DSA - ALBSNF(L,N)) *       
     &                                    ( 1. - EXP(-MASKD*SNOW(L)) )  
          ENDDO                                                         
        ENDDO                                                           
                                                                        
! Copy albedo to all bands                                              
        DO BAND=2,4                                                     
          DO N=1,NTYPE                                                  
            DO J=1,TILE_PTS(N)                                          
              L = TILE_INDEX(J,N)                                       
              ALB_TYPE(L,N,BAND) = ALB_TYPE(L,N,1)                      
            ENDDO                                                       
          ENDDO                                                         
        ENDDO                                                           
                                                                        
      ENDIF       ! Spectral or non-spectral albedo schemes             
                                                                        
!---------------------------------------------------------------------- 
! Calculate GBM surface albedo                                          
!---------------------------------------------------------------------- 
                                                                        
      DO BAND=1,4                                                       
        DO I=1,P_FIELD                                                  
          LAND_ALBEDO(I,BAND) = 0.                                      
        ENDDO                                                           
        DO N=1,NTYPE                                                    
          DO J=1,TILE_PTS(N)                                            
            L = TILE_INDEX(J,N)                                         
            I = LAND_INDEX(L)                                           
            LAND_ALBEDO(I,BAND) = LAND_ALBEDO(I,BAND) +                 
     &                            FRAC(L,N)*ALB_TYPE(L,N,BAND)          
          ENDDO                                                         
        ENDDO                                                           
      ENDDO                                                             
                                                                        
!---------------------------------------------------------------------- 
! Copy albedos as required for aggregate or distinct tiles              
!---------------------------------------------------------------------- 
                                                                        
      IF (NTILES.EQ.1) THEN                                             
        DO BAND=1,4                                                     
          DO L=LAND1,LAND1+LAND_PTS-1                                   
            I = LAND_INDEX(L)                                           
            ALB_TILE(L,1,BAND) = LAND_ALBEDO(I,BAND)                    
          ENDDO                                                         
        ENDDO                                                           
      ELSE                                                              
        DO BAND=1,4                                                     
          DO N=1,NTYPE                                                  
            DO J=1,TILE_PTS(N)                                          
              L = TILE_INDEX(J,N)                                       
              ALB_TILE(L,N,BAND) = ALB_TYPE(L,N,BAND)                   
            ENDDO                                                       
          ENDDO                                                         
        ENDDO                                                           
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
*DECK ALBPFT
C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 2000, METEOROLOGICAL OFFICE, All Rights Reserved.
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
! Routine to calculate the spectral albedo of the land surface using
! the two stream approach of Sellers, 1995.
!
! Written by Peter Cox (Mar 1999)
! Adapted for MOSES 2.2 by Richard Essery (March 2000)
!
C**********************************************************************
      SUBROUTINE ALBPFT (P_FIELD,LAND_FIELD,LAND_INDEX,TILE_INDEX,
     &                   TILE_PTS,ALBSOIL,COSZ,LAI,ALB_TYPE)

      IMPLICIT NONE

*CALL NSTYPES

      INTEGER
     & P_FIELD                     ! IN Total number of grid points.
     &,LAND_FIELD                  ! IN Number of land points.
     &,LAND_INDEX(P_FIELD)         ! IN Index of land points.
     &,TILE_PTS(NTYPE)             ! IN Number of land points which
!                                  !    include the surface type.
     &,TILE_INDEX(LAND_FIELD,NTYPE)! IN Indices of land points which
!                                  !    include the surface type.

      REAL
     & ALBSOIL(LAND_FIELD)         ! IN Soil albedo.
     &,COSZ(P_FIELD)               ! IN Cosine of the zenith angle.
     &,LAI(LAND_FIELD,NPFT)        ! IN Leaf area index.

      REAL
     & ALB_TYPE(LAND_FIELD,NTYPE,4)!OUT Albedos for surface types.
!                                  !     (*,*,1) - Direct beam visible
!                                  !     (*,*,2) - Diffuse visible
!                                  !     (*,*,3) - Direct beam near-IR
!                                  !     (*,*,4) - Diffuse near-IR

      REAL
     & ALBUDIF(LAND_FIELD,2)      ! Diffuse albedo of the underlying
!                                 ! surface.
     &,ALBUDIR(LAND_FIELD,2)      ! Direct albedo of the underlying
!                                 ! surface.
     &,ALPL(2)                    ! Leaf reflection coefficient.
     &,BETADIR                    ! Upscatter parameter for direct beam.
     &,BETADIF                    ! Upscatter parameter for diffuse beam
     &,COSZM                      ! Mean value of COSZ.
     &,K                          ! Optical depth per unit leaf area.
     &,G                          ! Relative projected leaf area in
                                  ! direction cosz.
     &,OM(2)                      ! Leaf scattering coefficient.
     &,SALB                       ! Single scattering albedo.
     &,SQCOST                     ! Cosine squared of the mean leaf angl
                                  ! to the horizontal.
     &,TAUL(2)                    ! Leaf transmission coefficient.
     &,B,C,CA,D,F,H,U1            ! Miscellaneous variables from
     &,P1,P2,P3,P4,D1             ! Sellers (1985).
     &,H1,H2,H3,H7,H8             !
     &,S1,S2,SIG                  !

      INTEGER
     & BAND,I,J,L,N               ! Loop counters.

*CALL TRIF

      DO L=1,LAND_FIELD
        ALBUDIF(L,1) = ALBSOIL(L)
        ALBUDIF(L,2) = ALBSOIL(L)
        ALBUDIR(L,1) = ALBSOIL(L)
        ALBUDIR(L,2) = ALBSOIL(L)
      ENDDO

      DO N=1,NPFT

        OM(1) = OMEGA(N)
        OM(2) = OMNIR(N)
        ALPL(1) = ALPAR(N)
        ALPL(2) = ALNIR(N)

        DO BAND=1,2  ! Visible and near-IR bands
          TAUL(BAND) = OM(BAND) - ALPL(BAND)
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            I = LAND_INDEX(L)
            IF (ORIENT(N).EQ.0) THEN
              SQCOST = 1./3.
              G = 0.5
              COSZM = 1.0
              SALB = 0.5*OM(BAND)
              IF (COSZ(I).GT.0.01)
     &          SALB = 0.5*OM(BAND) *
     &                 ( 1. - COSZ(I)*LOG((COSZ(I)+1.)/COSZ(I)) )
            ELSEIF (ORIENT(N).EQ.1) THEN
              SQCOST = 1.
              G = COSZ(I)
              COSZM = 1.
              SALB = OM(BAND)/4.
            ENDIF
            K = G / 0.01
            IF (COSZ(I).GT.0.01) K = G / COSZ(I)
            BETADIR = (1. + COSZM*K)/(OM(BAND)*COSZM*K)*SALB
            C = 0.5*( ALPL(BAND) + TAUL(BAND) +
     &               (ALPL(BAND) - TAUL(BAND))*SQCOST )
            BETADIF = C / OM(BAND)
            B = 1. - (1. - BETADIF)*OM(BAND)
            D = OM(BAND)*COSZM*K*BETADIR
            F = OM(BAND)*COSZM*K*(1. - BETADIR)
            H = SQRT(B*B - C*C) / COSZM
            SIG = (COSZM*K)**2 + C*C - B*B
            U1 = B - C/ALBUDIF(L,BAND)
            CA = C*ALBUDIR(L,BAND)/ALBUDIF(L,BAND)
            S1 = EXP(-H*LAI(L,N))
            S2 = EXP(-K*LAI(L,N))
            P1 = B + COSZM*H
            P2 = B - COSZM*H
            P3 = B + COSZM*K
            P4 = B - COSZM*K
            D1 = P1*(U1 - COSZM*H)/S1 - P2*(U1 + COSZM*H)*S1
            H1 = -D*P4 - C*F
            H2 = ( (D - P3*H1/SIG) * (U1 - COSZM*H) / S1 -
     &             (D - CA - (U1 + COSZM*K)*H1/SIG)*P2*S2 ) / D1
            H3 = - ( (D - P3*H1/SIG) * (U1 + COSZM*H)*S1 -
     &               (D - CA - (U1 + COSZM*K)*H1/SIG)*P1*S2 ) / D1
            H7 = (C/D1)*(U1 - COSZM*H) / S1
            H8 = - (C/D1)*(U1 + COSZM*H) * S1
            ALB_TYPE(L,N,2*BAND-1) = H1/SIG + H2 + H3   ! Direct beam
            ALB_TYPE(L,N,2*BAND) = H7 + H8              ! Diffuse
          ENDDO
        ENDDO

      ENDDO

      RETURN
      END
*DECK ALBSNOW
C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 2000, METEOROLOGICAL OFFICE, All Rights Reserved.
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
! Routine to calculate spectral snow albedos for MOSES II, based on
! the Marshall (1989) parametrization of the Wiscombe and Warren (1980)
! model. Influence of contaminants in snow has not been included - see
! UM vn4.5 deck FTSA1A.
!
! Richard Essery (March 2000)
!
C *********************************************************************
      SUBROUTINE ALBSNOW (P_FIELD,LAND_FIELD,LAND1,LAND_PTS,LAND_INDEX,
     &                    NTILES,TILE_INDEX,TILE_PTS,
     &                    COSZ,RGRAIN,SNOW_TILE,ALB_SNOW)

      IMPLICIT NONE

*CALL NSTYPES

      INTEGER
     & P_FIELD                     ! IN Total number of grid points.
     &,LAND_FIELD                  ! IN Number of land points.
     &,LAND1                       ! IN First land point to be processed
     &,LAND_PTS                    ! IN No of land pts to be processed.
     &,NTILES                      ! IN Number of surface tiles.
     &,LAND_INDEX(P_FIELD)         ! IN Index of land points.
     &,TILE_PTS(NTYPE)             ! IN Number tile points.
     &,TILE_INDEX(LAND_FIELD,NTYPE)! IN Index of tile points.

      REAL
     & COSZ(P_FIELD)               ! IN Zenith cosine.
     &,RGRAIN(LAND_FIELD,NTILES)   ! IN Snow grain size (microns).
     &,SNOW_TILE(LAND_FIELD,NTILES)! IN Lying snow (kg/m2).

      REAL
     & ALB_SNOW(LAND_FIELD,NTYPE,4)! OUT Snow albedo.
!                                  !     (*,*,1) - Direct beam visible
!                                  !     (*,*,2) - Diffuse visible
!                                  !     (*,*,3) - Direct beam near-IR
!                                  !     (*,*,4) - Diffuse near-IR

      REAL
     & REFF                        ! Zenith effective grain size.
     &,R0                          ! Grain size for fresh snow (microns)
      PARAMETER ( R0 = 50. )

      INTEGER
     & BAND,I,J,L,N                ! Loop counters.

      REAL
     & AMAX(2)                     ! Maximum albedo for fresh snow
!                 VIS     NIR
      DATA AMAX / 0.98,   0.7   /

      DO N=1,NTILES
        DO L=LAND1,LAND1+LAND_PTS-1
          I = LAND_INDEX(L)
          IF (SNOW_TILE(L,N) .GT. 0.) THEN
            REFF = RGRAIN(L,N) * ( 1. + 0.77*(COSZ(I)-0.65) )**2
            ALB_SNOW(L,N,1) = AMAX(1) - 0.002*(SQRT(REFF) - SQRT(R0))
            ALB_SNOW(L,N,2) = AMAX(1) -
     &                        0.002*(SQRT(RGRAIN(L,N)) - SQRT(R0))
            ALB_SNOW(L,N,3) = AMAX(2) - 0.09*ALOG(REFF/R0)
            ALB_SNOW(L,N,4) = AMAX(2) - 0.09*ALOG(RGRAIN(L,N)/R0)
          ENDIF
        ENDDO
      ENDDO

      IF (NTILES .EQ. 1) THEN
        DO BAND=1,4
          DO N=2,NTYPE
            DO J=1,TILE_PTS(N)
              L = TILE_INDEX(J,N)
              ALB_SNOW(L,N,BAND) = ALB_SNOW(L,1,BAND)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END
*DECK SICEHT7A                                                          
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
C *****************************COPYRIGHT******************************  
C (c) CROWN COPYRIGHT 2000, METEOROLOGICAL OFFICE, All Rights Reserved. 
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
!!!  SUBROUTINE SICE_HTF----------------------------------------------- 
!!!                                                                     
!!!  Purpose: Updates sea-ice surface layer temperature.                
!!!                                                                     
!!!  Model            Modification history                              
!!! version  date                                                       
!!!                                                                     
!!! New deck for MOSES II. R. Essery, 4/4/00.                           
!!!                                                                     
!!!  Note: At present the formulation is so simple as to make this      
!!!        routine fairly trivial; but in future the formulation may    
!!!        be revised so as to make a separate routine more obviously   
!!!        worthwhile.                                                  
!!!                                                                     
!!!  Programming standard: Unified Model Documentation Paper No.4       
!!!                        version no.2, dated 18/1/90.                 
!!!                                                                     
!!!  System component covered: P241                                     
!!!                                                                     
!!!  Documentation: ??                                                  
!!!                                                                     
                                                                        
! Arguments:---------------------------------------------------------   
      SUBROUTINE SICE_HTF (                                             
     & POINTS,P_FIELD,P1,FLANDG,SIMLT,                               
     & DI,ICE_FRACTION,SURF_HT_FLUX_SICE,
     & TSTAR_SEA,TSTAR_SICE,TIMESTEP,                 
     & TI,SICE_MLT_HTF,SEA_ICE_HTF,                                     
     & LTIMER)                                                          
                                                                        
      IMPLICIT NONE                                                     
                                                                        
      LOGICAL LTIMER                                                    
                                                                        
      INTEGER                                                           
     & POINTS               ! IN No of gridpoints to be processed.      
     &,P_FIELD              ! IN Total Number of points on p-grid       
     &,P1                   ! IN First point of p grid to be processed  
                                                                        
      LOGICAL                                                           
     & LAND_MASK(P_FIELD)   ! IN Land mask (T for land, F for sea).     
     &,SIMLT                ! IN STASH flag for sea-ice melting ht flux.
                                                                        
      REAL                                                              
     & FLANDG(P_FIELD)      ! IN Land fraction.                 
     &,DI(P_FIELD)          ! IN "Equivalent thickness" of sea-ice (m). 
     &,ICE_FRACTION(P_FIELD)! IN Fraction of gridbox covered by sea-ice.
     &,SURF_HT_FLUX_SICE(P_FIELD)                               
!                           ! IN Net downward heat flux at      
!                           !    sea-ice surface W/m2           
     &,TSTAR_SICE(P_FIELD)  ! IN Sea-ice surface                
!                           !    temperature (K).               
     &,TSTAR_SEA(P_FIELD)   ! IN Open sea surface               
!                           !    temperature (K).               
     &,TIMESTEP             ! IN Timestep (s).                          
                                                                        
      REAL                                                              
     & TI(P_FIELD)          ! INOUT Sea-ice surface layer temperature(K)
!                           !       Set to TSTAR for unfrozen sea,      
!                           !       missing data for land.              
     &,SICE_MLT_HTF(P_FIELD)! INOUT Heat flux due to melting of sea-ice 
!                           !       (W/m2).                             
     &,SEA_ICE_HTF(P_FIELD) ! OUT Heat flux through sea-ice (W per sq m,
!                                 positive downwards).                  
                                                                        
      EXTERNAL TIMER                                                    
                                                                        
!  Common and local physical constants.                                 
*CALL C_0_DG_C                                                          
*CALL C_KAPPAI                                                          
*CALL C_SICEHC                                                          
                                                                        
      REAL                                                              
     & ASURF(P_FIELD)       ! Reciprocal areal heat capacity of         
!                             sea-ice surface layer (Km2/J).            
      INTEGER I             ! Loop Counter; horizontal field index.     
                                                                        
      IF (LTIMER) THEN                                                  
        CALL TIMER('SICEHTF ',3)                                        
      ENDIF                                                             
                                                                        
      DO I=P1,P1+POINTS-1                                               
        IF (FLANDG(I).EQ.1.0) THEN                                    
          SEA_ICE_HTF(I)=0.0                                            
          TI(I) = 1.0E30                                                
        ELSE IF (ICE_FRACTION(I).LE.0.0) THEN                           
          SEA_ICE_HTF(I)=0.0                                            
          TI(I) = TFS                                                 
          TSTAR_SICE(I) = TFS                                         
        ELSE                                                            
          ASURF(I) = AI / ICE_FRACTION(I)                               
          SEA_ICE_HTF(I) = ICE_FRACTION(I)*KAPPAI*(TI(I) - TFS)/DI(I)   
          TI(I) = TI(I) + ASURF(I)*TIMESTEP*                            
     &                    (SURF_HT_FLUX_SICE(I) - SEA_ICE_HTF(I))       
          if (ti(I).lt.150) then
            write(6,*)"sicehtf: capping stupid cold sea-ice temp",ti(i)
            write(6,*)"sicehtf: frac, equiv. thick"
            write(6,*)"sicehtf:",ICE_FRACTION(I),DI(I)
            ti(i)=150.
          endif
          IF ( TI(I) .GT. TM ) THEN                                     
            IF (SIMLT) SICE_MLT_HTF(I) = SICE_MLT_HTF(I) +              
     &                                  (TI(I) - TM)/(ASURF(I)*TIMESTEP)
            TI(I) = TM                                                  
          ENDIF                                                         
        ENDIF                                                           
      ENDDO                                                             
                                                                        
      IF (LTIMER) THEN                                                  
        CALL TIMER('SICEHTF ',4)                                        
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
*ENDIF                                                                  
*DECK SOILEV7A
*IF DEF,A03_7A,OR,DEF,A03_8A
C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 2000, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!**********************************************************************
! Subroutine to adjust canopy conductance and soil moisture extraction
! for soil evaporation beneath vegetation.
!
! Richard Essery. 6/4/00.
!
!**********************************************************************

      SUBROUTINE SOIL_EVAP (NPNTS,NSHYD,TILE_PTS,TILE_INDEX,
     &                      GSOIL,LAI,GS,WT_EXT)

      IMPLICIT NONE

      INTEGER
     & NPNTS                ! IN Number of gridpoints.
     &,NSHYD                ! IN Number of soil moisture layers.
     &,TILE_PTS             ! IN Number of points containing the
!                           !    given surface type.
     &,TILE_INDEX(NPNTS)    ! IN Indices on the land grid of the
!                           !    points containing the given
!                           !    surface type.

       REAL
     & GSOIL(NPNTS)         ! IN Soil surface conductance (m/s).
     &,LAI(NPNTS)           ! IN Leaf area index.

      REAL
     & GS(NPNTS)            ! INOUT Surface conductance (m/s).
     &,WT_EXT(NPNTS,NSHYD)  ! INOUT Fraction of evapotranspiration
!                           !       extracted from each soil layer.

       REAL
     & FSOIL(NPNTS)         ! Fraction of ground below canopy
!                           ! contributing to evaporation.

       INTEGER
     & J,K,L                ! Loop indices

        DO J=1,TILE_PTS
          L=TILE_INDEX(J)
          FSOIL(L) = EXP(-0.5*LAI(L))
        ENDDO

        DO K=2,NSHYD
          DO J=1,TILE_PTS
            L=TILE_INDEX(J)
            WT_EXT(L,K) = GS(L)*WT_EXT(L,K)/(GS(L) + FSOIL(L)*GSOIL(L))
          ENDDO
        ENDDO

        DO J=1,TILE_PTS
          L=TILE_INDEX(J)
          WT_EXT(L,1) = (GS(L)*WT_EXT(L,1) + FSOIL(L)*GSOIL(L))
     &                   / (GS(L) + FSOIL(L)*GSOIL(L))
          GS(L) = GS(L) + FSOIL(L)*GSOIL(L)
        ENDDO

      RETURN
      END
*ENDIF
*DECK BDY_EXPL7A
*IF DEF,A03_7A
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
!!!  SUBROUTINE BDY_EXPL-----------------------------------------------
!!!
!!!  Purpose: Calculate explicit boundary layer fluxes of heat, moisture
!!!           and momentum. Also calculates boundary layer exchange
!!!           coefficients required for implicit update of boundary
!!!           layer fluxes
!!!
!!!
!!! F.Hewer     <- programmer of some or all of previous code or changes
!!! C.Wilson    <- programmer of some or all of previous code or changes
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   4.3  7/2/97     New deck. S Jackson
!!!   4.4 25/6/97     Modified for MOSES II tile model. R Essery
!!!   4.4 25/6/97     Move grid definitions up to BL_INTCT.  R.A.Betts
!!!  4.5    Jul. 98  Kill the IBM specific lines. (JCThil)
!!!   4.5  7/5/98     Set TSTAR, SNOW_SURF_HTF and SOIL_SURF_HTF to 0
!!!                   at all land points, to avoid problems of
!!!                   non-initialised data.  R.A.Betts
!!!   4.5 21/5/98     Add optional error check for negative surface
!!!                   temperature.  R.A.Betts
!!!
!!!  Programming standard: Unified Model Documentation Paper No 4,
!!!                        Version ?, dated ?.
!!!
!!!  System component covered: P24.
!!!
!!!  Project task:
!!!
!!!  Documentation: UMDP 24.
!!!
!!!---------------------------------------------------------------------

!    Arguments :-
      SUBROUTINE BDY_EXPL (

! IN values defining field dimensions and subset to be processed :
     & P_FIELD,U_FIELD,ROW_LENGTH,
     & N_P_ROWS,N_U_ROWS,P_POINTS,P1,U_POINTS,U1,

! IN values defining vertical grid of model atmosphere :
     & BL_LEVELS,P_LEVELS,AK,BK,AKH,BKH,DELTA_AK,DELTA_BK,
     & EXNER,

! IN sea/sea-ice data :
     & U_0,V_0,

! IN cloud data :
     & CF,QCF,QCL,CCA,CCB,CCT,

! IN everything not covered so far :
     & PSTAR,RAD_HR,RADHR_DIM1,
     & FB_SURF,U_S,T1_SD,Q1_SD,TV1_SD,
     & H_BLEND_OROG,Z0M_EFF,
     & TIMESTEP,L_BL_LSPICE,L_MOM,

! INOUT data :
     & Q,T,U,V,ZH,

! OUT Diagnostic not requiring STASH flags :
     & QW,TL,FQW,FTL,
     & RHOKH,RHOKM_UV,
     & TAUX,TAUY,ZHT,
     & BL_TYPE_1,BL_TYPE_2,BL_TYPE_3,BL_TYPE_4,BL_TYPE_5,BL_TYPE_6,

! OUT data required for tracer mixing :
     & NRML,

! OUT data required for 4D-VAR :
     & RHO_KM,

! OUT data required elsewhere in UM system :
     & DTRDZ,RDZ,RDZUV,
     & DU,DV,CT_CTQ,DQW,DTL,CQ_CM,

! LOGICAL LTIMER
     & LTIMER
     & )

      IMPLICIT NONE

!  Inputs :-

! (a) Defining horizontal grid and subset thereof to be processed.
!    Checked for consistency.

      INTEGER
     & P_FIELD                     ! IN No. of P-points in whole grid
!                                  !    (for dimensioning only).
     &,U_FIELD                     ! IN No. of UV-points in whole grid.
     &,ROW_LENGTH                  ! IN No. of points in one row.
     &,N_P_ROWS   ! IN No of P-rows being processed.
     &,N_U_ROWS   ! IN No of UV-rows being processed.
     &,P_POINTS   ! IN No of P-points being processed.
     &,P1         ! IN First P-point to be processed.
     &,U_POINTS   ! IN No of UV-points being processed.
     &,U1         ! IN First UV-point to be processed.

! (b) Defining vertical grid of model atmosphere.

      INTEGER
     & BL_LEVELS                   ! IN Max. no. of "boundary" levels
!                                  !    allowed. Assumed <= 30 for dim-
!                                  !    ensioning GAMMA in common deck
!                                  !    C_GAMMA used in SF_EXCH and KMKH
     &,P_LEVELS                    ! IN Total no. of vertical levels in
!                                  !    the model atmosphere.
      REAL
     & AK(P_LEVELS)                ! IN Hybrid 'A' for all levels.
     &,BK(P_LEVELS)                ! IN Hybrid 'B' for all levels.
     &,AKH(P_LEVELS+1)             ! IN Hybrid 'A' for layer interfaces.
     &,BKH(P_LEVELS+1)             ! IN Hybrid 'B' for layer interfaces.
     &,DELTA_AK(P_LEVELS)          ! IN Difference of hybrid 'A' across
!                                  !    layers (K-1/2 to K+1/2).
!                                  !    NB: Upper minus lower.
     &,DELTA_BK(P_LEVELS)          ! IN Difference of hybrid 'B' across
!                                  !     layers (K-1/2 to K+1/2).
!                                  !     NB: Upper minus lower.
     &,EXNER(P_FIELD,BL_LEVELS+1)  ! IN Exner function.  EXNER(,K) is
!                                  !    value for LOWER BOUNDARY of
!                                  !    level K.

! (d) Sea/sea-ice data.

      REAL
     & U_0(U_FIELD)                ! IN W'ly component of surface
!                                  !    current (m/s).
     &,V_0(U_FIELD)                ! IN S'ly component of surface
!                                  !    current (m/s).

! (e) Cloud data.

      REAL
     & CF(P_FIELD,BL_LEVELS)       ! IN Cloud fraction (decimal).
     &,QCF(P_FIELD,BL_LEVELS)      ! IN Cloud ice (kg per kg air)
     &,QCL(P_FIELD,BL_LEVELS)      ! IN Cloud liquid water (kg
!                                  !    per kg air).
     &,CCA(P_FIELD)                ! IN Convective Cloud Amount
!                                  !    (decimal)

      INTEGER
     & CCB(P_FIELD)                ! IN Convective Cloud Base
     &,CCT(P_FIELD)                ! IN Convective Cloud Top

! (f) Atmospheric + any other data not covered so far, incl control.

      REAL
     & PSTAR(P_FIELD)              ! IN Surface pressure (Pascals).
     &,TIMESTEP                    ! IN Timestep (seconds).
     &,H_BLEND_OROG(P_FIELD)       ! IN Blending height used as part of
!                                  !    effective roughness scheme
     &,Z0M_EFF(P_FIELD)            ! IN Effective grid-box roughness
!                                  !    length for momentum
     &,Q(P_FIELD,BL_LEVELS)        ! IN Specific humidity ( kg/kg air).
     &,T(P_FIELD,BL_LEVELS)        ! IN Atmospheric temperature (K).
     &,U(U_FIELD,BL_LEVELS)        ! IN W'ly wind component (m/s)
     &,V(U_FIELD,BL_LEVELS)        ! IN S'ly wind component (m/s)
!                                  !       length for momentum (m).

      LOGICAL
     & LTIMER                      ! IN Logical switch for TIMER diags
     &,L_BL_LSPICE                 ! IN Use if 3A large scale precip
     &,L_MOM                       ! IN Switch for convective momentum
!                                  !    transport.

!  Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-
!
      REAL
     & FQW(P_FIELD,BL_LEVELS)      ! OUT Moisture flux between layers
!                                  !     (kg per square metre per sec).
!                                  !     FQW(,1) is total water flux
!                                  !     from surface, 'E'.
     &,FTL(P_FIELD,BL_LEVELS)      ! OUT FTL(,K) contains net turbulent
!                                  !     sensible heat flux into layer K
!                                  !     from below; so FTL(,1) is the
!                                  !     surface sensible heat, H.(W/m2)
     &,RHOKH(P_FIELD,BL_LEVELS)    ! OUT Exchange coeffs for moisture.
     &,RHOKM_UV(U_FIELD,BL_LEVELS) ! OUT Exchange coefficients for
!                                  !     momentum (on UV-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,TAUX(U_FIELD,BL_LEVELS)     ! OUT W'ly component of surface wind
!                                  !     stress (N/sq m). (On UV-grid
!                                  !     with first and last rows
!                                  !     undefined or, at present,
!                                  !     set to missing data
     &,TAUY(U_FIELD,BL_LEVELS)     ! OUT S'ly component of surface wind
!                                  !     stress (N/sq m).  On UV-grid;
!                                  !     comments as per TAUX.
     &,RHO_KM(P_FIELD,2:BL_LEVELS) ! OUT Air density * turbulent mixing
!                                  !     coefficient for momentum before
!                                  !     interpolation.
     &,QW(P_FIELD,BL_LEVELS)       ! OUT Total water content, but
!                                  !     replaced by specific humidity
!                                  !     in LS_CLD.
     &,DTRDZ(P_FIELD,BL_LEVELS)    ! OUT -g.dt/dp for model layers.
     &,RDZ(P_FIELD,BL_LEVELS)      ! OUT RDZ(,1) is the reciprocal of
!                                  !     the height of level 1, i.e. of
!                                  !     the middle of layer 1.  For
!                                  !     K > 1, RDZ(,K) is the
!                                  !     reciprocal of the vertical
!                                  !     distance from level K-1 to
!                                  !     level K.
     &,RDZUV(U_FIELD,BL_LEVELS)    ! OUT RDZ (K > 1) on UV-grid.
!                                  !     Comments as per RHOKM (RDZUV).
     &,TL(P_FIELD,BL_LEVELS)       ! OUT Ice/liquid water temperature,
!                                  !     but replaced by T in LS_CLD.

       REAL
     & CT_CTQ(P_FIELD,BL_LEVELS)   ! OUT Coefficient in T and q
!                                        tri-diagonal implicit matrix
     &,CQ_CM(U_FIELD,BL_LEVELS)    ! OUT Coefficient in U and V
!                                        tri-diagonal implicit matrix
     &,DQW(P_FIELD,BL_LEVELS)      ! OUT BL increment to q field
     &,DTL(P_FIELD,BL_LEVELS)      ! OUT BL increment to T field
     &,DU(U_FIELD,BL_LEVELS)       ! OUT BL increment to u wind field
     &,DV(U_FIELD,BL_LEVELS)       ! OUT BL increment to v wind field

      INTEGER
     & NRML(P_FIELD)               ! OUT Number of model layers in the
!                                  !     Rapidly Mixing Layer; set to
!                                  !     zero in SF_EXCH for MOSES II.

!-2 Genuinely output, needed by other atmospheric routines :-

      REAL
     & ZH(P_FIELD)                 ! OUT Height above surface of top of
!                                  !     boundary layer (metres).

      INTEGER
     & RADHR_DIM1            ! DUMMY Used in 6A boundary layer scheme

      REAL
     & RAD_HR(RADHR_DIM1,BL_LEVELS)
!                            ! DUMMY Used in 6A boundary layer scheme
     &,FB_SURF(P_FIELD)      ! DUMMY Used in 6A boundary layer scheme
     &,U_S(P_FIELD)          ! DUMMY Used in 6A boundary layer scheme
     &,T1_SD(P_FIELD)        ! DUMMY Used in 6A boundary layer scheme
     &,Q1_SD(P_FIELD)        ! DUMMY Used in 6A boundary layer scheme
     &,TV1_SD(P_FIELD)       ! DUMMY Used in 6A boundary layer scheme
     &,ZHT(P_FIELD)          ! DUMMY Used in 6A boundary layer scheme
     &,BL_TYPE_1(P_FIELD)    ! DUMMY Used in 6A boundary layer scheme
     &,BL_TYPE_2(P_FIELD)    ! DUMMY Used in 6A boundary layer scheme
     &,BL_TYPE_3(P_FIELD)    ! DUMMY Used in 6A boundary layer scheme
     &,BL_TYPE_4(P_FIELD)    ! DUMMY Used in 6A boundary layer scheme
     &,BL_TYPE_5(P_FIELD)    ! DUMMY Used in 6A boundary layer scheme
     &,BL_TYPE_6(P_FIELD)    ! DUMMY Used in 6A boundary layer scheme


!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL Z,HEAT_CON,SF_EXCH,BOUY_TQ,BTQ_INT,
     & KMKH,EX_FLUX_TQ,EX_FLUX_UV,IM_CAL_TQ,SICE_HTF,SF_EVAP,SF_MELT,
     & IM_CAL_UV,SCREEN_TQ
      EXTERNAL TIMER
*IF -DEF,SCMA
      EXTERNAL UV_TO_P,P_TO_UV
*ENDIF

!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

*CALL C_R_CP
*CALL C_G
*CALL C_LHEAT
*CALL C_GAMMA
*CALL SOIL_THICK
*IF DEF,MPP
! MPP Common block
*CALL PARVARS
*ENDIF

! Derived local parameters.

      REAL LCRCP,LS,LSRCP

      PARAMETER (
     & LCRCP=LC/CP           ! Evaporation-to-dT conversion factor.
     &,LS=LF+LC              ! Latent heat of sublimation.
     &,LSRCP=LS/CP           ! Sublimation-to-dT conversion factor.
     &  )

!-----------------------------------------------------------------------

!  Workspace :-

      REAL
     & BF(P_FIELD,BL_LEVELS)    ! A buoyancy parameter (beta F tilde)
     &,BQ(P_FIELD,BL_LEVELS)    ! A buoyancy parameter (beta q tilde).
     &,BT(P_FIELD,BL_LEVELS)    ! A buoyancy parameter (beta T tilde).
     &,BT_CLD(P_FIELD,BL_LEVELS)
!                               ! A buoyancy parameter for cloudy air
!                               ! on p,T,q-levels (full levels).
     &,BQ_CLD(P_FIELD,BL_LEVELS)! A buoyancy parameter for cloudy air
!                               ! on p,T,q-levels (full levels).
     &,A_QS(P_FIELD,BL_LEVELS)  ! Saturated lapse rate factor
!                               ! on p,T,q-levels (full levels).
     &,A_DQSDT(P_FIELD,BL_LEVELS)
!                               ! Saturated lapse rate factor
!                               ! on p,T,q-levels (full levels).
     &,DQSDT(P_FIELD,BL_LEVELS) ! Derivative of q_SAT w.r.t. T
     &,DELTAP(P_FIELD,BL_LEVELS)! Difference in pressure between levels
     &,DELTAP_UV(P_FIELD,BL_LEVELS)
!                               ! Difference in pressure between levels
!                               ! on UV points
     &,DTRDZ_UV(U_FIELD,BL_LEVELS)
!                               ! -g.dt/dp for model wind layers.
     &,DZL(P_FIELD,BL_LEVELS)   ! DZL(,K) is depth in m of layer
!                               ! K, i.e. distance from boundary
!                               ! K-1/2 to boundary K+1/2.
     &,DU_NT(U_FIELD,BL_LEVELS) ! non-turbulent inc. to u wind field
     &,DV_NT(U_FIELD,BL_LEVELS) ! non-turbulent inc. to v wind field
     &,DTL_NT(P_FIELD,BL_LEVELS)! non-turbulent inc. to TL field
     &,DQW_NT(P_FIELD,BL_LEVELS)! non-turbulent inc. to QW field
     &,P(P_FIELD,BL_LEVELS)     ! Pressure at model levels
     &,RHO(P_FIELD,BL_LEVELS)   ! Density of model layer
     &,RHOKM(P_FIELD,BL_LEVELS) ! Exchange coefficients for
!                               ! momentum on P-grid
     &,TV(P_FIELD,BL_LEVELS)    ! Virtual temp
     &,U_P(P_FIELD,BL_LEVELS)   ! U on P-grid.
     &,V_P(P_FIELD,BL_LEVELS)   ! V on P-grid.
     &,ZLB(P_FIELD,0:BL_LEVELS) ! ZLB(,K) is the height of the
!                               ! upper boundary of layer K
!                               ! ( = 0.0 for "K=0").


!  Local scalars :-

      REAL
     & WK         ! LOCAL 0.5 * DZL(I,K) * RDZ(I,K)
     &,WKM1       ! LOCAL 0.5 * DZL(I,K-1) * RDZ(I,K)

      INTEGER
     & I          ! LOCAL Loop counter (horizontal field index).
     &,K          ! LOCAL Loop counter (vertical level index).

      IF (LTIMER) THEN
        CALL TIMER('BDYLAYR ',3)
      ENDIF

!-----------------------------------------------------------------------
!! 1.1 Initialise ZLB(,0) (to zero, of course, this being the height
!!     of the surface above the surface).
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1
        ZLB(I,0)=0.0
      ENDDO

!-----------------------------------------------------------------------
!! 1.2 Calculate layer depths and heights, and construct wind fields on
!!     P-grid.  This involves calling subroutines Z and UV_TO_P.
!!     Virtual temperature is also calculated, as a by-product.
!-----------------------------------------------------------------------
!  NB RDZ  TEMPORARILY used to return DELTA_Z_LOWER, the lower half
!     layer thickness

      DO K=1,BL_LEVELS
        CALL Z(P_POINTS,EXNER(P1,K),EXNER(P1,K+1),PSTAR(P1),
     &    AKH(K),BKH(K),Q(P1,K),QCF(P1,K),
     &    QCL(P1,K),T(P1,K),ZLB(P1,K-1),TV(P1,K),
     &    ZLB(P1,K),DZL(P1,K),RDZ(P1,K),LTIMER)

*IF -DEF,SCMA
        CALL UV_TO_P(U(U1,K),U_P(P1,K),
     &               U_POINTS,P_POINTS,ROW_LENGTH,N_U_ROWS)
        CALL UV_TO_P(V(U1,K),V_P(P1,K),
     &               U_POINTS,P_POINTS,ROW_LENGTH,N_U_ROWS)

*IF DEF,MPP
! DZL can contain incorrect data in halos, so call SWAPBOUNDS.
      CALL SWAPBOUNDS(DZL(P1,1),ROW_LENGTH,N_U_ROWS,1,0,BL_LEVELS)

*ENDIF
! du_nt 'borrowed to store dzl on uv grid
        CALL P_TO_UV (DZL(P1,K),DU_NT(U1+ROW_LENGTH,K),
     &     P_POINTS,U_POINTS,ROW_LENGTH,N_P_ROWS)

*ELSE
      DO I = U1, U1-1+U_POINTS
        U_P(i,K) = U(i,K)
        V_P(i,K) = V(i,K)
      END DO
*ENDIF
      ENDDO

! set pressure array.
      DO K=1,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          P(I,K) = AK(K) + BK(K)*PSTAR(I)

! These will be used in new dynamics scheme - currently unused
          DTL_NT(I,K)=0.0
          DQW_NT(I,K)=0.0

        ENDDO

      ENDDO  ! end of loop over bl_levels

      DO K=BL_LEVELS,2,-1

        DO I=P1,P1+P_POINTS-1
          RDZ(I,K)=1.0/(RDZ(I,K)+(DZL(I,K-1)-RDZ(I,K-1)))
          DELTAP(I,K)=DELTA_AK(K) + PSTAR(I)*DELTA_BK(K)

          DTRDZ(I,K) = -G * TIMESTEP/ DELTAP(I,K)
!     &                  (DELTA_AK(K) + PSTAR(I)*DELTA_BK(K))
        ENDDO
      ENDDO

      DO I=P1,P1+P_POINTS-1
        RDZ(I,1)=1.0/RDZ(I,1)
        DELTAP(I,1)=DELTA_AK(1) + PSTAR(I)*DELTA_BK(1)
        DTRDZ(I,1) = -G * TIMESTEP/DELTAP(I,1)
!     &                  (DELTA_AK(1) + PSTAR(I)*DELTA_BK(1))
      ENDDO

      DO K=1,BL_LEVELS

! Calculate RDZUV here

        IF(K.GE.2)THEN
*IF -DEF,SCMA

          DO I=U1+ROW_LENGTH,U1-ROW_LENGTH+U_POINTS-1
            RDZUV(I,K) = 2.0 / ( DU_NT(I,K) + DU_NT(I,K-1) )
          ENDDO

!-----------------------------------------------------------------------
! 1.3 Set first and last rows to "missing data indicator"
!-----------------------------------------------------------------------

*IF DEF,MPP
      IF (attop) THEN
*ENDIF
        DO I=U1,U1+ROW_LENGTH-1
          RDZUV(I,K) = 1.0E30
        ENDDO
*IF DEF,MPP
      ENDIF

      IF (atbase) THEN
*ENDIF
        DO I= U1+(N_U_ROWS-1)*ROW_LENGTH, U1 + N_U_ROWS*ROW_LENGTH-1
          RDZUV(I,K) = 1.0E30
        ENDDO
*IF DEF,MPP
      ENDIF
*ENDIF

*ELSE
      DO I = U1, U1-1+U_POINTS
        RDZUV(i,K) = 2.0 / ( DZL(i,K) + DZL(i,K-1) )
      ENDDO
*ENDIF
        ENDIF   ! K .ge. 2

! Calculate DTRDZ_UV here.

*IF -DEF,SCMA
!        CALL P_TO_UV (DTRDZ(P1,K),DTRDZ_UV(U1+ROW_LENGTH,K),
!     &     P_POINTS,U_POINTS,ROW_LENGTH,N_P_ROWS)

        CALL P_TO_UV (DELTAP(P1,K),DELTAP_UV(U1+ROW_LENGTH,K),
     &     P_POINTS,U_POINTS,ROW_LENGTH,N_P_ROWS)

        DO I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1
          DTRDZ_UV(I,K) = -G * TIMESTEP / DELTAP_UV(I,K)
        ENDDO

*ELSE
      DO I = P1, P1-1+P_POINTS
        DTRDZ_UV(i,K) = DTRDZ(i,K)
      ENDDO
*ENDIF

      ENDDO ! loop over bl_levels

! "borrowed" du_nt reset to zero
! Non turbulent increments for new dynamics scheme (currently not used)
        DO K=1,BL_LEVELS
          DO I=1,U_FIELD
            DU_NT(I,K) =0.0
            DV_NT(I,K) =0.0
          ENDDO
        ENDDO

!-----------------------------------------------------------------------
!! Calculate total water content, QW and Liquid water temperature, TL
!-----------------------------------------------------------------------
      DO K=1,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          QW(I,K) = Q(I,K) + QCL(I,K) + QCF(I,K)              ! P243.10
          TL(I,K) = T(I,K) - LCRCP*QCL(I,K) - LSRCP*QCF(I,K)  ! P243.9
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!! Calculate buoyancy parameters BT and BQ.
!-----------------------------------------------------------------------
      CALL BOUY_TQ (
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,P,CF,T,TL,Q,QCF,QCL
     &,BT,BQ,BF,BT_CLD,BQ_CLD,A_QS,A_DQSDT,DQSDT
     &,L_BL_LSPICE,LTIMER
     & )

!-----------------------------------------------------------------------
!! 5.  Turbulent exchange coefficients and "explicit" fluxes between
!!     model layers in the boundary layer (P243b, routine KMKH).
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!!      Interpolate BT and BQ to interface between layers.
!-----------------------------------------------------------------------

      CALL BTQ_INT (
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,BQ,BT,BF,DZL,RDZ,QW,QCF,TL
     &,L_BL_LSPICE,LTIMER
     &  )

!-----------------------------------------------------------------------
!! 5.3  Calculate the diffusion coefficients Km and Kh.
!-----------------------------------------------------------------------

! Repeat of KMKH calculation, could be passed in from KMKH.

      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          WKM1 = 0.5 * DZL(I,K-1) * RDZ(I,K)
          WK = 0.5 * DZL(I,K) * RDZ(I,K)

! Calculate rho at K-1/2, from P243.111 :-
          RHO(I,K) =
     &     ( AKH(K) + BKH(K)*PSTAR(I) )    ! Pressure at K-1/2, P243.112
     &     /                               ! divided by ...
     &     ( R *                           ! R times ...
     &     ( TV(I,K-1)*WK + TV(I,K)*WKM1 ) ! TV at K-1/2, from P243.113
     &     )
        ENDDO
      ENDDO

      CALL KMKH (
     & P_FIELD,P1,P_POINTS,BL_LEVELS,
     & TIMESTEP,P,CCA,BT,BQ,BF,CF,DZL,DTRDZ,
     & RDZ,U_P,V_P,FTL,FQW,
     & RHO,Z0M_EFF,ZLB(1,0),H_BLEND_OROG,
     & QW,QCF,RHOKM,RHO_KM(1,2),RHOKH,TL,ZH,
     & CCB,CCT,L_MOM,
     & NRML,L_BL_LSPICE,LTIMER
     & )

!-----------------------------------------------------------------------
!! 5.4 Interpolate RHOKM's and CDR10M to uv points ready for the
!!     calculation of the explcit fluxes TAU_X and TAU_Y at levels
!!     above the surface.
!-----------------------------------------------------------------------

*IF DEF,MPP
! RHOKM(*,1) contains duff data in halos. The P_TO_UV can interpolate
! this into the real data, so first we must update east/west halos

      CALL SWAPBOUNDS(RHOKM(1,2),ROW_LENGTH,
     &                U_FIELD/ROW_LENGTH,1,1,BL_LEVELS-1)
*ENDIF

      DO K=2,BL_LEVELS

*IF -DEF,SCMA
        CALL P_TO_UV (RHOKM(P1,K),RHOKM_UV(U1+ROW_LENGTH,K),
     &     P_POINTS,U_POINTS,ROW_LENGTH,N_P_ROWS)
*IF DEF,MPP
      IF (attop) THEN
*ENDIF
        DO I=U1,U1+ROW_LENGTH-1
          RHOKM_UV(I,K) = 1.0E30
        ENDDO
*IF DEF,MPP
      ENDIF

      IF (atbase) THEN
*ENDIF
        DO I= U1+(N_U_ROWS-1)*ROW_LENGTH, U1+N_U_ROWS*ROW_LENGTH-1
          RHOKM_UV(I,K) = 1.0E30
        ENDDO
*IF DEF,MPP
      ENDIF
*ENDIF

*ELSE
      DO I = P1, P1-1+P_POINTS
        RHOKM_UV(i,K) = RHOKM(i,K)
      ENDDO
*ENDIF
      ENDDO ! loop over bl_levels

!-----------------------------------------------------------------------
!! 5.5 Calculation of explicit fluxes of T,Q
!-----------------------------------------------------------------------

      CALL EX_FLUX_TQ (
     &  P_POINTS,P_FIELD,P1,BL_LEVELS
     &, TL,QW,RDZ,FTL,FQW,RHOKH
     &, LTIMER
     &  )

!-----------------------------------------------------------------------
!! 5.6 Calculation of explicit fluxes of U and V.
!-----------------------------------------------------------------------

      CALL EX_FLUX_UV ( ! For U
     &  U_POINTS,U_FIELD,ROW_LENGTH,BL_LEVELS,U1
     &, U,U_0,RDZUV(1,2),RHOKM_UV,TAUX
     &, LTIMER
     &  )

      CALL EX_FLUX_UV ( ! For V
     &  U_POINTS,U_FIELD,ROW_LENGTH,BL_LEVELS,U1
     &, V,V_0,RDZUV(1,2),RHOKM_UV,TAUY
     &, LTIMER
     &  )

*IF -DEF,SCMA
!-----------------------------------------------------------------------
!! Set first and last rows to "missing data indicator"
!-----------------------------------------------------------------------
      DO K=2,BL_LEVELS
*IF DEF,MPP
      IF (attop) THEN
*ENDIF
        DO I=U1,U1+ROW_LENGTH-1
          TAUX(I,K)=1.E30
          TAUY(I,K)=1.E30
        ENDDO
*IF DEF,MPP
      ENDIF

      IF (atbase) THEN
*ENDIF
        DO I= U1 + (N_U_ROWS-1)*ROW_LENGTH, U1 + N_U_ROWS*ROW_LENGTH -1
          TAUX(I,K)=1.E30
          TAUY(I,K)=1.E30
        ENDDO
*IF DEF,MPP
      ENDIF
*ENDIF
      ENDDO
*ENDIF

!-----------------------------------------------------------------------
!! 6.  "Implicit" calculation of increments for TL and QW
!-----------------------------------------------------------------------

      CALL IM_BL_PT1 (
     & P_FIELD,P1,U_FIELD,U1
     &,P_POINTS,U_POINTS,ROW_LENGTH,BL_LEVELS
     &,DTRDZ,DTRDZ_UV,RHOKH(1,2),RHOKM_UV(1,2)
     &,RDZ,RDZUV(1,2),GAMMA
     &,DQW_NT,DTL_NT,DU_NT,DV_NT
     &,FQW,FTL,TAUX,TAUY
     &,CT_CTQ,DQW,DTL,CQ_CM,DU,DV
     &,LTIMER
     &)

      IF (LTIMER) THEN
        CALL TIMER('BDYLAYR ',4)
      ENDIF

      RETURN
      END
*ENDIF
*DECK BDY_EXPL8A
*IF DEF,A03_8A
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
!!!  SUBROUTINE BDY_EXPL-----------------------------------------------
!!!
!!!  Purpose: Calculate explicit boundary layer fluxes of heat, moisture
!!!           and momentum. Also calculates boundary layer exchange
!!!           coefficients required for implicit update of boundary
!!!           layer fluxes
!!!
!!!
!!! F.Hewer     <- programmer of some or all of previous code or changes
!!! C.Wilson    <- programmer of some or all of previous code or changes
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   4.3  7/2/97     New deck. S Jackson
!!!   4.4 25/6/97     Modified for MOSES II tile model. R Essery
!!!   4.4 25/6/97     Move grid definitions up to BL_INTCT.  R.A.Betts
!!!  4.5    Jul. 98  Kill the IBM specific lines. (JCThil)
!!!   4.5  7/5/98     Set TSTAR, SNOW_SURF_HTF and SOIL_SURF_HTF to 0
!!!                   at all land points, to avoid problems of
!!!                   non-initialised data.  R.A.Betts
!!!   4.5 21/5/98     Add optional error check for negative surface
!!!                   temperature.  R.A.Betts
!!!
!!!  Programming standard: Unified Model Documentation Paper No 4,
!!!                        Version ?, dated ?.
!!!
!!!  System component covered: P24.
!!!
!!!  Project task:
!!!
!!!  Documentation: UMDP 24.
!!!
!!!---------------------------------------------------------------------

!    Arguments :-
      SUBROUTINE BDY_EXPL (

! IN values defining field dimensions and subset to be processed :
     & P_FIELD,U_FIELD,ROW_LENGTH,
     & N_P_ROWS,N_U_ROWS,P_POINTS,P1,U_POINTS,U1,

! IN values defining vertical grid of model atmosphere :
     & BL_LEVELS,P_LEVELS,AK,BK,AKH,BKH,DELTA_AK,DELTA_BK,
     & EXNER,

! IN sea/sea-ice data :
     & U_0,V_0,

! IN cloud data :
     & CF,QCF,QCL,CCA,CCB,CCT,

! IN everything not covered so far :
     & PSTAR,RAD_HR,RADHR_DIM1,
     & FB_SURF,U_S,T1_SD,Q1_SD,TV1_SD,
     & H_BLEND_OROG,Z0M_EFF_GB,
     & TIMESTEP,L_BL_LSPICE,L_MOM,

! INOUT data :
     & Q,T,U,V,ZH,

! OUT Diagnostic not requiring STASH flags :
     & QW,TL,FQW,FTL,
     & RHOKH,RHOKM_UV,
     & TAUX,TAUY,ZHT,
     & BL_TYPE_1,BL_TYPE_2,BL_TYPE_3,BL_TYPE_4,BL_TYPE_5,BL_TYPE_6,

! OUT data required for tracer mixing :
     & NRML,

! OUT data required for 4D-VAR :
     & RHO_KM,

! OUT data required elsewhere in UM system :
     & DTRDZ,RDZ,RDZUV,
     & DU,DV,CT_CTQ,DQW,DTL,CQ_CM,

! LOGICAL LTIMER
     & LTIMER
     & )

      IMPLICIT NONE

!  Inputs :-

! (a) Defining horizontal grid and subset thereof to be processed.

      INTEGER
     & P_FIELD                     ! IN No. of P-points in whole grid
!                                     (for dimensioning only).
     &,P1                          ! IN First point to be processed in
!                                       P-grid.
     &,RADHR_DIM1                  ! IN Dimension of Radiative heating
!                                  !    rate (P_FIELD but used for
!                                  !    dynamic allocation)
     &,U_FIELD                     ! IN No. of UV-points in whole grid.
!                                     (Checked for consistency with
!                                     P_FIELD and P_ROWS; there must
!                                     be 1 less UV than P row.)
     &,U1                          ! IN First point to be processed in
!                                       U_V-grid.
     &,P_POINTS                    ! IN Number of P-grid points to be
!                                       processed.
     &,U_POINTS                    ! IN Number of U_V-grid points.
     &,N_P_ROWS                    ! IN No of P-rows being processed.
     &,N_U_ROWS                    ! IN No of UV-rows being processed.
     &,ROW_LENGTH                  ! IN No. of points in one row.
!                                     (Checked for consistency with
!                                     P_FIELD and N_ROWS.)

! (b) Defining vertical grid of model atmosphere.

      INTEGER
     & BL_LEVELS                   ! IN Max. no. of "boundary" levels
!                                     allowed.Assumed <= 30 for dim-
!                                     sioning of GAMMA in common deck
!                                     C_GAMMA used in SF_EXCH and KMKH
     &,P_LEVELS                    ! IN Total no. of vertical levels in
!                                       the model atmosphere.
      REAL
     & AK(P_LEVELS)                ! IN Hybrid 'A' for all levels.
     &,BK(P_LEVELS)                ! IN Hybrid 'B' for all levels.
     &,AKH(P_LEVELS+1)             ! IN Hybrid 'A' for layer interfaces.
     &,BKH(P_LEVELS+1)             ! IN Hybrid 'B' for layer interfaces.
     &,DELTA_AK(P_LEVELS)          ! IN Difference of hybrid 'A' across
!                                     layers (K-1/2 to K+1/2).
!                                     NB: Upper minus lower.
     &,DELTA_BK(P_LEVELS)          ! IN Difference of hybrid 'B' across
!                                     layers (K-1/2 to K+1/2).
!                                     NB: Upper minus lower.
     &,EXNER(P_FIELD,BL_LEVELS+1)  ! IN Exner function.  EXNER(,K) is
!                                     value for LOWER BOUNDARY of
!                                     level K.

! (c) Soil/vegetation/land surface parameters (mostly constant).

      LOGICAL
     & L_BL_LSPICE                 ! IN True if 3A large-scale ppn
!                                       scheme is used.
     &,L_MOM                       ! IN Switch for convective momentum
!                                  !    transport.

! (d) Sea/sea-ice data.

      REAL
     & U_0(U_FIELD)                ! IN W'ly component of surface
!                                     current (m/s).
     &,V_0(U_FIELD)                ! IN S'ly component of surface
!                                     current (m/s).

! (e) Cloud data.

      REAL
     & CF(P_FIELD,BL_LEVELS)       ! IN Cloud fraction (decimal).
     &,QCF(P_FIELD,BL_LEVELS)      ! IN Cloud ice (kg per kg air)
     &,QCL(P_FIELD,BL_LEVELS)      ! IN Cloud liquid water (kg
!                                     per kg air).
     &,CCA(P_FIELD)                ! IN Convective Cloud Amount
!                                     (decimal)

      INTEGER
     & CCB(P_FIELD)                ! IN Convective Cloud Base
     &,CCT(P_FIELD)                ! IN Convective Cloud Top

! (f) Atmospheric + any other data not covered so far, incl control.

      REAL
     & PSTAR(P_FIELD)              ! IN Surface pressure (Pascals).
     &,RAD_HR(RADHR_DIM1,BL_LEVELS)! IN Radiative heating rate (K/s).
     &,FB_SURF(P_FIELD)            ! IN Surface flux buoyancy over
!                                  !  density (m^2/s^3)
     &,U_S(P_FIELD)                ! IN Surface friction velocity (m/s)
     &,T1_SD(P_FIELD)              ! IN Standard deviation of turbulent
!                                     fluctuations of layer 1
!                                     temperature; for use in initiating
!                                     convection.
     &,Q1_SD(P_FIELD)              ! IN Standard deviation of turbulent
!                                     fluctuations of layer 1 humidity;
!                                     for use in initiating convection.
     &,TV1_SD(P_FIELD)             ! IN Standard deviation of turbulent
!                                     fluctuations of surface layer
!                                  !  virtual temperature (K).
     &,H_BLEND_OROG(P_FIELD)       ! IN Blending height used as part of
!                                     effective roughness scheme
     &,Z0M_EFF_GB(P_FIELD)         ! IN Effective grid-box roughness
!                                     length for momentum
     &,TIMESTEP                    ! IN Timestep (seconds).

      LOGICAL LTIMER               ! Logical switch for TIMER diags

!  In/outs :-

      REAL
     & Q(P_FIELD,BL_LEVELS)        ! INOUT Input:specific humidity
!                                      ( kg/kg air).
!                                      Output:total water content
!                                      (Q)(kg/Kg air).
     &,T(P_FIELD,BL_LEVELS)        ! INOUT Input:atmospheric temp(K)
!                                      Output:liquid/frozen water
!                                      temperature (TL) (K)
     &,U(U_FIELD,BL_LEVELS)        ! INOUT W'ly wind component (m/s)
     &,V(U_FIELD,BL_LEVELS)        ! INOUT S'ly wind component (m/s)
     &,ZH(P_FIELD)                 ! INOUT Height above surface of top
!                                      of boundary layer (metres).

!  Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-
!
      REAL
     & FQW(P_FIELD,BL_LEVELS)      ! OUT Moisture flux between layers
!                                     (kg per square metre per sec).
!                                     FQW(,1) is total water flux
!                                     from surface, 'E'.
     &,FTL(P_FIELD,BL_LEVELS)      ! OUT FTL(,K) contains net turbulent
!                                     sensible heat flux into layer K
!                                     from below; so FTL(,1) is the
!                                     surface sensible heat, H. (W/m2)
     &,RHOKH(P_FIELD,BL_LEVELS)    ! OUT Exchange coeffs for moisture.
     &,RHOKM_UV(U_FIELD,BL_LEVELS) ! OUT Exchange coefficients for
!                                     momentum (on UV-grid, with 1st
!                                     and last rows undefined (or, at
!                                     present, set to "missing data"))
     &,TAUX(U_FIELD,BL_LEVELS)     ! OUT W'ly component of surface wind
!                                     stress (N/sq m).(On UV-grid with
!                                     first and last rows undefined or
!                                     at present, set to missing data
     &,TAUY(U_FIELD,BL_LEVELS)     ! OUT S'ly component of surface wind
!                                     stress (N/sq m).  On UV-grid;
!                                     comments as per TAUX.
     &,ZHT(P_FIELD)                ! OUT Height below which there may be
!                                  !     turbulent mixing (m).
     &,BL_TYPE_1(P_FIELD)          ! OUT Indicator set to 1.0 if stable
!                                  !     b.l. diagnosed, 0.0 otherwise.
     &,BL_TYPE_2(P_FIELD)          ! OUT Indicator set to 1.0 if Sc over
!                                  !     stable surface layer diagnosed,
!                                  !     0.0 otherwise.
     &,BL_TYPE_3(P_FIELD)          ! OUT Indicator set to 1.0 if well
!                                  !     mixed b.l. diagnosed,
!                                  !     0.0 otherwise.
     &,BL_TYPE_4(P_FIELD)          ! OUT Indicator set to 1.0 if
!                                  !     decoupled Sc layer (not over
!                                  !     cumulus) diagnosed,
!                                  !     0.0 otherwise.
     &,BL_TYPE_5(P_FIELD)          ! OUT Indicator set to 1.0 if
!                                  !     decoupled Sc layer over cumulus
!                                  !     diagnosed, 0.0 otherwise.
     &,BL_TYPE_6(P_FIELD)          ! OUT Indicator set to 1.0 if a
!                                  !     cumulus capped b.l. diagnosed,
!                                  !     0.0 otherwise.
     &,RHO_KM(P_FIELD,2:BL_LEVELS) ! OUT Air density * turbulent mixing
!                                     coefficient for momentum before
!                                     interpolation.

!-2 Output needed by implicit surface and boundary layer routines :-

      REAL
     & DTRDZ(P_FIELD,BL_LEVELS)    ! OUT -g.dt/dp for model layers.
     &,QW(P_FIELD,BL_LEVELS)       ! OUT Total water content, but
!                                     replaced by specific humidity
!                                     in LS_CLD.
     &,TL(P_FIELD,BL_LEVELS)       ! OUT Ice/liquid water temperature,
!                                     but replaced by T in LS_CLD.
     &,CT_CTQ(P_FIELD,BL_LEVELS)   ! OUT Coefficient in T and q
!                                        tri-diagonal implicit matrix
     &,CQ_CM(U_FIELD,BL_LEVELS)    ! OUT Coefficient in U and V
!                                        tri-diagonal implicit matrix
     &,DQW(P_FIELD,BL_LEVELS)      ! OUT BL increment to q field
     &,DTL(P_FIELD,BL_LEVELS)      ! OUT BL increment to T field
     &,DU(U_FIELD,BL_LEVELS)       ! OUT BL increment to u wind field
     &,DV(U_FIELD,BL_LEVELS)       ! OUT BL increment to v wind field
     &,RDZ(P_FIELD,BL_LEVELS)      ! OUT RDZ(,1) is the reciprocal of
!                                     the height of level 1, i.e. of
!                                     the middle of layer 1.  For K > 1,
!                                     RDZ(,K) is the reciprocal
!                                     of the vertical distance
!                                     from level K-1 to level K.
     &,RDZUV(U_FIELD,BL_LEVELS)    ! OUT  RDZ (K > 1) on UV-grid.
!                                     Comments as per RHOKM (RDZUV).


      INTEGER
     & NRML(P_FIELD)      ! DUMMY Used in 7A boundary layer scheme



!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL Z,HEAT_CON,SMC_ROOT,SF_EXCH,BOUY_TQ,BTQ_INT,
     & KMKH,EX_FLUX_TQ,EX_FLUX_UV,IM_CAL_TQ,SICE_HTF,SF_EVAP,
     & IM_CAL_UV
      EXTERNAL TIMER
*IF -DEF,SCMA
      EXTERNAL UV_TO_P,P_TO_UV
*ENDIF

!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

*CALL C_R_CP
*CALL C_G
*CALL C_LHEAT
*CALL C_GAMMA
*IF DEF,MPP
! MPP Common block
*CALL PARVARS
*ENDIF

! Derived local parameters.

      REAL LCRCP,LS,LSRCP

      PARAMETER (
     & LCRCP=LC/CP           ! Evaporation-to-dT conversion factor.
     &,LS=LF+LC              ! Latent heat of sublimation.
     &,LSRCP=LS/CP           ! Sublimation-to-dT conversion factor.
     &  )

!-----------------------------------------------------------------------

!  Workspace :-

      REAL
     & A_DQSDT(P_FIELD,BL_LEVELS)
!                               ! Saturated lapse rate factor
!                               ! on p,T,q-levels (full levels).
     &,A_DQSDTM(P_FIELD,BL_LEVELS)
!                               ! Saturated lapse rate factor
!                               ! on intermediate levels (half levels).
     &,A_QS(P_FIELD,BL_LEVELS)  ! Saturated lapse rate factor
!                               ! on p,T,q-levels (full levels).
     &,A_QSM(P_FIELD,BL_LEVELS)
!                               ! Saturated lapse rate factor
!                               ! on intermediate levels (half levels).
     &,BF(P_FIELD,BL_LEVELS)    ! A buoyancy parameter (beta F tilde)
     &,BQ(P_FIELD,BL_LEVELS)    ! A buoyancy parameter for clear air
!                               ! on p,T,q-levels (full levels).
     &,BQ_CLD(P_FIELD,BL_LEVELS)! A buoyancy parameter for cloudy air
!                               ! on p,T,q-levels (full levels).
     &,BQM(P_FIELD,BL_LEVELS)   ! A buoyancy parameter for clear air
!                               ! on intermediate levels (half levels).
     &,BQM_CLD(P_FIELD,BL_LEVELS)
!                               ! A buoyancy parameter for cloudy air
!                               ! on intermediate levels (half levels).
     &,BT(P_FIELD,BL_LEVELS)    ! A buoyancy parameter for clear air
!                               ! on p,T,q-levels (full levels).
     &,BT_CLD(P_FIELD,BL_LEVELS)
!                               ! A buoyancy parameter for cloudy air
!                               ! on p,T,q-levels (full levels).
     &,BTM(P_FIELD,BL_LEVELS)   ! A buoyancy parameter for clear air
!                               ! on intermediate levels (half levels).
     &,BTM_CLD(P_FIELD,BL_LEVELS)
!                               ! A buoyancy parameter for cloudy air
!                               ! on intermediate levels (half levels).
     &,DB(P_FIELD,2:BL_LEVELS)
!                               ! Buoyancy jump across layer interface.
     &,DELTAP(P_FIELD,BL_LEVELS)! Difference in pressure between levels
     &,DELTAP_UV(P_FIELD,BL_LEVELS)
!                                 Difference in pressure between levels
!                                 on UV points
     &,DQSDT(P_FIELD,BL_LEVELS) ! Derivative of q_SAT w.r.t. T
     &,DTRDZ_UV(U_FIELD,BL_LEVELS)
!                                 -g.dt/dp for model wind layers.
     &,DZL(P_FIELD,BL_LEVELS)   ! DZL(,K) is depth in m of layer
!                                 K, i.e. distance from boundary
!                                 K-1/2 to boundary K+1/2.
     &,DU_NT(U_FIELD,BL_LEVELS) ! non-turbulent inc. to u wind field
     &,DV_NT(U_FIELD,BL_LEVELS) ! non-turbulent inc. to v wind field
     &,DTL_NT(P_FIELD,BL_LEVELS)! non-turbulent inc. to TL field
     &,DQW_NT(P_FIELD,BL_LEVELS)! non-turbulent inc. to QW field
!
     &,GRAD_Q_ADJ(P_FIELD)      ! Humidity gradient adjustment
!                                 for non-local mixing in unstable
!                                 turbulent boundary layer.
     &,GRAD_T_ADJ(P_FIELD)      ! Temperature gradient adjustment
!                                 for non-local mixing in unstable
!                                 turbulent boundary layer.
     &,P(P_FIELD,BL_LEVELS)     ! P(*,K) is pressure at full level k.
     &,P_HALF(P_FIELD,BL_LEVELS)! P_HALF(*,K) is pressure at half
!                               ! level k-1/2.
     &,Z_FULL(P_FIELD,BL_LEVELS)! Z_FULL(*,K) is height of full level k.
     &,Z_HALF(P_FIELD,BL_LEVELS)! Z_HALF(*,K) is height of half level
!                               ! k-1/2.
     &,Z_UV(P_FIELD,BL_LEVELS)  ! Z_UV(*,K) is height of half level
!                               ! k-1/2.
     &,Z_TQ(P_FIELD,BL_LEVELS)  ! Z_TQ(*,K) is height of half level
!                               ! k+1/2.
     &,RHO_FULL(P_FIELD,BL_LEVELS)
!                               ! RHO_FULL(*,K) is the density at full
!                               ! model level k.
     &,RHO_HALF(P_FIELD,BL_LEVELS)
!                               ! RHO_HALF(*,K) is the density at half
!                               ! level k-1/2.
     &,RHO_UV(P_FIELD,BL_LEVELS)
!                               ! RHO_UV(*,K) is the density at half
!                               ! level k-1/2.
     &,RHO_TQ(P_FIELD,BL_LEVELS)
!                               ! RHO_TQ(*,K) is the density at half
!                               ! level k+1/2.
     &,RHOKHZ(P_FIELD,2:BL_LEVELS)
!                               ! Non-local turbulent mixing
!                                 coefficient for heat and moisture.
     &,RHOKH_TOP(P_FIELD,2:BL_LEVELS)
!                               ! Non-local turbulent mixing coefficient
!                               ! for top-down mixing of heat and
!                               ! moisture.
     &,RHOKM(P_FIELD,BL_LEVELS) ! Turbulent mixing coefficient for
!                                 momentum on P-grid.
     &,RHOKMZ(P_FIELD,2:BL_LEVELS)
!                               ! Non-local turbulent mixing
!                                 coefficient for momentum.
     &,RHOKM_TOP(P_FIELD,2:BL_LEVELS)
!                               ! Non-local turbulent mixing coefficient
!                               ! for top-down mixing of momentum.
     &,TV(P_FIELD,BL_LEVELS)    ! Virtual temp
     &,U_P(P_FIELD,BL_LEVELS)   ! U on P-grid.
     &,V_P(P_FIELD,BL_LEVELS)   ! V on P-grid.
     &,ZLB(P_FIELD,0:BL_LEVELS) ! ZLB(,K) is the height of the
!                                 upper boundary of layer K
!                                 ( = 0.0 for "K=0").
       REAL
     & Z_LCL(P_FIELD)           ! Height of lifting condensation level.
!
      INTEGER
     & NTML(P_FIELD)            ! Number of model levels in the
!                                 turbulently mixed layer.
     &,NTDSC(P_FIELD)           ! Top level for turbulent mixing in
!                               ! cloud layer.
      LOGICAL
     & CUMULUS(P_FIELD)         ! Logical switch for cumulus in the b.l.
     &,UNSTABLE(P_FIELD)        ! Logical switch for unstable
!                                 surface layer.
     &,DSC(P_FIELD)             ! Flag set if decoupled stratocumulus
!                               ! layer found.

!  Local scalars :-

      REAL
     & WK         ! LOCAL 0.5 * DZL(I,K) * RDZ(I,K)
     &,WKM1       ! LOCAL 0.5 * DZL(I,K-1) * RDZ(I,K)

      INTEGER
     & I          ! LOCAL Loop counter (horizontal field index).
     &,K          ! LOCAL Loop counter (vertical level index).

      IF (LTIMER) THEN
        CALL TIMER('BDYLAYR ',3)
      ENDIF
!-----------------------------------------------------------------------
!! 1.  Perform calculations in what the documentation describes as
!!     subroutine Z_DZ.  In fact, a separate subroutine isn't used.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!! 1.1 Initialise ZLB(,0) (to zero, of course, this being the height
!!     of the surface above the surface).
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1
        ZLB(I,0)=0.0
      ENDDO

!-----------------------------------------------------------------------
!! 1.2 Calculate layer depths and heights, and construct wind fields on
!!     P-grid.  This involves calling subroutines Z and UV_TO_P.
!!     Virtual temperature is also calculated, as a by-product.
!-----------------------------------------------------------------------

!  NB RDZ  TEMPORARILY used to return DELTA_Z_LOWER, the lower half
!     layer thickness

      DO K=1,BL_LEVELS
        CALL Z(P_POINTS,EXNER(P1,K),EXNER(P1,K+1),PSTAR(P1),
     &    AKH(K),BKH(K),Q(P1,K),QCF(P1,K),
     &    QCL(P1,K),T(P1,K),ZLB(P1,K-1),TV(P1,K),
     &    ZLB(P1,K),DZL(P1,K),RDZ(P1,K),LTIMER)
      ENDDO
      DO K=1,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          Z_FULL(I,K) = ZLB(I,K) - 0.5 * DZL(I,K)
          Z_HALF(I,K) = ZLB(I,K-1)
          Z_UV(I,K) = ZLB(I,K-1)
          Z_TQ(I,K) = ZLB(I,K)
       ENDDO
      ENDDO
      DO K=1,BL_LEVELS

*IF -DEF,SCMA
        CALL UV_TO_P(U(U1,K),U_P(P1,K),
     &               U_POINTS,P_POINTS,ROW_LENGTH,N_U_ROWS)
        CALL UV_TO_P(V(U1,K),V_P(P1,K),
     &               U_POINTS,P_POINTS,ROW_LENGTH,N_U_ROWS)


! du_nt 'borrowed to store dzl on uv grid
        CALL P_TO_UV (DZL(P1,K),DU_NT(U1+ROW_LENGTH,K),
     &     P_POINTS,U_POINTS,ROW_LENGTH,N_P_ROWS)

*ELSE
      DO I = P1, P1-1+P_POINTS
        U_P(i,K) = U(i,K)
        V_P(i,K) = V(i,K)
      ENDDO
*ENDIF
      ENDDO

! set pressure array.
      DO K=1,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          P(I,K) = AK(K) + BK(K)*PSTAR(I)
          P_HALF(I,K) = AKH(K) + BKH(K)*PSTAR(I)

! These will be used in new dynamics scheme - currently unused
          DTL_NT(I,K)=0.0
          DQW_NT(I,K)=0.0

        ENDDO

      ENDDO  ! end of loop over bl_levels

      DO K=BL_LEVELS,2,-1

        DO I=P1,P1+P_POINTS-1
          RDZ(I,K)=1.0/(RDZ(I,K)+(DZL(I,K-1)-RDZ(I,K-1)))
          DELTAP(I,K)=DELTA_AK(K) + PSTAR(I)*DELTA_BK(K)

          DTRDZ(I,K) = -G * TIMESTEP/ DELTAP(I,K)
!     &                  (DELTA_AK(K) + PSTAR(I)*DELTA_BK(K))
        ENDDO
      ENDDO

      DO I=P1,P1+P_POINTS-1
        RDZ(I,1)=1.0/RDZ(I,1)

        DELTAP(I,1)=DELTA_AK(1) + PSTAR(I)*DELTA_BK(1)
        DTRDZ(I,1) = -G * TIMESTEP/DELTAP(I,1)
!     &                  (DELTA_AK(1) + PSTAR(I)*DELTA_BK(1))
      ENDDO

      DO K=1,BL_LEVELS


! Calculate RDZUV here

        IF(K.GE.2)THEN
*IF -DEF,SCMA

          DO I=U1+ROW_LENGTH,U1-ROW_LENGTH+U_POINTS-1
            RDZUV(I,K) = 2.0 / ( DU_NT(I,K) + DU_NT(I,K-1) )
          ENDDO

!-----------------------------------------------------------------------
! 1.3 Set first and last rows to "missing data indicator"
!-----------------------------------------------------------------------

*IF DEF,MPP
      IF (attop) THEN
*ENDIF
        DO I=U1,U1+ROW_LENGTH-1
          RDZUV(I,K) = 1.0E30
        ENDDO
*IF DEF,MPP
      ENDIF

      IF (atbase) THEN
*ENDIF
        DO I= U1+(N_U_ROWS-1)*ROW_LENGTH, U1 + N_U_ROWS*ROW_LENGTH-1
          RDZUV(I,K) = 1.0E30
        ENDDO
*IF DEF,MPP
      ENDIF
*ENDIF


*ELSE
      DO I = P1, P1-1+P_POINTS
        RDZUV(i,K) = 2.0 / ( DZL(i,K) + DZL(i,K-1) )
      ENDDO
*ENDIF
        ENDIF   ! K .ge. 2

! Calculate DTRDZ_UV here.

*IF -DEF,SCMA
!        CALL P_TO_UV (DTRDZ(P1,K),DTRDZ_UV(U1+ROW_LENGTH,K),
!     &     P_POINTS,U_POINTS,ROW_LENGTH,N_P_ROWS)

        CALL P_TO_UV (DELTAP(P1,K),DELTAP_UV(U1+ROW_LENGTH,K),
     &     P_POINTS,U_POINTS,ROW_LENGTH,N_P_ROWS)

        DO I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1
          DTRDZ_UV(I,K) = -G * TIMESTEP / DELTAP_UV(I,K)
        ENDDO

*ELSE
      DO I = P1, P1-1+P_POINTS
        DTRDZ_UV(i,K) = DTRDZ(i,K)
      ENDDO
*ENDIF

      ENDDO ! loop over bl_levels

! "borrowed" du_nt reset to zero
! Non turbulent increments for new dynamics scheme (currently not used)
        DO K=1,BL_LEVELS
          DO I=1,U_FIELD
            DU_NT(I,K) =0.0
            DV_NT(I,K) =0.0
          ENDDO
        ENDDO

!-----------------------------------------------------------------------
!! Calculate total water content, QW and Liquid water temperature, TL
!-----------------------------------------------------------------------

      DO K=1,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          QW(I,K) = Q(I,K) + QCL(I,K) + QCF(I,K)              ! P243.10
          TL(I,K) = T(I,K) - LCRCP*QCL(I,K) - LSRCP*QCF(I,K)  ! P243.9
        ENDDO
      ENDDO


!-----------------------------------------------------------------------
!! 5.  Turbulent exchange coefficients and "explicit" fluxes between
!!     model layers in the boundary layer (P243b, routine KMKH).
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!! 5.1  Calculate bouyancy parameters BT and BQ.
!-----------------------------------------------------------------------

      CALL BOUY_TQ (
     & P_FIELD,P1
     &,P_POINTS,BL_LEVELS
     &,P,CF,T,TL,Q,QCF,QCL
     &,BT,BQ,BF,BT_CLD,BQ_CLD,A_QS,A_DQSDT,DQSDT
     &,L_BL_LSPICE,LTIMER
     &  )


!-----------------------------------------------------------------------
!! 5.2  Interpolate BT and BQ to half levels.
!-----------------------------------------------------------------------

      CALL BTQ_INT (
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,DZL,RDZ,BQ,BT,BQ_CLD,BT_CLD,A_QS,A_DQSDT
     &,BQM,BTM,BQM_CLD,BTM_CLD,A_QSM,A_DQSDTM
     &,LTIMER
     &  )


!-----------------------------------------------------------------------
!! 5.3  Calculate the diffusion coefficients Km and Kh.
!-----------------------------------------------------------------------

      DO K=1,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          RHO_FULL(I,K) =
     &     ( AK(K) + BK(K)*PSTAR(I) )      ! Pressure at K
     &     /                               ! divided by ...
     &     ( R * TV(I,K) )                 ! R times TV at K
        ENDDO
      ENDDO
      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          WKM1 = 0.5 * DZL(I,K-1) * RDZ(I,K)
          WK = 0.5 * DZL(I,K) * RDZ(I,K)
          RHO_HALF(I,K) = WK*RHO_FULL(I,K-1) + WKM1*RHO_FULL(I,K)
        ENDDO
      ENDDO
      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          RHO_UV(I,K) = RHO_HALF(I,K)
        ENDDO
      ENDDO
      DO K=1,BL_LEVELS-1
        DO I=P1,P1+P_POINTS-1
          RHO_TQ(I,K) = RHO_HALF(I,K+1)
        ENDDO
      ENDDO
      DO I=P1,P1+P_POINTS-1
        RHO_HALF(I,1) = RHO_FULL(I,1)
        RHO_UV(I,1) = RHO_FULL(I,1)
        RHO_TQ(I,BL_LEVELS) = RHO_FULL(I,BL_LEVELS)
      ENDDO

      CALL KMKHZ (
     & P_FIELD,P1,P_POINTS,BL_LEVELS,
     & P,P_HALF,T,Q,QCL,QCF,BT,BQ,CF,DZL,
     & RDZ,DELTAP,FTL,FQW,
     & Z0M_EFF_GB,Z_FULL,Z_HALF,Z_UV,Z_TQ,U_S,FB_SURF,
     & QW,RHOKMZ(1,2),DB(1,2),RHOKHZ(1,2),TL,ZH,TV1_SD,T1_SD,Q1_SD,
     & NTML,GRAD_T_ADJ,GRAD_Q_ADJ,
     & BTM,BQM,DQSDT,BTM_CLD,BQM_CLD,A_QSM,A_DQSDTM,RHO_TQ,RHO_UV,
     & RAD_HR,RADHR_DIM1,CUMULUS,Z_LCL,RHOKM_TOP(1,2),RHOKH_TOP(1,2),
     & ZHT,BL_TYPE_1,BL_TYPE_2,BL_TYPE_3,BL_TYPE_4,BL_TYPE_5,BL_TYPE_6,
     & UNSTABLE,NTDSC,DSC,
     & LTIMER
     & )

      CALL EX_COEF (
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,CCB,CCT,NTML,L_MOM
     &,CCA,DZL,RDZ,DB(1,2),U_P,V_P
     &,RHO_HALF,ZH,Z_HALF,Z0M_EFF_GB,H_BLEND_OROG
     &,CUMULUS,Z_LCL
     &,RHOKM,RHOKH
     &,LTIMER
     & )

      CALL KMKH (
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,RHOKM,RHO_KM(1,2),RHOKH
     &,RHOKMZ(1,2),RHOKHZ(1,2)
     &,NTML,CUMULUS,RHOKM_TOP(1,2),RHOKH_TOP(1,2)
     &,UNSTABLE,NTDSC,DSC
     &,LTIMER
     & )

!
!-----------------------------------------------------------------------
!! 5.4 Interpolate RHOKM's and CDR10M to uv points ready for the
!!     calculation of the explcit fluxes TAU_X and TAU_Y at levels
!!     above the surface.
!-----------------------------------------------------------------------

*IF DEF,MPP
! RHOKM(*,1) contains duff data in halos. The P_TO_UV can interpolate
! this into the real data, so first we must update east/west halos

      CALL SWAPBOUNDS(RHOKM(P1,1),ROW_LENGTH,N_U_ROWS,1,0,1)
      CALL SWAPBOUNDS(RHOKM(1,2),ROW_LENGTH,
     &                U_FIELD/ROW_LENGTH,1,1,BL_LEVELS-1)
*ENDIF

      DO K=2,BL_LEVELS

*IF -DEF,SCMA
        CALL P_TO_UV (RHOKM(P1,K),RHOKM_UV(U1+ROW_LENGTH,K),
     &     P_POINTS,U_POINTS,ROW_LENGTH,N_P_ROWS)
*IF DEF,MPP
      IF (attop) THEN
*ENDIF
        DO I=U1,U1+ROW_LENGTH-1
          RHOKM_UV(I,K) = 1.0E30
        ENDDO
*IF DEF,MPP
      ENDIF

      IF (atbase) THEN
*ENDIF
        DO I= U1+(N_U_ROWS-1)*ROW_LENGTH, U1+N_U_ROWS*ROW_LENGTH-1
          RHOKM_UV(I,K) = 1.0E30
        ENDDO
*IF DEF,MPP
      ENDIF
*ENDIF

*ELSE
      DO I = P1, P1-1+P_POINTS
        RHOKM_UV(i,K) = RHOKM(i,K)
      ENDDO
*ENDIF
      ENDDO ! loop over bl_levels

      IF (L_BL_LSPICE) THEN

        DO K = 1,BL_LEVELS
          DO I = P1,P1+P_POINTS-1
            QW(I,K) = Q(I,K) + QCL(I,K)
            TL(I,K) = T(I,K) - LCRCP * QCL(I,K)
          ENDDO
        ENDDO

      ENDIF

!-----------------------------------------------------------------------
!! 5.5 Calculation of explicit fluxes of T,Q
!-----------------------------------------------------------------------


      CALL EX_FLUX_TQ (
     &  P_POINTS,P_FIELD,P1,BL_LEVELS
     &, TL,QW,RDZ,FTL,FQW,RHOKH
     &, RHOKHZ(1,2)
     &, GRAD_T_ADJ,GRAD_Q_ADJ
     &, NTML
     &, LTIMER
     &  )

!-----------------------------------------------------------------------
!! 5.6 Calculation of explicit fluxes of U and V.
!-----------------------------------------------------------------------


      CALL EX_FLUX_UV ( ! For U
     &  U_POINTS,U_FIELD,ROW_LENGTH,BL_LEVELS,U1
     &, U,U_0,RDZUV(1,2),RHOKM_UV,TAUX
     &, LTIMER
     &  )


      CALL EX_FLUX_UV ( ! For V
     &  U_POINTS,U_FIELD,ROW_LENGTH,BL_LEVELS,U1
     &, V,V_0,RDZUV(1,2),RHOKM_UV,TAUY
     &, LTIMER
     &  )


*IF -DEF,SCMA
!-----------------------------------------------------------------------
!! Set first and last rows to "missing data indicator"
!-----------------------------------------------------------------------
      DO K=2,BL_LEVELS
*IF DEF,MPP
      IF (attop) THEN
*ENDIF
        DO I=U1,U1+ROW_LENGTH-1
          TAUX(I,K)=1.E30
          TAUY(I,K)=1.E30
        ENDDO
*IF DEF,MPP
      ENDIF

      IF (atbase) THEN
*ENDIF
        DO I= U1 + (N_U_ROWS-1)*ROW_LENGTH, U1 + N_U_ROWS*ROW_LENGTH -1
          TAUX(I,K)=1.E30
          TAUY(I,K)=1.E30
        ENDDO
*IF DEF,MPP
      ENDIF
*ENDIF
      ENDDO
*ENDIF


!-----------------------------------------------------------------------
!! 6.  "Implicit" calculation of increments for TL and QW
!-----------------------------------------------------------------------

      CALL IM_BL_PT1 (
     & P_FIELD,P1,U_FIELD,U1
     &,P_POINTS,U_POINTS,ROW_LENGTH,BL_LEVELS
     &,DTRDZ,DTRDZ_UV,RHOKH(1,2),RHOKM_UV(1,2)
     &,RDZ,RDZUV(1,2),GAMMA
     &,DQW_NT,DTL_NT,DU_NT,DV_NT
     &,FQW,FTL,TAUX,TAUY
     &,CT_CTQ,DQW,DTL,CQ_CM,DU,DV
     &,LTIMER
     &)


      IF (LTIMER) THEN
        CALL TIMER('BDYLAYR ',4)
      ENDIF

      RETURN
      END
*ENDIF
*DECK BDY_IMPL8A
*IF DEF,A03_7A,OR,DEF,A03_8A
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
!!!  SUBROUTINE BDY_IMPL-----------------------------------------------
!!!
!!!  Purpose: Calculate implicit correction to boundary layer fluxes of
!!!           heat, moisture and momentum.
!!!
!!!
!!! F.Hewer     <- programmer of some or all of previous code or changes
!!! C.Wilson    <- programmer of some or all of previous code or changes
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   4.3  7/2/97     New deck. S Jackson
!!!   4.4 25/6/97     Modified for MOSES II tile model. R Essery
!!!   4.4 25/6/97     Move grid definitions up to BL_INTCT.  R.A.Betts
!!!  4.5    Jul. 98  Kill the IBM specific lines. (JCThil)
!!!   4.5  7/5/98     Set TSTAR, SNOW_SURF_HTF and SOIL_SURF_HTF to 0
!!!                   at all land points, to avoid problems of
!!!                   non-initialised data.  R.A.Betts
!!!   4.5 21/5/98     Add optional error check for negative surface
!!!                   temperature.  R.A.Betts
!!!
!!!  Programming standard: Unified Model Documentation Paper No 4,
!!!                        Version ?, dated ?.
!!!
!!!  System component covered: P24.
!!!
!!!  Project task:
!!!
!!!  Documentation: UMDP 24.
!!!
!!!---------------------------------------------------------------------

!    Arguments :-
      SUBROUTINE BDY_IMPL (

! IN values defining field dimensions and subset to be processed :
     & P_FIELD,P1,U_FIELD,U1,P_POINTS,U_POINTS,ROW_LENGTH,

! IN values defining vertical grid of model atmosphere :
     & BL_LEVELS,

! IN data :
     & RHOKH,RHOKM_UV,RDZ,RDZUV,

! INOUT data :
     & Q,T,U,V,QW,TL,FQW,FTL,TAUX,TAUY,
     & DU,DV,CT_CTQ,DQW,DTL,CQ_CM,

! LOGICAL LTIMER
     & LTIMER
     & )

      IMPLICIT NONE

!  Inputs :-

! (a) Defining horizontal grid and subset thereof to be processed.

      INTEGER
     & P_FIELD                     ! IN No. of P-points in whole grid
!                                     (for dimensioning only).
     &,P1                          ! IN First point to be processed in
!                                       P-grid.
     &,U_FIELD                     ! IN No. of UV-points in whole grid.
!                                     (Checked for consistency with
!                                     P_FIELD and P_ROWS; there must
!                                     be 1 less UV than P row.)
     &,U1                          ! IN First point to be processed in
!                                       U_V-grid.
     &,P_POINTS                    ! IN Number of P-grid points to be
!                                       processed.
     &,U_POINTS                    ! IN Number of U_V-grid points.
     &,ROW_LENGTH                  ! IN No. of points in one row.
!                                     (Checked for consistency with
!                                     P_FIELD and N_ROWS.)

! (b) Defining vertical grid of model atmosphere.

      INTEGER
     & BL_LEVELS                   ! IN Max. no. of "boundary" levels
!                                     allowed.Assumed <= 30 for dim-
!                                     sioning of GAMMA in common deck
!                                     C_GAMMA used in SF_EXCH and KMKH

!  In :-

      REAL
     & RHOKH(P_FIELD,BL_LEVELS)    ! IN Exchange coeffs for moisture.
     &,RHOKM_UV(U_FIELD,BL_LEVELS) ! IN Exchange coefficients for
!                                     momentum (on UV-grid, with 1st
!                                     and last rows undefined (or, at
!                                     present, set to "missing data"))
     &,RDZ(P_FIELD,BL_LEVELS)      ! IN RDZ(,1) is the reciprocal of the
!                                     height of level 1, i.e. of the
!                                     middle of layer 1.  For K > 1,
!                                     RDZ(,K) is the reciprocal
!                                     of the vertical distance
!                                     from level K-1 to level K.
     &,RDZUV(U_FIELD,BL_LEVELS)    ! IN  RDZ (K > 1) on UV-grid.
!                                     Comments as per RHOKM (RDZUV).

      LOGICAL LTIMER               ! Logical switch for TIMER diags

!  In/outs :-

      REAL
     & Q(P_FIELD,BL_LEVELS)        ! INOUT Input:specific humidity
!                                      ( kg/kg air).
!                                      Output:total water content
!                                      (Q)(kg/Kg air).
     &,T(P_FIELD,BL_LEVELS)        ! INOUT Input:atmospheric temp(K)
!                                      Output:liquid/frozen water
!                                      temperature (TL) (K)
     &,U(U_FIELD,BL_LEVELS)        ! INOUT W'ly wind component (m/s)
     &,V(U_FIELD,BL_LEVELS)        ! INOUT S'ly wind component (m/s)
     &,QW(P_FIELD,BL_LEVELS)       ! INOUT Total water content, but
!                                      replaced by specific humidity
!                                      in LS_CLD.
     &,TL(P_FIELD,BL_LEVELS)       ! INOUT Ice/liquid water temperature,
!                                      but replaced by T in LS_CLD.
     &,FQW(P_FIELD,BL_LEVELS)      ! INOUT Moisture flux between layers
!                                      (kg per square metre per sec).
!                                      FQW(,1) is total water flux
!                                      from surface, 'E'.
     &,FTL(P_FIELD,BL_LEVELS)      ! INOUT FTL(,K) contains net
!                                      turbulent sensible heat flux into
!                                      layer K from below; so FTL(,1) is
!                                      the surface sensible heat, H.
!                                      (W/m2)
     &,TAUX(U_FIELD,BL_LEVELS)     ! INOUT W'ly component of surface
!                                      wind stress (N/sq m).(On UV-grid
!                                      with first and last rows
!                                      undefined or at present, set to
!                                      missing data
     &,TAUY(U_FIELD,BL_LEVELS)     ! INOUT S'ly component of surface
!                                      wind stress (N/sq m).  On
!                                      UV-grid; comments as per TAUX.
     &,CT_CTQ(P_FIELD,BL_LEVELS)   ! INOUT Coefficient in T and q
!                                          tri-diagonal implicit matrix
     &,CQ_CM(U_FIELD,BL_LEVELS)    ! INOUT Coefficient in U and V
!                                          tri-diagonal implicit matrix
     &,DQW(P_FIELD,BL_LEVELS)      ! INOUT BL increment to q field
     &,DTL(P_FIELD,BL_LEVELS)      ! INOUT BL increment to T field
     &,DU(U_FIELD,BL_LEVELS)       ! INOUT BL increment to u wind field
     &,DV(U_FIELD,BL_LEVELS)       ! INOUT BL increment to v wind field

!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

*CALL C_R_CP
*CALL C_GAMMA
*IF DEF,MPP
! MPP Common block
*CALL PARVARS
*ENDIF

!  Local scalars :-

      INTEGER
     & I          ! LOCAL Loop counter (horizontal field index).
     &,K          ! LOCAL Loop counter (vertical level index).

      IF (LTIMER) THEN
        CALL TIMER('BDYLAYR ',3)
      ENDIF

      CALL IM_BL_PT2 (
     & P_FIELD,P1,U_FIELD,U1
     &,P_POINTS,U_POINTS,ROW_LENGTH,BL_LEVELS
     &,RHOKH(1,2),RHOKM_UV(1,2)
     &,RDZ,RDZUV(1,2),GAMMA
     &,CT_CTQ,DQW,DTL,CQ_CM,DU,DV
     &,FQW,FTL,TAUX,TAUY
     &,QW,TL,U,V
     &,LTIMER
     &)


!-----------------------------------------------------------------------
!! 6.1 Convert FTL to sensible heat flux in Watts per square metre.
!      Also, IMPL_CAL only updates FTL_TILE(*,1) and FQW_TILE(*,1)
!      over sea points, so copy this to remaining tiles
!-----------------------------------------------------------------------

      DO K=2,BL_LEVELS
Cfpp$ Select(CONCUR)
        DO  I=P1,P1+P_POINTS-1
          FTL(I,K) = FTL(I,K)*CP
        ENDDO
      ENDDO

!7.1 Copy T and Q from workspace to INOUT space.

      DO K=1,BL_LEVELS
Cfpp$  Select(CONCUR)
        DO I=P1,P1+P_POINTS-1
          T(I,K)=TL(I,K)
          Q(I,K)=QW(I,K)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! 10 Set RHOKH, the coefficients required for tracer mixing.
!    Required 5B and after due to change in contents of RHOKH in rest
!    of routine.
!-----------------------------------------------------------------------

      DO K = 2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          RHOKH(I,K) = GAMMA(K)*RHOKH(I,K)*RDZ(I,K)
        ENDDO
      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('BDYLAYR ',4)
      ENDIF

      RETURN
      END
*ENDIF
*DECK IMBLPT18A
*IF DEF,A03_7A,OR,DEF,A03_8A
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
!!!  SUBROUTINE IM_BL_PT1 ----------------------------------------------
!!!
!!!  Purpose: Calculate increments for
!!!           T and Q in the boundary layer, using an
!!!           implicit numerical scheme.  The tridiagonal matrices are
!!!           inverted using simple Gaussian elimination.
!!!
!!!
!!!  Model           Modification history
!!! version  Date
CLL  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!!!
!!!  JJ  <- Programmers of some or all of previous code or changes
!!!
!!!
!!!  Programming standard: UM Documentation Paper No 4, Version 2,
!!!                        dated 18/1/90
!!!
!!!  System component covered: P244
!!!
!!!  Project task: P24
!!!
!!!  Documentation: UM Documentation Paper No 24.
!!!
!!!---------------------------------------------------------------------
!!  Arguments :-
      SUBROUTINE IM_BL_PT1 (
     & P_FIELD,P1,U_V_FIELD,U1_V1
     &,P_POINTS,U_V_POINTS,ROW_LENGTH,BL_LEVELS
     &,DTRDZ,DTRDZ_U_V,RHOKH,RHOKM_U_V
     &,RDZ,RDZ_U_V,GAMMA
     &,DQW_NT,DTL_NT,DU_NT,DV_NT
     &,FQW,FTL,TAU_X,TAU_Y
     &,CT_CTQ,DQW,DTL,CQ_CM,DU,DV
     &,LTIMER
     &)

      IMPLICIT NONE

      LOGICAL LTIMER

      INTEGER
     & P_FIELD                     ! IN No. of points in P-grid.
     &,P1                          ! IN First point to be processed in
!                                       P-grid.
     &,U_V_FIELD                   ! IN No. of points in U_V-grid.
     &,U1_V1                       ! IN First point to be processed in
!                                       U_V-grid.
     &,P_POINTS                    ! IN Number of P-grid points to be
!                                       processed.
     &,U_V_POINTS                  ! IN Number of U_V-grid points.
     &,ROW_LENGTH                  ! IN No. of points in latitude row.
     &,BL_LEVELS                   ! IN No. of atmospheric levels for
!                                       which boundary layer fluxes are
!                                       calculated.

      REAL
     & DTRDZ(P_FIELD,BL_LEVELS)    ! IN dz for bottom BL_LEVELS
     &,DTRDZ_U_V(U_V_FIELD,BL_LEVELS)
!                                  ! IN -g.dt/dp for model wind layers
     &,RHOKH(P_FIELD,2:BL_LEVELS)  ! IN Exchange coeff for FTL above
!                                       surface.
     &,RHOKM_U_V(U_V_FIELD,2:BL_LEVELS)
!                                  ! IN Exchange coefficients for
!                                       momentum, on UV-grid with
!                                       first and last rows ignored.
!                                       for K>=2 (from KMKH).
     &,RDZ(P_FIELD,BL_LEVELS)      ! IN 1./dz
     &,RDZ_U_V(U_V_FIELD,2:BL_LEVELS)
!                                  ! IN Reciprocal of the vertical
!                                       distance from level K-1 to
!                                       level K. (K > 1) on wind levels
     &,GAMMA(BL_LEVELS)            ! IN Implicit weighting.
     &,DQW_NT(P_FIELD,BL_LEVELS)   ! IN Non-turbulent increment for QW.
     &,DTL_NT(P_FIELD,BL_LEVELS)   ! IN Non-turbulent increment for TL.
     &,DU_NT(U_V_FIELD,BL_LEVELS)
!                                  ! IN u non-turbulent increments.
     &,DV_NT(U_V_FIELD,BL_LEVELS)
!                                  ! IN v non-turbulent increments.
     &,FQW(P_FIELD,BL_LEVELS)      ! IN Flux of QW (ie., for surface,
!                                       total evaporation). Kg/sq m/s
     &,FTL(P_FIELD,BL_LEVELS)      ! IN Flux of TL (ie., for surface,
!                                       H/Cp where H is sensible heat
!                                       in W per sq m).
     &,TAU_X(U_V_FIELD,BL_LEVELS)  ! IN x-component of turbulent
!                                       stress at levels k-1/2;
!                                       eg. TAUX(,1) is surface stress.
!                                       UV-grid, 1st and last rows set
!                                       to "missing data". (N/sq m)
!                                       IN as "explicit" fluxes from
!                                       ex_flux_uv, OUT as "implicit
     &,TAU_Y(U_V_FIELD,BL_LEVELS)  ! IN y-component of turbulent
!                                       stress at levels k-1/2;
!                                       eg. TAUX(,1) is surface stress.
!                                       UV-grid, 1st and last rows set
!                                       to "missing data". (N/sq m)
!                                       IN as "explicit" fluxes from
!                                       ex_flux_uv, OUT as "implicit


      REAL
     & CT_CTQ(P_FIELD,BL_LEVELS)   ! OUT Coefficient in T and q
!                                        tri-diagonal implicit matrix
     &,CQ_CM(U_V_FIELD,BL_LEVELS)  ! OUT Coefficient in U and V
!                                        tri-diagonal implicit matrix
     &,DQW(P_FIELD,BL_LEVELS)      ! OUT BL increment to q field
     &,DTL(P_FIELD,BL_LEVELS)      ! OUT BL increment to T field
     &,DU(U_V_FIELD,BL_LEVELS)     ! OUT BL increment to u wind field
     &,DV(U_V_FIELD,BL_LEVELS)     ! OUT BL increment to v wind field

!  External references :-
      EXTERNAL TIMER

!  Local scalars :-
      REAL
     & AT       ! Matrix element in "T" row in eqn P244.79.
     &,RBT      ! Reciprocal of BT' (eqns P244.107, 110, 113).
     &,AM       ! Matrix element in eqn P244.80.
     &,RBM      ! Reciprocal of BM(') (eqns P244.81, 85, 89).

      INTEGER
     & BLM1     ! BL_LEVELS minus 1.
     &,I        ! Loop counter (horizontal field index).
     &,K        ! Loop counter (vertical index).





!-----------------------------------------------------------------------
!!  0.  Check that the scalars input to define the grid are consistent.
!       See comments to routine SF_EXCH for details.
!-----------------------------------------------------------------------

      IF (LTIMER) THEN
        CALL TIMER('IM_BL_PT1 ',3)
      ENDIF

      BLM1 = BL_LEVELS-1


      DO I=P1,P1+P_POINTS-1
! Include non-turbulent increments.
        DQW(I,BL_LEVELS) = DTRDZ(I,BL_LEVELS) * FQW(I,BL_LEVELS)
     &                     +DQW_NT(I,BL_LEVELS)
        DTL(I,BL_LEVELS) = DTRDZ(I,BL_LEVELS) * FTL(I,BL_LEVELS)
     &                     +DTL_NT(I,BL_LEVELS)

        CT_CTQ(I,BL_LEVELS) = -DTRDZ(I,BL_LEVELS) *
     &         GAMMA(BL_LEVELS)*RHOKH(I,BL_LEVELS)*
     &          RDZ(I,BL_LEVELS)

        RBT = 1.0 / ( 1.0 - CT_CTQ(I,BL_LEVELS) )

        DQW(I,BL_LEVELS) = RBT * DQW(I,BL_LEVELS)
        DTL(I,BL_LEVELS) = RBT * DTL(I,BL_LEVELS)

        CT_CTQ(I,BL_LEVELS) = RBT * CT_CTQ(I,BL_LEVELS)         ! P244.1
      ENDDO


      DO K=BLM1,2,-1
        DO I=P1,P1+P_POINTS-1

            DQW(I,K) = -DTRDZ(I,K) * ( FQW(I,K+1) - FQW(I,K) )
     &                + DQW_NT(I,K)
            DTL(I,K) = -DTRDZ(I,K) * ( FTL(I,K+1) - FTL(I,K) )
     &                + DTL_NT(I,K)

            AT = -DTRDZ(I,K) * GAMMA(K+1)*RHOKH(I,K+1)*RDZ(I,K+1)

            CT_CTQ(I,K) = -DTRDZ(I,K) * GAMMA(K)*RHOKH(I,K)*RDZ(I,K)

            RBT = 1.0 / ( 1.0 - CT_CTQ(I,K) -
     &                             AT*( 1.0 + CT_CTQ(I,K+1) ) )

            DQW(I,K) = RBT * (DQW(I,K) - AT*DQW(I,K+1) )
            DTL(I,K) = RBT * (DTL(I,K) - AT*DTL(I,K+1) )

            CT_CTQ(I,K) = RBT * CT_CTQ(I,K)                     ! P244.1
        ENDDO ! P_points
      ENDDO !blm1,2,-1

!-----------------------------------------------------------------------
!! 3.3 Bottom model layer QW row of matrix equation.
!-----------------------------------------------------------------------

       DO I=P1,P1+P_POINTS-1

            DQW(I,1) = -DTRDZ(I,1) * FQW(I,2) + DQW_NT(I,1)
            DTL(I,1) = -DTRDZ(I,1) * FTL(I,2) + DTL_NT(I,1)

            AT = -DTRDZ(I,K) * GAMMA(2)*RHOKH(I,2)*RDZ(I,2)

            RBT = 1.0 / ( 1.0 - AT*( 1.0 + CT_CTQ(I,2) ) )

            DQW(I,1) = RBT * (DQW(I,1) - AT*DQW(I,2) )
            DTL(I,1) = RBT * (DTL(I,1) - AT*DTL(I,2) )
!
! Now set CT_CTQ(1) to be BETA
            CT_CTQ(I,1) = -DTRDZ(I,1) * RBT

      ENDDO ! P_points


*IF -DEF,SCMA
        DO I=U1_V1+ROW_LENGTH,U1_V1+U_V_POINTS-ROW_LENGTH-1
*ELSE
        DO I=1,U_V_POINTS
*ENDIF

        DU(I,BL_LEVELS) = -DTRDZ_U_V(I,BL_LEVELS) *
     &                     TAU_X(I,BL_LEVELS)
        DV(I,BL_LEVELS) = -DTRDZ_U_V(I,BL_LEVELS) *
     &                     TAU_Y(I,BL_LEVELS)

! addition of non-turbulent increments
        DU(I,BL_LEVELS) = DU(I,BL_LEVELS)
     &                       + DU_NT(I,BL_LEVELS)
        DV(I,BL_LEVELS) = DV(I,BL_LEVELS)
     &                       + DV_NT(I,BL_LEVELS)

        CQ_CM(I,BL_LEVELS) = -DTRDZ_U_V(I,BL_LEVELS) *
     &        GAMMA(BL_LEVELS)*
     &        RHOKM_U_V(I,BL_LEVELS)*RDZ_U_V(I,BL_LEVELS)

        RBM = 1.0 / ( 1.0 - CQ_CM(I,BL_LEVELS) )

        DU(I,BL_LEVELS) = RBM * DU(I,BL_LEVELS)
        DV(I,BL_LEVELS) = RBM * DV(I,BL_LEVELS)

        CQ_CM(I,BL_LEVELS) = RBM * CQ_CM(I,BL_LEVELS)
      ENDDO


      DO K=BLM1,2,-1

*IF -DEF,SCMA
        DO I=U1_V1+ROW_LENGTH,U1_V1+U_V_POINTS-ROW_LENGTH-1
*ELSE
        DO I=1,U_V_POINTS
*ENDIF

          DU(I,K) = DTRDZ_U_V(I,K) *
     &                   ( TAU_X(I,K+1) - TAU_X(I,K) )
          DV(I,K) = DTRDZ_U_V(I,K) *
     &                   ( TAU_Y(I,K+1) - TAU_Y(I,K) )
! addition of non-turbulent increments
          DU(I,K) = DU(I,K) + DU_NT(I,K)
          DV(I,K) = DV(I,K) + DV_NT(I,K)

          AM = -DTRDZ_U_V(I,K) * GAMMA(K+1)*RHOKM_U_V(I,K+1)*
     &                    RDZ_U_V(I,K+1)

          CQ_CM(I,K) = -DTRDZ_U_V(I,K) * GAMMA(K)*RHOKM_U_V(I,K)*
     &          RDZ_U_V(I,K)

          RBM = 1.0 / ( 1.0 - CQ_CM(I,K) -AM*( 1.0 + CQ_CM(I,K+1) ) )

          DU(I,K) = RBM * ( DU(I,K) - AM*DU(I,K+1) )
          DV(I,K) = RBM * ( DV(I,K) - AM*DV(I,K+1) )

          CQ_CM(I,K) = RBM * CQ_CM(I,K)
        ENDDO !loop over u_v_points
      ENDDO ! loop over 2,BLM1


*IF -DEF,SCMA
        DO I=U1_V1+ROW_LENGTH,U1_V1+U_V_POINTS-ROW_LENGTH-1
*ELSE
        DO I=1,U_V_POINTS
*ENDIF

        DU(I,1) = DTRDZ_U_V(I,1) * TAU_X(I,2)
        DV(I,1) = DTRDZ_U_V(I,1) * TAU_Y(I,2)

! addition of non-turbulent increments
        DU(I,1) = DU(I,1) + DU_NT(I,1)
        DV(I,1) = DV(I,1) + DV_NT(I,1)

        AM = -DTRDZ_U_V(I,1) * GAMMA(2)*RHOKM_U_V(I,2)
     &               *RDZ_U_V(I,2)

        RBM = 1.0 / ( 1.0 - AM *( 1.0 + CQ_CM(I,2) ) )

        DU(I,1) = RBM * ( DU(I,1) - AM*DU(I,2) )
        DV(I,1) = RBM * ( DV(I,1) - AM*DV(I,2) )

        CQ_CM(I,1) = DTRDZ_U_V(I,1) * RBM
      ENDDO ! loop over U_V_POINTS

      IF (LTIMER) THEN
        CALL TIMER('IM_BL_PT1 ',4)
      ENDIF

      RETURN
      END
*ENDIF
*DECK IMBLPT28A
*IF DEF,A03_7A,OR,DEF,A03_8A
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
!!!  SUBROUTINE IM_BL_PT2 ----------------------------------------------
!!!
!!!  Purpose: Calculate increments for
!!!           T and Q in the boundary layer, using an
!!!           implicit numerical scheme.  The tridiagonal matrices are
!!!           inverted using simple Gaussian elimination.
!!!
!!!
!!!  Model           Modification history
!!! version  Date
CLL  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!!!
!!!  JJ  <- Programmers of some or all of previous code or changes
!!!
!!!
!!!  Programming standard: UM Documentation Paper No 4, Version 2,
!!!                        dated 18/1/90
!!!
!!!  System component covered: P244
!!!
!!!  Project task: P24
!!!
!!!  Documentation: UM Documentation Paper No 24.
!!!
!!!---------------------------------------------------------------------
!!  Arguments :-
      SUBROUTINE IM_BL_PT2 (
     & P_FIELD,P1,U_V_FIELD,U1_V1
     &,P_POINTS,U_V_POINTS,ROW_LENGTH,BL_LEVELS
     &,RHOKH,RHOKM_U_V
     &,RDZ,RDZ_U_V,GAMMA
     &,CT_CTQ,DQW,DTL,CQ_CM,DU,DV
     &,FQW,FTL,TAU_X,TAU_Y
     &,QW,TL,U,V
     &,LTIMER
     &)

      IMPLICIT NONE

      LOGICAL LTIMER

      INTEGER
     & P_FIELD                     ! IN No. of points in P-grid.
     &,P1                          ! IN First point to be processed in
!                                       P-grid.
     &,U_V_FIELD                   ! IN No. of points in U_V-grid.
     &,U1_V1                       ! IN First point to be processed in
!                                       U_V-grid.
     &,P_POINTS                    ! IN Number of P-grid points to be
!                                       processed.
     &,U_V_POINTS                  ! IN Number of U_V-grid points.
     &,ROW_LENGTH                  ! IN No. of points in latitude row.
     &,BL_LEVELS                   ! IN No. of atmospheric levels for
!                                       which boundary layer fluxes are
!                                       calculated.

      REAL
     & RHOKH(P_FIELD,2:BL_LEVELS)  ! IN Exchange coeff for FTL above
!                                       surface.
     &,RHOKM_U_V(U_V_FIELD,2:BL_LEVELS)
!                                  ! IN Exchange coefficients for
!                                       momentum, on UV-grid with
!                                       first and last rows ignored.
!                                       for K>=2 (from KMKH).
     &,RDZ(P_FIELD,BL_LEVELS)      ! IN 1./dz
     &,RDZ_U_V(U_V_FIELD,2:BL_LEVELS)
!                                  ! IN Reciprocal of the vertical
!                                       distance from level K-1 to
!                                       level K. (K > 1) on wind levels
     &,GAMMA(BL_LEVELS)            ! IN Implicit weighting.


      REAL
     & CT_CTQ(P_FIELD,BL_LEVELS)   ! INOUT Coefficient in T and q
!                                          tri-diagonal implicit matrix
     &,CQ_CM(U_V_FIELD,BL_LEVELS)  ! INOUT Coefficient in U and V
!                                          tri-diagonal implicit matrix
     &,DQW(P_FIELD,BL_LEVELS)      ! INOUT BL increment to q field
     &,DTL(P_FIELD,BL_LEVELS)      ! INOUT BL increment to T field
     &,DU(U_V_FIELD,BL_LEVELS)     ! INOUT BL increment to u wind field
     &,DV(U_V_FIELD,BL_LEVELS)     ! INOUT BL increment to v wind field
     &,FQW(P_FIELD,BL_LEVELS)      ! INOUT Flux of QW (ie., for surface,
!                                          total evaporation). Kg/sq m/s
     &,FTL(P_FIELD,BL_LEVELS)      ! INOUT Flux of TL (ie., for surface,
!                                          H/Cp where H is sensible heat
!                                          in W per sq m).
     &,TAU_X(U_V_FIELD,BL_LEVELS)  ! INOUT x-component of turbulent
!                                          stress at levels k-1/2;
!                                          eg. TAUX(,1) is surface stres
!                                          UV-grid, 1st and last rows se
!                                          to "missing data". (N/sq m)
!                                          IN as "explicit" fluxes from
!                                          ex_flux_uv, OUT as "implicit
     &,TAU_Y(U_V_FIELD,BL_LEVELS)  ! INOUT y-component of turbulent
!                                          stress at levels k-1/2;
!                                          eg. TAUX(,1) is surface stres
!                                          UV-grid, 1st and last rows se
!                                          to "missing data". (N/sq m)
!                                          IN as "explicit" fluxes from
!                                          ex_flux_uv, OUT as "implicit
     &,QW(P_FIELD,BL_LEVELS)       ! INOUT Total water content (kg per
!                                          kg air).  From P243.
     &,TL(P_FIELD,BL_LEVELS)       ! INOUT Liquid/frozen water
!                                          temperature (K).  From P243.
     &,U(U_V_FIELD,BL_LEVELS)      ! INOUT delta (U) elements of
!                                          vector on RHS, then LHS, of
!                                          eqn P244.80.
     &,V(U_V_FIELD,BL_LEVELS)      ! INOUT delta (V) elements of
!                                          vector on RHS, then LHS, of
!                                          eqn P244.80.


*CALL C_R_CP

!  External references :-
      EXTERNAL TIMER

!  Local scalars :-
      INTEGER
     & I        ! Loop counter (horizontal field index).
     &,K        ! Loop counter (vertical index).




!-----------------------------------------------------------------------
!!  0.  Check that the scalars input to define the grid are consistent.
!       See comments to routine SF_EXCH for details.
!-----------------------------------------------------------------------

      IF (LTIMER) THEN
        CALL TIMER('IM_BL_PT2 ',3)
      ENDIF

*IF -DEF,SCMA
      DO I=U1_V1+ROW_LENGTH,U1_V1+U_V_POINTS-ROW_LENGTH-1
*ELSE
      DO I=1,U_V_POINTS
*ENDIF

        DU(I,1) = DU(I,1) - CQ_CM(I,1)*TAU_X(I,1)
        U(I,1) = U(I,1) + DU(I,1)
        DV(I,1) = DV(I,1) - CQ_CM(I,1)*TAU_Y(I,1)
        V(I,1) = V(I,1) + DV(I,1)

      ENDDO


      DO K=2,BL_LEVELS

*IF -DEF,SCMA
        DO I=U1_V1+ROW_LENGTH,U1_V1+U_V_POINTS-ROW_LENGTH-1
*ELSE
        DO I=1,U_V_POINTS
*ENDIF

          DU(I,K) = DU(I,K) - CQ_CM(I,K)*DU(I,K-1)
          U(I,K) = U(I,K) + DU(I,K)
          DV(I,K) = DV(I,K) - CQ_CM(I,K)*DV(I,K-1)
          V(I,K) = V(I,K) + DV(I,K)

        ENDDO
      ENDDO


      DO I=P1,P1+P_POINTS-1
        DTL(I,1) = DTL(I,1) - CT_CTQ(I,1)*FTL(I,1)/CP
        TL(I,1) = TL(I,1) + DTL(I,1)
        DQW(I,1) = DQW(I,1) - CT_CTQ(I,1)*FQW(I,1)
        QW(I,1) = QW(I,1) + DQW(I,1)
      ENDDO !p_points


      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1

          DTL(I,K) = DTL(I,K) - CT_CTQ(I,K)*DTL(I,K-1)
          TL(I,K) = TL(I,K) + DTL(I,K)
          DQW(I,K) = DQW(I,K) - CT_CTQ(I,K)*DQW(I,K-1)
          QW(I,K) = QW(I,K) + DQW(I,K)

        ENDDO !p_points
      ENDDO !bl_levels


      DO K=2,BL_LEVELS
*IF -DEF,SCMA
        DO I=U1_V1+ROW_LENGTH,U1_V1+U_V_POINTS-ROW_LENGTH-1
*ELSE
        DO I=1,U_V_POINTS
*ENDIF

          TAU_X(I,K) = TAU_X(I,K) +
     &    GAMMA(K) * RHOKM_U_V(I,K) * RDZ_U_V(I,K)
     &                        *( DU(I,K) - DU(I,K-1) )

          TAU_Y(I,K) = TAU_Y(I,K) +
     &    GAMMA(K) * RHOKM_U_V(I,K) * RDZ_U_V(I,K)
     &                        *( DV(I,K) - DV(I,K-1) )

        ENDDO !u_v_points
      ENDDO ! bl_levels


      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1

!  Calculate and store fluxes due to local mixing.
!  FTL(local mixing) stored in array AT,
!  FQW(local mixing) stored in array AQ_AM.

          FTL(I,K) = FTL(I,K) - GAMMA(K)*RHOKH(I,K)*RDZ(I,K)
     &                            * ( DTL(I,K) - DTL(I,K-1) )
          FQW(I,K) = FQW(I,K) - GAMMA(K)*RHOKH(I,K)*RDZ(I,K)
     &                            * ( DQW(I,K) - DQW(I,K-1) )

        ENDDO ! p_points
      ENDDO ! bl_levels

      IF (LTIMER) THEN
        CALL TIMER('IM_BL_PT2 ',4)
      ENDIF

      RETURN
      END
*ENDIF
*DECK IMSFPT8A                                                          
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
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
!!!  SUBROUTINE IM_SF_PT ---------------------------------------------- 
!!!                                                                     
!!!  Purpose: Calculate increments for                                  
!!!           T and Q in the boundary layer, using an                   
!!!           implicit numerical scheme.  The tridiagonal matrices are  
!!!           inverted using simple Gaussian elimination.               
!!!                                                                     
!!!                                                                     
!!!  Model           Modification history                               
!!! version  Date                                                       
CLL  4.5    Jul. 98  Kill the IBM specific lines (JCThil)               
!!!                                                                     
!!!  JJ  <- Programmers of some or all of previous code or changes      
!!!                                                                     
!!!                                                                     
!!!  Programming standard: UM Documentation Paper No 4, Version 2,      
!!!                        dated 18/1/90                                
!!!                                                                     
!!!  System component covered: P244                                     
!!!                                                                     
!!!  Project task: P24                                                  
!!!                                                                     
!!!  Documentation: UM Documentation Paper No 24.                       
!!!                                                                     
!!!---------------------------------------------------------------------
!!  Arguments :-                                                        
      SUBROUTINE IM_SF_PT (                                             
     & P_FIELD,P1,U_V_FIELD,U1_V1                                       
     &,P_POINTS,U_V_POINTS,ROW_LENGTH,LAND_FIELD                        
     &,LAND_INDEX,NTILES,TILE_INDEX,TILE_PTS                            
     &,FLANDG,TILE_FRAC,SNOW_TILE,ICE_FRACT                             
     &,GAMMA_1,ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE                      
     &,RESFT,RHOKPM,RHOKPM_SICE                                         
     &,RHOKM_1_U_V,RHOKH_1,RHOKH1_SICE                                  
     &,CT_CTQ_1,DQW_1,DTL_1,CQ_CM_1,DU_1,DV_1
     &,FLANDG_UV                                                
     &,FQW_GB,FTL_GB
     &,TAUX_1,TAUX_LAND,TAUX_SSI,TAUY_1,TAUY_LAND,TAUY_SSI              
     &,FQW_TILE,FTL_TILE,FQW_ICE,FTL_ICE,E_SEA,H_SEA                    
     &,LTIMER                                                           
     &)                                                                 
                                                                        
                                                                        
      IMPLICIT NONE                                                     
                                                                        
      LOGICAL LTIMER                                                    
                                                                        
      INTEGER                                                           
     & P_FIELD                     ! IN No. of points in P-grid.        
     &,P1                          ! IN First point to be processed in  
!                                       P-grid.                         
     &,U_V_FIELD                   ! IN No. of points in U_V-grid.      
     &,U1_V1                       ! IN First point to be processed in  
!                                       U_V-grid.                       
     &,P_POINTS                    ! IN Number of P-grid points to be   
!                                       processed.                      
     &,U_V_POINTS                  ! IN Number of U_V-grid points.      
     &,ROW_LENGTH                  ! IN No. of points in latitude row.  
     &,LAND_FIELD                  ! IN Total number of land points.    
     &,LAND_INDEX(P_FIELD)         ! IN LAND_INDEX(I)=J => the Jth      
!                                       point in P_FIELD is the Ith     
!                                       land point.                     
     &,NTILES                      ! IN Number of land surface tiles.   
     &,TILE_INDEX(LAND_FIELD,NTILES)!IN Index of tile points.           
     &,TILE_PTS(NTILES)            ! IN Number of tiles.                
                                                                        
                                                                        
      REAL                                                              
     & FLANDG(P_FIELD)             ! IN Land fraction                   
     &,TILE_FRAC(LAND_FIELD,NTILES)! IN Tile fraction                   
     &,SNOW_TILE(LAND_FIELD,NTILES)! IN Lying snow on land tiles (kg/m2)
     &,ICE_FRACT(P_FIELD)          ! IN Fraction of grid-box which is   
!                                       sea-ice (decimal fraction).     
     &,GAMMA_1                     ! IN Implicit weighting.             
     &,ALPHA1(LAND_FIELD,NTILES)   ! IN Gradient of saturated specific  
!                                       humidity with respect to        
!                                       temperature between the bottom  
!                                       model layer and the surface.    
     &,ALPHA1_SICE(P_FIELD)        ! IN ALPHA1 for sea-ice              
     &,ASHTF(P_FIELD)              ! IN Coefficient to calculate surface
!                                       heat flux into soil or sea-ice  
!                                       (W/m2/K).                       
                                                                        
     &,ASHTF_TILE(LAND_FIELD,NTILES)!IN Coefficient to calculate heat   
!                                  !    flux into land tiles (W/m2/K).  
     &,RESFT(LAND_FIELD,NTILES)    ! IN Total resistance factor         
     &,RHOKPM(LAND_FIELD,NTILES)   ! IN Surface exchange coeff for tiles
     &,RHOKPM_SICE(P_FIELD)        ! IN Sea-ice surface exchange coeff. 
     &,RHOKM_1_U_V(U_V_FIELD)      ! IN Level 1 exchange coefficient for
!                                       momentum                        
     &,RHOKH_1(LAND_FIELD,NTILES)  ! IN Surface exchange coeffs for FTL 
                                                                        
     &,RHOKH1_SICE(P_FIELD)        ! IN Sea and sea-ice surface exchange
     &,CT_CTQ_1(P_FIELD)           ! IN Coefficient in T and q          
!                                        tri-diagonal implicit matrix   
     &,CQ_CM_1(U_V_FIELD)          ! IN Coefficient in U and V          
!                                        tri-diagonal implicit matrix   
     &,DQW_1(P_FIELD)              ! IN Level 1 increment to q field    
     &,DTL_1(P_FIELD)              ! IN Level 1 increment to T field    
     &,DU_1(U_V_FIELD)             ! IN Level 1 increment to u wind     
!                                       field                           
     &,DV_1(U_V_FIELD)             ! IN Level 1 increment to v wind     
!                                       field                           
     &,FLANDG_UV(U_V_FIELD)        ! IN Land fraction on UV grid.       
                                                                       
                                                                        
      REAL                                                              
     & FQW_GB(P_FIELD)             ! INOUT Grid-box value of QW flux at 
!                                          Kg/sq m/s                    
     &,FTL_GB(P_FIELD)             ! INOUT Grid-box value of TL flux at 
!                                          i.e. H/Cp where H is sensible
!                                          in W per sq m).              
     &,TAUX_1(U_V_FIELD)           ! OUT   x-component of turbulent     
!                                          stress at surface.           
     &,TAUX_LAND(U_V_FIELD)        ! INOUT x-component of turbulent     
!                                  !       stress at land surface.      
     &,TAUX_SSI(U_V_FIELD)         ! INOUT x-component of turbulent     
!                                  !       stress at sea surface.       
     &,TAUY_1(U_V_FIELD)           ! OUT   y-component of turbulent     
!                                          stress at surface.           
     &,TAUY_LAND(U_V_FIELD)        ! INOUT y-component of turbulent     
!                                  !       stress at land surface.      
     &,TAUY_SSI(U_V_FIELD)         ! INOUT y-component of turbulent     
!                                  !       stress at sea surface.       
     &,FQW_TILE(LAND_FIELD,NTILES) ! INOUT Tile flux of QW. Kg/sq m/s   
     &,FTL_TILE(LAND_FIELD,NTILES) ! INOUT Tile flux of TL              
     &,E_SEA(P_FIELD)              ! INOUT Evaporation from sea times   
!                                          leads fraction (kg/m2/s).    
!                                          Zero over land.              
     &,H_SEA(P_FIELD)              ! INOUT Surface sensible heat flux ov
!                                          sea times leads fraction (W/m
!                                          Zero over land.              
                                                                        
                                                                        
                                                                        
!  External references :-                                               
      EXTERNAL TIMER                                                    
                                                                        
                                                                        
!  Local and other symbolic constants :-                                
*CALL C_LHEAT                                                           
*CALL C_R_CP                                                            
                                                                        
                                                                        
      REAL LS                                                           
      PARAMETER (                                                       
     & LS=LC+LF     ! Latent heat of sublimation (J per kg).            
     &)                                                                 
                                                                        
! Workspace :-                                                          
      REAL                                                              
     & FQW_ICE(P_FIELD)            ! "Explicit" surface flux of QW for  
!                                     sea-ice fraction of gridsquare.   
     &,FTL_ICE(P_FIELD)            ! "Explicit" surface flux of TL for  
!                                     sea-ice fraction of gridsquare.   
     &,LAT_HT   ! Latent heat of evaporation for snow-free land         
!               ! or sublimation for snow-covered land and ice.         
     &,APART(P_FIELD,2)            ! Tempary array                      
     &,BPART(P_FIELD,2)            ! Tempary array                      
     &,RECIP(P_FIELD)              ! Tempary array                      
     &,FTL_LAND(P_FIELD)           ! Tempary array                   
     &,FQW_LAND(P_FIELD)           ! Tempary array                   
                                                                        
!  Local scalars :-                                                     
      INTEGER                                                           
     & I        ! Loop counter (horizontal field index).                
     &,J        ! Loop counter (tile index).                            
     &,L        ! Loop counter (horizontal land index).                 
     &,N        ! Loop counter (tile counter).                          
                                                                        
      REAL                                                              
     & FTL_OLD  ! Used to hold current value of FTL_GB before updating  
                                                                        
                                                                        
!-----------------------------------------------------------------------
!!  0.  Check that the scalars input to define the grid are consistent. 
!       See comments to routine SF_EXCH for details.                    
!-----------------------------------------------------------------------
                                                                        
      IF (LTIMER) THEN                                                  
        CALL TIMER('IM_SF_PT ',3)                                       
      ENDIF                                                             
                                                                        
! Initialise APART and BPART to zero                                    
      DO I=P1,P1+P_POINTS-1                                             
        APART(I,1)=0.0                                                  
        APART(I,2)=0.0                                                  
        BPART(I,1)=0.0                                                  
        BPART(I,2)=0.0                                                  
        FTL_LAND(I)=0.0                                               
        FQW_LAND(I)=0.0                                               
      ENDDO                                                             
                                                                        
                                                                        
! Land tiles                                                            
      DO N=1,NTILES                                                     
        DO J=1,TILE_PTS(N)                                              
          L = TILE_INDEX(J,N)                                           
          I = LAND_INDEX(L)                                             
          LAT_HT = LC                                                   
          IF (SNOW_TILE(L,N).GT.0.) LAT_HT = LS                         
                                                                        
          APART(I,1)=APART(I,1) - TILE_FRAC(L,N) *                      
     &               GAMMA_1 * RHOKPM(L,N) *                            
     &            ( LAT_HT*RESFT(L,N)*RHOKH_1(L,N)*ALPHA1(L,N) +        
     &                         ASHTF_TILE(L,N) )                        
          APART(I,2)=APART(I,2) + TILE_FRAC(L,N) *                      
     &               GAMMA_1 * RHOKPM(L,N) *                            
     &               LAT_HT*RESFT(L,N)*RHOKH_1(L,N)                     
          BPART(I,1)=BPART(I,1) + TILE_FRAC(L,N) *                      
     &               GAMMA_1 * RESFT(L,N)*RHOKPM(L,N) *                 
     &               CP*RHOKH_1(L,N)*ALPHA1(L,N)                        
          BPART(I,2)=BPART(I,2) - TILE_FRAC(L,N) *                      
     &               GAMMA_1 * RESFT(L,N)*RHOKPM(L,N) *                 
     &               ( CP*RHOKH_1(L,N) + ASHTF_TILE(L,N) )              
                                                                        
        ENDDO                                                           
      ENDDO                                                             
                                                                        
                                                                        
                                                                        
! Sea points                                                            
      DO I=P1,P1+P_POINTS-1                                             
                                                                        
        IF(FLANDG(I).LT.1.0 .AND. ICE_FRACT(I).GT.0.0) THEN         
! Sea ice point                                                         
          APART(I,1)=FLANDG(I)*APART(I,1)                         
     &       - GAMMA_1 * (1.0-FLANDG(I)) * ICE_FRACT(I)             
     &      * RHOKPM_SICE(I) *                                        
     &      ( LS*RHOKH1_SICE(I)*ALPHA1_SICE(I) + ASHTF(I) )       
     &       - GAMMA_1 * (1.0-FLANDG(I)) * ( 1.0 - ICE_FRACT(I) )   
     &      * RHOKH1_SICE(I)                                          
                                                                        
          APART(I,2)=FLANDG(I)*APART(I,2)                         
     &       + GAMMA_1 * (1.0-FLANDG(I)) *ICE_FRACT(I)              
     &       * RHOKPM_SICE(I) * LS*RHOKH1_SICE(I)                   
                                                                        
          BPART(I,1)=FLANDG(I)*BPART(I,1)                         
     &       + GAMMA_1 * ICE_FRACT(I) * ( 1.0 - FLANDG(I) )         
     &       * RHOKPM_SICE(I) *CP*RHOKH1_SICE(I)*ALPHA1_SICE(I)   
                                                                        
          BPART(I,2)=FLANDG(I)*BPART(I,2)                         
     &       - GAMMA_1 * ICE_FRACT(I) * ( 1.0 - FLANDG(I) )         
     &       * RHOKPM_SICE(I) * ( CP*RHOKH1_SICE(I) + ASHTF(I) )  
     &       - GAMMA_1 * ( 1.0 - ICE_FRACT(I) )                       
     &       * ( 1.0 - FLANDG(I) ) * RHOKH1_SICE(I)                 
                                                                        
        ELSEIF(FLANDG(I).LT.1.0 .AND. .NOT.ICE_FRACT(I).GT.0.0) THEN   
! Ordinary sea point                                                    
          APART(I,1)= FLANDG(I)*APART(I,1)                        
     &       - GAMMA_1 * ( 1.0 - FLANDG(I) ) * RHOKH1_SICE(I)       
          APART(I,2)= FLANDG(I)*APART(I,2)                        
                                                                        
          BPART(I,1)= FLANDG(I)*BPART(I,1)                        
          BPART(I,2)= FLANDG(I)*BPART(I,2)                        
     &       - GAMMA_1 * ( 1.0 - FLANDG(I) ) * RHOKH1_SICE(I)       
                                                                        
        ENDIF                                                           
      ENDDO                                                             
                                                                        
                                                                        
                                                                        
! Calculate grid-box fluxes of heat and moisture                        
      DO I=P1,P1+P_POINTS-1                                             
        RECIP(I)=( 1.0 + CT_CTQ_1(I)*APART(I,1) ) *                     
     &           ( 1.0 + CT_CTQ_1(I)*BPART(I,2) ) -                     
     &             CT_CTQ_1(I)*APART(I,2)*CT_CTQ_1(I)*BPART(I,1)        
                                                                        
        FTL_OLD=FTL_GB(I)                                               
                                                                        
        FTL_GB(I) = ( ( 1.0 + CT_CTQ_1(I)*BPART(I,2) ) * ( FTL_OLD +    
     &                APART(I,1)*DTL_1(I) + APART(I,2)*DQW_1(I)) -      
     &                  CT_CTQ_1(I)*APART(I,2) * ( FQW_GB(I) +          
     &                  BPART(I,1)*DTL_1(I) + BPART(I,2)*DQW_1(I)) ) /  
     &                  RECIP(I)                                        
                                                                        
        FQW_GB(I) = ( ( 1.0 + CT_CTQ_1(I)*APART(I,1) ) * ( FQW_GB(I) +  
     &                  BPART(I,1)*DTL_1(I) + BPART(I,2)*DQW_1(I)) -    
     &                  CT_CTQ_1(I)*BPART(I,1) * ( FTL_OLD +            
     &                APART(I,1)*DTL_1(I) + APART(I,2)*DQW_1(I)) ) /    
     &                  RECIP(I)                                        
                                                                        
      ENDDO                                                             
                                                                        
                                                                        
! Make implicit correction to tile fluxes                               
                                                                        
! Land tiles                                                            
      DO N=1,NTILES                                                     
        DO J=1,TILE_PTS(N)                                              
          L = TILE_INDEX(J,N)                                           
          I = LAND_INDEX(L)                                             
          LAT_HT = LC                                                   
          IF (SNOW_TILE(L,N).GT.0.) LAT_HT = LS                         
                                                                        
          FTL_TILE(L,N)=FTL_TILE(L,N) -                                 
     &               GAMMA_1 * RHOKPM(L,N) *                            
     &            ( LAT_HT*RESFT(L,N)*RHOKH_1(L,N)*ALPHA1(L,N) +        
     &                       ASHTF_TILE(L,N) ) *                        
     &         ( DTL_1(I) - CT_CTQ_1(I)*FTL_GB(I) ) +                   
     &               GAMMA_1 * RHOKPM(L,N) *                            
     &               LAT_HT*RESFT(L,N)*RHOKH_1(L,N) *                   
     &         ( DQW_1(I) - CT_CTQ_1(I)*FQW_GB(I) )                     
                                                                        
          FQW_TILE(L,N)=FQW_TILE(L,N) +                                 
     &               GAMMA_1 * RESFT(L,N)*RHOKPM(L,N) *                 
     &               CP*RHOKH_1(L,N)*ALPHA1(L,N) *                      
     &         ( DTL_1(I) - CT_CTQ_1(I)*FTL_GB(I) ) -                   
     &               GAMMA_1 * RESFT(L,N)*RHOKPM(L,N) *                 
     &               ( CP*RHOKH_1(L,N) + ASHTF_TILE(L,N) ) *            
     &         ( DQW_1(I) - CT_CTQ_1(I)*FQW_GB(I) )                     
                                                                        
          FQW_LAND(I)=FQW_LAND(I)+FQW_TILE(L,N)*TILE_FRAC(L,N)      
          FTL_LAND(I)=FTL_LAND(I)+FTL_TILE(L,N)*TILE_FRAC(L,N)      
                                                                        
        ENDDO                                                           
      ENDDO                                                             
                                                                        
                                                                        
                                                                        
! Sea points                                                            
      DO I=P1,P1+P_POINTS-1                                             
                                                                        
        IF(FLANDG(I).LT.1.0 .AND. ICE_FRACT(I).GT.0.0) THEN         
! Sea ice point                                                         
          H_SEA(I)=H_SEA(I) - GAMMA_1 * (1.0 - ICE_FRACT(I)) * CP *     
     &       RHOKH1_SICE(I) * ( DTL_1(I) - CT_CTQ_1(I) * FTL_GB(I) )    
          E_SEA(I)=E_SEA(I) - GAMMA_1 * (1.0 - ICE_FRACT(I)) *          
     &       RHOKH1_SICE(I) * ( DQW_1(I) - CT_CTQ_1(I) * FQW_GB(I) )    
          FTL_ICE(I)=(FTL_GB(I)                                     
     &       - FTL_LAND(I)*FLANDG(I))/(1.-FLANDG(I))              
     &       - H_SEA(I)/CP                                            
          FQW_ICE(I)=(FQW_GB(I)                                     
     &       - FQW_LAND(I)*FLANDG(I))/(1.-FLANDG(I))              
     &        - E_SEA(I)                                              
                                                                        
        ELSEIF(FLANDG(I).LT.1.0 .AND. .NOT.ICE_FRACT(I).GT.0.0) THEN    
! Ordinary sea point                                                    
          H_SEA(I)=CP * (FTL_GB(I)                                  
     &       - FTL_LAND(I)*FLANDG(I))/(1.-FLANDG(I))              
          E_SEA(I)=(FQW_GB(I)                                       
     &       - FQW_LAND(I)*FLANDG(I))/(1.-FLANDG(I))              
          FTL_ICE(I)=0.0                                                
          FQW_ICE(I)=0.0                                                
                                                                        
        ENDIF                                                           
      ENDDO                                                             
                                                                        
                                                                        
*IF -DEF,SCMA                                                           
      DO I=U1_V1+ROW_LENGTH,U1_V1+U_V_POINTS-ROW_LENGTH-1             
*ELSE                                                                   
      DO I=1,U_V_POINTS                                               
*ENDIF                                                                  
                                                                        
        IF(FLANDG_UV(I).GT.0.0)THEN                                    
          TAUX_LAND(I) = ( TAUX_LAND(I) +                           
     &                 GAMMA_1*RHOKM_1_U_V(I)*DU_1(I) ) /             
     &                ( 1.0 + GAMMA_1*RHOKM_1_U_V(I)*CQ_CM_1(I) )   
        ELSE                                                            
          TAUX_LAND(I) = 0.0                                          
        ENDIF                                                           
                                                                        
        IF(FLANDG_UV(I).LT.1.0)THEN                                    
          TAUX_SSI(I) = ( TAUX_SSI(I) +                             
     &                 GAMMA_1*RHOKM_1_U_V(I)*DU_1(I) ) /             
     &                ( 1.0 + GAMMA_1*RHOKM_1_U_V(I)*CQ_CM_1(I) )   
        ELSE                                                            
          TAUX_SSI(I) = 0.0                                           
        ENDIF                                                           
                                                                        
        TAUX_1(I) = FLANDG_UV(I)*TAUX_LAND(I)                      
     &                + ( 1.0-FLANDG_UV(I))*TAUX_SSI(I)       

        IF(FLANDG_UV(I).GT.0.0)THEN                                    
          TAUY_LAND(I) = ( TAUY_LAND(I) +                           
     &                 GAMMA_1*RHOKM_1_U_V(I)*DV_1(I) ) /             
     &                ( 1.0 + GAMMA_1*RHOKM_1_U_V(I)*CQ_CM_1(I) )   
        ELSE                                                            
          TAUY_LAND(I) = 0.0                                          
        ENDIF                                                           
                                                                        
        IF(FLANDG_UV(I).LT.1.0)THEN                                    
          TAUY_SSI(I) = ( TAUY_SSI(I) +                             
     &                 GAMMA_1*RHOKM_1_U_V(I)*DV_1(I) ) /             
     &                ( 1.0 + GAMMA_1*RHOKM_1_U_V(I)*CQ_CM_1(I) )   
        ELSE                                                            
          TAUY_SSI(I) = 0.0                                           
        ENDIF                                                           
                                                                        
        TAUY_1(I) = FLANDG_UV(I)*TAUY_LAND(I)                      
     &                + ( 1.0-FLANDG_UV(I))*TAUY_SSI(I)              

      ENDDO  !u_v_points                                                
                                                                        
      IF (LTIMER) THEN                                                  
        CALL TIMER('IM_SF_PT ',4)                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
*ENDIF                                                                  
*DECK SFEXPL8A                                                          
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
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
!!!  SUBROUTINE SF_EXPL------------------------------------------------ 
!!!                                                                     
!!!  Purpose: Calculate explicit surface fluxes of heat, moisture and   
!!!           momentum. Also calculates surface exchange coefficients   
!!!           required for implicit update of surface fluxes and surface
!!!           information required by the explicit boundary layer routin
!!!                                                                     
!!!                                                                     
!!! F.Hewer     <- programmer of some or all of previous code or changes
!!! C.Wilson    <- programmer of some or all of previous code or changes
!!!                                                                     
!!!  Model            Modification history:                             
!!! version  Date                                                       
!!!                                                                     
!!!                                                                     
!!!  Programming standard: Unified Model Documentation Paper No 4,      
!!!                        Version ?, dated ?.                          
!!!                                                                     
!!!  System component covered: P24.                                     
!!!                                                                     
!!!  Project task:                                                      
!!!                                                                     
!!!  Documentation: UMDP 24.                                            
!!!                                                                     
!!!---------------------------------------------------------------------
                                                                        
!    Arguments :-                                                       
      SUBROUTINE SF_EXPL (                                              
                                                                        
! IN values defining field dimensions and subset to be processed :      
     & P_FIELD,U_FIELD,LAND_FIELD,ROW_LENGTH,                           
     & N_P_ROWS,N_U_ROWS,P_POINTS,P1,LAND1,LAND_PTS,U_POINTS,U1,        
                                                                        
! IN values defining vertical grid of model atmosphere :                
     & AK_1,BK_1,AKH_1,BKH_1,DELTA_AK_1,DELTA_BK_1,                     
     & EXNER,                                                           
                                                                        
! IN soil/vegetation/land surface data :                                
     & LAND_INDEX,                                                      
     & LAND_MASK,L_Z0_OROG,                                             
     & NTILES,TILE_INDEX,TILE_PTS,SM_LEVELS,                            
     & CANHC_TILE,CANOPY,CATCH,FLAKE,GC,HCON,HO2R2_OROG,                
     & FLAND,FLANDG,                                                    
     & SNOW_TILE,SIL_OROG_LAND,SMVCST,STHF,STHU,                        
     & TILE_FRAC,VFRAC_TILE,Z0_TILE,                                    
                                                                        
! IN sea/sea-ice data :                                                 
     & ICE_FRACT,U_0,V_0,                                               
                                                                        
! IN cloud data :                                                       
     & CF_1,QCF_1,QCL_1,                                                
                                                                        
! IN everything not covered so far :                                    
     & PSTAR,LW_DOWN,RAD_SICE,SW_TILE,TIMESTEP,
     & VSHR,VSHR_LAND,VSHR_SSI,ZH,                 
     & Q_1,T_1,T_SOIL,TI,                                               
     & TSTAR,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,
     & TSTAR_TILE,U_1,V_1,                                        
     & L_BL_LSPICE,                                                     
                                                                        
! IN STASH flags :-                                                     
     & SFME,SQ1P5,ST1P5,SU10,SV10,                                      
                                                                        
! INOUT data :                                                          
     & Z0MSEA,                                                          
                                                                        
! OUT Diagnostic not requiring STASH flags :                            
     & CD,CH,E_SEA,QW_1,TL_1,FQW_1,                                     
     & FTL_1,FTL_TILE,LE_TILE,H_SEA,RADNET_SICE,RADNET_TILE,            
     & RHOKM_UV_1,RIB,RIB_TILE,TAUX_1,TAUY_1,                           
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           
                                                                        
! OUT diagnostic requiring STASH flags :                                
     & FME,                                                             
                                                                        
! OUT data required for tracer mixing :                                 
     & RHO_ARESIST,ARESIST,RESIST_B,                                    
     & RHO_ARESIST_TILE,ARESIST_TILE,RESIST_B_TILE,                     
     & NRML,                                                            
                                                                        
! OUT data required for 4D-VAR :                                        
     & RHO_CD_MODV1,                                                    
                                                                        
! OUT data required elsewhere in UM system :                            
     & FB_SURF,U_S,T1_SD,Q1_SD,TV1_SD,                                  
                                                                        
! OUT data required elsewhere in boundary layer or surface code         
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,DTRDZ_1,FQW_TILE,            
     & FQW_ICE,FTL_ICE,TSTAR_TILE_OLD,FRACA,RESFS,RESFT,                
     & RHOKH,RHOKH_TILE,RHOKH_SICE,RHOKPM,RHOKPM_SICE,                  
     & Z1,H_BLEND_OROG,Z0H,Z0H_TILE,Z0M,Z0M_TILE,Z0M_EFF,               
     & CDR10M_UV,CHR1P5M,CHR1P5M_SICE,                                  
     & FLANDG_UV,                                               
                                                                        
! LOGICAL LTIMER                                                        
     & LTIMER                                                           
     & )                                                                
                                                                        
      IMPLICIT NONE                                                     
                                                                        
!  Inputs :-                                                            
                                                                        
! (a) Defining horizontal grid and subset thereof to be processed.      
!    Checked for consistency.                                           
                                                                        
      INTEGER                                                           
     & P_FIELD                     ! IN No. of P-points in whole grid   
!                                  !    (for dimensioning only).        
     &,U_FIELD                     ! IN No. of UV-points in whole grid. 
     &,LAND_FIELD                  ! IN No.of land points in whole grid.
     &,ROW_LENGTH                  ! IN No. of points in one row.       
     &,N_P_ROWS   ! IN No of P-rows being processed.                    
     &,N_U_ROWS   ! IN No of UV-rows being processed.                   
     &,P_POINTS   ! IN No of P-points being processed.                  
     &,P1         ! IN First P-point to be processed.                   
     &,LAND1      ! IN First land-point to be processed.                
!                 !       1 <= LAND1 <= LAND_FIELD                      
     &,LAND_PTS   ! IN No of land points being processed.               
     &,U_POINTS   ! IN No of UV-points being processed.                 
     &,U1         ! IN First UV-point to be processed.                  
                                                                        
! (b) Defining vertical grid of model atmosphere.                       
                                                                        
      REAL                                                              
     & AK_1                        ! IN Hybrid 'A' for all levels.      
     &,BK_1                        ! IN Hybrid 'B' for all levels.      
     &,AKH_1                       ! IN Hybrid 'A' for layer interfaces.
     &,BKH_1                       ! IN Hybrid 'B' for layer interfaces.
     &,EXNER(P_FIELD,2)            ! IN Exner function.  EXNER(,K) is   
!                                  !    value for LOWER BOUNDARY of     
!                                  !    level K.                        
     &,DELTA_AK_1               ! IN Difference of hybrid 'A' across    
!                               !    layers (K-1/2 to K+1/2).           
!                               !    NB: Upper minus lower.             
     &,DELTA_BK_1               ! IN Difference of hybrid 'B' across    
!                               !     layers (K-1/2 to K+1/2).          
!                               !     NB: Upper minus lower.            
                                                                        
! (c) Soil/vegetation/land surface parameters (mostly constant).        
                                                                        
      LOGICAL                                                           
     & LAND_MASK(P_FIELD)          ! IN T if land, F elsewhere.         
     &,L_Z0_OROG                   ! IN T to use orog.roughness         
!                                  !    treatment in SF_EXCH            
                                                                        
      INTEGER                                                           
     & LAND_INDEX(P_FIELD)         ! IN LAND_INDEX(I)=J => the Jth      
!                                  !    point in P_FIELD is the Ith     
!                                  !    land point.                     
                                                                        
      INTEGER                                                           
     & SM_LEVELS                   ! IN No. of soil moisture levels     
     &,NTILES                      ! IN No. of land-surface tiles       
     &,TILE_INDEX(LAND_FIELD,NTILES)!IN Index of tile points            
     &,TILE_PTS(NTILES)            ! IN Number of tile points           
                                                                        
      REAL                                                              
     & CANHC_TILE(LAND_FIELD,NTILES)!IN Areal heat capacity of canopy   
!                                  !    for land tiles (J/K/m2).        
     &,CANOPY(LAND_FIELD,NTILES)   ! IN Surface/canopy water for        
!                                  !    snow-free land tiles (kg/m2)    
     &,CATCH(LAND_FIELD,NTILES)    ! IN Surface/canopy water capacity   
!                                  !    of snow-free land tiles (kg/m2).
     &,FLAKE(LAND_FIELD,NTILES)    ! IN Lake fraction.                  
     &,GC(LAND_FIELD,NTILES)       ! IN "Stomatal" conductance to       
!                                  !     evaporation for land tiles     
!                                  !     (m/s).                         
     &,HCON(LAND_FIELD)            ! IN Soil thermal conductivity       
!                                  !    (W/m/K).                        
     &,SNOW_TILE(LAND_FIELD,NTILES)! IN Lying snow on tiles (kg/m2)     
     &,SMVCST(LAND_FIELD)          ! IN Volumetric saturation point     
!                                  !    (m3/m3 of soil).                
     &,STHF(LAND_FIELD,SM_LEVELS)  ! IN Frozen soil moisture content of 
!                                  !    each layer as a fraction of     
!                                  !    saturation.                     
     &,STHU(LAND_FIELD,SM_LEVELS)  ! IN Unfrozen soil moisture content  
!                                  !    of each layer as a fraction of  
!                                  !    saturation.                     
     &,TILE_FRAC(LAND_FIELD,NTILES)! IN Tile fractions including        
!                                  ! snow cover in the ice tile.        
     &,VFRAC_TILE(LAND_FIELD,NTILES)!IN Fractional canopy coverage for  
!                                  !    land tiles.                     
     &,Z0_TILE(LAND_FIELD,NTILES)  ! IN Tile roughness lengths (m).     
     &,SIL_OROG_LAND(LAND_FIELD)   ! IN Silhouette area of unresolved   
!                                  !    orography per unit horizontal   
!                                  !    area on land points only.       
     &,HO2R2_OROG(LAND_FIELD)      ! IN Standard Deviation of orography.
!                                  !    equivilent to peak to trough    
!                                  !    height of unresolved orography  
!                                  !    divided by 2SQRT(2) on land     
!                                  !    points only (m)                 
     &,FLAND(LAND_FIELD)           ! IN Land fraction on land tiles.    
     &,FLANDG(P_FIELD)             ! IN Land fraction on all tiles.     
                                                                        
! (d) Sea/sea-ice data.                                                 
                                                                        
      REAL                                                              
     & ICE_FRACT(P_FIELD)          ! IN Fraction of gridbox covered by  
!                                  !     sea-ice (decimal fraction).    
     &,U_0(U_FIELD)                ! IN W'ly component of surface       
!                                  !    current (m/s).                  
     &,V_0(U_FIELD)                ! IN S'ly component of surface       
!                                  !    current (m/s).                  
                                                                        
! (e) Cloud data.                                                       
                                                                        
      REAL                                                              
     & CF_1(P_FIELD)               ! IN Cloud fraction (decimal).       
     &,QCF_1(P_FIELD)              ! IN Cloud ice (kg per kg air)       
     &,QCL_1(P_FIELD)              ! IN Cloud liquid water (kg          
!                                  !    per kg air).                    
                                                                        
! (f) Atmospheric + any other data not covered so far, incl control.    
                                                                        
      REAL                                                              
     & PSTAR(P_FIELD)              ! IN Surface pressure (Pascals).     
     &,LW_DOWN(P_FIELD)            ! IN Surface downward LW radiation   
!                                  !    (W/m2).                         
     &,RAD_SICE(P_FIELD)           ! IN Surface net shortwave and       
!                                  !    downward LWradiation for        
!                                  !    sea-ice (W/sq m).               
     &,SW_TILE(LAND_FIELD,NTILES)  ! IN Surface net SW radiation on     
!                                  !    land tiles (W/m2).              
     &,TIMESTEP                    ! IN Timestep (seconds).             
     &,VSHR(P_FIELD)               ! IN Magnitude of surface-to-lowest  
!                                  !    atm level wind shear (m per s). 
     &,VSHR_LAND(P_FIELD)          ! IN Magnitude of surface-to-lowest 
!                                  !    atm level wind shear (m per s).
     &,VSHR_SSI(P_FIELD)           ! IN Magnitude of surface-to-lowest 
!                                  !    atm level wind shear (m per s).
     &,ZH(P_FIELD)                 ! IN Height above surface of top of  
!                                  !    boundary layer (metres).        
     &,Q_1(P_FIELD)                ! IN  Specific humidity ( kg/kg air).
     &,T_1(P_FIELD)                ! IN  Atmospheric temperature (K).   
     &,T_SOIL(LAND_FIELD,SM_LEVELS)! IN Soil temperatures (K).          
     &,TI(P_FIELD)                 ! IN Sea-ice surface layer           
!                                  !    temperature (K).                
     &,TSTAR(P_FIELD)              ! IN GBM surface temperature (K).    
     &,TSTAR_LAND(P_FIELD)         ! IN Land mean surface temperature (K
     &,TSTAR_SEA(P_FIELD)          ! IN Open sea surface temperature (K)
     &,TSTAR_SICE(P_FIELD)         ! IN Sea-ice surface temperature (K).
     &,TSTAR_SSI(P_FIELD)          ! IN mean sea surface temperature (K)
     &,TSTAR_TILE(LAND_FIELD,NTILES)!IN Surface tile temperatures       
     &,U_1(U_FIELD)                ! IN W'ly wind component (m/s)       
     &,V_1(U_FIELD)                ! IN S'ly wind component (m/s)       
                                                                        
      LOGICAL                                                           
     & LTIMER                      ! IN Logical switch for TIMER diags  
     &,L_BL_LSPICE                 ! IN Use if 3A large scale precip    
                                                                        
!  STASH flags :-                                                       
                                                                        
      LOGICAL                                                           
     & SFME    ! IN Flag for FME (q.v.).                                
     &,SQ1P5   ! IN Flag for Q1P5M (q.v.)                               
     &,ST1P5   ! IN Flag for T1P5M (q.v.)                               
     &,SU10    ! IN Flag for U10M (q.v.)                                
     &,SV10    ! IN Flag for V10M (q.v.)                                
                                                                        
!  In/outs :-                                                           
                                                                        
      REAL                                                              
     & Z0MSEA(P_FIELD)             ! INOUT Sea-surface roughness        
!                                  !       length for momentum (m).     
                                                                        
!  Outputs :-                                                           
!-1 Diagnostic (or effectively so - includes coupled model requisites):-
                                                                        
!  (a) Calculated anyway (use STASH space from higher level) :-         
!                                                                       
      REAL                                                              
     & CD(P_FIELD)                 ! OUT Turbulent surface exchange     
!                                  !     (bulk transfer) coefficient for
!                                  !     momentum.                      
     &,CH(P_FIELD)                 ! OUT Turbulent surface exchange     
!                                  !     (bulk transfer) coefficient for
!                                  !     heat and/or moisture.          
     &,E_SEA(P_FIELD)              ! OUT Evaporation from sea times     
!                                  !     leads fraction. Zero over land.
!                                  !     (kg per square metre per sec). 
     &,QW_1(P_FIELD)               ! OUT Total water content            
     &,TL_1(P_FIELD)               ! OUT Ice/liquid water temperature   
     &,FQW_1(P_FIELD)              ! OUT Moisture flux between layers   
!                                  !     (kg per square metre per sec). 
!                                  !     FQW(,1) is total water flux    
!                                  !     from surface, 'E'.             
     &,FTL_1(P_FIELD)              ! OUT FTL(,K) contains net turbulent 
!                                  !     sensible heat flux into layer K
!                                  !     from below; so FTL(,1) is the  
!                                  !     surface sensible heat, H.(W/m2)
     &,FTL_TILE(LAND_FIELD,NTILES) ! OUT Surface FTL for land tiles     
     &,LE_TILE(LAND_FIELD,NTILES)  ! OUT Surface latent heat flux for   
!                                  !     land tiles                     
     &,H_SEA(P_FIELD)              ! OUT Surface sensible heat flux over
!                                  !     sea times leads fraction (W/m2)
     &,RADNET_SICE(P_FIELD)        ! OUT Surface net radiation on       
!                                  !     sea-ice (W/m2)                 
     &,RADNET_TILE(LAND_FIELD,NTILES)                                   
!                                  ! OUT Surface net radiation on       
!                                  !     land tiles (W/m2)              
     &,RHOKM_LAND(P_FIELD)         ! OUT Exchange coefficients for      
!                                  !     momentum on P-grid             
     &,RHOKM_SSI(P_FIELD)          ! OUT Exchange coefficients for      
!                                  !     momentum on P-grid             
     &,RHOKM_UV_LAND(U_FIELD)      ! OUT Exchange coefficients for   
!                                  !     momentum (on UV-grid, with 1st 
!                                  !     and last rows undefined or, at 
!                                  !     present, set to "missing data")
     &,RHOKM_UV_SSI(U_FIELD)       ! OUT Exchange coefficients for    
!                                  !     momentum (on UV-grid, with 1st 
!                                  !     and last rows undefined or, at 
!                                  !     present, set to "missing data")
     &,FLANDG_UV(U_FIELD)          ! OUT Land frac (on UV-grid, with 1st
!                                  !     and last rows undefined or, at 
!                                  !     present, set to "missing data")
     &,RHOKM_UV_1(U_FIELD)         ! OUT Exchange coefficients for      
!                                  !     momentum (on UV-grid, with 1st 
!                                  !     and last rows undefined or, at 
!                                  !     present, set to "missing data")
     &,RIB(P_FIELD)                ! OUT Mean bulk Richardson number for
!                                  !     lowest layer.                  
     &,RIB_TILE(LAND_FIELD,NTILES) ! OUT RIB for land tiles.            
     &,TAUX_LAND(U_FIELD)          ! OUT W'ly component of sfc wind     
!                                  !     stress (N/sq m). (On UV-grid   
!                                  !     with first and last rows       
!                                  !     undefined or, at present,      
!                                  !     set to missing data            
     &,TAUX_SSI(U_FIELD)           ! OUT W'ly component of sfc wind     
!                                  !     stress (N/sq m). (On UV-grid   
!                                  !     with first and last rows       
!                                  !     undefined or, at present,      
!                                  !     set to missing data            
     &,TAUX_1(U_FIELD)             ! OUT W'ly component of surface wind 
!                                  !     stress (N/sq m). (On UV-grid   
!                                  !     with first and last rows       
!                                  !     undefined or, at present,      
!                                  !     set to missing data            
     &,TAUY_LAND(U_FIELD)          ! OUT S'ly component of sfc wind     
!                                  !     stress (N/sq m).  On UV-grid;  
!                                  !     comments as per TAUY.          
     &,TAUY_SSI(U_FIELD)           ! OUT S'ly component of sfc wind     
!                                  !     stress (N/sq m).  On UV-grid;  
!                                  !     comments as per TAUY.          
     &,TAUY_1(U_FIELD)             ! OUT S'ly component of surface wind 
!                                  !     stress (N/sq m).  On UV-grid;  
!                                  !     comments as per TAUX.          
     &,RHO_CD_MODV1(P_FIELD)       ! OUT Surface air density * drag coef
!                                  !     *mod(v1 - v0) before interp    
     &,RHO_ARESIST(P_FIELD)        ! OUT RHOSTAR*CD_STD*VSHR for Sulphur
!                                  !     cycle                          
     &,ARESIST(P_FIELD)            ! OUT 1/(CD_STD*VSHR) for Sulphur    
!                                  !     cycle                          
     &,RESIST_B(P_FIELD)           ! OUT (1/CH-1/(CD_STD)/VSHR for      
!                                  !     Sulphur cycle                  
     &,RHO_ARESIST_TILE(LAND_FIELD,NTILES)                              
!                                  ! OUT RHOSTAR*CD_STD*VSHR on land    
!                                  !     tiles                          
     &,ARESIST_TILE(LAND_FIELD,NTILES)                                  
!                                  ! OUT 1/(CD_STD*VSHR) on land tiles  
     &,RESIST_B_TILE(LAND_FIELD,NTILES)                                 
!                                  ! OUT (1/CH-1/CD_STD)/VSHR on land   
!                                  !     tiles                          
                                                                        
      INTEGER                                                           
     & NRML(P_FIELD)               ! OUT Number of model layers in the  
!                                  !     Rapidly Mixing Layer; set to   
!                                  !     zero in SF_EXCH for MOSES II.  
                                                                        
!  (b) Not passed between lower-level routines (not in workspace at this
!      level) :-                                                        
                                                                        
      REAL                                                              
     & FME(P_FIELD)                ! OUT Wind mixing "power" (W/m2).    
                                                                        
!-2 Genuinely output, needed by other atmospheric routines :-           
                                                                        
      REAL                                                              
     & FB_SURF(P_FIELD)            ! OUT Surface flux buoyancy over     
!                                  !     density (m^2/s^3)              
     &,U_S(P_FIELD)                ! OUT Surface friction velocity (m/s)
     &,T1_SD(P_FIELD)              ! OUT Standard deviation of turbulent
!                                  !     fluctuations of layer 1 temp;  
!                                  !     used in initiating convection. 
     &,Q1_SD(P_FIELD)              ! OUT Standard deviation of turbulent
!                                  !     flucs of layer 1 humidity;     
!                                  !     used in initiating convection. 
     &,TV1_SD(P_FIELD)             ! OUT Standard deviation of turbulent
!                                  !     fluctuations of surface layer  
!                                  !     virtual temperature (K).       
                                                                        
      REAL                                                              
     & ALPHA1(LAND_FIELD,NTILES)   ! OUT Mean gradient of saturated     
!                                  !     specific humidity with respect 
!                                  !     to temperature between the     
!                                  !     bottom model layer and tile    
!                                  !     surfaces                       
     &,ALPHA1_SICE(P_FIELD)        ! OUT ALPHA1 for sea-ice.            
     &,ASHTF(P_FIELD)              ! OUT Coefficient to calculate       
!                                  !     surface heat flux into soil or 
!                                  !     sea-ice.                       
     &,ASHTF_TILE(LAND_FIELD,NTILES)!OUT Coefficient to calculate       
!                                  !     surface heat flux into land    
!                                  !     tiles.                         
     &,DTRDZ_1(P_FIELD)            ! OUT -g.dt/dp for model layers.     
     &,FQW_TILE(LAND_FIELD,NTILES) ! OUT Surface FQW for land tiles     
     &,FQW_ICE(P_FIELD)            ! OUT Surface FQW for sea-ice        
     &,FTL_ICE(P_FIELD)            ! OUT Surface FTL for sea-ice        
     &,TSTAR_TILE_OLD(LAND_FIELD,NTILES)                                
!                                  ! OUT Tile surface temperatures at   
!                                  !     beginning of timestep.         
     &,FRACA(LAND_FIELD,NTILES)    ! OUT Fraction of surface moisture   
!                                  !     flux with only aerodynamic     
!                                  !     resistance for snow-free land  
!                                  !     tiles.                         
     &,RESFS(LAND_FIELD,NTILES)    ! OUT Combined soil, stomatal        
!                                  !     and aerodynamic resistance     
!                                  !     factor for fraction (1-FRACA)  
!                                  !     of snow-free land tiles.       
     &,RESFT(LAND_FIELD,NTILES)    ! OUT Total resistance factor.       
!                                  !     FRACA+(1-FRACA)*RESFS for      
!                                  !     snow-free land, 1 for snow.    
     &,RHOKH(P_FIELD)              ! OUT Grid-box surface exchange      
!                                  !     coefficients                   
     &,RHOKH_TILE(LAND_FIELD,NTILES)                                    
!                                  ! OUT Surface exchange coefficients  
!                                  !     for land tiles                 
     &,RHOKH_SICE(P_FIELD)         ! OUT Surface exchange coefficients  
!                                  !     for sea and sea-ice            
     &,RHOKPM(LAND_FIELD,NTILES)   ! OUT Land surface exchange coeff.   
     &,RHOKPM_SICE(P_FIELD)        ! OUT Sea-ice surface exchange coeff.
     &,Z1(P_FIELD)                 ! OUT Height of lowest level (i.e.   
!                                  !     height of middle of lowest     
!                                  !     layer).                        
     &,H_BLEND_OROG(P_FIELD)       ! OUT Blending height used as part of
!                                  !     effective roughness scheme     
     &,Z0H(P_FIELD)                ! OUT Roughness length for heat and  
!                                  !     moisture (m).                  
     &,Z0H_TILE(LAND_FIELD,NTILES)                                      
!                                  ! OUT Tile roughness lengths for heat
!                                  !     and moisture (m).              
     &,Z0M(P_FIELD)                ! OUT Roughness length for           
!                                  !     momentum (m).                  
     &,Z0M_TILE(LAND_FIELD,NTILES)                                      
!                                  ! OUT Tile roughness lengths for     
!                                  !     momentum.                      
     &,Z0M_EFF(P_FIELD)            ! OUT Effective grid-box roughness   
!                                  !     length for momentum            
     &,CDR10M_UV(U_FIELD)          ! OUT Ratio of CD's reqd for         
!                                  !     calculation of 10 m wind. On   
!                                  !     UV-grid; comments as per RHOKM.
     &,CHR1P5M(LAND_FIELD,NTILES)  ! OUT Ratio of coefffs for           
!                                  !     calculation of 1.5m temp for   
!                                  !     land tiles.                    
     &,CHR1P5M_SICE(P_FIELD)       ! OUT CHR1P5M for sea and sea-ice    
!                                  !     (leads ignored).               
                                                                        
                                                                        
!---------------------------------------------------------------------  
!  External routines called :-                                          
                                                                        
      EXTERNAL Z,HEAT_CON,SF_EXCH,BOUY_TQ,BTQ_INT,                      
     & KMKH,EX_FLUX_TQ,EX_FLUX_UV,IM_CAL_TQ,SICE_HTF,SF_EVAP,SF_MELT,   
     & IM_CAL_UV,SCREEN_TQ                                              
      EXTERNAL TIMER                                                    
*IF -DEF,SCMA                                                           
      EXTERNAL UV_TO_P,P_TO_UV                                          
*ENDIF                                                                  
                                                                        
!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-        
                                                                        
*CALL C_R_CP                                                            
*CALL C_G                                                               
*CALL C_LHEAT                                                           
*CALL C_GAMMA                                                           
*CALL CSIGMA                                                            
*CALL SOIL_THICK                                                        
*IF DEF,MPP                                                             
! MPP Common block                                                      
*CALL PARVARS                                                           
*ENDIF                                                                  
                                                                        
! Derived local parameters.                                             
                                                                        
      REAL LCRCP,LS,LSRCP                                               
                                                                        
      PARAMETER (                                                       
     & LCRCP=LC/CP           ! Evaporation-to-dT conversion factor.     
     &,LS=LF+LC              ! Latent heat of sublimation.              
     &,LSRCP=LS/CP           ! Sublimation-to-dT conversion factor.     
     &  )                                                               
                                                                        
!-----------------------------------------------------------------------
                                                                        
!  Workspace :-                                                         
                                                                        
      REAL                                                              
     & BF_1(P_FIELD)            ! A buoyancy parameter (beta F tilde)   
     &,BQ_1(P_FIELD)            ! A buoyancy parameter (beta q tilde).  
     &,BT_1(P_FIELD)            ! A buoyancy parameter (beta T tilde).  
     &,BT_CLD_1(P_FIELD)        ! A buoyancy parameter for cloudy air   
     &,BQ_CLD_1(P_FIELD)        ! A buoyancy parameter for cloudy air   
     &,A_QS_1(P_FIELD)          ! Saturated lapse rate factor           
     &,A_DQSDT_1(P_FIELD)       ! Saturated lapse rate factor           
     &,DQSDT_1(P_FIELD)         ! Derivative of q_SAT w.r.t. T          
     &,DZL_1(P_FIELD)           ! DZL(,K) is depth in m of layer        
!                               ! K, i.e. distance from boundary        
!                               ! K-1/2 to boundary K+1/2.              
     &,HCONS(LAND_FIELD)        ! Soil thermal conductivity including   
!                               ! the effects of water and ice (W/m2)   
     &,P_1(P_FIELD)             ! Pressure at model levels              
     &,RDZ_1(P_FIELD)           ! RDZ(,1) is the reciprocal of the      
!                               ! height of level 1, i.e. of the        
!                               ! middle of layer 1.  For K > 1,        
!                               ! RDZ(,K) is the reciprocal             
!                               ! of the vertical distance              
!                               ! from level K-1 to level K.            
     &,RHOKM_1(P_FIELD)         ! Exchange coefficients for             
!                               ! momentum on P-grid                    
     &,DELTAP_1(P_FIELD)        ! Difference in pressure between levels 
     &,TV_1(P_FIELD)            ! Virtual temp                          
     &,ZLB(P_FIELD,0:1)         ! ZLB(,K) is the height of the          
!                               ! upper boundary of layer K             
!                               ! ( = 0.0 for "K=0").                   
     &,CDR10M(P_FIELD)          ! Ratio of CD's reqd for calculation    
!                               ! of 10 m wind. On P-grid               
     &,FLANDG_TMP(P_FIELD)      ! Land fraction on P-grid               
                                                                        
                                                                        
!  Local scalars :-                                                     
                                                                        
      INTEGER                                                           
     & I,J,L      ! LOCAL Loop counter (horizontal field index).        
     &,N          ! LOCAL Loop counter (tile index).                    
                                                                        
      IF (LTIMER) THEN                                                  
        CALL TIMER('BDYLAYR ',3)                                        
      ENDIF                                                             
                                                                       
!-----------------------------------------------------------------------
!! 1.  Perform calculations in what the documentation describes as      
!!     subroutine Z_DZ.  In fact, a separate subroutine isn't used.     
!-----------------------------------------------------------------------
                                                                        
!-----------------------------------------------------------------------
!! 1.1 Initialise ZLB(,0) (to zero, of course, this being the height    
!!     of the surface above the surface).                               
!-----------------------------------------------------------------------
                                                                        
      DO I=P1,P1+P_POINTS-1                                             
        ZLB(I,0)=0.0                                                    
      ENDDO                                                             
                                                                        
!-----------------------------------------------------------------------
!! 1.2 Calculate layer depths and heights, and construct wind fields on 
!!     P-grid.  This involves calling subroutines Z and UV_TO_P.        
!!     Virtual temperature is also calculated, as a by-product.         
!-----------------------------------------------------------------------
!  NB RDZ  TEMPORARILY used to return DELTA_Z_LOWER, the lower half     
!     layer thickness                                                   
                                                                        
      CALL Z(P_POINTS,EXNER(P1,1),EXNER(P1,2),PSTAR(P1),                
     &  AKH_1,BKH_1,Q_1(P1),QCF_1(P1),                                  
     &  QCL_1(P1),T_1(P1),ZLB(P1,0),TV_1(P1),                           
     &  ZLB(P1,1),DZL_1(P1),RDZ_1(P1),LTIMER)                           
                                                                        
                                                                        
! set pressure array.                                                   
      DO I=P1,P1+P_POINTS-1                                             
        P_1(I) = AK_1 + BK_1*PSTAR(I)                                   
      ENDDO                                                             
                                                                        
      DO I=P1,P1+P_POINTS-1                                             
        Z1(I)=RDZ_1(I)                                                  
        RDZ_1(I)=1.0/RDZ_1(I)                                           
        DELTAP_1(I)=DELTA_AK_1 + PSTAR(I)*DELTA_BK_1                    
        DTRDZ_1(I) = -G * TIMESTEP/DELTAP_1(I)                          
!     &                  (DELTA_AK_1 + PSTAR(I)*DELTA_BK_1)             
      ENDDO                                                             
                                                                        
      IF (LAND_FIELD.GT.0) THEN    ! Omit if no land points             
                                                                        
!-----------------------------------------------------------------------
! Calculate the thermal conductivity of the top soil layer.             
!-----------------------------------------------------------------------
        CALL HEAT_CON (LAND_FIELD,HCON,STHU,STHF,SMVCST,HCONS,LTIMER)   
                                                                        
      ENDIF                     ! End test on land points               
                                                                        
!-----------------------------------------------------------------------
!! Calculate total water content, QW and Liquid water temperature, TL   
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1                                             
        QW_1(I) = Q_1(I) + QCL_1(I) + QCF_1(I)              ! P243.10   
        TL_1(I) = T_1(I) - LCRCP*QCL_1(I) - LSRCP*QCF_1(I)  ! P243.9    
      ENDDO                                                             
                                                                        
!-----------------------------------------------------------------------
!! Calculate net radiation on land tiles and sea-ice                    
!-----------------------------------------------------------------------
      DO N=1,NTILES                                                     
        DO L=1,LAND_FIELD                                               
          RADNET_TILE(L,N) = 0.                                         
          LE_TILE(L,N) = 0.                                             
        ENDDO                                                           
      ENDDO                                                             
      DO N=1,NTILES                                                     
        DO J=1,TILE_PTS(N)                                              
          L = TILE_INDEX(J,N)                                           
          I = LAND_INDEX(L)                                             
          RADNET_TILE(L,N) = SW_TILE(L,N) +                             
     &                       LW_DOWN(I) - SBCON*T_SOIL(L,1)**4          
          TSTAR_TILE_OLD(L,N) = TSTAR_TILE(L,N)                         
        ENDDO                                                           
      ENDDO                                                             
                                                                        
      DO I=P1,P1+P_POINTS-1                                             
        RADNET_SICE(I) = 0.                                             
        IF (FLANDG(I).LT.1.0 .AND. ICE_FRACT(I).GT.0.)                 
     &    RADNET_SICE(I) = RAD_SICE(I) - ICE_FRACT(I)*SBCON*TI(I)**4    
      ENDDO                                                             
                                                                        
!-----------------------------------------------------------------------
!! Calculate buoyancy parameters BT and BQ.                             
!-----------------------------------------------------------------------
      CALL BOUY_TQ (                                                    
     & P_FIELD,P1,P_POINTS,1                                            
     &,P_1,CF_1,T_1,TL_1,Q_1,QCF_1,QCL_1                                
     &,BT_1,BQ_1,BF_1,BT_CLD_1,BQ_CLD_1                                 
     &,A_QS_1,A_DQSDT_1,DQSDT_1                                         
     &,L_BL_LSPICE,LTIMER                                               
     & )                                                                
                                                                        
!-----------------------------------------------------------------------
!! 4.  Surface turbulent exchange coefficients and "explicit" fluxes    
!!     (P243a, routine SF_EXCH).                                        
!!     Wind mixing "power" and some values required for other, later,   
!!     diagnostic calculations, are also evaluated if requested.        
!-----------------------------------------------------------------------
                                                                        
      CALL SF_EXCH (                                                    
     & P_POINTS,P_FIELD,P1,LAND1,LAND_PTS,LAND_FIELD,NTILES,LAND_INDEX, 
     & TILE_INDEX,TILE_PTS,FLAND,FLANDG,                                
     & BQ_1,BT_1,CANHC_TILE,CANOPY,CATCH,DZSOIL(1),FLAKE,GC,HCONS,      
     & HO2R2_OROG,ICE_FRACT,SNOW_TILE,PSTAR,QW_1,RADNET_SICE,           
     & RADNET_TILE,SIL_OROG_LAND,SMVCST,TILE_FRAC,TIMESTEP,             
     & T_1,Q_1,QCF_1,QCL_1,                                             
     & TL_1,TI,T_SOIL(1,1),
     & TSTAR_TILE,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,
     & VFRAC_TILE,VSHR_LAND,VSHR_SSI,                 
     & ZH,Z0_TILE,Z1,Z1,                                           
     & LAND_MASK,SU10,SV10,SQ1P5,ST1P5,SFME,LTIMER,L_Z0_OROG,Z0MSEA,    
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,CD,CH,CDR10M,CHR1P5M,        
     & CHR1P5M_SICE,E_SEA,FME,FQW_1,FQW_TILE,FQW_ICE,                   
     & FTL_1,FTL_TILE,FTL_ICE,FRACA,H_BLEND_OROG,H_SEA,                 
     & RESFS,RESFT,RIB,RIB_TILE,                                        
     & FB_SURF,U_S,Q1_SD,T1_SD,TV1_SD,Z0M_EFF,                          
     & Z0H,Z0H_TILE,Z0M,Z0M_TILE,RHO_ARESIST,ARESIST,RESIST_B,          
     & RHO_ARESIST_TILE,ARESIST_TILE,RESIST_B_TILE,                     
     & RHO_CD_MODV1,RHOKH_TILE,RHOKH_SICE,RHOKM_1,RHOKM_LAND,RHOKM_SSI,
     & RHOKPM,RHOKPM_SICE,   
     & NRML                                                             
     & )                                                                
                                                                        
*IF DEF,MPP                                                             
! RHOKM(*,1) contains duff data in halos. The P_TO_UV can interpolate   
! this into the real data, so first we must update east/west halos      
                                                                        
! Call fixed to swap north/south halos as well 
      CALL SWAPBOUNDS(RHOKM_1(U1),ROW_LENGTH,N_U_ROWS,1,1,1)
      CALL SWAPBOUNDS(RHOKM_LAND(U1),ROW_LENGTH,N_U_ROWS,1,1,1)
      CALL SWAPBOUNDS(RHOKM_SSI(U1),ROW_LENGTH,N_U_ROWS,1,1,1)          
*ENDIF                                                                  

      DO I=P1,P1-1+P_POINTS                                          
        FLANDG_TMP(I)=FLANDG(I)                                         
      ENDDO                                                             
                                                                        
*IF -DEF,SCMA                                                           
      CALL P_TO_UV (RHOKM_1(P1),RHOKM_UV_1(U1+ROW_LENGTH),              
     &   P_POINTS,U_POINTS,ROW_LENGTH,N_P_ROWS)

      CALL P_TO_UV_LAND (RHOKM_LAND(P1),RHOKM_UV_LAND(U1+ROW_LENGTH),   
     &   P_POINTS,U_POINTS,FLANDG(P1),ROW_LENGTH,N_P_ROWS)
      CALL P_TO_UV_SEA (RHOKM_SSI(P1),RHOKM_UV_SSI(U1+ROW_LENGTH),      
     &   P_POINTS,U_POINTS,FLANDG(P1),ROW_LENGTH,N_P_ROWS)

      CALL P_TO_UV (FLANDG_TMP(P1),FLANDG_UV(U1+ROW_LENGTH),            
     &   P_POINTS,U_POINTS,ROW_LENGTH,N_P_ROWS)                         
*IF DEF,MPP                                                             
      IF (attop) THEN                                                   
*ENDIF                                                                  
      DO I=U1,U1+ROW_LENGTH-1                                           
        RHOKM_UV_1(I) = 1.0E30
        RHOKM_UV_LAND(I) = 1.0E30
        RHOKM_UV_SSI(I) = 1.0E30
        FLANDG_UV(I) = 1.0E30                              
      ENDDO                                                             
*IF DEF,MPP                                                             
      ENDIF                                                             
                                                                        
      IF (atbase) THEN                                                  
*ENDIF                                                                  
      DO I= U1+(N_U_ROWS-1)*ROW_LENGTH, U1+N_U_ROWS*ROW_LENGTH-1        
        RHOKM_UV_1(I) = 1.0E30                                          
        RHOKM_UV_LAND(I) = 1.0E30
        RHOKM_UV_SSI(I) = 1.0E30
        FLANDG_UV(I) = 1.0E30
      ENDDO                                                             
*IF DEF,MPP                                                             
      ENDIF                                                             
*ENDIF                                                                  
                                                                        
*ELSE                                                                   
      DO I = P1, P1-1+P_POINTS                                          
        RHOKM_UV_1(I) = RHOKM_1(I)
        RHOKM_UV_LAND(I) = RHOKM_LAND(I)
        RHOKM_UV_SSI(I) = RHOKM_SSI(I)                                  
        FLANDG_UV(I) = FLANDG(I)
      ENDDO                    
*ENDIF                                                                  
                                                                        
*IF -DEF,SCMA                                                           
        DO I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1                     
*ELSE                                                                   
        DO I=1,U_POINTS                                                 
*ENDIF                                                                  
          TAUX_LAND(I) = RHOKM_UV_LAND(I) * U_1(I)                  
          TAUX_SSI(I) = RHOKM_UV_SSI(I) * ( U_1(I) - U_0(I) )           
          TAUX_1(I) = FLANDG_UV(I)*TAUX_LAND(I)                     
     &            +(1.-FLANDG_UV(I))*TAUX_SSI(I)

          TAUY_LAND(I) = RHOKM_UV_LAND(I) * V_1(I)                  
          TAUY_SSI(I) = RHOKM_UV_SSI(I) * ( V_1(I) - V_0(I) )           
          TAUY_1(I) = FLANDG_UV(I)*TAUY_LAND(I)                     
     &            +(1.-FLANDG_UV(I))*TAUY_SSI(I)                        
        ENDDO                                                           
                                                                        
*IF -DEF,SCMA                                                           
!-----------------------------------------------------------------------
!! Set first and last rows to "missing data indicator"                  
!-----------------------------------------------------------------------
*IF DEF,MPP                                                             
      IF (attop) THEN                                                   
*ENDIF                                                                  
      DO I=U1,U1+ROW_LENGTH-1 
        TAUX_LAND(I)=1.E30
        TAUX_SSI(I)=1.E30                                        
        TAUX_1(I)=1.E30
        
        TAUY_LAND(I)=1.E30
        TAUY_SSI(I) =1.E30                                              
        TAUY_1(I)=1.E30                                                 
      ENDDO                                                             
*IF DEF,MPP                                                             
      ENDIF                                                             
                                                                        
      IF (atbase) THEN                                                  
*ENDIF                                                                  
      DO I= U1 + (N_U_ROWS-1)*ROW_LENGTH, U1 + N_U_ROWS*ROW_LENGTH -1   
        TAUX_LAND(I)=1.E30
        TAUX_SSI(I)=1.E30                                        
        TAUX_1(I)=1.E30
        
        TAUY_LAND(I)=1.E30
        TAUY_SSI(I) =1.E30                                              
        TAUY_1(I)=1.E30                                                 
      ENDDO                                                             
*IF DEF,MPP                                                             
      ENDIF                                                             
*ENDIF                                                                  
*ENDIF                                                                  
                                                                        
        IF (SU10. OR. SV10)THEN                                         

*IF DEF,MPP
! Fix diags by swapping halos
        CALL SWAPBOUNDS(CDR10M(U1),ROW_LENGTH,N_U_ROWS,1,1,1)
*ENDIF

*IF -DEF,SCMA                                                           
        CALL P_TO_UV (CDR10M(P1),CDR10M_UV(U1+ROW_LENGTH),P_POINTS,     
     &     U_POINTS,ROW_LENGTH,N_P_ROWS)                                
!-----------------------------------------------------------------------
!! Set first and last rows to "missing data indicator"                  
!-----------------------------------------------------------------------
*IF DEF,MPP                                                             
        IF (attop) THEN                                                 
*ENDIF                                                                  
          DO I=U1,U1+ROW_LENGTH-1                                       
            CDR10M_UV(I) = 1.0E30                                       
          ENDDO                                                         
*IF DEF,MPP                                                             
        ENDIF                                                           
                                                                        
        IF (atbase) THEN                                                
*ENDIF                                                                  
          DO I= U1+(N_U_ROWS-1)*ROW_LENGTH, U1+N_U_ROWS*ROW_LENGTH-1    
            CDR10M_UV(I) = 1.0E30                                       
          ENDDO                                                         
*IF DEF,MPP                                                             
        ENDIF                                                           
*ENDIF                                                                  
                                                                        
*ELSE                                                                   
      DO I = P1, P1-1+P_POINTS                                          
        CDR10M_UV(I) = CDR10M(I)                                        
      ENDDO                                                             
*ENDIF                                                                  
        ENDIF                                                           
                                                                        
                                                                        
!-----------------------------------------------------------------------
!! Set grid-box surface exchange coefficients                           
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1                                             
        IF( FLANDG(I).LT.1.0 ) THEN                                   
          RHOKH(I) = (1.0 - FLANDG(I))*RHOKH_SICE(I)                    
        ELSE                                                            
          RHOKH(I) = 0.0                                                
        ENDIF                                                           
      ENDDO                                                             
                                                                        
      DO N=1,NTILES                                                     
        DO J=1,TILE_PTS(N)                                              
          L = TILE_INDEX(J,N)                                           
          I = LAND_INDEX(L)                                             
          RHOKH(I) = RHOKH(I) 
     &      + FLANDG(I)*TILE_FRAC(L,N)*RHOKH_TILE(L,N)          
        ENDDO                                                           
      ENDDO                                                             
                                                                        
                                                                        
      IF (LTIMER) THEN                                                  
        CALL TIMER('BDYLAYR ',4)                                        
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
*ENDIF                                                                  
*DECK SFIMPL8A                                                          
*IF DEF,A03_7A,OR,DEF,A03_8A                                            
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
!!!  SUBROUTINE SF_IMPL-----------------------------------------------  
!!!                                                                     
!!!  Purpose: Calculate implicit correction to surface fluxes of heat,  
!!!           moisture and momentum. Also calculates screen level       
!!!           temperature and humidity as well as 10 m winds.           
!!!                                                                     
!!!                                                                     
!!! F.Hewer     <- programmer of some or all of previous code or changes
!!! C.Wilson    <- programmer of some or all of previous code or changes
!!!                                                                     
!!!  Model            Modification history:                             
!!! version  Date                                                       
!!!                                                                     
!!!   4.3  7/2/97     New deck. S Jackson                               
!!!   4.4 25/6/97     Modified for MOSES II tile model. R Essery        
!!!   4.4 25/6/97     Move grid definitions up to BL_INTCT.  R.A.Betts  
!!!  4.5    Jul. 98  Kill the IBM specific lines. (JCThil)              
!!!   4.5  7/5/98     Set TSTAR, SNOW_SURF_HTF and SOIL_SURF_HTF to 0   
!!!                   at all land points, to avoid problems of          
!!!                   non-initialised data.  R.A.Betts                  
!!!   4.5 21/5/98     Add optional error check for negative surface     
!!!                   temperature.  R.A.Betts                           
!!!                                                                     
!!!  Programming standard: Unified Model Documentation Paper No 4,      
!!!                        Version ?, dated ?.                          
!!!                                                                     
!!!  System component covered: P24.                                     
!!!                                                                     
!!!  Project task:                                                      
!!!                                                                     
!!!  Documentation: UMDP 24.                                            
!!!                                                                     
!!!---------------------------------------------------------------------
                                                                        
!    Arguments :-                                                       
      SUBROUTINE SF_IMPL (                                              
                                                                        
! IN values defining field dimensions and subset to be processed :      
     & P_FIELD,U_FIELD,LAND_FIELD,ROW_LENGTH,                           
     & P_POINTS,P1,LAND1,LAND_PTS,U_POINTS,U1,                          
                                                                        
! IN soil/vegetation/land surface data :                                
     & LAND_INDEX,LAND_MASK,                                            
     & NTILES,TILE_INDEX,TILE_PTS,SM_LEVELS,                            
     & CANHC_TILE,CANOPY,FLAKE,SMC,                                     
     & TILE_FRAC,WT_EXT_TILE,                                           
     & FLAND,FLANDG,                                                    
                                                                        
! IN sea/sea-ice data :                                                 
     & DI,ICE_FRACT,U_0,V_0,                                            
                                                                        
! IN everything not covered so far :                                    
     & PSTAR,LW_DOWN,RAD_SICE,SW_TILE,TIMESTEP,                         
     & T_SOIL,QW_1,TL_1,U_1,V_1,RHOKM_UV_1,                             
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,                             
     & DTRDZ_1,DU_1,DV_1,FQW_TILE,FQW_ICE,FTL_ICE,TSTAR_TILE_OLD,       
     & FRACA,RESFS,RESFT,RHOKH_TILE,RHOKH_SICE,RHOKPM,RHOKPM_SICE,      
     & Z1,Z0H,Z0H_TILE,Z0M,Z0M_TILE,CDR10M_UV,CHR1P5M,CHR1P5M_SICE,     
     & CT_CTQ_1,DQW_1,DTL_1,CQ_CM_1,                                    
     & L_NEG_TSTAR,                                                     
     & FLANDG_UV,                                               
                                                                        
! IN STASH flags :-                                                     
     & SIMLT,SMLT,SLH,SQ1P5,ST1P5,SU10,SV10,                            
                                                                        
! INOUT data :                                                          
     & TI,TSTAR,                                                        
     & TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,TSTAR_SSI,                       
     & TSTAR_TILE,SNOW_TILE,                                            
     & LE_TILE,RADNET_SICE,RADNET_TILE,                                 
     & E_SEA,FQW_1,FTL_1,FTL_TILE,H_SEA,OLR,TAUX_1,TAUY_1,              
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           
                                                                        
! OUT Diagnostic not requiring STASH flags :                            
     & ECAN,EI_TILE,ESOIL_TILE,                                         
     & SEA_ICE_HTF,SURF_HT_FLUX,SURF_HT_FLUX_LAND,SURF_HT_FLUX_SICE,    
                                                                        
! OUT diagnostic requiring STASH flags :                                
     & SICE_MLT_HTF,SNOMLT_SURF_HTF,LATENT_HEAT,                        
     & Q1P5M,Q1P5M_TILE,T1P5M,T1P5M_TILE,U10M,V10M,                     
                                                                        
! OUT data required elsewhere in UM system :                            
     & ECAN_TILE,EI,ESOIL,EXT,SNOWMELT,MELT_TILE,                       
     & ERROR,                                                           
                                                                        
! LOGICAL LTIMER                                                        
     & LTIMER                                                           
     & )                                                                
                                                                        
      IMPLICIT NONE                                                     
                                                                        
!  Inputs :-                                                            
                                                                        
! (a) Defining horizontal grid and subset thereof to be processed.      
!    Checked for consistency.                                           
                                                                        
      INTEGER                                                           
     & P_FIELD                     ! IN No. of P-points in whole grid   
!                                  !    (for dimensioning only).        
     &,U_FIELD                     ! IN No. of UV-points in whole grid. 
     &,LAND_FIELD                  ! IN No.of land points in whole grid.
     &,ROW_LENGTH                  ! IN No. of points in one row.       
     &,P_POINTS   ! IN No of P-points being processed.                  
     &,P1         ! IN First P-point to be processed.                   
     &,LAND1      ! IN First land-point to be processed.                
!                 !       1 <= LAND1 <= LAND_FIELD                      
     &,LAND_PTS   ! IN No of land points being processed.               
     &,U_POINTS   ! IN No of UV-points being processed.                 
     &,U1         ! IN First UV-point to be processed.                  
                                                                        
! (c) Soil/vegetation/land surface parameters (mostly constant).        
                                                                        
      LOGICAL                                                           
     & LAND_MASK(P_FIELD)          ! IN T if land, F elsewhere.         
                                                                        
      INTEGER                                                           
     & LAND_INDEX(P_FIELD)         ! IN LAND_INDEX(I)=J => the Jth      
!                                  !    point in P_FIELD is the Ith     
!                                  !    land point.                     
                                                                        
      INTEGER                                                           
     & SM_LEVELS                   ! IN No. of soil moisture levels     
     &,NTILES                      ! IN No. of land tiles               
     &,TILE_INDEX(LAND_FIELD,NTILES)!IN Index of tile points            
     &,TILE_PTS(NTILES)            ! IN Number of tile points           
                                                                        
      REAL                                                              
     & CANHC_TILE(LAND_FIELD,NTILES)!IN Areal heat capacity of canopy   
!                                  !    for land tiles (J/K/m2).        
     &,CANOPY(LAND_FIELD,NTILES)   ! IN Surface/canopy water for        
!                                  !    snow-free land tiles (kg/m2)    
     &,FLAKE(LAND_FIELD,NTILES)    ! IN Lake fraction.                  
     &,SMC(LAND_FIELD)             ! IN Available soil moisture (kg/m2).
     &,TILE_FRAC(LAND_FIELD,NTILES)! IN Tile fractions including        
!                                  ! snow cover in the ice tile.        
     &,WT_EXT_TILE(LAND_FIELD,SM_LEVELS,NTILES)                         
!                                  ! IN Fraction of evapotranspiration  
!                                  !    extracted from each soil layer  
!                                  !    by each tile.                   
     &,FLAND(LAND_FIELD)           ! IN Land fraction on land pts.      
     &,FLANDG(P_FIELD)             ! IN Land fraction on all pts.       
     &,FLANDG_UV(U_FIELD)          ! IN Land fraction on UV grid.       
                                                                        
! (d) Sea/sea-ice data.                                                 
                                                                        
      REAL                                                              
     & DI(P_FIELD)                 ! IN "Equivalent thickness" of       
!                                  !     sea-ice(m).                    
     &,ICE_FRACT(P_FIELD)          ! IN Fraction of gridbox covered by  
!                                  !     sea-ice (decimal fraction).    
     &,U_0(U_FIELD)                ! IN W'ly component of surface       
!                                  !    current (m/s).                  
     &,V_0(U_FIELD)                ! IN S'ly component of surface       
!                                  !    current (m/s).                  
                                                                        
! (f) Atmospheric + any other data not covered so far, incl control.    
                                                                        
      REAL                                                              
     & PSTAR(P_FIELD)              ! IN Surface pressure (Pascals).     
     &,LW_DOWN(P_FIELD)            ! IN Surface downward LW radiation   
!                                  !    (W/m2).                         
     &,RAD_SICE(P_FIELD)           ! IN Surface net SW and downward LW  
!                                  !    radiation for sea-ice (W/sq m). 
     &,SW_TILE(LAND_FIELD,NTILES)  ! IN Surface net SW radiation on     
!                                  !    land tiles (W/m2).              
     &,TIMESTEP                    ! IN Timestep (seconds).             
     &,T_SOIL(LAND_FIELD,SM_LEVELS)! IN Soil temperatures (K).          
     &,QW_1(P_FIELD)               ! IN Total water content             
     &,TL_1(P_FIELD)               ! IN Ice/liquid water temperature    
     &,U_1(U_FIELD)                ! IN W'ly wind component (m/s)       
     &,V_1(U_FIELD)                ! IN S'ly wind component (m/s)       
     &,RHOKM_UV_1(U_FIELD)         ! IN Exchange coefficients for       
!                                  !    momentum (on UV-grid, with 1st  
!                                  !    and last rows undefined or, at  
!                                  !    present, set to "missing data") 
                                                                        
      REAL                                                              
     & ALPHA1(LAND_FIELD,NTILES)   ! IN Mean gradient of saturated      
!                                  !    specific humidity with respect  
!                                  !    to temperature between the      
!                                  !    bottom model layer and tile     
!                                  !    surfaces                        
     &,ALPHA1_SICE(P_FIELD)        ! IN ALPHA1 for sea-ice.             
     &,ASHTF(P_FIELD)              ! IN Coefficient to calculate surface
!                                  !    heat flux into soil or sea-ice. 
     &,ASHTF_TILE(LAND_FIELD,NTILES)!IN Coefficient to calculate        
!                                  !    surface heat flux into land     
!                                  !    tiles.                          
     &,DTRDZ_1(P_FIELD)            ! IN -g.dt/dp for model layers.      
     &,FQW_TILE(LAND_FIELD,NTILES) ! IN Surface FQW for land tiles      
     &,FQW_ICE(P_FIELD)            ! IN Surface FQW for sea-ice         
     &,FTL_ICE(P_FIELD)            ! IN Surface FTL for sea-ice         
     &,TSTAR_TILE_OLD(LAND_FIELD,NTILES)                                
!                                  ! IN Tile surface temperatures at    
!                                  !    beginning of timestep.          
     &,FRACA(LAND_FIELD,NTILES)    ! IN Fraction of surface moisture    
!                                  !    flux with only aerodynamic      
!                                  !    resistance for snow-free land   
!                                  !    tiles.                          
     &,RESFS(LAND_FIELD,NTILES)    ! IN Combined soil, stomatal         
!                                  !    and aerodynamic resistance      
!                                  !    factor for fraction (1-FRACA) of
!                                  !    snow-free land tiles.           
     &,RESFT(LAND_FIELD,NTILES)    ! IN Total resistance factor.        
!                                  !    FRACA+(1-FRACA)*RESFS for       
!                                  !    snow-free land, 1 for snow.     
     &,RHOKH_TILE(LAND_FIELD,NTILES)!IN Surface exchange coefficients   
!                                  !    for land tiles                  
     &,RHOKH_SICE(P_FIELD)         ! IN Surface exchange coefficients   
!                                  !    for sea and sea-ice             
     &,RHOKPM(LAND_FIELD,NTILES)   ! IN Land surface exchange coeff.    
     &,RHOKPM_SICE(P_FIELD)        ! IN Sea-ice surface exchange coeff. 
                                                                        
       REAL                                                             
     & Z1(P_FIELD)                 ! IN Height of lowest level (i.e.    
!                                  !    height of middle of lowest      
!                                  !    layer).                         
     &,Z0H(P_FIELD)                ! IN Roughness length for heat and   
!                                  !    moisture (m).                   
     &,Z0H_TILE(LAND_FIELD,NTILES) ! IN Tile roughness lengths for heat 
!                                  !    and moisture (m).               
     &,Z0M(P_FIELD)                ! IN Roughness length for momentum   
!                                  !    (m).                            
     &,Z0M_TILE(LAND_FIELD,NTILES) ! IN Tile roughness lengths for      
!                                  !    momentum.                       
     &,CDR10M_UV(U_FIELD)          ! IN Ratio of CD's reqd for          
!                                  !    calculation of 10 m wind. On    
!                                  !    UV-grid; comments as per RHOKM. 
     &,CHR1P5M(LAND_FIELD,NTILES)  ! IN Ratio of coefffs for calculation
!                                  !    of 1.5m temp for land tiles.    
     &,CHR1P5M_SICE(P_FIELD)       ! IN CHR1P5M for sea and sea-ice     
!                                  !    (leads ignored).                
     &,CT_CTQ_1(P_FIELD)           ! IN Coefficient in T and q          
!                                       tri-diagonal implicit matrix    
     &,CQ_CM_1(U_FIELD)            ! IN Coefficient in U and V          
!                                       tri-diagonal implicit matrix    
     &,DQW_1(P_FIELD)              ! IN Level 1 increment to q field    
     &,DTL_1(P_FIELD)              ! IN Level 1 increment to T field    
     &,DU_1(U_FIELD)               ! IN Level 1 increment to u wind     
!                                       field                           
     &,DV_1(U_FIELD)               ! IN Level 1 increment to v wind     
!                                       field                           
                                                                        
      LOGICAL                                                           
     & LTIMER                      ! IN Logical switch for TIMER diags  
     &,L_NEG_TSTAR                ! IN Switch for -ve TSTAR error check 
                                                                        
!  STASH flags :-                                                       
                                                                        
      LOGICAL                                                           
     & SIMLT   ! IN Flag for SICE_MLT_HTF (q.v.)                        
     &,SMLT    ! IN Flag for SNOMLT_SURF_HTF (q.v.)                     
     &,SLH     ! IN Flag for LATENT_HEAT (q.v.)                         
     &,SQ1P5   ! IN Flag for Q1P5M (q.v.)                               
     &,ST1P5   ! IN Flag for T1P5M (q.v.)                               
     &,SU10    ! IN Flag for U10M (q.v.)                                
     &,SV10    ! IN Flag for V10M (q.v.)                                
                                                                        
!  In/outs :-                                                           
                                                                        
      REAL                                                              
     & TI(P_FIELD)                 ! INOUT Sea-ice surface layer        
!                                  !       temperature (K).             
     &,TSTAR(P_FIELD)              ! OUT   GBM surface temperature (K). 
     &,TSTAR_LAND(P_FIELD)         ! OUT   Land mean sfc temperature (K)
     &,TSTAR_SEA(P_FIELD)          ! IN    Open sea sfc temperature (K).
     &,TSTAR_SICE(P_FIELD)         ! OUT   Sea-ice sfc temperature (K). 
     &,TSTAR_SSI(P_FIELD)          ! INOUT Sea mean sfc temperature (K).
     &,TSTAR_TILE(LAND_FIELD,NTILES)!INOUT Surface tile temperatures    
     &,SNOW_TILE(LAND_FIELD,NTILES)! INOUT Lying snow on tiles (kg/m2)  
     &,LE_TILE(LAND_FIELD,NTILES)  ! INOUT Surface latent heat flux for 
!                                  !     land tiles                     
     &,RADNET_SICE(P_FIELD)        ! INOUT Surface net radiation on     
!                                  !       sea-ice (W/m2)               
     &,RADNET_TILE(LAND_FIELD,NTILES)                                   
!                                  ! INOUT Surface net radiation on     
!                                  !     land tiles (W/m2)              
     &,E_SEA(P_FIELD)              ! INOUT Evaporation from sea times   
!                                  !       leads fraction. Zero over    
!                                  !       land. (kg per square metre   
!                                  !       per sec).                    
     &,FQW_1(P_FIELD)              ! INOUT Moisture flux between layers 
!                                  !       (kg per square metre per sec)
!                                  !       FQW(,1) is total water flux  
!                                  !       from surface, 'E'.           
     &,FTL_1(P_FIELD)              ! INOUT FTL(,K) contains net         
!                                  !       turbulent sensible heat flux 
!                                  !       into layer K from below; so  
!                                  !       FTL(,1) is the surface       
!                                  !       sensible heat, H.(W/m2)      
     &,FTL_TILE(LAND_FIELD,NTILES) ! INOUT Surface FTL for land tiles   
     &,H_SEA(P_FIELD)              ! INOUT Surface sensible heat flux   
!                                  !       over sea times leads fraction
!                                  !       (W/m2)                       
     &,OLR(P_FIELD)                ! IN    TOA - surface upward LW on   
!                                  !       last radiation timestep      
!                                  ! OUT   Corrected TOA outward LW     
     &,TAUX_1(U_FIELD)             ! OUT   W'ly component of surface    
!                                  !       wind stress (N/sq m). (On    
!                                  !       UV-grid with first and last  
!                                  !       rows undefined or, at        
!                                  !       present, set to missing data 
     &,TAUX_LAND(U_FIELD)          ! INOUT W'ly component of surface    
!                                  !       wind stress over land        
!                                  !       (N/sq m). (On                
!                                  !       UV-grid with first and last  
!                                  !       rows undefined or, at        
!                                  !       present, set to missing data 
     &,TAUX_SSI(U_FIELD)           ! INOUT W'ly component of surface    
!                                  !       wind stress over mean sea    
!                                  !       (N/sq m). (On                
!                                  !       UV-grid with first and last  
!                                  !       rows undefined or, at        
!                                  !       present, set to missing data 
     &,TAUY_1(U_FIELD)             ! OUT   S'ly component of surface    
!                                  !       wind stress (N/sq m).  On    
!                                  !       UV-grid; comments as per TAUX
     &,TAUY_LAND(U_FIELD)          ! INOUT S'ly component of land sfc   
!                                  !       wind stress (N/sq m).  On    
!                                  !       UV-grid; comments as per TAUX
     &,TAUY_SSI(U_FIELD)           ! INOUT S'ly compt of mean sea sfc   
!                                  !       wind stress (N/sq m).  On    
!                                  !       UV-grid; comments as per TAUX
                                                                        
!  Outputs :-                                                           
!-1 Diagnostic (or effectively so - includes coupled model requisites):-
                                                                        
!  (a) Calculated anyway (use STASH space from higher level) :-         
!                                                                       
      REAL                                                              
     & ECAN(P_FIELD)               ! OUT Gridbox mean evaporation from  
!                                  !     canopy/surface store (kg/m2/s).
!                                  !     Zero over sea.                 
     &,ESOIL_TILE(LAND_FIELD,NTILES)                                    
!                                  ! OUT ESOIL for snow-free land tiles 
     &,SEA_ICE_HTF(P_FIELD)        ! OUT Heat flux through sea-ice      
!                                  !     (W/m2, positive downwards).    
     &,SURF_HT_FLUX(P_FIELD)       ! OUT Net downward heat flux at      
!                                  !     surface over land and sea-ice  
!                                  !     fraction of gridbox (W/m2).    
     &,SURF_HT_FLUX_LAND(P_FIELD)  ! OUT Net downward heat flux at      
!                                  !     surface over land              
!                                  !     fraction of gridbox (W/m2).    
     &,SURF_HT_FLUX_SICE(P_FIELD)  ! OUT Net downward heat flux at      
!                                  !     surface over sea-ice           
!                                  !     fraction of gridbox (W/m2).    
                                                                        
!  (b) Not passed between lower-level routines (not in workspace at this
!      level) :-                                                        
                                                                        
      REAL                                                              
     & SICE_MLT_HTF(P_FIELD)       ! OUT Heat flux due to melting of    
!                                  !     sea-ice (Watts per sq metre).  
     &,SNOMLT_SURF_HTF(P_FIELD)    ! OUT Heat flux required for surface 
!                                  !     melting of snow (W/m2).        
     &,LATENT_HEAT(P_FIELD)        ! OUT Surface latent heat flux, +ve  
!                                  !     upwards (Watts per sq m).      
     &,Q1P5M(P_FIELD)              ! OUT Q at 1.5 m (kg water / kg air).
     &,Q1P5M_TILE(LAND_FIELD,NTILES)!OUT Q1P5M over land tiles.         
     &,T1P5M(P_FIELD)              ! OUT T at 1.5 m (K).                
     &,U10M(U_FIELD)               ! OUT U at 10 m (m per s).           
     &,T1P5M_TILE(LAND_FIELD,NTILES)!OUT T1P5M over land tiles.         
     &,V10M(U_FIELD)               ! OUT V at 10 m (m per s).           
                                                                        
!-2 Genuinely output, needed by other atmospheric routines :-           
                                                                        
      REAL                                                              
     & EI(P_FIELD)                 ! OUT Sublimation from lying snow or 
!                                  !     sea-ice (kg/m2/s).             
     &,EI_LAND(P_FIELD)            ! OUT Sublimation from lying snow    
!                                  !     (kg/m2/s).                     
     &,EI_SICE(P_FIELD)            ! OUT Sublimation from sea-ice       
!                                  !     (kg/m2/s).                     
     &,EI_TILE(LAND_FIELD,NTILES)  ! OUT EI for land tiles.             
     &,ECAN_TILE(LAND_FIELD,NTILES)! OUT ECAN for snow-free land tiles  
     &,ESOIL(P_FIELD)              ! OUT Surface evapotranspiration     
!                                  !     from soil moisture store       
!                                  !     (kg/m2/s).                     
     &,EXT(LAND_FIELD,SM_LEVELS)   ! OUT Extraction of water from each  
!                                  !     soil layer (kg/m2/s).          
     &,SNOWMELT(P_FIELD)           ! OUT Snowmelt (kg/m2/s).            
     &,MELT_TILE(LAND_FIELD,NTILES)! OUT Snowmelt on land tiles (kg/m2/s
                                                                        
      INTEGER                                                           
     & ERROR          ! OUT 0 - AOK;                                    
!                     !     1 to 7  - bad grid definition detected;     
                                                                        
!---------------------------------------------------------------------  
!  External routines called :-                                          
                                                                        
      EXTERNAL Z,HEAT_CON,SF_EXCH,BOUY_TQ,BTQ_INT,                      
     & KMKH,EX_FLUX_TQ,EX_FLUX_UV,IM_CAL_TQ,SICE_HTF,SF_EVAP,SF_MELT,   
     & IM_CAL_UV,SCREEN_TQ                                              
      EXTERNAL TIMER                                                    
*IF -DEF,SCMA                                                           
      EXTERNAL UV_TO_P,P_TO_UV                                          
*ENDIF                                                                  
                                                                        
!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-        
                                                                        
*CALL C_0_DG_C                                                          
*CALL C_R_CP                                                            
*CALL C_G                                                               
*CALL C_LHEAT                                                           
*CALL C_GAMMA                                                           
*CALL CSIGMA                                                            
*CALL SOIL_THICK                                                        
*IF DEF,MPP                                                             
! MPP Common block                                                      
*CALL PARVARS                                                           
*ENDIF                                                                  
                                                                        
! Derived local parameters.                                             
                                                                        
      REAL LS                                                           
                                                                        
      PARAMETER (                                                       
     & LS=LF+LC              ! Latent heat of sublimation.              
     &  )                                                               
                                                                        
!-----------------------------------------------------------------------
                                                                        
!  Workspace :-                                                         
                                                                        
      REAL                                                              
     & ELAKE_TILE(LAND_FIELD,NTILES) ! Lake evaporation.                
     &,QIM_1(P_FIELD)  ! Implicit value of first model level humidity   
     &,TIM_1(P_FIELD)  ! Implicit value of first model level temperature
     &,TSTAR_RAD4(P_FIELD)! Effective surface radiative temperature for 
!                         ! land and sea-ice                            
                                                                        
!  Local scalars :-                                                     
                                                                        
      REAL                                                              
     & LAT_HT     ! Latent heat of evaporation for snow-free land       
!                 ! or sublimation for snow-covered land and ice.       
                                                                        
      INTEGER                                                           
     & I,J,L      ! LOCAL Loop counter (horizontal field index).        
     &,N          ! LOCAL Loop counter (tile index).                    
                                                                        
      IF (LTIMER) THEN                                                  
        CALL TIMER('BDYLAYR ',3)                                        
      ENDIF                                                             
      ERROR = 0                                                         
                                                                        
      CALL IM_SF_PT (                                                   
     & P_FIELD,P1,U_FIELD,U1                                            
     &,P_POINTS,U_POINTS,ROW_LENGTH,LAND_FIELD                          
     &,LAND_INDEX,NTILES,TILE_INDEX,TILE_PTS                            
     &,FLANDG,TILE_FRAC,SNOW_TILE,ICE_FRACT                             
     &,GAMMA(1),ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE                     
     &,RESFT,RHOKPM,RHOKPM_SICE                                         
     &,RHOKM_UV_1,RHOKH_TILE,RHOKH_SICE                                 
     &,CT_CTQ_1,DQW_1,DTL_1,CQ_CM_1,DU_1,DV_1
     &,FLANDG_UV                                                
     &,FQW_1,FTL_1                                                      
     &,TAUX_1,TAUX_LAND,TAUX_SSI,TAUY_1,TAUY_LAND,TAUY_SSI              
     &,FQW_TILE,FTL_TILE,FQW_ICE,FTL_ICE,E_SEA,H_SEA                    
     &,LTIMER                                                           
     &)                                                                 
                                                                        
!-----------------------------------------------------------------------
!! 6.1 Convert FTL to sensible heat flux in Watts per square metre.     
!-----------------------------------------------------------------------
                                                                        
Cfpp$ Select(CONCUR)                                                    
      DO  I=P1,P1+P_POINTS-1                                            
        FTL_1(I) = FTL_1(I)*CP                                          
      ENDDO                                                             
                                                                        
      DO I=P1,P1+P_POINTS-1                                             
        FTL_ICE(I) = CP*FTL_ICE(I)                                      
      ENDDO                                                             
                                                                        
      DO N=1,NTILES                                                     
        DO J=1,TILE_PTS(N)                                              
          L = TILE_INDEX(J,N)                                           
          FTL_TILE(L,N) = CP*FTL_TILE(L,N)                              
        ENDDO                                                           
      ENDDO                                                             
                                                                        
!-----------------------------------------------------------------------
! Diagnose the GBM surface temperature for points with sea-ice          
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1                                             
        IF ( FLANDG(I).LT.1.0 .and. ICE_FRACT(I).GT.0. ) THEN       
          SURF_HT_FLUX_SICE(I) = RADNET_SICE(I) - LS*FQW_ICE(I) - 
     &                        FTL_ICE(I)                              
                                                                        
          TSTAR_SSI(I) = (1. - ICE_FRACT(I))*TSTAR_SEA(I) +       
     &                   ICE_FRACT(I)*TI(I) +                       
     &                    SURF_HT_FLUX_SICE(I) / ASHTF(I)           
        ENDIF                                                           
      ENDDO                                                             
                                                                        
!-----------------------------------------------------------------------
! Optional error check : test for negative top soil layer temperature   
!-----------------------------------------------------------------------
      IF (L_NEG_TSTAR) THEN                                             
        DO L=LAND1,LAND1+LAND_PTS-1                                     
          IF (T_SOIL(L,1).LT.0) THEN                                    
            ERROR = 1                                                   
            WRITE(6,*) '*** ERROR DETECTED BY ROUTINE BDY_LAYR ***'     
            WRITE(6,*) 'NEGATIVE TEMPERATURE IN TOP SOIL LAYER AT '     
            WRITE(6,*) 'LAND POINT ',L                                  
          ENDIF                                                         
        ENDDO                                                           
      ENDIF                                                             
                                                                        
!-----------------------------------------------------------------------
!!   Diagnose the land surface temperature                              
!-----------------------------------------------------------------------
                                                                        
      DO N=1,NTILES                                                     
        DO L=LAND1,LAND1+LAND_PTS-1                                     
          TSTAR_TILE(L,N) = T_SOIL(L,1)                                 
          IF (SNOW_TILE(L,N).GT.0.)                                     
     &      TSTAR_TILE(L,N) =  MIN( T_SOIL(L,1), TM )                   
        ENDDO                                                           
      ENDDO                                                             
                                                                        
      DO N=1,NTILES                                                     
        DO J=1,TILE_PTS(N)                                              
          L = TILE_INDEX(J,N)                                           
          LAT_HT = LC                                                   
          IF (SNOW_TILE(L,N).GT.0.) LAT_HT = LS                         
          TSTAR_TILE(L,N) = T_SOIL(L,1) + ( RADNET_TILE(L,N)            
     &                          - LAT_HT*FQW_TILE(L,N) - FTL_TILE(L,N)  
     &                          + (CANHC_TILE(L,N)/TIMESTEP) *          
     &                            (TSTAR_TILE_OLD(L,N) - T_SOIL(L,1)) ) 
     &                                                 / ASHTF_TILE(L,N)
        ENDDO                                                           
      ENDDO                                                             
                                                                        
!-----------------------------------------------------------------------
!! 7.  Surface evaporation components and updating of surface           
!!     temperature (P245, routine SF_EVAP).                             
!-----------------------------------------------------------------------
      CALL SF_EVAP (                                                    
     & P_POINTS,P_FIELD,P1,LAND1,LAND_PTS,LAND_FIELD,NTILES,            
     & LAND_INDEX,TILE_INDEX,TILE_PTS,SM_LEVELS,LTIMER,FLAND,           
     & ASHTF_TILE,CANOPY,DTRDZ_1,FLAKE,FRACA,SNOW_TILE,RESFS,           
     & RESFT,RHOKH_TILE,TILE_FRAC,SMC,WT_EXT_TILE,TIMESTEP,             
     & FQW_1,FQW_TILE,FTL_1,FTL_TILE,TSTAR_TILE,                        
     & ECAN,ECAN_TILE,ELAKE_TILE,ESOIL,ESOIL_TILE,EI_TILE,EXT           
     & )                                                                
                                                                        
!-----------------------------------------------------------------------
!!     Surface melting of sea-ice and snow on land tiles.               
!-----------------------------------------------------------------------
                                                                        
      CALL SF_MELT (                                                    
     & P_POINTS,P_FIELD,P1,LAND_FIELD,NTILES,LAND_INDEX,                
     & TILE_INDEX,TILE_PTS,LAND_MASK,LTIMER,SIMLT,SMLT,FLANDG,          
     & ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE,DTRDZ_1,ICE_FRACT,           
     & RHOKH_TILE,RHOKH_SICE,TILE_FRAC,TIMESTEP,                        
     & EI_TILE,FQW_1,FQW_ICE,FTL_1,FTL_ICE,FTL_TILE,                    
     & TSTAR_SEA,TSTAR_SSI,TSTAR_TILE,SNOW_TILE,                        
     & EI_LAND,EI_SICE,                                                 
     & SICE_MLT_HTF,SNOMLT_SURF_HTF,SNOWMELT,MELT_TILE                  
     & )                                                                
                                                                        
*IF -DEF,SCMA                                                           
      DO I=P1,P1+P_POINTS-1                                             
*ELSE                                                                   
      DO I=1,P_POINTS                                                   
*ENDIF                                                                  
        QIM_1(I)=QW_1(I) + DQW_1(I)-CT_CTQ_1(I)*FQW_1(I)                
        TIM_1(I)=TL_1(I) + DTL_1(I)-CT_CTQ_1(I)*FTL_1(I)/CP             
      ENDDO                                                             
                                                                        
!-----------------------------------------------------------------------
!!     Specific humidity and temperature at 1.5 metres.                 
!-----------------------------------------------------------------------
      CALL SCREEN_TQ (                                                  
     & P_POINTS,P_FIELD,P1,LAND1,LAND_PTS,LAND_FIELD,NTILES,            
     & LAND_INDEX,TILE_INDEX,TILE_PTS,FLANDG,                        
     & SQ1P5,ST1P5,CHR1P5M,CHR1P5M_SICE,PSTAR,QIM_1,RESFT,              
     & TILE_FRAC,TIM_1,TSTAR_SSI,TSTAR_TILE,                            
     & Z0H,Z0H_TILE,Z0M,Z0M_TILE,Z1,                                    
     & Q1P5M,Q1P5M_TILE,T1P5M,T1P5M_TILE                                
     & )                                                                
                                                                        
!-----------------------------------------------------------------------
!!     Gridbox-mean surface temperature and net surface heat fluxes     
!-----------------------------------------------------------------------
      DO I=P1,P1+P_POINTS-1                                             
        SURF_HT_FLUX_LAND(I) = 0.                                     
        SURF_HT_FLUX_SICE(I) = 0.                                     
        TSTAR_RAD4(I) = 0.                                            
        IF (FLANDG(I).LT.1.0 .AND. ICE_FRACT(I).GT.0.) THEN         
          TSTAR_SICE(I) = (TSTAR_SSI(I) -                           
     &              (1.-ICE_FRACT(I))*TSTAR_SEA(I))/ICE_FRACT(I)  
          TSTAR_RAD4(I) = (1.0-FLANDG(I))                           
     &              *ICE_FRACT(I)*TSTAR_SICE(I)**4                  
          RADNET_SICE(I) = RAD_SICE(I) -                            
     &                     ICE_FRACT(I)*SBCON*TSTAR_SICE(I)**4      
          SURF_HT_FLUX_SICE(I) = RADNET_SICE(I) - LS*FQW_ICE(I) - 
     &                        FTL_ICE(I)                              
        ENDIF                                                           
      ENDDO                                                             
                                                                        
      DO L=1,LAND_FIELD                                                 
        I = LAND_INDEX(L)                                               
        TSTAR_LAND(I) = 0.                                              
      ENDDO                                                             
                                                                        
      DO N=1,NTILES                                                     
        DO J=1,TILE_PTS(N)                                              
          L = TILE_INDEX(J,N)                                           
          I = LAND_INDEX(L)                                             
          RADNET_TILE(L,N) = SW_TILE(L,N) +                             
     &                       LW_DOWN(I) - SBCON*TSTAR_TILE(L,N)**4      
          LE_TILE(L,N) = LC*ECAN_TILE(L,N) + LC*ESOIL_TILE(L,N) +       
     &                   LC*ELAKE_TILE(L,N) + LS*EI_TILE(L,N)           
          SURF_HT_FLUX_LAND(I) = SURF_HT_FLUX_LAND(I)               
     &                      + TILE_FRAC(L,N) *                          
     &                      ( RADNET_TILE(L,N) - FTL_TILE(L,N) -        
     &                        LE_TILE(L,N) - LF*MELT_TILE(L,N) -        
     &                       (CANHC_TILE(L,N)/TIMESTEP) *               
     &                       (TSTAR_TILE(L,N) - TSTAR_TILE_OLD(L,N)) )  
          TSTAR_LAND(I) = TSTAR_LAND(I)                             
     &               + TILE_FRAC(L,N)*TSTAR_TILE(L,N)                   
          TSTAR_RAD4(I) = TSTAR_RAD4(I) + FLANDG(I)*              
     &                    TILE_FRAC(L,N)*TSTAR_TILE(L,N)**4             
        ENDDO                                                           
      ENDDO                                                             
                                                                        
! TOA outward LW radiation after boundary layer                         
      DO I=P1,P1+P_POINTS-1                                             
        OLR(I) = OLR(I) + SBCON*TSTAR_RAD4(I)                           
      ENDDO                                                             
                                                                        
!-----------------------------------------------------------------------
! Optional error check : test for negative surface temperature          
!-----------------------------------------------------------------------
      IF (L_NEG_TSTAR) THEN                                             
        DO L=LAND1,LAND1+LAND_PTS-1                                     
          I = LAND_INDEX(L)                                             
          IF (TSTAR_LAND(I).LT.0) THEN                                  
            ERROR = 1                                                   
            WRITE(6,*) '*** ERROR DETECTED BY ROUTINE BDY_LAYR ***'     
            WRITE(6,*) 'NEGATIVE SURFACE TEMPERATURE AT LAND POINT ',L  
          ENDIF                                                         
        ENDDO                                                           
      ENDIF                                                             
                                                                        
!-----------------------------------------------------------------------
! Update sea-ice surface layer temperature.                             
!-----------------------------------------------------------------------
                                                                        
      CALL SICE_HTF(                                                    
     & P_POINTS,P_FIELD,P1,FLANDG,SIMLT,                             
     & DI,ICE_FRACT,SURF_HT_FLUX_SICE,                                  
     & TSTAR_SEA,TSTAR_SICE,TIMESTEP,                                   
     & TI,SICE_MLT_HTF,SEA_ICE_HTF,                                     
     & LTIMER)                                                          
                                                                        
!---------------------------------------------------------------------- 
!! 8.1 Update U_V.                                                      
!---------------------------------------------------------------------- 
                                                                        
! U component of 10m wind                                               
      IF (SU10)THEN                                                     
*IF -DEF,SCMA                                                           
        DO I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1                     
*ELSE                                                                   
        DO I=1,U_POINTS                                                 
*ENDIF                                                                  
          U10M(I) = (U_1(I) + (DU_1(I) - CQ_CM_1(I)*TAUX_1(I)) -        
     &                      U_0(I))*CDR10M_UV(I) + U_0(I)               
        ENDDO                                                           
      ENDIF                                                             
                                                                        
! V component of 10m wind                                               
      IF (SV10)THEN                                                     
*IF -DEF,SCMA                                                           
        DO I=U1+ROW_LENGTH,U1+U_POINTS-ROW_LENGTH-1                     
*ELSE                                                                   
        DO I=1,U_POINTS                                                 
*ENDIF                                                                  
          V10M(I) = (V_1(I) + (DV_1(I) - CQ_CM_1(I)*TAUY_1(I)) -        
     &                      V_0(I))*CDR10M_UV(I) + V_0(I)               
        ENDDO                                                           
      ENDIF                                                             
                                                                        
!-----------------------------------------------------------------------
!! 9.  Calculate surface latent heat flux.                              
!-----------------------------------------------------------------------
                                                                        
      IF (SLH) THEN                                                     
        DO I=P1,P1+P_POINTS-1                                           
          LATENT_HEAT(I) = LC*FQW_1(I)                              
     &     + LF*(FLANDG(I)*EI_LAND(I)+(1.-FLANDG(I))*EI_SICE(I))
        ENDDO                                                           
      ENDIF                                                             
                                                                        

      DO I=P1,P1+P_POINTS-1                                           
        TSTAR(I)=FLANDG(I)*TSTAR_LAND(I)                      
     &    +(1.-FLANDG(I))*TSTAR_SSI(I)                          
                                                                        
        EI(I)=FLANDG(I)*EI_LAND(I)                            
     &    +(1.-FLANDG(I))*EI_SICE(I)                            
                                                                        
        SURF_HT_FLUX(I)=FLANDG(I)*SURF_HT_FLUX_LAND(I)        
     &    +(1.-FLANDG(I))*SURF_HT_FLUX_SICE(I)                  
      ENDDO                                                             


      IF (LTIMER) THEN                                                  
        CALL TIMER('BDYLAYR ',4)                                        
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
*ENDIF                                                                  
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  namelist_famous_land_060912.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/=============================================
*/Mod to enable a perturbed physics ensemble
*/of terrestrial carbon cycle parameters to
*/be performed. It is based on mods from 
*/Ben Booth but has been adapted to take 
*/account of the necessary different coding 
*/standard needed for FAMOUS (UM version 4.5.3)  
*/versus HadCM3 (Um version 4.7). Comments
*/within this mod refer to the line below the
*/comment.
*/
*/Jonny Williams, April 2011
*/
*/=============================================
*/
*ID ID_LAND_CC
*/
*DECLARE TRIF
*/
*/Delete variable definitions of F0,LAI_MIN,NL0,
*/R_GROW, TLOW and TUPP
*D TRIF.30
*D TRIF.43
*D TRIF.44
*D TRIF.45
*D TRIF.54
*D TRIF.59
*D TRIF.60
*D TRIF.61
*D TRIF.62
*/
*/Delete the numerical definition of F0
*D TRIF.76
*/Delete the numerical definition of LAI_MIN
*D MOSES_CTILE.4912
*/Delete the numerical definition of NL0
*D MOSES_CTILE.4913
*/Delete the numerical definition of R_GROW
*D TRIF.91
*/Delete the numerical definition of TLOW
*D ABX1F405.1738
*/Delete the numerical definition of TUPP
*D ABX1F405.1739
*/
*/Define the common deck of variables to be perturbed.
*COMDECK C_LAND_CC
      REAL F0(NPFT),LAI_MIN(NPFT),NL0(NPFT),R_GROW(NPFT),
     +TLOW(NPFT),TUPP(NPFT),Q10,V_CRIT_ALPHA,KAPS
C     
      COMMON /COMMON_LAND_CC/ F0,LAI_MIN,NL0,R_GROW,
     +TLOW,TUPP,Q10,V_CRIT_ALPHA,KAPS
*/    
*/Some variables in the EDDY namelist in CNTLOCN are being
*/changed also so read them in too. 
*COMDECK C_EDDY
      REAL FNUB_SI,KAPPA0_SI,DKAPPA_DZ_SI,
     +FNU0_SI,STABLM_SI,AHI1_SI,AHI2_SI,AHI3_SI,SLOPE_MAX,
     +ATHKDF1_SI,ATHKDF2_SI,ATHKDF3_SI,AM0_SI,AM1_SI,
     +AH_SI,BM_SI,CRIT_RI,MAX_QLARGE_DEPTH  
C
      COMMON /COMMON_EDDY/ FNUB_SI,KAPPA0_SI,DKAPPA_DZ_SI,
     +FNU0_SI,STABLM_SI,AHI1_SI,AHI2_SI,AHI3_SI,SLOPE_MAX,
     +ATHKDF1_SI,ATHKDF2_SI,ATHKDF3_SI,AM0_SI,AM1_SI,
     +AH_SI,BM_SI,CRIT_RI,MAX_QLARGE_DEPTH      
*/
*/Globally recall the variables which were deleted
*/from the TRIF deck
*I TRIF.63
*CALL C_LAND_CC
C
*/
*DECLARE MICROB7A
*I MICROB7A.56
*/NSTYPES defines the number of Plant Functional Types (PFTs)
*/as well as various other surface-specific quantities
*CALL NSTYPES
C
*/This calls the previously defined common block
*CALL C_LAND_CC
C
*/The following 6 lines Remove Q10 and KAPS 
*/from the microbiology deck. In Ben Booth's original mod
*/it is only Q10 which is deleted in this way. I can't
*/work out why this is but this solution is necessary 
*/due to the inclusion of the definition of KAPS
*/within the common deck
*D MICROB7A.63
*D MICROB7A.64
*D MICROB7A.65
*D MICROB7A.66
*D MICROB7A.67
*D MICROB7A.68
*/	    
*/Reads in control information from NAMELIST, and
*/copies into dump header/history file for reference.
*DECLARE READLSA1
*I READLSA1.88
C
*CALL NSTYPES
C
*CALL C_LAND_CC
C
      NAMELIST /LAND_CC/ F0,LAI_MIN,NL0,R_GROW,
     +TLOW,TUPP,Q10,V_CRIT_ALPHA,KAPS
C 
*CALL C_EDDY
C
      NAMELIST /EDDY/ FNUB_SI,KAPPA0_SI,DKAPPA_DZ_SI,
     +FNU0_SI,STABLM_SI,AHI1_SI,AHI2_SI,AHI3_SI,SLOPE_MAX,
     +ATHKDF1_SI,ATHKDF2_SI,ATHKDF3_SI,AM0_SI,AM1_SI,
     +AH_SI,BM_SI,CRIT_RI,MAX_QLARGE_DEPTH 
C     
*/
*B READLSA1.132
      WRITE(6,*) 'Print variables from namelists'
      REWIND 5
      READ(5,LAND_CC)
      REWIND 5
      READ(5,EDDY)
      REWIND 5
      WRITE(6,16) 'F0:       ',F0(1),F0(3),F0(3),F0(4),F0(5)
      WRITE(6,16) 'LAI_MIN:  ',LAI_MIN(1),LAI_MIN(2),LAI_MIN(3)
     +                        ,LAI_MIN(4),LAI_MIN(5)
      WRITE(6,16) 'NL0:      ',NL0(1),NL0(2),NL0(3),NL0(4),NL0(5)
      WRITE(6,16) 'R_GROW:   ',R_GROW(1),R_GROW(2),R_GROW(3)
     +                        ,R_GROW(4),R_GROW(5)
      WRITE(6,16) 'TLOW:     ',TLOW(1),TLOW(2),TLOW(3)
     +                        ,TLOW(4),TLOW(5)
      WRITE(6,16) 'TUPP:     ',TUPP(1),TUPP(2),TUPP(3)
     +                        ,TUPP(4),TUPP(5)
      WRITE(6,22) 'Q10:      ',Q10
      WRITE(6,22) 'V_CRIT_ALPHA:   ',V_CRIT_ALPHA
      WRITE(6,21) 'KAPS:     ',KAPS
      WRITE(6,23) 'RHCRIT:     ',RHCRIT(1),RHCRIT(2)
     +                             ,RHCRIT(3),RHCRIT(4)
     +                             ,RHCRIT(5),RHCRIT(6)
     +                             ,RHCRIT(7),RHCRIT(8)
     +                             ,RHCRIT(9),RHCRIT(10)
     +                             ,RHCRIT(11)
      WRITE(6,22) 'VF1:     ',VF1
      WRITE(6,21) 'CT:     ',CT
      WRITE(6,21) 'CW_SEA:     ',CW_SEA
      WRITE(6,21) 'CW_LAND:     ',CW_LAND
      WRITE(6,21) 'KAY_GWAVE:     ',KAY_GWAVE
      WRITE(6,21) 'KAY_LEE_GWAVE:     ',KAY_LEE_GWAVE
      WRITE(6,21) 'Z0FSEA:     ',Z0FSEA
      WRITE(6,22) 'ALPHAM:     ',ALPHAM
      WRITE(6,24) 'DIFF_COEFF:     ',DIFF_COEFF(1)
     +                              ,DIFF_COEFF(2)
     +                              ,DIFF_COEFF(3)
     +                              ,DIFF_COEFF(4)
     +                              ,DIFF_COEFF(5)
     +                              ,DIFF_COEFF(6)
     +                              ,DIFF_COEFF(7)
     +                              ,DIFF_COEFF(8)
     +                              ,DIFF_COEFF(9)
     +                              ,DIFF_COEFF(10)
     +                              ,DIFF_COEFF(11)      
      WRITE(6,24) 'DIFF_COEFF_Q:   ',DIFF_COEFF_Q(1)
     +                              ,DIFF_COEFF_Q(2)
     +                              ,DIFF_COEFF_Q(3)
     +                              ,DIFF_COEFF_Q(4)
     +                              ,DIFF_COEFF_Q(5)
     +                              ,DIFF_COEFF_Q(6)
     +                              ,DIFF_COEFF_Q(7)
     +                              ,DIFF_COEFF_Q(8)
     +                              ,DIFF_COEFF_Q(9)
     +                              ,DIFF_COEFF_Q(10)
     +                              ,DIFF_COEFF_Q(11)
      WRITE(6,21) 'FNUB_SI:     ',FNUB_SI	
      WRITE(6,21) 'KAPPA0_SI:     ',KAPPA0_SI	
      WRITE(6,21) 'DKAPPA_DZ_SI:     ',DKAPPA_DZ_SI	
      WRITE(6,21) 'FNU0_SI:     ',FNU0_SI	
      WRITE(6,21) 'STABLM_SI:     ',STABLM_SI      
      WRITE(6,21) 'AHI1_SI:     ',AHI1_SI	
      WRITE(6,21) 'AHI2_SI:     ',AHI2_SI	
      WRITE(6,21) 'AHI3_SI:     ',AHI3_SI	
      WRITE(6,21) 'SLOPE_MAX:     ',SLOPE_MAX
      WRITE(6,21) 'ATHKDF1_SI:     ',ATHKDF1_SI	
      WRITE(6,21) 'ATHKDF2_SI:     ',ATHKDF2_SI	
      WRITE(6,21) 'ATHKDF3_SI:     ',ATHKDF3_SI	
      WRITE(6,21) 'AM0_SI:     ',AM0_SI	
      WRITE(6,21) 'AM1_SI:     ',AM1_SI	
      WRITE(6,21) 'AH_SI:     ',AH_SI	
      WRITE(6,21) 'BM_SI:     ',BM_SI	
      WRITE(6,21) 'CRIT_RI:     ',CRIT_RI	
      WRITE(6,21) 'MAX_QLARGE_DEPTH:     ',MAX_QLARGE_DEPTH      
C
16    format(A13,5F9.3)
21    format(A13,e12.4)
22    format(A20,F5.3) 
23    format(A13,11F9.3)
24    format(A13,11e12.4)
C
*/
*/Introduce V_CRIT_ALPHA, V_CRIT_NEW into the plant physiology deck
*DECLARE PHYSIO7A
*/Introduce the perturbed volumetric soil moisture 
*/concentration above which stomata are not sensitive 
*/to soil water (m3 H2O/m3 soil).
*/
*D MOSES_CTILE.2393 
*B MOSES_CTILE.2394
     &,                   V_CRIT_NEW,V_SAT,V_WILT,WIND,Z0_TILE,Z1
*D PHYSIO7A.80
*B PHYSIO7A.81
     &,V_CRIT_NEW(LAND_FIELD)
*I PHYSIO7A.135
*CALL C_LAND_CC
*I PHYSIO7A.142
*/
*/Define the new critical soil moisture concentration 
      DO L=1,LAND_FIELD
         V_CRIT_NEW(L) = V_WILT(L) + 
     &      V_CRIT_ALPHA * (V_SAT(L) - V_WILT(L))
      ENDDO
*/Remove reference to old V_CRIT and replace with 
*/reference to V_CRIT_NEW
*D MOSES_CTILE.2435
*D MOSES_CTILE.2436
        IF (V_CRIT_NEW(L).GT.0.)  
     &    GSOIL(L) = 
     &    GS_NVG(SOIL-NPFT)*(STHU(L,1)*V_SAT(L)/V_CRIT_NEW(L))**2
*/Remove reference to old V_CRIT and replace with 
*/reference to V_CRIT_NEW
*D MOSES_CTILE.2449
     &,               F_ROOT,STHU,V_CRIT_NEW,V_SAT,V_WILT
*/
*/Change calls to V_CRIT to V_CRIT_NEW in SMC_EXT
*/subroutine for completeness
*/
*DECLARE SMCEXT7A
*/
*D MOSES_CTILE.4214
*B SMCEXT7A.25
     &,                   F_ROOT,STHU,V_CRIT_NEW,V_SAT,V_WILT
*/
*D SMCEXT7A.71
*B SMCEXT7A.72
     &,V_CRIT_NEW(NPNTS)
*/
*D SMCEXT7A.120
*B SMCEXT7A.121
     &               /(V_CRIT_NEW(I)-V_WILT(I))     
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  namelist_famous_ocean_cc.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  oco2_cap_060912.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/NOT GREAT - cap pco2 and ocean co2 flux calcs that can (independently)
*/blow up easily and screw up ocean TCO2, and, in coupled carbon runs,
*/everything else besides. In practice doesn't appear to affect even some
*/pretty big perturbations to a modern climate, but far from ideal
*/
*DECLARE FLUX_CO2
*B FLUX_CO2.197
      DO I = IFROM_CYC, ITO_CYC
        if (co2_flux(i).ne.co2_flux(i)) then
!          write(6,*)"f_co2: NaN!",piston_vel(I),PCO2(I)
          co2_flux(i)=0.
        endif
        if (abs(co2_flux(i)).gt.30) then
!          write(6,*)"f_co2 cap",piston_vel(I),PCO2(I)
          co2_flux(i)=sign(30.,co2_flux(i))
        end if
      ENDDO
*DECLARE PPCO2
*B PPCO2.139
      if (PCO2(I).gt.1000) then
!        write(6,*)"ppco2 blown up:",PCO2(I)
        PCO2(I)=1000.
      endif
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  ozone_mod-3levsane
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID CONSTANT O3
*DECLARE RAD_CTL1
*I ADB2F404.937
      DO ROW=1,P_ROWS
        DO I=1,ROW_LENGTH
          POINT=I+(ROW-1)*ROW_LENGTH
          DO LEVEL=1,OZONE_LEVELS
            IF(LEVEL.LT.TRINDX(POINT)) OZONE_1(POINT,LEVEL)=2.0E-8
            IF(LEVEL.EQ.TRINDX(POINT)) OZONE_1(POINT,LEVEL)=1.0E-7
            IF(LEVEL.GT.TRINDX(POINT)) OZONE_1(POINT,LEVEL)=2.0E-6

            IF(LEVEL.EQ.OZONE_LEVELS)  OZONE_1(POINT,LEVEL)=6.0E-6
          ENDDO
        ENDDO
      ENDDO
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  rayleigh_fric.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*DECLARE CRUNTIMC
*// declare array, put into common block
*B CRUNTIMC.52
     & ,RAY_TAU(MAX_P_LEVELS)
*B CRUNTIMC.63
     & ,RAY_TAU
*DECLARE READLSA1
*// make part of namelist, initialise. Will be read in with namelist
*I AJX3F405.151
     &  ,RAY_TAU
*I ADR1F305.178
        RAY_TAU(LEVEL) = 0.
*DECLARE ATMDYN1
*// common block declared here already, pass to difctl as arg
*I APB0F401.78
     &          RAY_TAU,
*DECLARE DIFCTL1A
*// pick up array, use where it has valid values
*I APB0F401.1446
     &                   RAY_TAU,
*I DIFCTL1A.74
     &,RAY_TAU(P_LEVELS)
*I DIFCTL1A.337
! Hack in simple Rayleigh friction to supplment GWD in FAMOUS
          if (abs(ray_tau(k)-0.) .gt. 1e-6) then 
            U(I,K)=U(I,K)-(U(I,K)*ADVECTION_TIMESTEP
     &                           /RAY_TAU(k)/86400.)
            V(I,K)=V(I,K)-(V(I,K)*ADVECTION_TIMESTEP
     &                           /RAY_TAU(K)/86400.)
          end if
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  snowonice_moses2.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*DECLARE FTSA1A
*B ADB1F400.10
     &     S_SEA,
*I AWI2F400.5
     &     ,S_SEA(L1)                         ! Snow amount on sea-ice
*D AJG1F405.62
               if (s_sea(j).gt.0.0) then
*D AJG1F405.69
     &           +(snow_albedo-alphab)*(1.0-exp(-maskd*s_sea(j)))
*DECLARE RAD_CTL1
*B ADB1F400.89
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
*I AJS1F400.175
     &  JSNODEP_SEA,
