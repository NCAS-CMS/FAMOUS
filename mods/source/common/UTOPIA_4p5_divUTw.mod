*ID UTOPIA
*/
*/ "backported" for 4.5 FAMOUS from Dave Storkey's UTOPIA_basic_bit_ys_4p7.m*
*/ All the sequence numbers have changed in the header section, and a 
*/ version of TOTVEL has been scavenged from a 5.3 mod
*/ r.s.smith@reading.ac.uk
*/
*/ <start of header UTOPIA_basic_bit_ys_4p7.mh>
*/
*/ Mod to implement UTOPIA advection scheme (with only first order
*/ cross terms).
*/ Replaces ADVSRCE.dk with ADVECT.dk
*/ Also implements ULTIMATE flux limiter.
*/
*/ Written for UM version 4.5
*/ J+1 variables needed in blokinit, (j=j_2-1): DYU(J-1),DYT2R(J)
*/ J+1 variables needed in tracer, : DYT2R(J+1),CSTR(J+2)
*/
*/ D. Storkey
*/ 22/3/99
*/
*/ Updated 25/01/01 - I had some problems with bit-reproducibility
*/ on different numbers of PE's - but only after more than a day or
*/ so. I've had to slightly change the initialisation of FVST in
*/ blokinit, and I'm pretty sure it compares over more than 2 days!
*/ Malcolm Roberts
*/******************************************************************
*DECLARE COCAWRKA
*/******************************************************************
*/
*/ *D COCAWRKA.6 ,COCAWRKA.7 - these lines inactive. Guess...
*D COCAWRKA.8
*D OOM3F405.21
     &,UCLIN,VCLIN,USAV,VSAV,RHON,RHOS,FUWTP,FUW,FUWT,FVNTP,FVN
     &,FVNT,FVSU,FVSTP,FVST,FMM,FM,FMP,GM,FMPP,UCLINP,VCLINP,D2U,D2V
     &,SSFUBPP,SSFVBPP            
*/
*/ *D COCAWRKA.9  - this line inactive. Guess...
*D ORH3F403.259
     &,UOVER,UDIF,UUNDER,VOVER,VDIF,VUNDER,WTP,W,WT,TDIF  
*/
*/******************************************************************
*DECLARE COCTWRKA
*/******************************************************************
*/
*/ *D COCTWRKA.20 ,COCTWRKA.21  - seems like not the right line!
*D COCTWRKA.18
*D COCTWRKA.19
     * RHON (IMT,KM),RHOS (IMT,KM),
     * FUWTP (IMTP1,KM),FUW  (IMTP1,KM),FUWT  (IMTP1,KM),
     * FVNTP (IMT,KM),FVN  (IMT,KM),FVNT  (IMT,KM), 
     * FVSU (IMT,KM),FVSTP (IMT,KM),FVST (IMT,KM), 
*/
*/ *D COCTWRKA.30  - again, seems wrong
*D ORH1F304.71
     * WTP(IMT,KMP1),W(IMT,KMP1),WT(IMT,KMP1),
     * TDIF(IMT,KMP2,NTMIN2)     
*/
*DECLARE ARGOCTOP
*/ *I ARGOCTOP.44  - this line inactive already. Guess...
*I OOM3F405.32
     &, O_SPCON(jocp_dyujm),O_SPCON(jocp_dyt2rjm)
     &, O_SPCON(jocp_cstrjp),O_SPCON(jocp_cstrjm),O_SPCON(jocp_cstjm)
*DECLARE TYPOCDPT
*/ *I TYPOCDPT.98  - wrong number?
*I OOM3F405.41
     &, jocp_dyujm
     &, jocp_dyt2rjm
     &, jocp_cstrjp
     &, jocp_cstrjm
     &, jocp_cstjm
*DECLARE COMOCDPT
*/ *I COMOCDPT.97  - wrong number
*I OOM3F405.50
     &, jocp_dyujm
     &, jocp_dyt2rjm
     &, jocp_cstrjp
     &, jocp_cstrjm
     &, jocp_cstjm
*DECLARE ARGOCONE
*/ *D ARGOCONE.8  - line numbering seems unlikely
*D OOM3F405.73
     & CSRJP,DYU2RJP,CSTJP,DYTRJP,CSJM,DYURJM,DYUJM,DYT2RJM,CSTRJP,
     & CSTRJM,CSTJM,
*DECLARE TYPOCONE
*/ *I TYPOCONE.88 - line doesn't exist anymore?
*I OOM3F405.79
     &,      DYUJM,   ! DYU for J-1 outside halo
     &      DYT2RJM,   ! DYT2R for J-1 outside halo
     &      CSTRJP,   ! CSTR for J+2 outside halo
     &      CSTRJM,   ! CSTR for J-1 outside halo
     &      CSTJM    ! CST  for J-1 outside halo
*/
*/ <end of header UTOPIA_basic_bit_ys_4p7.mh>
*/ <start of mod UTOPIA_basic_bit_ys_4p7.mf77>
*/
*DECLARE OCNARYPT
*I OOM3F405.64
       jocp_dyujm=icount
       icount=icount+1
       jocp_dyt2rjm=icount
       icount=icount+1
       jocp_cstrjp=icount
       icount=icount+1
       jocp_cstrjm=icount
       icount=icount+1
       jocp_cstjm=icount
       icount=icount+1
*DECLARE OSETCON
*I OOM3F405.209
C Send value of DYU(J_JMT) forward to next PE, to be DYUJM 
      CALL GC_GSYNC(O_NPROC,INFO)  
      IF (J_PE_JFINP1.GE.0) THEN   
c       ! We must send row J_JMT-1    
         PE_RECV=J_PE_JFINP1                     
         CALL GC_RSEND(5027,1,PE_RECV,INFO,DYUJM,DYU(J_JMT-1))  
      ENDIF                                      
                                                  
      CALL GC_GSYNC(O_NPROC,INFO)                
                                                
      IF (J_PE_JSTM1.GE.0) THEN                
c       ! We're expecting to receive a message:  
         PE_SEND = J_PE_JSTM1                     
         CALL GC_RRECV(5027,1,PE_SEND,INFO,DYUJM,DYU) 
      ENDIF                                       
                                                   
      CALL GC_GSYNC(O_NPROC,INFO)                  
                                                    
C Send value of DYT2R(J_JMT) forward to next PE, to be DYT2RJM 
      CALL GC_GSYNC(O_NPROC,INFO)  
      IF (J_PE_JFINP1.GE.0) THEN   
c       ! We must send row J_JMT-1    
         PE_RECV=J_PE_JFINP1                     
         CALL GC_RSEND(5028,1,PE_RECV,INFO,DYT2RJM,DYT2R(J_JMT-1))  
      ENDIF                                      
                                                  
      CALL GC_GSYNC(O_NPROC,INFO)                
                                                
      IF (J_PE_JSTM1.GE.0) THEN                
c       ! We're expecting to receive a message:  
         PE_SEND = J_PE_JSTM1                     
         CALL GC_RRECV(5028,1,PE_SEND,INFO,DYT2RJM,DYT2R) 
      ENDIF                                       
                                                   
      CALL GC_GSYNC(O_NPROC,INFO)   
               
C Send value of CSTR(J_1+1) back to previous PE, to be CSTRJP 
      CALL GC_GSYNC(O_NPROC,INFO)               
      IF (J_PE_JSTM1.GE.0) THEN               
c       ! We must send row J_1+1             
         PE_RECV=J_PE_JSTM1                  
         CALL GC_RSEND(5029,1,PE_RECV,INFO,CSTRJP,CSTR(J_1+1))
      ENDIF                                           
                                                     
      CALL GC_GSYNC(O_NPROC,INFO)                     
                                                       
      IF (J_PE_JFINP1.GE.0) THEN                     
c       ! We're expecting to receive a message:     
         PE_SEND = J_PE_JFINP1                       
         CALL GC_RRECV(5029,1,PE_SEND,INFO,CSTRJP,CSTR)   
      ENDIF                                            
                                                     
      CALL GC_GSYNC(O_NPROC,INFO)                        

C Send value of CSTR(J_JMT) forward to next PE, to be CSTRJM 
      CALL GC_GSYNC(O_NPROC,INFO)  
      IF (J_PE_JFINP1.GE.0) THEN   
c       ! We must send row J_JMT-1    
         PE_RECV=J_PE_JFINP1                     
         CALL GC_RSEND(5030,1,PE_RECV,INFO,CSTRJM,CSTR(J_JMT-1))  
      ENDIF                                      
                                                  
      CALL GC_GSYNC(O_NPROC,INFO)                
                                                
      IF (J_PE_JSTM1.GE.0) THEN                
c       ! We're expecting to receive a message:  
         PE_SEND = J_PE_JSTM1                     
         CALL GC_RRECV(5030,1,PE_SEND,INFO,CSTRJM,CSTR) 
      ENDIF                                       
                                                   
      CALL GC_GSYNC(O_NPROC,INFO)   
               
C Send value of CST(J_JMT) forward to next PE, to be CSTJM 
      CALL GC_GSYNC(O_NPROC,INFO)  
      IF (J_PE_JFINP1.GE.0) THEN   
c       ! We must send row J_JMT-1    
         PE_RECV=J_PE_JFINP1                     
         CALL GC_RSEND(5031,1,PE_RECV,INFO,CSTJM,CST(J_JMT-1))  
      ENDIF                                      
                                                  
      CALL GC_GSYNC(O_NPROC,INFO)                
                                                
      IF (J_PE_JSTM1.GE.0) THEN                
c       ! We're expecting to receive a message:  
         PE_SEND = J_PE_JSTM1                     
         CALL GC_RRECV(5031,1,PE_SEND,INFO,CSTJM,CST) 
      ENDIF                                       
                                                   
      CALL GC_GSYNC(O_NPROC,INFO)   

*I ORH4F402.25
      CALL SWAPBOUNDS(DYT2R,1,JMT,O_EW_HALO,O_NS_HALO,1)
*/******************************************************************
*/ this is done in the rho_av mod for pressure averaging
*DECLARE ROWCALC
*/******************************************************************
*I ORH1F305.2318
C-----------------------------------------------------------------------
C Copy velocities at base of model levels to STASH workspace. WT is for
C the present row and is already calculated, either by BLOKINIT or by the
C execution of TRACER for the previous row.
C-----------------------------------------------------------------------
C
      IF (sf201_30) THEN
        DO K=2,KM
          DO I=IFROM_CYC,ICOL_CYC
            sw201_30(I,J,K-1)=WT(I,K)*FM(I,K)+(1.0-FM(I,K))*RMDI
          END DO
          IF (L_OCYCLIC) THEN
              sw201_30(1,J,K-1)=WT(IMTM1,K)*FM(IMTM1,K)+
     *                                    (1.0-FM(IMTM1,K))*RMDI
          ENDIF
        END DO
      END IF
*D ROWCALC.857,ROWCALC.877
*/ *D ROWCALC.870
*/             sw201_30(I,J,K-1)=WT(I,K)*FM(I,K)+(1.0-FM(I,K))*RMDI
*/ *D ORH1F305.2354
*/               sw201_30(1,J,K-1)=WT(IMTM1,K)*FM(IMTM1,K)+
*I ORH1F305.2456
          WT(1  ,K)=WT(IMTM1,K)
          WT(IMT,K)=WT(2    ,K)
*D ORH4F404.85
       CALL MATRIX(WT,IMT,ISTRT,ISTOP,0,KMP1,SCL,0,0)
*/******************************************************************
*DECLARE TRACER
*/******************************************************************
*I OJG2F401.19
       LOGICAL l_xadv,l_yadv,l_zadv   ! advection diags from 
*/
*/
*/ Delete some commented out lines!
*D OSY1F405.86,87
*/
*/ Remove section that calculates FVST for the (globally) southernmost
*/ row. This now gets done in BLOKINIT.
*/
*D OSI1F405.20,OSI1F405.21
C 
C Section to calculate FVST for row J+J_OFFSET=2 is now in BLOKINIT.
C
*/
*/ Change so that code calculates FUWTP and FVNTP rather than FUW and 
*/ FVN. Quantities will be passed back as slab calculation progresses.
*/ Note also setting quantities near the EW grid boundaries more carefully,
*/ so the flux limiter works properly.
*/
*D OSY1F405.84,OSY1F405.85
        FXA=DYT2R(J+1)
*D TRACER.310,311
      DO 691 I=2,IMTP1 
        FUWTP(I,K)=(UP(I-1,K)*DYU (J+1)+U(I-1,K)*DYU (J  ))*FXA  
*D TRACER.315
      IF(L_OCYCLIC) THEN
        FUWTP(1,K)=(UP(IMTM2,K)*DYU (J+1)+U(IMTM2,K)*DYU (J  ))*FXA  
      ELSE
        FUWTP(1,K)=0.0
      ENDIF 
*D OSY1F405.88
      FXB=(CSTR(J+1)*CS(J+1))
*D OSY1F405.91
               FVNTP(I,K)=((VP(I,K)*DXU(I))+(VP(I-1,K)*DXU(I-1)))*FXB 
*D OSY1F405.94
           IF(L_OCYCLIC) THEN
             FVNTP(1,K)=((VP(1,K)*DXU(1))+(VP(IMTM2,K)*DXU(IMTM2)))*FXB 
     &                                   *DXT2R(1)  
           ELSE
             FVNTP(1,K)=0.0
           ENDIF            
*/
*/ Calculate WTP rather than W (via continuity) and make sure it's  
*/ calculated in column IMT rather than just being set to zero.
*/
*D ORL1F404.648
          WTP(I,KMP1) = FX
*D ORL1F404.651
          WTP(I,1) = FX 
*D ORH1F305.821,OSY1F405.97
             DO I=1,IMT
        WTP(I,K+KOFF)=DZ(K)
     &             *( (FUWTP(I+1,K)-FUWTP (I,K))*(DXTR(I)*CSTR(J+1))
     &               +(FVNTP(I  ,K)-FVSTP (I,K))*DYTR(J+1))  
*D ORL1F404.665
*D ORL1F404.674
             WTP(I,K-1)=WTP(I,K)-WTP(I,K-1) 
*D ORH1F305.838
               WTP(I,K+1)=WTP(I,K)+WTP(I,K+1)
*/
*/ Put the code that sets W=0 on bathymetry into the "calculation of W"
*/ section rather than the "set boundary conditions for tracer diffusion"
*/ section!
*/
*B TRACER.363

! Eliminate rounding error by setting WTP=0 at the bottom.

      DO I=1,IMT                                   
        KZ = KMTP(I)
        IF(KZ.GT.0) WTP(I,KZ+1)=0.
      ENDDO
*/
*D TRACER.369
*D TRACER.384
*D TRACER.394
*/
*/ Need to change the call to ADV_SOURCE (now called "ADVECT"), and 
*/ change the argument list. Also make sure the TA array is initialised.
*/
*D OSY1F405.104,124
      DO K=1,KM
        DO I=1,IMT
          TA(I,K,M) = 0.0
        ENDDO
      ENDDO

      if (M.eq.1) then
        l_xadv=(sf_dt(1).or.sf_dt(12))
        l_yadv=(sf_dt(2).or.sf_dt(12))
        l_zadv=(sf_dt(3).or.sf_dt(12))
      else if (M.eq.2) then
        l_xadv=(sf_ds(1).or.sf_ds(11))
        l_yadv=(sf_ds(2).or.sf_ds(11))
        l_zadv=(sf_ds(3).or.sf_ds(11))
      endif

*/      
      CALL ADVECT(                                                 
     & O_ADVECT_SCHEME(1,M),                                 
     & J,                                                    
     & IMT,J_JMT,KM,                                         
     & TA(1,1,M),                                            
     & WDTXADV,WDTYADV,TEMPA,WDTZADV,                        
     & T(1,1,M),TB(1,1,M),                                   
     & TM(1,1,M),TBM(1,1,M),TP(1,1,M),TBP(1,1,M),            
     & TPP(1,1,M),TBPP(1,1,M),                                        
     & FUWT,FUWTP,FVNT,FVNTP,FVST,FVSTP,WT,WTP,
     & FLUXST(1,1,M),FLUXNT(1,1,M),                             
     & TIMESTEP(1),DXTR,DYTR,DXUR,DYUR,DZ2R,DZZ2R,CSTR,
     & KMTJM,KMT,KMTP,KMTPP,                                          
     & L_OIMPADDF,                                                    
     & L_OFREESFC,                                                    
     & L_BOOTSTRAP,                                                   
     & L_OCYCLIC,                                                     
     & J_OFFSET,imout,jmout,imout_hud,jmout_hud,ATTEND,HUDTEND,       
     & NMEDLEV,m,NT,L_OMEDADV,L_OHUDOUT,SF_DT(15),SF_DS(14),WDTMED,   
! the following line used to be 2 seperate ones
     & WDSMED,CS,l_xadv,l_yadv,l_zadv
     & )                                                              

! Set cyclic boundary conditions on TA at this stage, since this
! isn't done within ADVECT.

      IF (L_OCYCLIC) THEN
      DO K=1,KM
        TA(1  ,K,M)=TA(IMTM1,K,M)
        TA(IMT,K,M)=TA(2    ,K,M)
      ENDDO
      ENDIF

*/  
*/ Pass lots more variables from slab to slab.
*/
*D OSY1F405.141
*IF DEF,MPP
      IF (J.LT.J_JMT) THEN
        FXA = (CST(J+1)*CSTR(J+2))
      ELSE
        FXA = (CST(J+1)*CSTRJP)
      ENDIF
*ELSE
      FXA = (CST(J+1)*CSTR(J+2))
*ENDIF
      FXB = (CST(J)*CSTR(J+1))
*D TRACER.1503,1505
      DO K=1,KM                       
        DO I=1,IMT                                
          FUWT(I,K) = FUWTP(I,K)
          FVST(I,K) = FVSTP(I,K)
          FVNT(I,K) = FVNTP(I,K)
          FVSTP(I,K) = FVNTP(I,K)*FXA
          WT(I,K)   = WTP(I,K)
*D OSY1F405.143
            FLUXST(I,K,M)=FLUXNT(I,K,M)*FXB
*D TRACER.1507
        ENDDO ! over i
        FUWT(IMTP1,K) = FUWTP(IMTP1,K)
      ENDDO ! over k

      DO I=1,IMT
        WT(I,KMP1)   = WTP(I,KMP1)
      ENDDO
*/
*/ Change "W" to "WT" in calculation of BUOY.
*/
*D TRACER.1167
          BUOY=BUOY-FX*DZZ(K)*WT(I,K)*(RHOS(I,K-1)+RHOS(I,K))  
*/ change advection diagnostics
*D TRACER.1123,TRACER.1124
     *                (WT(I,K+1)*(T(I,K  ,M)+T(I,KP1,M))
     *                -WT(I,K  )*(T(I,KM1,M)+T(I,K  ,M)))*DZ2R(K)
*D TRACER.1128
     *                (WT(I,K+1)*(TB(I,K  ,M)+TB(I,KP1,M)+
*D TRACER.1130
     *                -WT(I,K  )*(TB(I,KM1,M)+TB(I,K  ,M)+ 
*/
*/
*/******************************************************************
*DECLARE BLOKINIT
*/******************************************************************
*/ As the logical structure of BLOKINIT hardly merits the name, 
*/ I'm going to feel free to muck about with it for the sake of not
*/ having too much duplication of code.
*I ORH7F402.317
        REAL FXA1
	REAL DUMMYA(IMT_FSF)
        LOGICAL l_xadv,l_yadv,l_zadv   ! advection diagnostics
*/
*B BLOKINIT.96
     &     KZ,                ! Number of sea levels at point
*B OOM3F405.305
     &,    FXC               ! Temporary value 
*I BLOKINIT.129
     &,    SFUM(IMT)        ! SFUM at J - 1
     &,    SFVM(IMT)        ! SFVM at J - 1
*B ORH6F404.1110
C calculate value of UM = total velocity at j-1 for Utopia scheme
*IF DEF,MPP
        CALL TOTVEL(
*CALL ARGSIZE
     & P(1,J),P(1,J+1),DYU2RJ,DXU2R,HRJ,CSRJ,SFUM,SFVM
     &,DUMMYA,DUMMYA,.false.,L_OCYCLIC)
*ELSE
        CALL TOTVEL(
*CALL ARGSIZE
     & P(1,J),P(1,J+1),DYU2R(J),DXU2R,HR(1,j),CSR(J),SFUM,SFVM
     &,DUMMYA,DUMMYA,.false.,L_OCYCLIC)
*ENDIF
*B ORH6F404.1133
               SFUM(IMT)=SFUM(2)  
               SFVM(IMT)=SFVM(2)
*B ORH6F404.1148
                     UM(I,K)=UM(I,K)+SFUM(I) 
                     VM(I,K)=VM(I,K)+SFVM(I) 
*/
*/ Delete these lines which initialize FLUXST for the southernmost
*/ row. This will now get done in the main section to calculate advective
*/ quantities.
*/ 
*/*D OSY1F405.23,29
*/
*/ End the IF(JST.GT.1) block. The section to calculate the advective
*/ quantities will get done for all blocks with IF tests round certain
*/ lines as appropriate.
*/
*B ORH6F404.1163
      ENDIF  !  jst.gt.1 (see line BLOKINIT.642 if you're confused).
*/
*/ This line (ORH6F404.1164) is actually incorrect in the code, as the 
*/ row-by-row calculation starts at J_2, not J_1. It didn't matter before
*/ since BLOKINIT didn't deal with the southernmost block, and for all 
*/ other blocks, J_1 = J_2. But now BLOKINIT does deal with the southernmost
*/ block, so it does matter.
*/
*/*D ORH6F404.1164
*/         J=J_2-1 

*/
*/ Change comment.
*/
*D BLOKINIT.803
C  BOOTSTRAP SLAB-BY-SLAB CALCULATION OF ADVECTIVE VELOCITIES 
C  AND FLUXES
*D BLOKINIT.805
*/
*B ORH0F405.5
     &     KOFF,             ! offset used in divergence calculation    
*D OSY1F405.31,ORH6F404.1172
*IF DEF,MPP
         J=J_2-1 
*ELSE
         J=JST-1
*ENDIF
c        FXA=(CSTR(J)*CSJM) 
        FXA=(CSTRJM*CSJM)
        FXA1=(CSTJM*CSTR(J)) 
        FXB=(CSTR(J)*CS(J)) 
        FXC=(CSTR(J+1)*CS(J+1)) 
        IF(JST.GT.1) THEN
          DO K=1,KM
            DO I=2,IMT
              FUWT(I,K)=(U(I-1,K)*DYU (J)+UM(I-1,K)*DYUJM)
     &                                                   *DYT2R(J)
              FVST(I,K)=(((VM(I,K)*DXU(I))+(VM(I-1,K)*DXU(I-1)))*FXA   
     &                                   *DXT2R(I))*FXA1        
              FVNT(I,K)=((V(I,K)*DXU(I))+(V(I-1,K)*DXU(I-1)))*FXB
     &                                   *DXT2R(I)        
              FVSTP(I,K)=FVNT(I,K)*(CST(J)*CSTR(J+1))
              FUWTP(I,K)=(UP(I-1,K)*DYU (J+1)+U(I-1,K)*DYU (J  ))
     &                                                   *DYT2R(J+1)
              FVNTP(I,K)=((VP(I,K)*DXU(I))+(VP(I-1,K)*DXU(I-1)))*FXC
     &                                 *DXT2R(I)        
            ENDDO ! on i

            IF(L_OCYCLIC) THEN
              FUWT(1,K)=(U(IMTM2,K)*DYU (J)+UM(IMTM2,K)*DYUJM)
     &                                                     *DYT2R(J)
              FVST(1,K)=(((VM(1,K)*DXU(1))+(VM(IMTM2,K)*DXU(IMTM2)))*FXA   
     &                                     *DXT2R(1))*FXA1    
              FVNT(1,K)=((V(1,K)*DXU(1))+(V(IMTM2,K)*DXU(IMTM2)))*FXB
     &                                     *DXT2R(1)        
              FVSTP(1,K)=FVNT(1,K)*(CST(J)*CSTR(J+1))
              FUWTP(1,K)=(UP(IMTM2,K)*DYU (J+1)+U(IMTM2,K)*DYU (J  ))
     &                                                   *DYT2R(J+1)
              FVNTP(1,K)=((VP(1,K)*DXU(1))+(VP(IMTM2,K)*DXU(IMTM2)))*FXC
     &                                 *DXT2R(1)        
            ELSE ! if not l_ocyclic
              FUWT(1,K)=0.0 
              FUWTP(1,K)=0.0 
              FVST(1,K)=0.0
              FVSTP(1,K)=0.0 
              FVNT(1,K)=0.0 
              FVNTP(1,K)=0.0 
            ENDIF ! l_ocyclic

            FUWT(IMTP1,K)=(U(IMT,K)*DYU (J)+UM(IMT,K)*DYUJM)
     &                                                   *DYT2R(J)
            FUWTP(IMTP1,K)=(UP(IMT,K)*DYU (J+1)+U(IMT,K)*DYU (J  ))
     &                                                   *DYT2R(J+1)
          ENDDO ! on k

        ELSE   !  if(jst.eq.1) then...

          DO K=1,KM
            DO I=2,IMT
              FUWT(I,K)=0.0
              FVST(I,K)=0.0
              FVNT(I,K)=((V(I,K)*DXU(I))+(V(I-1,K)*DXU(I-1)))*FXB
     &                                   *DXT2R(I)        
              FVSTP(I,K)=FVNT(I,K)*(CST(J)*CSTR(J+1))
              FUWTP(I,K)=(UP(I-1,K)*DYU (J+1)+U(I-1,K)*DYU (J  ))
     &                                                   *DYT2R(J+1)
              FVNTP(I,K)=((VP(I,K)*DXU(I))+(VP(I-1,K)*DXU(I-1)))*FXC
     &                                 *DXT2R(I)        
            ENDDO ! on i

            IF(L_OCYCLIC) THEN
              FUWT(1,K)=0.0
              FVST(1,K)=0.0
              FVNT(1,K)=((V(1,K)*DXU(1))+(V(IMTM2,K)*DXU(IMTM2)))*FXB
     &                                   *DXT2R(1)        
              FVSTP(1,K)=FVNT(1,K)*(CST(J)*CSTR(J+1))
              FUWTP(1,K)=(UP(IMTM2,K)*DYU (J+1)+U(IMTM2,K)*DYU (J  ))
     &                                                   *DYT2R(J+1)
              FVNTP(1,K)=((VP(1,K)*DXU(1))+(VP(IMTM2,K)*DXU(IMTM2)))*FXC
     &                                 *DXT2R(1)        
            ELSE ! if not l_ocyclic
              FUWT(1,K)=0.0
              FVST(1,K)=0.0
              FVNT(1,K)=0.0
              FVSTP(1,K)=0.0
              FUWTP(1,K)=0.0 
              FVNTP(1,K)=0.0 
            ENDIF ! l_ocyclic

            FUWT(IMTP1,K)=0.0
            FUWTP(IMTP1,K)=(UP(IMT,K)*DYU (J+1)+U(IMT,K)*DYU (J  ))
     &                                                   *DYT2R(J+1)
          ENDDO ! on k

        ENDIF  !  jst.gt.1
*/
*B OSY1F405.34
      IF(JST.GT.1) THEN
C                                                                       
C---------------------------------------------------------------------    
C  COMPUTE VERTICAL VELOCITY WT IN T COLUMNS                                 
C---------------------------------------------------------------------    
C                                                                         
C  1st, rearrange the calculation of the vertical velocity to allow for   
C  the free surface solution. This method integrates upwards from the     
C  bottom level to calculate the free surface vertical velocity rather    
C  than integrating downwards from the rigid lid boundary condition.      
C                                                                         
      FX = 0.0                                                            
      DO I=1,IMT                                                          
        IF (L_OFREESFC) THEN                                              
          WT(I,KMP1) = FX                                                  
        ELSE                                                              
                                                                          
          WT(I,1) = FX                                                     
        ENDIF         ! L_OFREESFC                                        
      ENDDO           ! over i                                            
     
                                                                       
C                                                                      
C  2ND, CALCULATE THE DIVERGENCE (DEPTH WEIGHTED) AT EACH LEVEL        
C                                                                      
      IF (L_OFREESFC) THEN                                             
        KOFF = 0                                                       
      ELSE                                                             
        KOFF = 1                                                       
      ENDIF                                                            
                                                                       
      DO K=1,KM                                                        
             DO I=1,IMT
        WT(I,K+KOFF)=DZ(K)*((FUWT(I+1,K)-FUWT (I,K))*(DXTR(I)*CSTR(J))    
     &                    +(FVNT(I  ,K)-FVST(I,K))*DYTR(J))             
             ENDDO ! over I                                            
      ENDDO ! over K                                                   
C
C  3RD, INTEGRATE DOWNWARDS (RIGID LID SOLUTION)                       
C                 UPWARDS   (FREE SURFACE SOLUTION)                    
C                                                                      
      IF (L_OFREESFC) THEN                                             
                                                                        
        DO K=KMP1,2,-1                                           
          DO I=1,IMT                                             
             WT(I,K-1)=WT(I,K)-WT(I,K-1)                            
          ENDDO   ! over i                                       
        ENDDO ! over K                                           
      ELSE                                                           
      DO K=1,KM                                            
            DO I=1,IMT                                     
               WT(I,K+1)=WT(I,K)+WT(I,K+1)                    
            ENDDO                                          
C                                                          
      ENDDO ! over K                                       
                                                           
      ENDIF      ! L_OFREESFC                              

! Eliminate rounding error by setting WT=0 at the bottom.

      DO I=1,IMT                                   
        KZ = KMT(I)
        IF(KZ.GT.0) WT(I,KZ+1)=0.
      ENDDO

      ENDIF      ! jst.gt.1
C                                                                       
C---------------------------------------------------------------------    
C  COMPUTE VERTICAL VELOCITY WTP IN T COLUMNS                                 
C---------------------------------------------------------------------    
C                                                                         
C  1st, rearrange the calculation of the vertical velocity to allow for   
C  the free surface solution. This method integrates upwards from the     
C  bottom level to calculate the free surface vertical velocity rather    
C  than integrating downwards from the rigid lid boundary condition.      
C                                                                         
      FX = 0.0                                                            
      DO I=1,IMT                                                          
        IF (L_OFREESFC) THEN                                              
          WTP(I,KMP1) = FX                                                  
        ELSE                                                              
                                                                          
          WTP(I,1) = FX                                                     
        ENDIF         ! L_OFREESFC                                        
      ENDDO           ! over i                                            
     
                                                                       
C                                                                      
C  2ND, CALCULATE THE DIVERGENCE (DEPTH WEIGHTED) AT EACH LEVEL        
C                                                                      
      IF (L_OFREESFC) THEN                                             
        KOFF = 0                                                       
      ELSE                                                             
        KOFF = 1                                                       
      ENDIF                                                            
                                                                       
      DO K=1,KM                                                        
             DO I=1,IMT
        WTP(I,K+KOFF)=DZ(K)
     &            *( (FUWTP(I+1,K)-FUWTP (I,K))*(DXTR(I)*CSTR(J+1))    
     &              +(FVNTP(I  ,K)-FVSTP(I,K))*DYTR(J+1) )             
             ENDDO ! over I                                            
      ENDDO ! over K                                                   
C                                                                      
C  3RD, INTEGRATE DOWNWARDS (RIGID LID SOLUTION)                       
C                 UPWARDS   (FREE SURFACE SOLUTION)                    
C                                                                      
      IF (L_OFREESFC) THEN                                             
                                                                        
        DO K=KMP1,2,-1                                           
          DO I=1,IMT                                             
             WTP(I,K-1)=WTP(I,K)-WTP(I,K-1)                            
          ENDDO   ! over i                                       
        ENDDO ! over K                                           
      ELSE                                                           
      DO K=1,KM                                            
            DO I=1,IMT                                     
               WTP(I,K+1)=WTP(I,K)+WTP(I,K+1)                    
            ENDDO                                          
C                                                          
      ENDDO ! over K                                       
                                                           
      ENDIF      ! L_OFREESFC                              

! Eliminate rounding error by setting WTP=0 at the bottom.

      DO I=1,IMT                                   
        KZ = KMTP(I)
        IF(KZ.GT.0) WTP(I,KZ+1)=0.
      ENDDO
*/
*/ Need to change the call to ADV_SOURCE (now called "ADVECT"), and 
*/ change the argument list. Also correct the preceding comment.
*/
*D OSY1F405.37,40
C  ADVECT to calculate FLUXNT for the halo row, and this then becomes 
C  FLUXST for the first row of the block.
*B OSY1F405.42
      IF(JST.GT.1) THEN
*D OSY1F405.48,68
        l_xadv=.false.
        l_yadv=.false.
        l_zadv=.false.

      CALL ADVECT(                                                 
     & O_ADVECT_SCHEME(1,M),                                 
     & J,                                                    
     & IMT,J_JMT,KM,                                         
     & TEMPA,                                 ! Dummy return   
     & TEMPA,TEMPA,TEMPA,TEMPA,               ! variables
     & T(1,1,M),TB(1,1,M),                                   
     & TM(1,1,M),TBM(1,1,M),TP(1,1,M),TBP(1,1,M),            
     & TPP(1,1,M),TBPP(1,1,M),                                        
     & FUWT,FUWTP,FVNT,FVNTP,FVST,FVSTP,WT,WTP,    
     & FLUXST(1,1,M),FLUXNT(1,1,M),                             
     & C2DTTS,DXTR,DYTR,DXUR,DYUR,DZ2R,DZZ2R,CSTR,
     & KMTJM,KMT,KMTP,KMTPP,                                          
     & L_OIMPADDF,                                                    
     & L_OFREESFC,                                                    
     & L_BOOTSTRAP,                                                   
     & L_OCYCLIC,                                                     
     & J_OFFSET,imout,jmout,imout_hud,jmout_hud,TEMPTEND,TEMPTEND,       
     & NMEDLEV,m,NT,L_OMEDADV,L_OHUDOUT,.FALSE.,.FALSE.,TEMPMED,   
! the following used to be two seperate lines
     & TEMPMED,CS,l_xadv,l_yadv,l_zadv                                                
     & )                                                              
*B BLOKINIT.816
      ELSE   !  if(jst.eq.1) then...
        DO M=1,NT                               
          DO K=1,KM                   
            DO I=1,IMT                
              FLUXNT(I,K,M)=0.        
            ENDDO                     
          ENDDO                       
        ENDDO                         
      ENDIF  !  jst.gt.1
*/
*D OSY1F405.70
*IF DEF,MPP
      IF (J.LT.J_JMT) THEN
      FXA = (CST(J+1)*CSTR(J+2))
      ELSE
      FXA = (CST(J+1)*CSTRJP)
      ENDIF
*ELSE
      FXA = (CST(J+1)*CSTR(J+2))
*ENDIF
      FXB = (CST(J)*CSTR(J+1))
*D ORH6F404.1174,1178
         DO K=1,KM           
            DO I=1,IMT       
               FUWT(I,K)=FUWTP(I,K)
               FVST(I,K)=FVSTP(I,K)
               FVNT(I,K)=FVNTP(I,K)
               FVSTP(I,K)=FVNTP(I,K)*FXA
               WT(I,K)=WTP(I,K)
               DO M=1,NT                              
                 FLUXST(I,K,M)=FLUXNT(I,K,M)*FXB   
               ENDDO                    
            ENDDO    ! Over I
            FUWT(IMTP1,K) = FUWTP(IMTP1,K)
         ENDDO       ! Over K 

         DO I=1,IMT
            WT(I,KMP1)=WTP(I,KMP1)
         ENDDO
*IF DEF,MPP
         J=J_1-1 
*ELSE
         J=JST-1
*ENDIF
*/
*/ Restore the IF(JST.GT.1) statement that I interrupted just before 
*/ this section.
*/
*B ORH6F404.1179
      IF(JST.GT.1) THEN
*/
*/ the above calls for a TOTVEL subroutine: the following is adapted
*/ from ORHAF503.mf77
*/
*I BLOKINIT.866
      SUBROUTINE TOTVEL(
!
!LL
!LL   Author: M J Roberts
!LL
!LL   Description: This brings the calculation of the barotropic
!LL                part of the ocean velocity to one subroutine,
!LL                rather than having it inline in the code in
!LL                BLOKINIT and CLINIC, which can give non-bit
!LL                reproducibility. Included CYCLIC and free
!LL                surface conditions - R. Hill.
!LL
!LL   External documentation:
!LL
!LL   Implemented at UM vn 4.5 30th Nov 1998
!LL
!LL   Modification History  :
!LL   Version   Date     Comment & Name
!LL   ------- --------  --------------------------------------------
!LL
!LL   Subroutine dependencies.
!LL   Called by BLOKINIT and CLINIC
!LL
!LLEND-----------------------------------------------------------------

*CALL ARGSIZE
     &  PB1,PB2,DYU2R,DXU2R,HR,CSR,SFU,SFV
     &, UBTBBC,VBTBBC,L_OFREESFC,L_OCYCLIC)

      IMPLICIT NONE

*CALL OARRYSIZ
*CALL PARPARM
*CALL TYPSIZE
      LOGICAL L_OFREESFC,L_OCYCLIC  ! IN

      REAL UBTBBC(IMT_FSF),VBTBBC(IMT_FSF) ! IN

      INTEGER I
      REAL SFU(IMT),SFV(IMT),PB1(IMT_STREAM),PB2(IMT_STREAM),
     &     DIAG1,DIAG2,HR(IMT),CSR,DYU2R,DXU2R(IMT)

      IF (L_OFREESFC) THEN
        DO I=1,IMTM1
          SFU(I) = UBTBBC(I)
          SFV(I) = VBTBBC(I)
        ENDDO       ! over i
      ELSE
        DO I=1,IMT-1
          DIAG1=PB2(I+1) - PB1(I)
          DIAG2=PB2(I)   - PB1(I+1)
          SFU(I)=-(DIAG1+DIAG2)*DYU2R*HR(I)
          SFV(I)= (DIAG1-DIAG2)*DXU2R(I )*HR(I)*CSR
        ENDDO  ! Over I
      ENDIF

      IF (L_OCYCLIC) THEN
! SET CYCLIC BOUNDARY CONDITIONS
        SFU (IMT)=SFU (2)
        SFV (IMT)=SFV (2)
      ELSE
        SFU (IMT)=0.0
        SFV (IMT)=0.0
      ENDIF

      RETURN
      END
*/******************************************************************
*DECK -ADVSRCE
*/******************************************************************
*/
*/******************************************************************
*DECK ADVECT
*/******************************************************************
*/
      real function LIMITER(
     & f_adv,
     & fu,fc,fd,
     & rcfl,fcurvn,
     & c_flag
     & )
 
      implicit none
                                                                               
!------------------------------------------------------------------- 
!  Declare argument list.                                            
!------------------------------------------------------------------- 

c      real limiter
      real f_adv
      real fc, fu, fd
      real rcfl,fcurvn
      logical c_flag

!------------------------------------------------------------------- 
!  Declare local variables.                                            
!------------------------------------------------------------------- 

      real fnc,fnf,fnfm,fdel,rfdel,fnmax,fnc_in_range,fnf_in_range
      real limited_f,courant_in_range

*CALL COCTOL

!===================================================================
! Begin executable code.                        
!===================================================================

            c_flag = .false.

            fdel = fd - fu

C !!!! Do something about TOL_SMALL.

!  Y.S. 2/11/05 commented out     rfdel = 1.0 / (fdel + tol_small) 
!   and added:

            rfdel = 1.0 / (fdel + SIGN(tol_small,fdel))
            
	    fnf = (f_adv - fu) * rfdel
            fnc = (fc - fu) * rfdel
            fnmax = fnc * rcfl
            fnc_in_range = sign (0.5, abs(fdel) - abs(fcurvn) )
            fnf_in_range = max (fnc, min(fnf, fnmax, 1.0) )

! We must have S >= 1:
c            courant_in_range = sign(0.5, rcfl - 1.0)
c            if(courant_in_range .eq. -0.5) c_flag = .true.

            limiter = fc*(0.5-fnc_in_range) +   
     &                  (fnf_in_range*fdel+fu)*(0.5+fnc_in_range)

      return
      end

!===================================================================
!===================================================================

      subroutine ADVECT(                                          
     & scheme,                                                       
     & j,                                                             
     & i_max,j_max,k_max,                                             
     & source,                                                        
     & x_div,y_div,h_div,z_div,                                       
     & field,field_b,                                                 
     & field_m,field_bm,                                              
     & field_p,field_bp,                                              
     & field_pp,field_bpp,                                            
     & fuw,fuw_p,fvn,fvn_p,fvs,fvs_p,w,w_p,         
     & fluxs,fluxn,                                          
     & dt,dxtr,dytr,dxur,dyur,dz2r,dzz2r,recip_cos,     
     & levm,lev,levp,levpp,                                           
     & l_oimpaddf,                                                    
     & l_ofreesfc,                                                    
     & l_bootstrap,                                                   
     & l_ocyclic,                                                     
     & j_offset,imout,jmout,imout_hud,jmout_hud,attend,hudtend,       
     & nmedlev,m,nt,l_omedadv,l_ohudout,sf_dtmed,sf_dsmed,wdtmed, 
! the following used to be two seperate lines  
     & wdsmed,cosu,l_xadv,l_yadv,l_zadv                                                 
     & )                                                              

      implicit none

!------------------------------------------------------------------- 
!  Declare argument list.                                            
!------------------------------------------------------------------- 
      integer                                                         
     & j,                      ! in: Local (to this PE) value of j.   
     & j_offset,               ! local offset from global row         
     & i_max,j_max,k_max,      ! in: Local (to this PE) array         
     &                         !     dimensions.                      
     & scheme(2)               ! in: scheme(1) is the choice of       
                               !     scheme in the horizontal and     
                               !     scheme(2) in the vertical.       
                               !     =0 for upstream differencing.    
                               !     =1 for centred differencing      
                               !     =2 for QUICKEST.                    
      INTEGER                                                         
     &    imout(4),jmout(4),imout_hud(4),jmout_hud(4),NMEDLEV,NT,m    
                                                                      
       INTEGER                                                        
     & levm(i_max),   ! in: Number of levels at (i,j-1).                
     & lev(i_max),    ! in: Number of levels at (i,j).                
     & levp(i_max),   ! in: Number of levels at (i,j+1).                
     & levpp(i_max)   ! in: Number of levels at (i,j+2).                
                                                                     
      real                                                           
     & source(i_max,k_max),    ! out: Advection term.                
                                                                     
     & x_div(i_max,k_max),     ! out: STASH diagnostics: x- and y-   
     & y_div(i_max,k_max),     !      divergences, horizontal        
     & h_div(i_max,k_max+1),   !      divergence and z-divergence.   
     & z_div(i_max,k_max),     !                                     
                                                                     
     & field(i_max,k_max),     !                                     
     & field_b(i_max,k_max),   !                                     
     & field_m(i_max,k_max),   !                                     
     & field_bm(i_max,k_max),  ! in: The advected field. The b,m,    
     & field_p(i_max,k_max),   !     and p suffixes have the usual   
     & field_bp(i_max,k_max),  !     Cox code meanings. 
     & field_pp(i_max,k_max),  !                                     
     & field_bpp(i_max,k_max), !                                     
                                                                     
     & fuw(i_max+1,k_max),     !                                     
     & fuw_p(i_max+1,k_max),   !                                     
     & fvn(i_max,k_max),       ! in: Cell face velocities.           
     & fvn_p(i_max,k_max),     ! in: Cell face velocities.           
     & fvs(i_max,k_max),       !                                     
     & fvs_p(i_max,k_max),     !                                     
     & w(i_max,k_max+1),       !                                     
     & w_p(i_max,k_max+1),     !                                     
                                                                     
     & fluxs(i_max,k_max),     ! in:  South-face fluxes.              
     & fluxn(i_max,k_max),     ! out: North-face fluxes.             
     & flux_bot(i_max,k_max+1),! out: Bottom-face fluxes, for STASH  
     &                         !      diagnostics for biology.       
                                                                     
     & dt,                        ! timestep
     & dxtr(i_max),dytr(j_max),   !                                     
     & dxur(i_max),dyur(j_max),   ! in: Reciprocal grid spacings.
     & dz2r(k_max),dzz2r(k_max+1),! 
                                                                     
     & recip_cos(j_max)        ! in: Recip. cosine scaling factors.   
                                                                     
      real                                          
     &     attend(k_max,nt,4)                       
     &,    hudtend(k_max,nt,4)                      
     &,    wdtmed(i_max,k_max)                      
     &,    wdsmed(i_max,k_max)                      
     &,    cosu(j_max)                              
                                                   
      logical                                                        
     & l_oimpaddf,             ! in: Control for Crank-Nicholson.    
     & l_ofreesfc,             ! in: Control for free surface.       
     & l_bootstrap,            ! in: Control for row bootstrap       
c                              !     calculation.                    
     & l_ocyclic               ! in: cyclic e-w condition              
     &, l_omedadv              ! advective med outflow                  
     &, l_ohudout              ! advective hudson bay outflow           
     &, sf_dtmed,sf_dsmed      ! stash flags for med & hud diagnostics  
      logical
     & l_xadv,l_yadv,l_zadv    ! advection diagnostics

!------------------------------------------------------------------    
!  Declare local variables.                                            
!------------------------------------------------------------------    

      integer i,k
      integer ku,icu,kcu
      integer levu,levu1,levu2

      integer signu,signv,signw   
      integer upos,uneg,vpos,vneg  
      integer wpos,wneg            
      integer ufpos, ufneg
      integer vfpos, vfneg
      integer wfpos, wfneg

      real R6   ! reciprocal of 6
      real fmdi ! upper bound on missing value

      real abscn, signcn, cnpos, cnneg
      real fc, fu, fd, fcu
      real uface, vface, wface

      real dzr(k_max),dzzr(k_max+1)
      real dxr0, dyr0, dzr0

      real fcu1, fcu2
      real fgt1, fgt2

      real cflcoef, fcurvn

      real f_adv(i_max,k_max+1)! face value for advection
      real flux(i_max,k_max+1) ! advected flux through face - used
                               ! successively for west and top faces

      real rcfl(i_max,k_max)   ! recip. of sum of outflow courant 
                               ! numbers (="S").
      real rcflp(i_max,k_max)  ! as above for row J+1.
                               ! 
      real rcfl_in             ! value of rcfl passed to limiter.
      real sum_cfl_out         ! 1/S 
     &,fluxmed(i_max,k_max,nt)                
     &,fluxhud(i_max,k_max,nt)                                             
     &,dtwork(i_max,k_max)                                           
     &,dswork(i_max,k_max)                                           

      real temp

      logical c_flag

      real limiter

      external limiter

*CALL COCTOL
*CALL OTIMER

!===================================================================
! Begin executable code.                        
!===================================================================

      IF (L_OTIMER) CALL TIMER('ADVECT_UT',103)
      
!
! Initialise the advection diagnostic arrays.
!
       do k=1,k_max
        do i=1,i_max
          x_div(i,k) = 0.0
          y_div(i,k) = 0.0
          z_div(i,k) = 0.0
          h_div(i,k) = 0.0
        enddo
       enddo


*/ ys, 4/11/05
*/ added the initialisation of source(i,k) 
*/       
       do k=1,k_max
        do i=1,i_max
          source(i,k) = 0.0
        enddo
       enddo


*/ys, 4/11/05
*/ LSETV is a function for initialising variables
*/ that worked on the t3e using benchlib
*/ but is not available on the NEC
*/      if (l_xadv.or.l_yadv.or.l_zadv) then
*/       CALL LSETV(x_div,1,i_max*k_max,0.0)
*/       CALL LSETV(y_div,1,i_max*k_max,0.0)
*/       CALL LSETV(z_div,1,i_max*k_max,0.0)
*/      endif
*/      CALL LSETV(h_div,1,i_max*k_max,0.0)
*/      CALL LSETV(source,1,i_max*k_max,0.0)
*/

      do k=1,k_max                                                    
        dzr(k)  = 2.0*dz2r(k)
        dzzr(k) = 2.0*dzz2r(k)
      enddo 

      dzzr(k_max+1) = 2.0*dzz2r(k_max+1)

      if(scheme(1).eq.3 .and. scheme(2).eq.3) then
!
! Initialise the RCFL array for the flux limiter.
!
        do k=1,k_max                                                    
          do i=1,i_max                                                  

c !!!!!!!!!!!!! NB. Should I be using dxUr etc. here rather than dxTr ??

            sum_cfl_out = dt * (
     &    ( max(0.0,fuw(i+1,k) + abs( min(0.0,fuw(i,k))))) * dxtr(i)
     &  + ( max(0.0,fvn(i,k)   + abs( min(0.0,fvs(i,k))))) * dytr(j)
     &  + ( max(0.0,w(i,k)     + abs( min(0.0,w(i,k+1))))) * dzr(k))

            rcfl(i,k) = 1.0 / (sum_cfl_out + tol_small)

            sum_cfl_out = dt * (
     &    (max(0.0,fuw_p(i+1,k)+abs( min(0.0,fuw_p(i,k)))))*dxtr(i)
     &  + (max(0.0,fvn_p(i,k)  +abs( min(0.0,fvs_p(i,k)))))
     &     *dytr(min(j+1,j_max))
     &  + (max(0.0,w_p(i,k)    +abs( min(0.0,w_p(i,k+1)))))*dzr(k))

            rcflp(i,k) = 1.0 / (sum_cfl_out + tol_small)

          enddo ! over i                                                
        enddo   ! over k                                                

        end if ! scheme.eq.3

      if(.NOT.l_bootstrap) then                                         
C initialise the local outflow variables to zero        
       do k=1,k_max                                     
         do i=1,i_max                                   
           fluxmed(i,k,m)=0.                            
           fluxhud(i,k,m)=0.                            
         enddo                                         
       enddo                                           
!-------------------------------------------------------------------    
!  1st, compute flux through west face of T box                        
!-------------------------------------------------------------------    

      if(scheme(1).eq.0) then                                           
!                                                                        
!  first order upstream differencing                        
!                                                                        
        do k=1,k_max                                                    
          do i=2,i_max                                                  
            signu=sign(0.5,fuw(i,k))                                    
            upos=0.5+signu                                              
            uneg=0.5-signu                                              
            f_adv(i,k)=field_b(i-1,k)*upos+field_b(i,k)*uneg        
          enddo ! over i                                                
          f_adv(1,k)=0.0
        enddo   ! over k                                                
                                                                        
      else if (scheme(1).eq.1) then                  
!                                                                        
!  centred differencing
!                                                                        
        do k=1,k_max                                                    
          do i=2,i_max                                                  
            f_adv(i,k) = 0.5*(field(i-1,k)+field(i,k))             
          enddo ! over i                                                
          f_adv(1,k)=0.0
        enddo   ! over k                                                
                                                                        
      else if(scheme(1).ge.2) then                                         
!                                                                        
!  UTOPIA
!

       R6 = 1./6. 

        do k=1,k_max                                                  
          do i=2,i_max

! Calculate stencil points and 3 components of the velocity at the cell face.

! the normal velocity is in the x direction at west face
            abscn = abs(fuw(i,k)) * dt * dxur(i-1) * recip_cos(j)
            cnpos = nint(0.5 + sign(0.5,fuw(i,k)))             
            cnneg = 1-cnpos
            fd=field_b(i  ,k)*cnpos+field_b(i-1,k)*cnneg
            fc=field_b(i-1,k)*cnpos+field_b(i  ,k)*cnneg       
!  Because there is only one EW wrap-round column in the model, the    
!  upstream point when i=2 and i=i_max is treated as a special case.          
            if(i.eq.2) then
              if(l_ocyclic) then 
!	       
! Yvonne, 26/01/04 changed i_max-1 to i-max-2 in the following 2 lines
! since a bug means the UTOPIA scheme doesn't conserve tracer
! as wrap point 
! See /u/um1/vn5.4/mods/source/osp0504/osp0f504.mf77
!            
                fu=field_b(i_max-2,k)*cnpos+field_b(i+1,k)*cnneg        
                levu=lev(i_max-2)*cnpos+lev(i+1)*cnneg
              else
                fu=field_b(i-1,k)*cnpos+field_b(i+1,k)*cnneg        
                levu=lev(i-1)*cnpos+lev(i+1)*cnneg
              endif
            else if(i.eq.i_max) then
              if(l_ocyclic) then              
                fu=field_b(i-2,k)*cnpos+field_b(3,k)*cnneg        
                levu=lev(i-2)*cnpos + lev(3)*cnneg
              else
                fu=field_b(i-2,k)*cnpos+field_b(i,k)*cnneg        
                levu=lev(i-2)*cnpos+lev(i)*cnneg
              endif
            else
              fu=field_b(i-2,k)*cnpos+field_b(i+1,k)*cnneg        
              levu=lev(i-2)*cnpos+lev(i+1)*cnneg
            endif             
            if (levu .lt. k) fu = fc     

            vface = 0.25 * (   fvn(i-1,k) + fvn(i,k) 
     &                       + fvs(i-1,k) + fvs(i,k)  )
            vfpos = nint( 0.5 + sign (0.5, vface) )
            vfneg = 1-vfpos
            icu = (i-1)*cnpos + i*cnneg 
            fcu1 = field_bm(icu,k)*vfpos + field_bp(icu,k)*vfneg
            levu1 = levm(icu)*vfpos + levp(icu)*vfneg
            if (levu1 .lt. k) fcu1 = fc  

!  Note that w is positive for upward motion; wface is set 
!  to -w to be consistent with level indexing
            wface = -0.25 * (  w(i-1,k) + w(i-1,k+1)
     &                       + w(i  ,k) + w(i  ,k+1) )      
            wfpos = nint( 0.5 + sign (0.5, wface) )
            wfneg = 1-wfpos
            kcu = (k+1)*wfneg + max(1,(k-1))*wfpos
            if (kcu.le.k_max) then
              fcu2 = field_b(icu,kcu)
            else
              fcu2 = fc
            end if
            levu2 = lev(icu)
            if (levu2 .lt. kcu) fcu2 = fc     

! Calculate 1D curv term 
            cflcoef = R6 * ( 1.0 - abscn * abscn )
            fcurvn = fd - 2.0*fc + fu

! First gradt term (in y direction)   
! Yvonne included write statements, 19/01/04
!
!         WRITE(6,*) 'dyr0 before ',dyr0
	 
            dyr0 = dyur(j-1)*vfpos + dyur(j)*vfneg

!         WRITE(6,*) 'dyr0 after ',dyr0
        
	    fgt1 = abs(vface) * 0.5 * dyr0 * dt * (fc-fcu1)

! Second gradt term (in z direction) 
            dzr0 = dzzr(k)*wfpos + dzzr(k+1)*wfneg
            fgt2 = abs(wface) * 0.5 * dzr0 * dt * (fc-fcu2)

! Final expression                                      
            f_adv(i,k) = 0.5*(fd+fc) - 0.5*abscn*(fd-fc)
     &                     - cflcoef*fcurvn 
     &                     - fgt1 - fgt2         

! Use flux limiter if required.
      if(scheme(1).eq.3 .and. scheme(2).eq.3) then
              rcfl_in = rcfl(i-1,k) * cnpos 
     &                + rcfl(i,k)   * cnneg
              f_adv(i,k) = limiter(f_adv(i,k),fu,fc,fd,rcfl_in,
     &                             fcurvn,c_flag)
c              if(c_flag) then 
c                write(6,*) 'LIMITER: courant overflow at i=',i,
c     &                      'j=',j,'k=',k
c              endif
            endif ! scheme.eq.3

          enddo ! over i                                              
          f_adv(1,k)=0.0
        enddo   ! over k                       
      endif   ! scheme(1).eq.0                                          

      do k=1,k_max                                                      
        do i=2,i_max                                                    
          flux(i,k)=fuw(i,k)*f_adv(i,k)                           
        enddo ! over i                                                  
        flux(1,k)=0.0
      enddo   ! over k                                                  
                                                                        
C need to calculate the zonal flux divergence for the Med outflow,    
C which looks at non-adjacent points. Also, in order to separate      
C the rate of change advection and Med outflow diagnostics, the       
C appropriate TEMPA values are zeroed.                                
       if (l_omedadv) then                                            
                                                                      
       if (j+j_offset.eq.jmout(1)) then                               
        do k=1,nmedlev                                                
         fluxmed(imout(1)+1,k,m)=fuw(imout(1)+1,k)*attend(k,m,1)      
         flux(imout(1)+1,k)=0.                                        
        enddo                                                         
       endif                                                          
                                                                      
       if (j+j_offset.eq.jmout(2)) then                               
        do k=1,nmedlev                                                
         fluxmed(imout(2)+1,k,m)=fuw(imout(2)+1,k)*attend(k,m,2)      
         flux(imout(2)+1,k)=0.                                        
        enddo                                                         
       endif                                                          
                                                                      
       if (j+j_offset.eq.jmout(3)) then                               
        do k=1,nmedlev                                                
         fluxmed(imout(3),k,m)=fuw(imout(3),k)*attend(k,m,3)          
         flux(imout(3),k)=0.                                          
        enddo                                                         
       endif                                                          
                                                                      
       if (j+j_offset.eq.jmout(4)) then                               
        do k=1,nmedlev                                                
         fluxmed(imout(4),k,m)=fuw(imout(4),k)*attend(k,m,4)          
         flux(imout(4),k)=0.                                          
        enddo                                                         
       endif                                                          
                                                                      
C need to calculate the zonal flux divergence for the Hudson Bay      
C outflow, which looks at non-adjacent points. Also, in order to      
C separate the rate of change advection and Med outflow diagnostics,  
C the appropriate TEMPA values are zeroed.                            
                                                                      
      if (l_ohudout) then                                             
       if (j+j_offset.eq.jmout_hud(1)) then                           
        do k=1,k_max                                                  
         fluxhud(imout_hud(1)+1,k,m)=fuw(imout_hud(1)+1,k)            
     &                   *hudtend(k,m,1)                              
         flux(imout_hud(1)+1,k)=0.                                    
        enddo                                                         
       endif                                                          
                                                                      
       if (j+j_offset.eq.jmout_hud(2)) then                           
        do k=1,k_max                                                  
         fluxhud(imout_hud(2)+1,k,m)=fuw(imout_hud(2)+1,k)            
     &                   *hudtend(k,m,2)                              
         flux(imout_hud(2)+1,k)=0.                                    
        enddo                                                         
       endif                                                          
                                                                      
       if (j+j_offset.eq.jmout_hud(3)) then                           
        do k=1,k_max                                                  
         fluxhud(imout_hud(3),k,m)=fuw(imout_hud(3),k)                
     &                   *hudtend(k,m,3)                              
         flux(imout_hud(3),k)=0.                                      
        enddo                                                         
       endif                                                          
                                                                      
       if (j+j_offset.eq.jmout_hud(4)) then                           
        do k=1,k_max                                                  
         fluxhud(imout_hud(4),k,m)=fuw(imout_hud(4),k)                
     &                   *hudtend(k,m,4)                              
         flux(imout_hud(4),k)=0.                                      
        enddo                                                         
       endif                                                          
       endif  ! l_ohudout                                             
                                                                      
       endif  ! l_omedadv                                             

!-------------------------------------------------------------------    
!  2nd, compute zonal flux divergence                              
!-------------------------------------------------------------------    
!                                                                       
      do k=1,k_max                                                      
        do i=2,i_max-1   
          source(i,k)=(flux(i,k)-flux(i+1,k))*dxtr(i)*recip_cos(j)  
        end do
      end do

!  Need an if test round all this diagnostics only stuff !!!!!   
      if (l_xadv) then
      do k=1,k_max                                                      
        do i=1,i_max-1
c use if diagnostics required are u.grad(T)    
c         x_div(i,k)=(f_adv(i,k)-f_adv(i+1,k))          
c    &                                       *dxtr(i)*recip_cos(j) 
c         x_div(i,k)=0.5*x_div(i,k)*(fuw(i,k)+fuw(i+1,k))          
c use if diagnostics required are div(uT)                               
           x_div(i,k)=(flux(i,k)-flux(i+1,k))*dxtr(i)*recip_cos(j) 
        enddo ! over i                                                  
      enddo ! over k                                                    
      endif !l_xadv

      if (l_omedadv) then                                  
c     mediterranean outflow diagnostics:                   
        if (sf_dtmed.and.(m.eq.1)) then                    
          do k=1,k_max                                     
            do i=1,i_max                                   
              dtwork(i,k)=source(i,k)                      
            enddo                                          
          enddo                                            
        endif                                              
                                                           
        if (sf_dsmed.and.(m.eq.2)) then                    
          do k=1,k_max                                     
           do i=1,i_max                                    
             dswork(i,k)=source(i,k)                       
           enddo                                           
          enddo                                            
        endif                                              
                                                                       
C                                                                      
C add in Mediterranean zonal flux divergence                           
        do k=1,k_max                                                   
          do i=1,i_max-1                                               
           source(i,k)=source(i,k)+(fluxmed(i,k,m)-fluxmed(i+1,k,m))   
     &                                    *dxtr(i)                      
          enddo                                                        
          source(i_max,k)=0.0                                          
        enddo                                                          
        if (l_ohudout) then                                           
c add in hudson bay zonal flux divergence                             
        do k=1,k_max                                                  
          do i=1,i_max-1                                              
           source(i,k)=source(i,k)+(fluxhud(i,k,m)-fluxhud(i+1,k,m))  
     &                                    *dxtr(i)                     
          enddo                                                       
          source(i_max,k)=0.0                                         
        enddo                                                         
        endif ! l_ohudout                                             
                                                                      
C     Mediterranean outflow diagnostics:                              
        if (sf_dtmed.and.(m.eq.1)) then                               
          do k=1,k_max                                                
           do i=1,i_max                                               
             wdtmed(i,k)=source(i,k)-dtwork(i,k)                      
           enddo                                                      
          enddo                                                       
        endif                                                         
                                                                      
        if (sf_dsmed.and.(m.eq.2)) then                               
          do k=1,k_max                                                
           do i=1,i_max                                               
             wdsmed(i,k)=source(i,k)-dswork(i,k)                      
           enddo                                                      
          enddo                                                       
        endif                                                         
                                                                      
      endif  ! l_omedadv                                              

      endif ! l_bootstrap

!-------------------------------------------------------------------    
!  3rd, compute flux through north face of T box                        
!-------------------------------------------------------------------    
!                                                                       
!  Note have to use fluxn, fluxs variables in the meridional            
!  direction.                                                           
!                                                                       
      if(scheme(1).eq.0) then                                           
!                                                                        
!  first order upstream differencing 
!                                                                        
        do k=1,k_max                                                    
          do i=1,i_max                                                  
            signv=sign(0.5,fvn(i,k))                                    
            vpos=0.5+signv                                              
            vneg=0.5-signv                                              
            f_adv(i,k)=field_b(i,k)*vpos+field_bp(i,k)*vneg        
          enddo ! over i                                                
        enddo   ! over k                                                
                                                                        
      else if (scheme(1).eq.1) then                       
!                                                                        
!  centred differencing
!                                                                        
        do k=1,k_max                                                    
          do i=1,i_max                                                  
            f_adv(i,k) = 0.5*(field(i,k)+field_p(i,k))             
          enddo ! over i                                                
        enddo   ! over k                                                
                                                                        
      else if(scheme(1).ge.2) then                                         
!                                                                        
!  UTOPIA
!
        R6 = 1./6. 

        do k=1,k_max                                                  
          do i=2,i_max-1 

! Calculate stencil points and 3 components of the velocity at the cell face.

! the normal velocity is in the y direction at north face
            abscn = abs(fvn(i,k)) * dt * dyur(j)      
            cnpos = nint(0.5 + sign(0.5,fvn(i,k)))   
            cnneg = 1-cnpos
            fd=field_bp(i,k)*cnpos+field_b(  i,k)*cnneg
            fc=field_b( i,k)*cnpos+field_bp( i,k)*cnneg       

! Yvonne put in conditional WRITE statements, 21/01/04
!   
!           IF (field_bpp(i,k).NE.0) THEN   
!             WRITE(6,*) 'field_bpp(i,k)  ',field_bpp(i,k)
!           ENDIF

            fu=field_bm(i,k)*cnpos+field_bpp(i,k)*cnneg 
	           
            levu = levm(i)*cnpos + levpp(i)*cnneg
            if (levu .lt. k) fu = fc      

!  Note that w is positive for upward motion; wface is set 
!  to -w to be consistent with level indexing
            wface = -0.25 * (  w(i,k)   + w(i,k+1)
     &                       + w_p(i,k) + w_p(i,k+1) )      
            wfpos = nint( 0.5 + sign (0.5, wface) )
            wfneg = 1-wfpos
            kcu = (k+1)*wfneg + max(1,(k-1))*wfpos
            if (kcu.le.k_max) then
              fcu1 = field_b(i,kcu)*cnpos + field_bp(i,kcu)*cnneg
            else
              fcu1 = fc
            end if
            levu1 = lev(i)*cnpos + levp(i)*cnneg
            if (levu1 .lt. kcu) fcu1 = fc     

! Yvonne put in conditional WRITE statements, 19/01/04
!   
!           IF (uface.NE.0) THEN   
!            WRITE(6,*) 'uface before  ',uface
!            ENDIF
            uface = 0.25 * (   fuw(i,k)   + fuw(i+1,k) 
     &                       + fuw_p(i,k) + fuw_p(i+1,k)  )
!   
!            IF (uface.NE.0) THEN   
!             WRITE(6,*) 'uface after  ',uface
!           ENDIF
	 
            ufpos = nint( 0.5 + sign (0.5, uface) )
            ufneg = 1-ufpos
            icu = (i+1)*ufneg + (i-1)*ufpos
            fcu2 = field_b(icu,k)*cnpos + field_bp(icu,k)*cnneg
            levu2 = lev(icu)*cnpos + levp(icu)*cnneg
            if (levu2 .lt. k) fcu2 = fc  

! Calculate 1D curv term 
            cflcoef = R6 * ( 1.0 - abscn * abscn )
            fcurvn = fd - 2.0*fc + fu

! First gradt term (in z direction)   
            if (kcu.le.k_max) then
!              fgt1 = abs(wface) * 0.25 
!     &                 * (dzr(k)+dzr(kcu)) * dt * (fc-fcu1)
              fgt1 = abs(wface) * 0.25  ! FHL fix from Vn 6.1
     &            * (dzr(k)+dzr(MIN(k_max,kcu))) * dt * (fc-fcu1)

            else
              fgt1 = 0
            end if
! Second gradt term (in x direction) 
            fgt2 = abs(uface) * 0.25 
     &                 * (dxtr(i)+dxtr(icu)) * dt * (fc-fcu2)

! Final expression                                         
            f_adv(i,k) = 0.5*(fd+fc) - 0.5*abscn*(fd-fc)
     &                     - cflcoef*fcurvn 
     &                     - fgt1 - fgt2         

! Use flux limiter if required.
            if(scheme(1).eq.3 .and. scheme(2).eq.3) then
              rcfl_in = rcfl(i,k)  * cnpos
     &                + rcflp(i,k) * cnneg
              f_adv(i,k) = limiter(f_adv(i,k),fu,fc,fd,rcfl_in,
     &                             fcurvn,c_flag)
c              if(c_flag) then 
c                write(6,*) 'LIMITER: courant overflow at i=',i,
c     &                      'j=',j,'k=',k
c              endif
            endif ! scheme.eq.3

          enddo ! over i                                              
        enddo   ! over k                                      
      endif   ! scheme(1).eq.0                                          
                                                                        
      do k=1,k_max                                                      
        do i=2,i_max-1                                                    
          fluxn(i,k)=fvn(i,k)*f_adv(i,k)                           
        enddo ! over i                                                  
! Make sure the whole fluxn array is properly initialised.
        fluxn(1,k)=0.0
        fluxn(i_max,k)=0.0
      enddo   ! over k                                                  
                                                                        
      if(.NOT.l_bootstrap) then                                         
!-------------------------------------------------------------------    
!  4th, compute meridional flux divergence                              
!-------------------------------------------------------------------    
!                                                                       
      do k=1,k_max                                                      
        do i=2,i_max-1   
          source(i,k)=source(i,k)+(fluxs(i,k)-fluxn(i,k))*dytr(j)  
        end do
      end do
                 
!  Need an if test round all this diagnostics only stuff !!!!!   
      if (l_yadv) then
      do k=1,k_max                                                      
        do i=1,i_max                           
c use if diagnostics required are u.grad(T)    
!!!!!!!!!!!!!!   is this still valid ?? !!!!!!!!!!!!!!                         
c         if(fvs(i,k).ne.0.0) then                                      
c           y_div(i,k)=(fluxs(i,k)/fvs(i,k)-f_adv(i,k))*dytr(j)     
c         else                                                          
c           y_div(i,k)=-f_adv(i,k)*dytr(j)                          
c         endif                                                         
c         y_div(i,k)=0.5*y_div(i,k)*(fvs(i,k)+fvn(i,k)) 
c use if diagnostics required are div(uT)                               
           y_div(i,k)=(fluxs(i,k)-fluxn(i,k))*dytr(j)                    
c          source(i,k)=source(i,k)+y_div(i,k)                           
c          h_div(i,k)=source(i,k)                                        
        enddo ! over i                                                  
      enddo ! over k                                                    
      endif ! l_yadv

!-------------------------------------------------------------------    
!  5th, compute flux through top face of T box                        
!-------------------------------------------------------------------    

      if(scheme(2).eq.0) then                                           
!                                                                        
!  first order upstream differencing
!                                                                        
        do k=2,k_max                                                    
          do i=1,i_max                                                  
            signw=sign(0.5,w(i,k))                                    
            wpos=0.5+signw                                              
            wneg=0.5-signw                                              
            f_adv(i,k)=field_b(i,k)*wpos+field_b(i,k-1)*wneg        
          enddo ! over i                                                
        enddo   ! over k                                                
                                                                        
      else if (scheme(2).eq.1) then 
!                                                                        
!  centred differencing                                   
!                                                                        
        do k=2,k_max                                                    
          do i=1,i_max                                                  
            f_adv(i,k) = 0.5*(field(i,k)+field(i,k-1))             
          enddo ! over i                                                
        enddo   ! over k                                                
                                                                        
      else if(scheme(2).ge.2) then                                         
!                                                                        
!  UTOPIA
!
        R6 = 1./6. 

        do k=2,k_max                                                  
          do i=2,i_max-1

! Calculate stencil points and 3 components of the velocity at the cell face.

! the normal velocity is in the z direction at top face
            abscn = abs(w(i,k)) * dt * dzzr(k)                             
            cnpos = nint(0.5 + sign(0.5,w(i,k))) 
            cnneg = 1-cnpos
            fd=field_b(i,k-1)*cnpos+field_b(i,k  )*cnneg
            fc=field_b(i,k  )*cnpos+field_b(i,k-1)*cnneg       
            if(k.eq.2) then 
              ku=(k+1)*cnpos + (k-1)*cnneg
            else if(k.eq.k_max) then
              ku=(k)*cnpos + (k-2)*cnneg
            else 
              ku=(k+1)*cnpos + (k-2)*cnneg
            endif
            fu=field_b(i,ku)
            if(lev(i).lt.ku) fu = fc

            uface = 0.25 * (  fuw(i  ,k  ) + fuw(i  ,k-1)
     &                      + fuw(i+1,k  ) + fuw(i+1,k-1) )      
            ufpos = nint( 0.5 + sign (0.5, uface) )
            ufneg = 1-ufpos
            icu = (i-1)*ufpos + (i+1)*ufneg
            kcu = (k  )*cnpos + max(1,(k-1))*cnneg
            fcu1 = field_b(icu,kcu)
            levu1 = lev(icu)
            if(levu1.lt.kcu) fcu1 = fc

            vface = 0.25 * (   fvn(i,k) + fvn(i,k-1) 
     &                       + fvs(i,k) + fvs(i,k-1)  )
            vfpos = nint( 0.5 + sign (0.5, vface) )
            vfneg = 1-vfpos
            fcu2 = field_bm(i,kcu)*vfpos + field_bp(i,kcu)*vfneg
            levu2 = levm(i)*vfpos + levp(i)*vfneg
            if(levu2.lt.kcu) fcu2 = fc

! Calculate 1D curv term 
            cflcoef = R6 * ( 1.0 - abscn * abscn )
            fcurvn = fd - 2.0*fc + fu

! First gradt term (in x direction)   
! Yvonne included write statements, 15/01/04
!
!        WRITE(6,*) 'dxr0 before ',dxr0
	 
            dxr0 = dxur(i-1)*ufpos + dxur(i)*ufneg
	    
!        WRITE(6,*) 'dxr0 after  ',dxr0
	
            fgt1 = abs(uface) * 0.5 * dxr0 * dt * (fc-fcu1)

! Second gradt term (in y direction)   
            dyr0 = dyur(j-1)*vfpos + dyur(j)*vfneg
            fgt2 = abs(vface) * 0.5 * dyr0 * dt * (fc-fcu2)

! Final expression                                              
            f_adv(i,k) = 0.5*(fd+fc) - 0.5*abscn*(fd-fc)
     &                     - cflcoef*fcurvn 
     &                     - fgt1 - fgt2         

! Use flux limiter if required.
            if(scheme(1).eq.3 .and. scheme(2).eq.3) then
              rcfl_in = rcfl(i,k-1) * cnpos
     &                + rcfl(i,k)   * cnneg
              f_adv(i,k) = limiter(f_adv(i,k),fu,fc,fd,rcfl_in,
     &                             fcurvn,c_flag)
c              if(c_flag) then 
c                write(6,*) 'LIMITER: courant overflow at i=',i,
c     &                      'j=',j,'k=',k
c              endif
            endif ! scheme.eq.3

          enddo ! over i                                              
        enddo   ! over k                    
      endif   ! scheme(1).eq.0                                          
                                                                        
      do k=2,k_max                                                      
        do i=1,i_max                                                    
          flux(i,k)=w(i,k)*f_adv(i,k)                           
        enddo ! over i                                                  
      enddo   ! over k                                                  
!                                                                      
! The following calculation for the flux at the surface for the free   
! surface solution follows the method used by Killworth and used in    
! the MOMA code. It is not second order accurate.                      
!                                                                   
      if (l_ofreesfc) then                                             
        do i=1,i_max                                                   
          f_adv(i,1)=field(i,1)                                   
          f_adv(i,k_max+1)=0.0                                    
          flux(i,1)=w(i,1)*field(i,1)                                  
          flux(i,k_max+1)=0.0                                          
        enddo                                                          
      else                                                             
        do i=1,i_max                                                   
          f_adv(i,1)=0.0                                          
          f_adv(i,k_max+1)=0.0                                    
          flux(i,1)=0.0                                                
          flux(i,k_max+1)=0.0                                          
        enddo                                                          
      endif     ! l_ofreesfc                                           
                                                                       
!-------------------------------------------------------------------    
!  6th, compute vertical flux divergence                              
!-------------------------------------------------------------------    
!                                                                       
      do k=1,k_max                                                     
        do i=2,i_max-1   
          source(i,k)=source(i,k)+(flux(i,k+1)-flux(i,k))*dzr(k)  
        end do
      end do
                 
!  Need an if test round all this diagnostics only stuff !!!!!   
      if (l_zadv) then
      do k=1,k_max                                                      
        do i=1,i_max                           
c use if diagnostics required are u.grad(T)    
!!!!!!!!!!!!!!   is this still valid ?? !!!!!!!!!!!!!!                         
c         z_div(i,k)=(f_adv(i,k+1)-f_adv(i,k))*dzr(k)   
c         z_div(i,k)=0.5*z_div(i,k)*(w(i,k)+w(i,k+1))             
c use if diagnostics required are div(uT)                               
           z_div(i,k)=(flux(i,k+1)-flux(i,k))*dzr(k)                    
c          source(i,k)=source(i,k)+y_div(i,k)                           
c          h_div(i,k)=source(i,k)                                        
        enddo ! over i                                                  
      enddo ! over k                                                    
      endif ! l_zadv

      endif ! l_bootstrap

      IF (L_OTIMER) CALL TIMER('ADVECT_UT',104)


      return
      end
*/ <end of mod UTOPIA_basic_bit_ys_4p7.mf77>
