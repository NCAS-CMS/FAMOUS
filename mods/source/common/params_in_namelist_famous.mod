*ID SLAB
*/
*/ This is a FAMOUS mod. 
*/ 
*/ copied from David Sexton's version, but changed to include
*/ Z0MSEA_VAL from ctile_9_new, because mods clashed. Should now be used in
*/ conjunction with ctile_9_new_param.
*/  CDJ 3/9/02
*/  
*/-------------------------------------------------------------------------
*/  
*/ CLIMATEPREDICTION.COM - Slab model mod -
*/ This HadAM3 vn4.5 mod enables to transfer the hardwired
*/ precipitation parameter "VF1" and ..."CT" to a namelist file
*/ which is related to the ...
*/      Pierre-philippe MATHIEU, 12 january 2001
*/      mathieu@met.rdg.ac.uk
*/
*/ Adjusted to read SLBC21 from the CNTLATM namelist with the RUNCNST
*/ variables so that cp.com only needs to change CNTLATM for each
*/ participant.
*/      Dave Stainforth. 22/1/1
*/
*/ More parameters added to SLBC21. SLBC21 is now store in its own
*/ COMDECK C_SLBC21 declared above. The variables are removed from the
*/ original COMDECKs and replaced by commented out lines that set the
*/ variables up as REAL variables rather than PARAMETERS. This makes it 
*/ easy to make a mod which changes the PARAMETERs to variables - 
*/ necessary to obtain bit comparison with this mod, as PARAMETER and
*/ REALs compile differently.
*/
*/      David Sexton. 30/07/02
*/ 
*/ Achieved bit comparison at end of atmosphere timestep. Needed 
*/ -g option on compiler for LSPFOR2B, LSPFOR2D, and LSCLD1A to do this.
*/ 
*/      David Sexton. 13/08/02
*/ 
*/ EACF now read is an array with a value for each level. *CALLs to
*/ CMAXSIZE added where needed.
*/
*/
*/      David Sexton. 13/08/02
*/
*/
*COMDECK C_SLBC21
C CMAXSIZE must be called before *CALLing this comdeck.
C Declare variables for use in climateprediction.com. DS
      REAL VF1,CT,ENTCOEF,CAPE_TS,G0,CHARNOCK,Z0FSEA
      INTEGER R_LAYERS(4)  ! Number of soil layers from which  
                           ! water can be extracted
      REAL ICE_SIZE        ! New variable which can be set in CNTLATM
                           ! for ice particle size			   
      REAL EACF(MAX_P_LEVELS) ! New variable which can be set in CNTLATM
                           ! for EACF scheme			   
      REAL ASYM_LAMBDA     ! New variable which can be set in CNTLATM
                           ! for asymptotic length scale lambda			   
      COMMON /COMMON_SLBC21/ VF1,CT,ENTCOEF,CAPE_TS,G0,CHARNOCK,Z0FSEA,
     +                       ICE_SIZE,EACF,ASYM_LAMBDA,R_LAYERS    
C
*DECLARE READLSA1
*I READLSA1.80
*CALL C_SLBC21
      NAMELIST /SLBC21/ VF1,CT,ENTCOEF,CAPE_TS,G0,CHARNOCK,Z0FSEA,
     +                  ICE_SIZE,EACF,ASYM_LAMBDA,R_LAYERS
      
*I READLSA1.131
C Read in variables for use in climateprediction.com. DS
      WRITE(6,*) 'DS in READSA1: read RUNCNST' !ppm
      REWIND 5            !ppm
      READ(5,SLBC21)      !ppm
      WRITE(6,*) 'DS in READSA1: '
      WRITE(6,10) 'VF1:      ',VF1
      WRITE(6,10) 'CT:       ',CT
      WRITE(6,10) 'ENTCOEF:  ',ENTCOEF
      WRITE(6,13) 'CAPE_TS:  ',CAPE_TS
      WRITE(6,10) 'G0:       ',G0
      WRITE(6,10) 'CHARNOCK: ',CHARNOCK
      WRITE(6,10) 'Z0FSEA:   ',z0FSEA
      WRITE(6,14) 'ICE_SIZE: ',ICE_SIZE
      WRITE(6,10) 'ASYM_LAM: ',ASYM_LAMBDA
      WRITE(6,11) 'R_LAYERS: ',R_LAYERS(1),R_LAYERS(2),
     +             R_LAYERS(3),R_LAYERS(4)
      WRITE(6,*) 'EACF at 1st, 6th, 7th and 8th levels'
      WRITE(6,12) EACF(1),EACF(6),EACF(7), EACF(8)
      WRITE(6,11) 'R_LAYERS: ',R_LAYERS(1),R_LAYERS(2),
     +             R_LAYERS(3),R_LAYERS(4)
10    format(A10,F9.5)
11    format(A10,4I5)
12    format(4F8.4)
13    format(A10,F7.1)
14    format(A10,F10.8)
*/
*/ This is based on PPM's and DS's original code to change VF1 and CT
*/
*DECLARE C_LSPFRM
*D CR190293.6
C      REAL CT,VF1
      REAL CA,PCF,CFMIN,FPL,VFPOWER
*D CR190293.11
C     +CT=1.0E-4,        ! Reciprocal time constant for conversion of 
*D CR190293.18
C     +VF1=1.0,          ! Tunable speed in Heymsfield formula (m/s). 
*D ARN1F304.6
C      REAL CT,VF1
      REAL CA,PCF,CFMIN,FPL,VFPOWER
*D ARN1F304.11
C     +CT=1.0E-4,        ! Reciprocal time constant for conversion of
*D ARN1F304.18
C     +VF1=1.0,          ! Tunable speed in Heymsfield formula (m/s).
*D AYY2F400.28
C      REAL CT,VF1
      REAL CA,PCF,CFMIN,FPL,VFPOWER
*D AYY2F400.33
C     +CT=1.0E-4,        ! Reciprocal time constant for conversion of
*D AYY2F400.35
C     +VF1=1.0,          ! Tunable speed in Heymsfield formula (m/s).
*/
*DECLARE LSPFOR2B
*I LSPFOR2B.78
*CALL CMAXSIZE
*CALL C_SLBC21
*/
*DECLARE LSPFOR2D
*I LSPFOR2D.102
*CALL CMAXSIZE
*CALL C_SLBC21
*/
*/
*/
*/ Added by David Sexton to include ENTCOEF 30/07/02
*/ Adding commented out lines to set REAL ENTCOEF=3.0
*/
*DECLARE ENTCNST
*D API2F400.404
*D API2F400.407
C      REAL ENTCOEF
C      PARAMETER(ENTCOEF=3.0)
      PARAMETER (SH_FAC=1.0)    
*/
*DECLARE CONVEC3A
*I CONVEC3A.68
*CALL CMAXSIZE    
*CALL C_SLBC21
*/
*DECLARE CONVEC3C
*I CONVEC3C.96    
*CALL CMAXSIZE    
*CALL C_SLBC21
*/
*DECLARE LAYCN1A
*I LAYCN1A.39    
*CALL CMAXSIZE    
*CALL C_SLBC21
*/
*DECLARE LAYCN3C
*I LAYCN3C.59    
*CALL CMAXSIZE    
*CALL C_SLBC21
*/
*DECLARE LAYERD2A
*I LAYERD2A.43    
*CALL CMAXSIZE    
*CALL C_SLBC21
*/
*DECLARE LAYERD2C
*I ADR1F405.55    
*CALL CMAXSIZE    
*CALL C_SLBC21
*/
*DECLARE LAYERD3C
*I ADR1F405.64    
*CALL CMAXSIZE    
*CALL C_SLBC21
*/
*/
*/
*/
*/ Added by David Sexton to include CAPE_TS 30/07/02
*/
*DECLARE CAPECNST
*D CAPECNST.2,CAPECNST.6
C      REAL CAPE_TS      !  TIMESCALE FOR DESTRUCTION OF CONVECTIVE
C                        !  AVAILABLE POTENTIAL ENERGY BY CONVECTION
C                       !  WHEN A CAPE CLOSURE TO THE CONVECTION
C                        !  SCHEME IS EMPLOYED (S)        
C      PARAMETER(CAPE_TS = 7200.0                                   )
C
*/
*/ Would have to declare following code if CAPE_TS was being included
*/ in namelist and VF1 and CT were not. But as they are, comment out
*/ lines to avoid COMMON blocks in C_SLBC21 being declared twice in CONVECT.
*/
*/*DECLARE CONVEC3A
*/*I CONVEC3A.68
*/*CALL CMAXSIZE
*/*CALL C_SLBC21
*/*/
*/*DECLARE CONVEC3C
*/*I CONVEC3C.96
*/*CALL CMAXSIZE
*/*CALL C_SLBC21
*/
*/
*/
*/
*/ Added by David Sexton to include G0 31/07/02
*/ As DM and DH in original EXCOEF code are PARAMETERS that depend
*/ on G0 remove these from PARAMETER list and set up elsewhere in
*/ code.
*/
*DECLARE EXCOEF3A
*D EXCOEF3A.100
C      REAL G0
C      PARAMETER(G0=10.0,DH=G0/EH,DM=G0/EM)
      REAL EH,EM,DH,DM,LAMBDA_MIN,A_LAMBDA                             
*D EXCOEF3A.104,EXCOEF3A.106 
C     +,G0=10.0                  ! Used in stability function calcs. 
C     +,DH=G0/EH                 ! Used in calc of stability function FH.  
C     +,DM=G0/EM                 ! Used in calc of stability function FM.
*I EXCOEF3A.99
*CALL CMAXSIZE
*CALL C_SLBC21
*I EXCOEF3A.154
      DH=G0/EH                 ! Used in calc of stability function FH.
      DM=G0/EM                 ! Used in calc of stability function FM.
*/
*/
*DECLARE EXCOEF6A
*D EXCOEF6A.119   
C      REAL G0
C      PARAMETER(G0=10.0,DH=G0/EH,DM=G0/EM)
      REAL EH,EM,DH,DM,LAMBDA_MIN,A_LAMBDA                             
*D EXCOEF6A.123,EXCOEF6A.125 
C     +,G0=10.0                  ! Used in stability function calcs. 
C     +,DH=G0/EH               ! Used in calc of stability function FH.
C     +,DM=G0/EM               ! Used in calc of stability function FM.  
*I EXCOEF6A.117
*CALL CMAXSIZE
*CALL C_SLBC21
*I EXCOEF6A.172
      DH=G0/EH                 ! Used in calc of stability function FH.
      DM=G0/EM                 ! Used in calc of stability function FM.
*/
*/
*DECLARE EXCOEF7A
*D EXCOEF7A.112   
C      REAL G0
C      PARAMETER(G0=10.0,DH=G0/EH,DM=G0/EM)
      REAL EH,EM,DH,DM,LAMBDA_MIN,A_LAMBDA
*D EXCOEF7A.116,EXCOEF7A.118 
C     +,G0=10.0                ! Used in stability function calcs. 
C     +,DH=G0/EH               ! Used in calc of stability function FH.
C     +,DM=G0/EM               ! Used in calc of stability function FM.
*I EXCOEF7A.110
*CALL CMAXSIZE
*CALL C_SLBC21
*I EXCOEF7A.163
      DH=G0/EH                 ! Used in calc of stability function FH.
      DM=G0/EM                 ! Used in calc of stability function FM.
*/ The following code is added to check for bit comparison and can be omitted later.
*/*DECLARE EXCOEF3A
*/*I EXCOEF3A.109   
*/      REAL G0_2,DH_1,DM_1
*/      PARAMETER(G0_2=10.0,DH_1=G0_2/EH,DM_1=G0_2/EM)
*/*I EXCOEF3A.155
*/      DH=DH_1                 ! Used in calc of stability function FH.
*/      DM=DM_1                 ! Used in calc of stability function FM.
*/*/
*/*DECLARE EXCOEF6A
*/*I EXCOEF6A.129   
*/      REAL G0_2,DH_1,DM_1
*/      PARAMETER(G0_2=10.0,DH_1=G0_2/EH,DM_1=G0_2/EM)
*/*I EXCOEF6A.173
*/      DH=DH_1                 ! Used in calc of stability function FH.
*/      DM=DM_1                 ! Used in calc of stability function FM.
*/*/
*/*DECLARE EXCOEF7A
*/*I EXCOEF7A.122   
*/      REAL G0_2,DH_1,DM_1
*/      PARAMETER(G0_2=10.0,DH_1=G0_2/EH,DM_1=G0_2/EM)
*/*I EXCOEF7A.164
*/      DH=DH_1                 ! Used in calc of stability function FH.
*/      DM=DM_1                 ! Used in calc of stability function FM.
*/
*/
*/
*/
*/
*/ Added by David Sexton to include CHARNOCK 31/07/02
*/
*DECLARE C_CHARNK
*D C_CHARNK.5,C_CHARNK.7
C      REAL CHARNOCK
C      PARAMETER(CHARNOCK=0.012)
*D ARN1F404.41,ARN1F404.43
C      REAL CHARNOCK
C      PARAMETER(CHARNOCK=0.011)
*/
*DECLARE SFEXCH3A
*I SFEXCH3A.326
*CALL CMAXSIZE   
*CALL C_SLBC21
*/
*DECLARE SFEXCH5A
*I SFEXCH5A.401   
*CALL CMAXSIZE   
*CALL C_SLBC21
*/
*DECLARE SFEXCH6A
*I SFEXCH6A.402   
*CALL CMAXSIZE   
*CALL C_SLBC21
*/
*DECLARE SFEXCH7A
*I SFEXCH7A.282   
*CALL CMAXSIZE   
*CALL C_SLBC21
*/
*/
*/
*/
*/
*/ Added by David Sexton to include R_LAYERS 31/07/02
*/ Based on Pete COx's original modset, R_LAYERS is not
*/ changed in BDYLYR5A.
*/* DECLARE BDYLYR5A
*/*I BDYLYR5A.705
*/*CALL CMAXSIZE   
*/*CALL C_SLBC21
*/*D BDYLYR5A.707
*/C      INTEGER R_LAYERS(4)
*/C      DATA R_LAYERS/   4,   4,   3,   3/
*/*D BDYLYR5A.712
*/*/
*DECLARE SMROOT5A
*D SMROOT5A.110,SMROOT5A.111
C      INTEGER R_LAYERS(4)
C      DATA R_LAYERS/   4,   4,   3,   3/
*D SMROOT5A.116
*I SMROOT5A.109
*CALL CMAXSIZE   
*CALL C_SLBC21
*/
*/
*/
*/
*/
*/
*/ Added by David Sexton to include Z0FSEA 31/07/02
*/
*/ changed to include Z0MSEA_VAL from ctile_9_new, because mods clashed
*/  CDJ 3/9/02
*/
*DECLARE C_ROUGH
*D C_ROUGH.11,C_ROUGH.14
      REAL Z0MSEA_VAL,Z0HSEA,Z0MIZ,Z0SICE
C      REAL Z0FSEA
C      PARAMETER(Z0FSEA=1.3E-3)
      PARAMETER(Z0HSEA = 1.0E-4,    
     &          Z0MSEA_VAL = 10.0E-4,
*/
*/ Would have to declare following code if CAPE_TS was being included
*/ in namelist and VF1 and CT were not. But as they are, comment out
*/ lines to avoid COMMON blocks in C_SLBC21 being declared twice in CONVECT.
*/
*/*DECLARE SFEXCH3A
*/*I SFEXCH3A.330
*/*CALL CMAXSIZE   
*/*CALL C_SLBC21
*/
*/*DECLARE SFEXCH5A
*/*I SFEXCH5A.405   
*/*CALL CMAXSIZE   
*/*CALL C_SLBC21
*/
*/*DECLARE SFEXCH6A
*/*I SFEXCH6A.405   
*/*CALL CMAXSIZE   
*/*CALL C_SLBC21
*/
*/*DECLARE SFEXCH7A
*/*I SFEXCH7A.285   
*/*CALL CMAXSIZE   
*/*CALL C_SLBC21
*/
*/ The following DECLARE block should be included if using BL scheme 6A.
*/
*/*DECLARE SFLBES5B
*/*I SFLBES5B.110   
*/*CALL CMAXSIZE   
*/*CALL C_SLBC21
*/
*DECLARE SFRUGH5A
*I SFRUGH5A.112   
*CALL CMAXSIZE   
*CALL C_SLBC21
*/
*DECLARE SFRUGH6A
*I SFRUGH6A.114   
*CALL CMAXSIZE   
*CALL C_SLBC21
*/
*/
*/
*/
*/
*/
*/ Modset to change the default size of ice crystals in HadAM3 at
*/ version 4.5
*/               J. M. Edwards 07/01/02
*/
*/ Change the effective radius of ice spheres from 30 microns to
*/ 25 microns
*/
*/ adapted from ice_size 
*/ Gareth S. Jones 17/1/2002
*/
*/ Added by David Sexton  31/07/02
*/
*/
*DECLARE FILL3A
*I FILL3A.787
*CALL CMAXSIZE
*CALL C_SLBC21
*I FILL3A.788
C     ICE_SIZE=30.E-6    ! the default value for ice particle size
*D ADB2F404.89
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)=ICE_SIZE
*D ADB2F404.122
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)=ICE_SIZE
*/
*/
*/
*/
*/
*/ Added by David Sexton to include aysmptotic length scale lambda 31/07/02
*/ This section assumes that C_SLBC21 has already been added to EXCOEF*.dk for G0
*/
*DECLARE EXCOEF3A
*D EXCOEF3A.177
        LAMBDAM = MAX ( LAMBDA_MIN , ASYM_LAMBDA*ZH(I) )   
*I EXCOEF3A.150
C     ASYM_LAMBDA=0.15
*/
*DECLARE EXCOEF6A
*D ARN0F405.803
        LAMBDAM = MAX ( LAMBDA_MIN , ASYM_LAMBDA*ZH_LOCAL(I) )
*I EXCOEF6A.168
C     ASYM_LAMBDA=0.15
*/
*DECLARE EXCOEF7A
*D EXCOEF7A.179
        LAMBDAM = MAX ( LAMBDA_MIN , ASYM_LAMBDA*ZH(I) )
*I EXCOEF7A.160
C     ASYM_LAMBDA=0.15
*/
*/
*/
*/
*/ 
*/ Based on Mod to include the empirically adjusted cloud fraction at vn4.5.
*/ CF=0.6 when (QT/QSAT)=1 for all clouds.
*/ I cannot guarantee at present that this will be compatible with
*/ versions 2D or 2E of section 4.
*/ 
*/                                                  S. Cusack June 2002
*/
*/ Included EACF 31/07/02. David Sexton
*/
*/ To include EACF, we test if EACF=0.5. If so, use old code, if not then
*/ then use adjusted QCN. Could have tried to have one bit of code that 
*/ does it all but tried this and it lost bit comparability.
*/ 
*/      David Sexton. 13/08/02
*/
*/
*DECLARE LSCLD1A
*/
*D AYY2F400.127
     &                  INDEX,QC_POINTS,POINTS,K)
*D AYY2F400.131
     &,INDEX,POINTS,POINTS_F,K)
*I AAD2F404.211
      INTEGER K
*I LSCLD1A.339
     &,QN_ADJ              ! Adjusted value of QN
*/
*/*D AAD2F404.217
*/*D AAD2F404.218
*/*D AAD2F404.219
*/*D AAD2F404.220
*I LSCLD1A.323
C Extra CNTLATM namelist COMDECK
*CALL CMAXSIZE
*CALL C_SLBC21
      REAL RK1, RK2,QN2
*I LSCLD1A.119
C Set up comment for EACF default. DS
C     EACF=0.5
*D LSCLD1A.361,LSCLD1A.386
*I AYY1F401.25
!
! if EACF is effectively 0.5 do old way, else use adjusted QCN
        IF (ABS(EACF(K)-0.5).LT.1.0E-12) THEN
          IF (RHCRIT .LT. 1.) THEN                                      
!-----------------------------------------------------------------------
!L 2. Calculate cloud fraction C, BS (ie. sigma*sqrt(6), where sigma is 
!L    as in P292.14) and normalised cloud water QCN=qc/BS, using eqs    
!L    P292.15 & 16 if RHcrit < 1.                                       
!  N.B. QN (input) is initially in QCN                                   
!  N.B. QN does not depend on AL and so CF and QCN can be calculated     
!       outside the iteration (which is performed in LS_CLD_C).          
!       QN is > -1 for all points processed so CF > 0.                     
!----------------------------------------------------------------------- 
!                                                                        
            BS(I) = (1.0 - RHCRIT) * AL * QSL_F(II)    ! P292.14         
            IF (QCN(I) .LE. 0.) THEN                                     
              CF_F(II)=0.5*(1.+QCN(I))*(1.+QCN(I))                       
              QCN(I)=(1.+QCN(I))*(1.+QCN(I))*(1.+QCN(I))/6.              
            ELSEIF (QCN(I) .LT. 1.) THEN                                  
              CF_F(II)=1.-0.5*(1.-QCN(I))*(1.-QCN(I))                    
              QCN(I)=QCN(I) + (1.-QCN(I))*(1.-QCN(I))*(1.-QCN(I))/6.    
            ELSE ! QN .GE. 1                                            
              CF_F(II)=1.                                                
            ENDIF ! Tests on QN                                         
          ELSE ! i.e. if RHcrit =1                                       
C----------------------------------------------------------------------- 
!L 3.a Set the cloud fraction to 1 if RHcrit = 1.                        
C      For the case RHcrit =1, QN is > 0 for all points processed         
C      so CF =1.                                                           
C----------------------------------------------------------------------- 
            BS(I) = AL                                                     
            CF_F(II) = 1.                                                  
          ENDIF ! Test on RHCRIT                                          
!	  
	ELSE
!	
	  IF (EACF(K).GT.0.5) THEN
	    QN2=1.0 - SQRT( 2.0 * (1.0 - EACF(K)) )
	    RK2=QN2 / ( 1.0 + QN2 )
            RK1=1.0 - RK2
            QN_ADJ=(QCN(I)+RK2)/RK1
	  ELSE
	    QN2=SQRT( 2.0 * EACF(K) ) -1.0
	    RK2=QN2 /( 1.0 + QN2 )
            RK1=1.0 - RK2
            QN_ADJ=(QCN(I)+RK2)/RK1
	  ENDIF
!
	  IF (RHCRIT .LT. 1.) THEN
            BS(I) = (1.0 - RHCRIT) * AL * QSL_F(II)    ! P292.14         
            IF (QN_ADJ .LE. 0.) THEN
              CF_F(II) = 0.5 * (1.0 + QN_ADJ) * (1.0 + QN_ADJ)
            ELSEIF (QN_ADJ .LT. 1.) THEN
              CF_F(II) = 1. - 0.5 * (1.0 - QN_ADJ) * (1.0 - QN_ADJ)
            ELSE ! QN .GE. 1
              CF_F(II) = 1.
            END IF ! QCN_if 
          ELSE ! i.e. if RHcrit = 1
            BS(I) = AL                                                     
            CF_F(II) = 1.
          END IF ! Rhcrit_if1
	  
          IF (QCN(I) .LE. 0.) THEN                                     
            QCN(I)=(1.+QCN(I))*(1.+QCN(I))*(1.+QCN(I))/6.              
          ELSEIF (QCN(I) .LT. 1.) THEN                                  
            QCN(I)=QCN(I) + (1.-QCN(I))*(1.-QCN(I))*(1.-QCN(I))/6.    
          ENDIF ! Tests on QN                                         
	ENDIF
*/
*/ The following lines are included in the mod lux32bit.mod.jul3
*/ which should already be included in any 32 bit compilation.
*/ D. Stainforth 12/4/01
*/
*/*DECLARE CONVEC3C
*/*D CONVEC3C.764
*/        PARAMETER (SAFETY_MARGIN = 1.0E-44 )      
*/
*/*DECLARE GRB1F405
*/*D GRB1F405.83
*/        PARAMETER (SAFETY_MARGIN = 1.0E-44 )      
