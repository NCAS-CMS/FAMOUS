*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  FHL_bugs.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*DECLARE BL_CTL1
*D ACN1F405.42
            IF (ABS(D1(J_CO2_EMITS+I-1)) <=  ABS(0.999*RMDI)) THEN
*D ACN1F405.48
            IF (ABS(D1(J_CO2FLUX+I-1)) <=  ABS(0.999*RMDI)) THEN
*D ACN1F405.53
            IF (ABS(LAND_CO2(I)) <=  ABS(0.999*RMDI)) THEN
*DECLARE TRANO2A1
*D TRANO2A1.378 
        IF (ABS(UOUT(I,J)) <=  ABS(0.999*RMDI)) THEN
*D CCN1F405.337
            IF (ABS(CO2FLUXOUT(i,j))  <=  ABS(0.999*RMDI)) THEN
*DECLARE MICROB7A
*B MICROB7A.55
     &,STH_NEW                   ! WORK defines new regime for
!                                ! soil respiration dependance
!                                ! on moisture
*I MICROB7A.68
      STH_NEW = 0.0             ! FHL: NEC QUMP-like behaviour
*D MICROB7A.78,86
          IF (STH_SOIL(L) .LE. STH_WILT) THEN
            FSTH = 0.2
          ELSEIF (STH_SOIL(L) .GT. STH_WILT .AND.
     &            STH_SOIL(L) .LE. STH_NEW) THEN
            FSTH = 0.2 + 0.8 * ((STH_SOIL(L) - STH_WILT)
     &                        / (STH_NEW - STH_WILT))
          ELSEIF (STH_SOIL(L) .GT. STH_NEW .AND.
     &            STH_SOIL(L) .LE. STH_OPT) THEN
            FSTH = 1.0
          ELSEIF (STH_SOIL(L) .GT. STH_OPT) THEN
            FSTH = 1 - 0.8 * (STH_SOIL(L) - STH_OPT)
          ENDIF
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  LANDTOGLOBAL
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/-------------------------------------------------------------------
*/-------------------------------------------------------------------
*/-------------------------------------------------------------------
*/ This is a FAMOUS mod. 
*/
*/ Add new subroutine:
*DECK LAND2GLO
CLL Subroutine LAND_TO_GLOBAL ---------------------------------------
CLL
CLL  Purpose: convert real land-only field to global field
CLL
      SUBROUTINE LAND_TO_GLOBAL(LAND,FIELD_L,FIELD_G,LAND_PTS,P_POINTS)
C
      IMPLICIT NONE
C
      INTEGER
     & LAND_PTS                    !IN Number of land points to be
C                                  !   processed
     &,P_POINTS                    !IN Number of global points to be
C                                  !   processed
C
      REAL
     & FIELD_L(P_POINTS)           !IN land-only field
     & ,FIELD_G(P_POINTS)          !OUT Global field
C
      LOGICAL
     &  LAND(P_POINTS)             !IN Land mask
C
      INTEGER
     & I,L                         ! LOCAL Loop counters
C
      L=0
      DO I=1,P_POINTS
        IF(LAND(I))THEN
          L=L+1
          FIELD_G(I)=FIELD_L(L)
        ELSE
          FIELD_G(I)=0.0
        ENDIF
      ENDDO
C
      RETURN
      END
*/
*/---------------------------------------------------------------------
*/
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  PTOUV_LAND
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*DECK PTULAND
*/ 
*/ This is a FAMOUS mod. 
*/ 
*IF DEF,C90_1A,OR,DEF,C90_2A,OR,DEF,C90_2B
C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1996, METEOROLOGICAL OFFICE, All Rights Reserved.
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
CLL  SUBROUTINE P_TO_UV---------------------------------------------
CLL
CLL  Purpose:   Interpolates a horizontal land field from pressure to
CLL             wind points on an Arakawa B grid. Under UPDATE
CLL             identifier GLOBAL the data is assumed periodic along
CLL             rows. Otherwise, the last value on each row is set 
CLL             equal to the penultimate value on each row.Output array
CLL             contains one less row than the input array.
CLL
CLL  Not suitable for single column use.
CLL
CLL  Adapted from PTOUV by Nic Gedney (09/99)
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL  Programming standard:  Unified Model Documentation Paper No 3
CLL
CLL
CLLEND-------------------------------------------------------------
 
C
C*L  Arguments:----------------------------------------------------
      SUBROUTINE P_TO_UV_LAND
     1  (P_DATA,U_DATA,P_FIELD,U_FIELD,FLANDG,ROW_LENGTH,ROWS)
 
      IMPLICIT NONE
 
      INTEGER
     &  ROWS               !IN    Number of rows to be updated.
     &, ROW_LENGTH         !IN    Number of points per row
     &, P_FIELD            !IN    Number of points in input field
     &, U_FIELD            !IN    Number of points in output field
 
      REAL
     & P_DATA(P_FIELD)     !INOUT Data on p points
     &,U_DATA(U_FIELD)     !OUT   Data on uv points
     &,FLANDG(P_FIELD)      !IN    Fraction of land in grid box.
C*---------------------------------------------------------------------
 
C*L  Local arrays:-----------------------------------------------------
C    None
C*---------------------------------------------------------------------
*IF DEF,MPP
! Parameters and Common blocks
*CALL PARVARS
*ENDIF
 
C*L  External subroutine calls:----------------------------------------
C    None
C*---------------------------------------------------------------------
 
C----------------------------------------------------------------------
C    Define local variables
C----------------------------------------------------------------------
      INTEGER
     &  U_POINTS      !     Number of values at u points
     &,I              !     Horizontal loop indices
 
      REAL
     &  TOTFRAC       !     Fractional weighting
C---------------------------------------------------------------------
CL    1.     Initialise local constants
C---------------------------------------------------------------------
 
      U_POINTS      =  ROW_LENGTH * (ROWS-1)
 
C---------------------------------------------------------------------
CL    2.     Calculate horizontal average at u points
C---------------------------------------------------------------------
 
      DO 200 I=1,U_POINTS-1
       TOTFRAC=FLANDG(I)+FLANDG(I+1)+FLANDG(I+ROW_LENGTH) +
     &         FLANDG(I+1+ROW_LENGTH)
       IF(TOTFRAC.GT.0.0)THEN
         U_DATA(I)=FLANDG(I)*P_DATA(I)+FLANDG(I+1)*P_DATA(I+1) +
     &           FLANDG(I+ROW_LENGTH)*P_DATA(I+ROW_LENGTH) +
     &           FLANDG(I+1+ROW_LENGTH)*P_DATA(I+1+ROW_LENGTH)
         U_DATA(I)=U_DATA(I)/TOTFRAC
       ELSE
         U_DATA(I)=0.0
       ENDIF
200   CONTINUE
 
C  End points
 
*IF DEF,GLOBAL
*IF -DEF,MPP
C Cyclic wrap around
      DO 201 I=ROW_LENGTH,U_POINTS,ROW_LENGTH
       TOTFRAC=FLANDG(I)+FLANDG(I+1-ROW_LENGTH)+FLANDG(I+ROW_LENGTH) +
     &           FLANDG(I+1)
       IF(TOTFRAC.GT.0.0)THEN
         U_DATA(I)=FLANDG(I)*P_DATA(I) +
     &           FLANDG(I+1-ROW_LENGTH)*P_DATA(I+1-ROW_LENGTH) +
     &           FLANDG(I+ROW_LENGTH)*P_DATA(I+ROW_LENGTH) +
     &           FLANDG(I+1)*P_DATA(I+1)
         U_DATA(I)=U_DATA(I)/TOTFRAC
       ELSE
         U_DATA(I)=0.0
       ENDIF
201   CONTINUE
*ELSE
! Cyclic wrap around already taken account of via halo
*ENDIF
*ELSE
C Set last values on row
*IF -DEF,MPP
      DO 201 I=ROW_LENGTH,U_POINTS,ROW_LENGTH
       U_DATA(I)=U_DATA(I-1)
201   CONTINUE
*ELSE
       IF (atright) THEN
! Set the last values on row
         DO I=ROW_LENGTH,U_POINTS,ROW_LENGTH
           U_DATA(I-Offx)=U_DATA(I-Offx-1)
         ENDDO
       ENDIF
 
*ENDIF
 
*ENDIF
*IF DEF,MPP
! Set a sensible number in the bottom right corner
      U_DATA(U_POINTS)=U_DATA(U_POINTS-1)
 
*ENDIF
 
      RETURN
      END
*ENDIF
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  PTOUV_SEA
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*DECK PTUSEA
*/
*/ This is a FAMOUS mod. 
*/
*IF DEF,C90_1A,OR,DEF,C90_2A,OR,DEF,C90_2B
C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1996, METEOROLOGICAL OFFICE, All Rights Reserved.
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
CLL  SUBROUTINE P_TO_UV---------------------------------------------
CLL
CLL  Purpose:   Interpolates a horizontal sea field from pressure to
CLL             wind points on an Arakawa B grid. Under UPDATE
CLL             identifier GLOBAL the data is assumed periodic along
CLL             rows. Otherwise, the last value on each row is set 
CLL             equal to the penultimate value on each row.Output array
CLL             contains one less row than the input array.
CLL
CLL  Not suitable for single column use.
CLL
CLL  Adapted from PTOUV by Nic Gedney (09/99)
 
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL  Programming standard:  Unified Model Documentation Paper No 3
CLL
CLL
CLLEND-------------------------------------------------------------
 
C
C*L  Arguments:----------------------------------------------------
      SUBROUTINE P_TO_UV_SEA
     1  (P_DATA,U_DATA,P_FIELD,U_FIELD,FLANDG,ROW_LENGTH,ROWS)
 
      IMPLICIT NONE
 
      INTEGER
     &  ROWS               !IN    Number of rows to be updated.
     &, ROW_LENGTH         !IN    Number of points per row
     &, P_FIELD            !IN    Number of points in input field
     &, U_FIELD            !IN    Number of points in output field
 
      REAL
     & P_DATA(P_FIELD)     !INOUT Data on p points
     &,U_DATA(U_FIELD)     !OUT   Data on uv points
     &,FLANDG(P_FIELD)      !IN    Fraction of land in grid box.
C*---------------------------------------------------------------------
 
C*L  Local arrays:-----------------------------------------------------
C    None
C*---------------------------------------------------------------------
*IF DEF,MPP
! Parameters and Common blocks
*CALL PARVARS
*ENDIF
 
C*L  External subroutine calls:----------------------------------------
C    None
C*---------------------------------------------------------------------
 
C----------------------------------------------------------------------
C    Define local variables
C----------------------------------------------------------------------
      INTEGER
     &  U_POINTS      !     Number of values at u points
     &,I              !     Horizontal loop indices
 
      REAL
     &  TOTFRAC       !     Fractional weighting
C---------------------------------------------------------------------
CL    1.     Initialise local constants
C---------------------------------------------------------------------
 
      U_POINTS      =  ROW_LENGTH * (ROWS-1)
 
C---------------------------------------------------------------------
CL    2.     Calculate horizontal average at u points
C---------------------------------------------------------------------
 
      DO 200 I=1,U_POINTS-1
       TOTFRAC=4.-(FLANDG(I)+FLANDG(I+1)+FLANDG(I+ROW_LENGTH) +
     &         FLANDG(I+1+ROW_LENGTH))
       IF(TOTFRAC.GT.0.0)THEN
         U_DATA(I)=
     &           (1.-FLANDG(I))*P_DATA(I)+(1.-FLANDG(I+1))*P_DATA(I+1)+
     &           (1.-FLANDG(I+ROW_LENGTH))*P_DATA(I+ROW_LENGTH) +
     &           (1.-FLANDG(I+1+ROW_LENGTH))*P_DATA(I+1+ROW_LENGTH)
         U_DATA(I)=U_DATA(I)/TOTFRAC
       ELSE
         U_DATA(I)=0.0
       ENDIF
200   CONTINUE
 
C  End points
 
*IF DEF,GLOBAL
*IF -DEF,MPP
C Cyclic wrap around
      DO 201 I=ROW_LENGTH,U_POINTS,ROW_LENGTH
       TOTFRAC=4.-(FLANDG(I)+FLANDG(I+1-ROW_LENGTH) +
     &           FLANDG(I+ROW_LENGTH) + FLANDG(I+1))
       IF(TOTFRAC.GT.0.0)THEN
         U_DATA(I)=(1.-FLANDG(I))*P_DATA(I) +
     &           (1.-FLANDG(I+1-ROW_LENGTH))*P_DATA(I+1-ROW_LENGTH) +
     &           (1.-FLANDG(I+ROW_LENGTH))*P_DATA(I+ROW_LENGTH) +
     &          (1.- FLANDG(I+1))*P_DATA(I+1)
         U_DATA(I)=U_DATA(I)/TOTFRAC
       ELSE
         U_DATA(I)=0.0
       ENDIF
201   CONTINUE
*ELSE
! Cyclic wrap around already taken account of via halo
*ENDIF
*ELSE
C Set last values on row
*IF -DEF,MPP
      DO 201 I=ROW_LENGTH,U_POINTS,ROW_LENGTH
       U_DATA(I)=U_DATA(I-1)
201   CONTINUE
*ELSE
       IF (atright) THEN
! Set the last values on row
         DO I=ROW_LENGTH,U_POINTS,ROW_LENGTH
           U_DATA(I-Offx)=U_DATA(I-Offx-1)
         ENDDO
       ENDIF
 
*ENDIF
 
*ENDIF
*IF DEF,MPP
! Set a sensible number in the bottom right corner
      U_DATA(U_POINTS)=U_DATA(U_POINTS-1)
 
*ENDIF
 
      RETURN
      END
*ENDIF
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  UTOPIA_4p5_divUTw.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  a2o_specifiedCO2.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*IDENT CO2FIX
*/
*/ Ensure ocean biogechemistry sees the same OC2 levels as 
*/ the atmosphere for runs with a fully coupled carbon cycle. 
*/ 
*/ 
*DECLARE SWAPA2O2
*B SWAPA2O2.186
*CALL CRUNTIMC
*CALL UMSCALAR
*B SWAPA2O2.608
c**********************************************************************
c****PASS THE SPATIALLY CONSTANT ATMOSPHERIC pCO2 INTO THE OCEAN*******
c****(values copied across from TRANA2O)*******************************
      IF (.NOT. L_CO2_INTERACTIVE) THEN
        PCO2_ATM_0=CO2_MMR*1e6/1.5194
      END IF
c**********************************************************************

*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  abortfix.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID ABORTFIX
*DECLARE UMSHELL1
*/
*/  Do not call the MPP timer if there has been an error.
*/  This is because the call can hang in situations where only one PE
*/  has returned to UMSHELL because of an error code, while all the other
*/  PEs are executing a different communications call.
*/
*/  Also, flush the output buffer before doing a potential abort;
*/  much improves the usefulness of the leave file.
*/
*/  Alan Iwi -- 4/2/02
*/
*D GSM1F401.25
*IF DEF,MPP,AND,DEF,C97_3A
      if (iCode .eq. 0) then
        CALL TIMER('UM_SHELL',2)
      else
        print *,'Model aborted; timings not available'
      endif
*ELSE
      CALL TIMER('UM_SHELL',2)
*ENDIF      
*/
*B UMSHELL1.101
      call flush(6)
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  alk_2.mod_4.7.mh
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID ALK
*/
*/  set rain_ratio to 0.0070 from 0.01 on advice from Ian Tots.
*/  this will result in spun-up alkalinity nearer to GEOSECS
*/
*/  used for runs post abnd... which had psmax=0.6
*/
*DECLARE OBIOCONST
*D OBIOCONST.87
      PARAMETER (rain_ratio = 0.0070) ! carbon export as calcite, as a
*/
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  ask1f406
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID ASK1F406
*/ U.M. 4.6 unix / source code change form / header   version 04/01/99
*/Instructions: see http://fr0800/umdoc/hegui/t3e4.6.html#chgfinst
*/ 
*/SOC:-> prevent bit-negative arguments to SQRT in new deck at 4.5
*/ 
*/->: fuller description of mod.: purpose, relevant configurations, 
*/-> applicable previous releases, dependencies
*/   At vn4.5 a new parametrization of RHcrit was built into the UM.
*/   The subroutine which calculates the value of RHcrit uses the SQRT
*/   function four times. The argument to the SQRT function should never
*/   be less than zero, but at the bit level, on very rare occasions,the
*/   argument may become negative and cause a model crash. This mod
*/   ensures the argument to the SQRT calls are always greater than or
*/   equal to zero.
*/   THIS MOD IS ONLY REQUIRED FOR 4.5 JOBS WHICH USE THE RHcrit 
*/   PARAMETRIZATION
*/ 
*/ Has an entry been lodged in the Problem Reporting System? Y - number 410
*/
*/ THIS CODE IS INTENDED FOR INCLUSION IN THE 4.6 BUILD      Y
*/ .....................................................................
*/   Author[s]:-> E-mail:-> scusack@meto.gov.uk 
*/ Reviewer[s]:-> E-mail:-> abushell@meto.gov.uk
*/
*/    "I have checked this change. When provided, the advance design 
*/  specification was agreed and adequate, and the new code conforms to
*/  Unified Model standards."
*/
*/  DESIGN SPEC. WAS REVIEWED ON: ......      REVIEWER[S] SIGNATURES
*/                                            ----------------------
*/    DATE CODE REVIEWED: ......
*/  .....................................................................
*/
*/  WILL CHANGES AFFECT ANCILLARY FILES?         N
*/  ARE ANY CHANGES TO STASHMASTER FILES NEEDED? N 
*/  USER INTERFACE ACTION REQUIRED?              N
*/ 
*/  TESTED IN CONFIGURATIONS:-> CLIMATE
*/  TESTS RUN BY [PERSON]:-> S Cusack
*/ 
*/  WILL THE CHANGES SLOW DOWN THE MODEL?        N
*/  -> Further details
*/  CHANGES WILL INCREASE MEMORY CONSUMPTION?    N   
*/  -> Further details
*/
*/ | Re-start dumps bit compare with those created without the change 
*/ V MARK [Y| ] BELOW; leave rest of lines untouched.
*/    
*/   Control Code    loses bit comparison
*/   Atmosphere (assuming same science options chosen)   loses b.c.
*/   Ocean       loses bit comparison
*/   Wave        loses bit comparison
*/   Reconfiguration   loses bit comparison
*/   Diagnostics      lose bit comparison
*/
*/                    ***** Bit-Comparable (unless failure is occurring)
*/
*/ For Y2K compliance checking:  
*/ DOES THIS CHANGE INTERACT WITH DATE CALCULATIONS IN ANY WAY?  N 
*/ 
*/  SECTIONS (TO BE) CHANGED: Section 9, version 2B
*/
*/  SECTIONS (TO BE) DELETED? None
*/
*/  NEW SECTIONS?  Fill in form http://www-hc/~hadmk/STASHmaster_change.html,
*/  and give section numbers below:
*/  
*/  *DEFS ADDED OR REMOVED: None
*/
*/  **Existing** decks being changed [with *I, *D, *B directives]
*/ -> CALRHC2B
*/
*/  Decks being created or purged [with *DECK, *COMDECK, *PURGEDK]
*/ *......K  Deck name   Section#.vr
*/ -> None
*/ ......................................................................
*/ ANY REFERENCES TO EXTERNAL DOCUMENTS-> instead of design spec.
*/  ...OR ... ADVANCE DESIGN SPECIFICATION (optional) 
*/ ->    NONE
*/
*//--------------
*DECLARE CALRHC2B
*//--------------
*I CALRHC2B.59
!  4.6    06/01/99   Ensure arguments to SQRT function never negative.
!                                                    S. Cusack
*I CALRHC2B.266
! 
          TOT_VAR=ABS(TOT_VAR)
*I CALRHC2B.327
            TOT_VAR=ABS(TOT_VAR)
*I CALRHC2B.386
! 
            TOT_VAR=ABS(TOT_VAR)
*I CALRHC2B.443
              TOT_VAR=ABS(TOT_VAR)
*/
*//////////////////////////////////////////////////////////////////////// 
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  ask6f406
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID ASK6F406
*/ U.M. 4.6 unix / source code change form / header   version 06/01/99
*/ CODE WRITERS MUST READ THE ACCOMPANYING INSTRUCTIONS FOR THIS BUILD:
*/  - See http://fr0800/umdoc/hegui/t3e4.6.html#chgfinst
*/ 
*/SOC:-> Fix a bug in tracer advection.
*/->: A mod introduced at vn4.5 to optimise tracer advection used an unset
*/->  array in TRAVAD1A.dk - -this is now set.
*/
*/ Has an entry been lodged in the Problem Reporting System? [Y|N]      
*/
*/ THIS CODE IS INTENDED FOR INCLUSION IN THE 4.6 BUILD      [Y|N]
*/ .....................................................................
*/   Author[s]:-> E-mail:-> scusack@meto.gov.uk 
*/ Reviewer[s]:-> E-mail:-> pmburton@meto.gov.uk
*/
*/    "I have checked this change. When provided, the advance design 
*/  specification was agreed and adequate, and the new code conforms to
*/  Unified Model standards."
*/
*/  DESIGN SPEC. WAS REVIEWED ON: ......      REVIEWER[S] SIGNATURES
*/                                            ----------------------
*/    DATE CODE REVIEWED: ......
*/  .....................................................................
*/
*/  WILL CHANGES AFFECT ANCILLARY FILES?         [N]
*/  ARE ANY CHANGES TO STASHMASTER FILES NEEDED? [N] 
*/  USER INTERFACE ACTION REQUIRED?              [N]
*/ 
*/  TESTED IN CONFIGURATIONS:-> mesoscale
*/  TESTS RUN BY [PERSON]:-> S. Cusack
*/ 
*/  WILL THE CHANGES SLOW DOWN THE MODEL?        [N]
*/  -> Further details
*/  CHANGES WILL INCREASE MEMORY CONSUMPTION?    [N]   
*/  -> Further details
*/
*/ | Re-start dumps bit compare with those created without the change 
*/ V MARK [Y| ] BELOW; leave rest of lines untouched.
*/
*/   Control Code    loses bit comparison
*/   Atmosphere (assuming same science options chosen)   loses b.c.
*/   Ocean       loses bit comparison
*/   Wave        loses bit comparison
*/   Reconfiguration   loses bit comparison
*/   Diagnostics      lose bit comparison
*/ For Y2K compliance checking:  
*/ DOES THIS CHANGE INTERACT WITH DATE CALCULATIONS IN ANY WAY? [N]   
*/ 
*/  SECTIONS (TO BE) CHANGED: 
*/
*/  SECTIONS (TO BE) DELETED? 
*/
*/  NEW SECTIONS?  Fill in form http://www-hc/~hadmk/STASHmaster_change.html,
*/  and give section numbers below:
*/  
*/  *DEFS ADDED OR REMOVED: 
*/
*/  **Existing** decks being changed [with *I, *D, *B directives]
*/ ->
*/
*/  Decks being created or purged [with *DECK, *COMDECK, *PURGEDK]
*/ *......K  Deck name   Section#.vr
*/ -> 
*/ ......................................................................
*/ ANY REFERENCES TO EXTERNAL DOCUMENTS-> instead of design spec.
*/  ...OR ... ADVANCE DESIGN SPECIFICATION (optional) 
*/ ->    
*//////////////////////////////////////////////////////////////////////// 
*DECLARE ATMDYN1
*I GPB7F405.18
         DO I=1,P_FIELD
           WORK2(I)=0.0
         END DO
*I GPB7F405.37
         DO I=1,P_FIELD
           WORK2(I)=0.0
         END DO
*I GPB7F405.53
         DO I=1,P_FIELD
           WORK2(I)=0.0
         END DO
*I GPB7F405.70
         DO I=1,P_FIELD
           WORK2(I)=0.0
         END DO
*I AWO2F405.29
!
         DO I=1,P_FIELD
           WORK2(I)=0.0
         END DO
*I GPB7F405.87
         DO I=1,P_FIELD
           WORK2(I)=0.0
         END DO
*I GPB7F405.103
         DO I=1,P_FIELD
           WORK2(I)=0.0
         END DO
*I GPB7F405.120
         DO I=1,P_FIELD
           WORK2(I)=0.0
         END DO
*I GPB7F405.135
         DO I=1,P_FIELD
           WORK2(I)=0.0
         END DO
*I AWO2F405.105
         DO I=1,P_FIELD
           WORK2(I)=0.0
         END DO
*I AWO2F405.180
         DO I=1,P_FIELD
           WORK2(I)=0.0
         END DO
*I AWO2F405.255
         DO I=1,P_FIELD
           WORK2(I)=0.0
         END DO
*I ACN2F405.168
!
         DO I=1,P_FIELD
           WORK2(I)=0.0
         END DO
*/ END OF MOD
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  asm1f406
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*IDENT ASM1F406
*/ U.M. 4.6 unix / source code change form / header   version 08/12/98
*/Instructions: see http://fr0800/umdoc/hegui/t3e4.6.html#chgfinst
*/ 
*/ SOC: Correct overwriting of U if 4,218 selected but L_MURK false
*/
*/ Prevents serious overwriting of U prognostic if 4,218 selected
*/ but total aerosol fields not selected. Advised for use in ALL
*/ atmosphere runs.
*/ 18/03/99
*/
*/ Has an entry been lodged in the Problem Reporting System? [Y]
*/
*/ THIS CODE IS INTENDED FOR INCLUSION IN THE 4.6 BUILD      [Y]
*/ .....................................................................
*/   Author[s]:-> Steve Mullerworth E-mail:-> sdmullerworth@meto.gov.uk 
*/ Reviewer[s]:-> E-mail:-> paclark@meto.gov.uk
*/
*/    "I have checked this change. When provided, the advance design 
*/  specification was agreed and adequate, and the new code conforms to
*/  Unified Model standards."
*/
*/  DESIGN SPEC. WAS REVIEWED ON: ......      REVIEWER[S] SIGNATURES
*/                                            ----------------------
*/    DATE CODE REVIEWED: ......
*/  .....................................................................
*/ ANY REFERENCES TO EXTERNAL DOCUMENTS-> instead of design spec.
*/
*/  WILL CHANGES AFFECT ANCILLARY FILES?         [N]
*/  ARE ANY CHANGES TO STASHMASTER FILES NEEDED? [N] 
*/  USER INTERFACE ACTION REQUIRED?              [N]
*/ 
*/  TESTED IN CONFIGURATIONS:-> HadCM3
*/  TESTS RUN BY [PERSON]:-> Steve Mullerworth
*/ 
*/  WILL THE CHANGES SLOW DOWN THE MODEL?        [N]
*/  -> Further details
*/  CHANGES WILL INCREASE MEMORY CONSUMPTION?    [N]   
*/  -> Further details
*/
*/ | Forecast dumps bit compare with those created without the change 
*/ V MARK [Y| ] BELOW; leave rest of lines untouched.
*/
*/   Control Code    loses bit comparison
*/   Atmosphere (assuming same science options chosen)   loses b.c.
*/   Ocean       loses bit comparison
*/   Wave        loses bit comparison
*/   Reconfiguration   loses bit comparison
*/   Diagnostics      lose bit comparison
*/ For Y2K compliance checking:  
*/ DOES THIS CHANGE INTERACT WITH DATE CALCULATIONS IN ANY WAY? [N]   
*/ 
*/  SECTIONS (TO BE) CHANGED:
*/
*/  SECTIONS (TO BE) DELETED? 
*/
*/  NEW SECTIONS?  Fill in form http://www-hc/~hadmk/STASHmaster_change.html,
*/  and give section numbers below:
*/  
*/  *DEFS ADDED OR REMOVED: 
*/ ......................................................................
*/                ADVANCE DESIGN SPECIFICATION (optional) 
*/ ->    
*//////////////////////////////////////////////////////////////////////// 
*DECLARE VISBTY1A
*I APC0F405.119   
!LL  4.6  18/03/99  Correct overwriting of U field that occurs if 4,218
!LL                 selected but L_MURK false. S.D.Mullerworth
*I APC0F405.234   
      IF (L_MURK) THEN
*D APC0F405.236,VISBTY1A.87   
        DO Point=1,P_FIELD                                                
          AEROSOL(Point)=MAX(AEROSOL(Point),AERO0)                        
        ENDDO
      ENDIF
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  atmstep_flush.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID FLUSH
*/
*/  Flush the buffers on unit 6 after printing the "ATMOS TIMESTEP" line.
*/  This means that the output file of a running job shows more accurately
*/  which timestep the job has reached.
*/
*DECLARE ATMSTEP1
*I ARB0F400.45
      Call flush(6)
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  boundsfix_famous.mod_xcpsa
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*IDENT BOUNDSFAMOUS
*/
*/ Fixes necessary for bounds checking FAMOUS. 
*/
*/ These fixes are all from Paul Valdes's coupled_bugs4.
*/ 	 Annette Osprey 19-March-2008
*/ 
*/ bounds_fix_famous2.mod :
*/ Fixed so that pp_lblev set correctly for ocean.
*/       Annette Osprey 06-January-2009
*/
*/
*/ This next one arose from famous and is to do with an invalid
*/ bound on an if statement array index. (PJV)
*/
*DECLARE PPHEAD1A
*I PPHEAD1A.126
     *,  ilevelfix    ! used for array bound problem in IF statement
*D GPB0F403.29,34
           if (lvcode.eq.ppx_half_level) then
              ilevelfix=LevIndex
              if (ilevelfix.gt.t_levels+1) ilevelfix=t_levels+1
              if (bkh(ilevelfix).eq.1.0) then
                 pp_lblev=ppx_meto8_surf
              else
                 pp_lblev=level+0.00001
              end if
           else
              pp_lblev=level+0.00001
           end if
*/
*/
*/   Another bounds problem, found for famous.
*/   This has been solved by adhoc fix but could do
*/   with locating real problem. (PJV)
*/
*DECLARE STZONM1A
*I STZONM1A.81
      INTEGER JLOOP
*D GPB0F403.2853
      JLOOP=n_rows_full_zonal_data
      if (jloop.gt.global_yend-global_ystart+1) 
     :    jloop=global_yend-global_ystart+1
      do j=1,jloop
*/
*/
*/   The CO2 related code calls a subroutine with one argument
*/   D(JA_CO2) (or equivalent flux) but JA_CO2 is 0
*/   if CO2 flag if false. This causes a bounds error, although
*/   is not important since the array is never used within the
*/   subroutine. Solution is a but cludgy but works.
*/
*DECLARE SWAPO2A1
*D CCN1F405.412
     & O_CO2FLUX,D1(MAX(JA_CO2FLUX,1)),
     : CO2_ICOLS,CO2_JROWS,CO2_IMT,CO2_JMT,
*DECLARE SWAPA2O1
*D CCN1F405.372   
     & ATMCO2, D1(MAX(JO_CO2,1)), 
     : CO2_ICOLS, CO2_JROWS, CO2_IMT, CO2_JMT,
*/
*/
*/  Bounds problem when calculating flow between med and Atlantic
*/  Standard version of low resolution model only uses one point
*/  on each side of straits (standard HadCM3 uses average of two
*/  points). The second pair of points are set to coordinate zero
*/  and this causes a bounds error (but code appears to be 
*/  otherwise OK). Add if statement, so that code is only executed 
*/  if using all 4 points.
*/
*/  Also, last line of this routine needs to be deleted to make
*/  that update can be used in all model configurations (otherwise
*/  when using atmosphere only, subroutine ends up being 1 byte long)
*/
*DECLARE READ_REM
*I READ_REM.165
      if (jmsend(3).gt.0) then
*I READ_REM.219
      end if
*I READ_REM.222
      if (jmsend(4).gt.0) then
*I READ_REM.275
      end if
*D READ_REM.283
*/
*/
*/   Another bounds problem, found using famous when adding
*/   standard BRIDGE STASH output. Like the bug above, the
*/   solution is a pragmatic fix without understanding cause
*/   of problem
*/
*DECLARE INITCTL1
*D GSS1F305.459
           if (-stlist(st_input_bottom,ii).le.NUM_LEVEL_LISTS) then
              NUM_LEVELS=STASH_LEVELS(1,-STLIST(st_input_bottom,II))
           else 
              num_levels=1 
           end if
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  boundsfix_hadocc.mod_xdbua
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*IDENT BOUNDSHADOCC
*/
*/ Bounds checking fix for HadOCC
*/
*DECLARE BIOLOGY
*B OJP0F404.233
            if(kmt(i).gt.0) then
*I OJP0F404.237
            end if
*B OJP0F404.250 
        if(kmt(i).gt.0) then
*B OJP0F404.281
        end if !kmt gt 0
*DECLARE ROW_CTL
*I ORH3F405.14
        IF(SI_BIO_LOCAL(ITEM).lt.1) SI_BIO_LOCAL(ITEM)=1
            
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  boundsfix_nonmpp.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*IDENT BOUNDSNONMPP
*/
*/ Fixes necessary for bounds checking non-MPP code. 
*/
*/ These fixes are all from Paul Valdes's coupled_bugs4.
*/ 	Annette Osprey 11-Sep-2006
*/
*/
*/ First problem was that someone forgot to caveat last call 
*/ for when not running in MPP. (PJV)
*/
*DECLARE CNVSTOP
*B CNVSTOP.173
*IF DEF,MPP
*I CNVSTOP.174
*ENDIF
*/
*/
*/   This is a bounds problem assoicated with the calculation of 
*/   mountain torque. It only effects the non-MPP version. Fairly
*/   confident about the solution. (PJV)
*/
*DECLARE DYNDIA1A
*D GSM1F405.732
        DO i=FIRST_FLD_PT,LAST_U_FLD_PT-1
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  boundsfix_vn4.5.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*IDENT BOUNDSFIX
*/
*/ General fixes necessary for bounds checking vn4.5 code.
*/ Tested on HadCM3 and FAMOUS
*/ 
*/ These fixes are all from Paul Valdes's coupled_bugs4.
*/ 	Annette Osprey 11-Sep-2006
*/
*/
*/   For some unknown reason, with some combinations of stash/input 
*/   dumps, the reconfiguration gives a variable which has undefined
*/   (i.e. 0) stashcode. This causes a bounds error. Not too sure
*/   about solution since the new if statements means that no space
*/   is allocated for these mysterious stash variables which may be
*/   OK but ...... (PJV)
*/
*DECLARE INANCA1A
*I GRB4F305.223
        if (stashancil(i).gt.0) then
*I GRB4F305.232
        end if
*I INANCA1A.655
        if (stashancil(i).gt.0) then
*I GDR8F400.95
        end if
*I ADR1F304.115
        if (fileancil(i).gt.0) then
*I INANCA1A.701
        end if
*/
*/
*/   This is a real bounds error problem but don't fully
*/   understand. In PP files, level dependent constants are
*/   defined with a second dimension 4 within subroutine
*/   initpp but the calling argument suggests that they 
*/   should be dimensioned 1 for ocean pp files. Simple
*/   correction stops bounds error. (PJV)
*/
*DECLARE INITPP1A
*D INITPP1A.98
      PP_FIXHD(112)=
     : MIN(LEN2_LEVDEPC,PP_LEN2_LEVDEPC)
*D INITPP1A.136
      DO 5 II=1,LEN1_LEVDEPC*
     :    MIN(LEN2_LEVDEPC,PP_LEN2_LEVDEPC)
*/
*/
*/   The FFT routine FOURIER passes array arguments to FTRANS in the 
*/   f77 style by the starting array element only. FTRANS accesses 
*/   beyond the dimensions of the dummy array (which causes bounds 
*/   error), but refers to element of the actual array in FOURIER 
*/   (presumably at the appropriate place). 
*/   So just define all as assumed-size arrays in FTRANS. (AO)
*/
*DECLARE FOURIE3A
*D PXORDER.17,21
      REAL A(*),     ! First real input vector
     &     B(*),
     &     C(*),     ! First real output vector
     &     D(*),
     &     TRIGS(*)  ! Precalculated list of sines & cosines
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  diagsw_inlansea.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/
*/ if L_INLANSEA is used to cap SSS values to be between
*/ SALLOW and SALUP, the diagnostic DIAGSW (2,30,280) is 
*/ modified to reflect the virtual salinity flux that would
*/ be required to make this change. DIAGSW has the same sense
*/ as P-E however, such that a negative change in SSS (fresher)
*/ requires a positive change in DIAGSW (more water). This
*/ appears to have been originally coded with the opposite sign
*/
*/ r.s.smith@reading.ac.uk 24/09/15
*/ 
*DECLARE TRACER
*D OJL1F405.31
     &          -1e9*dz(1)*(sallow-TA(I,1,2))/(C2DTTS*100.0)
*D OJL1F405.36
     &          -1e9*dz(1)*(salup-TA(I,1,2))/(C2DTTS*100.0)
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  dummy.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  fixfill3a.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID FIXFILL3A
*/
*/  fix a multiple initialisation in deck FILL3A
*/
*/  Having initialised the whole array to .FALSE. we are not
*/  allowed to initialise individual elements to .TRUE. using
*/  DATA statements.
*/
*/  Most compilers treat this as a warning, but the Intel compiler
*/  treats it as an error.
*/
*/  Just set the array elements using assignments instead...
*/
*/  Alan Iwi 26/11/01
*/
*DECLARE FILL3A
*D ADB2F404.235,239
*I ADB2F404.247
      L_IN_CLIMAT(IP_WATER_SOLUBLE) = .TRUE.
      L_IN_CLIMAT(IP_DUST_LIKE) = .TRUE.
      L_IN_CLIMAT(IP_OCEANIC) = .TRUE.
      L_IN_CLIMAT(IP_SOOT) = .TRUE.
      L_IN_CLIMAT(IP_SULPHURIC) = .TRUE.
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  fixmeanctl.mf77
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID FIXMEANCTL
*/
*/  There are some lines in mean_ctl which potentially do a divide-by zero.
*/
*/  Some of the time, they should't actually get executed, but still cause
*/  a problem because of a bug in Portland compiler optimisations.
*/
*/  This mod changes comparisons of form:
*/
*/                 i_hour .eq. (24/dumps_per_day)
*/
*/  to an equivalent form which relies on integer multiplication rather than
*/  integer division
*/
*DECLARE MEANCTL1
*D GMG1F404.152
     &     i_hour*dumps_per_day <= 24
     &     .and.(i_hour+1)*dumps_per_day > 24)) then
*D GMG1F404.202
     &     i_hour*dumps_per_day <= 24
     &     .and.(i_hour+1)*dumps_per_day > 24)) then
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  fixsolang.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID FIXSOLANG
*/
*/ Fix a problem which was leading to anomalous clear-sky SW fluxes
*/ at certain points on the equator. (vn4.5) 
*/     -- Alan Iwi, 5/12/02
*/      
*DECLARE SOLANG1A
*I SOLANG1A.70
      REAL EPS
      PARAMETER (EPS=1E-3) ! small number for floating-point comparisons      
*D SOLANG1A.134,137
          IF (DIFTIM.LT.0.) THEN
            IF (ABS(DIFTIM-(DTRAD-TWOPI)).LT.EPS) THEN
              !
              ! This is the genuine case where the sun sets and rises
              ! within the timestep.
              !
              DIFSIN = DIFSIN + 2. * SQRT(1.-COSHLD**2)
              DIFTIM = DIFTIM + 2. * HLD
            ELSEIF (DIFTIM.GT.-EPS) THEN
              !
              ! This is an artefact which can arise where the timestep
              ! ends exactly at sunrise (OMEGAE.eq.0).  This should give 
              ! OMEGA1.eq.-HLD, OMEGA2.eq.-HLD, hence DIFTIM.eq.0, but
              ! under compiler optimisations can sometimes give a tiny
              ! negative value.  (Seen with Intel compiler on linux.)
              ! Set the marker to indicate no sunshine during timestep.
              !
              DIFTIM=0.
            ELSE
              !
              ! If DIFTIM is negative for any other reason, 
              ! something's awry.
              !
              CALL ABORT
     $          ('Unexpected value for DIFTIM in SOLANG; please report')
            ENDIF
          ENDIF
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  fixspin3a.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID FIXSPIN3A
*/
*/  fix a multiple initialisation in deck SPIN3A
*/
*/  Most compilers treat this as a warning, but the Intel compiler
*/  treats it as an error.
*/
*/  Just set the array elements using assignments instead...
*/
*/  Lois Steenman-Clark 18.08.05
*/
*DECLARE SPIN3A
*D SPIN3A.1129
*I SPIN3A.1130
      N_SCALE_VARIABLE(IP_SCALE_POWER_LAW)=2
      N_SCALE_VARIABLE(IP_SCALE_FNC_NULL)=0
      N_SCALE_VARIABLE(IP_SCALE_POWER_QUAD)=3
      N_SCALE_VARIABLE(IP_SCALE_DOPPLER_QUAD)=4
*/
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  fixstdia.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID FIXSTDIA
*/
*/  Change the way that deck ST_DIA11 does the selection of levels
*/  for diagnostics.  This method will choose the *closest* of the
*/  available levels.
*/
*/  The code being replaced insists on an *exact* match; this involves
*/  a floating point equality test -- and if it evaluates false, the
*/  level index is uninitialised and is later used as an array offset,
*/  segfaulting the program.
*/
*DECLARE ST_DIA11      
*D ST_DIA11.338,ST_DIA11.345
      Call level_indices(STASH_LEVELS(1,NI),
     $     T_P_LEVS,T_PRESS,T2_IND)
*D ST_DIA11.366,ST_DIA11.373
      Call level_indices(STASH_LEVELS(1,NI),
     $     UCOMP_P_LEVS,UCOMP_PRESS,U2_IND)
*D ST_DIA11.394,ST_DIA11.401
      Call level_indices(STASH_LEVELS(1,NI),
     $     VCOMP_P_LEVS,VCOMP_PRESS,V2_IND)
*D ST_DIA11.226,ST_DIA11.238
      Call level_indices2(STASH_LEVELS(1,NI),
     $     UCOMP_P_LEVS,UCOMP_PRESS,VCOMP_P_LEVS,VCOMP_PRESS,
     $     UV_IND)
*D ST_DIA11.272,ST_DIA11.284
      Call level_indices2(STASH_LEVELS(1,NI),
     $     UCOMP_P_LEVS,UCOMP_PRESS,T_P_LEVS,T_PRESS,
     $     UT_IND)
*D ST_DIA11.305,ST_DIA11.317
      Call level_indices2(STASH_LEVELS(1,NI),
     $     VCOMP_P_LEVS,VCOMP_PRESS,T_P_LEVS,T_PRESS,
     $     VT_IND)
*D ST_DIA11.435,ST_DIA11.447
      Call level_indices2(STASH_LEVELS(1,NI),
     $     W_P_LEVS,W_PRESS,T_P_LEVS,T_PRESS,
     $     WT_IND)
*D ST_DIA11.468,ST_DIA11.480
      Call level_indices2(STASH_LEVELS(1,NI),
     $     W_P_LEVS,W_PRESS,UCOMP_P_LEVS,UCOMP_PRESS,
     $     WU_IND)
*D ST_DIA11.501,ST_DIA11.513
      Call level_indices2(STASH_LEVELS(1,NI),
     $     W_P_LEVS,W_PRESS,VCOMP_P_LEVS,VCOMP_PRESS,
     $     WV_IND)
*D ST_DIA11.548,ST_DIA11.560
      Call level_indices2(STASH_LEVELS(1,NI),
     $     Q_P_LEVS,Q_PRESS,UCOMP_P_LEVS,UCOMP_PRESS,
     $     QU_IND)
*D ST_DIA11.581,ST_DIA11.593
      Call level_indices2(STASH_LEVELS(1,NI),
     $     Q_P_LEVS,Q_PRESS,VCOMP_P_LEVS,VCOMP_PRESS,
     $     QV_IND)
*D ARS1F404.210,ARS1F404.222
      Call level_indices2(STASH_LEVELS(1,NI),
     $     Q_P_LEVS,Q_PRESS,W_P_LEVS,W_PRESS,
     $     QW_IND)
*D ARS1F404.264,ARS1F404.276
      Call level_indices2(STASH_LEVELS(1,NI),
     $     Z_P_LEVS,Z_PRESS,UCOMP_P_LEVS,UCOMP_PRESS,
     $     UZ_IND)
*D ARS1F404.293,ARS1F404.305
      Call level_indices2(STASH_LEVELS(1,NI),
     $     Z_P_LEVS,Z_PRESS,VCOMP_P_LEVS,VCOMP_PRESS,
     $     VZ_IND)
*B ST_DIA11.736
!
      CONTAINS
!
!---------------------------------------
      Subroutine level_indices2
     $     (stlev, npres1, press1, npres2, press2, levels)
!
! wrapper for level_indices where indices from two sets of levels 
! are to be extracted on the stash levels
!      
      implicit none
      integer,intent(in)::stlev(*),npres1,npres2
      real,intent(in)::press1(npres1),press2(npres2)
      integer,intent(out)::levels(stlev(1)*2)
!
      Call level_indices(stlev,npres1,press1,levels(1))
      Call level_indices(stlev,npres2,press2,levels(1+stlev(1)))
      End subroutine level_indices2
!     
!---------------------------------------
      Subroutine level_indices(stlev, npres, press, levels)
      implicit none
      integer,intent(in)::npres
      real,intent(in)::press(npres) ! pressures we have; sorted descending
!
      integer,intent(in)::stlev(*)
             ! pressures we want; sorted descending
             !  STLEV(1) is number of pressures
             !  STLEV(2:STLEV(1)+1) is integers representing pressures
!
      integer,intent(out)::levels(stlev(1))
          ! returns indices of pressures in PRESS array closest to 
          ! those specified by STLEV
!---
      integer::i,index
      real::required
      real,parameter::multiplier=2e-3
          ! factor of 1e-3 is to convert STLEV integers into pascals
          ! factor of 2 is so we can compare directly with sum of 2
          !     pressures rather than explicitly computing the mean.
!
      index=1
      do i=1,stlev(1)
         required=stlev(i+1)*multiplier
         findmatch: do while (index < npres)
           if (required > press(index)+press(index+1)) exit findmatch
           index=index+1
         enddo findmatch
         levels(i)=index
      enddo
!
!      PRINT *,'DEBUG level_indices...'
!      PRINT *,'INPUTS:'
!      PRINT *,npres,press
!      PRINT *,stlev(1:1+stlev(1))
!      PRINT *,'OUTPUT:'
!      PRINT *,levels
!      PRINT *,'============'
!      CALL FLUSH(6)
!      
      End subroutine level_indices
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  ftj1f405
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID FTJ1F405
*/ User-diagnostic code in RAD_CTL1 to provide ozone as a full-field
*/ diagnostic after LW radiation (radiation timesteps only).
*/ The stashcode is 2,260 and 3D workspace in STASHWORK is required 
*/ only if the diagnostic is requested as output.
*/ 
*/ This mod should work against vn4.3 or vn4.4 or 4.5.
*/   Author:   Tim Johns  30 Jan 1998.
*/   Reviewer: William Ingram  ??
*/ 
*/ The user STASHmaster file to go with this is:
*/ 
*/ H1| SUBMODEL_NUMBER=1
*/ H2| SUBMODEL_NAME=ATMOS
*/ H3| UM_VERSION=4.4
*/ #
*/ #|Model |Sectn | Item |Name                                |
*/ #|Space |Point | Time | Grid |LevelT|LevelF|LevelL|PseudT|PseudF|PseudL|LevCom|
*/ #| Option Codes         | Version Mask         |
*/ #|DataT |DumpP | PC1  PC2  PC3  PC4  PC5  PC6  PC7  PC8  PC9  PCA |
*/ #|Rotate| PPFC | USER | LBVC | BLEV | TLEV |RBLEVV| CFLL | CFFF |
*/ #
*/ #===============================================================================
*/ #
*/ 1|    1 |    2 |  260 |OZONE CONCENTRATION AFTER LW        |
*/ 2|    0 |    0 |    2 |    1 |    1 |    1 |    2 |    0 |    0 |    0 |    1 |
*/ 3| 00000000000000000000 | 00000000000000000101 |
*/ 4|    1 |    2 | -99  -99  -99  -99  -99  -99  -99  -99  -99  -99 |
*/ 5|    0 |  453 |    0 |    9 |    0 |    0 |    0 |    0 |    0 |
*/ #
*/ .. NB: it sets upper and lower level range to theta levels
*/        rather than available ozone levels, which didn't work.
*/
*DECLARE RAD_CTL1
*I RAD_CTL1.1004

        IF(SF(260,2)) THEN

CL 2.4.3.1  Full field ozone (expanded from zonal as necessary)

      CALL COPYDIAG_3D (STASHWORK(SI(260,2,im_index)),
     &                  OZONE_1(1,1),
     &         FIRST_POINT,LAST_POINT,P_FIELD,ROW_LENGTH,OZONE_LEVELS,
     &         STLIST(1,STINDEX(1,260,2,im_index)),LEN_STLIST,
     &         STASH_LEVELS,NUM_STASH_LEVELS+1,
     &         im_ident,2,260,
*CALL ARGPPX
     &         ICODE,CMESSAGE) 

          IF(ICODE.GT.0) THEN
            RETURN
          END IF

        END IF
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  gbc0f406
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID GBC0F406
*/ U.M. 4.6 unix / source code change form / header   version 06/01/99
*/ CODE WRITERS MUST READ THE ACCOMPANYING INSTRUCTIONS FOR THIS BUILD:
*/  - See http://fr0800/umdoc/hegui/t3e4.6.html#chgfinst
*/ 
*/SOC: Fix TIMER so that CPU Times are correct with SHMEM_NAM
*/Applies also to vn4.4
*/
*/ Has an entry been lodged in the Problem Reporting System? [Y]      
*/
*/ THIS CODE IS INTENDED FOR INCLUSION IN THE 4.6 BUILD      [Y]
*/ .....................................................................
*/   Author[s]:-> E-mail:-> bcarruthers@meto.gov.uk 
*/ Reviewer[s]:-> E-mail:-> 
*/
*/    "I have checked this change. When provided, the advance design 
*/  specification was agreed and adequate, and the new code conforms to
*/  Unified Model standards."
*/
*/  DESIGN SPEC. WAS REVIEWED ON: ......      REVIEWER[S] SIGNATURES
*/                                            ----------------------
*/    DATE CODE REVIEWED: ......
*/  .....................................................................
*/
*/  WILL CHANGES AFFECT ANCILLARY FILES?         [N]
*/  ARE ANY CHANGES TO STASHMASTER FILES NEEDED? [N] 
*/  USER INTERFACE ACTION REQUIRED?              [N]
*/ 
*/  TESTED IN CONFIGURATIONS: Global, Meso-scale
*/  TESTS RUN BY [PERSON]: Bob Carruthers
*/ 
*/  WILL THE CHANGES SLOW DOWN THE MODEL?        [N]
*/  -> Further details
*/  CHANGES WILL INCREASE MEMORY CONSUMPTION?    [N]   
*/  -> Further details
*/
*/ | Re-start dumps bit compare with those created without the change 
*/ V MARK [Y| ] BELOW; leave rest of lines untouched.
*/
*/   Control Code    loses bit comparison
*/   Atmosphere (assuming same science options chosen)   loses b.c.
*/   Ocean       loses bit comparison
*/   Wave        loses bit comparison
*/   Reconfiguration   loses bit comparison
*/   Diagnostics      lose bit comparison
*/ For Y2K compliance checking:  
*/ DOES THIS CHANGE INTERACT WITH DATE CALCULATIONS IN ANY WAY? [N]   
*/ 
*/  SECTIONS (TO BE) CHANGED:
*/
*/  SECTIONS (TO BE) DELETED? 
*/
*/  NEW SECTIONS?  Fill in form http://www-hc/~hadmk/STASHmaster_change.html,
*/  and give section numbers below:
*/  
*/  *DEFS ADDED OR REMOVED: 
*/
*/  **Existing** decks being changed [with *I, *D, *B directives]
*/
*/  TIMER3A
*/
*/  Decks being created or purged [with *DECK, *COMDECK, *PURGEDK]
*/ *......K  Deck name   Section#.vr
*/ -> 
*/ ......................................................................
*/ ANY REFERENCES TO EXTERNAL DOCUMENTS-> instead of design spec.
*/  ...OR ... ADVANCE DESIGN SPECIFICATION (optional) 
*/ ->    
*//////////////////////////////////////////////////////////////////////// 
*DC TIMER3A
*I GPB3F403.124
!LL   4.6       06/04/99  Moved the initial synchronisation in
!LL                       TIMER_OUTPUT to make sure that the
!LL                       data is not overwritten, if PE 0 is getting
!LL                       from the other PE's.
!LL                         Author: Bob Carruthers
*I TIMER3A.913
!
! Ensure that we do not change the output data area on the 
! T3E before PE 0 has had a chance to get the data, if that is the
! protocol being used
!
          CALL GC_GSYNC (nproc,info)
*D TIMER3A.933,TIMER3A.934

*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  gdr2f406
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID GDR2F406
*/ U.M. 4.6 unix / source code change form / header   version 06/01/99
*/ CODE WRITERS MUST READ THE ACCOMPANYING INSTRUCTIONS FOR THIS BUILD:
*/  - See http://fr0800/umdoc/hegui/t3e4.6.html#chgfinst
*/ 
*/SOC:-> Remove flush of Unit 6 in U_MODEL.
*/
*/->: fuller description of mod.: purpose, relevant configurations, 
*/-> applicable previous releases, dependencies
*/
*/ At 4.5, code was added to flush Unit 6 at the end if each timestep.
*/ It increases the timings in UM runs so is being removed.
*/ This modset is recommended for use in all 4.5 runs.
*/
*/ Has an entry been lodged in the Problem Reporting System? [Y]      
*/
*/ THIS CODE IS INTENDED FOR INCLUSION IN THE 4.6 BUILD      [Y]
*/ .....................................................................
*/   Author[s]: Dave Robinson E-mail: drobinson@meto.gov.uk 
*/ Reviewer[s]: Steve Mullerworth E-mail: sdmullerworth@meto.gov.uk
*/
*/    "I have checked this change. When provided, the advance design 
*/  specification was agreed and adequate, and the new code conforms to
*/  Unified Model standards."
*/
*/  DESIGN SPEC. WAS REVIEWED ON: ......      REVIEWER[S] SIGNATURES
*/                                            ----------------------
*/    DATE CODE REVIEWED: ......
*/  .....................................................................
*/
*/  WILL CHANGES AFFECT ANCILLARY FILES?         [N]
*/  ARE ANY CHANGES TO STASHMASTER FILES NEEDED? [N] 
*/  USER INTERFACE ACTION REQUIRED?              [N]
*/ 
*/  TESTED IN CONFIGURATIONS:-> Global/Climate
*/  TESTS RUN BY [PERSON]:-> D. Robinson/S Mullerworth
*/ 
*/  WILL THE CHANGES SLOW DOWN THE MODEL?        [N]
*/  -> Further details
*/  CHANGES WILL INCREASE MEMORY CONSUMPTION?    [N]   
*/  -> Further details
*/
*/ | Re-start dumps bit compare with those created without the change 
*/ V MARK [Y| ] BELOW; leave rest of lines untouched.
*/
*/   Control Code    loses bit comparison
*/   Atmosphere (assuming same science options chosen)   loses b.c.
*/   Ocean       loses bit comparison
*/   Wave        loses bit comparison
*/   Reconfiguration   loses bit comparison
*/   Diagnostics      lose bit comparison
*/ For Y2K compliance checking:  
*/ DOES THIS CHANGE INTERACT WITH DATE CALCULATIONS IN ANY WAY? [N]   
*/ 
*/  SECTIONS (TO BE) CHANGED:
*/
*/  SECTIONS (TO BE) DELETED? 
*/
*/  NEW SECTIONS?  Fill in form http://www-hc/~hadmk/STASHmaster_change.html,
*/  and give section numbers below:
*/  
*/  *DEFS ADDED OR REMOVED: 
*/
*/  **Existing** decks being changed [with *I, *D, *B directives]
*/
*/  U_MODEL1
*/
*/  Decks being created or purged [with *DECK, *COMDECK, *PURGEDK]
*/ *......K  Deck name   Section#.vr
*/
*/  None 
*/
*/ ......................................................................
*/ ANY REFERENCES TO EXTERNAL DOCUMENTS-> instead of design spec.
*/  ...OR ... ADVANCE DESIGN SPECIFICATION (optional) 
*/ ->    
*//////////////////////////////////////////////////////////////////////// 
*DECLARE U_MODEL1
*B U_MODEL1.14 
!LL  4.6  14/07/99  Remove flush of Unit 6 at end of timestep.
!LL                 D. Robinson 
*D GPB0F405.15,18
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  gdr3f406
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID GDR3F406
*/
*/ U.M. 4.6 unix / source code change form / header   version 06/01/99
*/ CODE WRITERS MUST READ THE ACCOMPANYING INSTRUCTIONS FOR THIS BUILD:
*/  - See http://fr0800/umdoc/hegui/t3e4.6.html#chgfinst
*/ 
*/SOC:-> Correct filenames of reinitialised boundary datasets.
*/
*/->: fuller description of mod.: purpose, relevant configurations, 
*/-> applicable previous releases, dependencies
*/
*/ The filenaming convention for boundary datasets was changed at
*/ Vn 4.5 (in GDR2F405) and the filenaming for re-initialised boundary
*/ datasets was set up incorrectly. This modset corrects this.
*/ Required in atmos or ocean runs setting up reinitialised
*/ boundary datasets. Not required in operational suite. 
*/
*/ Has an entry been lodged in the Problem Reporting System? [Y]      
*/
*/ THIS CODE IS INTENDED FOR INCLUSION IN THE 4.6 BUILD      [Y]
*/ .....................................................................
*/   Author[s]: Dave Robinson E-mail:-> drobinson@meto.gov.uk 
*/ Reviewer[s]:-> E-mail:-> @meto.gov.uk
*/
*/    "I have checked this change. When provided, the advance design 
*/  specification was agreed and adequate, and the new code conforms to
*/  Unified Model standards."
*/
*/  DESIGN SPEC. WAS REVIEWED ON: ......      REVIEWER[S] SIGNATURES
*/                                            ----------------------
*/    DATE CODE REVIEWED: ......
*/  .....................................................................
*/
*/  WILL CHANGES AFFECT ANCILLARY FILES?         [N]
*/  ARE ANY CHANGES TO STASHMASTER FILES NEEDED? [N] 
*/  USER INTERFACE ACTION REQUIRED?              [N]
*/ 
*/  TESTED IN CONFIGURATIONS:->
*/  TESTS RUN BY [PERSON]:->
*/ 
*/  WILL THE CHANGES SLOW DOWN THE MODEL?        [N]
*/  -> Further details
*/  CHANGES WILL INCREASE MEMORY CONSUMPTION?    [N]   
*/  -> Further details
*/
*/ | Re-start dumps bit compare with those created without the change 
*/ V MARK [Y| ] BELOW; leave rest of lines untouched.
*/
*/   Control Code    loses bit comparison
*/   Atmosphere (assuming same science options chosen)   loses b.c.
*/   Ocean       loses bit comparison
*/   Wave        loses bit comparison
*/   Reconfiguration   loses bit comparison
*/   Diagnostics      lose bit comparison
*/ For Y2K compliance checking:  
*/ DOES THIS CHANGE INTERACT WITH DATE CALCULATIONS IN ANY WAY? [Y|N]   
*/ 
*/  SECTIONS (TO BE) CHANGED:
*/
*/  SECTIONS (TO BE) DELETED? 
*/
*/  NEW SECTIONS?  Fill in form http://www-hc/~hadmk/STASHmaster_change.html,
*/  and give section numbers below:
*/  
*/  *DEFS ADDED OR REMOVED: 
*/
*/  **Existing** decks being changed [with *I, *D, *B directives]
*/ 
*/  GET_NAM2 
*/
*/  Decks being created or purged [with *DECK, *COMDECK, *PURGEDK]
*/ *......K  Deck name   Section#.vr
*/
*/  None 
*/ ......................................................................
*/ ANY REFERENCES TO EXTERNAL DOCUMENTS-> instead of design spec.
*/  ...OR ... ADVANCE DESIGN SPECIFICATION (optional) 
*/ 
*/  None
*/    
*////////////////////////////////////////////////////////////////////////
*/
*DECLARE GET_NAM2
*B TJ080294.479
!LL  4.6  19/07/99  Correct new filenaming convention for reinitialised
!LL                 boundary files. D. Robinson.
*D DR240293.600,GET_NAM2.247
                  FILETYPE_2=LETTER_3
*D GET_NAM2.255
                  FILETYPE_2=LETTER_3
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  gnuml45.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*IDENT GNUML45
*/
*/ Update to set gnu to FKAPB (i.e. to set vert tracer diffusivity 
*/ to its background value) everywhere.
*/
*/ Written for UM vn 3.1              RAW 09/11/94
*/ Converted for vn 3.2               RAW 3/3/95
*/ Should work OK for vn4.3           JRP 27/6/97
*/
*/ Re-written to replace all the executable code in VERTCOFT.
*/	Annette Osprey 15/09/06
*/ 
*/ This is a FAMOUS mod. 
*/ 
*/
*DECLARE VERTCOFT
*D VERTCOFT.60,VERTCOFT.157
C
      DO K=1,KM
        DO I=1,IMT
          gnu(I,K)=KAPPA_B_SI(K)*1.E4*FM(I,K) ! use background value
        END DO 
      END DO
C*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  gps0f406
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID GPS0F406
*/ U.M. 4.6 unix / source code change form / header   version 04/01/99
*/Instructions: see http://fr0800/umdoc/hegui/t3e4.6.html#chgfinst
*/ 
*/SOC: Fix STASH accumulations/means with times list.
*/ 
*/ Sets st_end_time_code correctly in PRELIM1 for time-processing
*/ diagnostics (*not* time series though) that use a times list
*/ rather than time frequencies.
*/
*/ Has an entry been lodged in the Problem Reporting System? [Y]      
*/
*/ THIS CODE IS INTENDED FOR INCLUSION IN THE 4.6 BUILD      [Y]
*/ .....................................................................
*/   Author[s]: Paul Selwood         E-mail: pmselwood@meto.gov.uk 
*/ Reviewer[s]:-> E-mail:-> @meto.gov.uk
*/
*/    "I have checked this change. When provided, the advance design 
*/  specification was agreed and adequate, and the new code conforms to
*/  Unified Model standards."
*/
*/  DESIGN SPEC. WAS REVIEWED ON: ......      REVIEWER[S] SIGNATURES
*/                                            ----------------------
*/    DATE CODE REVIEWED: ......
*/  .....................................................................
*/
*/  WILL CHANGES AFFECT ANCILLARY FILES?         [N]
*/  ARE ANY CHANGES TO STASHMASTER FILES NEEDED? [N] 
*/  USER INTERFACE ACTION REQUIRED?              [N]
*/ 
*/  TESTED IN CONFIGURATIONS:->
*/  TESTS RUN BY [PERSON]:->
*/ 
*/  WILL THE CHANGES SLOW DOWN THE MODEL?        [N]
*/  -> Further details
*/  CHANGES WILL INCREASE MEMORY CONSUMPTION?    [N]   
*/  -> Further details
*/
*/ | Re-start dumps bit compare with those created without the change 
*/ V MARK [Y| ] BELOW; leave rest of lines untouched.
*/
*/   Control Code    loses bit comparison
*/   Atmosphere (assuming same science options chosen)   loses b.c.
*/   Ocean       loses bit comparison
*/   Wave        loses bit comparison
*/   Reconfiguration   loses bit comparison
*/   Diagnostics      lose bit comparison
*/ For Y2K compliance checking:  
*/ DOES THIS CHANGE INTERACT WITH DATE CALCULATIONS IN ANY WAY? [N]   
*/ 
*/  SECTIONS (TO BE) CHANGED:
*/
*/  SECTIONS (TO BE) DELETED? 
*/
*/  NEW SECTIONS?  Fill in form http://www-hc/~hadmk/STASHmaster_change.html,
*/  and give section numbers below:
*/  
*/  *DEFS ADDED OR REMOVED: 
*/
*/  **Existing** decks being changed [with *I, *D, *B directives]
*/  PRELIM1
*/
*/  Decks being created or purged [with *DECK, *COMDECK, *PURGEDK]
*/ *......K  Deck name   Section#.vr
*/ -> 
*/ ......................................................................
*/ ANY REFERENCES TO EXTERNAL DOCUMENTS-> instead of design spec.
*/  ...OR ... ADVANCE DESIGN SPECIFICATION (optional) 
*/ ->    
*//////////////////////////////////////////////////////////////////////// 
*//-------------
*DECLARE PRELIM1
*//-------------
*I ABX1F405.2
!   4.6    5/1/99      Correct st_end_time_code when using time
!                      processing                    P.Selwood
*I PRELIM1.92
      INTEGER      IMAX          ! to find max of times-table
      INTEGER      ITIMLST       ! column of times-table
*I ABX1F405.42
          END IF                                                        
                                                                        
    ! For the NRECS item an end_time_code needs to be set if we
    ! are dealing with a times table rather than  regular diagn.
    ! This should be the maximum timestep in the time list. The list
    ! should be ready sorted (and thus maximum is last member) but
    ! will run through and find maximum to be on the safe side.

          IMAX = 0
          ITIMLST = -LIST_S(st_freq_code,NRECS+1)

          IF (ITIMLST .GT. 0) THEN      ! List *not* regular
            DO I = 1, ITIM_T(ITIMLST)
              IF (IMAX .LT. ITIM_S(I, ITIMLST)) THEN
                IMAX = ITIM_S(I, ITIMLST)
              END IF
            END DO

            LIST_S(st_end_time_code,NRECS) = IMAX

*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  gsm1f406
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID GSM1F406
*/ U.M. 4.6 unix / source code change form / header   version 08/12/98
*/Instructions: see http://fr0800/umdoc/hegui/t3e4.6.html#chgfinst
*/ 
*/ SOC: Change TIMER call to prevent barrier hang.
*/
*/ Advised to use with Hydrology 5A (Section 8 Moses I, HadCM3 option)
*/ Strictly only required if using TIMER (LTIMER = .TRUE. in CNTLALL
*/ and CONTCNTL) but has no side effects if not.
*/ Prevents failure when certain PEs have no landpoints (eg on 6x4 pes)
*/ 26/11/98
*/
*/ Has an entry been lodged in the Problem Reporting System? [Y]
*/
*/ THIS CODE IS INTENDED FOR INCLUSION IN THE 4.6 BUILD      [Y]
*/ .....................................................................
*/   Author[s]:-> Steve Mullerworth E-mail:-> sdmullerworth@meto.gov.uk 
*/ Reviewer[s]:-> E-mail:-> @meto.gov.uk
*/
*/    "I have checked this change. When provided, the advance design 
*/  specification was agreed and adequate, and the new code conforms to
*/  Unified Model standards."
*/
*/  DESIGN SPEC. WAS REVIEWED ON: ......      REVIEWER[S] SIGNATURES
*/                                            ----------------------
*/    DATE CODE REVIEWED: ......
*/  .....................................................................
*/ ANY REFERENCES TO EXTERNAL DOCUMENTS-> instead of design spec.
*/
*/  WILL CHANGES AFFECT ANCILLARY FILES?         [N]
*/  ARE ANY CHANGES TO STASHMASTER FILES NEEDED? [N] 
*/  USER INTERFACE ACTION REQUIRED?              [N]
*/ 
*/  TESTED IN CONFIGURATIONS:-> HadCM3
*/  TESTS RUN BY [PERSON]:-> Cyndy Bunton
*/ 
*/  WILL THE CHANGES SLOW DOWN THE MODEL?        [N]
*/  -> Further details
*/  CHANGES WILL INCREASE MEMORY CONSUMPTION?    [N]   
*/  -> Further details
*/
*/ | Forecast dumps bit compare with those created without the change 
*/ V MARK [Y| ] BELOW; leave rest of lines untouched.
*/
*/   Control Code    loses bit comparison
*/   Atmosphere (assuming same science options chosen)   loses b.c.
*/   Ocean       loses bit comparison
*/   Wave        loses bit comparison
*/   Reconfiguration   loses bit comparison
*/   Diagnostics      lose bit comparison
*/ For Y2K compliance checking:  
*/ DOES THIS CHANGE INTERACT WITH DATE CALCULATIONS IN ANY WAY? [N]   
*/ 
*/  SECTIONS (TO BE) CHANGED:
*/
*/  SECTIONS (TO BE) DELETED? 
*/
*/  NEW SECTIONS?  Fill in form http://www-hc/~hadmk/STASHmaster_change.html,
*/  and give section numbers below:
*/  
*/  *DEFS ADDED OR REMOVED: 
*/ ......................................................................
*/                ADVANCE DESIGN SPECIFICATION (optional) 
*/ ->    
*//////////////////////////////////////////////////////////////////////// 
*DECLARE HEATCP5A
*I ACB2F405.2         
!  4.5   11/12/98    TIMER calls changed. S.D.Mullerworth
*D HEATCP5A.113
        CALL TIMER('HEATCAP ',103)
*D HEATCP5A.175
        CALL TIMER('HEATCAP ',104)
*DECLARE SFSNOW2A
*I ARE2F404.490
!  4.5   11/12/98    TIMER calls changed. S.D.Mullerworth
*D AJS1F401.1519
        CALL TIMER('SFSNOW  ',103)
*D AJS1F401.1522
        CALL TIMER('SFSNOW  ',104)
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  gsm2f406
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID GSM2F406
*/ U.M. 4.6 unix / source code change form / header   version 08/12/98
*/Instructions: see http://fr0800/umdoc/hegui/t3e4.6.html#chgfinst
*/ 
*/ SOC: Remove most of LSMASK=.TRUE. lines from TSTMSK
*/
*/ This modset can be used in all configurations. It rejects requests
*/ for some unavailable STASH items where before, the request was not
*/ rejected.
*/ TSTMSK initialises LSMASK to TRUE at top, so no other LSMASK=.TRUE.
*/ are required. Some that were removed were bugs that reversed a
*/ preceding setting to false (eg GSS3F401.1462 and ACN2F405.14)
*/ others that were removed were superfluous and confusing - people
*/ making changes were copying the logic and introducing bugs.
*/ 15/01/99
*/
*/ Has an entry been lodged in the Problem Reporting System? [Y]
*/
*/ THIS CODE IS INTENDED FOR INCLUSION IN THE 4.6 BUILD      [Y]
*/ .....................................................................
*/   Author[s]:-> Steve Mullerworth E-mail:-> sdmullerworth@meto.gov.uk 
*/ Reviewer[s]:-> E-mail:-> @meto.gov.uk
*/
*/    "I have checked this change. When provided, the advance design 
*/  specification was agreed and adequate, and the new code conforms to
*/  Unified Model standards."
*/
*/  DESIGN SPEC. WAS REVIEWED ON: ......      REVIEWER[S] SIGNATURES
*/                                            ----------------------
*/    DATE CODE REVIEWED: ......
*/  .....................................................................
*/ ANY REFERENCES TO EXTERNAL DOCUMENTS-> instead of design spec.
*/
*/  WILL CHANGES AFFECT ANCILLARY FILES?         [N]
*/  ARE ANY CHANGES TO STASHMASTER FILES NEEDED? [N] 
*/  USER INTERFACE ACTION REQUIRED?              [N]
*/ 
*/  TESTED IN CONFIGURATIONS:-> HadCM3
*/  TESTS RUN BY [PERSON]:-> S.D.Mullerworth
*/ 
*/  WILL THE CHANGES SLOW DOWN THE MODEL?        [N]
*/  -> Further details
*/  CHANGES WILL INCREASE MEMORY CONSUMPTION?    [N]   
*/  -> Further details
*/
*/ | Forecast dumps bit compare with those created without the change 
*/ V MARK [Y| ] BELOW; leave rest of lines untouched.
*/
*/   Control Code    loses bit comparison
*/   Atmosphere (assuming same science options chosen)   loses b.c.
*/   Ocean       loses bit comparison
*/   Wave        loses bit comparison
*/   Reconfiguration   loses bit comparison
*/   Diagnostics      lose bit comparison
*/ For Y2K compliance checking:  
*/ DOES THIS CHANGE INTERACT WITH DATE CALCULATIONS IN ANY WAY? [N]   
*/ 
*/  SECTIONS (TO BE) CHANGED:
*/
*/  SECTIONS (TO BE) DELETED? 
*/
*/  NEW SECTIONS?  Fill in form http://www-hc/~hadmk/STASHmaster_change.html,
*/  and give section numbers below:
*/  
*/  *DEFS ADDED OR REMOVED: 
*/ ......................................................................
*/                ADVANCE DESIGN SPECIFICATION (optional) 
*/ ->    
*//////////////////////////////////////////////////////////////////////// 
*DECLARE TSTMSK1
*I OOM1F405.12    
!  4.6  15/01/98  Get rid of all the LMASK=.TRUE. after the first.
!                 S.D.Mullerworth
*D GSS3F401.1419,GSS3F401.1424 
              END IF                                                    
*D GSS3F401.1434,GSS3F401.1436 
        IF (SUM_IOPN.NE.0) THEN                                         
*D GSS3F401.1443,GSS3F401.1444 
*D GSS3F401.1449,GSS3F401.1451 
        IF (SUM_IOPN.NE.0) THEN                                         
*D GSS3F401.1456,GSS3F401.1457 
*D GSS3F401.1461,GSS3F401.1462 
*D ACN2F405.13,ACN2F405.14   
*D GSS3F401.1469,GSS3F401.1470 
*D GSS3F401.1476,GSS3F401.1477 
*D ASK1F405.36,ASK1F405.37   
*D GSS3F401.1481,GSS3F401.1483 
        IF (SUM_IOPN.NE.0) THEN                   
*D GSS3F401.1487,GSS3F401.1488 
*D AWO0F405.59,AWO0F405.60   
*D AWO0F405.66,AWO0F405.67   
*D GSS3F401.1493,GSS3F401.1495 
        IF (SUM_IOPN.NE.0) THEN      
*D GSS3F401.1502,GSS3F401.1503 
*D GSS3F401.1507,GSS3F401.1509 
        IF (SUM_IOPN.NE.0) THEN      
*D GSS3F401.1597,GSS3F401.1601 
            END IF            
*D GSS3F401.1671,GSS3F401.1675 
              END IF                                                    
*D GSS3F401.1679,GSS3F401.1682 
        IF (SUM_IOPN.NE.0) THEN                                         
          N2N1=MOD( IOPN(1),100)                                          
*D GSS3F401.1738,GSS3F401.1739 
*D GSS3F401.1744,GSS3F401.1746 
        IF (SUM_IOPN.NE.0) THEN                                         
*D GSS3F401.1763,GSS3F401.1764 
*D GSS3F401.1773,GSS3F401.1776 
          IF (SUM_IOPN.NE.0) THEN                                
            N2N1=MOD( IOPN(1),100)                                        
*D GSS3F401.1834,GSS3F401.1835 
*D GSS3F401.1846,GSS3F401.1848 
        IF (SUM_IOPN.NE.0) THEN                                         
*D GSS3F401.1864,GSS3F401.1865 
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  gsm4f406.PUM
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*IDENT GSM4F406
*/ U.M. 4.6 unix / source code change form / header   version 08/12/98
*/Instructions: see http://fr0800/umdoc/hegui/t3e4.6.html#chgfinst
*/ 
*/ SOC: Add oneover_v call to CONVEC3C for T3E
*/
*/ Optimise calculation of reciprocal of PSTAR
*/ S.D.Mullerworth 19/02/99
*/
*/ Has an entry been lodged in the Problem Reporting System? [N]
*/
*/ THIS CODE IS INTENDED FOR INCLUSION IN THE 4.6 BUILD      [Y]
*/ .....................................................................
*/   Author[s]:-> Steve Mullerworth E-mail:-> sdmullerworth@meto.gov.uk 
*/ Reviewer[s]:-> E-mail:-> @meto.gov.uk
*/
*/    "I have checked this change. When provided, the advance design 
*/  specification was agreed and adequate, and the new code conforms to
*/  Unified Model standards."
*/
*/  DESIGN SPEC. WAS REVIEWED ON: ......      REVIEWER[S] SIGNATURES
*/                                            ----------------------
*/    DATE CODE REVIEWED: ......
*/  .....................................................................
*/ ANY REFERENCES TO EXTERNAL DOCUMENTS-> instead of design spec.
*/
*/  WILL CHANGES AFFECT ANCILLARY FILES?         [N]
*/  ARE ANY CHANGES TO STASHMASTER FILES NEEDED? [N] 
*/  USER INTERFACE ACTION REQUIRED?              [N]
*/ 
*/  TESTED IN CONFIGURATIONS:-> HadCM3,HadAM3,HadSM3
*/  TESTS RUN BY [PERSON]:-> S.D.Mullerworth
*/ 
*/  WILL THE CHANGES SLOW DOWN THE MODEL?        [N]
*/  -> Further details
*/  CHANGES WILL INCREASE MEMORY CONSUMPTION?    [N]   
*/  -> Further details
*/
*/ | Forecast dumps bit compare with those created without the change 
*/ V MARK [Y| ] BELOW; leave rest of lines untouched.
*/
*/   Control Code    loses bit comparison
*/ Y Atmosphere (assuming same science options chosen)   loses b.c.
*/   Ocean       loses bit comparison
*/   Wave        loses bit comparison
*/   Reconfiguration   loses bit comparison
*/   Diagnostics      lose bit comparison
*/ For Y2K compliance checking:  
*/ DOES THIS CHANGE INTERACT WITH DATE CALCULATIONS IN ANY WAY? [N]   
*/ 
*/  SECTIONS (TO BE) CHANGED:
*/  A05_3C
*/  SECTIONS (TO BE) DELETED? 
*/
*/  NEW SECTIONS?  Fill in form http://www-hc/~hadmk/STASHmaster_change.html,
*/  and give section numbers below:
*/  
*/  *DEFS ADDED OR REMOVED: T3E added to CONVEC3C
*/ ......................................................................
*/                ADVANCE DESIGN SPECIFICATION (optional) 
*/ ->    
*//////////////////////////////////////////////////////////////////////// 
*DECLARE CONVEC3C
*I ADR1F405.28    
!LL  4.6   19/02/99  Replace reciprocal with oneover_v. S.D.Mullerworth
*D CONVEC3C.796
     *         DQS_DTH,COR_ENGY,DD_CALL,CALC_3D_CCA,oneover_v
*I CONVEC3C.814   
*IF DEF,VECTLIB
      CALL oneover_v(NPNTS,PSTAR,RECIP_PSTAR)
*ELSE
*I CONVEC3C.817   
*ENDIF
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  gsm5f406
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID GSM5F406
*/ U.M. 4.6 unix / source code change form / header   version 08/12/98
*/Instructions: see http://fr0800/umdoc/hegui/t3e4.6.html#chgfinst
*/ 
*/ SOC: Correct sign on OFFSET_DUMPSim
*/
*/ Only affects runs using climate meaning whose mean reference time
*/ does not coincide with a period 1 mean period. Without this mod,
*/ the history file gets written at the wrong times which makes the
*/ run unrestartable, and each automatically resubmitted chunk may
*/ run for the wrong period.
*/ 01/03/99
*/
*/ Has an entry been lodged in the Problem Reporting System? [Y]
*/
*/ THIS CODE IS INTENDED FOR INCLUSION IN THE 4.6 BUILD      [Y]
*/ .....................................................................
*/   Author[s]:-> Steve Mullerworth E-mail:-> sdmullerworth@meto.gov.uk 
*/ Reviewer[s]:-> E-mail:-> @meto.gov.uk
*/
*/    "I have checked this change. When provided, the advance design 
*/  specification was agreed and adequate, and the new code conforms to
*/  Unified Model standards."
*/
*/  DESIGN SPEC. WAS REVIEWED ON: ......      REVIEWER[S] SIGNATURES
*/                                            ----------------------
*/    DATE CODE REVIEWED: ......
*/  .....................................................................
*/ ANY REFERENCES TO EXTERNAL DOCUMENTS-> instead of design spec.
*/
*/  WILL CHANGES AFFECT ANCILLARY FILES?         [N]
*/  ARE ANY CHANGES TO STASHMASTER FILES NEEDED? [N] 
*/  USER INTERFACE ACTION REQUIRED?              [N]
*/ 
*/  TESTED IN CONFIGURATIONS:-> HadAM3
*/  TESTS RUN BY [PERSON]:-> Steve Mullerworth
*/ 
*/  WILL THE CHANGES SLOW DOWN THE MODEL?        [N]
*/  -> Further details
*/  CHANGES WILL INCREASE MEMORY CONSUMPTION?    [N]   
*/  -> Further details
*/
*/ | Forecast dumps bit compare with those created without the change 
*/ V MARK [Y| ] BELOW; leave rest of lines untouched.
*/
*/   Control Code    loses bit comparison
*/   Atmosphere (assuming same science options chosen)   loses b.c.
*/   Ocean       loses bit comparison
*/   Wave        loses bit comparison
*/   Reconfiguration   loses bit comparison
*/   Diagnostics      lose bit comparison
*/ For Y2K compliance checking:  
*/ DOES THIS CHANGE INTERACT WITH DATE CALCULATIONS IN ANY WAY? [N]   
*/ 
*/  SECTIONS (TO BE) CHANGED:
*/
*/  SECTIONS (TO BE) DELETED? 
*/
*/  NEW SECTIONS?  Fill in form http://www-hc/~hadmk/STASHmaster_change.html,
*/  and give section numbers below:
*/  
*/  *DEFS ADDED OR REMOVED: 
*/ ......................................................................
*/                ADVANCE DESIGN SPECIFICATION (optional) 
*/ ->    
*//////////////////////////////////////////////////////////////////////// 
*DECLARE SETTSCT1
*I GMB1F405.540   
!LL  4.6  01/03/99  Correct sign on OFFSET_DUMPSim. S.D.Mullerworth
*D GRB1F305.592
          LHISTORY=  (MOD(STEP+OFFSET_DUMPS*DUMPFREQ,
*DC INITTIM1
*I GSM2F405.15    
!LL  4.6  01/03/99  Correct sign on OFFSET_DUMPSim. S.D.Mullerworth
*D GJC0F405.26
      IF (TARGET_END_STEP .GE. (-OFFSET_DUMPSim(im) *                   
*D GGH1F401.22
     &      MOD(TARGET_END_STEP + (OFFSET_DUMPSim(im)*DUMPFREQim(im)),
*D GGH1F401.26
     &      MOD(TARGET_END_STEP + (OFFSET_DUMPSim(im)*DUMPFREQim(im)),
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  gsm8f406
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*IDENT GSM8F406
*/ U.M. 4.6 unix / source code change form / header   version 08/12/98
*/Instructions: see http://fr0800/umdoc/hegui/t3e4.6.html#chgfinst
*/ 
*/ SOC: Reduce MEANDIAG buffer size - for very high res runs
*/
*/ MEANDIAG buf array is dimensioned for maximum size of all levels
*/ of a field when it only deals with one level at a time. This is
*/ corrected by using ST_LIST(st_dump_level_output_length,*) instead of 
*/ st_dump_output_length.
*/   Used in, eg, 1/3 degree ocean where a 20 level ocean field is
*/ about 8MW.
*/
*/ Has an entry been lodged in the Problem Reporting System? [Y]
*/
*/ THIS CODE IS INTENDED FOR INCLUSION IN THE 4.6 BUILD      [Y]
*/ .....................................................................
*/   Author[s]:-> Steve Mullerworth E-mail:-> sdmullerworth@meto.gov.uk 
*/ Reviewer[s]:-> E-mail:-> @meto.gov.uk
*/
*/    "I have checked this change. When provided, the advance design 
*/  specification was agreed and adequate, and the new code conforms to
*/  Unified Model standards."
*/
*/  DESIGN SPEC. WAS REVIEWED ON: ......      REVIEWER[S] SIGNATURES
*/                                            ----------------------
*/    DATE CODE REVIEWED: ......
*/  .....................................................................
*/ ANY REFERENCES TO EXTERNAL DOCUMENTS-> instead of design spec.
*/
*/  WILL CHANGES AFFECT ANCILLARY FILES?         [N]
*/  ARE ANY CHANGES TO STASHMASTER FILES NEEDED? [N] 
*/  USER INTERFACE ACTION REQUIRED?              [N]
*/ 
*/  TESTED IN CONFIGURATIONS:-> HadAM3
*/  TESTS RUN BY [PERSON]:-> Steve Mullerworth
*/ 
*/  WILL THE CHANGES SLOW DOWN THE MODEL?        [N]
*/  -> Further details
*/  CHANGES WILL INCREASE MEMORY CONSUMPTION?    [N]   
*/  -> Further details
*/
*/ | Forecast dumps bit compare with those created without the change 
*/ V MARK [Y| ] BELOW; leave rest of lines untouched.
*/
*/   Control Code    loses bit comparison
*/   Atmosphere (assuming same science options chosen)   loses b.c.
*/   Ocean       loses bit comparison
*/   Wave        loses bit comparison
*/   Reconfiguration   loses bit comparison
*/   Diagnostics      lose bit comparison
*/ For Y2K compliance checking:  
*/ DOES THIS CHANGE INTERACT WITH DATE CALCULATIONS IN ANY WAY? [N]   
*/ 
*/  SECTIONS (TO BE) CHANGED:
*/
*/  SECTIONS (TO BE) DELETED? 
*/
*/  NEW SECTIONS?  Fill in form http://www-hc/~hadmk/STASHmaster_change.html,
*/  and give section numbers below:
*/  
*/  *DEFS ADDED OR REMOVED: 
*/ ......................................................................
*/                ADVANCE DESIGN SPECIFICATION (optional) 
*/ ->    
*//////////////////////////////////////////////////////////////////////// 
*DECLARE MEANCTL1
*B MEANCTL1.14
!LL  4.5+ 30/09/99  Correct oversizing of MAXSIZE. S.D.Mullerworth
*D GPB2F405.97
     &                      STLIST(st_dump_level_output_length,IE))
*D GPB2F405.114
     &                      STLIST(st_dump_level_output_length,IE))

*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  gsm9f406
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*IDENT GSM9F406
*/ U.M. 4.6 unix / source code change form / header   version 08/12/98
*/Instructions: see http://fr0800/umdoc/hegui/t3e4.6.html#chgfinst
*/ 
*/ SOC: Add declaration for functions to COEX1A
*/
*/ There is no "IMPLICIT NONE" in COEX1A which in some cases results
*/ in CRI2IBM returning a real value to an integer.
*/
*/ Has an entry been lodged in the Problem Reporting System? [Y]
*/
*/ THIS CODE IS INTENDED FOR INCLUSION IN THE 4.6 BUILD      [Y]
*/ .....................................................................
*/   Author[s]:-> Steve Mullerworth E-mail:-> sdmullerworth@meto.gov.uk 
*/ Reviewer[s]:-> E-mail:-> @meto.gov.uk
*/
*/    "I have checked this change. When provided, the advance design 
*/  specification was agreed and adequate, and the new code conforms to
*/  Unified Model standards."
*/
*/  DESIGN SPEC. WAS REVIEWED ON: ......      REVIEWER[S] SIGNATURES
*/                                            ----------------------
*/    DATE CODE REVIEWED: ......
*/  .....................................................................
*/ ANY REFERENCES TO EXTERNAL DOCUMENTS-> instead of design spec.
*/
*/  WILL CHANGES AFFECT ANCILLARY FILES?         [N]
*/  ARE ANY CHANGES TO STASHMASTER FILES NEEDED? [N] 
*/  USER INTERFACE ACTION REQUIRED?              [N]
*/ 
*/  TESTED IN CONFIGURATIONS:-> HadAM3
*/  TESTS RUN BY [PERSON]:-> David Hassel
*/ 
*/  WILL THE CHANGES SLOW DOWN THE MODEL?        [N]
*/  -> Further details
*/  CHANGES WILL INCREASE MEMORY CONSUMPTION?    [N]   
*/  -> Further details
*/
*/ | Forecast dumps bit compare with those created without the change 
*/ V MARK [Y| ] BELOW; leave rest of lines untouched.
*/
*/   Control Code    loses bit comparison
*/   Atmosphere (assuming same science options chosen)   loses b.c.
*/   Ocean       loses bit comparison
*/   Wave        loses bit comparison
*/   Reconfiguration   loses bit comparison
*/   Diagnostics      lose bit comparison
*/ For Y2K compliance checking:  
*/ DOES THIS CHANGE INTERACT WITH DATE CALCULATIONS IN ANY WAY? [N]   
*/ 
*/  SECTIONS (TO BE) CHANGED:
*/
*/  SECTIONS (TO BE) DELETED? 
*/
*/  NEW SECTIONS?  Fill in form http://www-hc/~hadmk/STASHmaster_change.html,
*/  and give section numbers below:
*/  
*/  *DEFS ADDED OR REMOVED: 
*/ ......................................................................
*/                ADVANCE DESIGN SPECIFICATION (optional) 
*/ ->    
*//////////////////////////////////////////////////////////////////////// 
*DECLARE COEX1A
*I UIE2F403.3
!LL   4.6  29/11/99 Add declaration for CRI2IBM. S.D.Mullerworth
*I COEX1A.28    
      INTEGER CRI2IBM,IBM2CRI
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  h3dbqlim
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/
*/  copied from ~t20db/um_mods1/h3dbqlim  18/3/96
*/  unchanged at vn4.4                    13/01/98
*/  unchanged at vn4.5                    15/01/99
*/ Lower limit on specific humidity.
*/
*DECLARE FILL3A
*D FILL3A.150
      PARAMETER(H2OLMN=2.5E-6)
*D FILL3A.206
            GAS_MIX_RATIO(L, I, IUMP_H2O)
     &         =MAX(H2O(LG, NLEVS+1-I), H2OLMN)
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  inittime_info
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID EXTRA INFO
*DECLARE INITTIM1
*I INITTIM1.130
      write(6,*)' Model basis time = '
      write(6,*)(model_basis_time(i),i=1,6)
      write(6,*)' Atmos. dump headers 28,33= '
      write(6,*)(a_fixhd(i),i=28,33)
*I INITTIM1.143
      write(6,*)' Model basis time = '
      write(6,*)(model_basis_time(i),i=1,6)
      write(6,*)' Ocean dump headers 28,33= '
      write(6,*)(o_fixhd(i),i=28,33)
*I INITTIM1.162
      write(6,*)' Model data time = '
      write(6,*)(model_data_time(i),i=1,6)
      write(6,*)' Atmos. dump headers 21,26 = '
      write(6,*)(a_fixhd(i),i=21,36)
*I INITTIM1.175
      write(6,*)' Model data time = '
      write(6,*)(model_data_time(i),i=1,6)
      write(6,*)' ocean dump headers 21,26= '
      write(6,*)(o_fixhd(i),i=21,26)
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  linuxf_mpp.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID GCB3F404
*/ AV - one of Bob's mods for SGI. MPI must be initialised before
*/      opening any files. On the T3E, this doesn't matter.
*/
*/ AH - mod to delete lines relating to openwing + sending of 
*/      wakeup call to channel 8
*/      Also to delete closing of channel 8.  With these lines 
*/      intact the Linux mpi code crashes after the last timestep
*/      with a p4 error.
*/
*DC UMSHELL1
*D GPB0F305.5,GPB0F305.15
*D GRR2F305.293,GRR2F305.297

*I GPB0F305.163

!   Open file for UNIT 5 before initialisation of model. All runtime
!   control variables subsequently read in from UNIT 5 by namelist.
      CALL GET_FILE(5,FILENAME,80,ICODE)
      OPEN(5,FILE=FILENAME,IOSTAT=ISTATUS)
      IF(ISTATUS.NE.0) THEN
        ICODE=500
        WRITE(6,*) ' ERROR OPENING FILE ON UNIT 5'
        WRITE(6,*) ' FILENAME =',FILENAME
        WRITE(6,*) ' IOSTAT =',ISTATUS
        GOTO 999
      END IF

CL------------------------------------------------------------------
CL 0.1 Get submodel/internal model components of model run.
CL
      ICODE=0
      CALL UM_Submodel_Init(ICODE)
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  long_output_names.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*IDENT RENAME
*/ modified version of Julia Tindall's rename_output.f mod to 
*/ cope with long filenames in paleoruns. As Absolute Long 
*/ format, but uses the year in an (i9.9) format in the 
*/ filename, with a + or - appended for past/future runs. The 
*/ middle character is changed from @ to # to indicate the new 
*/ format. Needs to be used with modified qx execs for full
*/ continuation run capabilities.
*/ - robin smith + jonathan gregory, 2007
*DECLARE DUMPCTL1
*D DUMPCTL1.80
      CHARACTER*21 DUMPNAME
*DECLARE MEANCTL1
*D MEANCTL1.85
*/ change size of files
      CHARACTER*21 PSNAME_READ,PSNAME_WRITE,PSNAME_DELETE
*D GKR1F404.356
 890      FORMAT('%%% ',A21,' DELETE')                                     
*DECLARE U_MODEL1
*D U_MODEL1.68
      CHARACTER*21 PPNAME       ! Work - Dummy PP filename   
*DECLARE INITIAL1
*D INITIAL1.93
*/ change size of file
      CHARACTER*21 PPNAME       ! Dummy PP file name returned by PPCTL  
*DECLARE MEANDIA2
*D MEANDIA2.105
*/ change size of file
      CHARACTER*21
*D MEANDIA2.507,MEANDIA2.511
 410  FORMAT("%%% ",A21," ARCHIVE PPCHART")
 420  FORMAT("%%% ",A21," ARCHIVE PPNOCHART")
 430  FORMAT("%%% ",A21," PLOTONLY PPCHART")
 440  FORMAT("%%% ",A21," DELETE")
 450  FORMAT("%%% ",A21," HPSEND")                                         
*DECLARE PPCTL2
*D PPCTL2.46
*/ change size of file
      CHARACTER*21 PPNAME       ! OUT - PP filename generated by GET_NAME 
*D PPCTL2.74
      CHARACTER*21 OLDPPFILE        ! Previous PPfile on given unit
*D PPCTL2.154,PPCTL2.155
      STRING(18:38)=PPNAME                       
      STRING(39:80)='                                       '
*D PPCTL2.281 
      OLDPPFILE=STRING(IPOS+1:IPOS+21)
*D PPCTL2.330,PPCTL2.335 
 100  FORMAT('%%% ',A21,' ARCHIVE PPCHART')
 110  FORMAT('%%% ',A21,' ARCHIVE PPNOCHART')
 120  FORMAT('%%% ',A21,' PLOTONLY PPCHART')
 130  FORMAT('%%% ',A21,' ARCHIVE BNDY')
 140  FORMAT('%%% ',A21,' HPSEND')
 150  FORMAT('%%% ',A21,' DELETE')
*D PPCTL2.372,PPCTL2.373
      STRING(18:38)=PPNAME        
      STRING(39:80)='                                       '   
*D PPCTL2.482,PPCTL2.483
      STRING(18:38)=PPNAME                 
      STRING(39:80)='                                       '      
*DECLARE GET_NAM2
*I GET_NAM2.71
      INTEGER absyear
*I GET_NAM2.86
      CHARACTER*1 PorM          ! Character plus/minus identifier
*D GET_NAM2.53
*/ change size of filename
      CHARACTER*21  FILENAME    ! OUT - Generated file name
*I GET_NAM2.178
      ELSE IF (TIME_CONVENTION.EQ.'Absolute_verylong') THEN
          SEPARATOR='#'
          STYLE='B'
*I GET_NAM2.377
      ELSE IF (TIME_CONVENTION.EQ.'Absolute_verylong') THEN
            absyear=abs(YYYY)
            write(FILENAME(10:18),'(i9.9)')absyear
            
            FILENAME(19:19)=M                                  
            FILENAME(20:20)=D   
            PorM                = '+'
            IF (YYYY.LT.0) PorM = '-'
            FILENAME(21:21)=PorM
*I GET_NAM2.437    
      ELSE IF (TIME_CONVENTION.EQ.'Absolute_verylong') THEN
            absyear=abs(YYYY)
            write(FILENAME(10:18),'(i9.9)')absyear

            FILENAME(19:19)=M                                  
            FILENAME(20:20)=D   
            PorM                   = '+'
            IF (YYYY.LT.0) PorM    = '-'
            FILENAME(21:21)=PorM

            IF (MEANLEV.eq.0) THEN                                     
              FILENAME(19:20) = MONTH_2CHAR(MM)                
            ELSE ! means date routine called is at beginning of next     
C                   period                                               
             MON=MM-1                                                   
             IF (MON.EQ.0) THEN                                    
               MON = 12                                            
             ENDIF                                                 
             IF (FILETYPE_2.eq.'m') THEN                           
               FILENAME(19:20) = MONTH_2CHAR(MON)                  
             ELSE                                                  
               FILENAME(19:20) = SEASON_2CHAR(MON)                 
             ENDIF                                                 
            ENDIF                                                  
*DECLARE SEC2TIM1
*I SEC2TIM1.81 
            IF (I_DAY.GE.0) THEN
*I SEC2TIM1.85 
            ELSE
               I_YEAR=((I_DAY+1)/360)-1
               I_DAY=I_DAY - (I_YEAR * 360) 
               I_MONTH = MOD(I_DAY/30,12)+1
               I_DAY   = MOD(I_DAY,30)+1
               I_DAY_NUMBER = I_DAY+30*(I_MONTH-1)
            ENDIF
*DECLARE STWORK1A
*/ change field length for long run
*D STWORK1A.266
            CHARACTER*21                          
*D STWORK1A.1143
            PPNAME=STRING(JJ+1:JJ+21)  
*DECLARE CHISTG
*D CHISTG.20
      CHARACTER*21 END_DUMPim(N_INTERNAL_MODEL_MAX)!most recent dumpname
*D GGH3F400.1 
      CHARACTER*21 SAFEDMPim(N_INTERNAL_MODEL_MAX)
*D GGH3F400.3 
      CHARACTER*21 NEWSAFEim(N_INTERNAL_MODEL_MAX)
*D GKR1F404.245
      CHARACTER*21 LASTATMim(N_INTERNAL_MODEL_MAX) ! Keep name of last
*D GKR1F404.248 
      CHARACTER*21 CURRATMim(N_INTERNAL_MODEL_MAX) ! Keep name of
*D GKR1F404.251
      CHARACTER*21 LASTDMPim(N_INTERNAL_MODEL_MAX) ! Keep name of last
*DECLARE CHISTO
*D CHISTO.20,CHISTO.21
      CHARACTER*21 RUN_COMPCODE        ! Run completion code
      CHARACTER*21 RUN_LAST_MEAN       ! Last mean dump created by run
*DECLARE INITCHS1
*D GRB1F305.97,GLW4F403.6
      RUN_COMPCODE='                     '
      RUN_LAST_MEAN='                     '
*D GRB1F305.120
*D GLW4F403.18,GKR1F404.259
      END_DUMPim(i)='                     '
      SAFEDMPim(i)='                     '
      NEWSAFEim(i)='                     '
      LASTATMim(i)='                     '
      CURRATMim(i)='                     '
      LASTDMPim(i)='                     '
c
*DECLARE DUMPCTL1
*D GKR1F404.242
 620  FORMAT('%%% ',A21,' ARCHIVE DUMP')
*D DUMPCTL1.372
 610    FORMAT('%%% ',A21,' DELETE')
*D DUMPCTL1.270
      ARESTART(18:38)=DUMPNAME
*D GDR3F305.27
      AOMEAN(18:38) = DUMPNAME
*D DUMPCTL1.320
      ORESTART(18:38)=DUMPNAME
*D WRB1F401.141
      WRESTART(18:38)=DUMPNAME
c
*DECLARE PRINTHS1
*D WRB1F401.704,WRB1F401.705
     *'Last Restart dump(s) written      : ',A21,3(4X,A21)/
     *'Current Restart dump(s) name      : ',A21,3(4X,A21)/
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  meadlengths.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID meadlengths
*/ work around problems in calculating input and output stash lengths
*/ for use with ocean MEAD diagnostics
*DC ADDRLN1
*D ORH5F403.296
      ELSE IF (IGPL.EQ.41 .OR. IGPL.EQ.43) THEN
*D ORH5F403.300
      ELSE IF (IGPL.EQ.42 .OR. IGPL.EQ.44) THEN
*DC OUTPTL1
*I GPB1F402.516
        ! for data on ocean zonal grid, subdomain returned by
        ! GLOBAL_TO_LOCAL_SUBDOMAIN does not include the halos, but needs to,
        ! so force the right answer
        IF (IGP.eq.43 .or. IGP.eq.44) THEN
          local_north=1
          local_south=g_lasize(2,mype)
        ENDIF
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  medout44.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/ This modset is for version 4.4 of the Unified Model
*/ This is a FAMOUS mod. 
*/ **********************************************************************
*/ Modification set - MEDOUT44
*/ Author - Jonathan Palmer, 9th February 1998
*/ **********************************************************************
*/
*/ UM44 checks for the 3.75*2.5 degree ocean grid when setting
*/ up the parameters for the Mediterranean outflow parametrisation.
*/ The parameters set, however, a re harwired for HadCM2; they are
*/ timestep independent because HadCM2 mixes fully across the Straits
*/ of Gibraltar at every tstep. This modset makes the low res med
*/ outlow look like the hi-res verion ie HadCM3, by making tendfrc
*/ suitably grid and timestep dependent...
*/
*IDENT MEDOUT44
*/
*DECLARE OSETCON
*D OJG1F404.26
             tendfrc=9.6e-5 * DTTS/3600
*D OJG1F404.38
! Make tendfrc suitably dependent on timestep. There is one
! low res grid box to 6 hi-res boxes. But we are only mixing
! one pair here, instead of 2 pairs in the hi-res code above.
! Hence 9.6e-5*2/6 preserves mixing volume per second
             tendfrc=9.6e-5 * 2 / 6 * DTTS/3600
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  model_fix.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID MODELFIX
*/
*/ Delete barrier call in TIMER3A routine, can get deadlock otherwise.
*/
*DECLARE TIMER3A
*D GPB8F405.15
C   Skipping sync
C        CALL GC_GSYNC(nproc,info)
*/
*/ Initialise O_ATMCO2 array correctly
*/
*DECLARE SWAPA2O2
*D CCN1F405.160
*I GRR0F403.229
      DO I=1,CO2_DIMA
      O_ATMCO2(I) = RMDI
      ENDDO ! I
*/
*/ This code will work on any MPP machine, therefore change 
*/ *IF  DEF,MPP,AND,DEF,T3E to *IF  DEF,MPP.
*/
*DECLARE SWPLND1A
*D SWPLND1A.3
*IF  DEF,MPP
*/
*/ If running coupled model then J will be declared twice
*/
*DECLARE INTFCTL1
*I GMB1F405.146
      INTEGER J                   ! Ocean: For setting LBC_UNIT_NO_O
*D UDG1F305.178
*D GMB1F405.154
*/
*/ If using VECTLIB delete non-fortran line
*/
*DECLARE TRSFC3A
*DELETE PXVECTLB.149
*/
*/ On some machines (e.g. hp,ibm) abort doesn't take any arguments
*/ this causes a compilation error. If you declare abort as external
*/ compiler doesn't give an error but abort is still called, without
*/ writing the error message.
*/
*DECLARE FIELDOP1
*D FIELDOP1.105
      External readff,setpos,ioerror,fieldop_main,abort
*/
*DECLARE FIELDCOS
*D FIELDCOS.35
      EXTERNAL READFF,SETPOS,IOERROR,READ_WRITE,ABORT
*/
*DECLARE DUMMYVEG
*I DUMMYVEG.30

      EXTERNAL ABORT
*/
*DECLARE RDLSM1A
*B RDLSM1A.60

      EXTERNAL ABORT
*/
*DECLARE SOLANG1A
*B SOLANG1A.71

      EXTERNAL ABORT
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  no_gcom_setopt.mf77
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID GCOMSETOPT
*/ remove gc_setopt to allow gcom 7.3+
*/ SSW 11/09
*DECLARE FILTER1A
*D GPB2F404.122
*D GPB2F404.228
*DECLARE GATFLD1A
*D GPB2F402.83
*DECLARE GATZON1A
*D GATZON1A.164
*DECLARE GENGAT1A
*D GENGAT1A.112
*DECLARE GENSCT1A
*D GENSCT1A.109
*DECLARE GTALBC1A
*D GPB2F402.343
*D GPB2F402.346
*D GPB2F402.350
*DECLARE GTOLBC1A
*D GTOLBC1A.490
*D GTOLBC1A.604
*DECLARE MUSPAC1A
*D GPB0F403.1789
*D GPB0F403.1797
*DECLARE SCALBC1A
*D GPB2F402.352
*D GPB2F402.355
*D GPB2F402.359
*DECLARE SCOLBC1A
*D SCOLBC1A.494
*D SCOLBC1A.609
*DECLARE SCTFLD1A
*D GPB2F402.252
*DECLARE SCTZON1A
*D SCTZON1A.138
*DECLARE STGFLD1A
*D GPB0F403.3535
*DECLARE STMERM1A
*D GPB0F403.1557
*DECLARE STMERM1A
*D GPB0F403.1567
*D GPB0F403.1607
*DECLARE STSFLD1A
*D STSFLD1A.397
*DECLARE STZONM1A
*D GPB0F403.2831
*D GPB0F403.2841
*D GPB0F403.2881
*DECLARE SWPBND1A
*D GPB2F402.364
*D GPB2F402.368
*D GPB2F402.372
*D GPB2F402.376
*D GPB2F402.380
*D GPB2F402.384
*D GPB2F402.388
*D GPB2F402.392
*D GPB3F403.314
*D GPB3F403.340
*D GPB3F403.363
*D GPB3F403.388
*DECLARE UMSHELL1
*D GBCAF404.51,GBCAF404.60
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  nomsrest
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID NOMSREST
*/
*/  Modset to stop the total mass of the atmosphere being re-set for 
*/    use next time to the value it's just been corrected from !
*/  Just deletes one unthinking line, & alters a comment to match.
*/                                                        WJI 3/9/99
*DECLARE ATMPHY1
*D ATMPHY1.489
C SWAP MODIFIED TOTAL DRY ENERGY OF ATMOSPHERE
*D ATMPHY1.491
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  ofilter_mpp.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID OFILTMPP
*/
*/ Replace shmem load balancing code for ocean filtering with 
*/ gcom version. Based on NEC benchmarking code.
*/
*/ Need to add code for free surface filtering.
*/
*DECLARE ARGOCFIL
*D ORH1F405.54
     &,SLAV_CNT_T,SLAV_CNT_F,MIN_K_U,MIN_K_T,MAX_K_U,MAX_K_T,MYSLAVE,
*/
*DECLARE TYPOCFIL
*I ORH1F403.2
       INTEGER MIN_K_U(JMT,0:O_NPROC-1)
     &,        MIN_K_T(JMT,0:O_NPROC-1)
     &,        MAX_K_U(JMT,0:O_NPROC-1)
     &,        MAX_K_T(JMT,0:O_NPROC-1)
       LOGICAL MYSLAVE(JMT,0:O_NPROC-1)
*/
*DECLARE TYPOCDPT
*I ORH1F405.35
     &, JP_MIKU,JP_MIKT,JP_MAKU,JP_MAKT,JP_MYSL
*/
*DECLARE COMOCDPT
*I ORH1F405.37
     &, JP_MIKU,JP_MIKT,JP_MAKU,JP_MAKT,JP_MYSL
*/
*DECLARE OCNARYPT
*I ORH1F405.30
      JP_MIKU=icount
      icount=icount+JMT*O_NPROC
      JP_MIKT=icount
      icount=icount+JMT*O_NPROC
      JP_MAKU=icount
      icount=icount+JMT*O_NPROC
      JP_MAKT=icount
      icount=icount+JMT*O_NPROC
      JP_MYSL=icount
      icount=icount+JMT*O_NPROC
*/
*DECLARE ARGOCTOP
*I ORH1F405.43
     &,O_SPCON(JP_MIKU),O_SPCON(JP_MIKT),O_SPCON(JP_MAKU)
     &,O_SPCON(JP_MAKT),O_SPCON(JP_MYSL)
*/
*DECLARE CALCFVN
*I ORH6F402.86
     &,O_NPROC
*I PXORDER.7
     &,O_NPROC           ! IN O_NPROC for ocean
*/
*DECLARE BLOKINIT
*I ORH6F402.85
     &,O_NPROC
*/
*DECLARE CTODUMP
*I ORH6F402.83
     & O_NPROC,
*I ORH6F402.84
     &,O_NPROC       ! IN O_NPROC for ocean
*/
*DECLARE OCNFRST1
*I ORH6F402.82
     & O_NPROC,
*/
*DECLARE OFILTR
*I ORH6F401.21
*CALL COCNINDX
*D ORH1F403.147
*/
*DECLARE OSETCON
*D ORH1F405.77
*IF DEF,MPP
*/
*DECLARE DECMFLTR
*I ORH1F405.97
!     03.07.03     J.Cole      Rewrite load balancing to use two-sided
!                              GCOM communications. Based on NEC 
!                              benchmarking code.
*I DECMFLTR.37
*CALL COCNINDX
*D DECMFLTR.41
*I ORH1F405.133

      INTEGER KLAST
      INTEGER IP0
      INTEGER MASTERCOUNT, TOTALWORK
      INTEGER NPROC_AV, ASSIGNEDCOUNT, IPX, ICOUNT

      INTEGER NWORKLOAD(0:O_NPROC-1)
      INTEGER SLAVELIST(0:O_NPROC-1,JMT_GLOBAL)   ! Master -> Slave1 
                                                  !        -> Slave2..
      INTEGER MYMASTER(0:O_NPROC-1,JMT_GLOBAL)    
      INTEGER ISLAVECOUNT(0:O_NPROC-1)

      LOGICAL MASTERLOAD(0:O_NPROC-1,JMT_GLOBAL)  ! TRUE: I am a master
      LOGICAL ASSIGNED(0:O_NPROC-1)

      REAL XMAX, WORKAVERAGE
      REAL SLAVECOUNT(0:O_NPROC-1)
*I ORH1F405.170


      DO J=1,JMT_GLOBAL
         DO PE=0,O_NPROC-1
            MASTERLOAD(PE,J) = .FALSE.
            SLAVELIST(PE,J) = -1
            MYMASTER(PE,J) = -1
         ENDDO
      ENDDO
      DO PE=0,O_NPROC-1
         DO J=1,JMT
            MIN_K_U(J,PE) = -1
            MAX_K_U(J,PE) = -2
            MIN_K_T(J,PE) = -1
            MAX_K_T(J,PE) = -2
            MYSLAVE(J,PE) = .FALSE.
         ENDDO
      ENDDO
*I ORH1F405.183
            NWORKLOAD(IPROC) = 0
*I ORH1F405.219
                           NWORKLOAD(IPROC) = NWORKLOAD(IPROC) + SEG_LEN
*I ORH1F405.257
                           NWORKLOAD(IPROC) = NWORKLOAD(IPROC) + SEG_LEN
*D ORH1F405.277,ORH1F405.459
         MASTERCOUNT=0
         TOTALWORK=0
         DO IPROC=0, O_NPROC - 1   ! For all PEs
            MASTERLOAD(IPROC,JIND) = NWORKLOAD(IPROC) .GT. 0
            ASSIGNED(IPROC)=.FALSE.
            IF (NWORKLOAD(IPROC) .GT. 0)  THEN
               WRITE(6,*)'JIND,IPROC,NWORKLOAD=',
     &                    JIND,IPROC,NWORKLOAD(IPROC)
               MASTERCOUNT=MASTERCOUNT+1
               TOTALWORK=TOTALWORK+NWORKLOAD(IPROC)
            ENDIF
         ENDDO
         
         IF (TOTALWORK .GT. 0) THEN
            ASSIGNEDCOUNT=0
            NPROC_AV=O_NPROC
            WORKAVERAGE=TOTALWORK/REAL(O_NPROC)
            WRITE(6,*) 'TOTALWORK,Workaverage=',TOTALWORK,Workaverage
            DO IPROC=0, O_NPROC - 1 
               SLAVECOUNT(IPROC)=0.0
               ISLAVECOUNT(IPROC)=0
               IF (MASTERLOAD(IPROC,JIND)) THEN
                  IF (NWORKLOAD(IPROC)/WORKAVERAGE .LT. 1.0) THEN
                     SLAVECOUNT(IPROC)=1.0
                     ISLAVECOUNT(IPROC)=1
                     NPROC_AV = NPROC_AV - 1
                     TOTALWORK=TOTALWORK-NWORKLOAD(IPROC)
                     WORKAVERAGE=TOTALWORK/NPROC_AV 
                     WRITE(6,*) 'TOTALWORK,Workaverage=',
     &                           TOTALWORK,Workaverage
                     ASSIGNEDCOUNT=ASSIGNEDCOUNT+1
                  ENDIF
               ENDIF
            ENDDO
            DO IPROC=0, O_NPROC - 1 
               IF (MASTERLOAD(IPROC,JIND) .AND. 
     &             ISLAVECOUNT(IPROC) .EQ. 0) THEN
                  SLAVECOUNT(IPROC)=NWORKLOAD(IPROC)/WORKAVERAGE
                  ISLAVECOUNT(IPROC)=SLAVECOUNT(IPROC)
                  ASSIGNEDCOUNT=ASSIGNEDCOUNT+ISLAVECOUNT(IPROC)
               ENDIF
            ENDDO

            IF (ASSIGNEDCOUNT .GT. O_NPROC) STOP 'ASSIGNEDCOUNT'
            DO ICOUNT=ASSIGNEDCOUNT,O_NPROC-1
               XMAX=-1.
               DO IPROC=0, O_NPROC - 1
                  IF (SLAVECOUNT(IPROC) .GT. 0) THEN
                     IF (SLAVECOUNT(IPROC)-ISLAVECOUNT(IPROC) .GT. 
     &                   XMAX) THEN
                        XMAX=SLAVECOUNT(IPROC)-ISLAVECOUNT(IPROC)
                        IPX=IPROC
                     ENDIF
                  ENDIF
               ENDDO
               IF (XMAX .LT. 0) STOP 'XMAX'
               ISLAVECOUNT(IPX)=ISLAVECOUNT(IPX)+1 
               ASSIGNEDCOUNT=ASSIGNEDCOUNT+1
               SLAVECOUNT(IPX)=0
            ENDDO

            IF (ASSIGNEDCOUNT .NE. O_NPROC) STOP 'ASSIGNEDCOUNT'

            DO IPROC=0, O_NPROC - 1
               IF (ISLAVECOUNT(IPROC) .NE. 0) THEN
                  WRITE(6,*) 'IPROC,ISLAVECOUNT,NWORKLOAD/ISLAVECOUNT=',
     &                        IPROC,ISLAVECOUNT(IPROC),
     &                        NWORKLOAD(IPROC)/ISLAVECOUNT(IPROC)
               ELSE
                  WRITE(6,*) 'IPROC,ISLAVECOUNT=',
     &                        IPROC,ISLAVECOUNT(IPROC)
               ENDIF
            ENDDO

            ICOUNT=0
            DO IPROC=0, O_NPROC - 1
               IF (ISLAVECOUNT(IPROC) .GT. 0) THEN
                  IPX = IPROC
                  ICOUNT = ISLAVECOUNT(IPROC)-1
                  DO IP0=0, O_NPROC - 1
                     IF (ICOUNT .GT. 0) THEN
                        IF ( .NOT. ( MASTERLOAD(IP0,JIND) .OR. 
     &                               ASSIGNED(IP0) )) THEN
                           ASSIGNED(IP0)=.TRUE.
                           SLAVELIST(IPX,JIND)=IP0
                           MYMASTER(IP0,JIND)=IPROC
                           IPX=IP0
                           ICOUNT=ICOUNT-1
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO

            DO IPROC=0, O_NPROC - 1
               IF (MASTERLOAD(IPROC,JIND)) THEN
                  WRITE(6,*) 'Master=',IPROC
                  IPX=IPROC
                  DO IP0=2,ISLAVECOUNT(IPROC)
                     WRITE(6,*) 'Slave=',SLAVELIST(IPX,JIND)
                     IPX=SLAVELIST(IPX,JIND)
                     IF (IPROC .EQ. O_MYPE) MYSLAVE(JIND,IPX) = .TRUE.
                  ENDDO
               ENDIF
            ENDDO

!           We now have decided who works for whom.
!           Now we go over all masterprocessors and assign adjacent 
!           work to the slaves which is about size of the average.

            WRITE(6,*)'JIND,SEG_CNT_U = ',JIND,SEG_CNT_U
            WRITE(6,*)'JIND,SEG_CNT_T = ',JIND,SEG_CNT_T

            ! For each PE
            DO PE = 0, O_NPROC-1

               IF (MASTERLOAD(PE,JIND)) THEN          

                  MAX_WORK = 0
                  TOT_WORK_LEFT = NWORKLOAD(PE)
                  WORKAVERAGE = NWORKLOAD(PE)/ISLAVECOUNT(PE)
                  WRITE(6,*)'JIND,PE,NWORKLOAD,WORKAVERAGE = ',
     &                       JIND,PE,NWORKLOAD(PE),WORKAVERAGE
                  IPX = SLAVELIST(PE,JIND)
                  ICOUNT = ISLAVECOUNT(PE)
                  IF (ICOUNT .EQ. 1) IPX=PE
                  KLAST = -1

                  ! For each U/V segment
                  DO L = 1, SEG_CNT_U

                     IF (WORK_PE_U(L) .EQ. PE) THEN
                        ! If this is my segment

                        IF ( WORK_K_U(L) .NE. KLAST ) THEN
!                          Change slave only when K changes
                           KLAST=WORK_K_U(L)
                           IF ( MAX_WORK .GE. WORKAVERAGE ) THEN
                              WRITE(6,*)'IPX,WORKAVERAGE,MAX_WORK = ',
     &                                   IPX,WORKAVERAGE,MAX_WORK
                              IPX = SLAVELIST(IPX,JIND)
                              ICOUNT = ICOUNT-1
                              IF ( ICOUNT .LE. 1) IPX=PE
                              TOT_WORK_LEFT = TOT_WORK_LEFT - MAX_WORK
                              WORKAVERAGE = TOT_WORK_LEFT / 
     &                                      MAX(1,ICOUNT)
                              MAX_WORK=0
                           ENDIF 
                        ENDIF 

                        MAX_WORK = MAX_WORK + WORK_SIZE_U(L)

                        ! Assign to current slave
                        IF (IPX .EQ. O_MYPE) THEN
                           MAST_CNT_U(JIND) = MAST_CNT_U(JIND) + 1
                           MAST_PE_U(MAST_CNT_U(JIND),JIND)  = PE
                           MAST_K_U(MAST_CNT_U(JIND),JIND) = WORK_K_U(L)
                           MAST_SEG_U(MAST_CNT_U(JIND),JIND) = 
     &                        WORK_SEG_U(L)
                        ENDIF
                        IF (MIN_K_U(JIND,IPX) .LT. 0) 
     &                     MIN_K_U(JIND,IPX)=WORK_K_U(L)
                        MAX_K_U(JIND,IPX)=WORK_K_U(L)

                        IF (PE .EQ. O_MYPE) THEN
                           SLAV_CNT_U(JIND) = SLAV_CNT_U(JIND) + 1
                        ENDIF

                     ENDIF

                  ENDDO  ! Over L

                  KLAST=-1
                  DO L = 1, SEG_CNT_T

                     IF (WORK_PE_T(L) .EQ. PE) THEN
                        ! If this is my segment

                        IF ( WORK_K_T(L) .NE. KLAST ) THEN
!                          Change slave only when K changes
                           KLAST=WORK_K_T(L)
                           IF ( MAX_WORK .GE. WORKAVERAGE ) THEN
                              WRITE(6,*)'IPX,WORKAVERAGE,MAX_WORK = ',
     &                                   IPX,WORKAVERAGE,MAX_WORK
                              IPX = SLAVELIST(IPX,JIND)
                              ICOUNT = ICOUNT -1
                              IF ( ICOUNT .LE. 1) IPX=PE
                              TOT_WORK_LEFT = TOT_WORK_LEFT - MAX_WORK
                              WORKAVERAGE = TOT_WORK_LEFT / 
     &                                      MAX(1,ICOUNT)
                              MAX_WORK=0
                           ENDIF
                        ENDIF

                        ! Assign to current slave
                        MAX_WORK = MAX_WORK + WORK_SIZE_T(L)

                        IF (IPX .EQ. O_MYPE) THEN
                           MAST_CNT_T(JIND) = MAST_CNT_T(JIND) + 1
                           MAST_PE_T(MAST_CNT_T(JIND),JIND)  = PE
                           MAST_K_T(MAST_CNT_T(JIND),JIND) = WORK_K_T(L)
                           MAST_SEG_T(MAST_CNT_T(JIND),JIND) = 
     &                        WORK_SEG_T(L)
                        ENDIF
                        IF (MIN_K_T(JIND,IPX) .LT. 0) 
     &                     MIN_K_T(JIND,IPX)=WORK_K_T(L)
                        MAX_K_T(JIND,IPX)=WORK_K_T(L)

                        IF (PE .EQ. O_MYPE) THEN
                           SLAV_CNT_T(JIND) = SLAV_CNT_T(JIND) + 1
                        ENDIF 

                     ENDIF

                  ENDDO  ! Over L
                  WRITE(6,*)'IPX,WORKAVERAGE,MAX_WORK = ',
     &                       IPX,WORKAVERAGE,MAX_WORK

               ENDIF
            ENDDO  ! Over PE

            WRITE(6,*) 'JIND,SLAV_CNT_U,MAST_CNT_U=',
     &                  JIND,SLAV_CNT_U(JIND),MAST_CNT_U(JIND)
            WRITE(6,*) 'JIND,SLAV_CNT_T,MAST_CNT_T=',
     &                  JIND,SLAV_CNT_T(JIND),MAST_CNT_T(JIND)
            WRITE(6,*)
         ENDIF
*/
*DECLARE OFLTCN2A
*D OFLTCN2A.3
*IF -DEF,MPP
*/
*DECLARE OFLTCNTL
*D ORH1F405.522
*IF DEF,MPP
*I OFLTCNTL.45
!     03.07.03   J.Cole      Rewrite load balancing to use two-sided
!                            GCOM communications. Based on NEC 
!                            benchmarking code.
*D ORH1F405.527
*D OFLTCNTL.82,OFLTCNTL.85
       REAL FTARR(IMTIMT_FLT)
*I ORH1F405.535
     &,        ISTAT
*I ORH1F405.539
     &,     FIELD_TEMP(IMT*2,KM*2)
     &,     FIELD_FILT(IMT*2)
*D OFLTCNTL.87,OFLTCNTL.106
*IF DEF,DEBUG
       WRITE(6,*)'MAST_CNT_U(',J,') = ',MAST_CNT_U(J)
       WRITE(6,*)'SLAV_CNT_U(',J,') = ',SLAV_CNT_U(J)
       WRITE(6,*)'MAST_CNT_T(',J,') = ',MAST_CNT_T(J)
       WRITE(6,*)'SLAV_CNT_T(',J,') = ',SLAV_CNT_T(J)
*ENDIF
*D ORH1F405.562,ORH1F405.563
            ! I must set up my UA and VA fields in a temporary array
            ! ready to be send to the slaves
*D ORH1F405.570,ORH1F405.571
            ! Send U,V data to be filtered to the slave processors
            DO IPROC=0,NPROC-1
               IF (MYSLAVE(J,IPROC)) THEN
                  K=MIN_K_U(J,IPROC)
                  SIZEB=(MAX_K_U(J,IPROC)-K+1)*IMT*2
                  IF (SIZEB .GT. 0) THEN 
*IF DEF,DEBUG
                     WRITE(6,*)'SENDING U TO ',IPROC,SIZEB,
     &                                         K,MAX_K_U(J,IPROC)
*ENDIF
                     CALL GC_RSEND(100,SIZEB,IPROC,ISTAT,
     &                             FIELD_TEMP(1,K),FIELD_TEMP(1,K))
                     CALL GC_RSEND(101,1,IPROC,ISTAT,CS_FILT,CS(J))
                     CALL GC_RSEND(102,1,IPROC,ISTAT,PHI_FILT,PHI(J))
                  ENDIF
               ENDIF
            ENDDO
*I OFLTCNTL.190
         ! If I'm a slave for any processor other than myself then
         ! receive U,V data
         IF (MAST_CNT_U(J).GT.0 .AND. SLAV_CNT_U(J).EQ.0) THEN
            IPROC = MAST_PE_U(1,J)
            K=MIN_K_U(J,MYPE)
            SIZEB=(MAX_K_U(J,MYPE)-K+1)*IMT*2
*IF DEF,DEBUG
            WRITE(6,*)'RECEIVING U FROM ',IPROC,SIZEB,K,MAX_K_U(J,MYPE)
*ENDIF
            CALL GC_RRECV(100,SIZEB,IPROC,ISTAT,
     &                    FIELD_TEMP(1,K),FIELD_TEMP(1,K))
            CALL GC_RRECV(101,1,IPROC,ISTAT,CS_FILT,CS(J))
            CALL GC_RRECV(102,1,IPROC,ISTAT,PHI_FILT,PHI(J))
         ENDIF

         ! If I'm a slave for myself then I already have U,V data
         IF (MAST_CNT_U(J).GT.0 .AND. SLAV_CNT_U(J).GT.0) THEN
            CS_FILT = CS(J)
            PHI_FILT = PHI(J)
         ENDIF
*D ORH1F405.577
            ! temporary array ready to be send to the slaves.
*D ORH1F405.587,ORH1F405.588
            ! Send T,S data to be filtered to the slave processors
            DO IPROC=0,NPROC-1
               IF (MYSLAVE(J,IPROC)) THEN
                  K=MIN_K_T(J,IPROC)
                  SIZEB=(MAX_K_T(J,IPROC)-K+1)*IMT*2
                  IF (SIZEB .GT. 0) THEN 
*IF DEF,DEBUG
                     WRITE(6,*)'SENDING T TO ',IPROC,SIZEB,
     &                                         K,MAX_K_T(J,IPROC)
*ENDIF
                     CALL GC_RSEND(103,SIZEB,IPROC,ISTAT,
     &                            FIELD_TEMP(1,K+KM),FIELD_TEMP(1,K+KM))
                     CALL GC_RSEND(104,1,IPROC,ISTAT,CST_FILT,CST(J))
                     CALL GC_ISEND(105,1,IPROC,ISTAT,KMT1_FILT,KMT(1))
                  ENDIF
               ENDIF
            ENDDO
*I OFLTCNTL.209
         ! If I'm a slave for any processor other than myself then
         ! receive T,S data
         IF (MAST_CNT_T(J).GT.0 .AND. SLAV_CNT_T(J).EQ.0) THEN
            IPROC = MAST_PE_T(1,J)
            K=MIN_K_T(J,MYPE)
            SIZEB=(MAX_K_T(J,MYPE)-K+1)*IMT*2
*IF DEF,DEBUG
            WRITE(6,*)'RECEIVING T FROM ',IPROC,SIZEB,K,MAX_K_T(J,MYPE)
*ENDIF
            CALL GC_RRECV(103,SIZEB,IPROC,ISTAT,
     &                    FIELD_TEMP(1,K+KM),FIELD_TEMP(1,K+KM))
            CALL GC_RRECV(104,1,IPROC,ISTAT,CST_FILT,CST(J))
            CALL GC_IRECV(105,1,IPROC,ISTAT,KMT1_FILT,KMT(1))
         ENDIF

         ! If I'm a slave for myself then I already have T,S data
         IF (MAST_CNT_T(J).GT.0 .AND. SLAV_CNT_T(J).GT.0) THEN
            CST_FILT = CST(J)
            KMT1_FILT = KMT(J)
         ENDIF

*D ORH1F405.590
         CALL GC_GSYNC(O_NPROC,ISTAT)
*D ORH1F405.609,ORH1F405.613
*D ORH1F405.649,ORH1F405.653
*D ORH1F405.656,ORH1F405.660
*D ORH1F405.662,ORH1F405.676
           DO I=IS,IEA
               UDIF(I-ISM1,K)=-((FX*FIELD_TEMP(I,K))*SPSIN(I))-
     &                           FIELD_TEMP(I+IMT,K)*SPCOS(I)
               VDIF(I-ISM1,K)= ((FX*FIELD_TEMP(I,K))*SPCOS(I))-
     &                           FIELD_TEMP(I+IMT,K)*SPSIN(I)
           ENDDO

           IF (IE.GE.IMU)THEN
               DO I=2,IEB
                 UDIF(I+II,K)=((-FX*FIELD_TEMP(I,K))*SPSIN(I))-
     &                           FIELD_TEMP(I+IMT,K)*SPCOS(I)
                 VDIF(I+II,K)=(( FX*FIELD_TEMP(I,K))*SPCOS(I))-
     &                           FIELD_TEMP(I+IMT,K)*SPSIN(I)
               ENDDO
           ENDIF
*D ORH1F405.681,ORH1F405.688
           DO 720 I=IS,IEA
                    FIELD_TEMP(I,K)=FX*(-UDIF(I-ISM1 ,K)*SPSIN(I)
     &                  +VDIF(I-ISM1,K)*SPCOS(I))
                    FIELD_TEMP(I+IMT,K)=-UDIF(I-ISM1 ,K)*SPCOS(I)
     &                  -VDIF(I-ISM1,K)*SPSIN(I)
  720      CONTINUE

           IF(IE.GE.IMU) THEN
                    DO 722 I=2,IEB
                       FIELD_TEMP(I,K)=FX*(-UDIF(I+II,K)*SPSIN(I)
     &                    +VDIF(I+II,K)*SPCOS(I))
                       FIELD_TEMP(I+IMT,K)=-UDIF(I+II,K)*SPCOS(I)
     &                    -VDIF(I+II,K)*SPSIN(I)
  722               CONTINUE
           ENDIF
*D OFLTCNTL.324,ORH1F405.697
*D ORH1F405.716,OFLTCNTL.379
*I ORH1F405.734
           ISM1=IS-1
*D ORH1F405.742,ORH1F405.747
*D ORH1F405.750
*D ORH1F405.752,ORH1F405.753
              DO I=IS,IEA
                 FIELD_FILT(((MM-1-TLAST)*IMT)+I-ISM1) = 
     &           FIELD_TEMP(((MM-1-TLAST)*IMT)+I,K+KM)
              ENDDO
              IF (IE.GE.IMT) THEN
                 DO I=2,IEB
                    FIELD_FILT(((MM-1-TLAST)*IMT)+I+II) = 
     &              FIELD_TEMP(((MM-1-TLAST)*IMT)+I,K+KM)
                 ENDDO
              ENDIF
*D ORH1F405.758
     &           FTARR,FIELD_FILT(((MM-1-TLAST)*IMT)+1),IM,M,N,IDX)
*D ORH1F405.760,ORH1F405.768
              DO I=IS,IEA
                 FIELD_TEMP(((MM-1-TLAST)*IMT)+I,K+KM) = 
     &           FIELD_FILT(((MM-1-TLAST)*IMT)+I-ISM1)
              ENDDO
              IF (IE.GE.IMT) THEN
                 DO I=2,IEB
                    FIELD_TEMP(((MM-1-TLAST)*IMT)+I,K+KM) = 
     &              FIELD_FILT(((MM-1-TLAST)*IMT)+I+II)
                 ENDDO
              ENDIF
*D OFLTCNTL.541,OFLTCNTL.547
        ! Send filtered data from slaves to master processor
        IF (MAST_CNT_U(J) .GT. 0) THEN
           IPROC = MAST_PE_U(1,J)
           K=MAST_K_U(1,J)
           SIZEB=(MAST_K_U(MAST_CNT_U(J),J)-K+1)*IMT*2
           IF (IPROC .NE. MYPE) THEN
*IF DEF,DEBUG
              WRITE(6,*)'SENDING U TO ',IPROC,SIZEB,
     &                   K,MAST_K_U(MAST_CNT_U(J),J)
*ENDIF
              CALL GC_RSEND(106,SIZEB,IPROC,ISTAT,
     &                      FIELD_TEMP(1,K),FIELD_TEMP(1,K))
           ENDIF
        ENDIF
        IF (MAST_CNT_T(J) .GT. 0) THEN
           IPROC = MAST_PE_T(1,J)
           K=MAST_K_T(1,J)
           SIZEB=(MAST_K_T(MAST_CNT_T(J),J)-K+1)*IMT*2
           IF ( IPROC .NE. MYPE) THEN
*IF DEF,DEBUG
              WRITE(6,*)'SENDING T TO ',IPROC,SIZEB,
     &                   K,MAST_K_T(MAST_CNT_T(J),J)
*ENDIF
              CALL GC_RSEND(107,SIZEB,IPROC,ISTAT,
     &                      FIELD_TEMP(1,K+KM),FIELD_TEMP(1,K+KM))
           ENDIF
        ENDIF

        ! If master processor get filtered data from slaves
        DO IPROC=0,NPROC-1
           IF (MYSLAVE(J,IPROC) ) THEN
              K=MIN_K_U(J,IPROC)
              SIZEB=(MAX_K_U(J,IPROC)-K+1)*IMT*2
              IF (SIZEB .GT. 0) THEN
*IF DEF,DEBUG
                 WRITE(6,*)'RECEIVING U FROM ',IPROC,SIZEB,
     &                      K,MAX_K_U(J,IPROC)
*ENDIF
                 CALL GC_RRECV(106,SIZEB,IPROC,ISTAT,
     &                         FIELD_TEMP(1,K),FIELD_TEMP(1,K))
              ENDIF
           ENDIF
        ENDDO
        DO IPROC=0,NPROC-1
           IF ( MYSLAVE(J,IPROC) ) THEN
              K=MIN_K_T(J,IPROC)
              SIZEB=(MAX_K_T(J,IPROC)-K+1)*IMT*2
              IF ( SIZEB .GT. 0) THEN
*IF DEF,DEBUG
                 WRITE(6,*)'RECEIVING T FROM ',IPROC,SIZEB,
     &                      K,MAX_K_T(J,IPROC)
*ENDIF
                 CALL GC_RRECV(107,SIZEB,IPROC,ISTAT,
     &                         FIELD_TEMP(1,K+KM),FIELD_TEMP(1,K+KM))
              ENDIF
           ENDIF
        ENDDO

        ! Synchronize before further processing.
        CALL GC_GSYNC(O_NPROC,ISTAT)
*D ORH1F405.790
         CALL GC_GSYNC(O_NPROC,ISTAT)
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  orbital_parameters-6.1.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*IDENT CALC_SOLAR
*/
*/  calculate secular variations in orbital parameters for 
*/  paleo runs based on the routines of Berger78 (JAS 35).
*/  basically a backport of SOLPOS and ORBPRM from UM6.1
*/
*COMDECK CSOLAR
c
c common block to hold the orbital parameters so that they 
c can be read in, and then only have to be redone once a year
c
      REAL    SC,
     &        GAMMA,E,TAU0,SINOBL,
     &        E1, E2, E3, E4, DINY,
     &        SEC_VAR_FACTOR
      LOGICAL L_SEC_VAR,L_SEC_VAR_ONLINE,L_SEC_VAR_FILE
      CHARACTER*80 SEC_VAR_FILE
      INTEGER SEC_VAR_YEAR

      COMMON /COMSOLAR/ SC,
     &                  GAMMA,E,TAU0,SINOBL,
     &                  E1, E2, E3, E4, DINY,
     &                  L_SEC_VAR,L_SEC_VAR_ONLINE,
     &                  L_SEC_VAR_FILE,SEC_VAR_FILE,
     &                  SEC_VAR_YEAR, SEC_VAR_FACTOR

      NAMELIST /NLSTSOLAR/ SC,SEC_VAR_YEAR,
     &               L_SEC_VAR,L_SEC_VAR_ONLINE,L_SEC_VAR_FILE,
     &               SEC_VAR_FILE,SEC_VAR_FACTOR
*/
*DECLARE SWSC
*/  replace declaration of SC
*D SWSC.2,SWSC.3
*CALL CSOLAR
*/
*DECLARE READLSA1
*/ read in namelist values for orbital params, set up param. values
*/ for the start of the run
*/
*B READLSA1.34
*CALL CSOLAR
      REAL OBLQ                            ! Obliquity (local)
*B READLSA1.145
c set defaults for solar param namelist
      SC=1365.
      SEC_VAR_YEAR=2000
      L_SEC_VAR=.FALSE.
      L_SEC_VAR_ONLINE=.FALSE.
      L_SEC_VAR_FILE=.FALSE.
      SEC_VAR_FILE="dummy_filename"
      SEC_VAR_FACTOR=1.
c read in any user specified values
      READ(5,NLSTSOLAR)

c if using secular variations, take account of the acceleration factor
      if (L_SEC_VAR) SEC_VAR_YEAR=
     &   ( A_FIXHD(28)-MODEL_BASIS_TIME(1) )*SEC_VAR_FACTOR  + 
     &     MODEL_BASIS_TIME(1)

c initialise orbital params, whether fixed or varying
      CALL ORBPRM(L_SEC_VAR_ONLINE,L_SEC_VAR_FILE,
     &            SEC_VAR_FILE,SEC_VAR_YEAR,LCAL360,
     &            E,GAMMA,OBLQ,TAU0,DINY)

      SINOBL=SIN(OBLQ)

      write(6,*)"solar constants: year,varying?",SEC_VAR_YEAR,L_SEC_VAR
      write(6,*)"               : table of variations?",L_SEC_VAR_FILE
      write(6,*)"               : online variations?",L_SEC_VAR_ONLINE
      write(6,'(a,f8.1,f10.6,f10.6,f10.5,f10.6)')
     &          "                : values",SC,GAMMA,E,TAU0,SINOBL

      E1 = E * (2.-.25*E*E)
      E2 = 1.25 * E*E                ! Coefficients for 3.1.2        
      E3 = E*E*E * 13./12.
      E4=( (1.+E*E*.5)/(1.-E*E) )**2 ! Constant for 3.1.4         
*/
*DECLARE SOLPOS1A
*/ get rid of fixed parameter versions of gamma, e, tau and sinobl
*/ if we want varying parameters and it's a new year, re-call the
*/ orbprm routine. Use the 6.1 calculation to get new SINDEC/SCS
*/ 
*D GSS1F304.678
      SUBROUTINE SOLPOS (DAY, YEAR, SINDEC, SCS, LCAL360_IN)
*D GSS1F304.681        
      LOGICAL LCAL360_IN    !In, true if 360 day calendar in use.
*D SOLPOS1A.43,SOLPOS1A.79
      REAL M, V                            ! Mean & true anomaly
      REAL OBLQ                            ! Obliquity (local)
*CALL C_PI                                                    
      REAL TWOPI
      PARAMETER ( TWOPI = 2. * PI ) 
*CALL CSUBMODL
*CALL CMAXSIZE
*CALL CHSUNITS
*CALL CNTLALL
*CALL CTIME
*CALL CSOLAR

C============================================================
C
c get new versions of the orbital params if needed
c

      IF (L_SEC_VAR .AND. PREVIOUS_TIME(1).lt.I_YEAR) THEN

c work out "orbital year" using the acceleration factor
        SEC_VAR_YEAR= (I_YEAR-MODEL_BASIS_TIME(1))*SEC_VAR_FACTOR  + 
     &            MODEL_BASIS_TIME(1)

        CALL ORBPRM(L_SEC_VAR_ONLINE,L_SEC_VAR_FILE,
     &              SEC_VAR_FILE,SEC_VAR_YEAR,LCAL360,
     &              E,GAMMA,OBLQ,TAU0,DINY)

        SINOBL=SIN(OBLQ)

        write(6,'(a,i8,f10.6,f10.6,f10.5,f10.6)')
     &       "new solar constants for",SEC_VAR_YEAR, GAMMA,E,TAU0,SINOBL

        E1 = E * (2.-.25*E*E)
        E2 = 1.25 * E*E                ! Coefficients for 3.1.2        
        E3 = E*E*E * 13./12.
        E4=( (1.+E*E*.5)/(1.-E*E) )**2 ! Constant for 3.1.4         


      END IF

c     Calculate the mean anomaly at 12Z on the current day.
c     The 0.5 accounts for the time in days since mid-night.            
c     The references are to Smart 1944 (and UMDP23)                     
c     Eq 67 p. 113 and n=2pi/orbital period     (Eq 3.1.1)        
      IF (LCAL360) THEN
        M = (TWOPI / DINY)           * (FLOAT(DAY) - TAU0 - .5)
      ELSE
        M = (TWOPI / 365.2424)       * (FLOAT(DAY) - TAU0 - .5)
      END IF

c       True anomaly, equation 87 in Smart on p. 120 (UMDP23 Eq 3.1.2)  
      V  = M + E1*SIN(M) + E2*SIN(2.*M) + E3*SIN(3.*M)

c       Solar constant scaling factor (UMDP23 Eq 3.1.4)                 
      SCS = E4 * ( 1. + E * COS(V) ) **2

c       sin(solar declination) (UMDP23 Eq 3.1.5)                        
c       The solar declination is related to                             
c        the true longitude of the earth (lambda) by:                   
c        sindec = sin(obliquity) * sin(lambda)                          
c       Lambda is counted counterclockwise from the vernal equinox      
c        and is related to v (the true anomaly) through                 
c        lambda = v + (longitude of perihelion)                         

      SINDEC = SINOBL * SIN (V - GAMMA)
C
*DECK ORBPRM
c                                                                         
c+ Subroutine to calculate the parameters of the Earth's orbit.           
c                                                                         
c (slightly adapted to fit in FAMOUS UM4.5 solpos by r.s.smith 
c  can get E,OBLQ and LPH from precalculated table if req'd
c  NB: mod(TAU0,DINY), mod(GAMMA,TWOPI) loop added at end
c  to match Julia/Michel's offline lookup numbers. I think is
c  OK. The date of the vernal equinox is fixed for the use of
c  LCAL360 - it's more useful in paleoruns to have this as a
c  fixed point rather than use the "real" Gregorian DOY for it)
c
c  Purpose:                                                               
c  This routine returns the parameters of the Earth's orbit               
c   (the eccentricity, obliquity, supplement of the longitude of          
c    perihelion) and the time of the perihelion passage in days
c                                                                          
c  Method:                                                                 
c  For long runs there may be an interest in running with secular          
c   variations in the astronomy. The orbital constants have                
c   been derived from A. L. Berger 1978, J. Atm. Sci, Volume 35            
c   2362-2367. A copy of which can be found in the Met Office library.     
c  For short current runs, or long control runs it is preferrable          
c   not to allow the astronomy to vary, so fixed values are used.          
c                                                                          
c Current Owner of Code: J. M. Edwards                                     
c                                                                          
c History:                                                                 
c       Version         Date                    Comment                    
c       5.2             15/11/00                Original Code              
c                                               E. Ostrom                  
c                                                                          
c Description of Code:                                                     
c   FORTRAN90 complying with UMDP3 Version 7.2 from 5/2/98                 
c                                                                          
c- ---------------------------------------------------------------------   

      SUBROUTINE ORBPRM(L_SEC_VAR_ONLINE, L_SEC_VAR_FILE
     &                , SEC_VAR_FILE, YEAR, LCAL360
     &                , E, GAMMA, OBLQ, TAU0, DINY)                       
                                                                           
      IMPLICIT NONE                                                        
                                                                           
      INTEGER  YEAR       ! Calendar year                   
      LOGICAL  L_SEC_VAR_ONLINE ! Use a table of precalulated values
      LOGICAL  L_SEC_VAR_FILE   ! Calculate orbital values online

      LOGICAL  LCAL360    ! Use a calendar of 360 days      

      CHARACTER*80 SEC_VAR_FILE ! Location of orbital lookup table
                                                                           
c     Parameters of the Earth's orbit:                                     
c                                                                          
      REAL  E             ! Eccentricity of the orbit       
      REAL  GAMMA         ! Supplement of the longitude     
c                                        !  of the perihelion              
      REAL  OBLQ          ! Obliquity of the orbit          
      REAL  TAU0          ! Time of the perihelion          
c                                        !  passage in days                
      REAL  DINY          ! Length of the calendar year     
c                                        !  (in whole days)                
                                                                           
c     Local Variables for use within ORBPRM                                
c                                                                          
      REAL  YEAR_OFFSET                ! Offset of the year from the     
c                                        !  reference year when default    
c                                        !  values apply                   
      REAL  ECN_SN                     ! Eccentricity multiplied by      
c                                        !  the sine of the longitude      
c                                        !  of the perihelion              
      REAL  ECN_CN                     ! Eccentricity multiplied by      
c                                        !  the cosine of the longitude    
c                                        !  of the perihelion              
      REAL  LPH_FIXED_VE               ! Longitude of the perihelion     
c                                        !  relative to a fixed vernal     
c                                        !  equinox                        
      REAL  GN_PRCS                    ! General precession              
      REAL  DATE_VE                    ! Date of the vernal equinox      
c                                        !  in days into the year          
      REAL  NO_LEAP_DAYS               ! The number of leap days,        
c                                        !  used to calculate DATE_VE      
      REAL  MEAN_ANOM_VE               ! Mean anomaly at the vernal      
c                                        !  equinox                        
                                                                           
c     Synthetic constants                                                  
c                                                                          
      REAL  BETA                                                         
      REAL  EE1                                                          
      REAL  EE2                                                          
      REAL  EE3                                                          
                                                                           
      INTEGER  I                       ! Loop variable                   
                                                                           
c     Mathematical constants:                                              
*CALL C_PI
      REAL TWOPI 
      PARAMETER ( TWOPI = 2. * PI )
                                                                           
      REAL TropYearLength
      PARAMETER ( TropYearLength=365.2424)

c     Astronomical Parameters:                                             
c     Default values of the orbital elements
c      (currently those for the epoch J2000 which is 1.5d Jan. 2000):
c     The Eccentricity and Longitue of perhelion are recommended by NAS
c      see (http://ssd.jpl.nasa.gov/elem_planets.html)                     
c     The Obliquity value comes from the Astronomical Almanac for 1984
c      page S26 and is used on several webpages e.g.
c      nedwww.ipac.caltech.edu/help/calc_doc.txt
c      www.stargazing.net/kepler/astrovba2.html
c      http://edhs1.gsfc.nasa.gov/waisdata/docsw/txt/tp4450505.txt         
c
c     The data in the series expansions are adapted from Berger 1978.
c      Andre' L. Berger Journ. of Atm. Sci. Volume 35 p. 2362-2367,
c      and is available from the Met Office library.
c
c     ! Eccentricity of the orbit
      Real E_DFLT
      Parameter ( E_DFLT         = 1.6710222E-02)
c
c     ! Longitude of the perihelion in radians
      Real LPH_DFLT
      Parameter ( LPH_DFLT       = 102.94719*PI/180.0)
c
c     ! Obliquity of the orbit - corresponds to 23.43929111 degrees
      Real OBLQ_DFLT
      Parameter ( OBLQ_DFLT      = 0.409092804)
c
c     ! Reference year for setting the date of the vernal equinox
      Integer YEAR_REF_VE
      Parameter ( YEAR_REF_VE    = 2000)
c
c     ! Date of the vernal equinox in days after the start of the year
c     !  This date is for the year 2000.
      Real DATE_VE_DFLT
      Parameter ( DATE_VE_DFLT   = 79.3159)
c
c     The final parameter required is the time of the perihelion
c     passage, TAU0. For a pure Keplerian orbit, with a specified
c     eccentricity and longitude of the perihelion, this can be
c     deduced from the date of the vernal equinox (as is specified
c     in AMIP-2, for example). In practice it is somewhat more
c     complicated to calculate the time of the perihelion.
c     For simplicity, a mean value for the years 1995-2005 is used
c     here: note that the range of TAU0 in this period is from      
c     1.0 to 3.75 and that there is no simple relationship with leap
c     years.
c                                             
c     ! Time of the perihelion passage in days     
      Real TAU0_DFLT
      Parameter ( TAU0_DFLT      = 2.667)
c                                                                       
c     ------------------------------------------------------------------
c                                                                       
c     The parameters used to calculate secular variations of the orbital
c     elements are taken from A. L. Berger (1978), J. Atm. Sci., vol 35,
c     p. 2362.These have been converted so that:
c     amplitudes are in radians,                  
c     angular frequencies in radians per year and
c     phases in radians                                      
c     Phases have also been converted to be taken relative to
c     J2000 for consistency with the time for the default values of the
c     orbital parameters.
c     Berger's numbers (with time correction) differ slightly from the
c     default values above.
c
c     The obliquity and longitude of the perihelion have been adjusted
c     to agree with the NASA values for J2000, but adjustment of the
c     eccentricity to NASA values is not so easy and has not been done.
c
c
c     ! Reference year     YEAR_REF
      Integer  YEAR_REF
      Parameter  ( YEAR_REF       = 2000)
c
c   -----------------------------------------------------------------
c     Obliquity (Table 1) from the first 24 terms:
c     (enough for deviations of less than 0.002 degrees)
c
c     ! Constant term in the obliquity: from the Astr. Almanac for 1984
c     !  The following value corresponds to 23.320870 degrees at J2000
      Real OBLQ_CNST
      Parameter ( OBLQ_CNST      = 0.40702597)
c
c     ! Number of terms retained in the series for the obliquity
      Integer N_TERM_OBQ 
      Parameter (  N_TERM_OBQ     = 24)
c
c     ! Amplitude
      Real  A(N_TERM_OBQ)
c     ! Angular frequency                                                  
      Real  F(N_TERM_OBQ)                                                
c     ! Phase in the series                                                
      Real  D(N_TERM_OBQ)                                                
c                                                                          
c   -----------------------------------------------------------------      
c     Eccentricity and longitude of the fixed perihelion (Table 4):        
c                                                                          
c     ! Number of terms retained in the series for the                     
c     !  eccentricty and longitude of the perihelion                       
      Integer N_TERM_ECN_LPH
      Parameter  ( N_TERM_ECN_LPH = 19  )                          
c                                                                          
c     ! Amplitude                                                          
      Real  M(N_TERM_ECN_LPH)                                            
c     ! Angular frequency                                                  
      Real  G(N_TERM_ECN_LPH)                                            
c     ! Phase in the series                                                
      Real  B(N_TERM_ECN_LPH)                                            

c ---------------------------------------------------------------
c   some variables for using a lookup table
      INTEGER SEC_VARUNIT,ICODE
      CHARACTER DUMMY
      REAL TABLEYEAR,TABLEYEAR_REF
      REAL FILEYEAR1,FILEYEAR2,FILEYEAR3
      REAL FE1,FOBLQ1,FLPH1
      REAL FE2,FOBLQ2,FLPH2,FLPH3
      REAL LPH_MOVING_VE
      REAL DY1,DY2
c                                                                          
c   ------------------------------------------------------------------72   
c     General Precession (Table 5):                                        
c                                                                          
c     ! Linear rate of precession!                                         
c     ! The value corresponds to 50.439273 seconds per year -Berger 1979   
      Real LIN_RATE_GN_PRCS
      Parameter ( LIN_RATE_GN_PRCS  = 2.44536496E-04  )
c                                                                          
c     ! Constant offset to general precession (in seconds pre year),       
c     ! corrected for 50 years difference in reference time.               
      Real GN_PRCS_CNST
      Parameter ( GN_PRCS_CNST      = 7.14372244E-02  )
c                                                                          
c     ! Number of terms kept in the series for the general precession      
      Integer N_TERM_GN_PRCS
      Parameter ( N_TERM_GN_PRCS = 10  )
c                                                                          
c     ! Amplitude                                                          
      Real  C(N_TERM_GN_PRCS)                                            
c     ! Angular frequency                                                  
      Real  H(N_TERM_GN_PRCS)                                            
c     ! Phase in the series                                                
      Real  R(N_TERM_GN_PRCS)                                            
c                                                                          
c   -----------------------------------------------------------------      
c    Table 1                                                               
c                                                                          
      DATA A/                                                           &  
     &    -1.19372E-02, -4.15640E-03, -3.05103E-03, -2.00849E-03        &  
     &  , -1.51146E-03,  1.49778E-03, -7.88065E-04, -5.62917E-04        &  
     &  ,  4.90244E-04, -3.28170E-04,  1.20767E-04,  1.09471E-04        &  
     &  , -1.02587E-04, -7.58733E-05,  7.46128E-05,  7.11222E-05        &  
     &  , -5.68686E-05,  4.97904E-05,  3.14644E-05,  2.83616E-05        &  
     &  , -2.66163E-05, -2.63254E-05,  2.50164E-05,  2.46285E-05/          
      DATA F/                                                           &  
     &     1.5324946E-04,  1.5814864E-04,  1.1719011E-04                &  
     &  ,  1.5506174E-04,  2.1733392E-04,  1.5016256E-04                &  
     &  ,  2.1170962E-04,  1.5633636E-04,  1.4835028E-04                &  
     &  ,  2.0692488E-04,  2.1252514E-04,  2.2999289E-04                &  
     &  ,  3.0649899E-04,  3.1139817E-04,  4.8991877E-06                &  
     &  ,  3.6059331E-05,  2.7043965E-04,  1.8122966E-06                &  
     &  ,  6.4084427E-05,  3.0341210E-04,  3.0831127E-04                &  
     &  ,  3.7058338E-04,  2.2211866E-04,  4.0958519E-05/                  
      DATA D/                                                           &  
     &     4.4041E+00,  4.9093E+00,  2.2451E+00,  5.1167E+00            &  
     &  ,  2.7912E-01,  4.6115E+00,  5.3935E+00,  4.1966E+00            &  
     &  ,  3.8990E+00,  4.7014E+00,  5.5397E+00,  5.5896E+00            &  
     &  ,  2.5251E+00,  3.0303E+00,  5.0517E-01,  2.1589E+00            &  
     &  ,  3.6608E-01,  7.1253E-01,  2.1582E+00,  2.7325E+00            &  
     &  ,  3.2376E+00,  4.6833E+00,  9.7121E-01,  2.6640E+00/              
c                                                                          
c   -----------------------------------------------------------------      
c    Table 4                                                               
c                                                                          
      DATA M/                                                           &  
     &     1.8607980E-02,  1.6275220E-02, -1.3006600E-02                &  
     &  ,  9.8882900E-03, -3.3670000E-03,  3.3307700E-03                &  
     &  , -2.3540000E-03,  1.4001500E-03,  1.0070000E-03                &  
     &  ,  8.5700000E-04,  6.4990000E-04,  5.9900000E-04                &  
     &  ,  3.7800000E-04, -3.3700000E-04,  2.7600000E-04                &  
     &  ,  1.8200000E-04, -1.7400000E-04, -1.2400000E-04                &  
     &  ,  1.2500000E-05/                                                  
      DATA G/                                                           &  
     &     2.0397105E-05,  3.5614854E-05,  8.6574454E-05                &  
     &  ,  8.3487563E-05,  8.1675266E-05,  2.5205846E-05                &  
     &  ,  8.8386751E-05,  1.2710243E-04,  3.0830121E-05                &  
     &  ,  7.8588375E-05,  1.4860417E-05,  8.0400672E-05                &  
     &  ,  8.9661345E-05,  3.0014587E-05,  9.1473642E-05                &  
     &  ,  8.4481533E-05,  2.9990579E-05,  8.9290274E-05                &  
     &  ,  3.2378912E-06/                                                  
      DATA B/                                                           &  
     &     5.0053E-01,  3.3839E+00,  5.3852E+00,  5.5925E+00            &  
     &  ,  4.8800E+00,  1.5230E+00,  6.0977E+00,  2.2481E+00            &  
     &  ,  2.6918E+00,  5.0874E+00,  2.0054E+00,  5.8001E+00            &  
     &  ,  5.1778E+00,  2.5455E+00,  5.8903E+00,  2.6587E+00            &  
     &  ,  2.2151E+00,  3.6812E+00,  1.2585E+00/                           
c                                                                          
c   -----------------------------------------------------------------      
c    Table 5                                                               
c                                                                          
      DATA C/                                                           &  
     &     3.58327E-02,  1.23877E-02,  9.80662E-03, -9.56853E-03        &  
     &  ,  6.01280E-03,  4.62449E-03, -4.51725E-03,  4.22942E-03        &  
     &  ,  2.93967E-03, -2.40482E-03/                                      
      DATA H/                                                           &  
     &     1.5324946E-04,  1.5814864E-04,  1.1719011E-04,  3.0868911E-06&  
     &  ,  1.5506174E-04,  1.5217749E-05,  1.5016256E-04,  2.1733392E-04&  
     &  ,  4.8087409E-06,  1.8122966E-06/                                  
      DATA R/                                                           &  
     &     4.4041E+00,  4.9093E+00,  2.2451E+00,  6.0756E+00            &  
     &  ,  5.1167E+00,  2.8833E+00,  4.6115E+00,  2.7912E-01            &  
     &  ,  1.0225E+00,  7.1253E-01/
                                                                           
c     The length of the calendar year may be set for a 360-day calendar    
c      (as is often used in climate runs),                                 
c      or for a real Gregorian calendar which has 365 days in              
c      non-leap years and 366 in leap years.                               

                                                                           
      IF (LCAL360) THEN                                                    
                                                                           
        DINY=360.0                                                         
                                                                           
      ELSE                                                                 
c      Is this a leap year?                                                
        IF (mod(year,4)   .eq. 0 .AND.                                     
     &     (mod(year,400) .eq. 0 .OR. mod(year,100) .ne. 0)) then          
                                                                           
          DINY = 366.0                                                     
                                                                           c      Is this a normal year?                                              
        ELSE                                                               
                                                                           
          DINY = 365.0                                                     
                                                                           
        END IF                                                             
      END IF                                                               
                                                                           
      IF (L_SEC_VAR_FILE .OR. L_SEC_VAR_ONLINE) THEN
c     The orbital elements are normally set to default values, but         
c     secular variations may be required in some longer climate runs.      
c  
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     r.s.smith
c     Eccentricity, obliquity and the longitude of perihelion can
c     be calculated online for +/- 1Myr via Berger78 formula. For
c     older periods, precalculated tables such as Laskar04 can be used.
c     Assumes a certain format for the tables (note: dates must
c     go back in time, as in Laskar paleo tables):
c      KYEAR FROM
c      <offset year>        ECC   OBLQ     PERH_MOVING_VE
c       ---------------------------------------------------------
c      <kyear from offset> <ecc> <oblq>  <long.perih. rel. to v.e.>
ccccccccccccccccccccccccccccccccccccccccccccccccc

      IF (L_SEC_VAR_FILE) THEN !L_SEC_VAR_FILE, use precalculated values

       SEC_VARUNIT=615

       OPEN(SEC_VARUNIT,FILE=SEC_VAR_FILE,IOSTAT=ICODE)
       REWIND(SEC_VARUNIT)

       READ(SEC_VARUNIT,*,IOSTAT=ICODE)DUMMY
       READ(SEC_VARUNIT,*,IOSTAT=ICODE)TABLEYEAR_REF
       READ(SEC_VARUNIT,*,IOSTAT=ICODE)DUMMY

       TABLEYEAR = REAL( YEAR - TABLEYEAR_REF )                              

       do while (icode.eq.0)
         READ(SEC_VARUNIT,*,IOSTAT=ICODE)FILEYEAR1,FE1,FOBLQ1,FLPH1
         READ(SEC_VARUNIT,*,IOSTAT=ICODE)FILEYEAR2,FE2,FOBLQ2,FLPH2

         FILEYEAR1=FILEYEAR1*1000.
         FILEYEAR2=FILEYEAR2*1000.

         if (tableyear.le.FILEYEAR1 .AND. tableyear.ge.FILEYEAR2) 
     &       icode=9999

         BACKSPACE(SEC_VARUNIT)
       end do

       if (ICODE.ne.9999) then
         write(6,*)"orbital lookup: year not found/file error",
     &             TABLEYEAR,icode
         stop
       end if


c linearly interpolate between the years we've got 

         DY2=1-( (TABLEYEAR-FILEYEAR2)/(FILEYEAR1-FILEYEAR2) )
         DY1=1-( (FILEYEAR1-TABLEYEAR)/(FILEYEAR1-FILEYEAR2) )

         E            =   FE1*DY1 + FE2*DY2
         OBLQ         =FOBLQ1*DY1 + FOBLQ2*DY2

c CAREFUL HERE - doing mod(GAMMA,2PI) later goes wrong when LPH
c uses a linear interpolation as above due to LPH wrapping at 0,2PI
c in the table
c LPH decreases as we go back in time. If 2 is greater than 1,
c then we must have wrapped. Use the general gradient to get a
c -ve value to interpolate towards - the wrap on GAMMA later will
c sort the signs out in the end
         if (FLPH2 .gt. FLPH1) then
           READ(SEC_VARUNIT,*,IOSTAT=ICODE)DUMMY
           READ(SEC_VARUNIT,*,IOSTAT=ICODE)FILEYEAR3,FE2,FOBLQ2,FLPH3
           FILEYEAR3=FILEYEAR3*1000.

           FLPH2=(FLPH1-(FLPH2-FLPH3))*
     &           (FILEYEAR1-FILEYEAR2)/(FILEYEAR2-FILEYEAR3)

         end if

         LPH_MOVING_VE= FLPH1*DY1 + FLPH2*DY2

ccccccccccccccccccccccccccccccccccccccccccccccccc

      ELSE ! L_SEC_VAR_ONLINE, calculate online, not from table

                                                                           
        YEAR_OFFSET = REAL( YEAR - YEAR_REF )                              
        if (ABS(year_offset).GT.1e6) then
          write(6,*) year_offset,
     &  "ONLINE SOLUTION FOR SOLAR PARAMS ONLY VALID FOR +/- 1MYr BP"
          stop
        end if
                                                                           
c       Obliquity: (Equation 1 from Berger 1978)                           
                                                                           
        OBLQ = OBLQ_CNST                                                   
        DO I=1, N_TERM_OBQ                                                 
          OBLQ = OBLQ+A(I)*COS(F(I)*YEAR_OFFSET+D(I))                      
        END DO                                                             
                                                                           
c       Eccentricity: this is better computed from its components          
c       than directly from the series.(Equation (4) of Berger 1978).       
                                                                           
        ECN_SN = M(1) * SIN (G(1) * YEAR_OFFSET + B(1))                    
        ECN_CN = M(1) * COS (G(1) * YEAR_OFFSET + B(1))                    
                                                                           
        DO I=2, N_TERM_ECN_LPH                                             
          ECN_SN = ECN_SN + M(I) * SIN (G(I) * YEAR_OFFSET + B(I))         
          ECN_CN = ECN_CN + M(I) * COS (G(I) * YEAR_OFFSET + B(I))         
        END DO                                                             
        E = SQRT(ECN_SN*ECN_SN+ECN_CN*ECN_CN)                              
                                                                           
c       We now obtain the longitude of the perihelion relative to the      
c       fixed equinox.                                                     
                                                                           
        LPH_FIXED_VE = ATAN2 (ECN_SN,ECN_CN)                               
                                                                           
c       The longitude of perihelion and                                    
c        the supplement of the longitude of the perihelion relative to     
c        the actual vernal equinox requires the general precession.        
                                                                           
c      General Precession.                                                 
        GN_PRCS = LIN_RATE_GN_PRCS * YEAR_OFFSET + GN_PRCS_CNST            
        DO I=1, N_TERM_GN_PRCS                                             
          GN_PRCS = GN_PRCS + C(I) * SIN (H(I) * YEAR_OFFSET + R(I))       
        END DO                                                             
                                                                           
        LPH_MOVING_VE = LPH_FIXED_VE + GN_PRCS

      END IF  !L_SEC_VAR_FILE, lookup or online calc identical from here
ccccccccccccccccccccccccccccccccccccccccccccccccc

c      Supplement of the longitude of the perihelion                       
        GAMMA = PI - LPH_MOVING_VE
                                                                           
c       Time of perihelion: The time at which an object is at perihelion   
c        (its closest distance to the sun).                                
c       The time of perihelion is inferred from the date of                
c        the vernal equinox using the Gregorian calendar.                  
c                                                                          
c      Calculate the date of the vernal equinox.                           
c       First we need to:                                                  
c        Calculate the no of leap days between year & year_ref_ve.         
c        This needs to be corrected when using the Gregorian calendar.     
c         by adding (DINY-366.0) when the year_ref_ve is a leap year or    
c         by adding (DINY-365.0) when the year_ref_ve is a normal year.    
c        This correction is done when the DATE_VE is calculated below!     
c                                                                          
c        In the calculation of NO_LEAP_DAYS below, the divisions of type   
c         'YEAR'/x (where x is 4, 100 or 400) are integer computations.    
c         These integers are then subtracted and the resulting integer     
c         is then converted to a real.                                     
                                                                           
        NO_LEAP_DAYS = ( TropYearLength - 365.0)                           
     &    * REAL( YEAR     - YEAR_REF_VE     )                             
     &    - REAL( YEAR/4   - YEAR_REF_VE/4   )                             
     &    + REAL( YEAR/100 - YEAR_REF_VE/100 )                             
     &    - REAL( YEAR/400 - YEAR_REF_VE/400 )                             
                                                                           
c      Now we can calculate the date of the vernal equinox!                
c      Because the date of the vernal equinox is varying with the year,    
c      we have to keep track of its position in the sky.                   
c      In order to accomodate a time varying vernal equinox when using     
c      a 360-day year, we still have to calculate the difference in        
c      the vernal equinox depending on leap years, normal years and the    
c      difference between the length of the tropical year and the          
c      "normal" year and then we adjust this by multiplying the            
c      DATE_VE by 360/(length of tropical year).                           
c                                                                          
c      Is a 360 day calendar being used?                                   
                                                                           
        IF (LCAL360) THEN                                                  
c         DATE_VE = DATE_VE_DFLT + NO_LEAP_DAYS                            
c         DATE_VE = DATE_VE * DINY / TropYearLength                        

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c for paleo runs it's far more useful to use the VE as a fixed point
c to define your year's forcing around
          DATE_VE = DATE_VE_DFLT * DINY / TropYearLength                        
                                                                           
                                                                           
c      Is a 365 day calendar being used?                                   
        ELSE                                                               
                                                                           
c        Is the epoch reference year a leap year?                          
                                                                           
          IF (mod(YEAR_REF_VE,4)   .eq. 0 .AND.                            
     &       (mod(YEAR_REF_VE,400) .eq. 0 .OR.                             
     &        mod(YEAR_REF_VE,100) .ne. 0)) THEN                           
                                                                           
            DATE_VE = DATE_VE_DFLT + (NO_LEAP_DAYS + (DINY - 366.0))       
                                                                           
c        Is the epoch reference year a normal year?                        
                                                                           
          ELSE                                                             
                                                                           
            DATE_VE = DATE_VE_DFLT + (NO_LEAP_DAYS + (DINY - 365.0))       
                                                                           
          END IF                                                           
        END IF                                                             
                                                                           
        BETA = SQRT(1.0E+00-E*E)                                           
        EE1  = (0.5*E + 0.125*E*E*E)*(1.0 + BETA)                          
        EE2  = -0.25*E*E* (0.5 + BETA)                                     
        EE3  = 0.125*E*E*E*((1.0/3.0) + BETA)                              
        MEAN_ANOM_VE = GAMMA - 2.0E+00 * (                                 
     &      EE1 * SIN (GAMMA)                                              
     &    + EE2 * SIN (2.0 * GAMMA)                                        
     &    + EE3 * SIN (3.0 * GAMMA)                                        
     &    )                                                                
                                                                           
        TAU0 = DATE_VE - MEAN_ANOM_VE * TropYearLength/(TWOPI)             

      ELSE ! neither L_SEC_VAR_ONLINE or _FILE, just use default values

        E     = E_DFLT
        OBLQ  = OBLQ_DFLT
        GAMMA = PI - LPH_DFLT
        TAU0  = TAU0_DFLT

      END IF 

                                                                           
c     If using a 360-day calendar the time of the perihelion is            
c     adjusted.                                                            
      IF (LCAL360) THEN                                                    
        TAU0 = TAU0*(360.0/TropYearLength)+0.71                            
      ENDIF                                                               

ccccccccccccccccccccccccccccccccccccccccccccccccc
c r.s.smith - added to get numbers similar to those
c in Julia/Michel's offline lookup table for paleodates
ccccccccccccccccccccccccccccccccccccccccccccccccc
      do while (tau0 .lt.0)
        tau0=tau0+DINY
      end do
      do while (tau0 .gt.DINY)
        tau0=tau0-DINY
      end do
      do while (gamma .lt.0)
        gamma=gamma+2*PI
      end do
      do while (gamma .gt.2*PI)
        gamma=gamma-2*PI
      end do
c                                                                          
      RETURN                                                              
      END                                                                 
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  orh0f406
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID ORH0F406
*/ U.M. 4.6 unix / source code change form / header   version 06/01/99
*/ CODE WRITERS MUST READ THE ACCOMPANYING INSTRUCTIONS FOR THIS BUILD:
*/  - See http://fr0800/umdoc/hegui/t3e4.6.html#chgfinst
*/ 
*/SOC:-> Correct the unconditional use of ice model code
*/->: Protect calls to SWAPBOUNDS which are only relevant
*/    when the ice model is switched on AND certain diagnostics
*/    are requested. Also dimenstion associated array CARYSALT 
*/    conditionally to avoid wasting memory.
*/-> applicable to UM version 4.5
*/
*/ Has an entry been lodged in the Problem Reporting System? [Y|N]      
*/
*/ THIS CODE IS INTENDED FOR INCLUSION IN THE 4.6 BUILD      [Y]
*/ .....................................................................
*/   Author[s]:-> E-mail:-> rshill@meto.gov.uk 
*/ Reviewer[s]:-> E-mail:-> @meto.gov.uk
*/
*/    "I have checked this change. When provided, the advance design 
*/  specification was agreed and adequate, and the new code conforms to
*/  Unified Model standards."
*/
*/  DESIGN SPEC. WAS REVIEWED ON: ......      REVIEWER[S] SIGNATURES
*/                                            ----------------------
*/    DATE CODE REVIEWED: ......
*/  .....................................................................
*/
*/  WILL CHANGES AFFECT ANCILLARY FILES?         [N]
*/  ARE ANY CHANGES TO STASHMASTER FILES NEEDED? [N] 
*/  USER INTERFACE ACTION REQUIRED?              [N]
*/ 
*/  TESTED IN CONFIGURATIONS:-> Ocean (LAM atlantic)
*/  TESTS RUN BY [PERSON]:->
*/ 
*/  WILL THE CHANGES SLOW DOWN THE MODEL?        [N]
*/  -> Will speed up very slightly or at worst remain the same
*/     depending on presence of ice model.
*/  CHANGES WILL INCREASE MEMORY CONSUMPTION?    [N]   
*/  -> Will reduce memory consumption in the routine OSWAPDIAGS
*/     if no ice model.
*/
*/ | Re-start dumps bit compare with those created without the change 
*/ V MARK [Y| ] BELOW; leave rest of lines untouched.
*/
*/   Control Code    loses bit comparison
*/   Atmosphere (assuming same science options chosen)   loses b.c.
*/   Ocean       loses bit comparison
*/   Wave        loses bit comparison
*/   Reconfiguration   loses bit comparison
*/   Diagnostics      lose bit comparison
*/ For Y2K compliance checking:  
*/ DOES THIS CHANGE INTERACT WITH DATE CALCULATIONS IN ANY WAY? [N]   
*/ 
*/  SECTIONS (TO BE) CHANGED:
*/
*/  SECTIONS (TO BE) DELETED? 
*/
*/  NEW SECTIONS?  Fill in form http://www-hc/~hadmk/STASHmaster_change.html,
*/  and give section numbers below:
*/  
*/  *DEFS ADDED OR REMOVED: 
*/
*/  **Existing** decks being changed [with *I, *D, *B directives]
*/ -> OSWAPDIA
*/
*/  Decks being created or purged [with *DECK, *COMDECK, *PURGEDK]
*/ *......K  Deck name   Section#.vr
*/ -> 
*/ ......................................................................
*/ ANY REFERENCES TO EXTERNAL DOCUMENTS-> instead of design spec.
*/  ...OR ... ADVANCE DESIGN SPECIFICATION (optional) 
*/ ->    
*//////////////////////////////////////////////////////////////////////// 
*/
*/ Move unprotected ice model related SWAPBOUNDS calls to 
*/ a more appropriate position to ensure they are only
*/ called when needed. Also conditionally dimension
*/ CARYSALT which is only fully dimensioned when ice model
*/ is included.
*/ The net effect of all this is that timing figures
*/ from runs which dont need any of the diagnostics in
*/ OSWAPDIAG are not distorted by needless barriers
*/ and unwanted subroutine calls plus we save
*/ a bit of memory.
*/
*/ R. Hill   March 1999.
*/
*DECLARE OSWAPDIA
*D OOM1F405.718 
      REAL CARYSALT(IMT_ICE,JMT_ICE)
     &    ,OCEANHEATFLUX(IMT,JMT)
*D OOM1F405.741,OOM1F405.742
*I OSWAPDIA.250
         CALL SWAPBOUNDS(CARYSALT,IMT,JMT,O_EW_HALO,O_NS_HALO,1)             
         CALL SWAPBOUNDS(OCEANHEATFLUX,IMT,JMT,O_EW_HALO,O_NS_HALO,1)     
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  osy1f405_change
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID MJRADVCH
*DECLARE ADVSRCE
*D ADVSRCE.391,ADVSRCE.394
          x_div(i,k)=(flux(i,k)-flux(i+1,k))*dxr(i)*recip_cos(j)
          source(i,k)=x_div(i,k)                        
*D ADVSRCE.536,ADVSRCE.542
          y_div(i,k)=(fluxs(i,k)-fluxn(i,k))*dyr(j)  
          source(i,k)=source(i,k)+y_div(i,k)        
*D ADVSRCE.646,ADVSRCE.648
          z_div(i,k)=(flux(i,k+1)-flux(i,k))/dz(k)  
          source(i,k)=source(i,k)+z_div(i,k)       
  
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  params_in_namelist_famous.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  polescon
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID POLESCON
*/ 
*/ This is a FAMOUS mod. 
*/ 
*/ Rather than change POLAR & expect it not to matter for OMEGA_P, use
*/   new routine for all other calls ...
*/
*DECK NEWPOLAR
      SUBROUTINE NEWPOLAR(FIELD,
*CALL ARGFLDPT
     &                 FIELD_SIZE,ROW_LENGTH,N_LEVELS)
*COPY POLAR1A,APB2F401.105,109
*COPY POLAR1A,APB2F401.114,117
*CALL TYPFLDPT
*COPY POLAR1A,APB2F401.119,122
*IF DEF,MPP
      INTEGER info
*ENDIF
*COPY POLAR1A,APB2F401.131,141
*IF -DEF,MPP
*COPY POLAR1A,APB2F401.143,144
          MEAN_NP(K)=MEAN_NP(K)+FIELD(I+ROW_LENGTH,K)
          MEAN_SP(K)=MEAN_SP(K)+FIELD(I+FIELD_SIZE-2*ROW_LENGTH,K)
*COPY POLAR1A,APB2F401.147,152
          FIELD(I,K)=MEAN_NP(K)
          FIELD(I+FIELD_SIZE-ROW_LENGTH,K)=MEAN_SP(K)
*COPY POLAR1A,APB2F401.156,APB2F401.159
*ELSE
      IF (at_top_of_LPG) THEN
*IF DEF,REPROD
        CALL GCG_RVECSUMR(FIELD_SIZE,ROW_LENGTH-2*EW_Halo,
*ELSE
        CALL GCG_RVECSUMF(FIELD_SIZE,ROW_LENGTH-2*EW_Halo,
*ENDIF
     &                    2*ROW_LENGTH+1+EW_Halo,N_LEVELS,
     &                    FIELD,GC_ROW_GROUP,info,MEAN_NP)
*COPY POLAR1A,APB2F401.165,APB2F401.167
            FIELD(I,K)=MEAN_NP(K)
*COPY POLAR1A,APB2F401.169,APB2F401.173
*IF DEF,REPROD
        CALL GCG_RVECSUMR(FIELD_SIZE,ROW_LENGTH-2*EW_Halo,
*ELSE
        CALL GCG_RVECSUMF(FIELD_SIZE,ROW_LENGTH-2*EW_Halo,
*ENDIF
     &                    FIELD_SIZE-3*ROW_LENGTH+1+EW_Halo,N_LEVELS,
     &                    FIELD,GC_ROW_GROUP,info,MEAN_SP)
*COPY POLAR1A,APB2F401.177,APB2F401.179
            FIELD(I,K)=MEAN_SP(K)
*COPY POLAR1A,APB2F401.181,APB2F401.183
*ENDIF
*COPY POLAR1A,APB2F401.185,POLAR1A.66
*/
*DECLARE ATMPHY1
*D ATMPHY1.544,556
*D APB2F401.63
      CALL NEWPOLAR(D1(JTHETA(1)),
*D APB2F401.65,66
     &           P_FIELD,ROW_LENGTH,P_LEVELS)
*D ATMPHY1.573,595
*D APB2F401.70,71
      CALL NEWPOLAR(D1(JQ(1)),
*D APB2F401.73,74
     &           P_FIELD,ROW_LENGTH,Q_LEVELS)
*D APB2F401.78,79
      CALL NEWPOLAR(D1(JQCL(1)),
*D APB2F401.81,82
     &           P_FIELD,ROW_LENGTH,Q_LEVELS)
*D APB2F401.86,87
      CALL NEWPOLAR(D1(JQCF(1)),
*D APB2F401.89,90
     &           P_FIELD,ROW_LENGTH,Q_LEVELS)
*/
*/ Carry on using old sort of CALL from BL_CTL1 QTADV1A QTADV1C
*/   THADV1A THADV1C THADV1D THADV1E THQDIF1A THQDIF1B THQDIF1C
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  port_conv_f.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID PORTCONVF
*/--------------------------------------------------------
*/
*/ Backport vn5.5 mod gqz3505 and vn6.0 mod gqz1600 to vn4.5
*/ Original author frqz (Paul Dando)
*/
*/ Jeff Cole 14/11/03
*/
*/ Provides portable data conversion routines to replace
*/ the Cray-specific routines on non-Cray platforms
*/--------------------------------------------------------
*/
*/-----
*DECLARE COEX1A
*/-----
*D UIE2F403.4
*IF DEF,CRAY
      INTEGER :: CRI2IBM
      INTEGER :: IBM2CRI
*ELSE
      INTEGER :: IEEE2IBM
      INTEGER :: IBM2IEEE
*ENDIF
*I UIE2F403.3
!LL   5.5  28/02/03 Insert code for portable data conversion routines
!LL                 to replace Cray-specific CRI2IBM etc.     P.Dando
!     6.0  10/09/03 Conversion of portable data conversion routines
!                   (IEEE2IBM etc) into functions with error return
!                   codes matching those of CRAY routines.   P.Dando
!LL
*I COEX1A.67
*IF DEF,CRAY
*I UIE2F403.7
*ELSE
          IER=IEEE2IBM(2,1,ICOMP(1),32,ISC,1,64,32)
          IER=IEEE2IBM(2,1,ICOMP(2),0,IX,1,64,16)
          IER=IEEE2IBM(2,1,ICOMP(2),16,IY,1,64,16)
*ENDIF
*I COEX1A.100
*IF DEF,CRAY
*I UIE2F403.10
*ELSE
          IER=IBM2IEEE(2,1,ICOMP(1),32,ISC,1,64,32)
          IER=IBM2IEEE(2,1,ICOMP(2),0,IX,1,64,16)
          IER=IBM2IEEE(2,1,ICOMP(2),16,IY,1,64,16)
*ENDIF
*D UIE2F403.11
*IF DEF,CRAY
      INTEGER :: CRI2IBM
      INTEGER :: IBM2CRI
*ELSE
      INTEGER :: IEEE2IBM
      INTEGER :: IBM2IEEE
*ENDIF
*I COEX1A.199
*IF DEF,CRAY
*I UIE2F403.14
*ELSE
              IERR1(JJ)=IEEE2IBM(3,1, JCOMP(1,JJ),0, BASE(JJ),1,64,32)
              IERR2(JJ)=IEEE2IBM(2,1,JCOMP(1,JJ),32,IBIT(JJ),1,64,16)
              IERR3(JJ)=IEEE2IBM(2,1,JCOMP(1,JJ),48,NOP(JJ),1,64,16)
*ENDIF
*I COEX1A.246
*IF DEF,CRAY
*I COEX1A.247
*ELSE
              CALL MOVEBYTES(JCOMP(1,JJ),1,NOB,ICOMP(ICX),IST)
*ENDIF
*I COEX1A.252
*IF DEF,CRAY
*I UIE2F403.15
*ELSE
          IER2=IEEE2IBM(2,1,ICOMP(1),0,NUM,1,64,32)
*ENDIF
*I COEX1A.273
*IF DEF,CRAY
*I UIE2F403.18
*ELSE
                  IERR=IBM2IEEE(3,1,ICOMP(ICX),32,BASE(JJ),1,64,32)
                  IERR=IBM2IEEE(2,1,ICOMP(ICX+1),0,IBIT(JJ),1,64,16)
                  IERR=IBM2IEEE(2,1,ICOMP(ICX+1),16,NOP(JJ),1,64,16)
*ENDIF
*I COEX1A.280
*IF DEF,CRAY
*I UIE2F403.21
*ELSE
                  IERR=IBM2IEEE(3,1,ICOMP(ICX),0,BASE(JJ),1,64,32)
                  IERR=IBM2IEEE(2,1,ICOMP(ICX),32,IBIT(JJ),1,64,16)
                  IERR=IBM2IEEE(2,1,ICOMP(ICX),48,NOP(JJ),1,64,16)
*ENDIF
*I COEX1A.306
*IF DEF,CRAY
*I COEX1A.307
*ELSE
              CALL MOVEBYTES(ICOMP(ICX),IST,NOB,JCOMP(1,JJ),1)
*ENDIF
*I COEX1A.815
*IF DEF,CRAY
*I COEX1A.816
*ELSE
          CALL MOVEBITS(INUM,1,1,ICOMP(1),ISCOMP)
*ENDIF
*I COEX1A.824
*IF DEF,CRAY
*I COEX1A.825
*ELSE
      CALL MOVEBITS(INUM,ISNUM,NUM,ICOMP(1),ISCOMP)
*ENDIF
*I COEX1A.873
*IF DEF,CRAY
*I COEX1A.874
*ELSE
          CALL MOVEBITS(ICOMP,ISCOMP,NUM,INUM,ISNUM)
*ENDIF
*I COEX1A.880
*IF DEF,CRAY
*I COEX1A.881
*ELSE
          CALL MOVEBITS(ICOMP,ISCOMP,1,INUM,1)
*ENDIF
*I COEX1A.885
*IF DEF,CRAY
*I COEX1A.886
*ELSE
          CALL MOVEBITS(ICOMP,ISCOMP,NUM,INUM,ISNUM)
*ENDIF
*I COEX1A.901
*IF DEF,CRAY
*I COEX1A.902
*ELSE
          CALL MOVEBITS(ICOMP,ISCOMP,NUM,INUM,ISNUM)
*ENDIF
*I GBCQF405.111
*IF DEF,CRAY
*I GBCQF405.112
*ELSE
      INTEGER IEEE2IBM
*ENDIF
*I GBCQF405.123
*IF DEF,CRAY
*I GBCQF405.126
*ELSE
      IERR=IEEE2IBM(2,1,ICOMP(1,1),32,ISC,1,64,32)
      IERR=IEEE2IBM(2,1,ICOMP(2,1),0,IX,1,64,16)
      IERR=IEEE2IBM(2,1,ICOMP(2,1),16,IY,1,64,16)
*ENDIF
*I GBCQF405.281
*IF DEF,CRAY
*I GBCQF405.282
*ELSE
      INTEGER IEEE2IBM
*ENDIF
*I GBCQF405.362
*IF DEF,CRAY
*D GBCQF405.365
            IERR=CRI2IBM(2,1,JCOMP_PROC(1,JJ),48,NOP_PROC(JJ),1,64,16)
*ELSE
            IERR=IEEE2IBM(3,1,JCOMP_PROC(1,JJ),0, BASE_PROC,1,64,32)
            IERR=IEEE2IBM(2,1,JCOMP_PROC(1,JJ),32,IBIT_PROC,1,64,16)
            IERR=IEEE2IBM(2,1,JCOMP_PROC(1,JJ),48,NOP_PROC(JJ),
     &                       1,64,16)
            IERR=0
*ENDIF
*I GBCQF405.437
*IF DEF,CRAY
*I GBCQF405.438
*ELSE
            CALL MOVEBYTES(JCOMP(1,JJ,K),1,NOB,ICOMP(ICX,K),IST)
*ENDIF

*I GBCQF405.444
*IF DEF,CRAY
*I GBCQF405.445
*ELSE
          IERR=IEEE2IBM(2,1,ICOMP(1,K),0,NUM,1,64,32)
*ENDIF
*/-----
*DECLARE FIELDCOS
*/-----
*B FIELDCOS.21
!LL  5.5  28/02/03  Insert code for portable data conversion routines
!LL                 to replace Cray-specific CRI2IBM etc.     P.Dando
!    6.0  10/09/03  Conversion of portable data conversion routines
!                   (IEEE2IBM etc) into functions with error return
!                   codes matching those of CRAY routines.    P.Dando
*D UIE1F402.1,UIE1F402.2
      EXTERNAL READFF,INT_FROM_REAL,TIME2SEC,SEC2TIME
*IF DEF,CRAY
      EXTERNAL CRI2IBM
      INTEGER CRI2IBM
*ELSE
      INTEGER IEEE2IBM
*ENDIF
      INTEGER INT_FROM_REAL
*I PS050793.211
*IF DEF,CRAY
*I UIE1F402.4
*ELSE
      IER = IEEE2IBM(2,LEN_ILABEL,IBM_LABEL(IBM_ADDR),BIT_OFF,
     &               ILABEL,1,64,32)
*ENDIF
*I PS050793.216
*IF DEF,CRAY
*I UIE1F402.6
*ELSE
      IER = IEEE2IBM(3,LEN_RLABEL,IBM_LABEL(IBM_ADDR),BIT_OFF,
     &               RLABEL,1,64,32)
*ENDIF
*I URR2F405.11
*IF DEF,CRAY
*I UIE1F402.10
*ELSE
            IER = IEEE2IBM(3,NUM_VALUES-IEXTRAW,IBM_FIELD,
     &                  BIT_OFF,FIELD,1,64,32)
*ENDIF
*I PS050793.220
*IF DEF,CRAY
*I UIE1F402.12
*ELSE
          IER = IEEE2IBM(2,NUM_VALUES-IEXTRAW,IBM_FIELD,
     &               BIT_OFF,FIELD,1,64,32)
*ENDIF
*I PS050793.221
*IF DEF,CRAY
*I UIE1F402.14
*ELSE
          IER = IEEE2IBM(5,NUM_VALUES-IEXTRAW,IBM_FIELD,
     &               BIT_OFF,FIELD,1,64,32)
*ENDIF
*I FIELDCOS.560
*IF DEF,CRAY
*I UIE1F402.16
*ELSE
          IER=IEEE2IBM(2,1,IBM_FIELD(IBM_ADDR),BIT_OFF,
     &             FIELD(ADDR),1,64,32)
*ENDIF
*I FIELDCOS.577
*IF DEF,CRAY
*I UIE1F402.18
*ELSE
          IER=IEEE2IBM(3,DATA_VALUES,IBM_FIELD(IBM_ADDR),
     &      BIT_OFF,FIELD(ADDR),1,64,32)
*ENDIF
*I PS050793.229
*IF DEF,CRAY
*I UIE1F402.20
*ELSE
        IER = IEEE2IBM(2,LEN_ILABEL,IBM_LABEL(IBM_ADDR),
     &               BIT_OFF,ILABEL,1,64,32)
*ENDIF
*I PS050793.239
*IF DEF,CRAY
*I UIE1F402.22
*ELSE
        IER = IEEE2IBM(3,LEN_RLABEL,IBM_LABEL(IBM_ADDR),
     &               BIT_OFF,RLABEL,1,64,32)
*ENDIF
*D UIE1F402.39,UIE1F402.40
      EXTERNAL READFF,INT_FROM_REAL,TIME2SEC,SEC2TIME
*IF DEF,CRAY
      EXTERNAL CRI2IEG
      INTEGER CRI2IEG
*ELSE
      INTEGER IEEE2IEG
*ENDIF
      INTEGER INT_FROM_REAL
*I PS050793.439
*IF DEF,CRAY
*I UIE1F402.42
*ELSE
      IER=IEEE2IEG(2,LEN_ILABEL,IEEE_LABEL(IEEE_ADDR),
     &           BIT_OFF,ILABEL,1,64,32)
*ENDIF
*I PS050793.444
*IF DEF,CRAY
*I UIE1F402.44
*ELSE
      IER=IEEE2IEG(3,LEN_RLABEL,IEEE_LABEL(IEEE_ADDR),
     &           BIT_OFF,RLABEL,1,64,32)
*ENDIF
*I GIE0F403.192
*IF DEF,CRAY
*I UIE1F402.46
*ELSE
          IER = IEEE2IEG(5,NUM_VALUES-IEXTRAW,IEEE_FIELD,
     &             BIT_OFF,FIELD,1,64,32)
*ENDIF
*I APS2F304.33
*IF DEF,CRAY
*I UIE1F402.48
*ELSE
          IER = IEEE2IEG(3,NUM_VALUES-IEXTRAW,IEEE_FIELD
     &             ,BIT_OFF,FIELD,1,64,32)
*ENDIF
*I FIELDCOS.1266
*IF DEF,CRAY
*I UIE1F402.50
*ELSE
          IER = IEEE2IEG(2,NUM_VALUES-IEXTRAW,IEEE_FIELD,
     &             BIT_OFF,FIELD,1,64,32)
*ENDIF
*I FIELDCOS.1273
*IF DEF,CRAY
*I UIE1F402.52
*ELSE
          IER = IEEE2IEG(5,NUM_VALUES-IEXTRAW,IEEE_FIELD,
     &             BIT_OFF,FIELD,1,64,32)
*ENDIF
*I FIELDCOS.1306
*IF DEF,CRAY
*I UIE1F402.54
*ELSE
          IER = IEEE2IEG(2,1,IEEE_FIELD(IEEE_ADDR),BIT_OFF,
     &             FIELD(ADDR),1,64,32)
*ENDIF
*I FIELDCOS.1323
*IF DEF,CRAY
*I UIE1F402.56
*ELSE
          IER=IEEE2IEG(3,DATA_VALUES,IEEE_FIELD(IEEE_ADDR),
     &     BIT_OFF,FIELD(ADDR),1,64,32)
*ENDIF
*I PS050793.454
*IF DEF,CRAY
*I UIE1F402.58
*ELSE
        IER=IEEE2IEG(2,LEN_ILABEL,IEEE_LABEL(IEEE_ADDR),
     &           BIT_OFF,ILABEL,1,64,32)
*ENDIF
*I PS050793.462
*IF DEF,CRAY
*I UIE1F402.60
*ELSE
        IER=IEEE2IEG(3,LEN_RLABEL,IEEE_LABEL(IEEE_ADDR),
     &           BIT_OFF,RLABEL,1,64,32)
*ENDIF
*/-----
*DECLARE PPTOANC1
*/-----
*B PPTOANC1.143
!   5.5   28/02/03 Insert code for portable data conversion routines
!                  to replace Cray-specific CRI2IBM etc.     P.Dando
!   6.0   10/09/03 Conversion of portable data conversion routines
!                  (IEEE2IBM etc) into functions with error return
!                  codes matching those of CRAY routines.    P.Dando
*I PPTOANC1.2341
*IF DEF,CRAY
*I PPTOANC1.2343
*ELSE
      integer ibm2ieee
*ENDIF
*I PPTOANC1.2367
*IF DEF,CRAY
*I PPTOANC1.2368
*ELSE
        ierr=ibm2ieee(3,rows*columns,levels_in,0,levels_array,
     &                   1,64,32)
*ENDIF
*I PPTOANC1.2606
*IF DEF,CRAY
*I PPTOANC1.2607
*ELSE
      integer ibm2ieee
*ENDIF
*I PPTOANC1.2644
*IF DEF,CRAY
*I PPTOANC1.2645
*ELSE
        ierr=ibm2ieee(3,rows*columns,levels_in,0,levels_array,
     &                   1,64,32)
*ENDIF
*I PPTOANC1.2854
*IF DEF,CRAY
*I PPTOANC1.2855
*ELSE
      INTEGER IBM2IEEE
*ENDIF
*I PPTOANC1.2865
*IF DEF,CRAY
*I PPTOANC1.2866
*ELSE
        ierr=ibm2ieee(3,rows*columns,datain,0,data_field,1,64,32)
*ENDIF
*I PPTOANC1.3385
*IF DEF,CRAY
*I PPTOANC1.3386
*ELSE
      integer ibm2ieee
*ENDIF
*I PPTOANC1.3394
*IF DEF,CRAY
*I PPTOANC1.3399
*ELSE
C Convert Integer part of header (Words 1-45)
        ier = ibm2ieee(2,45,pp_buffer,0,pp_int,1,64,32)

C Convert Real part of header (Words 46-64)
        ier = ibm2ieee(3,19,pp_buffer(23),32,pp_real,1,64,32)
*ENDIF
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  self_shade_safe.mod_4.5
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*IDENT SELF_SHADE
*/
*/ CDJ 20/11/98
*/
*/ Self shading mod for use in HadCM3LC.
*/  copy of:
*/          ~t20ti/Modsets/UM44/self_shade2.mod
*/
*/
*/ Ian J. Totterdell 981111
*/ Mod to put effects of self-shading into the simple
*/  light model currently used in HadOCC.
*/
*/-----------------------------------------------------------
*/
*DECK SOLSETSS
      subroutine SOLSETSS (IMT,KM,NT,DELTA_SI,DZ,GRAV_SI,
     &  RSOL,RZ,TR,ZDZ,ZDZZ,KFIX,DECAY,DELPSF,DELPSL,SOL_PEN)
c
      implicit none
c
*CALL CNTLOCN
*CALL OARRYSIZ
*CALL OTRACPNT
c
      integer
     &  IMT
     &, KM
     &, NT
     &, KFIX
c
      real
     &  DELTA_SI
     &, DZ(KM)
     &, GRAV_SI
     &, RSOL
     &, RZ(KM)
     &, TR(IMT,KM,NT)
     &, ZDZ(KM)
     &, ZDZZ(KM)
     &, DECAY(KM)
     &, DELPSF
     &, DELPSL(IMT,0:KM)
     &, SOL_PEN(IMT,0:KM)
c
c  External subroutine called.
c
      intrinsic
     &  EXP
c
      integer
     &  i
     &, k
     &, k_temp
c
      real
     &  chla
     &, eta1b
     &, eta1c
     &, eta1e(KM)
     &, eta1l(KM)
     &, eta1s
     &, eta2b
     &, eta2c
     &, eta2e(KM)
     &, eta2l(KM)
     &, eta2s
     &, fxa
     &, grav
     &, pebwcb
     &, pebwcc
     &, rsolc
     &, sclfct
c
*CALL OBIOCONST
c
      fxa = c2n_p*c_mol_wt/( c2chl*chl2pig )
c
      eta1b = 1.0
      eta1c = 0.02
      eta2b = 0.059
      eta2c = 0.02
c
      eta1b = eta1b*0.01
      eta1c = eta1c*0.01
      eta2b = eta2b*0.01
      eta2c = eta2c*0.01
c
      grav = GRAV_SI*100.0
      rsolc = 1.0 - RSOL
c
      DELPSF = -0.5*grav*DZ(1)*DZ(1)
c
c  Calculate layer of maximum penetration (200.0 m)
c
      k_temp = 1
      do k=1,KM
         k_temp = k
         if ( ZDZ(k) .gt. 200.0e2 ) go to 5
       enddo
    5 continue
      KFIX = k_temp
c
c  Loop over longitude and depth, calculating ETA, SOL_PEN and 
c   DELPSL as we descend.
c
      do i=1,IMT
c
         SOL_PEN(i,0) = 1.0
c
         eta1s = 0.0
         eta2s = 0.0
c
         do k=1,KFIX
c
            if( TR(i,k,PHYTO_TRACER) .ge. 0.0 )then
               chla = fxa*TR(i,k,PHYTO_TRACER)
             else
               chla = 0.0
             endif
c
            eta1l(k) = eta1b + eta1c*chla
            eta1s = eta1s - eta1l(k)*DZ(k)
            if( eta1s .lt. -180.0 )then
               eta1e(k) = 0.0
             else
               eta1e(k) = EXP(eta1s)
             endif
c
            eta2l(k) = eta2b + eta2c*chla
            eta2s = eta2s - eta2l(k)*DZ(k)
            if( eta2s .lt. -180.0 )then
               eta2e(k) = 0.0
             else
               eta2e(k) = EXP(eta2s)
             endif
c
            SOL_PEN(i,k) = RSOL*eta1e(k) + rsolc*eta2e(k)
c
c  End of loop (k=1,KFIX).
c
          enddo
c
         do k=KFIX,KM
            SOL_PEN(i,k) = 0.0
          enddo
c
c  Now, still in loop over i, calculate DELPSL.
c
         sclfct = -grav*DZ(1)/( 1.0 - SOL_PEN(i,1) )
         DELPSL(i,0) = 0.0
c
         pebwcc = RSOL*( 1.0 - eta1e(1) )/eta1l(1) +
     &     rsolc*( 1.0 - eta2e(1) )/eta2l(1) - ZDZ(1)*SOL_PEN(i,1)
         pebwcb = ( 1.0 - SOL_PEN(i,1) )*ZDZZ(1)*(DZ(1)/RZ(1))
         DELPSL(i,1) = ( pebwcb - pebwcc )*sclfct
c
         do k=2,KFIX
c
            pebwcc = pebwcc + RSOL*( eta1e(k-1) - eta1e(k) )/eta1l(k)
     &        + rsolc*( eta2e(k-1) - eta2e(k) )/eta2l(k)
     &        + ZDZ(k-1)*SOL_PEN(i,(k-1)) - ZDZ(k)*SOL_PEN(i,k)
            pebwcb = pebwcb +
     &        ( SOL_PEN(i,(k-1)) - SOL_PEN(i,k) )*ZDZZ(k)*(DZ(k)/RZ(k))
            DELPSL(i,k) = ( pebwcb - pebwcc )*sclfct
c
          enddo
c
         do k=(KFIX+1),KM
            DELPSL(i,k) = DELPSL(i,KFIX)
          enddo
c
c  End of loop (i=1,IMT).
c
       enddo
c
c  Decay scale for wind mixing energy.
c
      do k=1,KM
         DECAY(k) = EXP(-DZ(k)/(DELTA_SI*100.0))
       enddo
c  
c  Return from SOLSETSS.
c
      return
      end
*C SOLSETSS
*/
*/-----------------------------------------------------------
*/
*/ Will need some mods in TRACER to call the subroutine
*/  and use the output correctly.
*/
*DECLARE TRACER
*B OJP0F404.795
         IF (L_OBIOLOGY) THEN
           call SOLSETSS (IMT,KM,NT,DELTA_SI,DZ,GRAV_SI,RSOL,GAMMA_DZ,
     &       T,ZDZ,ZDZZ,KFIX_BIO,DECAY,DELPSF,DELPSL_LOC,SOL_PEN_BIO)
           do K=0,KM
             do I=1,IMT
               SOL_PEN_LOC(I,K) = SOL_PEN_BIO(I,K)
             enddo
           enddo
           KFIX = KFIX_BIO
         ELSE
*I OJP0F404.801
         END IF
*C TRACER
*/ End of modset.
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  timerupd_new.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID TIMUPD
*/
*/ TIMEF is an SGI intrinsic function and a UM subroutine name so to avoid
*/ confusion remove all calls to TIMEF, replace with call to TIMEFN.
*/ TIMEFN is the same as TIMEF subroutine with machine specific
*/ calls removed and replaced with f90 intrinsic function system_clock.
*/
*DECLARE TIMER1A
*D PXTIMEFN.26,PXTIMEFN.28
*D PXTIMEFN.29,PXTIMEFN.32
      CALL TIMEFN(ELPSTART)
*D PXTIMEFN.33,PXTIMEFN.36
      CALL TIMEFN(ELPEND)
*D PXTIMEFN.37,PXTIMEFN.40
         CALL TIMEFN(ELPEND)
*D PXTIMEFN.41,PXTIMEFN.44
      CALL TIMEFN(ELPSTART)
*D PXTIMEFN.45,PXTIMEFN.48
      CALL TIMEFN(ELPEND)
*D PXTIMEFN.49,PXTIMEFN.52
      CALL TIMEFN(ELPSTART)
*/
*DECLARE TIMER3A
*D PXTIMEFN.53,PXTIMEFN.59
      CALL TIMEFN(temp)
*/
*/ Replace the wallclock timer with vn5.X version
*/
*DECLARE TIMEFN2A
*D PXTIMEFN.14,PXTIMEFN.25
!----------------------------------------------------------------------
!
!+ Gets the elapsed time from the system

! Function Interface:
      SUBROUTINE TIMEFN(Get_Wallclock_Time)

      IMPLICIT NONE
!
! Description:
!   The system function SYSTEM_CLOCK is used to return the numbers of
!   seconds which have elapsed.
!
! Current Code Owner: Anette Van der Wal
!
! IE, PB, Richard Barnes  <- programmer of some or all of previous
!                            code or changes
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 5.3       24/04/01  Complete re-write to simplify.  A Van der Wal
!
! Code Description:
!   Language: FORTRAN 77 plus some Fortran 90 features
!   This code is written to UMDP3 v6 programming standards
!
! System component covered:
! System Task:
!
!- End of header
      REAL Get_Wallclock_Time

! Local variables
      INTEGER, SAVE :: Start_Count=-1
      INTEGER, SAVE :: Old_Count=0
      REAL, SAVE    :: Rollover=0.0
      INTEGER       :: Count, Count_Rate, Count_Max, Elapsed_Count
      REAL, SAVE    :: OneOver_Count_Rate=0.0

! Intrinsic procedures called:
      INTRINSIC SYSTEM_CLOCK

      CALL SYSTEM_CLOCK(Count=Count,Count_Rate=Count_Rate,
     &                  Count_Max=Count_Max)

      IF ((Old_Count .LT. Start_Count) .AND.
     &    ((Count .LT. Old_Count) .OR. (Count .GT. Start_Count))) THEN
        Rollover=Rollover+(REAL(Count_Max)/REAL(Count_Rate))
      ENDIF

      IF (Start_Count .EQ. -1) THEN
        Start_Count=Count
        OneOver_Count_Rate=1.0/REAL(Count_Rate)
      ENDIF

      Elapsed_Count=Count-Start_Count
      IF (Elapsed_Count .LT. 0) Elapsed_Count=Elapsed_Count+Count_Max

      Get_Wallclock_Time = Rollover+
     &                    (REAL(Elapsed_Count)*OneOver_Count_Rate)
      Old_Count=Count

      RETURN
      END
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  tjnowrit
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID TJNOWRIT,DC=WRITDM1A
*/
*/ Do not write various statements when writing dumps (mod to 4.2)
*/
*/ Name(s) of Author(s):    Tim Johns
*/ Name(s) of Reviewer(s):  Not applicable.
*/
*/       - Modules affected:  WRITDM1A, WRITHE1A
*/
*/       - Depends on the following mods for update:  GDG0F401, UDR1F401
*/
*I GDG0F401.1596  
!     4.2  19/12/96  Delete various write statements.   Tim Johns
*D WRITDM1A.112,WRITDM1A.113  
!     WRITE(6,'(/,'' WRITING UNIFIED MODEL DUMP ON UNIT'',I3)')NFTOUT
!     WRITE(6,'('' #####################################'',/)')
*D WRITDM1A.208,WRITDM1A.214  
!      WRITE(6,'('' '')')
!      IF (FIXHD(5).GE.6 .AND. FIXHD(5).LE.9) THEN ! AC/VarObs/Cx/Cov
!        WRITE(6,'('' OBSERVATION DATA'')')
!      ELSE
!        WRITE(6,'('' MODEL DATA'')')
!      ENDIF
!      WRITE(6,'('' '',I8,'' words long'')')FIXHD(161)
*D WRITDM1A.218
!     WRITE(6,'('' '')')
*DC WRITHE1A
*I UDR1F401.2     
!     4.2    19/12/96  Delete various write statements.  Tim Johns
*D WRITHE1A.294,WRITHE1A.296  
!      WRITE(6,'('' '')')
!      WRITE(6,'('' LEVEL DEPENDENT CONSTANTS'')')
!      WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(111)*FIXHD(112)
*D WRITHE1A.339,WRITHE1A.341  
!      WRITE(6,'('' '')')
!      WRITE(6,'('' ROW DEPENDENT CONSTANTS'')')
!      WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(116)*FIXHD(117)
*D WRITHE1A.383,WRITHE1A.385  
!      WRITE(6,'('' '')')
!      WRITE(6,'('' COLUMN DEPENDENT CONSTANTS'')')
!      WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(121)*FIXHD(122)
*D WRITHE1A.427,WRITHE1A.429  
!      WRITE(6,'('' '')')
!      WRITE(6,'('' FIELD DEPENDENT CONSTANTS'')')
!      WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(126)*FIXHD(127)
*D WRITHE1A.471,WRITHE1A.473  
!      WRITE(6,'('' '')')
!      WRITE(6,'('' EXTRA CONSTANTS'')')
!      WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(131)
*D WRITHE1A.515,WRITHE1A.517  
!      WRITE(6,'('' '')')
!      WRITE(6,'('' TEMPORARY HISTORY BLOCK'')')
!      WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(136)
*D WRITHE1A.560,WRITHE1A.562  
!      WRITE(6,'('' '')')
!      WRITE(6,'('' COMPRESSED FIELD INDEX NO 1'')')
!      WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(141)
*D WRITHE1A.603,WRITHE1A.605  
!      WRITE(6,'('' '')')
!      WRITE(6,'('' COMPRESSED FIELD INDEX NO 2'')')
!      WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(143)
*D WRITHE1A.646,WRITHE1A.648  
!      WRITE(6,'('' '')')
!      WRITE(6,'('' COMPRESSED FIELD INDEX NO 3'')')
!      WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(145)
*D WRITHE1A.696,WRITHE1A.698  
!      WRITE(6,'('' '')')
!      WRITE(6,'('' LOOKUP TABLE'')')
!      WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(151)*FIXHD(152)
*/
*/       - End of mod TJNOWRIT
*/
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  tropin11.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID TROPIN11
*/
*DECLARE QCTROP1A
*D QCTROP1A.373
c set minimum tropopause pressure to the match level structure
      TRP_MAX=6000.

*DECLARE TROPIN1A
*I TROPIN1A.132
c carried over from original low-res ozone fix - 
c not actually a good idea, overestimates things at high lats.?
c     DTI=MAX_TROP_LEVEL

*B TROPIN1A.362
c with a 0.002 LAPSE_TROP TROPIN seems to overestimate tp by 1 level or so
c c.f qctrop diagnostics (although they are capped too)
c     DO I=1,L2
c       if (IT(I).gt.MIN_TROP_LEVEL) IT(I)=IT(I)-1
c     end do

*DECLARE C_LAPSE
*/ tuned from 0.002 to suit 11 level setup.
*D C_LAPSE.5
      PARAMETER(LAPSE_TROP=0.003) !  TROPOPAUSE LAPSE RATE

*DECLARE RAD_CTL1
*/
*/  Slacken criterion for highest possible tropopause so it works with
*/        the old set of 11 levels.
*/
*D AWI1F402.24
          IF ( AKH(LEVEL)/PREF+BKH(LEVEL) .LT. .07 ) MAX_TROP = LEVEL
*/
*DECLARE CONVEC3C
*/
*/  Truly evil workspace defn & use needs fudge - 38 probably excessive !
*/
*D CONVEC3C.535,536
      REAL WORK(NPNTS,38),    !  WORK SPACE
     *     WORK2(NPNTS,38)
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  tuning_2.mod_4.7.mf77
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID TUNING
*/
*/ as tuning.mod, but remove resetting of psmax=0.4
*/ using 0.6 is better.
*/
*/ 15/3/99. CDJ.
*/
*DECLARE BIOLOGY
*D OJP0F404.139,OJP0F404.143
        IF (K .GT. 6) THEN
C       Convert depth from cm to m by *100.0 on top of
C       this form of the Martin et al. expression...
          remin_rate=100.0 * 8.58/(ZDZ(K)-DZ(K)/2.0)
        ELSE
          remin_rate=0.1
        ENDIF
*/
*/-----------------------------------------------------------
*/
*DECLARE FLUX_CO2
*D OJP0F404.591
        piston_vel(I) = 0.31 * wind_speed * wind_speed
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  vflux_drift.mod_301015
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID VFLUX_CORR
*DECLARE CNTLOCN
*I OJL1F405.80
     &,L_VFLUX_CORR
*I OJL1F405.83
     &,L_VFLUX_CORR
*I OJL1F405.86
     &,L_VFLUX_CORR
*DECLARE TYPPTRO
*I TYPPTRO.66
     &, joc_vflux_mask
     &, joc_vflux_corr
     &, joc_vflux_Sacc
     &, joc_vflux_Aacc
     &, joc_vflux_Cacc
*I TYPPTRO.87
     &  joc_vflux_mask,joc_vflux_corr,
     &  joc_vflux_Sacc,joc_vflux_Aacc,joc_vflux_Cacc,
*DECLARE STOCNPT1
*B STOCNPT1.116
      joc_vflux_mask = SI(400,0,im_index)
      joc_vflux_corr = SI(401,0,im_index)
      joc_vflux_Sacc = SI(402,0,im_index)
      IF (L_OCARBON) THEN
        joc_vflux_Aacc = SI(403,0,im_index)
        joc_vflux_Cacc = SI(404,0,im_index)
      END IF
*DECLARE ROWCALC
*B ROWCALC.843
     &,D1(joc_vflux_corr)
*DECLARE TRACER
*B TRACER.31
     &,vflux_corr
*B RH011293.192
      REAL vflux_corr(IMT*jmt)
*I TRACER.862
c ***APPLY (LAST YEAR's) VFLUX CONSERVATION CORRECTION
      IF (L_VFLUX_CORR) THEN
      do i=1,imt
      do k=1,kmt(i)
        ta(i,k,S_TRACER)=   ta(i,k,S_TRACER)+   
     &                      vflux_corr(imt+S_TRACER)*c2dtts
        IF (L_OCARBON) THEN
        ta(i,k,ALK_TRACER)= ta(i,k,ALK_TRACER)+ 
     &                      vflux_corr(imt+ALK_TRACER)*c2dtts

        ta(i,k,TCO2_TRACER)=ta(i,k,TCO2_TRACER)+
     &                      vflux_corr(imt+TCO2_TRACER)*c2dtts
        END IF
      enddo
      enddo
      END IF
c ********
*DECLARE SFCADD
*I OJP0F404.956
c***include surface water flux adjustment in this calc now
      IF (L_FLUXCORR) THEN
        VTCO2_FLUX(i)=VTCO2_FLUX(i)
     &  + con_salt*fluxcorw(I)*TB(i,1,TCO2_TRACER)/C2DTTS
      END IF

*I OJP0F404.962
      IF (L_FLUXCORR) THEN
        VALK_FLUX(i)=VALK_FLUX(i)
     &  + con_salt*fluxcorw(I)*TB(i,1,ALK_TRACER)/C2DTTS
      END IF

c*********
*DECLARE OCNFRST1
*I OCNFRST1.173
        CALL VFLUX_CORR(
*CALL ARGSIZE
*CALL ARGD1
*CALL ARGSTS
*CALL ARGOINDX
*CALL ARGOCONE
*CALL ARGPTRO
     &           O_FLDDEPC,ICODE,CMESSAGE)
c ********

*I OCNFRST1.203
c
c ***SUBROUTINE VFLUX_CORR
c
c ***Gather fields to work out a correction such that the global virtual tracer
c ***fluxes at the surface (worked out with local tracer values) match the
c ***global water flux. 
c ***Once a year, calculate and store the correction

      SUBROUTINE VFLUX_CORR(
*CALL ARGSIZE
*CALL ARGD1
*CALL ARGSTS
*CALL ARGOINDX
*CALL ARGOCONE
*CALL ARGPTRO
     &           FKMP_G,ICODE,CMESSAGE)

      IMPLICIT NONE

      INTEGER
     &       ICODE        ! Return code : 0 Normal Exit
      CHARACTER*(80)
     &       CMESSAGE     ! Error message if return code >0

*CALL CSUBMODL
*CALL CMAXSIZE
*CALL TYPSIZE
*CALL TYPD1
*CALL TYPPTRO
*CALL TYPSTS
*CALL TYPOINDX
*CALL TYPOCONE

*CALL UMSCALAR
*CALL C_MDI
*CALL CHSUNITS
*CALL CNTLALL
*CALL CNTLOCN
*CALL PARPARM
*CALL PARCOMM
*CALL CTIME
*CALL COMOCFLW
*CALL C_SOILH
*CALL C_PERMA
*CALL OTRACPNT

      INTEGER I,J,K,I_G,INFO,ij,ij_ST,IMT_ST  ! counters. *_ST means -halo
      REAL DAYSYR

c need global, surface avg. values to convert water fluxes to tracer fluxes.
c salinity held in psu/1000
c ALK, TCO2 held in micromoles/litre.Surface avg. values from Sarmiento+Gruber book, WOCE data
      REAL REF_SAL, REF_ALK, REF_TCO2
      PARAMETER (REF_SAL=0.035, REF_ALK=2363., REF_TCO2=2075.)

      INTEGER VF_SNOAR,
     &        VF_ICEAR,
     &        VF_SFLUX,
     &        VF_AFLUX,
     &        VF_CFLUX

      REAL    WATER(imt*jmt), 
     &        S_ACC(imt*jmt),S_ACC_GLOBAL(imt*jmt_global),S_CORR,
     &        A_ACC(imt*jmt),A_ACC_GLOBAL(imt*jmt_global),A_CORR,
     &        C_ACC(imt*jmt),C_ACC_GLOBAL(imt*jmt_global),C_CORR

      REAL    VF_MASK(imt*jmt_global)

      REAL    DYT_G(jmt_global), 
     &        CST_G(jmt_global), 
     &        FKMP_G(imt,jmt_global),
     &        AREA_G, 
     &        VOL_G


c recover fluxes from stash - should have been specified by macro
c
      CALL FIND_PTR_ADD(32228,69,VF_ICEAR,  ! change in seaice H (thermd)
*CALL ARGSIZE
*CALL ARGSTS
     &           2,ICODE,CMESSAGE)

      CALL FIND_PTR_ADD(32229,69,VF_SNOAR,  ! change in snow on ice H (thermd)
*CALL ARGSIZE
*CALL ARGSTS
     &           2,ICODE,CMESSAGE)

      CALL FIND_PTR_ADD(30280,69,VF_SFLUX,    ! PMER salt flux applied
*CALL ARGSIZE
*CALL ARGSTS
     &           2,ICODE,CMESSAGE)

      IF (L_OCARBON) THEN
      CALL FIND_PTR_ADD(30293,69,VF_AFLUX,    ! PMER ALK flux applied
*CALL ARGSIZE
*CALL ARGSTS
     &           2,ICODE,CMESSAGE)

      CALL FIND_PTR_ADD(30292,69,VF_CFLUX,    ! PMER TCO2 flux applied
*CALL ARGSIZE
*CALL ARGSTS
     &           2,ICODE,CMESSAGE)
      END IF

      IF (L_OCYCLIC) THEN
        IMT_ST=IMTm2
      ELSE
        IMT_ST=IMT
      END IF

      DO J=1,jmt
      DO I=1,IMT
        ij=I+(J-1)*IMT
        water(ij)=0.
        s_acc(ij)=0.
        IF (L_OCARBON) THEN
          a_acc(ij)=0.
          c_acc(ij)=0.
        END IF
      end do
      end do
c
c this j counting assumes MPI version (not necessarily >1 PE, but an
c N,S halo per processor. At time of coding, only MPI version works
c anyway)     
c

      IF (MOD(stepim(o_im),2)) THEN

      DO J=2,jmt-1
      DO I=1,IMT_ST
        ij_ST=I+(J-1)*IMT_ST - 1
        ij   =I+(J-1)*IMT    - 1
 
c water flux used - "conserved"
      if (D1(joc_ple+ij).ne.RMDI) 
     &    water(ij+1)=water(ij+1)+D1(joc_ple+ij)

      if (D1(joc_river+ij).ne.RMDI) 
     &    water(ij+1)=water(ij+1)+D1(joc_river+ij)

      if (D1(VF_SNOAR+ij_st).ne.RMDI) 
     &    water(ij+1)=water(ij+1)-D1(VF_SNOAR+ij_st)*RHO_SNOW

      if (D1(VF_ICEAR+ij_st).ne.RMDI) 
     &    water(ij+1)=water(ij+1)-D1(VF_ICEAR+ij_st)*RHO_ICE

      if (L_FLUXCORR) then
        if (D1(joc_anom_salt+ij).ne.RMDI) 
     &        water(ij+1)=water(ij+1)+D1(joc_anom_salt+ij)
      end if

      water(ij+1)=water(ij+1)/RHO_WATER_SI
 
c salt flux actually applied
      if (D1(VF_SFLUX+ij_st).ne.RMDI) s_acc(ij+1)=s_acc(ij+1)+
     &    D1(VF_SFLUX+ij_st)*1e-9

      if (D1(joc_salinc+ij).ne.RMDI) s_acc(ij+1)=s_acc(ij+1)-
     &    D1(joc_salinc+ij)*0.01*DZ(1)

c biogeo fluxes 
       IF (L_OCARBON) THEN
         if (D1(VF_AFLUX+ij_st).ne.RMDI) a_acc(ij+1)=a_acc(ij+1)-
     &     D1(VF_AFLUX+ij_st)*0.01*DZ(1)

         if (D1(VF_CFLUX+ij_st).ne.RMDI) c_acc(ij+1)=c_acc(ij+1)-
     &     D1(VF_CFLUX+ij_st)*0.01*DZ(1)
       END IF

      END DO
      END DO

      DO i=imt,imt*(jmt-1),imt
        d1(joc_vflux_corr+i+NT+1)=d1(joc_vflux_corr+i+NT+1)+c2dtts
      END DO

      DO J=1,jmt
      DO I=1,IMT
        ij=I+(J-1)*IMT - 1
        D1(joc_vflux_Sacc+ij)=D1(joc_vflux_Sacc+ij) +
     &                      (s_acc(ij+1)-(water(ij+1)*REF_SAL) )*c2dtts
       IF (L_OCARBON) THEN
        D1(joc_vflux_Aacc+ij)=D1(joc_vflux_Aacc+ij) +
     &                      (a_acc(ij+1)-(water(ij+1)*REF_ALK) )*c2dtts

        D1(joc_vflux_Cacc+ij)=D1(joc_vflux_Cacc+ij) +
     &                      (c_acc(ij+1)-(water(ij+1)*REF_TCO2) )*c2dtts
       END IF
      END DO
      END DO

      END IF

c ***ONCE A YEAR, WORK OUT NEW RATE ADJUSTMENT
      IF (PREVIOUS_TIME(1).lt.I_YEAR) THEN

c global gather of *_ACC field, metrics
      CALL GATHER_FIELD(D1(joc_vflux_Sacc),s_acc_global,
     &                    imt,jmt,imt,jmt_global,
     &                    0,gc_all_proc_group,ICODE)
 
      IF (L_OCARBON) THEN
      CALL GATHER_FIELD(D1(joc_vflux_Aacc),a_acc_global,
     &                    imt,jmt,imt,jmt_global,
     &                    0,gc_all_proc_group,ICODE)

      CALL GATHER_FIELD(D1(joc_vflux_Cacc),c_acc_global,
     &                    imt,jmt,imt,jmt_global,
     &                    0,gc_all_proc_group,ICODE)
      END IF

      CALL GATHER_FIELD(D1(joc_vflux_mask),vf_mask,
     &                    imt,jmt,imt,jmt_global,
     &                    0,gc_all_proc_group,ICODE)

      CALL O_SMARTPASS(1,1,DYT(J_1),DYT_G
     &                ,jfin-jst+1,jmt_global,jst,2)

      CALL O_SMARTPASS(1,1,CST(J_1),CST_G
     &                ,jfin-jst+1,jmt_global,jst,2)

c work out correction on PE 0
      IF (MYPE.eq.0) THEN
        s_corr=0.
        a_corr=0.
        c_corr=0.
        area_g=0.
        vol_g=0.
        do j=1,jmt_global
        do i=1,imt_st
          i_g=i+(j-1)*imt
          if (vf_mask(i_g).gt.0) then
             s_corr=s_corr+s_acc_global(i_g)*
     &                CST_G(j)*DYT_G(j)*DXT(i)

      IF (L_OCARBON) THEN
             a_corr=a_corr+a_acc_global(i_g)*
     &                CST_G(j)*DYT_G(j)*DXT(i)

             c_corr=c_corr+c_acc_global(i_g)*
     &                CST_G(j)*DYT_G(j)*DXT(i)
      END IF

             area_g=area_g+CST_G(j)*DYT_G(j)*DXT(i)

            IF (FKMP_G(I,J).GT.0) THEN
                vol_g=vol_g+CST_G(j)*DYT_G(j)*DXT(i)*
     &                ZDZ(INT(FKMP_G(I,J)))*0.01
            END IF
          end if
        end do
        end do

c surface
c        s_corr=s_corr/area_g/(DZ(1)*0.01)/d1(joc_vflux_corr+imt+NT+1)
c        IF (L_OCARBON) THEN
c          a_corr=a_corr/area_g/(DZ(1)*0.01)/d1(joc_vflux_corr+imt+NT+1)
c          c_corr=c_corr/area_g/(DZ(1)*0.01)/d1(joc_vflux_corr+imt+NT+1)
c        END IF
c whole depth
         s_corr=s_corr/vol_g/d1(joc_vflux_corr+imt+NT+1)
         IF (L_OCARBON) THEN
           a_corr=a_corr/vol_g/d1(joc_vflux_corr+imt+NT+1)
           c_corr=c_corr/vol_g/d1(joc_vflux_corr+imt+NT+1)
         END IF

      END IF !MYPE

c broadcast correction to other PEs, store it and 0 the accumulation field
      call gc_rbcast(1,1,0,nproc,info,s_corr)
      call gc_rbcast(1,1,0,nproc,info,a_corr)
      call gc_rbcast(1,1,0,nproc,info,c_corr)

      do i=imt,imt*(jmt-1),imt
        d1(joc_vflux_corr+i+NT+1)=0.
        d1(joc_vflux_corr+i+S_TRACER-1)=s_corr
        IF (L_OCARBON) THEN
          d1(joc_vflux_corr+i+ALK_TRACER-1)=a_corr
          d1(joc_vflux_corr+i+TCO2_TRACER-1)=c_corr
        END IF
      end do

      do i=1,imt*jmt
        d1(joc_vflux_Sacc+I-1)=0.
        IF (L_OCARBON) THEN
          d1(joc_vflux_Aacc+I-1)=0.
          d1(joc_vflux_Cacc+I-1)=0.
        END IF
      end do

      write(6,'(a)')"****************** VFLUX ****************"
      write(6,'(a)')"Doing vflux correction"
      write(6,'(a,e12.5)')"Salinity correction   ",s_corr
       IF (L_OCARBON) THEN
      write(6,'(a,e12.5)')"Alkalinity correction ",a_corr
      write(6,'(a,e12.5)')"TCO2 correction       ",c_corr
       END IF
      write(6,'(a)')"*****************************************"

      END IF !LAST_TIMESTEP_OF_YEAR

      RETURN

      END SUBROUTINE VFLUX_CORR

      SUBROUTINE FIND_PTR_ADD (STASHCODE,STASHMACRO_TAG,ADDRESS,
C  
C from ~jeff/um/hurrikan/mods/oasis3.0_um4.5_v1.mod
C
*CALL ARGSIZE
*CALL ARGSTS
     &  INTERNAL_MODEL,ICODE,CMESSAGE)
C
      IMPLICIT NONE
C
*CALL CSUBMODL
*CALL TYPSIZE
*CALL TYPSTS
C
      INTEGER STASHCODE         ! IN  - STASH code
      INTEGER STASHMACRO_TAG    ! IN  - STASHmacro tag number
      INTEGER ADDRESS           ! OUT - Address in D1
      INTEGER INTERNAL_MODEL    ! IN  - internal_model id.
      INTEGER ICODE             ! OUT - Error return code
      CHARACTER*(*) CMESSAGE    ! OUT - Error return message
C
*CALL C_MDI
C
C  Subroutines called
C
      EXTERNAL FINDPTR
C
C     Local variables
C
      INTEGER
     &  SECTION,                ! STASH section number
     &  ITEM,                   ! STASH item number
     &  PROCESS_CODE,           ! processing code
     &  FREQ_CODE,              ! frequency code
     &  START,END,PERIOD,       ! start, end and period step
     &  GRIDPT_CODE,WEIGHT_CODE,! gridpt and weighting codes
     &  BOTTOM_LEVEL,TOP_LEVEL, ! bottom and top input level
     &  GRID_N,GRID_S,GRID_W,GRID_E  ! grid corner definitions
C
      PROCESS_CODE=IMDI
      FREQ_CODE=IMDI
      START=IMDI
      END=IMDI
      PERIOD=IMDI
      GRIDPT_CODE=IMDI
      WEIGHT_CODE=IMDI
      BOTTOM_LEVEL=IMDI
      TOP_LEVEL=IMDI
      GRID_N=IMDI
      GRID_S=IMDI
      GRID_E=IMDI
      GRID_W=IMDI
      SECTION = STASHCODE / 1000
      ITEM = STASHCODE - SECTION * 1000
C
      CALL FINDPTR(INTERNAL_MODEL,SECTION,ITEM,
     &             PROCESS_CODE,FREQ_CODE,START,END,PERIOD,
     &             GRIDPT_CODE,WEIGHT_CODE,
     &             BOTTOM_LEVEL,TOP_LEVEL,GRID_N,GRID_S,GRID_W,GRID_E,
     &             STASHMACRO_TAG,IMDI,ADDRESS,
*CALL ARGSIZE
*CALL ARGSTS
     &    ICODE,CMESSAGE)
C
      RETURN
      END
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/Orginal mod:  stash80.mod
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
*ID STASH80
*DECLARE VERSION
*D GSS1F400.1166
      CHARACTER*80 STASH_SET     !Names of stasets files
*DECLARE INACTR1
*D GSS1F400.960
      CHARACTER*80 DIR
*D GNF0F401.5
      CALL FORT_GET_ENV('STASETS_DIR',11,DIR,80,err)
*D INACTR1.237
          DO J=1,80
