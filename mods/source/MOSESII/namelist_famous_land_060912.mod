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
