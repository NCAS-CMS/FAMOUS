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
! Dummy version of routine INIT_VEG
      SUBROUTINE INIT_VEG(A_STEP,
C All sizes
C ======================== COMDECK ARGOCPAR ============================
C ========================= COMDECK ARGOCBAS ===========================
     * NT,IMT,JMT,KM,
C ===================== END OF COMDECK ARGOCBAS ========================
C
C ===================== END OF COMDECK ARGOCPAR ========================
C
! History:                                              
! Version  Date    Comment
!  4.2   11/10/96  Enable atmos-ocean coupling for MPP.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins                               
     &        D1_ADDR,D1,LD1,ID1, ! IN/OUT:Addressing of D1 & D1 array
     &        D1_A, D1_O,        ! IN/OUT: Copies of D1 (atmos/ocean)

C Dump headers
     &A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,
     &A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,
! PP lookup headers and Atmos stash array + index with lengths
     &A_LOOKUP,
     &A_MPP_LOOKUP,
     &a_ixsts, a_spsts,
! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  05/10/94  Add murk and user ancillary pointers. RTHBarnes
!  4.1  04/12/95  Add pointers JSTHU and JSTHF. J.Smith
!  4.1  26/04/96  Add pointers for Sulphur Cycle variables (12)  MJW
!  4.3   18/3/97  And for HadCM2 sulphate loading patterns.  Will Ingram
!  4.4   05/8/97  And for Conv. cloud amt on model levs. Julie Gregory
!  4.5  04/03/98   Remove pointer SO2_HILEM (add to CARGPT_ATMOS)
!                  Add 1 pointers for NH3 in S Cycle
!                  Add 3 pointers for Soot              M. Woodage
!  4.5  08/05/98   Add 16 new pointers for User Anc.    D. Goddard
!  4.5  13/05/98   Add pointer for RHcrit variable.     S. Cusack
!  4.5  15/07/98   Add pointers for new 3D CO2 array.   C.D.Jones
!  4.5  17/08/98   Remove JSOIL_FLDS and JVEG_FLDS      D. Robinson
C Argument list.
C Pointers for ATMOSPHERE model variables. Configuration dependent.
C
C Addresses in D1 array of primary variables and 'extra' space
C  variable (Exner pressures)
C  Array  variables (dimensions are resolution dependent)
     &       JU, JV, JTHETA, JQ, JQCL, JQCF, J_DEEP_SOIL_TEMP,  JSMCL,
     &       JOZONE, JTRACER, JP_EXNER,
     &       JSO4, JH2SO4, JSOOT, JMURK, JMURK_SOURCE,
     &       JUSER_MULT1, JUSER_MULT2, JUSER_MULT3, JUSER_MULT4,
     &       JUSER_MULT5, JUSER_MULT6, JUSER_MULT7, JUSER_MULT8,
     &       JUSER_MULT9, JUSER_MULT10, JUSER_MULT11, JUSER_MULT12,
     &       JUSER_MULT13, JUSER_MULT14, JUSER_MULT15, JUSER_MULT16,
     &       JUSER_MULT17, JUSER_MULT18, JUSER_MULT19, JUSER_MULT20,
     &       JSTHU, JSTHF,
     &       JSO2,JDMS,JSO4_AITKEN,JSO4_ACCU,JSO4_DISS,JH2O2,
     &       JSO2_NATEM,JOH,JHO2,JH2O2_LIMIT,JO3_CHEM,
     &       JHadCM2_SO4,JCCA,JRHC,JNH3,
     &       JSOOT_NEW,JSOOT_AGD,JSOOT_CLD,JCO2,
! History:
! Version  Date    Comment
!  3.4   18/05/94  Add new field sin_u_latitude. J F Thomson.
C Constants for dynamics routines
     &    AK_TO_THE_KAPPA, BK_TO_THE_KAPPA, AKH, BKH, AKH_TO_THE_KAPPA,
     &    BKH_TO_THE_KAPPA, COS_U_LATITUDE, COS_P_LATITUDE,
     &    SEC_U_LATITUDE, SEC_P_LATITUDE, TAN_U_LATITUDE, SIN_LONGITUDE,
     &    COS_LONGITUDE, TRUE_LONGITUDE, F1, F2, F3, F3_P,
     &    TRIGS, TWO_D_GRID_CORRECTION, ETA_MATRIX_INV,
C Constants for physics routines
     &    SOILB, LAND_LIST,
     &    SIN_U_LATITUDE,
     &              ICODE,CMESSAGE)

      INTEGER
     & A_STEP             ! IN Current timestep in atmosphere model

      WRITE (6,*) '**ERROR**: INIT_VEG has been called but is'
      WRITE (6,*) 'unavailable.  Either set L_VEG_FRACS to .FALSE. or'
      WRITE (6,*) 'select section A19_1A or A19_2A.'

      CALL ABORT('Run aborted in INIT_VEG')

      RETURN
      END
