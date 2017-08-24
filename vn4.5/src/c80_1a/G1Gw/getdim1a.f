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
CLL  Subroutine GET_DIM
CLL
CLL  Written by D Robinson 28/7/92
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
CLL
CLL  Project Task : P3
CLL
CLL  Purpose : Get dimensions of data set components
CLL            from fixed length header.
CLL
CLLEND------------------------------------------------------------------
      SUBROUTINE GET_DIM (FIXHD,LEN_FIXHD,
     +                    LEN_INTHD,LEN_REALHD,
     +                    LEN1_LEVDEPC,LEN2_LEVDEPC,
     +                    LEN1_ROWDEPC,LEN2_ROWDEPC,
     +                    LEN1_COLDEPC,LEN2_COLDEPC,
     +                    LEN1_FLDDEPC,LEN2_FLDDEPC,
     +                    LEN_EXTCNST,LEN_DUMPHIST,
     +                    LEN_CFI1,LEN_CFI2,LEN_CFI3,
     +                    LEN1_LOOKUP,LEN2_LOOKUP,
     +                    LEN_DATA)

CL----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER          !  Dimension of :-
     +  LEN_FIXHD      !  Fixed length header
     + ,LEN_INTHD      !  Integer header
     + ,LEN_REALHD     !  Real header
     + ,LEN1_LEVDEPC   !  Level dependent constants (1st)
     + ,LEN2_LEVDEPC   !  Level dependent constants (2nd)
     + ,LEN1_ROWDEPC   !  Rows  dependent constants (1st)
     + ,LEN2_ROWDEPC   !  Rows  dependent constants (2nd)
     + ,LEN1_COLDEPC   !  Col   dependent constants (1st)
     + ,LEN2_COLDEPC   !  Col   dependent constants (2nd)
     + ,LEN1_FLDDEPC   !  Field dependent constants (1st)
     + ,LEN2_FLDDEPC   !  Field dependent constants (2nd)
     + ,LEN_EXTCNST    !  Extra constants
     + ,LEN_DUMPHIST   !  Dump history
     + ,LEN_CFI1       !  Compressed field index 1
     + ,LEN_CFI2       !  Compressed field index 2
     + ,LEN_CFI3       !  Compressed field index 3
     + ,LEN1_LOOKUP    !  Look up table (1st)
     + ,LEN2_LOOKUP    !  Look up table (2nd)
     + ,LEN_DATA       !  Data section

      INTEGER  FIXHD(LEN_FIXHD)   ! IN  Fixed length header

      LEN_INTHD    = FIXHD(101)
      LEN_REALHD   = FIXHD(106)
      LEN1_LEVDEPC = FIXHD(111)
      LEN2_LEVDEPC = FIXHD(112)
      LEN1_ROWDEPC = FIXHD(116)
      LEN2_ROWDEPC = FIXHD(117)
      LEN1_COLDEPC = FIXHD(121)
      LEN2_COLDEPC = FIXHD(122)
      LEN1_FLDDEPC = FIXHD(126)
      LEN2_FLDDEPC = FIXHD(127)
      LEN_EXTCNST  = FIXHD(131)
      LEN_DUMPHIST = FIXHD(136)
      LEN_CFI1     = FIXHD(141)
      LEN_CFI2     = FIXHD(143)
      LEN_CFI3     = FIXHD(145)
      LEN1_LOOKUP  = FIXHD(151)
      LEN2_LOOKUP  = FIXHD(152)
      LEN_DATA     = FIXHD(161)

      RETURN
      END
