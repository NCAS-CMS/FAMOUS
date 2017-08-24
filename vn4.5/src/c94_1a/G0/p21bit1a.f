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
CLL  Function:  P21BITS ------------------------------------------------
CLL
CLL  Purpose: Returns the number of bits used for exponent in packing
CLL           and unpacking Cray 64 bit words to/from 32 bit words by
CLL           Cray routines PACK21 and EXPAND21.  As from vn2.6 this
CLL           becomes a function of dump/PP/ancillary file format
CLL           version number.
CLL           For the present P21BITS=6 whatever version number.
CLL           This may be reconsidered in future.
CLL
CLL  Author:   T.Johns            Date:           08 May 1992
CLL
CLL  Tested under compiler:   cf77
CLL  Tested under OS version: UNICOS 6.1.5a
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered: S4
CLL
CLL  Project task: S420
CLL
CLL  External documentation:
CLL    UMDP F3
CLL
C*L  Interface and arguments: ------------------------------------------
C
      INTEGER FUNCTION P21BITS(version)
C
      IMPLICIT NONE
C
      INTEGER
     &    version  ! IN  UM dump/pp/ancillary file format version
C
C* ---------------------------------------------------------------------
C
      P21BITS=6
C
      RETURN
      END
