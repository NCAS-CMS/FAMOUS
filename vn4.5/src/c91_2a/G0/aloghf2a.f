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
CLL Function ALOGHF
CLL
CLL Purpose:  Portable version of Cray library function to compute
CLL           the natural logarithm of ARG1.
CLL
CLL Tested under compiler:   fort77
CLL Tested under OS version: HP-UX A.08.07
CLL
CLL  Model            Modification history:
CLL version  Date
CLL  3.2   16/07/93   New deck. Tracey Smith.
CLL  3.3   22/09/93   Improved comments  Tracey Smith
CLL  4.0   23/03/95   Write out input value if <0. Tracey Smith
CLL
CLL Programming Standard: UM Doc Paper 3, version 5 (08/12/92)
CLL
      REAL FUNCTION ALOGHF(ARG1)
        IMPLICIT NONE
        REAL ARG1
        IF(ARG1.LE.0) THEN
          WRITE(6,*) 'Error: Attempt to LOG negative value ',ARG1
          CALL ABORT
        ELSE
          ALOGHF=ALOG(ARG1)
        END IF
      END
