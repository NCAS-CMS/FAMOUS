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
CLL  Routine: COLM  ----------------------------------------------------
CLL
CLL  Purpose: Service routine to calculate weighted column mean of a
CLL           3D field, as required in zonal mean subroutine.
CLL
CLL  Tested under compiler:   cf77
CLL  Tested under OS version: UNICOS 5.1
CLL
CLL  Author:  T.C.Johns
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL vn4.2  07/11/96:Allow this routine to be used in C93_2A (Farnon)
CLL
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered: D61
CLL
CLL  Project task: C61
CLL
CLL  External documentation:
CLL    Unified Model Doc Paper C61 - Printed diagnostics
CLL
CLL  -------------------------------------------------------------------
C*L  Interface and arguments: ------------------------------------------
C
      SUBROUTINE COLM(VAR,COLMN,P_MASS,NX,NY,NZ)
C
      IMPLICIT NONE
C
      INTEGER
     1    NX,                        !     - Longitude: p,u-rows
     2    NY,                        !     - Latitude: rowlength
     3    NZ,                        !     - Level:p,u-level
     4    I,
     5    J,
     6    K
      REAL
     1    VAR(NY,NX,NZ),             ! IN  - Variable for calculation
     2    COLMN(NY,NX),              ! OUT - Column mean of variable
     3    P_MASS(NY,NX,NZ)           ! IN  - Mass weighting (p,u-grid)
CL
C---
C    Subroutines called
C       NONE
C---
C
C    Local Variables
C
      REAL
     1    SUMCBOT(NY,NX)
C
CL  The calculation is self-explanatory
C
      DO 100 I=1,NX
        DO 110 J=1,NY
          SUMCBOT(J,I)=0.0
          COLMN(J,I)  =0.0
 110    CONTINUE
 100  CONTINUE
C
      DO 200 K=1,NZ
        DO 210 I=1,NX
          DO 220 J=1,NY
            COLMN(J,I)   = COLMN(J,I) + P_MASS(J,I,K)*VAR(J,I,K)
            SUMCBOT(J,I) = SUMCBOT(J,I) + P_MASS(J,I,K)
 220      CONTINUE
 210    CONTINUE
 200  CONTINUE
C
      DO 300 I=1,NX
        DO 310 J=1,NY
          COLMN(J,I) = COLMN(J,I)/SUMCBOT(J,I)
 310    CONTINUE
 300  CONTINUE
C
      RETURN
      END
