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
C
      SUBROUTINE GLBM(VAR,GVAR,START,END,AREA,MASK,GPTS,NX,NY,LOGIP)
C
C *** Subroutine calculates global or quarter globe means
C
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL vn4.2  07/11/96:Allow this routine to be used in C93_2A (Farnon)
CLL
CLL Logical components covered: D61
CLL
CLLEND --------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER
     1    NX,                    !     - Longitude:point on p,u-rows
     2    NY,                    !     - Longitude:pont on rowlength
     5    I,
     6    J,
     7    START,                 !     - Marker for start of 1/4 globe
     8    END                    !     - Marker for end of 1/4 globe
      REAL
     1    VAR(NY,NX),            ! IN  - Variable for calculation
     2    GVAR,                  ! IN  - Global mean of variable
     3    AREA(NY,NX),           ! IN  - Area weighting
     4    MASK(NY,NX)            ! IN  - Mask or mass weighting
      LOGICAL
     1    GPTS,                  !     - No of land points/row
     2    LOGIP                  !     - True if p-grid, false if u-grid
CL
C---
C    Subroutines called
C       NONE
C---
C
C    Local Variables
C
      REAL
     1    SUMGTOP,
     2    SUMGBOT
CL
C
CL 1. SUMGTOP= Sum of(AREA * MASK * VAR) for all pts
CL    SUMGBOT= Sum of(AREA * MASK) for all pts
CL    GVAR=SUMGTOP/SUMGBOT
CL
      SUMGTOP=0.0
      SUMGBOT=0.0
C
CL    If its a p grid, then multiply start and end rows by 0.5, so that
CL    rows are not counted twice when computing 1/4 globe means
C
      IF (LOGIP) THEN
        DO 10,I=START,END
          DO 20,J=1,NY
            IF ((I .NE. START) .AND. (I .NE. END)) THEN
               SUMGBOT = SUMGBOT + AREA(J,I) * MASK(J,I)
               SUMGTOP = SUMGTOP + (AREA(J,I) * MASK(J,I) * VAR(J,I))
            ELSE
               SUMGBOT = SUMGBOT + 0.5*(AREA(J,I) * MASK(J,I))
               SUMGTOP = SUMGTOP + 0.5*(AREA(J,I)*MASK(J,I)*VAR(J,I))
            END IF
20        CONTINUE
10      CONTINUE
      ELSE
        DO 30,I=START,END
          DO 40,J=1,NY
            SUMGBOT = SUMGBOT + AREA(J,I) * MASK(J,I)
            SUMGTOP = SUMGTOP + (AREA(J,I) * MASK(J,I) * VAR(J,I))
40        CONTINUE
30      CONTINUE
      END IF
C
      IF (GPTS) THEN               ! If logical pts/globe or 1/4 globe
        GVAR = SUMGTOP / SUMGBOT   !  Calculate global mean
      ELSE                         ! Else
        GVAR=0.0                   !  Set global mean=0
      ENDIF
      RETURN
      END
