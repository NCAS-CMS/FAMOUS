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
CLL Subroutine ZONM     ----------------------------------------------
CLL
CLL Purpose :
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL vn4.2  07/11/96:Allow this routine to be used in C93_2A (Farnon)
CLL
CLL Programming standard :
CLL
CLL Logical components covered :
CLL
CLL Project task :
CLL
CLL External documentation:
CLL
CLLEND -----------------------------------------------------------------
C
      SUBROUTINE ZONM(VAR,ZVAR,MASK,MASS,LGPTS,NY,NX)
C
      IMPLICIT NONE
C
      INTEGER
     1    NX,                        !     - Longitude: p,u-rows
     2    NY,                        !     - Latitude: rowlength
     4    I,
     5    J
      REAL
     1    VAR(NY,NX),                ! IN  - Variable for calculation
     2    ZVAR(NX),                  ! IN  - Zonal mean of variable
     3    MASK(NY,NX),               ! IN  - Mask (eg. Land/sea)
     3    MASS(NY,NX)                ! IN  - Mass weighting (p,u-grid)
      LOGICAL
     1    LGPTS(NX)                  !     - No of land points/row
CL
C---
C    Subroutines called
C       NONE
C---
C
C    Local Variables
C
      REAL
     1    SUMZTOP(NX),
     2    SUMZBOT(NX)
CL
C
CL 1. SUMZTOP= Sum of(MASS * MASK * VAR) for each I
CL    SUMZBOT= Sum of(MASS * VAR) for each I
CL    ZVAR= SUMZTOP/SUMZBOT
CL
      DO 10,I=1,NX
        SUMZTOP(I)=0.0
        SUMZBOT(I)=0.0
          DO 20,J=1,NY
            SUMZBOT(I) = SUMZBOT(I) + MASK(J,I) * MASS(J,I)
            SUMZTOP(I) = SUMZTOP(I) + (MASK(J,I) * MASS(J,I) * VAR(J,I))
20    CONTINUE
        IF (LGPTS(I)) THEN            ! If logical pts/row=true then
          ZVAR(I) = SUMZTOP(I)/SUMZBOT(I) !  Calculate zonal mean
        ELSE                          ! Else
          ZVAR(I)=0.0                 !  zonal mean=zero for that row
        ENDIF
10    CONTINUE
      RETURN
      END
C
CL======================================================================
CL 5.1 Subroutine to calculate the column mean of a variable
CL----------------------------------------------------------------------
