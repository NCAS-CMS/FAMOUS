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
CLL  Subroutine : DAY2CHAR  -------------------------------------------
CLL
CLL  Purpose: Convert days to the character to represent the period.
CLL
CLL  Tested under compiler:   cft77
CLL  Tested under OS version: UNICOS 6.1.5A
CLL
CLL  Author:   R A Stratton
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
CLL
CLL  Logical components covered: S51
CLL
CLL  Project task: S51
CLL
CLL  External documentation: UM documentation paper 7 - Filenaming
CLL                          conventions for the Unified Model
CLL
CLLEND -----------------------------------------------------------------
C*L  Interface and arguments: ------------------------------------------
C
      SUBROUTINE DAY2CHAR(NDAYS,DAYCHAR)
C
      IMPLICIT NONE
C
      INTEGER       NDAYS        ! IN  - number of days in period
      CHARACTER*1   DAYCHAR      ! OUT - character for period
C
C*----------------------------------------------------------------------
C  Common blocks
C
C
C  Local variables
C
C
C  IF period a multiple of years
C
        IF (MOD(NDAYS,360).EQ.0)THEN
         IF (NDAYS.EQ.360) THEN        ! 1 year mean
           DAYCHAR='y'
         ELSE IF (NDAYS.EQ.1800) THEN  ! 5 year means
           DAYCHAR='v'
         ELSE IF (NDAYS.EQ.3600) THEN  ! 10 year means
           DAYCHAR='x'
         ELSE IF (NDAYS.EQ.18000) THEN ! 50 year means
           DAYCHAR='l'
         ELSE IF (NDAYS.EQ.36000) THEN ! 100 year means
           DAYCHAR='u'
         ELSE IF (NDAYS.EQ.360000) THEN ! 1000 year means
           DAYCHAR='z'
         ELSE                           ! not a special period
           DAYCHAR='0'
         ENDIF
c
        ELSE
C periods less than one year
C
         IF (NDAYS.EQ.5) THEN          ! 5 days means
           DAYCHAR='p'
         ELSE IF (NDAYS.EQ.7) THEN     ! weekly means
           DAYCHAR='w'
         ELSE IF (NDAYS.EQ.10) THEN    ! 10 day means
           DAYCHAR='t'
         ELSE IF (NDAYS.EQ.14) THEN    ! fortnightly means
           DAYCHAR='r'
         ELSE IF (NDAYS.EQ.30) THEN    ! monthly means
           DAYCHAR='m'
         ELSE IF (NDAYS.EQ.90) THEN    ! seasonal means
           DAYCHAR='s'
         ELSE                          ! not a special period
           DAYCHAR='0'
         ENDIF
        ENDIF
C
      RETURN
CL----------------------------------------------------------------------
      END
