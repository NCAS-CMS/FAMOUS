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
!  Routine: GRIB_TIME_INT ------------------------------------------
!
! Description :Given 2 dates and times find the time interval
!              between them and the appropriate time unit to use for
!              encoding grib.
!
!    Author:   R A Stratton    Reviewer : D Goddard
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  3.4    12/10/94   Original code. R A Stratton
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: ?
! System Task: ?
!
! --------------------------------------------------------------------
!
      SUBROUTINE GRIB_TIME_INT (I_YR1,I_MON1,I_DAY1,I_HR1,I_MIN1,
     &                     I_YR2,I_MON2,I_DAY2,I_HR2,I_MIN2,
     &                     CAL360,TIME_INTERVAL,T_UNITS)
!
      IMPLICIT NONE

      INTEGER
     &     I_YR1,                 ! IN  - model time (years)
     &     I_MON1,                ! IN  - model time (months)
     &     I_DAY1,                ! IN  - model time (days)     T1
     &     I_HR1,                 ! IN  - model time (hours)
     &     I_MIN1,                ! IN  - model time (minutes)
     &     I_YR2,                 ! IN  - model time (years)
     &     I_MON2,                ! IN  - model time (months)
     &     I_DAY2,                ! IN  - model time (days)     T2
     &     I_HR2,                 ! IN  - model time (hours)
     &     I_MIN2,                ! IN  - model time (minutes)
     &     TIME_INTERVAL,         ! OUT - dt =t2-t1
     &     T_UNITS                ! OUT - 1 hours,2 days,3 months
!                                         4 years
      LOGICAL CAL360              ! IN - true if 360 day year

!*------------------------------------------------------------------
!
!  Common blocks
!
CLL  Comdeck: CDAYDATA -------------------------------------------------
CLL
CLL  Purpose: Constants needed by routines to calculate day/month/year
CLL           from incremental day number since calendar zero point,
CLL           and vice-versa, when using Gregorian calendar
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL  Logical components covered: S62
CLL
CLLEND ------------------------------------------------------------
      INTEGER DAYS_PER_4C,DAYS_PER_C,DAYS_PER_4Y,DAYS_PER_Y,
     1        DAYS_IN_MONTH(12),DAYS_TO_MONTH(12)
      DATA    DAYS_PER_4C,DAYS_PER_C,DAYS_PER_4Y,DAYS_PER_Y
     1      /    146097  ,   36524  ,    1461   ,   365   /
      DATA    DAYS_IN_MONTH
     1      /  31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
      DATA    DAYS_TO_MONTH
     1      /   0, 31, 59, 90,120,151,181,212,243,273,304,334 /
C
!
!  Local variables
!
      INTEGER
     &   YEAR1,YEAR2       ! years
     &,  DAY1,DAY2         ! number of days since calendar zero
     &,  MINS_PER_MONTH    ! number of minutes in a month
     &,  MINS_PER_YEAR     ! number of minutes in a year
     &,  MINS_PER_DECADE   ! number of minutes in a decade
     &,  MINS_PER_CENTURY  ! number of minutes in a century
     &,  TIME_INTMIN

!------------------------------------------------------------------
! 1. Add up days from time zero to specified time
!
!
! 1.1 30-day month (360 day year) calendar
!
      IF (CAL360) THEN
        DAY1 = 360*I_YR1+30*(I_MON1-1)+I_DAY1-1
        DAY2 = 360*I_YR2+30*(I_MON2-1)+I_DAY2-1
      ELSE
!
! 1.2 Gregorian calendar   - not tested available for completeness
!                            and possible future use.
!
!     If leap year and after 28 February, adjust day number by one
        IF (MOD(I_YR1,4).EQ.0   .AND.
     &      (MOD(I_YR1,400).EQ.0 .OR. MOD(I_YR1,100).NE.0) .AND.
     &                                            I_MON1.GT.2) THEN
          DAY1 = I_DAY1+1
        ELSE
          DAY1 = I_DAY1
        ENDIF
        IF (MOD(I_YR2,4).EQ.0   .AND.
     &      (MOD(I_YR2,400).EQ.0 .OR. MOD(I_YR2,100).NE.0) .AND.
     &                                           I_MON2.GT.2) THEN
          DAY2 = I_DAY2+1
        ELSE
          DAY2 = I_DAY2
        ENDIF
!     Add on days in the preceding months in the year
        DAY1 = DAY2 + DAYS_TO_MONTH(I_MON1) - 1
        DAY2 = DAY1 + DAYS_TO_MONTH(I_MON2) - 1
        YEAR1= I_YR1 - 1
        YEAR2= I_YR2 - 1
!     Add on days up to the specified year
        DAY1 =  DAY1+(YEAR1/400)*DAYS_PER_4C
        DAY2 =  DAY2+(YEAR2/400)*DAYS_PER_4C
        YEAR1 = YEAR1-(YEAR1/400)*400
        YEAR2 = YEAR2-(YEAR2/400)*400
        DAY1=  DAY1+(YEAR1/100)*DAYS_PER_C
        DAY2=  DAY2+(YEAR2/100)*DAYS_PER_C
        YEAR1 = YEAR1-(YEAR1/100)*100
        YEAR2 = YEAR2-(YEAR2/100)*100
        DAY1 =  DAY1+(YEAR1/4)*DAYS_PER_4Y
        DAY2 =  DAY2+(YEAR2/4)*DAYS_PER_4Y
        YEAR1 = YEAR1-(YEAR1/4)*4
        YEAR2 = YEAR2-(YEAR2/4)*4
        DAY1 =  DAY1+YEAR1*DAYS_PER_Y
        DAY2 =  DAY2+YEAR2*DAYS_PER_Y
      ENDIF
!--------------------------------------------------------------------
! 2. Convert days, hours and minutes to minutes since calendar zero,
!     and subtract times to get interval in minutes
!
      MINS_PER_MONTH=30*1440
      IF (CAL360)THEN
        MINS_PER_YEAR=360*MINS_PER_MONTH
        MINS_PER_DECADE=MINS_PER_YEAR*10
        MINS_PER_CENTURY=MINS_PER_YEAR*100
      ELSE
! Warning needs to be corrected if leap years
! Note not used at present as PP2GRIB not set up for meaning
! with real calendar years
!
        MINS_PER_YEAR=365*MINS_PER_MONTH ! problem if leap year
        MINS_PER_DECADE=MINS_PER_YEAR*10
        MINS_PER_CENTURY=MINS_PER_YEAR*100
      ENDIF
      TIME_INTMIN=1440*(DAY2-DAY1)+60*(I_HR2-I_HR1)+I_MIN2-I_MIN1


      IF (MOD(TIME_INTMIN,60).NE.0) THEN  ! assumes multiple of mins
        T_UNITS=0
        TIME_INTERVAL=TIME_INTMIN
      ELSE
       IF (MOD(TIME_INTMIN,1440).NE.0) THEN
        T_UNITS=1             ! HOURS
        TIME_INTERVAL=TIME_INTMIN/60
       ELSE
        IF (MOD(TIME_INTMIN,MINS_PER_MONTH).NE.0) THEN
          T_UNITS=2             ! DAYs
          TIME_INTERVAL=TIME_INTMIN/1440
        ELSE
          IF (MOD(TIME_INTMIN,MINS_PER_YEAR).NE.0) THEN
           T_UNITS=3             ! months
           TIME_INTERVAL=TIME_INTMIN/MINS_PER_MONTH
          ELSE
            IF (MOD(TIME_INTMIN,MINS_PER_DECADE).NE.0) THEN
              T_UNITS=4             ! years
              TIME_INTERVAL=TIME_INTMIN/MINS_PER_YEAR
            ELSE
              IF (MOD(TIME_INTMIN,MINS_PER_CENTURY).NE.0) THEN
                T_UNITS=5             ! decades
                TIME_INTERVAL=TIME_INTMIN/MINS_PER_DECADE
              ELSE
                T_UNITS=7             ! century
                TIME_INTERVAL=TIME_INTMIN/MINS_PER_CENTURY
! Note could be extended to cope with normals ?
!          T_UNITS=6       ! normal 30 years
!
!  Note cannot cope with intervals greater than 255 * time unit
!  ie 25500 years
!
              ENDIF
            ENDIF
          ENDIF
        ENDIF
       ENDIF
      ENDIF

      RETURN
!-----------------------------------------------------------------
      END
