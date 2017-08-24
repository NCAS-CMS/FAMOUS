C ******************************COPYRIGHT******************************
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
C
!LL  SUBROUTINE SETPERLEN---------------------------------------------
!LL
!LL  Purpose:
!LL           Return length of current meaning period using mean-level
!LL           (0, 1, 2 or 3) & current date. Days_in_month is declared
!LL           and provided by comdeck CDAYDATA.
!LL
!LL  Method: where this routine is called at the end of a month, the
!LL          month will already have been incremented, which is why
!LL          i_month-1 is used in setting the period length. Every
!LL          100 years is not leap unless year is divisible by 400.
!LL
!  Current Code Owner: Mark Gallani
!
!  History:
!  Version  Date     Comment
!  =======  ====     =======
!  4.4      15/01/97 Original code. (Mark Gallani)
!
!  Code description:
!   FORTRAN 77 + common extensions also in Fortran 90
!LL
!LL
!LL Programming standard :  UMDP 3 Version 7  (rev. 6/10/94)
!LL
!LL Logical components covered :
!LL
!LL Project task :
!LL
!LL External documentation:
!LL
!LL-----------------------------------------------------------------
!*L Arguments:------------------------------------------------------
      SUBROUTINE SETPERLEN (MEANLEV,I_MONTH,I_YEAR,PERIODLEN)

      IMPLICIT NONE

      INTEGER
     &  MEANLEV,         ! IN - Mean level indicator, e.g. 0, 1, 2 or 3
     &  I_MONTH,         ! IN - model time (month)
     &  I_YEAR,          ! IN - model time (year)
     &  PERIODLEN        ! OUT - length of current meaning period (days)

!-----------------------------------------------------------------------
! Workspace usage:------------------------------------------------------
! NONE
!-----------------------------------------------------------------------
! External subroutines called:------------------------------------------
! NONE
!*----------------------------------------------------------------------
! Define local variables:-----------------------------------------------

      LOGICAL L_LEAP     ! Leap year indicator

!-----------------------------------------------------------------------
!  Comdecks:
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
!-----------------------------------------------------------------------
! End of standard header info

      IF (mod(i_year,4) .eq. 0 .AND.          ! is this a leap year?
     &    (mod(i_year,400) .eq. 0 .OR. mod(i_year,100) .ne. 0)) then
        L_LEAP = .TRUE.
      ELSE
        L_LEAP = .FALSE.
      END IF

      IF (meanlev .eq. 0) then       ! instantaneous data (e.g. part-
        periodlen = 1                ! way through a monthly mean)
      ELSEIF (meanlev .eq. 1) then   ! end of monthly mean
        IF (L_LEAP .AND. (i_month .eq. 3)) then  ! Is it leap year Feb?
          periodlen = days_in_month(i_month-1) + 1
        ELSE IF (i_month .eq. 1) then  ! end of Dec, so can't use
          periodlen = 31               ! days_in_month(i_month-1)
        ELSE
          periodlen = days_in_month(i_month-1)
        END IF
      ELSE IF (meanlev .eq. 2) then  ! seasonal mean

!      find season length using current month as pointer
        IF (L_LEAP) then  ! do leap year seasons
          IF (i_month .eq. 5) then  ! season=FebMarApr
            periodlen = 90
          ELSE IF ((i_month .eq. 3) .or. (i_month .eq. 4) .or.
     &            (i_month .eq. 7) .or. (i_month .eq. 12)) then
            periodlen = 91          ! for DJF, JFM, AMJ or SON
          ELSE
            periodlen = 92
          END IF
        ELSE              ! do non-leap year seasons
          IF (i_month .eq. 5) then  ! season=FebMarApr
            periodlen = 89
          ELSE IF ((i_month .eq. 3) .or. (i_month .eq. 4)) then
            periodlen = 90  ! for DJF and JFM
          ELSE IF ((i_month .eq. 7) .or. (i_month .eq. 12)) then
            periodlen = 91  ! for AMJ and SON
          ELSE
            periodlen = 92  ! for all other seasons
          END IF
        END IF  ! end of IF test of L_LEAP, and end of seasons.
      ELSE IF (meanlev .eq. 3) then  ! annual mean

! Bear in mind period 3 may be 366 days if _previous_ year was leap, and
! may not always be 366 days even if current year is a leap year, since
! annual means are often not for calendar years
        IF (L_LEAP .AND. (i_month .ne. 2)) then
          periodlen = 366
        ELSE IF (mod(i_year-1,4) .eq. 0 .AND.      ! was last year leap?
     &    (mod(i_year-1,400) .eq. 0 .OR. mod(i_year-1,100) .ne. 0)) then
          IF (i_month .eq. 2) then
            periodlen = 366
          ENDIF
        ELSE
          periodlen = 365
        ENDIF
      ELSE                      ! meanlev has unexpected value
        periodlen = 1           ! so set weighting factor=1
        write(6,*)'SETPERLEN: MEANLEV not in allowed range of 0 to 3'
      END IF  ! end of IF tests on meanlev

      RETURN
      END

