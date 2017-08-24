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
CLL Subroutine SOLPOS   ----------------------------------------------
CLL
CLL Purpose :
CLL  Calculations of the earth's orbit described in the first page of
CLL  the "Calculation of incoming insolation" section of UMDP 23, i.e.
CLL  from the day of the year (and, in forecast mode, whether it is a
CLL  leap year) and the orbital "constants" (which vary over
CLL  "Milankovitch" timescales) it calculates the sin of the solar
CLL  declination and the inverse-square scaling factor for the solar
CLL  "constant".  It is thus intrinsically scalar.  The FORTRAN code
CLL  present depends on whether *DEF CAL360 is set during UPDATE: this
CLL  replaces the Julian calendar with the climate-mode 360-day calendar
CLL
CLL   Author:    William Ingram
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL   3.4    20/06/94 DEF CAL360 replaced by LOGICAL LCAL360;
CLL                   PARAMETER statements duplicated for 360 and
CLL                   365 day calendar.
CLL                                                S.J.Swarbrick
!LL   4.4    27/02/97 Testing for leap years modified to deal with
!LL                   no leap every 100y except for every 400y
!LL                   Author: M.Gallani
CLL
CLL Programming standard :
CLL    Written in FORTRAN 77, with the addition of "!" comments and
CLL  underscores in variable names.
CLL    Written to comply with 12/9/89 version of UMDP 4 (meteorological
CLL  standard).
CLL
CLL Logical components covered : P233
CLL
CLL Project task :
CLL
CLL External documentation: P23
CLL
CLLEND -----------------------------------------------------------------
C*L
      SUBROUTINE SOLPOS (DAY, YEAR, SINDEC, SCS, LCAL360)
C

      LOGICAL LCAL360    !In, true if 360 day calendar in use.
C
      INTEGER!, INTENT(IN) ::
     &     DAY,                            !  Day-number in the year
     &     YEAR                            !  Calendar year
      REAL!, INTENT(OUT) ::
     &     SINDEC,                         !  sin(solar declination)
     &     SCS                             !  solar constant scaling
C*                                                            factor
CL This routine has no dynamically allocated work areas and no
CL  significant structure.  It calls the intrinsic functions FLOAT, SIN
CL  & COS, but no user functions or subroutines.
CL
      REAL GAMMA, E, TAU0, SINOBL,         ! Basic orbital constants
     &     TAU1_360, TAU1_365, E1,E2,E3,E4,!Derived orbital contants
     &     TWOPI                           ! 2pi
      REAL DINY_360, DINY_365              ! Number of days in year
      REAL M, V                            ! Mean & true anomaly
C*L------------------COMDECK C_PI---------------------------------------
CLL
CLL 4.0 19/09/95  New value for PI. Old value incorrect
CLL               from 12th decimal place. D. Robinson
CLL
      REAL PI,PI_OVER_180,RECIP_PI_OVER_180

      PARAMETER(
     & PI=3.14159265358979323846, ! Pi
     & PI_OVER_180 =PI/180.0,   ! Conversion factor degrees to radians
     & RECIP_PI_OVER_180 = 180.0/PI ! Conversion factor radians to
     &                              ! degrees
     & )
C*----------------------------------------------------------------------

      PARAMETER ( TWOPI = 2. * PI )
      PARAMETER (GAMMA=1.352631, E=.0167,  ! Gamma, e
     &     TAU0 = 2.5,                     ! True date of perihelion
     &     SINOBL = .397789 )              ! Sin (obliquity)
      PARAMETER ( E1 = E * (2.-.25*E*E),
     &     E2 = 1.25 * E*E,                ! Coefficients for 3.1.2
     &     E3 = E*E*E * 13./12.,
     &     E4=( (1.+E*E*.5)/(1.-E*E) )**2 )! Constant for 3.1.4
      PARAMETER (DINY_360=360., TAU1_360=TAU0*DINY_360/365.25+0.71+.5)
C
      PARAMETER (TAU1_365=TAU0+.5)
C
      IF (.NOT. LCAL360) THEN
        IF (mod(year,4) .eq. 0 .AND.          ! is this a leap year?
     &    (mod(year,400) .eq. 0 .OR. mod(year,100) .ne. 0)) then
          DINY_365 = 366.
        ELSE
          DINY_365 = 365.
        END IF
      END IF
C
!     In forecast mode and in climate mode with real-year means, DINY
!     depends on whether it is a leap year, otherwise DINY_36x = 360.
!
C  TAU1 is modified so as to include the conversion of day-ordinal into
C  fractional-number-of-days-into-the-year-at-12-Z-on-this-day.
C
      IF (LCAL360) THEN
        M = TWOPI * (FLOAT(DAY)-TAU1_360) / DINY_360          ! Eq 3.1.1
      ELSE
        M = TWOPI * (FLOAT(DAY)-TAU1_365) / DINY_365          ! Eq 3.1.1
      END IF
      V = M + E1*SIN(M) + E2*SIN(2.*M) + E3*SIN(3.*M)         ! Eq 3.1.2
      SCS = E4 * ( 1. + E * COS(V) ) **2                      ! Eq 3.1.4
      SINDEC = SINOBL * SIN (V - GAMMA)                       ! Eq 3.1.6
      RETURN
      END
