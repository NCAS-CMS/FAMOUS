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
CLL  Subroutine SOLANG -----------------------------------------------
CLL
CLL Purpose :
CLL  Calculations of the earth's orbit described in the second page of
CLL  the "Calculation of incoming insolation" section of UMDP 23, i.e.
CLL  from the sin of the solar  declination, the position of each point
CLL  and the time limits it calculates how much sunlight, if any, it
CLL  receives.
CLL
CLL    Author:    William Ingram
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL Programming standard :
CLL    Written to comply with 12/9/89 version of UMDP 4 (meteorological
CLL  standard).
CLL    Written in FORTRAN 77 with the addition of "!" comments and
CLL  underscores in variable names.
CLL
CLL Logical components covered : P233
CLL
CLL Project task :
CLL
CLL External documentation: UMDP 23
CLL
CLLEND -----------------------------------------------------------------
C*L
      SUBROUTINE SOLANG (SINDEC, T, DT, SINLAT, LONGIT, K,
     &     LIT, COSZ)
      INTEGER!, INTENT(IN) ::
     &     K                          ! Number of points
      REAL!, INTENT(IN) ::
     &     SINDEC,                    ! Sin(solar declination)
     &     T, DT,                     ! Start time (GMT) & timestep
     &     SINLAT(K),                 ! sin(latitude) & longitude
     &     LONGIT(K)                  ! of each point
      REAL!, INTENT(OUT) ::
     &     LIT(K),                    ! Sunlit fraction of the timestep
     &     COSZ(K)                    ! Mean cos(solar zenith angle)
C                                     ! during the sunlit fraction
C*
CL This routine has no dynamically allocated work areas.  It calls the
CL intrinsic functions SQRT, ACOS & SIN, but no user functions or
CL subroutines.  The only structure is a loop over all the points to be
CL dealt with, with IF blocks nested inside to cover the various
CL possibilities.
      INTEGER J                       ! Loop counter over points
      REAL TWOPI,                     ! 2*pi
     &     S2R                        ! Seconds-to-radians converter
      REAL SINSIN,            ! Products of the sines and of the cosines
     &     COSCOS,            ! of solar declination and of latitude.
     &     HLD,               ! Half-length of the day in radians (equal
     &                        ! to the hour-angle of sunset, and minus
     &     COSHLD,            ! the hour-angle of sunrise) & its cosine.
     &     HAT,               ! Local hour angle at the start time.
     &     OMEGAB,            ! Beginning and end of the timestep and
     &     OMEGAE,            ! of the period over which cosz is
     &     OMEGA1,            ! integrated, and sunset - all measured in
     &     OMEGA2,            ! radians after local sunrise, not from
     &     OMEGAS,            ! local noon as the true hour angle is.
     &     DIFSIN,            ! A difference-of-sines intermediate value
     &     DIFTIM,            ! and the corresponding time period
     &     TRAD, DTRAD
C     ! These are the start-time and length of the timestep (T & DT)
C     ! converted to radians after midday GMT, or equivalently, hour
C     ! angle of the sun on the Greenwich meridian.
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

      PARAMETER ( TWOPI = 2. * PI, S2R = PI / 43200.)
C
      TRAD = T * S2R - PI
      DTRAD = DT * S2R
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO 100 J = 1, K                          ! Loop over points
       HLD = 0.                                ! Logically unnecessary
C statement without which the CRAY compiler will not vectorize this code
       SINSIN = SINDEC * SINLAT(J)
       COSCOS = SQRT( (1.-SINDEC**2) * (1.-SINLAT(J)**2) )
       COSHLD = SINSIN / COSCOS
       IF (COSHLD.LT.-1.) THEN                 ! Perpetual night
          LIT(J) = 0.
          COSZ(J) = 0.
        ELSE
          HAT = LONGIT(J) + TRAD               ! (3.2.2)
          IF (COSHLD.GT.1.) THEN               !   Perpetual day - hour
             OMEGA1 = HAT                      ! angles for (3.2.3) are
             OMEGA2 = HAT + DTRAD              ! start & end of timestep
           ELSE                                !   At this latitude some
C points are sunlit, some not.  Different ones need different treatment.
             HLD = ACOS(-COSHLD)               ! (3.2.4)
C The logic seems simplest if one takes all "times" - actually hour
C angles - relative to sunrise (or sunset), but they must be kept in the
C range 0 to 2pi for the tests on their orders to work.
             OMEGAB = HAT + HLD
             IF (OMEGAB.LT.0.)   OMEGAB = OMEGAB + TWOPI
             IF (OMEGAB.GE.TWOPI) OMEGAB = OMEGAB - TWOPI
             IF (OMEGAB.GE.TWOPI) OMEGAB = OMEGAB - TWOPI
C            !  Line repeated - otherwise could have failure if
C            !  longitudes W are > pi rather than < 0.
             OMEGAE = OMEGAB + DTRAD
             IF (OMEGAE.GT.TWOPI) OMEGAE = OMEGAE - TWOPI
             OMEGAS = 2. * HLD
C Now that the start-time, end-time and sunset are set in terms of hour
C angle, can set the two hour-angles for (3.2.3).  The simple cases are
C start-to-end-of-timestep, start-to-sunset, sunrise-to-end and sunrise-
C -to-sunset, but two other cases exist and need special treatment.
             IF (OMEGAB.LE.OMEGAS .OR. OMEGAB.LT.OMEGAE) THEN
                OMEGA1 = OMEGAB - HLD
              ELSE
                OMEGA1 = - HLD
             ENDIF
             IF (OMEGAE.LE.OMEGAS) THEN
                OMEGA2 = OMEGAE - HLD
              ELSE
                OMEGA2 = OMEGAS - HLD
             ENDIF
             IF (OMEGAE.GT.OMEGAB.AND.OMEGAB.GT.OMEGAS) OMEGA2=OMEGA1
C  Put in an arbitrary marker for the case when the sun does not rise
C  during the timestep (though it is up elsewhere at this latitude).
C  (Cannot set COSZ & LIT within the ELSE ( COSHLD < 1 ) block
C  because 3.2.3 is done outside this block.)
          ENDIF           ! This finishes the ELSE (perpetual day) block
          DIFSIN = SIN(OMEGA2) - SIN(OMEGA1)             ! Begin (3.2.3)
          DIFTIM = OMEGA2 - OMEGA1
C Next, deal with the case where the sun sets and then rises again
C within the timestep.  There the integration has actually been done
C backwards over the night, and the resulting negative DIFSIN and DIFTIM
C must be combined with positive values representing the whole of the
C timestep to get the right answer, which combines contributions from
C the two separate daylit periods.  A simple analytic expression for the
C total sun throughout the day is used.  (This could of course be used
C alone at points where the sun rises and then sets within the timestep)
          IF (DIFTIM.LT.0.) THEN
            DIFSIN = DIFSIN + 2. * SQRT(1.-COSHLD**2)
            DIFTIM = DIFTIM + 2. * HLD
          ENDIF
          IF (DIFTIM.EQ.0.) THEN
C Pick up the arbitrary marker for night points at a partly-lit latitude
             COSZ(J) = 0.
             LIT(J) = 0.
           ELSE
             COSZ(J) = DIFSIN*COSCOS/DIFTIM + SINSIN     ! (3.2.3)
             LIT(J) = DIFTIM / DTRAD
          ENDIF
       ENDIF            ! This finishes the ELSE (perpetual night) block
  100 CONTINUE
      RETURN
      END
