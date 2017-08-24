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
CLL  SUBROUTINE  LWPTSC
CLL
CLL    PURPOSE
CLL  It calculates scaled pathlengths of each gaseous absorber for each
CLL  layer and returns them in DPATH for use by LWMAST, which sums them
CLL  to get the total scaled pathlengths between each pair of layers so
CLL  that the gaseous transmissivities can be calculated.
CLL  Used in version 1A (gaseous effects treated as Slingo & Wilderspin,
CLL  1986) of the UM LW code.
CLL  If UPDATE *DEF CRAY is off, a version is produced which except
CLL  for the addition of ! comments is standard FORTRAN 77 (and which
CLL  sets the "vector length" to 1) but the standard version includes
CLL  CRAY automatic arrays also.
CLL          11/1/91 - modified so that the radiative
CLL  effectiveness of a given amount of CO2 or O3 contains a constant
CLL  term as well as the pressure-to-the-power-alpha term.  The latter,
CLL  based on the important absorption being in pressure-broadened
CLL  wings of strong lines, is all that is needed in the troposphere,
CLL  but at low pressures and low pathlengths temperature-broadened
CLL  wings and weak lines also contribute significantly.  The value of
CLL  the constants are based on single-column tests - altered on
CLL  14/5/91.
CLL
CLL                      Author: William Ingram 24 Oct 1990
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL   4.2    Sept.96  T3E migration: *DEF CRAY removed;
CLL                   *DEF T3E used for T3E library functions;
CLL                   dynamic allocation no longer *DEF controlled;
CLL                   cray HF functions replaced by T3E lib functions.
CLL                       S.J.Swarbrick
CLL
CLL  It conforms to standard A of UMDP 4 (version 2, 18/1/90), and
CLL  includes no 8X-deprecated features.  ,
CLL
CLL  It is part of component P232 (longwave radiation) which is in task
CLL  P23 (radiation).
CLL
CLL  External documentation is in UMDP 23.
C*L
      SUBROUTINE LWPTSC (H2O, CO2, O3, PSTAR, AC, BC, AB, BB, TAC,
     &     L2,                                                          
     &     NWET, NOZONE, NLEVS, L1,                            DPATH)
C*
      INTEGER NGASES
C Effective number of absorbing gases treated in the longwave
      PARAMETER (NGASES=4)    !  Standard set is water vapour line and
C continuum (which have to be treated separately because the pathlength
C scaling is different), ozone and carbon dioxide.
C*L
      INTEGER!, INTENT (IN) ::
     &     L2,                       ! Number of points to be treated   
     &     NWET,                     ! Number of levels with moisture -
C     ! above these zero is used for the continuum and a small constant
C     ! H2OLMN for line absorption (where zero would give trouble)
     &     NOZONE,                   ! Number of levels with ozone data
C     ! provided - below these the value in the lowest of them is used
     &     NLEVS,                    ! Number of levels
     &     L1                        ! First dimension of input arrays
C     !  (The different physical assumptions about water vapour and
C     !  ozone in levels where no data is provided means that separate
C     !  loops are used for levels with and without water vapour but
C     !  only the indexing needs changing for levels with and without
C     !  their own ozone data.)
      REAL!, INTENT(IN) ::
     &     H2O(L1,NWET), CO2,        ! Mass mixing ratio (mK in UMDP 23)
     &     O3(L1,NOZONE),            !             of each absorbing gas
     &     TAC(L1,NLEVS),            ! Mid-layer temperatures
     &     PSTAR(L1),                ! Surface pressure
     &     AC(NLEVS), BC(NLEVS),     ! A & B for layer centres and
     &     AB(NLEVS+1), BB(NLEVS+1)  !                       boundaries
      REAL!, INTENT(OUT) ::
     &     DPATH(L2,NGASES,NLEVS)
C     !  The scaled pathlengths are returned in DPATH, indexed by NGASES
C     !  1 is CO2, 2 is H2O line absorption, 3 is O3, 4 is H2O continuum
CL    !  LWPTSC has no EXTERNAL calls and no significant structure
CL    !  WORK is the only dynamically allocated array                   
      REAL WORK(L2,2,2)
C     !  WORK is used to hold powers of layer boundary pressures used
C     !  in 2.3.1 and passed from one level to the next to save
C     !  re-calculation.  (This does prevent autotasking over levels.)
      REAL S9,                       ! Pressure scaling normalization
     &     S4,                       !  constants, for alpha = 0.9 & 0.4
     &     DOPCO2,                   ! Constants allowing for Doppler
     &     DOPO3,                    !      broadening for CO2 and ozone
     &     CONCON1,                  ! Constants used to find the
     &     CONCON2,                  !  continuum pathlengths (Eq 2.3.8)
     &     CO2CON                    ! Constant and exponent used in
      REAL POWER,                    !      temperature scaling for CO2
     &     CO2PSC,                   ! Pressure-scaled CO2 pathlength
     &     DAB,                      ! Differences of A and B across
     &     DBB,                      !                   current layer
     &     PTOP,                     ! Pressure at top of curent layer
     &     DP,                       !        and its pressure thickness
     &     PH                        ! Pressure scaling with alpha=0.9,
C                                    !              used for H2O and CO2
      INTEGER LEVEL, J,              ! Loopers over levels & points
     &     ONETWO,                   ! Flipper
     &     OLEVEL                    ! Index for the ozone data to be
C                                    !        used in the current level
C*
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

      PARAMETER ( CONCON1 = 1.66 / (G*101300.) )
      PARAMETER ( DOPCO2 = 1.3E-3, DOPO3 = 3.E-2 )
CL    Note that these constants are related to the pressure p1 at which
CL    the "Lorentz" and "Doppler" terms are equally important:
CL
CL                                 (p1/p0) ** ALPHA
CL                       DOPxx  = ------------------
CL                                        G
CL
      REAL H2OLMN                    ! Minimum pathlength for H2O line
      PARAMETER ( H2OLMN = 1.E-10 )  !                        absorption
C     !  FORTRAN 77 will not allow the following constants to be
C     !  defined in a PARAMETER statement, but the CRAY compiler will
C     !  give the same effect as if they were.
      S4 = 101325.**(-0.4) / (G*1.4)
      S9 = 101325.**(-0.9) / (G*1.9)
      CONCON2 = EXP(-1800./296.) / (0.005*EPSILON)
      CO2CON = 2.0 / LOG(10.)
C
      DO 1 J=1, L2
       WORK(J,1,1) = ( PSTAR(J) ) ** 1.9
       WORK(J,1,2) = ( PSTAR(J) ) ** 1.4
C      !  These lines could of course have ( PSTAR(J) * BB(1) + AB(1) )
    1 CONTINUE
      ONETWO=1
C
      DO 2 LEVEL=1, NWET
       DAB = AB(LEVEL) - AB(LEVEL+1)
       DBB = BB(LEVEL) - BB(LEVEL+1)
       OLEVEL = MAX (1, LEVEL+NOZONE-NLEVS)
       DO 20 J=1, L2
        PTOP = PSTAR(J) * BB(LEVEL+1) + AB(LEVEL+1)
        DP = DAB + PSTAR(J) * DBB
C       !  First, the pressure scaling common to CO2 and H2O line:
        WORK(J,3-ONETWO,1) = PTOP ** 1.9
        PH = S9 * ( WORK(J,ONETWO,1) - WORK(J,3-ONETWO,1) )
C       !  CO2 has temperature scaling also: (2.3.1)-(2.3.3):
        CO2PSC = ( PH + DP * DOPCO2 ) * CO2
        POWER = 6.5 + CO2CON * LOG(CO2PSC)
        IF  ( POWER .LT. 0.0 )  POWER = 0.0
        DPATH(J,1,LEVEL) = CO2PSC * ( TAC(J,LEVEL) / 263. ) ** POWER
C       !  For H2O line absorption just apply (2.3.1):
        DPATH(J,2,LEVEL) = PH * H2O(J,LEVEL)
C       !  but re-set zero humidities to a small value:
        IF ( DPATH(J,2,LEVEL) .EQ. 0. )  DPATH(J,2,LEVEL) = H2OLMN
C       !  For ozone (2.3.1) only:
        WORK(J,3-ONETWO,2) = PTOP ** 1.4
        DPATH(J,3,LEVEL) = O3(J,OLEVEL) *
     &  ( S4 * ( WORK(J,ONETWO,2) - WORK(J,3-ONETWO,2) ) + DP * DOPO3 )
C       !  For the H2O continuum apply (2.3.8):
        DPATH(J,4,LEVEL) = CONCON1 * H2O(J,LEVEL) *
     &  ( PSTAR(J) * BC(LEVEL) + AC(LEVEL) ) * ( PSTAR(J) * DBB + DAB )
     &  * ( 1. + CONCON2 * H2O(J,LEVEL) * EXP ( 1800. / TAC(J,LEVEL) ) )
   20  CONTINUE
       ONETWO = 3 - ONETWO            !  Flip ONETWO so that the numbers
C      !   we are avoiding recalculating do not have to be copied either
    2 CONTINUE
C
C     !  Above the levels where moisture is calculated, put in constant
C     !  for the H2O pathlength but treat CO2 and O3 the same:
C
      DO 3 LEVEL=NWET+1, NLEVS
       DAB = AB(LEVEL) - AB(LEVEL+1)
       DBB = BB(LEVEL) - BB(LEVEL+1)
       OLEVEL = MAX (1, LEVEL+NOZONE-NLEVS)
       DO 30 J=1, L2
        PTOP = PSTAR(J) * BB(LEVEL+1) + AB(LEVEL+1)
        DP = DAB + PSTAR(J) * DBB
C       !  CO2 is scaled just as before:
        WORK(J,3-ONETWO,1) = PTOP ** 1.9
        PH = S9 * ( WORK(J,ONETWO,1) - WORK(J,3-ONETWO,1) )
        CO2PSC = ( PH + DP * DOPCO2 ) * CO2
        POWER = 6.5 + CO2CON * LOG(CO2PSC)
        IF  ( POWER .LT. 0.0 )  POWER = 0.0
        DPATH(J,1,LEVEL) = CO2PSC * ( TAC(J,LEVEL) / 263. ) ** POWER
        WORK(J,3-ONETWO,2) = PTOP ** 1.4
C       !  and so is O3:
        DPATH(J,3,LEVEL) = O3(J,OLEVEL) *
     &  ( S4 * ( WORK(J,ONETWO,2) - WORK(J,3-ONETWO,2) ) + DP * DOPO3 )
C       !  For H2O line absorption just put in a small constant:
        DPATH(J,2,LEVEL) = H2OLMN
C       !  For the H2O continuum can use zero since this goes into EXP
C       !  rather than LOG:
        DPATH(J,4,LEVEL) = 0.
   30  CONTINUE
       ONETWO = 3 - ONETWO            !  Flip ONETWO as before
    3 CONTINUE
C
      RETURN
      END
