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
CLL Subroutine SWCLOP   ----------------------------------------------
CLL
CLL Purpose :
CLL  It calculates cloud optical properties for use by SWMAST.
CLL  It is suitable for single-column use.
CLL
CLL         Author: William Ingram
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   4.0   4/7/95    Code added for safety & accuracy near the singular
CLL  case when the e-folding depths for direct & diffuse light are equal
CLL   4.2  22/10/96  The code for the singular case when the e-folding
CLL                  depths for direct & diffuse light are nearly equal
CLL                  rearranged to a form algebraically the same but
CLL                  avoiding numerical problems for very thick cloud
CLL                  which gave failure on the T3E (also slightly more
CLL                  efficient).                           W Ingram
CLL
CLL
CLL Programming standard :
CLL  Contains ! comments, but otherwise conforms to the FORTRAN 77
CLL  It conforms to programming standard A of UMDP 4 (version 2,
CLL  18/1/90).
CLL  standard with no features deprecated by 8X.
CLL
CLL Logical components covered : P234
CLL  (interaction of shortwave radiation with the atmosphere)
CLL
CLL Project task : P23 (radiation)
CLL
CLL External documentation:   UMDP 23, sub-section "cloud optical
CLL  properties" of shortwave section.
CLL
CLLEND -----------------------------------------------------------------
C*L
      SUBROUTINE SWCLOP (CW, RE, COSZ,
     &     L1, L2, NCLDS,                       REF, CTR)
      INTEGER NBANDS         ! Number of spectral bands in the shortwave
      PARAMETER (NBANDS=4)   ! This run uses the standard set of
C shortwave bands as described by Slingo (May 1989, J.Atmos.Sci., 46,
C 10, 1419-1427) or UMDP 23.
      INTEGER!, INTENT (IN)
     &     NCLDS,                      ! Number of clouds
     &     L1,                         ! First dimension of input arrays
     &     L2                          ! Number of points to be treated
      REAL!, INTENT (IN)
     &     CW(L1,NCLDS),               ! Cloud condensed water path
     &     RE(L1,NCLDS),               ! Effective cloud droplet radius
     &     COSZ(L1)                    ! Cos(solar zenith angle)
      REAL!, INTENT (OUT)
     &     REF(L2,NBANDS,NCLDS,2),     ! Cloud reflectivity and
     &     CTR(L2,NBANDS,NCLDS,2)      ! transmissivity for direct and
C                                      ! diffuse radiation respectively
C
C     !  SWCLOP has no dynamically allocated workspace, no EXTERNAL
C     !  calls and no significant structure - just nested loops.
C*
      REAL B0CON, U1                   ! Used in the quadrature
      PARAMETER (B0CON=.375, U1=2.)
C     !  These values are set for PIFM as defined by Zdunkowski & al
C     !  (1980, Beitraege, 53, 147-166).
      REAL SWCTOL                      ! Minimum value of D2 for which
C     !  the standard formulae for direct-beam optical properties will
C     !  be used - if D2 is too small, the singular solution for equal
C     !  direct & diffuse e-folding depths is used.
C     !  ("Too small" rather than "zero", because near the singular case
C     !  the general formulae become ill-conditioned and only the
C     !  asymptotic ones give sensible answers.)
      PARAMETER ( SWCTOL = .005 )
C     !  Surprisingly, the errors are not machine-dependent.
      REAL AFIT(NBANDS), BFIT(NBANDS), !  Used in (4.8.23)-(4.8.25)
     &     CFIT(NBANDS), DFIT(NBANDS), !  but called just a-f
     &     EFIT(NBANDS), FFIT(NBANDS)  !  in UMDP 23.
      REAL TAU, OMEGA, G, BETA0,       !  Intermediate quantities
     &     BETAMU, F, U2, ALPHA1,      !  as defined in UMDP 23
     &     ALPHA2, ALPHA3, ALPHA4, EPS,!
     &     M, E, H, A1, D1, D2, GAMMA1,!
     &     GAMMA2, RT, RM, C1N, P1, P2,!
     &     OM1MF, E2                   !  OMEGA*(1-F) & E**2
      INTEGER CLOUD, J, BAND           !  Loopers over clouds, points &
C                                                                bands
C     !  These are the values for the standard set of 4 bands.
C     !  The limits on Re they imply are that it must be between
C     !  0.344 microns (where OMEGA0 in band 1 becomes greater than 1)
C     !  and 37.6 microns (where G in band 3 becomes greater than 1).
      DATA AFIT  /  28.17,     26.82,     22.64,    12.81 /
      DATA BFIT  /  1.305E-3,  1.346E-3,  1.454E-3, 1.641E-3 /
      DATA CFIT  /  -56.17E-9, -6.938E-6, 463.6E-6, 0.2011 /
      DATA DFIT  /  0.1632,    23.49,     1238.,    7556. /
      DATA EFIT  /  0.8287,    0.7935,    0.7535,   0.8255 /
      DATA FFIT  /  2482.,     4226.,     6560.,    4353. /
C
      DO 100 CLOUD=1, NCLDS
       DO 100 BAND=1, NBANDS
Cfpp$  Select(CONCUR)
       DO 100 J=1, L2
C       !  First calculate basic optical parameters.
C       !  N.B.  The code does not check for OMEGA0 < 0 or G > 1 (which
C       !  will occur if Re is too large) or OMEGA0 > 1 (which will
C       !  occur if Re is too small), and these unphysical cases will
C       !  produce a negative argument for the square root.
        IF ( RE(J,CLOUD) .GT. 0. ) THEN ! (4.8.23) & (4.8.24)
           TAU = CW(J,CLOUD) * ( AFIT(BAND) + BFIT(BAND)/RE(J,CLOUD) )
           OMEGA = 1. - CFIT(BAND) - DFIT(BAND) * RE(J,CLOUD)
         ELSE                           ! RE may be unset where no cloud
           TAU = 1.
           OMEGA = .5
        ENDIF
C       !  The special-case code ends up with 0/0 if tau v.small
C       !  N.B. 32bit machines may need larger small value
      IF ( TAU .LE. 1.0E-4 ) THEN !
           CTR(J,BAND,CLOUD,1) = 1.
           CTR(J,BAND,CLOUD,2) = 1.
           REF(J,BAND,CLOUD,1) = 0.
           REF(J,BAND,CLOUD,2) = 0.
         ELSE
        G = EFIT(BAND) + FFIT(BAND) * RE(J,CLOUD)            ! (4.8.25)
C       !  Now apply the chosen quadrature - Zdunkovski's PIFM
        BETA0 = B0CON * (1.-G)                               ! (4.8.18)
        BETAMU = 0.5 - 0.75 * COSZ(J) * G/(1.+G)             ! (4.8.19)
        F  = G * G                                           ! (4.8.20)
        U2 = 2.                                              ! (4.8.22)
C       !  Now have quadrature-dept quantities.  (U1 is constant)
C       !  Next find alphas etc:
        ALPHA1 = U1 * ( 1. - OMEGA*(1.-BETA0) )              ! (4.8.14)
        ALPHA2 = U2 * OMEGA * BETA0                          ! (4.8.15)
        OM1MF  = OMEGA * (1.-F)
        ALPHA3 = OM1MF * BETAMU                              ! (4.8.16)
        ALPHA4 = OM1MF - ALPHA3                              ! (4.8.17)
        EPS = SQRT ( ALPHA1**2 - ALPHA2**2 )                 ! (4.8.13)
        M = ALPHA2 / (ALPHA1+EPS)                            ! (4.8.8)
        E = EXP ( - EPS * TAU )                              ! (4.8.7)
        E2 = E * E
        H  = 1. - OMEGA * F                                  ! (4.8.12)
C       ! Now can find a1, the contribution to the transmissivity for a
C       ! direct incident beam due to transmission as a direct beam.
        A1 = EXP ( - H * TAU / COSZ(J) )                     ! (4.8.1)
C       !  Now find transmissivities & reflectivities for diffuse
C       !  incident beam, a4 & a5:
        D1 = 1. - E2 * M * M                                 ! (4.8.6)
        CTR(J,BAND,CLOUD,2)  = E * (1.-M*M) / D1             ! (4.8.4)
        REF(J,BAND,CLOUD,2) = M * (1.-E2) / D1               ! (4.8.5)
C       ! Next find the gammas:
        D2 = ( COSZ(J)*EPS )**2 - H*H                        ! (4.8.11)
        P1 = COSZ(J) * ( ALPHA1*ALPHA3 + ALPHA2*ALPHA4 ) - H*ALPHA3
        P2 = COSZ(J) * ( ALPHA1*ALPHA4 + ALPHA2*ALPHA3 ) + H*ALPHA4
        CTR(J,BAND,CLOUD,1) = .5   ! Ensure valid REALs set for 2nd IF
        REF(J,BAND,CLOUD,1) = .5   ! block to test, even if 1st skipped.
        IF ( D2 .NE. 0. ) THEN
C         ! For direct incident beam, (a1+a2) & a3 (eqns 4.8.2-4.8.5):
          GAMMA1 = P1 / D2
          GAMMA2 = P2 / D2
          CTR(J,BAND,CLOUD,1)  = A1 * ( 1.-REF(J,BAND,CLOUD,2)*GAMMA1 )
     &                         - GAMMA2 * ( CTR(J,BAND,CLOUD,2) - A1 )
          REF(J,BAND,CLOUD,1) = GAMMA1 * ( 1.-CTR(J,BAND,CLOUD,2)*A1 )
     &                          - GAMMA2 * REF(J,BAND,CLOUD,2)
        ENDIF
        IF ( ABS(D2) .LT. SWCTOL .OR.
     & CTR(J,BAND,CLOUD,1) .GT. 1. .OR. CTR(J,BAND,CLOUD,1) .LT. 0. .OR.
     & REF(J,BAND,CLOUD,1) .GT. 1. .OR. REF(J,BAND,CLOUD,1) .LT. 0.)THEN
C         !  With the standard 4 bands and rE, this can only be
C         !   executed in band 4 (only a 30th of the total insolation.)
          GAMMA1 = .5 * P1 / EPS
          GAMMA2 = .5 * P2 / EPS
          RT = ALPHA2 * TAU / ( 1. + ALPHA3 / GAMMA1 )
          RM = M * ( 1. + M * RT ) / ( RT + M )
          C1N = GAMMA1 * TAU / ( RM - 1. )
          CTR(J,BAND,CLOUD,1)  = A1 +
     &      E * ( C1N * M * (1.-E2) + TAU * GAMMA2 ) / COSZ(J)
          REF(J,BAND,CLOUD,1) = C1N * ( E2 - RM ) / COSZ(J)
        ENDIF
        ENDIF
  100 CONTINUE
      RETURN
      END
