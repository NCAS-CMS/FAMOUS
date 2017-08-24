C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1996, METEOROLOGICAL OFFICE, All Rights Reserved.
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
CLL  Subroutine FTSA -------------------------------------------------
CLL
CLL  It calculates (true) surface albedos for P234.
CLL    This version (2B) can modify the surface albedo for the direct
CLL  beam to mimic the effect of a layer of anthropogenic sulphate
CLL  aerosol using the first-order approximation following Charlock
CLL  & al (1991), as part of "HADCM2 physics" (v Mitchell & al 1995).
CLL  Otherwise it matches the standard routine.
CLL                                            William Ingram 19/11/96
CLL  Suitable for single column model use.
CLL
CLL   Author: William Ingram
CLL
CLL  Model            Modification history:
CLL version  Date
CLL   4.2   19/11/96     Written William Ingram, reviewed Cath Senior.
CLL   4.3    18/3/97     Make sulphate calculations optional.  WJI
CLL
CLL  It conforms to programming standard A of UMDP 4, version 2.
CLL  It contains ! comments, but otherwise conforms to the FORTRAN 77
CLL  standard with no features deprecated by 8X.
CLL
CLL Logical components covered : P233
CLL   (ancillary calculations for the shortwave scheme)
CLL
CLL Project task : P23
CLL
CLL  Offline documentation is in UMDP 23, sections "True surface albedo
CLL  specification" and "Modifications to the radiation scheme to
CLL  accommodate the leads model"
CLLEND
C*L
      SUBROUTINE FTSA (
     &     LAND, AICE, TSTAR, SFA, MDSA, COSZ, S, SULPH,
     &     ALPHAC, ALPHAM, DTICE,
     &     SANAON, NLALBS, NSULPAT,
     &     L1, L2,
     &     SALI, SAOS, SANA)
!
      INTEGER !, INTENT(IN) ::
     &     L1,                            ! Full field dimension
     &     L2                             ! Number of points to treat
      LOGICAL !, INTENT(IN) ::
     &     LAND(L1)                       ! Land-sea mask (land .TRUE.)
     &     , SANAON                       ! Is SANA to be output ?
      REAL !, INTENT(IN) ::
     &     AICE(L1),                      ! Sea-ice fraction
     &     TSTAR(L1),                     ! Surface temperature
     &     SFA(L1),                       ! Snow-free surface albedo
     &     MDSA(L1),                      ! Cold deep-snow albedo
     &     ! (These two are alpha sub 0 & alpha sub S resp. in UMDP 23.)
     &     COSZ(L1),                      ! cos(solar zenith angle)
     &     S(L1),                         ! Snow amount (mass/area)
     &     SULPH(L1)                      ! sulphate loading pattern
      REAL!, INTENT(OUT)
     &     SALI(L1,NLALBS),               ! Surface Albedos for Land
     &     SAOS(L1,2),                    !  and Ice, and for Open Sea,
C     ! respectively, with zeroes for safety where no value applies
     &     SANA(L1,2)
C     ! Grid box mean albedo as it would be with no aerosol.
C
C     !  FTSA has no dynamically allocated workspace, no EXTERNAL calls
C     !  and no significant structure - just a single loop and some
C     !  IF blocks.
C*
      INTEGER J                           ! Loops over points
      REAL DSA,                           ! Deep-snow albedo (alphasubD)
     &     TICE                           ! Surface temperature for
C     !                      the sea-ice covered fraction of a grid-box.
      REAL DTLAND, KLAND, TCLAND, ADIFC, ALPHAC, ALPHAM, DTICE, TCICE,
     &     ICE1, ICE2,                    ! Local PARAMETERs
     &     MASKD                          ! Masking depth (S in 3.6.1)
     &     , ALPHA, BETA                  ! Mass scattering coefficient
C     ! and upward scattering fraction for the anthropogenic sulphate
     &     , AERCON                       ! Their product
      PARAMETER ( MASKD = 0.2 )
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C     ! Basic quantities for land CSSA calculations:
      PARAMETER ( DTLAND = 2.,            ! delta(T) in 3.6.2
     &     FCATM = 0.3 )                  ! Fraction by which deep-snow
C     ! albedo changes (from "cold" value towards snow-free value) at TM
C     ! From these, 2 constants precalculated for efficiency in 3.6.2:
      PARAMETER ( KLAND = 0.3/DTLAND,
     &     TCLAND = TM-DTLAND )
C
      PARAMETER ( ADIFC = 0.06 )          ! Surface albedo of ice-free
C                                         !    sea for the diffuse beam
      PARAMETER ( ALPHA = 8500., BETA = 0.29,
     &             AERCON = ALPHA*BETA )

C     ! derive 3 constants from the basic quantities (supplied in the
C     ! namelist RUNCNST) for sea-ice CSSA calculations:
      TCICE = TM - DTICE
      ICE1 = (ALPHAM-ALPHAC) / DTICE
      ICE2 = ALPHAM - TM*ICE1
C
C
      DO 100 J=1, L2
       IF (LAND(J)) THEN
C 3.6.2:
          IF ( TSTAR(J) .LT. TCLAND ) THEN
             DSA = MDSA(J)
           ELSE
             DSA= MDSA(J) + KLAND * (SFA(J)-MDSA(J)) * (TSTAR(J)-TCLAND)
          ENDIF
C 3.6.1:
          SALI(J,1) = SFA(J) + ( DSA-SFA(J) ) * (1. - EXP(-MASKD*S(J)) )
          SAOS(J,1) = SALI(J,1)
          SAOS(J,2) = SALI(J,1)
        ELSE
          SAOS(J,1) = 0.05 / ( 1.1 * COSZ(J)**1.4 + 0.15 )
          SAOS(J,2) = ADIFC
C         ! Note that the following will add in ICE1*(TSTAR-TFS) to CSSA
C         ! if AICE#0 when it should be - even if only very small: for
C         ! large enough TSTAR this will give very large surface heating
C         ! and even negative atmospheric heating.
          IF ( AICE(J) .EQ. 0. ) THEN
             SALI(J,1)=0.
           ELSE
C            !  Recover TICE from TSTAR:
             TICE = ( TSTAR(J) + (AICE(J)-1.) * TFS )  / AICE(J)
C            !  3.5.1:
             IF ( TICE .LT. TCICE ) THEN
                SALI(J,1) = ALPHAC
              ELSE
                SALI(J,1) = ICE1 * TICE + ICE2
             ENDIF
          ENDIF
       ENDIF

  100 CONTINUE

      IF ( NLALBS .EQ. 2 ) THEN         ! Sulphate aerosol is to be used

!      The standard no-aerosol calculations have now been done.  Next,
!      set the diffuse land/ice albedo to match the direct, calculate
!      the no-aerosol grid-box mean & then put the aerosol forcing in -
!      to the direct beam only - limiting the new albedo to 90% as for
!      low sun it could otherwise exceed 100%.

        IF ( SANAON ) THEN
          DO J=1, L2
            IF (LAND(J)) THEN        ! Compute grid box mean albedo
               SANA(J,1) = SALI(J,1)
               SANA(J,2) = SALI(J,1)
             ELSE
               SANA(J,1) =
     &                AICE(J) * SALI(J,1) + ( 1. - AICE(J) ) * SAOS(J,1)
               SANA(J,2) =
     &                AICE(J) * SALI(J,1) + ( 1. - AICE(J) ) * SAOS(J,2)
            ENDIF
          ENDDO
        ENDIF

        DO J=1, L2

          SALI(J,2) = SALI(J,1)

          IF ( COSZ(J) .GT. 0. ) THEN
            SAOS(J,1) = SAOS(J,1) +
     &         ( (1-SAOS(J,1)) **2 ) * AERCON * SULPH(J) / COSZ(J)
            SALI(J,1) = SALI(J,1) +
     &         ( (1-SALI(J,1)) **2 ) * AERCON * SULPH(J) / COSZ(J)
          ENDIF

          IF ( SAOS(J,1) .GT. 0.9 ) SAOS(J,1) = 0.9
          IF ( SALI(J,1) .GT. 0.9 ) SALI(J,1) = 0.9

        ENDDO
      ENDIF

      RETURN
      END
