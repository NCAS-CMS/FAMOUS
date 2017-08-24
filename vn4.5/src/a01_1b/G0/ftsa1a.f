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
CLL  Subroutine FTSA -------------------------------------------------
CLL
CLL  It calculates (true) surface albedos for P234.
CLL    Release 2.8 of the UM allows for separate surface
CLL  albedos for direct and diffuse light over sea, calculating the
CLL  former from the formula of Briegleb and Ramanathan (1982, J. Appl.
CLL  Met., 21, 1160-1171) and passes the albedos out in different form
CLL  Land & ice albedos are still the same for direct and diffuse beams.
CLL                                            William Ingram 25/9/92
CLL  Suitable for single column model use.
CLL
CLL   Author: William Ingram
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.1   25/1/93     alphaC and alphaM increased if version 2A of
CLL                                       Section 1 is selected.
CLL   3.4   09/09/94     alphaC,alphaM,dTice all set in U.I.  C D Hewitt
CLL   4.0   04/07/95     MASKD set locally.  William Ingram
!
!     4.0   01/09/95     Logical test for version of SW radiation
!                        introduced to allow for differing treatment
!                        of the direct beam by different versions of
!                        the radiation.
!                                           J. M. Edwards
CLL
CLL   4.1   25/01/96     Set SAOS to SALI over land because SAOS is
CLL                      used to find Net Surface SW in Band 1, which
CLL                      is now required over land for photosythesis
CLL                      modelling.
CLL                                         Richard Betts
CLL
CLL   4.4   17/09/97     Optional prognostic snow albedo scheme -
CLL                      separate albedos calculated for visible and
CLL                      near-infrared wavelengths, depending on snow
CLL                      grain size, soot content and zenith angle.
CLL                                         Richard Essery
CLL
CLL   4.4   18/09/97     SNOW_FRAC and TSTAR_SNOW are fractional
CLL                      coverage and surface temperature for snow on
CLL                      land points if MOSES II is selected, set to
CLL                      1 and TSTAR otherwise.
CLL                                         Richard Essery    
CLL   4.5   21.08.98     If l_ssice_albedo, the albedo of sea-ice is
CLL         altered by the presence of snow. Jonathan Gregory
CLL    
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
     &     LAND, AICE, TSTAR, TSTAR_SNOW, SNOW_FRAC, SFA, MDSA, COSZ,
     &     S, RGRAIN, SOOT,
     &     ALPHAM,ALPHAC,ALPHAB,DTICE,L_SSICE_ALBEDO,
     &     VERSION,
     &     L1, L2, L_SNOW_ALBEDO, SAL_DIM, SAL_VIS, SAL_NIR, SALI, SAOS)
!
      implicit none
!
      CHARACTER !, INTENT(IN) ::
     &     VERSION*3                      ! Version of radiation scheme
      INTEGER !, INTENT(IN) ::
     &     L1,                            ! Full field dimension
     &     L2                             ! Number of points to treat
     &    ,SAL_DIM                        ! Dimensions SAL_VIS and
C                                         ! SAL_NIR for prognostic snow
C                                         ! albedo
      LOGICAL !, INTENT(IN) ::
     &     LAND(L1)                       ! Land-sea mask (land .TRUE.)
     &    ,L_SNOW_ALBEDO                  ! Flag for prognostic snow 
C                                         ! albedo
!         Switch on the effect of snow on sea-ice albedo
     &    ,L_SSICE_ALBEDO
      REAL !, INTENT(IN) ::
     &     AICE(L1),                      ! Sea-ice fraction
     &     SNOW_FRAC(L1),                 ! Snow fraction
     &     TSTAR(L1),                     ! Surface temperature
     &     TSTAR_SNOW(L1),                ! Snow surface temperature (K)
     &     SFA(L1),                       ! Snow-free surface albedo
     &     MDSA(L1),                      ! Cold deep-snow albedo
     &     ! (These two are alpha sub 0 & alpha sub S resp. in UMDP 23.)
     &     COSZ(L1),                      ! cos(solar zenith angle)
     &     RGRAIN(L1),                    ! Snow grain size (microns)   
     &     SOOT(L1),                      ! Snow soot content (mass frac
     &     S(L1)                          ! Snow amount (mass/area)
!     Constants used to determine the albedo of sea-ice:
! Albedo of sea-ice at melting point (TM) if .not.l_ssice_albedo, or
! Albedo of snow on sea-ice at melting point (TM) if l_ssice_albedo
     &    ,ALPHAM  ! "M" for "melting"
! Albedo of sea-ice at and below TM-DTICE if .not.l_ssice_albedo, or
! Albedo of snow on sea-ice at and below TM-DTICE if l_ssice_albedo
     &    ,ALPHAC  ! "C" for "cold"
! Albedo of snow-free sea-ice if l_ssice_albedo
     &    ,ALPHAB  ! "B" for "bare"
! Temperature range in which albedo of sea-ice, if .not.l_ssice_albedo,
! or of snow on sea-ice, if l_ssice_albedo, varies between its limits
     &    ,DTICE
      REAL!, INTENT(OUT)
     &     SAL_VIS(SAL_DIM,2),            ! Visible albedo for land 
     &     SAL_NIR(SAL_DIM,2),            ! Near-IR albedo for land
     &     SALI(L1),                      ! Surface Albedos for Land
     &     SAOS(L1,2)                     !  and Ice, and for Open Sea,
C     ! respectively, with zeroes for safety where no value applies
C
C     !  FTSA has no dynamically allocated workspace, no EXTERNAL calls
C     !  and no significant structure - just one loop and some IF blocks
C*
      INTEGER J                           ! Loops over points
      REAL DSA,                           ! Deep-snow albedo (alphasubD)
     &     TICE                           ! Surface temperature for
C     !                      the sea-ice covered fraction of a grid-box.
! Temperature at which (snow on) sea-ice reaches its "cold" value
     &    ,TCICE     
! Slope and intercept of temperature dependence of the albedo of
! (snow on) sea-ice
     &    ,ICE1,ICE2
!     Local parameters
      REAL DTLAND, KLAND, TCLAND, ADIFC, FCATM,
     &     ALB_DIFF,                      ! Diffuse beam snow albedo
     &     ALB_DIR,                       ! Direct beam snow albedo
     &     R0,                            ! Grain size for fresh snow
C                                         ! (microns)
     &     REFF,                          ! Zenith effective grain size
     &     SIGMA,                         ! Scaled soot content
     &     SNOW_ALBEDO,                   ! Snow albedo            
     &     MASK,                          ! Snow masking factor
     &     MASKD                          ! Masking depth (S in 3.6.1)
C     Note that the same masking depth is always used, both for land,
C     regardless of vegetation cover, and for sea-ice. This assumption
C     of constancy may be doubtful.
      PARAMETER ( MASKD = 0.2, R0 = 50. )
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
C     ! derive 3 constants from the basic quantities (supplied in the
C     ! namelist RUNCNST) for sea-ice CSSA calculations:
      TCICE = TM - DTICE
      ICE1 = (ALPHAM-ALPHAC) / DTICE
      ICE2 = ALPHAM - TM*ICE1

C
C

      IF ( L_SNOW_ALBEDO ) THEN
C-----------------------------------------------------------------------
C Prognostic spectral snow albedos
C-----------------------------------------------------------------------
      DO J=1,L2
        SALI(J) = 0.    
       IF (LAND(J)) THEN

C  Increased deep-snow albedo required for land ice points
          IF ( MDSA(J) .GT. 0.78 ) MDSA(J) = 0.83

C  Effective grain size for zenith angle
          REFF = RGRAIN(J) * ( 1. + 0.77*(COSZ(J)-0.65) )**2

C  Visible albedos
          ALB_DIR = 1.2*MDSA(J) - 0.002*(SQRT(REFF)-SQRT(R0))
          ALB_DIFF = 1.2*MDSA(J) - 0.002*(SQRT(RGRAIN(J))-SQRT(R0))
          SIGMA = SOOT(J) * RGRAIN(J) / 0.0017
          IF ( SIGMA .GT. 1. ) THEN
            ALB_DIR = 0.07 + 0.5*(ALB_DIR - 0.07) / (SIGMA**0.46)
            ALB_DIFF = 0.07 + 0.5*(ALB_DIFF - 0.07) / (SIGMA**0.46)
          ELSE
            ALB_DIR = ALB_DIR - 0.5*(ALB_DIR - 0.07)*(SIGMA**0.6) 
            ALB_DIFF = ALB_DIFF - 0.5*(ALB_DIFF - 0.07)*(SIGMA**0.6) 
          ENDIF
          SAL_VIS(J,1) = ALB_DIR
          SAL_VIS(J,2) = ALB_DIFF

C  Near-IR albedos
          ALB_DIR = 0.8*MDSA(J)  -
     &                   0.18*( ALOG(SQRT(REFF)) - ALOG(SQRT(R0)) )
          ALB_DIFF = 0.8*MDSA(J)  -
     &                   0.18*( ALOG(SQRT(RGRAIN(J))) - ALOG(SQRT(R0)) )
          SIGMA = SOOT(J) * RGRAIN(J) / 0.004
          IF ( SIGMA .GT. 1. ) THEN
            ALB_DIR = 0.06 + 0.5*(ALB_DIR - 0.06) / (SIGMA**0.6)
            ALB_DIFF = 0.06 + 0.5*(ALB_DIFF - 0.06) / (SIGMA**0.6)
          ELSE 
            ALB_DIR = ALB_DIR - 0.5*(ALB_DIFF - 0.06)*(SIGMA**0.7)
            ALB_DIFF = ALB_DIFF - 0.5*(ALB_DIFF - 0.06)*(SIGMA**0.7)
          ENDIF
          SAL_NIR(J,1) = ALB_DIR
          SAL_NIR(J,2) = ALB_DIFF

C  Adjust albedos for snow depth and fraction
          MASK = 1.
          IF ( SNOW_FRAC(J).GT.0. ) MASK = EXP(-MASKD*S(J)/SNOW_FRAC(J))
          SAL_VIS(J,1) = (1. - SNOW_FRAC(J))*SFA(J) + SNOW_FRAC(J) *
     &                    (SFA(J) + (SAL_VIS(J,1) - SFA(J))*(1. - MASK))
          SAL_VIS(J,2) = (1. - SNOW_FRAC(J))*SFA(J) + SNOW_FRAC(J) *
     &                    (SFA(J) + (SAL_VIS(J,2) - SFA(J))*(1. - MASK))
          SAL_NIR(J,1) = (1. - SNOW_FRAC(J))*SFA(J) + SNOW_FRAC(J) *
     &                    (SFA(J) + (SAL_NIR(J,1) - SFA(J))*(1. - MASK))
          SAL_NIR(J,2) = (1. - SNOW_FRAC(J))*SFA(J) + SNOW_FRAC(J) *
     &                    (SFA(J) + (SAL_NIR(J,2) - SFA(J))*(1. - MASK))
        ENDIF
      ENDDO  

      ELSE
C-----------------------------------------------------------------------
C Diagnosed all-band snow albedo
C-----------------------------------------------------------------------
      DO J=1,L2
        IF (LAND(J)) THEN
          SALI(J) = SFA(J)
          IF (SNOW_FRAC(J) .GT. 0.0) THEN
            IF ( TSTAR_SNOW(J) .LT. TCLAND ) THEN
             DSA = MDSA(J)
           ELSE
              DSA = MDSA(J) + KLAND * (SFA(J) - MDSA(J)) *
     &                                          (TSTAR_SNOW(J) - TCLAND)
          ENDIF
            SNOW_ALBEDO = SFA(J) + (DSA-SFA(J)) *
     &                            ( 1. - EXP(-MASKD*S(J)/SNOW_FRAC(J)) )
            SALI(J) = (1. - SNOW_FRAC(J))*SFA(J) +
     &                                          SNOW_FRAC(J)*SNOW_ALBEDO
          ENDIF
          SAOS(J,1) = SALI(J)
          SAOS(J,2) = SALI(J)
        ENDIF
      ENDDO

      ENDIF  ! Prognostic or diagnostic snow albedo

      DO 100 J=1, L2
        IF ( .NOT. LAND(J) ) THEN
          IF (VERSION.NE.'03A') THEN
             SAOS(J,1) = 0.05 / ( 1.1 * COSZ(J)**1.4 + 0.15 )
          ELSE
             SAOS(J, 1)=0.026/(COSZ(J)**1.7+0.065)
     &          +0.15*(COSZ(J)-0.1)*(COSZ(J)-0.5)*(COSZ(J)-1.0)
          ENDIF
          SAOS(J,2) = ADIFC
C         ! Note that the following will add in ICE1*(TSTAR-TFS) to CSSA
C         ! if AICE#0 when it should be - even if only very small: for
C         ! large enough TSTAR this will give very large surface heating
C         ! and even negative atmospheric heating.  Check whether this
C         ! could occur.
          IF ( AICE(J) .EQ. 0. ) THEN
             SALI(J) = 0.
           ELSE
C            !  Recover TICE from TSTAR:
             TICE = ( TSTAR(J) + (AICE(J)-1.) * TFS )  / AICE(J)
             if (l_ssice_albedo) then
               if (s(j).gt.0.0) then
                 if (tice.gt.tcice) then
                   snow_albedo=ice2+ice1*tice
                 else
                   snow_albedo=alphac
                 endif
                 sali(j)=alphab
     &           +(snow_albedo-alphab)*(1.0-exp(-maskd*s(j)))
               else
                 sali(j)=alphab
               endif
             else
C            !  3.5.1:
             IF ( TICE .LT. TCICE ) THEN
                SALI(J) = ALPHAC
              ELSE
                SALI(J) = ICE1 * TICE + ICE2
             ENDIF
             endif
          ENDIF
       ENDIF
  100 CONTINUE
C
      RETURN
      END
