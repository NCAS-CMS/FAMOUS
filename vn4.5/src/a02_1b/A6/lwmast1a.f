CLL   (Comment here referred to the former deckname - deleted MJH)
CLL     Purpose:
CLL  It calculates net longwave fluxes (and optionally flux diagnostics)
CLL  from the Planck flux differences found by LWPLAN, transmissivities
CLL  found by LWTRAN, and cloud arrays filled by LWCLD.
CLL  If UPDATE *DEF CRAY is off, the code is standard FORTRAN 77 except
CLL  for having ! comments (it then sets the "vector length" to 1) but
CLL  otherwise it includes automatic arrays also.
CLL
CLL     Author: William Ingram
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   4.2    Sept.96  T3E migration: *DEF CRAY removed;
CLL                   *DEF T3E used for T3E library functions;
CLL                   dynamic allocation no longer *DEF controlled;
CLL                   cray HF functions replaced by T3E lib functions.
CLL                       S.J.Swarbrick
!     4.4    10/4/97  Pass logical through to LWCLD to indicate the
!                     prognostic cloud ice should be used. AC Bushell
CLL
CLL  It conforms with standard A of version 3 (07/9/90) of UMDP 4, and
CLL  contains no 8X-deprecated features.
CLL
CLL    Logical components covered: P232 , D23...
CLL  It is the top-level plug-compatible routine in component P232
CLL  (longwave radiation)
CLL
CLL   It also performs some of
CLL  the functions of D23 (radiation diagnostics).
CLL
CLL  System task P23 (radiation).
CLL
CLL  Offline documentation is in UMDP 23.
CLLEND ---------------------------------------------------------
C*L
      SUBROUTINE LWMAST (H2O, CO2, O3, TAC, PEXNER, TSTAR, PSTAR, AB,
     &     BB, AC, BC, AICE, LCA, LCCWC1, LCCWC2, CCA, CCCWP, CCB, CCT,
     &     LUT,
     &     CSOLRD, CSOLON, SFDN, SFDNON, CSSFDN, CSSDON,
     &     L_CLOUD_WATER_PARTITION,
     &     L2, NLEVS, NCLDS,                                            
     &     NWET, NOZONE, L1,                           SEAFX, FLUX)
C*
      EXTERNAL LWCLD, LWPLAN, LWPTSC, LWTRAN
      INTEGER NBANDS          ! Number of spectral bands in the longwave
      PARAMETER (NBANDS=6)    ! This run uses the standard set of
!  "ECMWF" bands as described by Morcrette et al (J.-J. Morcrette,
!   L.D. Smith & Y. Fouquart, 1986, Beitr. Phys. Atmosph., 59, 455-469)
      INTEGER NGASES
C Effective number of absorbing gases treated in the longwave
      PARAMETER (NGASES=12)
C     ! This set is for the "ECMWF-like" code.  Gases mostly have
C     ! different pressure and temperature scaling of their pathlengths
C     ! in each band, so that there are 6 absorber amounts for water
C     ! vapour line absorption, 1 for each of the foreign-broadened &
C     ! self-broadened water vapour continua, 2 for CO2 (for 3 different
C     ! bands) and 2 for ozone (with and without pressure broadening,
C     ! though only in one band).
      INTEGER NTRANS, NDIATR
      PARAMETER (NDIATR=0)
      PARAMETER (NTRANS=NBANDS+NDIATR)
C Number of transmissivities to be calculated - one for each band is
C needed to construct the actual fluxes, but we also allow for NDIATR
C for diagnostic uses, such as possible "narrow-band" flux diagnostics.
      INTEGER IT              ! Dimension of look-up tables
!  Dimension of Lookup Tables
      PARAMETER (IT=25)
C     !  Array dimensions must be constants in FORTRAN:
C*L
      INTEGER!, INTENT(IN) ::
     &     L2,                       ! Number of points to be treated
     &     NLEVS,                    ! Number of levels
     &     NCLDS,                    ! Number of possibly cloudy levels
     &     NWET,                     ! Number of levels with moisture
     &     NOZONE,                   ! Number of levels with ozone
     &     L1                        ! Full field dimension
      REAL!, INTENT(IN) ::
     &     H2O(L1,NWET), CO2,        ! Mixing ratios of the three
     &     O3(L1,NOZONE),            !               absorbing gases
     &     TAC(L1,NLEVS),            ! Temperature at layer centres
     &     PEXNER(L1,NLEVS+1),       ! Exner function @ layer boundaries
     &     TSTAR(L1), PSTAR(L1),     ! Surface temperature & pressure
     &     AC(NLEVS), BC(NLEVS),     ! A & B for layer centres and
     &     AB(NLEVS+1), BB(NLEVS+1), !                       boundaries
     &     LUT(IT,NTRANS,NGASES,2),  ! Look-up tables for LWTRAN
     &     AICE(L1),                 ! Sea-ice fraction
     &     LCCWC1(L1,1/(NCLDS+1)+NCLDS), LCCWC2(L1,1/(NCLDS+1)+NCLDS),
C     ! Layer cloud condensed water contents (specific contents, mass
C     ! per unit mass).  Only the sum of these two fields is used.
     &     LCA(L1,1/(NCLDS+1)+NCLDS),! Layer cloud fractional cover
     &     CCCWP(L1),                ! Convective cloud fractional cover
     &     CCA(L1)                   !          and condensed water path
C     ! (LWCLD describes precisely which these cloud quantities are.)
      INTEGER!, INTENT(IN) ::
     &     CCB(L1), CCT(L1)          ! Convective cloud base and top
      LOGICAL!, INTENT(IN)
     &     CSOLON                    ! Is CSOLRD wanted ?
     &     , SFDNON                  ! And is SFDN ?
     &     , CSSDON                  ! & is CSSFDN ?
C                                    ! (If so, SFDNON must be set too)
     &     ,L_CLOUD_WATER_PARTITION  ! Is cloud ice prognostic?
      REAL!, INTENT(OUT) ::
     &     FLUX(L1,NLEVS+1),         ! Net longwave flux (+ downwards)
     &     CSOLRD(L1),               ! diagnosed clear-sky OLR
     &     SFDN(L1),                 ! Diagnosed downward surface flux
     &     CSSFDN(L1),               ! and its clear-sky equivalent
     &     SEAFX(L1)                 ! Term calculated by LWPLAN & used
      ! by LWRAD to derive open-sea-only flux at sea-ice points.
C*
CL    !  After zeroing FLUX and SEAFX,
CL    !  LWMAST calls LWPLAN, LWPTSC & LWCLD to set up arrays.
CL    !  and initializes an array DNSRCE ("DO 10" loop).
CL    !  Then it adds in the half-layer terms for each layer to the
CL    !  fluxes at the boundaries of that layer ("DO 2" loop).
CL    !  Most of the code is inside the "DO 3" loop and calculates
CL    !  the contribution to the flux at every layer boundary from every
CL    !  other.  That loop supplies the lower level being treated.
CL    !  Loops inside it ("DO 30" for layers that may contain cloud
CL    !  and "DO 38" for the others) supply the upper layer.
      REAL                           ! WORKSPACE
     &     DB(L2,NLEVS,NBANDS,2),    ! Differences of the black-body
C     ! flux across bottom and top halves of layers,
C     ! DB(,LEVEL,,1) between half-level LEVEL-1/2 & full-level LEVEL,
C     ! DB(,LEVEL,,2) between full-level LEVEL & half-level LEVEL+1/2.
     &     TRANS(L2,NTRANS),         ! Transmissivities
     &     DPATH(L2,NGASES,NLEVS),   ! Scaled gas pathlengths for each
     &     PATH(L2,NGASES),          !        layer & current total path
     &     ECA(L2,NCLDS+1/(NCLDS+1),NBANDS),
C     ! Effective clear-sky fraction: 1-ECA is cloud amount*emissivity
     &     ECTA(L2,NCLDS,NBANDS),
     &     ECBA(L2,NCLDS,NBANDS),     ! Effective amount of
C     ! cloud (amount*emissivity) having its top or bottom in each layer
     &     EFFTRA, UPSRCE(L2,NBANDS),
C     ! EFFTRA is a temporary product of TRANS and a clear-sky term
     &     DNSRCE(L2,NBANDS,2:(NCLDS+2/(NCLDS+1))),
C     !       NCLDS+2/(NCLDS+1)=MAX(NCLDS,2) if NCLDS>=0
     &     CLRF(L2,NBANDS), CLDCLB
C     ! Because random overlap of clouds in different layers is assumed,
C     ! a dB source term for upward and downward flux at each layer
C     ! boundary, from so much cloud surface (top or base respectively)
C     ! and so much clear sky, can be pre-calculated.  These are UPSRCE
C     ! and DNSRCE respectively.  Then the contribution from one layer
C     ! to the flux at any other is given by the dB source term, the
C     ! gaseous transmissivity between the two layers, and the fraction
C     ! of the grid-box where the two layers' view of each other is not
C     ! blocked by intervening cloud, CLRF.  CLDCLB is defined so that
C     ! (1-CLDCLB) is the effective cloud cover crossing a layer
C     ! boundary.
      INTEGER FLSTBD(NGASES,2)       !  First and last band in which
C                                    !  each effective gas is active
      INTEGER BAND, GAS,             !  Loopers over band, absorbing
     &     LEVEL, LEVEL2, J,         !          gas, levels and points
     &     FSCLYR                    !  Start of the "DO 38" loop
CL
CL    ! SECTION 1
CL
CL    ! 1.1 Zero output space
CL
      DO 1 LEVEL=1, NLEVS+1
Cfpp$  Select(CONCUR)
       DO 1 J=1, L2
        FLUX(J,LEVEL) = 0.
    1 CONTINUE
CL
CL    ! 1.11 zero SEAFX
CL
      DO 11 J=1, L2
       SEAFX(J) = 0.
   11 CONTINUE
CL
CL    !  1.12  Zero CSOLRD:
CL
      IF ( CSOLON ) THEN
        DO J=1, L2
          CSOLRD(J) = 0.
        ENDDO
      ENDIF
CL
CL    !  and SFDN:
CL
      IF ( SFDNON ) THEN
        DO J=1, L2
          SFDN(J) = 0.
        ENDDO
      ENDIF
C
C     ! and CSSFDN:
C
      IF ( CSSDON ) THEN
        DO J=1, L2
          CSSFDN(J) = 0.
        ENDDO
      ENDIF
CL
CL    ! 1.2  Set up dB arrays from temperature arrays
CL
Cfpp$ Expand
      CALL LWPLAN (TAC, PEXNER, PSTAR, AB, BB, TSTAR, AICE,
     &     SFDN, SFDNON,
     &     L2, NLEVS, L1,                    SEAFX,  DB, DB(1,1,1,2))
C
CL
CL    ! 1.3 Set up arrays of scaled pathlengths (Eqs 2.3.1 to 2.3.10)
CL
Cfpp$ Expand
      CALL LWPTSC (H2O, CO2, O3, PSTAR, AC, BC, AB, BB, TAC,
     &     L2,                                                          
     &     NWET, NOZONE, NLEVS, L1,                          DPATH)
CL
CL    1.4 Set arrays of effective amount of clear sky, cloud base & top
CL
      IF ( NCLDS .GT. 0 ) THEN
Cfpp$ Expand
        CALL LWCLD (LCA, LCCWC1, LCCWC2, CCA, CCCWP, CCB, CCT, TAC,
     &     PSTAR, AB, BB, L_CLOUD_WATER_PARTITION, L1, NLEVS, NCLDS,
     &     L2,                                                          
     &     ECTA, ECBA,
     &     ECA)
      ENDIF
CL
CL    !  and set DNSRCE
CL
      DO 10 LEVEL=2, NCLDS
       DO 10 BAND=1, NBANDS
Cfpp$   Select(CONCUR)
        DO 10 J=1, L2
C        ! This division by CLDCLB, as well as the 2 later, would fail
C        ! if the inputs specified total cover by a black cloud more
C        ! than one layer thick - but not for the physically identical
C        ! case of total cover by black clouds in adjacent layers.
         CLDCLB = ECA(J,LEVEL,BAND) + ECBA(J,LEVEL,BAND)
         DNSRCE(J,BAND,LEVEL) =
     &      ECA(J,LEVEL,BAND) * DB(J,LEVEL,BAND,1) / CLDCLB +
     &         DB(J,LEVEL-1,BAND,2)
   10 CONTINUE
CL
CL    ! SECTION 2
CL
CL    !  Add in the "half-layer" contributions.
CL    !  Transmissivities are calculated from pathlengths which are
CL    !  a quarter those for the full layers (Eq 2.1.9):
CL
      DO 2 LEVEL=1, NLEVS
       DO 21 GAS=1, NGASES
Cfpp$   Select(CONCUR)
        DO 21 J=1, L2
         PATH(J,GAS) = .25 * DPATH(J,GAS,LEVEL)
   21  CONTINUE
C
       CALL LWTRAN (PATH, LUT, LUT(1,1,1,2), FLSTBD, KCONT,
     &     L2,                                                          
     &     TRANS)
       IF (LEVEL.LE.NCLDS) THEN
C         !  In levels low enough that cloud may occur, there is no
C         !  radiative flux in the part of the grid-box covered by the
C         !  equivalent black-body cloud.
          DO 22 BAND=1, NBANDS
Cfpp$      Select(CONCUR)
           DO 22 J=1, L2
            EFFTRA = TRANS(J,BAND) * ECA(J,LEVEL,BAND)
            FLUX(J,LEVEL) = FLUX(J,LEVEL) + EFFTRA * DB(J,LEVEL,BAND,1)
            FLUX(J,LEVEL+1) = FLUX(J,LEVEL+1) +
     &                          EFFTRA * DB(J,LEVEL,BAND,2)
   22     CONTINUE
        ELSE IF (LEVEL.LT.NLEVS) THEN
C         !  Further up, Eq 2.1.9 applies simply:
          DO 23 BAND=1, NBANDS
Cfpp$      Select(CONCUR)
           DO 23 J=1, L2
            FLUX(J,LEVEL) = FLUX(J,LEVEL) +
     &               TRANS(J,BAND) * DB(J,LEVEL,BAND,1)
            FLUX(J,LEVEL+1) = FLUX(J,LEVEL+1) +
     &               TRANS(J,BAND) * DB(J,LEVEL,BAND,2)
   23     CONTINUE
        ELSE !IF (LEVEL.EQ.NLEVS)
C         !  except right at the top, where the toa flux gets special
C         !  treatment (no transmissivity):
          DO 24 BAND=1, NBANDS
Cfpp$      Select(CONCUR)
           DO 24 J=1, L2
            FLUX(J,LEVEL) = FLUX(J,LEVEL) +
     &               TRANS(J,BAND) * DB(J,LEVEL,BAND,1)
            FLUX(J,NLEVS+1) = FLUX(J,NLEVS+1) + DB(J,NLEVS,BAND,2)
            IF ( CSOLON ) CSOLRD(J) = CSOLRD(J) - DB(J,NLEVS,BAND,2)
   24     CONTINUE
       ENDIF
       IF ( CSSDON .AND. LEVEL .EQ. 1 ) THEN
         DO 221 BAND=1, NBANDS
Cfpp$      Select(CONCUR)
           DO J=1, L2
             CSSFDN(J) = CSSFDN(J) + TRANS(J,BAND) * DB(J,LEVEL,BAND,1)
           ENDDO
  221    CONTINUE
       ENDIF
    2 CONTINUE
CL
CL    ! Separate DB for each half-layer are needed for the "half-layer"
CL    ! terms and the cloud boundary (and surface) source terms.  Now
CL    ! the "half-layer" terms have been dealt with, they can be
CL    ! combined above all clouds, to save calculations later in the
CL    ! "DO 36" loop.
CL
      DO 20 BAND=1, NBANDS
       DO 20 LEVEL=NCLDS+1, NLEVS-1
Cfpp$   Select(CONCUR)
        DO 20 J=1, L2
         DB(J,LEVEL,BAND,2) = DB(J,LEVEL,BAND,2) + DB(J,LEVEL+1,BAND,1)
   20 CONTINUE
CL
CL    ! SECTION 3
CL
CL    ! Now have the full-level terms to find, with a contribution from
CL    ! every layer boundary to the flux at every other layer boundary.
CL
C     ! The contributions from each of a pair of layers to the other are
C     ! added in at once so that transmissivities or pathlengths do not
C     ! have to be stored for all combinations nor calculated twice.
C     !   LEVEL and LEVEL2 are the lower and upper layer boundaries
C     ! respectively.  Layer boundaries are conventionally indexed by
C     ! half-integers in the documentation, and so for FORTRAN we must
C     ! add or subtract a half.
C     ! We currently add it, so the numerical value of LEVEL is the
C     ! number of the layer centre ABOVE it, consistent with the
C     ! indexing of most arrays (not ECTA, which is indexed by LEVEL-1)
C
C     ! The source term for the downward flux (i.e. from the upper layer
C     ! boundary, LEVEL2), DNSRCE, is needed on different passes through
C     ! the outer "DO 3" loop, so that the whole lot is precalculated
C     ! and stored, but the upward one, UPSRCE, is only needed for the
C     ! current value of LEVEL.  The calculations are complicated by the
C     ! need to take account of clouds more than one level thick.  This
C     ! means that over part of the grid-box neither the half-level dB
C     ! (for cloud top or bottom) nor the full-level one (for cloud-free
C     ! space) may need to be added in.  This is not easily allowed for
C     ! by CLRF alone.  The method adopted is to normalize UPSRCE and
C     ! DNSRCE to be the mean dB terms from the parts of the layer
C     ! boundary which do contribute, and initialize CLRF to CLDCLB
C     ! to take account of the effect of any cloud that does cross the
C     ! layer boundary at LEVEL.
      DO 3 LEVEL=1, NLEVS
C      ! Start by setting the pathlength to that for the layer above the
C      ! layer boundary where FLUX(,LEVEL) applies, and also initialize
C      ! overlap quantities.
       DO 31 GAS=1, NGASES
Cfpp$   Select(CONCUR)
        DO 31 J=1, L2
         PATH(J,GAS) = DPATH(J,GAS,LEVEL)
   31  CONTINUE
       DO 32 BAND=1, NBANDS
Cfpp$   Select(CONCUR)
        DO 32 J=1, L2
         IF (LEVEL.EQ.1) THEN
            UPSRCE(J,BAND) = DB(J,1,BAND,1)
            CLRF(J,BAND) = 1.
          ELSE IF (LEVEL.LE.(NCLDS+1)) THEN
            CLDCLB = ECA(J,LEVEL-1,BAND) + ECTA(J,LEVEL-1,BAND)
            UPSRCE(J,BAND) = DB(J,LEVEL,BAND,1)
     &             + ECA(J,LEVEL-1,BAND) * DB(J,LEVEL-1,BAND,2) / CLDCLB
            CLRF(J,BAND) = CLDCLB
          ELSE     ! IF LEVEL > NCLDS + 1
            UPSRCE(J,BAND) = DB(J,LEVEL-1,BAND,2)
            CLRF(J,BAND) = 1.
         ENDIF
   32  CONTINUE
       DO 30 LEVEL2=LEVEL+1, NCLDS
        CALL LWTRAN (PATH, LUT, LUT(1,1,1,2), FLSTBD, KCONT,
     &     L2,                                                          
     &     TRANS)
        IF ( CSSDON .AND. LEVEL .EQ. 1 ) THEN
          DO 363 BAND=1, NBANDS
Cfpp$       Select(CONCUR)
            DO J=1, L2
              CSSFDN(J) = CSSFDN(J) +
     &  TRANS(J,BAND) * ( DB(J,LEVEL2,BAND,1) + DB(J,LEVEL2-1,BAND,2) )
            ENDDO
  363     CONTINUE
        ENDIF
        DO 33 BAND=1, NBANDS
Cfpp$    Select(CONCUR)
         DO 33 J=1, L2
C         ! Since this code works up from LEVEL, a cloud starts having
C         ! an effect on the overlaps (other than through CLDCLB) when
C         ! its base is reached, and its effect is constant thereafter.
          CLDCLB = ECA(J,LEVEL2-1,BAND) + ECBA(J,LEVEL2-1,BAND)
          CLRF(J,BAND) = CLRF(J,BAND) *
     &                      ( 1. - ECBA(J,LEVEL2-1,BAND) / CLDCLB )
          EFFTRA = TRANS(J,BAND) * CLRF(J,BAND)
C
          FLUX(J,LEVEL) = FLUX(J,LEVEL) + EFFTRA * DNSRCE(J,BAND,LEVEL2)
          FLUX(J,LEVEL2) = FLUX(J,LEVEL2) + EFFTRA * UPSRCE(J,BAND)
C
   33   CONTINUE
C       ! Add in the next layer's contributions to the gas pathlengths.
        DO 34 GAS=1, NGASES
Cfpp$    Select(CONCUR)
         DO 34 J=1, L2
          PATH(J,GAS) = PATH(J,GAS) + DPATH(J,GAS,LEVEL2)
   34   CONTINUE
   30  CONTINUE
C      ! For layers above all cloud, use the last values of the
C      ! cloud overlap terms (or the initialized ones if LEVEL>=NCLDS)
C      !  - otherwise the physics is the same.
       FSCLYR=LEVEL2                    ! Next layer boundary to do is
C      ! MAX(NCLDS+1,LEVEL+1), where there cannot be a cloud term in the
C      ! downward source, but there may still be for CLRF
       IF (LEVEL.LE.NCLDS) THEN
         DO 35 BAND=1, NBANDS
Cfpp$     Select(CONCUR)
          DO 35 J=1, L2
           CLRF(J,BAND) = CLRF(J,BAND) * (1.-ECBA(J,NCLDS,BAND))
   35    CONTINUE
       ENDIF
       DO 38 LEVEL2=FSCLYR, NLEVS+1
        CALL LWTRAN (PATH, LUT, LUT(1,1,1,2), FLSTBD, KCONT,
     &     L2,                                                          
     &     TRANS)
        IF ( CSSDON .AND. LEVEL .EQ. 1 ) THEN
          DO 366 BAND=1, NBANDS
Cfpp$       Select(CONCUR)
            DO J=1, L2
              CSSFDN(J) =
     &            CSSFDN(J) + TRANS(J,BAND) * DB(J,LEVEL2-1,BAND,2)
              IF ( LEVEL2 .EQ. NCLDS+1 ) CSSFDN(J) =
     &            CSSFDN(J) + TRANS(J,BAND) * DB(J,LEVEL2,BAND,1)
            ENDDO
  366     CONTINUE
        ENDIF
        DO 36 BAND=1, NBANDS
Cfpp$    Select(CONCUR)
         DO 36 J=1, L2
C
          EFFTRA = TRANS(J,BAND) * CLRF(J,BAND)
          FLUX(J,LEVEL) = FLUX(J,LEVEL) + EFFTRA * DB(J,LEVEL2-1,BAND,2)
          IF ( LEVEL2 .EQ. NCLDS+1 )
     &      FLUX(J,LEVEL) = FLUX(J,LEVEL) + EFFTRA * DB(J,LEVEL2,BAND,1)
          FLUX(J,LEVEL2) = FLUX(J,LEVEL2) + EFFTRA * UPSRCE(J,BAND)
C
   36   CONTINUE
        IF (LEVEL2.LT.NLEVS+1) THEN
C         ! Add in the next contribution to the gas pathlengths.
          DO 37 GAS=1, NGASES
Cfpp$      Select(CONCUR)
           DO 37 J=1, L2
            PATH(J,GAS) = PATH(J,GAS) + DPATH(J,GAS,LEVEL2)
   37     CONTINUE
        ENDIF
   38  CONTINUE
CL     !  Put in the contributions to CSOLRD:
       IF ( CSOLON ) THEN
         DO 39 BAND=1, NBANDS
           DO J=1, L2
             IF ( LEVEL .LE. NCLDS+1 ) CSOLRD(J) = CSOLRD(J) -
     &                              TRANS(J,BAND) * DB(J,LEVEL,BAND,1)
             IF ( LEVEL .GT. 1 ) CSOLRD(J) = CSOLRD(J) -
     &                              TRANS(J,BAND) * DB(J,LEVEL-1,BAND,2)
           ENDDO
   39    CONTINUE
       ENDIF
    3 CONTINUE
C
C
CL    !  Change CSSFDN from the net downward flux which has been found
CL    !    so far to the downward flux wanted:
      IF ( CSSDON ) THEN
        DO J=1, L2
          CSSFDN(J) = SFDN(J) + CSSFDN(J)
        ENDDO
      ENDIF
C
CL    !  Change SFDN from the upward flux returned by LWPLAN to the
CL    !   downward flux wanted:
      IF ( SFDNON ) THEN
        DO J=1, L2
          SFDN(J) = SFDN(J) + FLUX(J,1)
        ENDDO
      ENDIF
      RETURN
      END
