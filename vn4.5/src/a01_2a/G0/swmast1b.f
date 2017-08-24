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
CLL Subroutine SWMAST   ----------------------------------------------
CLL
CLL Purpose :
CLL        The version with COMPATHS on.
CLL  It is the top-level plug-compatible routine in component P234
CLL  (interaction of shortwave radiation with the atmosphere), which is
CLL  in task P23 (radiation).
CLL  It acts as the master routine for P234, assembling the net solar
CLL  flux (normalized by the incoming insolation at the top of the
CLL  atmosphere) by considering the various beams and calling various
CLL  specialized routines.
CLL  Before SWMAST is called, SWLKIN (in deck SWTRAN) must be CALLed to
CLL  initialize LUT for SWTRAN.
CLL        Release 2.8 of the UM) uses different surface albedos
CLL    for direct and diffuse light, which in turn means that two
CLL    quantities that SWRAD used to calculate from FLUX and the surface
CLL    albedos now have to be found here - TDSS and DSFLUX.
CLL                                             William Ingram 25/9/92
CLL
CLL   Author: William Ingram
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   4.2    Sept.96  T3E migration: *DEF CRAY removed;
CLL                   *DEF T3E used for T3E library functions;
CLL                   dynamic allocation no longer *DEF controlled.
CLL                       S.J.Swarbrick
CLL
CLL Programming standard :
CLL  It conforms to standard A of UMDP 4 (version 3, 07/9/90), and
CLL  includes no features deprecated in 8X.
CLL  If *DEF CRAY is not set, the code is standard FORTRAN 77 except for
CLL  having ! comments (it then sets the "vector length" to be 1) but
CLL  otherwise it includes CRAY automatic arrays also.
CLL
CLL Logical components covered : P234
CLL
CLL Project task : P23
CLL
CLL External documentation: UMDP 23.
CLL
CLLEND -----------------------------------------------------------------
C*L
      SUBROUTINE SWMAST (H2O, CO2, O3, PSTAR, AB, BB, LCA, LCCWP,
     &     LRE, CCA, CCCWP, CRE, CCB, CCT, COSZ, TSA, DTSA, TRTAB,
     &     CSOSDI, CSOSON, NSSB1, NSS1ON, TDSS, TDSSON,
     &     CSSSD, CSSSDO, CSSSU, CSSSUO, LCAAR, LCAARO, LCAARL, LCAARB,
     &     LCAAF, LCAAFO, LCAAFL, LCAAFB, CCAAR, CCAARO, CCAARB, CCAAF,
     &     CCAAFO, CCAAFB,
     &     L2, NLEVS, NCLDS,                                            
     &     NWET, NOZONE, L1, L3,                 DSFLUX, FLUX)
C*
C     !  SWMAST has 4 EXTERNAL calls
      EXTERNAL SWTRAN, SWCLOP, SWMSAL, SWPTSC
      INTEGER NBANDS         ! Number of spectral bands in the shortwave
      PARAMETER (NBANDS=4)   ! This run uses the standard set of
C shortwave bands as described by Slingo (May 1989, J.Atmos.Sci., 46,
C 10, 1419-1427) or UMDP 23.
      INTEGER NGASES
C Number of absorbing gases treated in the shortwave
      PARAMETER (NGASES=3)    !  Standard set is water vapour, ozone
C                             !  and carbon dioxide.
      INTEGER NTRANS, NDIATR
      PARAMETER (NDIATR=0)
      PARAMETER (NTRANS=NBANDS+NDIATR)
C Number of transmissivities to be calculated - one for each band is
C needed to construct the actual fluxes, but we also allow for NDIATR
C for diagnostic uses, such as possible "narrow-band" flux diagnostics.
      INTEGER NLKUPS          ! DIMENSION OF LOOK-UP TABLES
      PARAMETER (NLKUPS=50)
      INTEGER!, INTENT (IN)
     &     L2,                         ! Number of points to be treated
     &     NLEVS,                      ! Number of levels
     &     NCLDS,                      ! Number of possibly cloudy ones
     &     NWET,                       ! Number of levels with moisture
     &     NOZONE,                     ! Number of levels with ozone
C     ! Need 0 =< NCLDS < NLEVS, 0 =< NWET =< NLEVS, 0 < NOZONE =< NLEVS
     &     L1,                         ! First dimension of input arrays
     &     L3,                         ! First dimension of flux output
     &     CCB(L1), CCT(L1)
C     ! Convective cloud base & top, counting down from the top, and in
C     ! terms of lowest and highest full layers occupied.  Thus
C     !  CCT(SW)=NLEVS+2-CCT(LW),  CCB(SW)=NLEVS+1-CCB(LW),
C     ! and a one-layer-thick con cloud has CCB=CCT.
      REAL!, INTENT (IN)
     &     H2O(L1,NWET), CO2,          ! Mass mixing ratio (mK in UMDP
     &     O3(L1,NOZONE),              !       23) of each absorbing gas
     &     COSZ(L1),                   ! Cos(solar zenith angle)
     &     PSTAR(L1),                  ! Surface pressure
     &     AB(NLEVS+1), BB(NLEVS+1),   ! As and Bs at layer boundaries
     &     LCA(L1,NLEVS-NCLDS+1:NLEVS),! Layer cloud amount, condensed
     &   LCCWP(L1,NLEVS-NCLDS+1:NLEVS),!     water path and effective
     &     LRE(L1,NLEVS-NCLDS+1:NLEVS),!     cloud droplet radius.
     &     CCA(L1),                    ! The same for convective cloud.
     &     CCCWP(L1),                  !
     &     CRE(L1),                    !
     &     TSA(L1,NBANDS,2),           ! True surface albedo - mean over
C     ! the whole grid-box for each band for direct & then diffuse light
     &     DTSA(L1,NBANDS,2),          ! True surface albedo -
C     !  different value for some specific part of the grid-box where
C     !  separate calculations are wanted.
     &     TRTAB(NLKUPS,NTRANS,NGASES,2)
C     !    Look-up tables of transmissivities for each gas and of
C     !    differences of their successive elements.
      LOGICAL!, INTENT(IN)
     &     CSOSON, NSS1ON              !  Are CSOSDI and NSSB1 wanted ?
     &     , CSSSDO, CSSSUO            !      & are CSSSD and CSSSU,
     &     , LCAARO, LCAAFO            !            LCAAR and LCAAF,
     &     , CCAARO, CCAAFO            !            CCAAR and CCAAF ?
     &     , LCAARL(NCLDS),  LCAARB(NBANDS), LCAAFL(NCLDS)
     &     , LCAAFB(NBANDS), CCAARB(NBANDS), CCAAFB(NBANDS)
C     !  If L/C AA R/L are wanted, on which (levels and) bands ?
C     !  (The levels are listed from the surface up in these.)
     &     , TDSSON                    !       & is TDSS ?
      REAL!, INTENT (OUT)
     &     FLUX(L3,0:NLEVS)            ! Net downward solar flux, as a
C                                      ! fraction of the incoming solar
     &     , DSFLUX(L3)                ! Net downward flux at the
C                                      !   surface where DTSA applies
     &     , CSOSDI(L1)                ! Diagnosed clear-sky outgoing SW
     &     , NSSB1(L1)                 !  and net surface flux in band 1
     &     , CSSSD(L1)                 ! Clear-sky total downward &
     &     , CSSSU(L1)                 !  upward SW flux at the surface
     &     , LCAAR(L3,*)               ! Layer/Convective Cloud Amount
     &     , LCAAF(L3,*)               !    * Albedo to diRect and
     &     , CCAAR(L3,*)               !    diFfuse light (set to zero
     &     , CCAAF(L3,*)               !    at night points)
C                                      ! for the area DTSA applies to
     &     , TDSS(L1)                  ! Total downward solar flux at
C     !        the surface (counting multiply reflected light multiply).
C
C*
CL    !  SWMAST has lots of dynamically allocated workspace:
CL    ! L2*
CL    ! ( NGASES*(NLEVS+3) + NBANDS*(4*NCLDS+10) + 3 )
      REAL PATH(L2,NGASES,2),          ! Scaled gas pathlengths for the
C     ! total paths to the current layer for direct and diffuse beams.
     &     GREY(L2,NBANDS,3),
C     !  Grey factor for each beam and band (fraction of the incoming
C     !  insolation in that band which would be in that beam at the
C     !  current level, allowing for clouds but not gaseous absorption).
C     !  The last dimension indexes the direct beam (1), total diffuse
C     !  light (2) and the light currently in convective cloud (3).
     &     RFGREY(L2,NBANDS),
     &     RFPATH(L2,NGASES),
C     !  Similarly, grey factors and pathlengths for whichever
C     !   reflected beam is currently being treated.
     &     DPATH(L2,NGASES,NLEVS),     !  Scaled absorber pathlengths
C     ! for each layer, crossed vertically.  Added up after multiplying
C     ! by terms allowing for the angular magnification, these give PATH
     &     GTRANS(L2,NBANDS),          !  Gaseous transmissivities
     &     CTRANS(L2,NBANDS,NLEVS-NCLDS+1-1/(NCLDS+1):NLEVS,2),
     &     REF(L2,NBANDS,NLEVS-NCLDS+1-1/(NCLDS+1):NLEVS,2),
C     ! Cloud transmissivity and reflectivity for direct and diffuse
C                                              radiation respectively.
     &     CCTRANS(L2,NBANDS,2),       ! The same for convective cloud
     &     CCREF(L2,NBANDS,2),         !                     only
     &     DIRFAC(L2),                 ! Magnification factor for the
C                                      !                     direct beam
     &     MODSA(L2,NBANDS,2),         ! Surface albedo TSA modified to
C                                      ! allow for multiple reflections
     &     TOWTGR(L2),                 ! Total grey factor for the
C     ! diffuse beam for all the "wet" bands (FSWTBD to NBANDS)
     &     NEWDIF(L2)                  ! Contribution to TOWTGR from
C                                      !      light newly made diffuse
      INTEGER J, BEAM,                 ! Loopers over points, beams,
     &     BAND, GAS, DIRDIF,          !   bands, gases, direct versus
     &     LEVEL, LEVEL2,              !   diffuse beam, levels, and
     &     INIT                        !   values being initialized
      REAL DIFFAC                      ! Diffusivity factor
      PARAMETER ( DIFFAC = 1.66 )
      REAL HBYA,                       ! Typical height/earth's radius,
     &     HBYAP1,                     ! and constants based on it used
     &     HBYAX2                      ! to find DIRFAC
      PARAMETER ( HBYA = .004, HBYAP1 = 1. + HBYA, HBYAX2 = 2. * HBYA )
      REAL RAYSCA,                     ! Fraction of incoming sunlight
C     ! assumed to be Rayleigh-scattered back to space immediately.
     &     ONELRS                      ! 1-RAYSCA
      PARAMETER ( RAYSCA = .03, ONELRS = 1. - RAYSCA )
      INTEGER FSWTBD          ! First "wet" band (water vapour absorbs)
      PARAMETER (FSWTBD=2)
      REAL FSCIEB(NBANDS)     ! fraction of solar constant in each band
      DATA FSCIEB / .429760, .326158, .180608, .033474 /
C     !  First value has 3% knocked off for Rayleigh scattering.
      REAL TTEC(NGASES,NTRANS+2)       ! Offsets & multipliers for use
C     ! finding the place in (D)TRTAB, and constant for finding
C     ! transmissivity for very small pathlengths, all used in SWTRAN.
      DATA ((TTEC(GAS,INIT), GAS=1, NGASES), INIT=NTRANS+1,NTRANS+2)
     &        / 23., 57., 11.4, 2.17145, 4.3429, .86858 /
C     !   Last three values are 5,10,2/log(10).   Why not put in as such
      REAL GREYNT,                     ! Net grey factor in current beam
     &     SIF,                        ! Surface incoming flux in 1 band
     &     MINPTH,                     ! Mininum pathlength catered for
C     ! in the look-up table for a particular absorber.
     &     DIFFTR,                     ! See "DO 6" loop
     &     NWDFCN,                     ! Contributions to NEWDIF due to
     &     NWDFLY,                     !     convective and layer cloud
     &     GRDFCL                      !  Temporary store for the grey
C     ! factor for diffuse light incident on the top of the current
C     ! layer in the clear part of the layer boundary, rather than where
C     ! convective cloud crosses the boundary, if it does.
      INTEGER LSTCLR,                  ! Lowest always clear layer, and
     &     FSTCLD,                     !    highest possibly cloudy one
     &     DIRECT,                     ! Subscript for PATH & GREY
     &     OFFSET,                     ! Index for cloud albedo*amount
C     ! diagnostics, which SWMAST returns (potentially) compressed,
C     ! allowing just the bands or level-and-band combinations wanted to
C     ! have space allocated by STASH and be set here.  Bands are in
C     ! standard order, and, following the UM standard, multi-level
C     ! data has the different levels for each band or other
C     ! "pseudo-dimension" together, running up from the surface.
     &     NBEAMS,                     !  Number of beams to deal with
C     ! on current pass through "DO 40" loop over potentially cloudy
C     ! layers - first just the direct beam, then a diffuse one too.
     &     CONCLD                      ! Subscript in GREY of the factor
C     ! for the beam inside convective cloud
      PARAMETER ( CONCLD = 3, DIRECT = 1 )
CL
CL    ! Section 1
CL      ~~~~~~~~~
CL    ! Various initialization etc. - setting up constants to address
CL    ! arrays, array TTEC to pass to SWTRAN, arrays of scaled
CL    ! pathlengths by CALLing SWPTSC, cloud optical properties using
CL    ! SWCLOP and thence modified surface albedo by CALLing SWMSAL.
CL
      LSTCLR = NLEVS - NCLDS
      FSTCLD = LSTCLR + 1
C
      DO 11 GAS=1, NGASES
       MINPTH = EXP ( (1.-TTEC(GAS,NTRANS+1)) / TTEC (GAS,NTRANS+2) )
       DO 11 INIT=1, NTRANS
        TTEC(GAS,INIT) = ( 1.-TRTAB(1,INIT,GAS,1) )/ MINPTH
   11 CONTINUE
C
CL    !  CALL SWPTSC to set up DPATH from the water vapour, carbon
CL    !  dioxide and ozone mixing ratios, and pressure information for
CL    !  the pressure scaling.
C
Cfpp$ Expand
      CALL SWPTSC (H2O, CO2, O3, PSTAR, AB, BB,
     &     L2,                                                          
     &     NLEVS, NWET, NOZONE,         L1, DPATH)
C
CL    !  Next, set up cloud-related quantities
C
      IF ( NCLDS.GT.0 ) THEN
C
CL       !  First CALL to SWCLOP is for layer cloud.
C        !   Cloud water pathlength (mass per unit area), effective
C        !   radius and solar zenith angle are used to calculate their
C        !   optical properties.
C
Cfpp$ Expand
         CALL SWCLOP (LCCWP, LRE, COSZ, L1, L2, NCLDS, REF, CTRANS)
C
CL       !  Multiplication by cloud cover gives the optical properties
CL       !   averaged over the grid-box (as far as layer cloud goes).
C
         DO 12 DIRDIF=1, 2
          DO 12 LEVEL=FSTCLD, NLEVS
           DO 12 BAND=1, NBANDS
Cfpp$       Select(CONCUR)
            DO 12 J=1, L2
             REF(J,BAND,LEVEL,DIRDIF) =
     &         REF(J,BAND,LEVEL,DIRDIF) * LCA(J,LEVEL)
   12    CONTINUE
C
CL       !  SWCLOP is then CALLed for convective cloud.
C        !  There is only one layer to deal with.
Cfpp$ Expand
         CALL SWCLOP (CCCWP, CRE, COSZ, L1, L2, 1, CCREF, CCTRANS)
C
CL       !  Then the CALL to SWMSAL.
C        !   This must come before the convective and layer cloud
C        !   reflectivities are combined, as the combination is done for
C        !   light coming down, and it would be different for light
C        !   coming up after surface reflection where there was non-zero
C        !   convective cloud more than one layer thick.
C
Cfpp$ Expand
         CALL SWMSAL (TSA, REF(1,1,FSTCLD,2), LCA, CCREF(1,1,2), CCA,
     &     CCB, LSTCLR,
     &     L2,                                                          
     &     L1, NBANDS, NCLDS,                          MODSA)
C
C        ! Diagnose cloud amounts * albedos if they are wanted
C
         IF ( LCAARO ) THEN
           OFFSET = 1
           DO BAND=1, NBANDS
             DO LEVEL=NLEVS, FSTCLD, -1
               IF ( LCAARL(NLEVS+1-LEVEL) .AND. LCAARB(BAND) ) THEN
Cfpp$            Select(CONCUR)
                 DO J=1, L2
                   LCAAR(J,OFFSET) = REF(J,BAND,LEVEL,1)
                 ENDDO
                 OFFSET = OFFSET + 1
               ENDIF
             ENDDO
           ENDDO
         ENDIF
         IF ( LCAAFO ) THEN
           OFFSET = 1
           DO BAND=1, NBANDS
             DO LEVEL=NLEVS, FSTCLD, -1
               IF ( LCAAFL(NLEVS+1-LEVEL) .AND. LCAAFB(BAND) ) THEN
Cfpp$            Select(CONCUR)
                 DO J=1, L2
                   LCAAF(J,OFFSET) = REF(J,BAND,LEVEL,2)
                 ENDDO
                 OFFSET = OFFSET + 1
               ENDIF
             ENDDO
           ENDDO
         ENDIF
         IF ( CCAARO ) THEN
           OFFSET = 1
           DO BAND=1, NBANDS
             IF ( CCAARB(BAND) ) THEN
Cfpp$          Select(CONCUR)
               DO J=1, L2
                 CCAAR(J,OFFSET) = CCREF(J,BAND,1) * CCA(J)
               ENDDO
               OFFSET = OFFSET + 1
             ENDIF
           ENDDO
         ENDIF
         IF ( CCAAFO ) THEN
           OFFSET = 1
           DO BAND=1, NBANDS
             IF ( CCAAFB(BAND) ) THEN
Cfpp$          Select(CONCUR)
               DO J=1, L2
                 CCAAF(J,OFFSET) = CCREF(J,BAND,2) * CCA(J)
               ENDDO
               OFFSET = OFFSET + 1
             ENDIF
           ENDDO
         ENDIF
C
CL       !  Finally combine the convective and layer cloud
CL       !   reflectivities to get the effective layer mean reflectivity
CL       !   Recall that the layer cloud cover and water path are
CL       !   deemed to describe the fraction of the grid-box outside
CL       !   the convective cloud.
C
         DO 13 DIRDIF=1, 2
          DO 13 BAND=1, NBANDS
Cfpp$      Select(CONCUR)
           DO 13 J=1, L2
            REF(J,BAND,CCT(J),DIRDIF) = CCREF(J,BAND,DIRDIF) * CCA(J) +
     &                     REF(J,BAND,CCT(J),DIRDIF) * ( 1. - CCA(J) )
   13    CONTINUE
C
       ELSE
C
CL       ! If there are no clouds to be treated, just copy the clear-sky
CL       !  surface albedos to be used as modified surface albedos, and
CL       !  leave the rest to the DO loop bounds.
CL       !  THIS MAY OR MAY NOT WORK - untested 29/10/90
C
         DO 14 DIRDIF=1, 2
         DO 14 BAND=1, NBANDS
Cfpp$     Select(CONCUR)
          DO 14 J=1, L2
           MODSA(J,BAND,DIRDIF) = TSA(J,BAND,DIRDIF)
   14    CONTINUE
      ENDIF
C
CL    !  Last bit of Section 1 sets up the magnification factor for the
CL    !    direct beam.
C
      DO 15 J=1, L2
       DIRFAC(J) = HBYAP1 / SQRT ( COSZ(J)**2 + HBYAX2 )
   15 CONTINUE
C
CL    !  Section 2
CL       ~~~~~~~~~
CL    !  Calculate the flux at the top of the atmosphere taking account
CL    !  of Rayleigh scattering, and initialize parts of FLUX and PATH.
C
C     !  Rayleigh scattering is represented by simply reflecting RAYSCA
C     !  of the incoming sunlight (in the shortest wavelength band)
C     !  before any interaction with the atmosphere.  This is done by
C     !  subtracting RAYSCA from FSCIEB(1) before inserting the value
C     !  in the code, and by the following code for the top of the
C     !  model, where FSCIEB is not automatically picked up via SWTRAN.
C     !
C     !  Obviously, if this is to be changed, consistency must be
C                                                           maintained.
      DO 21 J=1, L2
       FLUX(J,0) = ONELRS
   21 CONTINUE
C
      IF ( CSOSON ) THEN
        DO J=1, L2
          CSOSDI(J) = RAYSCA
        ENDDO
      ENDIF
      IF ( NSS1ON ) THEN
        DO J=1, L2
          NSSB1(J) = 0.
        ENDDO
      ENDIF
C
      DO 23 LEVEL=1, NLEVS
Cfpp$  Select(CONCUR)
       DO 23 J=1, L2
        FLUX(J,LEVEL) = 0.
   23 CONTINUE
      DO 24 GAS=1, NGASES
Cfpp$  Select(CONCUR)
       DO 24 J=1, L2
        PATH(J,GAS,DIRECT) = DIRFAC(J) * DPATH(J,GAS,1)
   24 CONTINUE
C
CL    !  Section 3
CL       ~~~~~~~~~
CL    !  For the remaining layers above the levels where cloud may occur
CL    !   calculations are very simple - just loop down accumulating the
CL    !   gaseous pathlengths for the direct beam, calculating gaseous
CL    !   transmissivities from them and adding these in without having
CL    !   to use grey factors.
C
      DO 3 LEVEL=2, LSTCLR
Cfpp$  Expand
       CALL SWTRAN (PATH, TTEC, TRTAB, TRTAB(1,1,1,2),
     &     L2,                                                          
     &     GTRANS)
       DO 32 BAND=1, NBANDS
Cfpp$   Select(CONCUR)
        DO 32 J=1, L2
         FLUX(J,LEVEL-1) = FLUX(J,LEVEL-1) + GTRANS(J,BAND)
   32  CONTINUE
       DO 34 GAS=1, NGASES
Cfpp$   Select(CONCUR)
        DO 34 J=1, L2
         PATH(J,GAS,DIRECT) =
     &         PATH(J,GAS,DIRECT) + DIRFAC(J) * DPATH(J,GAS,LEVEL)
   34  CONTINUE
    3 CONTINUE
C
CL    !  And, before starting the loop over cloudy layers, the code must
CL    !   initialize all the grey factors, and the diffuse pathlengths.
C
      DO 36 BAND=1, NBANDS
Cfpp$  Select(CONCUR)
       DO 36 J=1, L2
        GREY(J,BAND,DIRECT) = 1.
        GREY(J,BAND,CONCLD) = 0.
        GREY(J,BAND,2) = 0.
   36 CONTINUE
C
      DO 38 GAS=1, NGASES
       DO 38 J=1, L2
        PATH(J,GAS,2) = PATH(J,GAS,1)
   38 CONTINUE
C
CL    !  Section 4
CL       ~~~~~~~~~
CL    !  Start the loop over cloudy levels, which has to be long and
CL    !  complex to allow for all the permitted interactions, and so is
CL    !  split into three sections.  Section 4 calculates the flux terms
CL    !  for downward light at the top of layer LEVEL.
C
      NBEAMS = 1
C
      DO 4 LEVEL=FSTCLD, NLEVS
C      !  Inside the loop over the levels for which we are finding the
C      !   flux, loop over the "beams" impinging on the level from above
C      !   - i.e. the categories of light whose histories are different
C      !   enough for us to keep separate pathlengths.
       DO 40 BEAM=1, NBEAMS
        DIRDIF = BEAM
Cfpp$   Expand
        CALL SWTRAN (PATH(1,1,BEAM), TTEC, TRTAB, TRTAB(1,1,1,2),
     &     L2,                                                          
     &     GTRANS)
C       !  For each beam and band, add in its gaseous transmissivity
C       !   multiplied by the right grey factor.  Note that the grey
C       !   factors are defined for the total amount of light impinging
C       !   on the layer boundary from above, so that some manipulation
C       !   is needed to get GREYNT from the GREYs.
C
        DO 41 BAND=1, NBANDS
Cfpp$    Select(CONCUR)
         DO 41 J=1, L2
          GREYNT = (1.-REF(J,BAND,LEVEL,DIRDIF)) * GREY(J,BAND,BEAM)
C         ! The above line would be all that was needed if all the light
C         !   hitting the layer were liable to reflection, but some
C         !   diffuse light may cross it inside the convective cloud:
          IF ( BEAM .EQ. 2 )
     &       GREYNT = GREYNT + GREY(J,BAND,CONCLD) * REF(J,BAND,LEVEL,2)
          FLUX(J,LEVEL-1) = FLUX(J,LEVEL-1) + GTRANS(J,BAND) * GREYNT
   41   CONTINUE
C
CL      !  Section 5
CL        ~~~~~~~~~
CL      !  This section calculates the flux due to reflection from the
CL      !   clouds in layer LEVEL through higher layers up to space.
CL      !   Recall that this light is assumed to pass to space without
CL      !   interaction with clouds, only gaseous absorption occurring.
C
CL      !  First, set up RFGREY, the grey factor for the current beam:
        DO 53 BAND=1, NBANDS
Cfpp$    Select(CONCUR)
         DO 53 J=1, L2
C         ! (as in the "DO 43" loop, a little manipulation of the GREYs)
          GREYNT = GREY(J,BAND,BEAM)
          IF ( BEAM .EQ. 2 ) GREYNT = GREYNT - GREY(J,BAND,CONCLD)
          RFGREY(J,BAND) = REF(J,BAND,LEVEL,DIRDIF) * GREYNT
   53   CONTINUE
CL      !  and initialize RFPATH, its pathlength
        DO 55 GAS=1, NGASES
Cfpp$    Select(CONCUR)
         DO 55 J=1, L2
          RFPATH(J,GAS) =
     &      PATH(J,GAS,BEAM) + DIFFAC * DPATH(J,GAS,LEVEL-1)
   55   CONTINUE
CL      !  and then loop up to the top of the atmosphere, putting in
CL      !   each upward flux term with no further calculation of grey
CL      !   terms, just finding transmissivities:
        DO 50 LEVEL2=LEVEL-2, 0, -1
Cfpp$    Expand
         CALL SWTRAN (RFPATH, TTEC, TRTAB, TRTAB(1,1,1,2),
     &     L2,                                                          
     &     GTRANS)
         DO 57 BAND=1, NBANDS
Cfpp$     Select(CONCUR)
          DO 57 J=1, L2
           FLUX(J,LEVEL2) = FLUX(J,LEVEL2) -
     &                           RFGREY(J,BAND) * GTRANS(J,BAND)
   57    CONTINUE
CL      !  and incrementing RFPATH as each layer is crossed:
         IF (LEVEL2.NE.0) THEN
           DO 59 GAS=1, NGASES
Cfpp$       Select(CONCUR)
            DO 59 J=1, L2
             RFPATH(J,GAS) =
     &                     RFPATH(J,GAS) + DIFFAC * DPATH(J,GAS,LEVEL2)
   59      CONTINUE
         ENDIF
   50   CONTINUE
C
   40  CONTINUE                             ! End of the loop over BEAM
C
C      !  The first time through "DO 40" NBEAMS was 1: thereafter 2.
       NBEAMS = 2
C
C
CL     !  Section 6
CL        ~~~~~~~~~
CL     ! Sections 4 and 5 dealt with transmission to and reflection from
CL     !   layer LEVEL: Section 6 prepares to deal with transmission
CL     !   through it by finding the effects of the cloud(s) whose tops
CL     !   are in it on the grey terms for the next layer boundary, and
CL     !   setting up the pathlengths to be used at the next layer
CL                                                             boundary.
CL
CL     !  The code is a little involved, but the physics being
CL     !   implemented is straightforward enough.
CL     !  Without convective cloud, all that would happen would be that
CL     !   direct light would be transmitted through the cloud-free area
CL     !   as direct light and through the cloud as diffuse light, and
CL     !   both direct and diffuse light would be attenuated in the
CL     !   cloudy area according to the appropriate transmissivity.
CL     !  At convective cloud top both clouds' fractional cover and
CL     !   transmissivities must be considered, and the amount of light
CL     !   going into the convective cloud acounted for.  The latter
CL     !   will no longer be affected by clouds till it comes out of the
CL     !   convective cloud's base, when it can be combined with the
CL     !   rest of the diffuse light.  The convective cloud need not be
CL     !   explicitly considered to calculate the change in the other
CL     !   grey factors "beside" it (i.e. when considering layer
CL     !   boundaries which it crosses), as the layer cloud fractions
CL     !   are there the fractions in the convective-cloud-free area.
C
       DO 60 J=1, L2
        NEWDIF(J) = 0.
   60  CONTINUE
C
C      !  Rather than a line-by-line explanation of the code in this
C      !   loop, the physical explanation above and the definitions of
C      !   each quantity seem most useful.
C      !  DIFFTR is an effective grey transmissivity of the layer to
C      !   diffuse light: the fraction of that diffuse light to impinge
C      !   on the layer not in convective cloud which is transmitted not
C      !   in convective cloud, ignoring gaseous absorption.
C
       DO 6 BAND=1, NBANDS
Cfpp$   Select(CONCUR)
        DO 61 J=1, L2
         GRDFCL = GREY(J,BAND,2) - GREY(J,BAND,3)
         DIFFTR = 1. - LCA(J,LEVEL) * ( 1. - CTRANS(J,BAND,LEVEL,2) )
         NWDFCN = CCA(J) * CCTRANS(J,BAND,1) * GREY(J,BAND,DIRECT)
C        !  NWDFCN is only used if CCT=LEVEL, but is set unconditionally
C        !                                              for fpp's sake.
         IF ( CCT(J) .EQ. LEVEL ) THEN
           GREY(J,BAND,3) = NWDFCN + CCA(J) * CCTRANS(J,BAND,2) * GRDFCL
           DIFFTR = DIFFTR * ( 1. - CCA(J) )
           GREY(J,BAND,DIRECT) = GREY(J,BAND,DIRECT) * ( 1. - CCA(J) )
         ENDIF
         NWDFLY = LCA(J,LEVEL)*CTRANS(J,BAND,LEVEL,1)*GREY(J,BAND,1)
         GREY(J,BAND,2) = NWDFLY + GRDFCL * DIFFTR + GREY(J,BAND,3)
         IF ( BAND .GE. FSWTBD ) THEN
           NEWDIF(J) = NEWDIF(J) + NWDFLY * FSCIEB(BAND)
           IF ( CCT(J) .EQ. LEVEL )
     &                 NEWDIF(J) = NEWDIF(J) + NWDFCN * FSCIEB(BAND)
         ENDIF
         GREY(J,BAND,DIRECT) = GREY(J,BAND,DIRECT)*( 1. - LCA(J,LEVEL) )
         IF ( CCB(J) .EQ. LEVEL ) THEN
           GREY(J,BAND,CONCLD) = 0.
         ENDIF
   61   CONTINUE
    6  CONTINUE
C
C      !  NEWDIF now contains the increase in the sum of the diffuse
C      !   grey factors for the wet bands due to passing through the
C      !   layer: normalize by the new sum to get the fraction of the
C      !   diffuse light which has just become diffuse and so has a
C      !   pathlength so far equal to the direct beam's.
       DO 62 J=1, L2
        TOWTGR(J) = GREY(J,FSWTBD,2) * FSCIEB(FSWTBD)
   62  CONTINUE
       DO 63 BAND=FSWTBD+1, NBANDS
        DO 63 J=1, L2
         TOWTGR(J) = TOWTGR(J) + GREY(J,BAND,2) * FSCIEB(BAND)
   63  CONTINUE
       DO 65 J=1, L2
        IF ( TOWTGR(J) .NE. 0. )  NEWDIF(J) = NEWDIF(J) / TOWTGR(J)
   65  CONTINUE
C      !  Set up pathlengths for the bottom of the layer LEVEL.
       DO 66 GAS=1, NGASES
Cfpp$   Select(CONCUR)
        DO 66 J=1, L2
C        !  The diffuse beam's pathlengths are the average, according to
C        !   the split given by NEWDIF, of their old values and the
C        !   values for the direct beam, plus the contribution due to
C        !   crossing layer LEVEL diffusely:
         PATH(J,GAS,2) = NEWDIF(J) * PATH(J,GAS,DIRECT) +
     & ( 1. - NEWDIF(J) ) * PATH(J,GAS,2) + DIFFAC * DPATH(J,GAS,LEVEL)
C        !   while the direct beam again adds on the layer pathlengths
C        !   multiplied by the direct beam magnification factor:
         PATH(J,GAS,DIRECT) =
     &         PATH(J,GAS,DIRECT) + DIRFAC(J) * DPATH(J,GAS,LEVEL)
   66  CONTINUE
C
    4 CONTINUE                       !  End of the outer loop over LEVEL
C
CL    !  Section 8
CL       ~~~~~~~~~
CL    !  This section calculates the surface flux and the effects of the
CL    !   light reflected from the surface.  It is thus similar to
CL    !   Sections 4 and 5, but simpler in that the surface does not
CL    !   have fractional cover, transmissivity or a beam crossing it
CL     !   inside convective cloud.  However this physical simplicity
CL     !   is partly offset by the fact that extra outputs are
CL     !   calculated - DSFLUX and possibly surface diagnostics.
C
       DO J=1, L2
         DSFLUX(J) = 0.
       ENDDO
       IF ( TDSSON ) THEN
         DO J=1, L2
           TDSS(J) = 0.
         ENDDO
       ENDIF
C
C     !  First time round, use direct-light albedos:
      DIRDIF=1
C     ! The "DO 8" loop is similar to the "DO 40" loop.
      DO 8 BEAM=DIRECT, NBEAMS
Cfpp$  Expand
       CALL SWTRAN (PATH(1,1,BEAM), TTEC, TRTAB, TRTAB(1,1,1,2),
     &     L2,                                                          
     &     GTRANS)
C      ! The "DO 81" loops are similar to the "DO 41" loops, but rather
C      !   simpler, as there cannot be any convective cloud crossing
C      !   this layer boundary.
       DO 81 BAND=1, NBANDS
Cfpp$   Select(CONCUR)
        DO 81 J=1, L2
         SIF = GTRANS(J,BAND) * GREY(J,BAND,BEAM)
         FLUX(J,NLEVS) = FLUX(J,NLEVS) + SIF * (1.-MODSA(J,BAND,DIRDIF))
         DSFLUX(J)     = DSFLUX(J)     + SIF *
     &      (    1.  -  DTSA(J,BAND,DIRDIF)  +
     & ( TSA(J,BAND,DIRDIF) - MODSA(J,BAND,DIRDIF) ) *
     &        ( 1. - DTSA(J,BAND,2) ) / ( 1. - TSA(J,BAND,2) )      )
   81  CONTINUE
       IF ( TDSSON ) THEN
         IF ( BEAM .EQ. DIRECT ) THEN
            DO 811 BAND=1, NBANDS
              DO J=1, L2
                TDSS(J) = TDSS(J) + GTRANS(J,BAND) * GREY(J,BAND,BEAM) *
     &  ( 1. + (TSA(J,BAND,1)-MODSA(J,BAND,1)) / (1.-TSA(J,BAND,2)) )
              ENDDO
  811       CONTINUE
          ELSE              ! Diffuse beam
            DO 818 BAND=1, NBANDS
              DO J=1, L2
                TDSS(J) = TDSS(J) + GTRANS(J,BAND) * GREY(J,BAND,BEAM) *
     &          ( 1.-MODSA(J,BAND,DIRDIF) ) / ( 1.-TSA(J,BAND,DIRDIF) )
              ENDDO
  818       CONTINUE
         ENDIF
       ENDIF
       IF ( NSS1ON ) THEN
Cfpp$    Select(CONCUR)
         DO J=1, L2
           NSSB1(J) = NSSB1(J) +
     &          GTRANS(J,1) * GREY(J,1,BEAM) *
     &      (    1.  -  DTSA(J,1,DIRDIF)  +
     & ( TSA(J,1,DIRDIF) - MODSA(J,1,DIRDIF) ) *
     &        ( 1. - DTSA(J,1,2) ) / ( 1. - TSA(J,1,2) )          )
         ENDDO
       ENDIF
       IF ( CSSSDO .AND. BEAM .EQ. DIRECT ) THEN
         DO J=1, L2
           CSSSD(J) = GTRANS(J,1)
         ENDDO
         DO 812 BAND=2, NBANDS
Cfpp$      Select(CONCUR)
           DO J=1, L2
             CSSSD(J) = CSSSD(J) + GTRANS(J,BAND)
           ENDDO
  812    CONTINUE
       ENDIF
       IF ( CSSSUO .AND. BEAM .EQ. DIRECT ) THEN
         DO J=1, L2
           CSSSU(J) = GTRANS(J,1) * TSA(J,1,1)
         ENDDO
         DO 813 BAND=2, NBANDS
Cfpp$      Select(CONCUR)
           DO J=1, L2
             CSSSU(J) = CSSSU(J) + GTRANS(J,BAND) * TSA(J,BAND,1)
           ENDDO
  813    CONTINUE
       ENDIF
CL     !  As in "DO 53", set up RFGREY, grey factor for the current beam
       DO 83 BAND=1, NBANDS
Cfpp$   Select(CONCUR)
        DO 83 J=1, L2
         RFGREY(J,BAND) = MODSA(J,BAND,DIRDIF) * GREY(J,BAND,BEAM)
   83  CONTINUE
CL     !  and, as in "DO 55", initialize RFPATH, its pathlength
       DO 85 GAS=1, NGASES
Cfpp$   Select(CONCUR)
        DO 85 J=1, L2
         RFPATH(J,GAS) =
     &     PATH(J,GAS,BEAM) + DIFFAC * DPATH(J,GAS,NLEVS)
   85  CONTINUE
CL     !  and then, as in "DO 50" & "DO 57", loop up to the top of the
CL     !  atmosphere, putting in each upward flux term with no further
CL     !   calculation of grey terms, just finding transmissivities:
       DO 80 LEVEL2=NLEVS-1, 0, -1
Cfpp$   Expand
         CALL SWTRAN (RFPATH, TTEC, TRTAB, TRTAB(1,1,1,2),
     &     L2,                        
     &     GTRANS)
        DO 87 BAND=1, NBANDS
Cfpp$    Select(CONCUR)
         DO 87 J=1, L2
          FLUX(J,LEVEL2) = FLUX(J,LEVEL2) -
     &                            RFGREY(J,BAND) * GTRANS(J,BAND)
   87   CONTINUE
C       ! and, as in "DO 59", incrementing RFPATH for each layer crossed
        IF ( LEVEL2 .NE. 0 ) THEN
          DO 89 GAS=1, NGASES
Cfpp$      Select(CONCUR)
           DO 89 J=1, L2
            RFPATH(J,GAS) =
     &           RFPATH(J,GAS) + DIFFAC * DPATH(J,GAS,LEVEL2)
   89     CONTINUE
        ENDIF
   80  CONTINUE
      IF ( CSOSON .AND. BEAM .EQ. DIRECT ) THEN
        DO 92 BAND=1, NBANDS
          DO 93 J=1, L2
            CSOSDI(J) = CSOSDI(J) + GTRANS(J,BAND) * TSA(J,BAND,1)
   93     CONTINUE
   92   CONTINUE
      ENDIF
C     !  After the first time round, use diffuse-light albedos:
      DIRDIF = 2
    8 CONTINUE
C
      RETURN
      END
