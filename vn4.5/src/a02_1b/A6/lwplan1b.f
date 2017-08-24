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
CLL Subroutine  LWPLAN  ----------------------------------------------
CLL
CLL Purpose :
CLL           The version of routine LWPLAN used in
CLL  version 1B (gaseous effects treated as Morcrette et al, 1986) of
CLL  the UM LW code.
CLL  Version 3 using the ECMWF quintic fits for the ECMWF LW bands, part
CLL  of the alternative code giving ECMWF-like treatment of LW gaseous
CLL  transmissivities.  The changes are not great: temperatures are now
CLL  normalized to a reference temperature before using the fits, the
CLL  fits are quintic rather than cubic, the constants (in *COMDECK
CLL  LWPFCO) are different, and the indentation is changed.
CLL  Version 3 of LWPLAN was set up from version 2.2 to be part of
CLL  version 1B (ECMWF-like gaseous transmissivities) of the LW from
CLL  release 2.7 of the UM.                 William Ingram 22 June 1992
CLL  It calculates Planck ("black-body") fluxes in each longwave band as
CLL  a quintic function of atmospheric (layer centre and boundary) and
CLL  surface temperatures and returns their differences across
CLL  half-layers (and optionally the surface black-body flux) for use by
CLL  by LWMAST constructing longwave fluxes.
CLL
CLL     Author: William Ingram
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL Programming standard :
CLL  It conforms with standard A of version 3 (07/9/90) of UMDP 4, and
CLL  contains no 8X-deprecated features.
CLL  The code is standard FORTRAN 77 except for having ! comments.
CLL
CLL Logical components covered : P232, D23
CLL  P232 (longwave radiation) D23 (radiation diagnostics).
CLL
CLL Project task : P23
CLL
CLL External documentation:      UMDP 23.
CLL
CLLEND -----------------------------------------------------------------
C*L
      SUBROUTINE LWPLAN (TAC, PEXNER, PSTAR, AKH, BKH, TSTAR, AICE,
     &     SFUP, SFUPON,
     &     L2, NLEVS, L1,                     SEAFX,  BHDB, THDB)
C*
      INTEGER NBANDS          ! Number of spectral bands in the longwave
      PARAMETER (NBANDS=6)    ! This run uses the standard set of
!  "ECMWF" bands as described by Morcrette et al (J.-J. Morcrette,
!   L.D. Smith & Y. Fouquart, 1986, Beitr. Phys. Atmosph., 59, 455-469)
C*L
      INTEGER!, INTENT (IN) ::
     &     L1,                      !  First dimension of input arrays
     &     L2,                      !  Number of points to be treated
     &     NLEVS                    !  Number of levels
      REAL!, INTENT(IN) ::
     &     TAC(L1,NLEVS),           !  Atmospheric temperatures at layer
C                                   !                           centres
     &     PEXNER(L1,NLEVS+1),      !  Exner function @ layer boundaries
     &     PSTAR(L1),               !  Surface pressure
     &     AKH(NLEVS+1),            !  As at layer boundaries
     &     BKH(NLEVS+1),            !  Bs at layer boundaries
C     !  The above four are used to get temperature at layer boundaries
     &     TSTAR(L1),               !  Surface temperature
     &     AICE(L1)                 !  Sea-ice fraction
      LOGICAL!, INTENT(IN) ::
     &     SFUPON
      REAL!, INTENT(OUT) ::
     &     BHDB(L2,NLEVS,NBANDS),
     &     THDB(L2,NLEVS,NBANDS),
C     !  BHDB holds the difference of the Planck functions across the
C     !  bottom half of each layer and THDB across the top half.
     &     SEAFX(L1)                ! Difference between the open-sea
C     ! and sea-ice upward longwave fluxes at the surface.
     &     , SFUP(L1)               ! Upward surface flux
C*
CL    ! LWPLAN has no EXTERNAL calls and no dynamically allocated arrays
C
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_R_CP-------------------------------------
C R IS GAS CONSTANT FOR DRY AIR
C CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
C PREF IS REFERENCE SURFACE PRESSURE
      REAL R,CP,KAPPA,PREF  ! PREFTOKI

      PARAMETER(R=287.05,
     &          CP=1005.,
     &          KAPPA=R/CP,
     &          PREF=100000.)
C*----------------------------------------------------------------------

      REAL TRPLAN                   !  Reference temperature to which
C     !       temperatures are normalized before using the quintic fits.
      PARAMETER ( TRPLAN = 250. )
      REAL PFCO (6,NBANDS)          !  Coefficients of the quintics
      DATA PFCO /
     &    0.46430621E+02,   0.12928299E+03,   0.20732648E+03,
     &    0.31398411E+03,   0.18373177E+03,  -0.11412303E+03,
     &    0.73604774E+02,   0.27887914E+03,   0.27076947E+03,
     &   -0.57322111E+02,  -0.64742459E+02,   0.87238280E+02,
     &    0.37050866E+02,   0.20498759E+03,   0.37558029E+03,
     &    0.17401171E+03,  -0.13350302E+03,  -0.37651795E+02,
     &    0.14930141E+02,   0.89161160E+02,   0.17793062E+03,
     &    0.93433860E+02,  -0.70646020E+02,  -0.26373150E+02,
     &    0.40386780E+02,   0.10855270E+03,   0.50755010E+02,
     &   -0.31496190E+02,   0.12791300E+00,   0.18017770E+01,
     &    0.90811926E+01,   0.75073923E+02,   0.24654438E+03,
     &    0.39332612E+03,   0.29385281E+03,   0.89107921E+02 /
C
C
      INTEGER BAND, LEVEL, J        !  Loopers over band, level & point
      REAL TAB,                     !  Layer boundary temperature
     &     WTL, WTU,                !  and weights for its calculation
     &     PL,                      ! Pressure at upper bdry of current
     &     PU,                      !  layer and next one up
     &     PLM1,                    ! Pressure at lower bdry of layer
     &     BBC,                     !  Black-body flux at a model layer
     &     BBB,                     !             centre or boundary
     &     BBTFS                    !             or at TFS
      REAL TN,                      ! Simplify the calculation of
     &     PLANCK,                  !   Planck fluxes with statement
     &     T                        !   functions
C*L------------------ COMDECK P_EXNERC ---------------------------------
C statement function to define exner pressure at model levels
C from exner pressures at model half-levels
      REAL P_EXNER_C
      REAL R_P_EXNER_C
      REAL P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM
C consistent with geopotential see eqn 26 DOC PAPER 10
      P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)/
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )
      R_P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )/
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)

C*------------------- --------------------------------------------------

      TN(T) = T / TRPLAN - 1.
      PLANCK(T,BAND) =
     &      ( ( ( ( PFCO(6,BAND) * TN(T) + PFCO(5,BAND) ) * TN(T) +
     &              PFCO(4,BAND) ) * TN(T) + PFCO(3,BAND) ) * TN(T) +
     &              PFCO(2,BAND) ) * TN(T) + PFCO(1,BAND)
C
CL    !  Section 1
CL    !  ~~~~~~~~~
CL    !  First find Planck function for surface temperature:
C
      DO 1 BAND=1, NBANDS
        BBTFS = PLANCK(TFS,BAND)
Cfpp$   Select(CONCUR)
        DO 11 J=1, L2
          IF ( AICE(J) .GT. 0. ) THEN
             TAB = ( TSTAR(J) + (AICE(J)-1.) * TFS ) / AICE(J)
           ELSE
             TAB = TSTAR(J)
          ENDIF
          BBB = PLANCK(TAB,BAND)
          SEAFX(J) = SEAFX(J) + BBTFS - BBB
          IF ( AICE(J) .GT. 0. )
     &                BBB = AICE(J) * BBB + (1.-AICE(J)) * BBTFS
          BHDB(J,1,BAND) = BBB
          IF ( SFUPON ) SFUP(J) = SFUP(J) + BBB
   11   CONTINUE
    1 CONTINUE
C
CL    !  Section 2
CL    !  ~~~~~~~~~
CL    !  Loop over the model layers finding the dB across each half
CL    !     of each one :
C
      DO 2 LEVEL=1, NLEVS-1
C       !
C       !  First must find TAB for the top of the layer (same for each
C       !                    band) and store it in spare space in THDB:
Cfpp$   Select(CONCUR)
        DO 20 J=1, L2
          PU = PSTAR(J)*BKH(LEVEL+2) + AKH(LEVEL+2)
          PL = PSTAR(J)*BKH(LEVEL+1) + AKH(LEVEL+1)
          PLM1 = PSTAR(J)*BKH(LEVEL) + AKH(LEVEL)
          WTL = TAC(J,LEVEL+1)*( PEXNER(J,LEVEL+1) /
     &        P_EXNER_C(PEXNER(J,LEVEL+2),PEXNER(J,LEVEL+1),PU,PL,KAPPA)
     &         - 1.0 )
          WTU = TAC(J,LEVEL) * ( PEXNER(J,LEVEL) /
     &        P_EXNER_C(PEXNER(J,LEVEL+1),PEXNER(J,LEVEL),PL,PLM1,KAPPA)
     &         - 1.0 )
          TAB =
     &   ( WTL * TAC(J,LEVEL) + WTU * TAC(J,LEVEL+1) ) / ( WTL + WTU )
          THDB(J,NLEVS,1) = TAB
   20   CONTINUE
C       !
C       !  This loop finds the black-body fluxes at the middle & top of
C       !  the layer and sets BHDB & THDB for the layer, using the
C       !  black-body flux for the bottom of the layer which is already
C       !  in BHDB and leaving the value at the top in the next level
C       !                                                  up of BHDB :
        DO 21 BAND=1, NBANDS
Cfpp$   Select(CONCUR)
          DO 22 J=1, L2
            BBC = PLANCK(TAC(J,LEVEL),BAND)
            BHDB(J,LEVEL,BAND) = BBC - BHDB(J,LEVEL,BAND)
            TAB = THDB(J,NLEVS,1)
            BBB = PLANCK(TAB,BAND)
            THDB(J,LEVEL,BAND) = BBB - BBC
            BHDB(J,LEVEL+1,BAND) = BBB
   22     CONTINUE
   21   CONTINUE
    2 CONTINUE
C
CL    !  Section 3
CL    !  ~~~~~~~~~
CL    !  Finally, the top layer is treated slightly differently because
CL    !  the black-body flux at the top is zero:
C
      DO 3 BAND=1, NBANDS
Cfpp$   Select(CONCUR)
        DO 30 J=1, L2
          BBC = PLANCK(TAC(J,NLEVS),BAND)
          BHDB(J,NLEVS,BAND) = BBC - BHDB(J,NLEVS,BAND)
          THDB(J,NLEVS,BAND) = - BBC
   30   CONTINUE
    3 CONTINUE
C
      RETURN
      END
