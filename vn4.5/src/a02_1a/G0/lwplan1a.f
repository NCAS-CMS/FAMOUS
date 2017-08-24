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
CLL  Subroutine LWPLAN ----------------------------------------------
CLL
CLL  This version of routine LWPLAN used in
CLL  version 1A (gaseous effects treated as Slingo & Wilderspin, 1986)
CLL  of the UM LW code.
CLL  It calculates Planck ("black-body") fluxes in each longwave band as
CLL  a cubic function of atmospheric (layer centre and boundary) and
CLL  surface temperatures and returns their differences across
CLL  half-layers (and optionally the surface black-body flux) for use by
CLL  by LWMAST constructing longwave fluxes.
CLL
CLL    Author: William Ingram
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL Programming standard :
CLL  The code is standard FORTRAN 77 except for having ! comments.
CLL  It conforms with standard A of version 3 (07/9/90) of UMDP 4, and
CLL  contains no 8X-deprecated features.
CLL
CLL Logical components covered : P232, D23
CLL  It is part of component P232 (longwave radiation)
CLL  It also performs some of the functions of D23
CLL
CLL Project task : P23    (radiation)
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
! longwave bands as described by Slingo and Wilderspin
! (April 1986, Quart.J.R.Met.Soc., 112, 472, 371-386) or UMDP 23.
C*L
      INTEGER!, INTENT (IN) ::
     &     L1,                      !  First dimension of input arrays
     &     L2,                      !  Number of points to be treated
     &     NLEVS                    !  Number of levels
      REAL!, INTENT(IN) ::
     &     TAC(L1,NLEVS),           !  Atmospheric temperatures at layer
     &     PEXNER(L1,NLEVS+1),      !  centres and Exner function at
     &                              !  layer boundaries (to get T there)
     &     PSTAR(L1),               !  Surface pressure
     &     AKH(NLEVS+1),            !  AK at layer boundaries
     &     BKH(NLEVS+1),            !  BK at layer boundaries
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

C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

      REAL PFCO (4,NBANDS)          !  Coefficients of the cubics
      DATA PFCO /
     &   -0.131802E-5,  0.141539E-2, -0.288741E-1, -7.05873,
     &   -0.250228E-5,  0.309318E-2, -0.556455,     28.9946,
     &   -0.212247E-5,  0.527366E-2,  -1.33625,     94.6759,
     &    0.719088E-5, -0.139883E-2, -0.100461,     24.2084,
     &    0.662392E-5, -0.118008E-2,  -.107828,     22.2919,
     &    0.488197E-4, -0.282063E-1,   5.54521,    -368.734 /
C
C
      INTEGER BAND, LEVEL, J        !  Loopers over band, level & point
      REAL TAB,                     !  Layer boundary temperature
     &     WTL, WTU,                !  and weights for its calculation
     &     BBC,                     !  Black-body flux at a model layer
     &     BBB,                     !             centre or boundary
     &     BBTFS                    !             or at TFS

      REAL
     &    PU                        ! Pressure at upper bdry of layer+1
     &,   PL                        ! Pressure at upper bdry of layer
     &,   PLM1                      ! Pressure at lower bdry of layer
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


C
CL    !  Section 1
CL    !  ~~~~~~~~~
CL    !  First find Planck function for surface temperature:
C
      DO 1 BAND=1, NBANDS
       BBTFS = ( ( PFCO(1,BAND) * TFS + PFCO(2,BAND) )
     &                * TFS + PFCO(3,BAND) ) * TFS + PFCO(4,BAND)
Cfpp$  Select(CONCUR)
       DO 11 J=1, L2
        IF ( AICE(J) .GT. 0. ) THEN
           TAB = ( TSTAR(J) + (AICE(J)-1.) * TFS ) / AICE(J)
         ELSE
           TAB = TSTAR(J)
        ENDIF
        BBB = ( ( PFCO(1,BAND) * TAB + PFCO(2,BAND) )
     &                * TAB + PFCO(3,BAND) ) * TAB + PFCO(4,BAND)
        SEAFX(J) = SEAFX(J) + BBTFS - BBB
        IF ( AICE(J) .GT. 0. )
     &              BBB = AICE(J) * BBB + (1.-AICE(J)) * BBTFS
        BHDB(J,1,BAND) = BBB
        IF ( SFUPON ) SFUP(J) = SFUP(J) + BBB
   11  CONTINUE
    1 CONTINUE
C
CL    !  Section 2
CL    !  ~~~~~~~~~
CL    !  Loop over the model layers finding the dB across each half
CL    !     of each one :
C
      DO 2 LEVEL=1, NLEVS-1
C      !
C      !  First must find TAB for the top of the layer (same for each
C      !                     band) and store it in spare space in THDB:
Cfpp$  Select(CONCUR)
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
C
        TAB =
     &   ( WTL * TAC(J,LEVEL) + WTU * TAC(J,LEVEL+1) ) / ( WTL + WTU )
        THDB(J,NLEVS,1) = TAB
   20  CONTINUE
C      !
C      !  This loop finds the black-body fluxes at the middle & top of
C      !  the layer and sets BHDB & THDB for the layer, using the
C      !  black-body flux for the bottom of the layer which is already
C      !  in BHDB and leaving the value at the top in the next level up
C      !                                                     of BHDB :
       DO 21 BAND=1, NBANDS
Cfpp$  Select(CONCUR)
        DO 22 J=1, L2
         BBC = ( ( PFCO(1,BAND) * TAC(J,LEVEL) + PFCO(2,BAND) )
     &    * TAC(J,LEVEL) + PFCO(3,BAND) ) * TAC(J,LEVEL) + PFCO(4,BAND)
         BHDB(J,LEVEL,BAND) = BBC - BHDB(J,LEVEL,BAND)
         TAB = THDB(J,NLEVS,1)
         BBB = ( ( PFCO(1,BAND) * TAB + PFCO(2,BAND) )
     &                * TAB + PFCO(3,BAND) ) * TAB + PFCO(4,BAND)
         THDB(J,LEVEL,BAND) = BBB - BBC
         BHDB(J,LEVEL+1,BAND) = BBB
   22   CONTINUE
   21  CONTINUE
    2 CONTINUE
C
CL    !  Section 3
CL    !  ~~~~~~~~~
CL    !  Finally, the top layer is treated slightly differently because
CL    !  the black-body flux at the top is zero:
C
      DO 3 BAND=1, NBANDS
Cfpp$  Select(CONCUR)
       DO 30 J=1, L2
        BBC = ( ( PFCO(1,BAND) * TAC(J,NLEVS) + PFCO(2,BAND) )
     &    * TAC(J,NLEVS) + PFCO(3,BAND) ) * TAC(J,NLEVS) + PFCO(4,BAND)
        BHDB(J,NLEVS,BAND) = BBC - BHDB(J,NLEVS,BAND)
        THDB(J,NLEVS,BAND) = - BBC
   30  CONTINUE
    3 CONTINUE
C
      RETURN
      END
