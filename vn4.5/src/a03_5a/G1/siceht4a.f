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
C*LL  SUBROUTINE SICE_HTF-----------------------------------------------
CLL
CLL  Purpose: Calculates heat flux through sea-ice (+ve downwards).
CLL           Sea-ice leads heat flux calculated in P243 (SF_EXCH).
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL   3.4   06/06/94  DEF TIMER replaced by LOGICAL LTIMER
CLL                   Argument LTIMER added
CLL                                                 S.J.Swarbrick
CLL
CLL  Note: At present the formulation is so simple as to make this
CLL        routine fairly trivial; but in future the formulation may
CLL        be revised so as to make a separate routine more obviously
CLL        worthwhile.
CLL
CLL  Programming standard: Unified Model Documentation Paper No.4
CLL                        version no.2, dated 18/1/90.
CLL
CLL  System component covered: P241
CLL
CLL  Documentation: ??
CLL
CLL *********************************************
CLL Penman-Monteith model. RE 19/1/95
CLL *********************************************
CLL Updates surface layer temperature and diagnoses surface temperature
CLL for sea-ice.
CLL
C*
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE SICE_HTF(ASHTF,DI,ICE_FRACTION,SURF_HT_FLUX,TIMESTEP,
     +                    LAND_MASK,POINTS,TI,TSTAR,ASURF,SEA_ICE_HTF,
     +                    LTIMER)
      IMPLICIT NONE
      LOGICAL LTIMER
      INTEGER POINTS        ! IN No of gridpoints to be processed.
      REAL
     + ASHTF(POINTS)        ! IN Coefficient to calculate surface
C                           !    heat flux into sea-ice (W/m2/K).
     +,DI(POINTS)           ! IN "Equivalent thickness" of sea-ice (m).
     +,ICE_FRACTION(POINTS) ! IN Fraction of gridbox covered by sea-ice.
     +,SURF_HT_FLUX(POINTS) ! IN Net downward heat flux at surface
C                           !    (W/m2).
     +,TIMESTEP             ! IN Timestep (s).
      LOGICAL
     + LAND_MASK(POINTS)    ! IN Land mask (T for land, F for sea).
      REAL
     + TI(POINTS)           ! INOUT  Sea-ice surface layer temperature
C                           !        (K). Set to TSTAR for unfrozen sea,
C                           !        missing data for land.
     +,TSTAR(POINTS)        ! INOUT Gridbox mean surface temperature (K)
     +,ASURF(POINTS)        ! OUT Reciprocal areal heat capacity of
C                           !     sea-ice surface layer (Km2/J).
     +,SEA_ICE_HTF(POINTS)  ! OUT Heat flux through sea-ice (W per sq m,
C                           !     positive downwards).
C-----------------------------------------------------------------------
CL  No workspace or EXTERNAL routines required.
C-----------------------------------------------------------------------
      EXTERNAL TIMER
C*
C  Common and local physical constants.
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

      REAL KAPPAI
      PARAMETER (
     + KAPPAI=2.09          ! Thermal conductivity of sea-ice (W per
C                           ! m per K).
     +)
      REAL DE
      PARAMETER (
     + DE = 0.1             ! Effective thickness of sea-ice surface
C                           ! layer (m).
     +)
C-----------------------------------------------------------------------
C*L-----------COMDECK C_SICEHC FOR SUBROUTINE IMPL_CAL----------
C AI  = reciprocal effective areal heat capacity of sea-ice,
C          ( 1 / (J per sq m per K)).
      REAL AI

      PARAMETER(AI  = 4.8E-6)
C*----------------------------------------------------------------------
C  Define local scalar.
      INTEGER I             ! Loop Counter; horizontal field index.
C-----------------------------------------------------------------------
CL  No significant structure.
C-----------------------------------------------------------------------
      IF (LTIMER) THEN
      CALL TIMER('SICEHTF ',3)
      ENDIF
      DO 1 I=1,POINTS
        IF (LAND_MASK(I)) THEN
          SEA_ICE_HTF(I)=0.0
          TI(I) = 0.0
        ELSE IF (ICE_FRACTION(I).LE.0.0) THEN
          SEA_ICE_HTF(I)=0.0
          TI(I) = TSTAR(I)
        ELSE
          ASURF(I) = AI / ICE_FRACTION(I)
          SEA_ICE_HTF(I) = KAPPAI*ICE_FRACTION(I)*(TI(I) - TFS)/DI(I)
          TSTAR(I) = (1. - ICE_FRACTION(I))*TFS + ICE_FRACTION(I)*TI(I)
     &                  + SURF_HT_FLUX(I)/ASHTF(I)
          TI(I) = TI(I) + ASURF(I)*TIMESTEP*
     &                    (SURF_HT_FLUX(I) - SEA_ICE_HTF(I))
        ENDIF
    1 CONTINUE
      IF (LTIMER) THEN
      CALL TIMER('SICEHTF ',4)
      ENDIF
      RETURN
      END
