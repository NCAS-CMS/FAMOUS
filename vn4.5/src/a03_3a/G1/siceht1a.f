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
CLL   4.1   08/05/96  decks A03_2C and A03_3B removed
CLL                                     S D Jackson
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
C*
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE SICE_HTF (
     + DI,ICE_FRACTION,LAND_MASK,TSTAR,POINTS,SEA_ICE_HTF,LTIMER
     +)
      IMPLICIT NONE
      LOGICAL LTIMER
      INTEGER POINTS        ! IN No of gridpoints to be processed.
      REAL
     + DI(POINTS)           ! IN "Equivalent thickness" of sea-ice (m).
     +,ICE_FRACTION(POINTS) ! IN Fraction of gridbox covered by sea-ice.
     +,TSTAR(POINTS)        ! IN Gridbox mean surface temperature (K).
      LOGICAL
     + LAND_MASK(POINTS)    ! IN Land mask (T for land, F for sea).
      REAL
     + SEA_ICE_HTF(POINTS)  ! OUT Heat flux through sea-ice (W per sq m,
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
C  Define local scalar.
      INTEGER I             ! Loop Counter; horizontal field index.
C-----------------------------------------------------------------------
CL  No significant structure.
C-----------------------------------------------------------------------
      IF (LTIMER) THEN
      CALL TIMER('SICEHTF ',3)
      ENDIF
      DO 1 I=1,POINTS
        IF (.NOT.LAND_MASK(I) .AND. ICE_FRACTION(I).GT.0.0) THEN
          SEA_ICE_HTF(I)=KAPPAI*(TSTAR(I)-TFS)/DI(I)         ! Eq P241.3
        ELSE
          SEA_ICE_HTF(I)=0.0
        ENDIF
    1 CONTINUE
      IF (LTIMER) THEN
      CALL TIMER('SICEHTF ',4)
      ENDIF
      RETURN
      END
