C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!!!  SUBROUTINE SICE_HTF-----------------------------------------------
!!!
!!!  Purpose: Calculates heat flux through sea-ice (+ve downwards).
!!!           Sea-ice leads heat flux calculated in P243 (SF_EXCH).
!!!
!!!  Model            Modification history
!!! version  date
!!!
!!!  Note: At present the formulation is so simple as to make this
!!!        routine fairly trivial; but in future the formulation may
!!!        be revised so as to make a separate routine more obviously
!!!        worthwhile.
!!!
!!!  Programming standard: Unified Model Documentation Paper No.4
!!!                        version no.2, dated 18/1/90.
!!!
!!!  System component covered: P241
!!!
!!!  Documentation: ??
!!!
!!! *********************************************
!!! Penman-Monteith model. RE 19/1/95
!!! *********************************************
!!! Updates surface layer temperature and diagnoses surface temperature
!!! for sea-ice.


! Arguments:---------------------------------------------------------
      SUBROUTINE SICE_HTF (
     & ASHTF,DI,ICE_FRACTION,SURF_HT_FLUX,TIMESTEP
     &,LAND_MASK,P_FIELD,POINTS,P1,TI,TSTAR,ASURF,SEA_ICE_HTF
     &,LTIMER)
      IMPLICIT NONE

      LOGICAL LTIMER

      INTEGER
     & POINTS               ! IN No of gridpoints to be processed.
     &,P_FIELD              ! IN Total Number of points on p-grid
     &,P1                   ! IN First point of p grid to be processed

      REAL
     & ASHTF(P_FIELD)       ! IN Coefficient to calculate surface
!                                heat flux into sea-ice (W/m2/K).
     &,DI(P_FIELD)          ! IN "Equivalent thickness" of sea-ice (m).
     &,ICE_FRACTION(P_FIELD)! IN Fraction of gridbox covered by sea-ice.
     &,SURF_HT_FLUX(P_FIELD)! IN Net downward heat flux at surface W/m2

     &,TIMESTEP             ! IN Timestep (s).

      LOGICAL
     & LAND_MASK(P_FIELD)   ! IN Land mask (T for land, F for sea).

      REAL
     & TI(P_FIELD)          ! INOUT  Sea-ice surface layer temperature
!                              (K). Set to TSTAR for unfrozen sea,
!                               missing data for land.
     &,TSTAR(P_FIELD)       ! INOUT Gridbox mean surface temperature (K)
     &,ASURF(P_FIELD)       ! OUT Reciprocal areal heat capacity of
!                              sea-ice surface layer (Km2/J).
     &,SEA_ICE_HTF(P_FIELD) ! OUT Heat flux through sea-ice (W per sq m,
!                              positive downwards).

!-----------------------------------------------------------------------
!!  No workspace or EXTERNAL routines required.
!-----------------------------------------------------------------------

      EXTERNAL TIMER

!  Common and local physical constants.
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

!  Define local scalar.
      INTEGER I             ! Loop Counter; horizontal field index.
!-----------------------------------------------------------------------
!!  No significant structure.
!-----------------------------------------------------------------------

      IF (LTIMER) THEN
        CALL TIMER('SICEHTF ',3)
      ENDIF


      DO I=P1,P1+POINTS-1
        IF (LAND_MASK(I)) THEN
          SEA_ICE_HTF(I)=0.0
          TI(I) = 1.0E30
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
      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('SICEHTF ',4)
      ENDIF

      RETURN
      END
