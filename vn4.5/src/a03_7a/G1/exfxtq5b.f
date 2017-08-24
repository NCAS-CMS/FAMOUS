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
!!!
!!!  Purpose: Calculate explicit fluxes of TL and QT
!!!
!!!  Suitable for single column use
!!!
!!!  Model           Modification history
!!! version  Date
!!!  4.3    23/2/97   New deck
!!!
!!!  JJ  <- Programmers of some or all of previous code or changes
!!!
!!!
!!!  Programming standard: UM Documentation Paper No 4, Version 2,
!!!                        dated 18/1/90
!!!
!!!  System component covered: P244
!!!
!!!  Project task: P24
!!!
!!!  Documentation: UM Documentation Paper No 24.
!!!
!!!---------------------------------------------------------------------

! SUBROUTINE EX_FLUX_TQ

      SUBROUTINE EX_FLUX_TQ (
     &  P_POINTS
     &, P_FIELD
     &, P1
     &, BL_LEVELS
     &, TL
     &, QW
     &, RDZ
     &, FTL
     &, FQW
     &, RHOKH
     &, LTIMER
     &  )


      IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

      LOGICAL
     &  LTIMER          ! IN Flag for TIMER diagnostics

      INTEGER
     & P_FIELD                ! IN No. of P-grid points in whole field
     &,P1                     ! IN First P-grid point to be processed
     &,P_POINTS               ! IN No. of P-grid points to be processed
     &,BL_LEVELS              ! IN No. of atmospheric levels for which
!                                boundary layer fluxes are calculated.
!                                Assumed ! <=30 for dimensioning GAMMA()
!                                in common deck C_GAMMA

      REAL
     &  TL(P_FIELD, BL_LEVELS)   ! IN Liquid/frozen water temperature
!                                     (K).
     &, QW(P_FIELD, BL_LEVELS)   ! IN Total water content (kg/kg air)
     &, RHOKH(P_FIELD, BL_LEVELS)! IN Exchange coeffs for moisture.
     &, RDZ(P_FIELD, BL_LEVELS)  ! IN RDZ(,1) is the reciprocal of the
!                                   height of level 1, i.e. of the
!                                   middle of layer 1.  For K > 1,
!                                   RDZ(,K) is the reciprocal of the
!                                   vertical distance from level
!                                   K-1 to level K.



! ARGUMENTS WITH INTENT OUT. IE: OUTPUT VARIABLES.

      REAL
     &  FTL(P_FIELD, BL_LEVELS)  ! OUT FTL(,K) contains net turbulent
!                                   sensible heat flux into layer K
!                                   from below; so FTL(,1) is the
!                                   surface sensible heat, H. (W/m2)
     &, FQW(P_FIELD, BL_LEVELS)  ! OUT Moisture flux between layers
!                                   (kg per square metre per sec).
!                                   FQW(,1) is total water flux
!                                   from surface, 'E'.

! LOCAL VARIABLES.

      INTEGER
     &  L, K

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
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


      REAL
     & GRCP

      PARAMETER (
     & GRCP = G/CP
     & )

      EXTERNAL TIMER

!-----------------------------------------------------------------------

      IF (LTIMER) THEN
        CALL TIMER('EXFLUXTQ',3)
      ENDIF

      DO K=2,BL_LEVELS
!-----------------------------------------------------------------------
!! 1. "Explicit" fluxes of TL and QW, on P-grid.
!-----------------------------------------------------------------------
        DO L=P1,P1+P_POINTS-1
        FTL(L,K) = -RHOKH(L,K) *
     &      ( ( ( TL(L,K) - TL(L,K-1) ) * RDZ(L,K) ) + GRCP )
        FQW(L,K) = -RHOKH(L,K) * ( QW(L,K) - QW(L,K-1) ) * RDZ(L,K)
        ENDDO
      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('EXFLUXTQ',4)
      ENDIF

      RETURN
      END
