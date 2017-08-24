C *****************************COPYRIGHT******************************
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
!    SUBROUTINE HEAT_CON----------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE HEAT_CON(NPNTS,HCON,STHU,STHF,
     &                    V_SAT,HCONS
C LOGICAL LTIMER
     +,LTIMER
     +)

      IMPLICIT NONE
!
! Description:
!    Calculates the soil thermal conductivity including the
!    effects of water and ice using the method suggested by
!    Farouki (1981) (note error in Verseghy, 1991)
!                                                          (Cox, 6/95)
!
! Documentation : UM Documentation Paper 25
!
! Current Code Owner : David Gregory
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.1               New deck.    Peter Cox
!LL   4.5   18/06/98  Changed Timer calls to indicate non-barrier
!LL                                                   P.Burton
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!


! Subroutine arguments:
!   Scalar arguments with intent(IN) :
      INTEGER
     & NPNTS              ! IN Number of gridpoints

      REAL
     & HCON(NPNTS)        ! IN Dry soil thermal conductivity (W/m/K).
     &,STHU(NPNTS)        ! IN Fractional saturation of unfrozen water
!                         !    at layer boundaries.
     &,STHF(NPNTS)        ! IN Fractional saturation of frozen water
!                         !    at layer boundaries.
     &,V_SAT(NPNTS)       ! IN Volumetric soil moisture concentration
!                         !    at saturation (m3/m3 soil).
C
      LOGICAL LTIMER      ! Logical switch for TIMER diags

!   Array arguments with intent(OUT) :
      REAL
     & HCONS(NPNTS)       ! OUT The thermal conductivity between adjacen
!                         !     layers including effects of water and ic
!                         !     (W/m/K).
!
C RHO_WATER removed to avoid clash with declaration in C_DENSTY
C J.Smith 28/06/95
      REAL OMEGA1,RHO_SNOW,DEFF_SNOW,SNOW_HCON,SNOW_HCAP
      INTEGER PSOIL
      PARAMETER (
     + PSOIL=4                  ! No. of soil layers (must = NSOIL).
     +,OMEGA1=3.55088E-4        ! Tunable characteristic freq (rad/s).
     +,RHO_SNOW=250.0           ! Density of lying snow (kg per m**3).
     +,DEFF_SNOW=0.1            ! Depth of `effective' snow surface
C                               ! layer (m).
     +,SNOW_HCON=0.265          ! Thermal conductivity of lying snow
C                               ! (Watts per m per K).
     +,SNOW_HCAP=0.63E6         ! Thermal capacity of lying snow
C                               ! (J/K/m3)
     +)
! Local scalars:
      INTEGER
     & I,J                ! WORK Loop counter.

! Local arrays:
      REAL
     & HCSAT(NPNTS)       ! WORK The thermal conductivity of the
!                         !  saturated  soil at current ratio of ice to
!                         !      liquid water (W/m/K).
     &,STH(NPNTS)         ! WORK Fractional saturation of water
!                         !     (liquid+ice) at layer boundaries.
     &,THICE(NPNTS)       ! WORK The concentration of ice at saturation
!                         !      for the current mass fraction of liquid
!                         !      water (m3 H2O/m3 soil).
     &,THWAT(NPNTS)       ! WORK The concentration of liquid water at
!                         !      saturation for the current mass
!                         !  fraction of liquid water  (m3 H2O/m3 soil).
!-----------------------------------------------------------------------
! Local Parameters (Source: "The Frozen Earth" p.90)
!-----------------------------------------------------------------------
      REAL
     & HCAIR              ! Thermal conductivity of air (W/m/K).
     &,HCICE              ! Thermal conductivity of ice (W/m/K).
     &,HCWAT              ! Thermal conductivity of liquid water (W/m/K)
      PARAMETER (HCAIR=0.025,HCICE=2.24,HCWAT=0.56)

      IF (LTIMER) THEN
        CALL TIMER('HEATCON ',103)
      ENDIF

!----------------------------------------------------------------------
! Initialise all points
!----------------------------------------------------------------------
      DO I=1,NPNTS
        IF (V_SAT(I).GT.0.0) THEN ! Soil points
          HCONS(I)=HCON(I)
        ELSE ! Ice points
          HCONS(I)=SNOW_HCON
        ENDIF
      ENDDO

      DO I=1,NPNTS
!---------------------------------------------------------------
! Only do calculation for non land-ice pts
! V_SAT is set to zero for land-ice points
!---------------------------------------------------------------
        IF (V_SAT(I).GT.0.0) THEN

          IF (STHU(I).GT.0.0) THEN
            THWAT(I)=V_SAT(I)*STHU(I)/(STHU(I)+STHF(I))
          ELSE
            THWAT(I)=0.0
          ENDIF

          IF (STHF(I).GT.0.0) THEN
            THICE(I)=V_SAT(I)*STHF(I)/(STHU(I)+STHF(I))
          ELSE
            THICE(I)=0.0
          ENDIF

          STH(I)=STHU(I)+STHF(I)
          HCSAT(I)=HCON(I)*(HCWAT**THWAT(I))*(HCICE**THICE(I))
     &                   /(HCAIR**V_SAT(I))
          HCONS(I)=(HCSAT(I)-HCON(I))*STH(I)+HCON(I)
        ENDIF

      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('HEATCON ',104)
      ENDIF
      RETURN
      END
