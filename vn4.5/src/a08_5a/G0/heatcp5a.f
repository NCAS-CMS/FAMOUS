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
!    SUBROUTINE HEAT_CAP------------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE HEAT_CAP (NPNTS,SOIL_PTS,SOIL_INDEX,B,DZ,HCAP,
     &                     SATHH,SMCL,STHF,TSOIL,V_SAT,HCAPS
! LOGICAL LTIMER
     +,LTIMER
     +)

      IMPLICIT NONE
!
! Description:
!     Calculates the effective heat capacity of a soil layer
!     including the effects of ice, water and phase changes.
!                                                     (Cox, 6/95)
!
! Documentation : UM Documentation Paper 25
!
! Current Code Owner : David Gregory
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.1               New deck.    Peter Cox
!  4.5   26/5/98     Correction to stop overflow when SMCL very small.
!                    C.Bunton    
!
! Adds the vector shared, private directives.
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!

! Global variables:
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

! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER
     & NPNTS                ! IN Number of gridpoints.
     &,SOIL_PTS             ! IN Number of soil points.

!   Array arguments with intent(IN) :
      INTEGER
     & SOIL_INDEX(NPNTS)    ! IN Array of soil points.

      REAL
     & B(NPNTS)             ! IN Clapp-Hornberger exponent.
     &,DZ                   ! IN Thickness of the soil layer (m).
     &,HCAP(NPNTS)          ! IN Soil heat capacity (J/K/m3).
     &,SATHH(NPNTS)         ! IN Saturated soil water pressure (m).
     &,SMCL(NPNTS)          ! IN Soil moisture content of the
!                           !    layer (kg/m2).
     &,STHF(NPNTS)          ! IN Frozen soil moisture content of
!                           !    the layer as a fraction of
!                           !    saturation.
     &,TSOIL(NPNTS)         ! IN Layer temperature (K).
     &,V_SAT(NPNTS)         ! IN Volumetric soil moisture
!                           !    concentration at saturation
!                           !    (m3 H2O/m3 soil).
!
      LOGICAL LTIMER        ! Logical switch for TIMER diags

!   Array arguments with intent(OUT) :
      REAL
     & HCAPS(NPNTS)         ! OUT The total volumetric heat capacity
!                           !     (soil+water) of the layer (J/m3/K).

! Local scalars:
      INTEGER
     & I,J                  ! WORK Loop counter.

! Local arrays:
      REAL
     & DTHU(NPNTS)          ! WORK Rate of change of volumetric unfrozen
!                           !      soil moisture concentration with
!                           !      temperature (m3 liquid H2O/m3 soil/K)
     &,SMCLU(NPNTS)         ! WORK Unfrozen moisture content of each
!                           !      soil layer (kg/m2).
     &,SMCLSAT(NPNTS)       ! WORK The saturation moisture content of
!                           !      each layer (kg/m2).
     &,SMCLF(NPNTS)         ! WORK Frozen moisture content of each
!                           !      soil layer (kg/m2).
     &,TMAX(NPNTS)          ! WORK Temperature above which all water is
!                           !      unfrozen (Celsius)
     &,TSL(NPNTS)           ! WORK Soil layer temperature (Celsius)

C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
      REAL HCAPW,HCAPI,RHO_ICE,DPSIDT
      PARAMETER(
     + HCAPW=4180      ! Specific heat capacity of water (J/kg/K)
     +,HCAPI=2100      ! Specific heat capacity of ice (J/kg/K)
     +,RHO_ICE=917     ! Density of ice (kg/m3)
     +,DPSIDT=114.3    ! Rate of change of ice potential with temperatur
C                      ! RHO_ICE*LF/ZERODEGC*1/(RHO_WATER*G) (m/K)
     +)
C-----------------------------------------------------------------------
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------


      IF (LTIMER) THEN
        CALL TIMER('HEATCAP ',3)
      ENDIF
!----------------------------------------------------------------------
! Initialise all points
!----------------------------------------------------------------------
      DO I=1,NPNTS
        IF (V_SAT(I).GT.0.0) THEN ! Soil points
          HCAPS(I)=HCAP(I)
        ELSE ! Ice points
          HCAPS(I)=SNOW_HCAP
        ENDIF
      ENDDO

!-----------------------------------------------------------------------
! Calculate the temperature in Celsius and diagnose the saturation,
! frozen and unfrozen soil moisture contents of the layer
!-----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        TSL(I)=TSOIL(I)-ZERODEGC
        SMCLSAT(I)=RHO_WATER*DZ*V_SAT(I)
        SMCLF(I)=SMCLSAT(I)*STHF(I)
        SMCLU(I)=SMCL(I)-SMCLF(I)
      ENDDO


CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
!-----------------------------------------------------------------------
! Calculate TMAX, the temperature above which all soil water is
! unfrozen. Check that (SMCLSAT/SMCL)**B will not overflow when SMCL is
! very small. The function EPSILON  gives the number of type (real) of 
! SMCL that is negligeable compared to 1.                              
!-----------------------------------------------------------------------
          IF (SMCL(I) .GT. SMCLSAT(I)*(EPSILON(SMCL(I)) ) )THEN  
 
          TMAX(I)=-SATHH(I)/DPSIDT
     &           *(SMCLSAT(I)/SMCL(I))**(B(I))
        ELSE
          TMAX(I)=-ZERODEGC
        ENDIF

!-----------------------------------------------------------------------
! Calculate the rate of change of volumetric unfrozen soil moisture
! concentration with temperature
!-----------------------------------------------------------------------
        IF (TSL(I).GE.TMAX(I)) THEN
          DTHU(I)=0.0
        ELSE
          DTHU(I)=DPSIDT*SMCLSAT(I)
     &             /(B(I)*SATHH(I)*RHO_WATER*DZ)
     &             *(-DPSIDT*TSL(I)/SATHH(I))**(-1.0/B(I)-1.0)
        ENDIF

!-----------------------------------------------------------------------
! Calculate the effective heat capacity
!-----------------------------------------------------------------------
        HCAPS(I)=HCAP(I)+(HCAPW-HCAPI)*SMCLU(I)/DZ
     &          +HCAPI*SMCL(I)/DZ
     &          +RHO_WATER*DTHU(I)*((HCAPW-HCAPI)*TSL(I)+LF)
      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('HEATCAP ',4)
      ENDIF

      RETURN
      END
