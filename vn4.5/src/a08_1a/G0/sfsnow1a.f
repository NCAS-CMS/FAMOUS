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
CLL  SUBROUTINE SFSNOW ------------------------------------------------
CLL
CLL  Purpose:  Calculates the decrease/increase in snowdepth due to the
CLL            sublimation/deposition of lying snow; adds the large-
CLL            scale and convective snowfall to the snowdepth;
CLL            melts snow when the input surface temperature is above
CLL            the melting point of ice, and adjusts the surface
CLL            temperature to account for latent cooling thus caused.
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   4.2    Oct. 96  T3E migration: *DEF CRAY removed (was used to
CLL                    switch on dynamic allocation) 
CLL                                    S.J.Swarbrick
CLL   4.4    17/9/97  Updates snow grain size if required for 
CLL                   prognostic snow albedo        R. Essery
CLL   4.5    01/10/98 Removed old section-version defs. K Rogers
CLL
CLL  Programming standard: Unified Model Documentation Paper No.4
CLL                        version no. 2, dated 18/1/90.
CLL
CLL  Logical component covered: P251.
CLL
CLL  System task:
CLL
CLL  Documentation: um documentation paper no 25
CLLEND------------------------------------------------------------------
C
C*L  ARGUMENTS:---------------------------------------------------------
      SUBROUTINE SFSNOW(
     + CONV_SNOW,HCAP,HCON,LS_SNOW,SNOW_SUB,TIMESTEP,POINTS,
     + RGRAIN,L_SNOW_ALBEDO,LYING_SNOW,TSTAR,SNOWMELT
     +,SNOWMELT_HTF,STF_SNOWMELT_HTF
     +)
      IMPLICIT NONE
      INTEGER POINTS       ! IN Number of points to be processed.
      REAL
     + CONV_SNOW(POINTS)   ! IN Convective snowfall (kg per sq m per s)
     +,HCAP(POINTS)        ! IN Soil heat capacity (J per K per m**3).
     +,HCON(POINTS)        ! IN Soil thermal conductivity
C                          !    (J per m per s per K).
     +,LS_SNOW(POINTS)     ! IN Large-scale snowfall (kg per sq m per s)
     +,SNOW_SUB(POINTS)    ! IN Sublimation of lying snow (kg/sq m/s).
      REAL TIMESTEP        ! IN Timestep in seconds.
      REAL
     + LYING_SNOW(POINTS)   ! INOUT Snow on the ground (kg per sq m).
     +,RGRAIN(POINTS)       ! INOUT Snow grain size (microns).
     +,TSTAR(POINTS)        ! INOUT Surface temperature (K).
      REAL
     + SNOWMELT(POINTS)     ! OUT Snowmelt (kg per sq m per second).
     +,SNOWMELT_HTF(POINTS) ! OUT Snowmelt heat flux (Watts per sq m).
      LOGICAL STF_SNOWMELT_HTF ! IN stash flag for snow melt heat flux
     +,L_SNOW_ALBEDO        ! IN Flag for prognostic snow albedo   
C-----------------------------------------------------------------------
C*
C Define common then local physical constants --------------------------
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
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

      REAL OMEGA1
      PARAMETER (
     + OMEGA1=3.55088E-4   ! Characteristic freq (rad per sec; tunable).
     +)
C
C WORKSPACE USAGE-------------------------------------------------------
C  1 real work array of full length is required.
      REAL
     + ASOIL(POINTS)       ! WORK ASOIL(1) - may be replaced by i/p.
C NO EXTERNAL SUBROUTINES CALLED----------------------------------------
C*----------------------------------------------------------------------
C  Define local variables-----------------------------------------------
      REAL
     + SMELT_TO_DT         ! See comments for contents.
     +,R0                  ! Grain size for fresh snow (microns).
     +,RMAX                ! Maximum snow grain size (microns).
     +,RATE                ! Grain area growth rate (microns**2 / s).
     +,SNOWFALL            ! Snowfall in timestep (kg/m2).
      PARAMETER (R0 = 50., RMAX = 2000.)
      INTEGER I            ! Loop counter; horizontal field index.
C
      DO 1 I=1,POINTS
        ASOIL(I)=0.0
        SMELT_TO_DT=1.0
C
C-----------------------------------------------------------------------
CL 1. Alter snowdepth as a result of snowfall and turbulent mass
CL    transport (equations P251.1 and P251.2 combined).
C-----------------------------------------------------------------------
C
        LYING_SNOW(I) = MAX ( 0.0, LYING_SNOW(I)-TIMESTEP*SNOW_SUB(I) )
        LYING_SNOW(I) = LYING_SNOW(I) + TIMESTEP *
     &                                   ( LS_SNOW(I) + CONV_SNOW(I) )
C
C-----------------------------------------------------------------------
CL 2. Melt snow over land if TSTAR is above freezing.
CL    Snowmelt is calculated as kg per sq m PER TIMESTEP, for internal
CL    convenience/efficiency.  Bear this in mind when comparing code
CL    with equations in the documentation.
C-----------------------------------------------------------------------
CL 2.1 Calculate snowmelt according to equation P251.3.
C-----------------------------------------------------------------------
C
        IF(TSTAR(I).GT.TM .AND. LYING_SNOW(I).GT.0.0)THEN
C
C  See P242 documentation for definition of ASOIL, i.e. ASOIL(1), the
C  reciprocal areal heat capacity of the top soil layer (sq m K per J).
C
          ASOIL(I)=
     +             TIMESTEP*SQRT(OMEGA1/(2.0*HCON(I)*HCAP(I)))
C
C  SMELT_TO_DT is in Kelvin per (kg per sq m of snowmelt).
C  It is the factor to convert between snowmelt and temperature change.
C
          SMELT_TO_DT=ASOIL(I)*LF/TIMESTEP
          SNOWMELT(I)=MIN(LYING_SNOW(I),(TSTAR(I)-TM)/SMELT_TO_DT)
C
C-----------------------------------------------------------------------
CL 2.2 Adjust snowdepth and TSTAR.  See eqs. P251.4, P251.5.
C-----------------------------------------------------------------------
C
          LYING_SNOW(I)=LYING_SNOW(I)-SNOWMELT(I)
          TSTAR(I)=TSTAR(I)-SNOWMELT(I)*SMELT_TO_DT
C
C-----------------------------------------------------------------------
CL 2.3 Convert snowmelt to kg of water per sq metre second.
C      Calculate heat flux required to melt snow (defined to be GE 0).
C-----------------------------------------------------------------------
C
          SNOWMELT(I)=SNOWMELT(I)/TIMESTEP
          IF (STF_SNOWMELT_HTF) SNOWMELT_HTF(I)=LF*SNOWMELT(I)
C
C-----------------------------------------------------------------------
CL 3. Set snowmelt to zero if no snow, or below freezing, or not land.
C-----------------------------------------------------------------------
C
        ELSE  ! (If no snow, or below freezing.)
          SNOWMELT(I)=0.0
          IF (STF_SNOWMELT_HTF) SNOWMELT_HTF(I)=0.0
        ENDIF
    1 CONTINUE

!-----------------------------------------------------------------------
! Increment snow grain size used in albedo calculations
!-----------------------------------------------------------------------
      IF ( L_SNOW_ALBEDO ) THEN
        DO I=1,POINTS 
          IF ( LYING_SNOW(I) .GT. 0.) THEN
            SNOWFALL = TIMESTEP*(LS_SNOW(I) + CONV_SNOW(I))
            RATE = 0.6
            IF (TSTAR(I) .LT. TM) THEN
              IF (RGRAIN(I) .LT. 150.) THEN
                RATE = 0.06
              ELSE
                RATE = 0.23E6*EXP(-3.7E4/(8.13451*TSTAR(I)))
              ENDIF
            ENDIF
            RGRAIN(I) = SQRT( RGRAIN(I)**2 + (RATE/3.14159)*TIMESTEP )
     &                                   - (RGRAIN(I) - R0)*SNOWFALL/2.5
            RGRAIN(I) = MIN( RMAX, RGRAIN(I) )
            RGRAIN(I) = MAX( R0, RGRAIN(I) )
          ELSE
            RGRAIN(I) = R0
          ENDIF
        ENDDO
      ENDIF

      RETURN
      END
