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
CLL
CLL  4.4    17/9/97   Updates snow grain size if required for 
CLL                   prognostic snow albedo               R. Essery
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
     + ASOIL,CONV_SNOW,LS_SNOW,SNOW_SUB,TSTAR,TIMESTEP,POINTS,
     + LYING_SNOW,RGRAIN,L_SNOW_ALBEDO,TS1,SNOWMELT,
     + SNOMLT_SUB_HTF,STF_HF_SNOW_MELT,
     + LTIMER)
      IMPLICIT NONE
      INTEGER POINTS       ! IN Number of points to be processed.
      REAL
     + ASOIL(POINTS)       ! IN Reciprocal areal heat capacity of top
C                          !    soil layer (K m2 / J).
     +,CONV_SNOW(POINTS)   ! IN Convective snowfall (kg per sq m per s)
     +,LS_SNOW(POINTS)     ! IN Large-scale snowfall (kg per sq m per s)
     +,SNOW_SUB(POINTS)    ! IN Sublimation of lying snow (kg/sq m/s).
     +,TIMESTEP            ! IN Timestep in seconds.
     +,TSTAR(POINTS)       ! IN Surface temperature (K).   
      REAL
     + LYING_SNOW(POINTS)  ! INOUT Snow on the ground (kg per sq m).
     +,RGRAIN(POINTS)      ! INOUT Snow grain size (microns).    
     +,SNOWMELT(POINTS)    ! IN    Surface snowmelt (kg/m2/s).
C                          ! OUT   Total snowmelt (kg/m2/s).
     +,TS1(POINTS)         ! INOUT Surface layer temperature (K).
      REAL
     + SNOMLT_SUB_HTF(POINTS)! OUT Sub-surface snowmelt heat flux (W/m2)
      LOGICAL
     + STF_HF_SNOW_MELT    ! IN Stash flag for snow melt heat flux
     +,L_SNOW_ALBEDO       ! IN Flag for prognostic snow albedo   
      LOGICAL
     + LTIMER
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

C NO EXTERNAL SUBROUTINES CALLED----------------------------------------
C*----------------------------------------------------------------------
C  Define local variables-----------------------------------------------
      REAL
     + SNOMLT_SUB          ! Sub-surface snow melt.
     +,R0                  ! Grain size for fresh snow (microns).
     +,RMAX                ! Maximum snow grain size (microns).
     +,RATE                ! Grain area growth rate (microns**2 / s).
     +,SNOWFALL            ! Snowfall in timestep (kg/m2).
      PARAMETER (R0 = 50., RMAX = 2000.)    
      INTEGER I            ! Loop counter; horizontal field index.
C
      IF (LTIMER) THEN
        CALL TIMER('SFSNOW  ',3)
      ENDIF
      DO 1 I=1,POINTS
        SNOMLT_SUB = 0.0
C
C-----------------------------------------------------------------------
CL Alter snowdepth as a result of snowfall, turbulent mass transport and
CL surface melt.
C-----------------------------------------------------------------------
C
        LYING_SNOW(I) = TIMESTEP*(LS_SNOW(I) + CONV_SNOW(I)) +
     &      MAX( 0.0, LYING_SNOW(I)-TIMESTEP*(SNOW_SUB(I)+SNOWMELT(I)) )
C
C-----------------------------------------------------------------------
CL Melt snow over land if TS1 is above freezing.
CL Adjust TS1 accordingly.
C-----------------------------------------------------------------------
C
        IF (TS1(I).GT.TM .AND. LYING_SNOW(I).GT.0.0) THEN
          IF (ASOIL(I).GT.0.0) THEN
            SNOMLT_SUB = MIN( LYING_SNOW(I)/TIMESTEP,
     &                            (TS1(I) - TM)/(LF*ASOIL(I)*TIMESTEP) )
!-----------------------------------------------------------------------
! For N/S boundaries in LAMS and polar rows in global model
! ASOIL will not have been calculated for these rows so diagnostics
! are not valid so just set to zero
!-----------------------------------------------------------------------
          ELSE
            SNOMLT_SUB = 0.0
          ENDIF
          TS1(I) = TS1(I) - ASOIL(I)*TIMESTEP*LF*SNOMLT_SUB
          LYING_SNOW(I) = LYING_SNOW(I) - TIMESTEP*SNOMLT_SUB
          SNOWMELT(I) = SNOWMELT(I) + SNOMLT_SUB
        ENDIF
        IF (STF_HF_SNOW_MELT) SNOMLT_SUB_HTF(I) = LF*SNOMLT_SUB
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
      
      IF (LTIMER) THEN
        CALL TIMER('SFSNOW  ',4)
      ENDIF
      RETURN
      END
