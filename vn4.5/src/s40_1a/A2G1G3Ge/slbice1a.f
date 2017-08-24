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
C*LL
CLL   SUBROUTINE SLABICE
CLL   ------------------
CLL
CLL   THIS ROUTINE IS FOR USE WITH THE 'SLAB' OCEAN MODEL ONLY.
CLL   FOR COUPLED OCEAN-ATMOSPHERE RUNS THE ROUTINE ICEFLOE SHOULD
CLL   BE USED. (OBVIOUSLY MUCH OF THE BASIC PHYSICS IS THE SAME.)
CLL
CLL   THERMODYNAMIC SEA-ICE MODEL, BASED ON THE ZERO-LAYER MODEL OF
CLL   SEMTNER, A.J. (1976) : J.PHYS.OCEANOGR., 6, 379-389,
CLL   MODIFIED ALONG THE LINES SUGGESTED BY
CLL   GORDON, C., AND BOTTOMLEY, M. (1984) : DCTN 1.
CLL   INCLUDES A REPRESENTATION OF LEADS (WITHOUT DYNAMICS) DUE TO
CLL   HIBLER, W.D. (1979) : J.PHYS.OCEANOGR., 9, 815-846.
CLL
CLL   THIS ROUTINE FORMS PART OF SYSTEM COMPONENT P40.
CLL   IT CAN BE COMPILED BY CFT77, BUT DOES NOT CONFORM TO ANSI
CLL   FORTRAN77 STANDARDS, BECAUSE OF THE INLINE COMMENTS
CLL   AND USE OF ENDDO.
CLL   IT ADHERES TO THE STANDARDS OF DOCUMENTATION PAPER 3, VERSION 5.
CLL
CLL   ALL QUANTITIES IN THIS ROUTINE ARE IN S.I. UNITS UNLESS
CLL   OTHERWISE STATED.
CLL
CLL   CALLED BY: SLABCNTL
CLL
CLL   WRITTEN BY D.L.ROBERTS (12/11/90)
CLL   MODIFIED BY A.B.KEEN (02/02/93)
CLL   MODIFIED BY C.A.SENIOR (22/03/93)
CLL   MODIFIED BY A.B.KEEN (05/08/93) (WHITE ICE CODE CORRECTED)
CLL   MODIFIED BY C.A.SENIOR (25/02/94)
CLL   MODIFIED BY J.F.THOMSON (26/05/94)
CLL   Version  Description of change
CLL     4.0    Increased small value used as minimum ice depth to
CLL            prevent rounding errors to match that used in the
CLL            coupled model and correct code for dynamic sea ice.
CLL            J.F.Crossley
CLL   VERSION NUMBER 1.1
CLL   REVIEWER: W.INGRAM (01/03/93)
CLL   DOCUMENTATION: UM DOCUMENTATION PAPER 58; THE SLAB OCEAN MODEL
CLL
CLLEND---------------------------------------------------------------
C*L
      SUBROUTINE SLABICE(ICY,
     +                   LDSFLUX,
     +                   NEWICE,
     +                   OIFLUX,
     +                   ATMSFLUX,
     +                   ADJHCONV,
     +                   HICE,
     +                   HSNOW,
     +                   SNOWRATE,
     +                   SUBLIMZ,
     +                   AICE,
     +                   CARYHEAT,
     +                   TOPMELTZ,
     +                   BOTMELTZ,
     +                   SLABTEMP,
     +                   L1,L2,DT,DZ1,H0,
     +                   L_THERM,L_IDYN,L_IDRIF,
     +                   AICEMAX,
     +                   AICEMIN,
     +                   TCLIMC,
     +                   HCLIM,
     +                   CALIB)
C
      INTEGER L1,        ! IN SIZE OF DATA VECTORS
     + L2                ! IN AMOUNT OF DATA TO BE PROCESSED
C
      REAL ATMSFLUX(L1)  ! IN NET HEAT FLUX OVER LEADS (W M-2)
     +                   !    (FLUX INTO SLAB IN NEW THERMODYNAMICS)
     +,ADJHCONV(L1)      ! IN ADJUSTED HEAT CONVERGENCE W M-2
     +                   ! INOUT IF CALIBRATION
     +,LDSFLUX(L1)       ! IN NET HEAT FLUX OVER LEADS INTO ICE
     +                   !    (NEW THERMODYNAMICS)
     +,OIFLUX(L1)        ! IN OCEAN TO ICE HEAT FLUX
     +                   !    (NEW THERMODYNAMICS)
     +,SNOWRATE(L1)      ! IN RATE OF SNOWFALL, IN KG M-2 S-1.
     +,SUBLIMZ(L1)       ! IN RATE OF SUBLIMATION, IN KG M-2 S-1.
     +,TOPMELTZ(L1)      ! IN RATE OF MELTING OF SNOW IN W M-2.
     +                   !    (THIS CAN BE TRANSFERRED TO ICE.)
     +,BOTMELTZ(L1)      ! IN DIFFUSIVE HEAT FLUX THROUGH ICE. IF
     +                   !    THIS IS +VE, ICE MELTS AT THE BASE.
     +                   !    IF IT IS -VE, ICE ACCRETES THERE.
     +,CARYHEAT(L1)      ! IN ZERO EXCEPT AT POINTS WHERE ICE IS
     +                   !   ABOUT TO FORM, WHERE IT SHOULD BE THE
     +                   !   NEGATIVE HEAT FLUX USED TO CREATE THE ICE.
     +,TCLIMC(L1)        ! IN CLIMATOLOGICAL SEA SURFACE TEMPS K
     +,HCLIM(L1)         ! IN CLIMATOLOGICAL SEA-ICE EXTENTS M
      REAL DT            ! IN TIMESTEP FOR UPDATING THE SLAB OCEAN.
     +,DZ1               ! IN THICKNESS OF THE SLAB OCEAN IN METRES.
C
      REAL
     + H0        ! IN MINIMUM LOCAL DEPTH OF NEWLY-FORMED ICE, IN M.
     +,AICEMAX   ! IN MAX ICE CONCENTRATION ALLOWED IN AN ICY GRID BOX.
     +,AICEMIN   ! IN MIN ICE CONCENTRATION ALLOWED IN AN ICY GRID BOX.
C
      LOGICAL ICY(L1)    ! INOUT TRUE IF BOX CONTAINS ICE.
     +                   ! RESET FOR NEW ICE POINTS + MELTED POINTS
     +,NEWICE(L1)        ! IN TRUE IF ICE IS FORMING.
     +,CALIB             ! IN TRUE IF CALIBRATION EXPT.
     +,L_THERM           ! IN TRUE FOR NEW THERMODYNAMICS.
     +,L_IDYN            ! IN TRUE FOR CAV FLUID ICE DYNAMICS.
     +,L_IDRIF           ! IN TRUE FOR ICE DEPTH ADVECTION.
C
      REAL AICE(L1)      ! INOUT ICE CONCENTRATION.
     +,HICE(L1)          ! INOUT MEAN ICE DEPTH OVER WHOLE GRID BOX.
     +,HSNOW(L1)         ! INOUT SNOW DEPTH, NOT AVERAGED OVER GRID
     +                   !       BOX, JUST OVER THE ICE PORTION.
     +,SLABTEMP(L1)      ! INOUT TEMPERATURE OF THE SLAB OCEAN.
C*
C     Include COMDECKS
C
C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------
C*L---------------COMDECK C_SLAB----------------------------------------
C PARAMETERS REQUIRED BY SLAB OCEAN MODEL AND NOT DEFINED ELSEWHERE
C
C CONRATIO IS THE RATIO OF THERMAL CONDUCTIVITIES OF ICE AND SNOW
C (DIMENSIONLESS)
C RHOCP IS THE VOLUMETRIC HEAT CAPACITY OF SEA WATER (J K-1 M-3)
C RHOICE IS THE DENSITY OF ICE (KG M-3)
C RHOSNOW IS THE DENSITY OF SNOW (KG M-3)
C NB ** RHOSNOW is also defined in the common deck C_SOILH, which
C cannot be used in the slab routines as it contains a duplicate
C definition of RHO_WATER, which is also found in C_DENSTY **
C ** It should be noted that the value of RHOSNOW defined here matches
C    the value defined in C_SOIL_H, but differs from that currently
C    used in the ocean GCM (300 Kg m-3)
C
       REAL CONRATIO,RHOCP,RHOICE,RHOSNOW
C
       PARAMETER(CONRATIO=6.5656)
       PARAMETER(RHOCP=4.04E6)
       PARAMETER(RHOICE=900.0)
       PARAMETER(RHOSNOW=250.0)
C
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

C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
C
C     VARIABLES LOCAL TO THIS ROUTINE ARE NOW DEFINED.
C
      INTEGER J           ! LOOP COUNTERS
      REAL XSENERGY       ! TEMPORARY SCALAR VARIABLE FOR HOLDING
     +                    ! AMOUNTS OF SURPLUS ENERGY.
     +,TEST               ! ANOTHER TEMPORARY SCALAR VARIABLE.
     +,DELA               ! CHANGE IN AICE OVER ONE TIME STEP.
     +,DELH               ! CHANGE IN HICE OVER ONE TIME STEP.
     +,XSSUBLIM           ! REMAINING SUBLIMATION AFTER ALL
     +                    ! SNOW HAS BEEN MELTED
     +,LEADFLUX           ! FLUX OVER LEADS TO CHANGE ICE DEPTH
     +,WORKC              ! TEMPORARY STORE
     +,AICEREF            ! ORIGINAL VALUES OF AICE
C
C
      INTEGER L2BY2    ! HALF THE NUMBER OF DATA POINTS
      INTEGER L2BY2P1  ! HALF THE NUMBER OF DATA POINTS, PLUS ONE
      REAL H0R         ! THE RECIPROCAL OF H0.
     +,TFREEZE   ! THE FREEZING POINT OF SEAWATER
     +,DTBYRHOS  ! THE RATIO OF THE TIMESTEP TO THE DENSITY OF SNOW.
     +,DENRATIO  ! THE RATIO OF THE DENSITY OF ICE TO THAT OF WATER.
     +,DENRAT2   ! THE RATIO OF THE DENSITY OF SNOW TO THAT OF WATER.
     +,DENRAT3   ! THE RATIO OF THE DENSITY OF SNOW TO THAT OF ICE.
     +,DENRATM1  ! DENRATIO MINUS ONE.
C
      REAL QI    ! THE VOLUMETRIC HEAT OF FUSION OF ICE, IN J/M**3.
     +,QS        ! THE VOLUMETRIC HEAT OF FUSION OF SNOW, IN J/M**3.
     +,QIR       ! THE RECIPROCAL OF QI.
     +,QSR       ! THE RECIPROCAL OF QS.
     +,QIRDT     ! THE RATIO OF THE TIMESTEP TO THE VOLUMETRIC HEAT
     +           !   OF FUSION OF ICE.
     +,QIBYDT    ! THE RATIO OF THE VOLUMETRIC HEAT OF FUSION OF ICE
     +           !   TO THE TIMESTEP.
     +,DTTOFLUX   ! CONVERSION FACTOR FOR K TO WM-2
     +,FLUXTODT   ! CONVERSION FACTOR FOR WM-2 TO K
     +,HEATCON    ! TEMPORARY STORE FOR HEAT CONVERGENCE
     +,ONEEM4     ! SMALL +VE VALUE TO ELIMINATE ROUNDING PROBLEMS
C
C     SET VARIOUS CONSTANTS,CONVERTING TO S.I. UNITS WHERE NECESSARY.
C
      PARAMETER ( TFREEZE  = TFS - ZERODEGC)
      PARAMETER ( DENRATIO = RHOICE / RHOSEA)
      PARAMETER ( DENRAT2  = RHOSNOW / RHOSEA)
      PARAMETER ( DENRAT3  = RHOSNOW / RHOICE)
      PARAMETER ( DENRATM1 = DENRATIO - 1.0)
      PARAMETER      ( QI  = LF * RHOICE)
      PARAMETER      ( QS  = LF * RHOSNOW)
      PARAMETER      ( QIR = 1.0 / QI)
      PARAMETER      ( QSR = 1.0 / QS)
      PARAMETER   ( oneem4 = 1.0E-4 )
C
      L2BY2    = L2 / 2
      L2BY2P1  = L2BY2 + 1
      H0R      = 1.0 / H0
      QIRDT    = QIR * DT
      QIBYDT   = QI / DT
      DTTOFLUX = RHOCP * DZ1 / DT
      FLUXTODT = 1.0 / DTTOFLUX
      DTBYRHOS = DT / RHOSNOW
C
C
      DO J = 1,L2     ! Loop over all points
C
C     TRAP AREAL FRACTIONS GREATER THAN THE SPECIFIED MAXIMUM
C     THIS IS DONE NOW RATHER THAN IN TRANSA2S AS NORTHERN AND
C     SOUTHERN HEMISPHERE MAX ARE REQUIRED SEPERATELY
C
        IF ( AICE(J) .GT. AICEMAX ) AICE(J) = AICEMAX
C
C     ZERO SOME VARIABLES
C
        DELH = 0.0
        DELA = 0.0
        WORKC = 0.0
C
C       NOTE THAT ATMSFLUX,TOPMELT,BOTMELT AND SUBLIMZ ARE GRID_BOX
C       MEAN VALUES AND MUST BE ADJUSTED APPROPRIATELY
C
C       ATMSFLUX IS THE TOTAL DOWNWARD HEAT FLUX FROM ATMOSPHERE OVER
C       THE LEADS AND IS DIVIDED BY THE LEAD FRACTION TO
C       CONCENTRATE THE FLUX OVER THE LEAD AREA
C
C       TOPMELT,BOTMELT AND SUBLIMATION ARE DIVIDED BY ICE FRACTION
C       TO CONCENTRATE FLUX OVER ICE AREA
C
        IF ( ICY(J) ) THEN
C
C     UPDATE THE SNOW DEPTH. (NOTE THAT ANY NEGATIVE SUBLIMATION,
C     I.E. FROST, GETS USED HERE TO INCREASE SNOW DEPTH. IT DOES NOT
C     PENETRATE TO THE ICE SECTION OF THE CODE.)
C
          HSNOW(J) = HSNOW(J) + DTBYRHOS *
     +               ( SNOWRATE(J) - ( SUBLIMZ(J) / AICE(J) ) )
          IF ( HSNOW(J) .GE. 0.0 ) THEN
            XSSUBLIM = 0.0
          ELSE
            XSSUBLIM = - ( RHOSNOW / RHOICE ) * HSNOW(J)
            HSNOW(J) = 0.0
          ENDIF
C
C     XSSUBLIM NOW CONTAINS ANY SUBLIMATION UNUSED DUE TO LACK OF SNOW.
C     CONVERTED TO M OF ICE
C
          HSNOW(J) = HSNOW(J) - DT * QSR * TOPMELTZ(J) / AICE(J)
C
C     IF SNOW GOES -VE,USE EXCESS HEAT TO MELT ICE INSTEAD.
C
          IF ( HSNOW(J) .LT. 0.0 ) THEN
            XSENERGY = - ( QS / QI ) * HSNOW(J)
            HSNOW(J) = 0.0
          ELSE
            XSENERGY = 0.0
          ENDIF
C
C     CALCULATE CHANGE IN ICE DEPTH OVER THE ICE FRACTION.
C
          DELH = - XSENERGY - XSSUBLIM
     +           - QIRDT * ( BOTMELTZ(J) / AICE(J) )
          IF (l_therm) delh = delh - qirdt*(oiflux(j)+caryheat(j))
C
C     FORM WEIGHTED MEAN OF DEPTH CHANGES OVER ICE AND OVER LEADS.
C
          LEADFLUX = DT*( LF * SNOWRATE(J) -
     +            ( ATMSFLUX(J) / (1.0 - AICE(J) ) ) )
          if (l_therm) leadflux = dt * ( lf * snowrate(j)
     &      - ( ldsflux(j) / (1.0 - aice(j) ) ) )
     &      - dt * ( oiflux(j) + caryheat(j) )
          DELH = QIR * LEADFLUX * ( 1.0 - AICE(J) )
     +           + DELH * AICE(J)
C
C     IF CALIBRATION EXPERIMENT, CALCULATE THE HEAT CONVERGENCE
C     NEEDED TO FORCE THE ICE DEPTH TO BE THE CLIMATOLOGICAL VALUE
C     (only done here if no ice dynamics)
          if (.not.(l_idyn.or.l_idrif)) then
            IF ( CALIB ) THEN
              ADJHCONV(J) = QIBYDT * ( HICE(J) + DELH - HCLIM(J) )
            END IF
C
C     SUBTRACT HEAT CONVERGENCE TERM, WHICH OPERATES OVER THE ENTIRE
C     GRID BOX.
C
C     IN CALIBRATION EXPERIMENT THIS GIVES DELH = HCLIM - HICE
C
            DELH = DELH - QIRDT * ADJHCONV(J)
          ENDIF
C
C     TEST WHETHER ALL ICE IS GOING TO MELT.
C
          IF ( (HICE(J) + DELH) .LE. oneem4 ) THEN
C
C     NOW SEE WHETHER THERE IS SUFFICIENT SNOW TO MAINTAIN ICE COVER
C     BY CONVERSION OF SNOW INTO ICE. IF SO, ALL THE SNOW IS CONVERTED.
C     IT COULD BE ARGUED THAT ONLY ENOUGH SNOW TO MAINTAIN ICE COVER
C     SHOULD BE CONVERTED. HOWEVER, THE WHITE ICE CODE WOULD
C     THEN CONVERT ALL BUT A THIN DUSTING OF SNOW TO ICE IN ANY CASE.
C
            TEST = HICE(J) + DELH
     +                     + DENRAT3 * AICE(J) * HSNOW(J)
C
            IF ( TEST .GT. oneem4 ) THEN
C
              DELH = DELH + DENRAT3 * AICE(J) * HSNOW(J)
              HSNOW(J) = 0.0
C
            ELSE
C
              ICY(J)      = .FALSE.
              HICE(J)     = 0.0
              AICE(J)     = 0.0
              HSNOW(J)    = 0.0
              SLABTEMP(J) = SLABTEMP(J) - (TEST*QI)/(RHOCP*DZ1)
C
C      IN CALIBRATION MODE CATER FOR THE CASE WHERE THE MODEL HAS
C      SIMULATED ICE BUT THE CLIMATOLOGY INDICATES OPEN WATER
C      THIS OCCURS IN PRACTICE WHERE THE CLIMATOLOGY HAS JUST
C      CHANGED FROM ICE TO WATER BUT THE MODEL KEEPS ICE
C      Only if no ice dynamics included.
C
              IF ( CALIB ) THEN
                HEATCON     = ( TCLIMC(J) - SLABTEMP(J) ) * DTTOFLUX
                ADJHCONV(J) = ADJHCONV(J) + HEATCON
                SLABTEMP(J) = SLABTEMP(J) + HEATCON * FLUXTODT
              END IF
            ENDIF
C
C     IF ICE WILL REMAIN,CALCULATE THE CHANGE IN THE AREAL FRACTION.
C
          ELSE  ! For points where HICE > 0.0
            IF (LEADFLUX .GT. 0.0) THEN
              DELA = ( 1.0 - AICE(J) )* H0R * QIR * LEADFLUX
            ENDIF
            IF (DELH .LT. 0.0) THEN
              DELA = DELA + ( AICE(J) * DELH ) / ( 2.0 * HICE(J) )
            ENDIF
C
C         END OF BLOCK CHECKING IF ALL ICE MELTS
C
          ENDIF
C
C       END OF IF(ICY) BLOCK
C
        ENDIF
C
C
C     SPECIAL CODE FOR GRID BOXES WHERE ICE FORMATION IS JUST STARTING.
C
        IF ( NEWICE(J) ) THEN
          DELH = -QIRDT * CARYHEAT(J)
          DELA = H0R * DELH
        ENDIF
C
C     END OF SPECIAL CODE FOR NEWLY-ICED GRID BOXES.
C
C     RESET ICY TO INCLUDE NEWLY-ICED POINTS. THEN UPDATE ICE DEPTHS
C     AND AREAL FRACTIONS. (NEWLY-MELTED POINTS WON'T BE UPDATED.)
C
C     IN CALIBRATION EXPERIMENT HICE = HICE + ( HCLIM - HICE)
C                                    = HCLIM
C
        ICY(J) = ICY(J) .OR. ( NEWICE(J) )
C
        IF ( ICY(J) ) THEN
C
C     STORE THE INITIAL AREAL FRACTIONS FOR USE IN UPDATING
C     SNOW DEPTHS.
C
        AICEREF = AICE(J)
C
          HICE(J) = HICE(J) + DELH
          AICE(J) = AICE(J) + DELA
C
C     TRAP AREAL FRACTIONS GREATER THAN THE SPECIFIED MAXIMUM
C     OR LESS THAN THE SPECIFIED MINIMUM.
C
          IF ( AICE(J) .GT. AICEMAX ) AICE(J)=AICEMAX
C
          IF ( AICE(J) .LT. AICEMIN ) AICE(J)=AICEMIN
C
C     REDISTRIBUTE SNOW OVER THE ALTERED ICE AREA,USING THE INITIAL
C     AREAL FRACTIONS HELD IN AICEREF.
C
          HSNOW(J) = ( HSNOW(J) * AICEREF ) / AICE(J)
C
C     THE NEXT SECTION OF CODE DEALS WITH THE FORMATION OF 'WHITE'
C     ICE, WHICH OCCURS WHEN SNOW ACCUMULATES TO SUCH AN EXTENT THAT
C     THE SNOW/ICE INTERFACE SINKS BELOW THE WATERLINE. IT BEGINS BY
C     FINDING THE HEIGHT OF THE WATERLINE ABOVE THE TOP OF THE ICE.
C
          WORKC = ( ( HICE(J) / AICE(J) ) * DENRATM1
     +                        + DENRAT2 * HSNOW(J) ) / DENRAT3
C
C     IF THIS HEIGHT IS POSITIVE, IT INDICATES THE DEPTH OF SNOW
C     THAT WILL BE CONVERTED TO ICE.
C
          IF ( WORKC .GT. 0.0 ) THEN
            HSNOW(J) = HSNOW(J) - WORKC
            WORKC    = WORKC * AICE(J)
            HICE(J)  = HICE(J) + DENRAT3 * WORKC
          ENDIF
        ENDIF
C
      END DO    ! End of loop over all points
C
      RETURN
      END
