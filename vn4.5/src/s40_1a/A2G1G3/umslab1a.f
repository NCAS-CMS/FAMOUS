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
CLL   SUBROUTINE UMSLAB
CLL   -----------------
CLL
CLL   THIS ROUTINE IS FOR USE WITH THE 'SLAB' OCEAN MODEL ONLY.
CLL
CLL   THE TEMPERATURE OF THE MIXED LAYER OCEAN IS UPDATED
CLL   ACCORDING TO THE EQUATION
CLL
CLL   CW * RHOW * DTM      NET HEAT FLUX IN  + HEAT CONVERGENCE
CLL               ---  =  ------------------------------------
CLL               DT                         H
CLL
CLL   THIS ROUTINE FORMS PART OF SYSTEM COMPONENT P40.
CLL   IT CAN BE COMPILED BY CFT77, BUT DOES NOT CONFORM TO ANSI
CLL   FORTRAN77 STANDARDS, BECAUSE OF THE INLINE COMMENTS.
CLL   IT ADHERES TO THE STANDARDS OF DOCUMENTATION PAPER 3, VERSION 5.
CLL
CLL   ALL QUANTITIES IN THIS ROUTINE ARE IN S.I. UNITS UNLESS
CLL   OTHERWISE STATED.
CLL
CLL   DIFFUSION OF SLAB TEMPERATURES ADDED
CLL
CLL   CALLED BY: SLABCNTL
CLL
CLL   WRITTEN BY C.A.SENIOR (14/1/91)
CLL   MODIFIED BY A.B.KEEN (02/02/93)
CLL   MODIFIED BY C.A.SENIOR(22/03/93)
CLL   MODIFIED BY C.A.SENIOR(14/09/93)
CLL   MODIFIED BY J.F.THOMSON (23/05/94) Thermodynamics altered to
CLL                                      allow addition of dynamics.
CLL   Version   Description of change
CLL     4.0     Vertical SST advection included (R.Carnell) and
CLL             change order of heat convergence calculation for
CLL             ice dynamics. J.F.Crossley
!LL  4.4   04/08/97  Add missing ARGOINDX to various argument lists.
!LL                  D. Robinson.
!LL
CLL   VERSION NUMBER 1.1
CLL   REVIEWER: W.INGRAM (01/03/93)
CLL
CLLEND---------------------------------------------------------------
C*L
      SUBROUTINE UMSLAB(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     +                 ICY,
     +                 ATMSFLUX,
     +                 ADJHCONV,
     +                 SNOWRATE,
     +                 CARYHEAT,
     +                 SLABTEMP,
     +                 wtsfc,
     +                 wtbase,
     +                 LAND,
     +                 AICE,
     +                 NEWICE,
     +                 OIFLUX,
     +                 DTADV,DTDIFF,
     +                 L1,L2,DT,DZ1,u_field,
     +                 JROWS,JROWSM1,ICOLS,
     +                 TCLIMC,
     +                 HCLIM,
     +                 CALIB,
     +                 L_THERM,L_IDYN,L_IDRIF,
     +                 L_SLBADV,LGLOBAL,
     +                 EDDYDIFF,
     +                 DELTA_LONG,DELTA_LAT,BASE_LAT,
     +                 AH,
     +                 OPENSEA,
     +                 UCURRENT,
     +                 VCURRENT,
     +                 COS_P_LATITUDE,COS_U_LATITUDE,
     +                 SEC_P_LATITUDE,SIN_U_LATITUDE)
C
      IMPLICIT NONE
C
!========================== COMDECK TYPOINDX ===========================
!
! Description:
!
!       This comdeck contains all the indices and row-wise loop
!       control variables required by the ocean MPP code.
!
! History:
!
!=======================================================================

      INTEGER J_1     ! Local value of loop control for J = 1, n
     &,       J_2     !   "     "    "   "     "     "  J = 2, n
     &,       J_3     !   "     "    "   "     "     "  J = 3, n
     &,       J_JMT   !   "     "    "   "     "     "  J = n, JMT
     &,       J_JMTM1 !   "     "    "   "     "     "  J = n, JMTM1
     &,       J_JMTM2 !   "     "    "   "     "     "  J = n, JMTM2
     &,       J_JMTP1 !   "     "    "   "     "     "  J = n, JMTP1
     &,       JST     ! First row this process considers (no halo)
     &,       JFIN    ! Last   "    "     "     "        "
     &,       J_FROM_LOC
     &,       J_TO_LOC
     &,       JMT_GLOBAL       ! Global value of JMT
     &,       JMTM1_GLOBAL     ! Global value of JMT - 1
     &,       JMTM2_GLOBAL     ! Global value of JMT - 2
     &,       JMTP1_GLOBAL     ! Global value of JMT + 1
     &,       J_OFFSET         ! Global value of JST - 1
     &,       O_MYPE           ! MYPE for ocean arg lists
     &,       O_EW_HALO        ! EW HALO for ocean arg lists
     &,       O_NS_HALO        ! NS HALO for ocean arg lists
     &,       J_PE_JSTM1
     &,       J_PE_JSTM2
     &,       J_PE_JFINP1
     &,       J_PE_JFINP2
     &,       O_NPROC
     &,       imout(4),jmout(4)! i,j indices for pts in Med outflow
     &,       J_PE_IND_MED(4)  ! no for each PE in Med outflow
     &,       NMEDLEV          ! no of levels for Med outflow
     &,      lev_med   ! level at which deep advective Med outflow
c                        exits the Mediterranean
     &,      lev_hud   ! level at which deep advective flow
c                        enters the Hudson Bay
     &,      imout_hud(4)  ! zonal index for Hudson Bay outflow
     &,      jmout_hud(4)  ! merid index for Hudson Bay outflow
     &,      J_PE_IND_HUD(4)  ! PE's involved in HB outflow
     &,      med_topflow   ! last level for which there is inflow to 
C                          ! Mediterranean




      INTEGER L1,       ! IN SIZE OF DATA VECTORS
     + L2               ! IN AMOUNT OF DATA TO BE PROCESSED
     +,JROWS             ! IN NO OF ROWS N-S
     +,JROWSM1           ! IN NO OF ROWS N-S - 1
     +,ICOLS             ! IN NO OF COLUMNS E-W
     +,u_field           ! In points in u field
C
      REAL ATMSFLUX(L1) ! IN NET DOWNWARD SURFACE HEAT FLUX (W M-2)
     +,ADJHCONV(L1)     ! IN ADJUSTED HEAT CONVERGENCE IN W M-2.
     +,SNOWRATE(L1)     ! IN RATE OF SNOWFALL, IN KG M-2 S-1.
     +,TCLIMC(L1)       ! IN CLIMATOLOGICAL SEA SURFACE TEMPS C
     +,HCLIM(L1)        ! IN CLIMATOLOGICAL SEA-ICE EXTENTS M
C
      REAL
     + COS_P_LATITUDE(L1)! IN COS LATITUDE ON P GRID
     +,COS_U_LATITUDE(L1)! IN COS LATITUDE ON UV GRID
     +,SEC_P_LATITUDE(L1)! 1.0/COS LATITUDE ON P GRID
     +,SIN_U_LATITUDE(L1)! IN SIN LATITUDE ON UV GRID
C
      REAL DT           ! IN TIMESTEP FOR UPDATING THE SLAB OCEAN IN S.
     +,DZ1              ! IN THICKNESS OF THE SLAB OCEAN IN METRES.
     +,AH               ! DIFFUSION COEFFICENT
     +,EDDYDIFF         ! IN EDDY DIFFUSION COEFF FOR OIFLUX CALCULATION
     +,UCURRENT(ICOLS,JROWSM1) ! IN ZONAL SURFACE CURRENT (M/S)
     +,VCURRENT(ICOLS,JROWSM1) ! IN MERIDIONAL SURFACE CURRENT (M/S)
C
      REAL
     & DELTA_LONG             ! IN EW grid spacing (degrees)
     &,DELTA_LAT              ! IN NS grid spacing (degrees)
     &,BASE_LAT               ! IN latitude of first row (degrees)
C
      LOGICAL ICY(L1)   ! IN TRUE IF BOX CONTAINS SEA-ICE.
     +,CALIB            ! IN TRUE IF CALIBRATION EXPT
     +,OPENSEA(L1)      ! IN TRUE IF BOX CONTAINS OPEN SEA (NO ICE)
     +,LAND(L1)         ! IN TRUE AT LAND POINTS
     +,NEWICE(L1)       ! OUT TRUE IF ICE IS FORMING
     +,L_THERM          ! IN TRUE FOR COUPLED LIKE ICE THERMODYNAMICS
     +,L_IDYN           ! IN TRUE FOR CAV FLUID ICE DYNAMICS
     +,L_IDRIF          ! IN TRUE FOR ICE DEPTH ADVECTION
     +,L_SLBADV         ! IN TRUE FOR SLAB TEMP ADVECTION
     +,LGLOBAL         ! IN TRUE IF GLOBAL MODEL
C
      REAL SLABTEMP(L1)  ! INOUT TEMPERATURE OF THE SLAB OCEAN.
     +,wtsfc(l1)         ! OUT W x sfc slab temp
     +,wtbase(l1)        ! OUT w x base slab temp
     +,AICE(L1)          ! INOUT ICE FRACTION
     +,OIFLUX(L1)        ! OUT OCEAN TO ICE HEAT FLUX
     +,CARYHEAT(L1)      ! OUT ZERO EXCEPT AT POINTS WHERE ICE IS
     +                   !     ABOUT TO FORM, WHERE IT SHOULD BE THE
     +                   !     NEGATIVE HEAT FLUX
     +,DTADV(L1)         ! OUT HEATING RATE DUE TO ADVECTION K/S
     +,DTDIFF(L1)        ! OUT HEATING RATE DUE TO DIFFUSION K/S
C*
C     Include COMDECKS
C
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

C
C
C     VARIABLES LOCAL TO THIS ROUTINE ARE NOW DEFINED.
C
      INTEGER J           ! LOOP COUNTERS
C
      REAL TFREEZE          ! FREEZING POINT OF SEAWATER IN C
     +,SLABD(L1)            ! TEMPORARY STORE FOR SLABTEMP-TCLIMC
     +,TWORK(L1)            ! TEMPORARY STORE FOR SLABTEMP
     +,DTTOFLUX             ! CONVERSION FACTOR FROM K TO WM-2
     +,FLUXTODT             ! CONVERSION FACTOR FROM WM-2 TO K
     +,DEPTHTOK             ! CONVERSION FACTOR FROM M ICE TO K
     +,AHDT                 ! DIFFUSION COEFFICENT * TIMESTEP
C
      PARAMETER(TFREEZE=TFS-ZERODEGC)
C
C
C     SET VARIOUS CONSTANTS FOR USE IN CALCULATING THE HEAT
C     CONVERGENCE
C
      DTTOFLUX = RHOCP * DZ1 / DT
      FLUXTODT = 1.0 / DTTOFLUX
      DEPTHTOK = RHOICE * LF / ( RHOCP * DZ1 )
C
      DO J = 1,L2
C
C     1. SET TEMPERATURE OF ICE POINTS TO TFREEZE (if coupled model
C        ice thermodynamics is not required)
C        AND INITIALISE CARYHEAT AND ADJHCONV (IN CALIBRATION)
C        and initialise newice for coupled model ice thermo.
C
        IF ( ICY(J) .and. .not. l_therm) THEN
           SLABTEMP(J) = TFREEZE
        END IF
C
        IF ( CALIB ) THEN
           ADJHCONV(J) = 0.0
        END IF
C
        CARYHEAT(J) = 0.0
        newice(j)   = .false.
C
C     2. CALCULATE THE CHANGE IN SLAB TEMPERATURE
C        and ocean to ice heat flux
C
C     PERFORM OPERATIONS OVER SEA POINTS ONLY
C  (for coupled ice thermodynamics, sea points means sea and sea ice)
C
C First coupled thermodynamics
        IF (L_THERM) THEN
          IF ( .not. LAND(J) ) THEN
            oiflux(j)   = rhocp * eddydiff / (0.5*dz1)
     &                  * ( slabtemp(j) -tfreeze )
            if (opensea(j)) oiflux(j) = 0.0
            SLABTEMP(J) = SLABTEMP(J) + ( ATMSFLUX(J)
     &                    - LF*SNOWRATE(J)*(1.0-aice(j))
     &                    - oiflux(j) ) * fluxtodt
C
          ENDIF
C Then original thermodynamic treatment
         ELSE
          IF ( OPENSEA(J) ) THEN
            SLABTEMP(J) = SLABTEMP(J) + ( ATMSFLUX(J)
     +                    - LF*SNOWRATE(J) ) * FLUXTODT
          ENDIF
        ENDIF
C End J loop.
      END DO
C
C     3. CALCULATE ADVECTION OF SLAB OCEAN TEMPERATURE
C (and diagnose heating rate)
      IF (L_SLBADV) THEN
      do j=1,l1
        twork(j) = slabtemp(j)
      end do
        CALL SLAB_T_UV(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &  L1,L2,ICOLS,JROWS,JROWSM1,LAND
     & ,LGLOBAL,u_field,dz1
     & ,DELTA_LAT,DELTA_LONG,BASE_LAT,DT
     & ,COS_U_LATITUDE,COS_P_LATITUDE
     & ,SEC_P_LATITUDE,SIN_U_LATITUDE
     & ,UCURRENT,VCURRENT
     & ,SLABTEMP
     & ,OPENSEA
     & ,wtsfc
     & ,wtbase
     &         )
      do j=1,l1
        dtadv(j) = ( slabtemp(j)-twork(j) ) / dt
      end do
      ENDIF
C
C
C        USE TEMPORARY SPACE TO HOLD SLABTEMP-TCLIMC FOR USE
C        IN DIFFUSION SUBROUTINE
C
      DO J=1,L2
        IF ( .NOT. CALIB ) THEN
          SLABD(J) = SLABTEMP(J) - TCLIMC(J)
        ENDIF
      END DO
C
C
C     4. CALCULATE DEL2 DIFFUSION OF SLAB TEMPERATURES
C        (DIFFUSE SLABTEMP-TCLIMC - SO DONT CALL IN
C         CALIBRATION EXPERIMENT)
C
        IF ( .NOT. CALIB) THEN
         IF ( AH .GT. 0.0) THEN
          AHDT= AH * DT
          CALL SLABDIFF(SLABD,
     +                  OPENSEA,
     +                  L1,L2,
     +                  JROWS,ICOLS,
     +                  AHDT,
     +                  DELTA_LONG,DELTA_LAT,BASE_LAT,
     +                  COS_P_LATITUDE,COS_U_LATITUDE,SEC_P_LATITUDE)
         ENDIF
        ENDIF
C
C     RECALCULATE SLABTEMP IN NON-CALIBRATION EXPERIMENT
C     PUT DIFFUSION INCREMENT INTO NEW STORE TO ADD TO HEAT CONV
C
C     5. IF MIXED LAYER TEMPERATURE IS BELOW FREEZING THEN ICE FORMATION
C        IS JUST STARTING. CALCULATE THE HEAT DEFICIT THAT ARISES
C        FROM THE TEMPERATURE CHANGE AND PASS TO SLABICE THROUGH
C        CARYHEAT(J) FOR ICE FORMATION TO OCCUR. SET SLABTEMP TO
C        TFREEZE
C First coupled thermodynamics
C
      IF (l_therm) then
      DO J = 1,L2
        IF ( .NOT. CALIB) THEN
          DTDIFF(J)   = SLABD(J) - SLABTEMP(J) + TCLIMC(J)
          ADJHCONV(J) = ADJHCONV(J) + DTDIFF(J) * DTTOFLUX
          dtdiff(j) = dtdiff(j) / dt
        ENDIF
C
        IF ( .NOT.LAND(J) ) THEN
C
          IF ( SLABTEMP(J) .LT. TFREEZE ) THEN
             CARYHEAT(J) = SLABTEMP(J) - TFREEZE
             CARYHEAT(J) = (CARYHEAT(J) * RHOCP * DZ1) / DT
             SLABTEMP(J) = TFREEZE
             NEWICE(J) = (.NOT. ICY(J))
          ENDIF
          IF ( CALIB ) THEN
            IF ( l_idyn .or. l_idrif ) THEN
              ADJHCONV(J) = ( TCLIMC(J) -
     +                    SLABTEMP(J) ) * DTTOFLUX
            ELSE
              ADJHCONV(J) = ( ( TCLIMC(J) - DEPTHTOK * HCLIM(J) ) -
     +                    SLABTEMP(J) ) * DTTOFLUX
            END IF
          END IF
C
          SLABTEMP(J) = SLABTEMP (J) + ADJHCONV(J) * FLUXTODT
C
        ENDIF
C End J loop.
      END DO
C Then original thermodynamic treatment
      ELSE
      DO J = 1,L2
        IF ( .NOT. CALIB) THEN
          DTDIFF(J)   = SLABD(J) - SLABTEMP(J) + TCLIMC(J)
          ADJHCONV(J) = ADJHCONV(J) + DTDIFF(J) * DTTOFLUX
          dtdiff(j) = dtdiff(j) / dt
        ENDIF
C
        IF ( OPENSEA(J) ) THEN
C
          IF ( CALIB ) THEN
            ADJHCONV(J) = ( ( TCLIMC(J) - DEPTHTOK * HCLIM(J) ) -
     +                  SLABTEMP(J) ) * DTTOFLUX
          END IF
C
          SLABTEMP(J) = SLABTEMP (J) + ADJHCONV(J) * FLUXTODT
C
          IF ( SLABTEMP(J) .LT. TFREEZE ) THEN
             CARYHEAT(J) = SLABTEMP(J) - TFREEZE
             CARYHEAT(J) = (CARYHEAT(J) * RHOCP * DZ1) / DT
             SLABTEMP(J) = TFREEZE
             NEWICE(J) = (.NOT. ICY(J))
          ENDIF
        ENDIF
C End J loop.
      END DO
C End IF block for coupled thermodynamics.
      ENDIF
      RETURN
      END
