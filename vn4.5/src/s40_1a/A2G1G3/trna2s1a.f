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
CLL    SUBROUTINE TRANSA2S
CLL    -------------------
CLL
CLL   THIS ROUTINE TRANSFERS DATA NEEDED FOR
CLL   COUPLING FROM THE ATMOSPHERE TO THE SLAB OCEAN.
CLL   IT CAN BE COMPILED BY CFT77, BUT DOES NOT CONFORM TO THE
CLL   ANSI FORTRAN77 STANDARDS BECAUSE OF ENDDOs AND !COMMENTS
CLL   CALLED BY: SLABCNTL
CLL   VERSION NUMBER 1.1
CLL   WRITTEN BY D L ROBERTS (14/1/91)
CLL   MODIFIED BY A.B.KEEN (02/02/93)
CLL   MODIFIED BY C.A.SENIOR (22/03/93)
CLL   MODIFIED BY C.A.SENIOR (08/07/93)
CLL   MODIFIED BY C.A.SENIOR (25/02/94)
CLL   REVIEWED BY W.INGRAM (01/03/93)
CLL   FOLLOWS DOCUMENTATION PAPER 3, VERSION 5 FOR STANDARDS.
CLL   DOCUMENTATION: UM DOCUMENTATION PAPER 58; THE SLAB OCEAN MODEL
CLLEND
C*L
C-----------------------------------------------------------------
C
      SUBROUTINE TRANSA2S(L1,L2,
     +                   L_THERM,
     +                   LAND,
     +                   TSTARATM,
     +                   SLABTEMP,
     +                   HICEATM,
     +                   HICESLB,
     +                   HICEMIN,
     +                   HSNOWATM,
     +                   HSNOWSLB,
     +                   AICEATM,
     +                   AICESLB,
     +                   AICEMIN,
     +                   TCLIM,
     +                   TCLIMC,
     +                   HCLIM)
C
C
C     THE FLOW OF CONTROL IS STRAIGHTFORWARD.
C
      INTEGER
     + L1              ! IN SIZE OF DATA VECTORS
     +,L2              ! IN AMOUNT OF DATA TO BE PROCESSED
C
      LOGICAL LAND(L1)              ! IN ATMOSPHERE MODEL LAND-SEA
     +                              ! MASK (FALSE AT OCEAN POINTS)
     +,L_THERM                      ! IN TRUE FOR COUPLED MODEL TYPE
     +                              !    ICE THERMODYNAMICS
C
      REAL
     + TSTARATM(L1)    ! IN SURFACE TEMPERATURE OF ATMOSPHERE (K)
     +,SLABTEMP(L1)    ! INOUT SLAB OCEAN TEMP (C)
     +,HICEATM(L1)     ! IN EQUIVALENT ICE DEPTH FROM ATMOSPHERE
     +                 !    THIS IS THE DEPTH OF ICE THAT HAS THE
     +                 !    SAME THERMAL CONDUCTIVITY AS THE SEA-ICE
     +                 !    TOGETHER WITH THE SNOW
     +,HICESLB(L1)     ! OUT ICE DEPTH FOR USE IN SLAB (M)
     +,HSNOWATM(L1)    ! IN SNOW DEPTH FROM ATMOSPHERE (KG/M**2)
     +,HSNOWSLB(L1)    ! OUT SNOW DEPTH FOR USE IN SLAB (M)
     +,AICEATM(L1)     ! IN ICE CONCENTRATION FROM ATMOSPHERE
     +,AICESLB(L1)     ! OUT ICE CONCENTRATION FOR USE IN SLAB
     +,TCLIM(L1)       ! IN CLIMATOLOGICAL SST K
     +,TCLIMC(L1)      ! OUT CLIMATOLOGICAL SST C
     +,HCLIM(L1)       ! IN CLIMATOLOGICAL SEA-ICE DEPTH M
C
      REAL
     + AICEMIN         ! IN MINIMUM CONCENTRATION OF ICE IF ICE PRESENT
     +,HICEMIN         ! IN MINIMUM DEPTH OF ICE IF ICE PRESENT
     +                 ! PREVENTS SMALL ICE DEPTHS CAUSING FAILURE
C
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
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------
C
C     LOCAL VARIABLES
C
      INTEGER
     + J                   ! LOOP COUNTER
C
      REAL TFREEZE      ! FREEZING POINT OF SEA WATER IN C
     +,ONEEM8           ! SMALL +VE VALUE TO ELIMINATE ROUNDING
C
      PARAMETER(TFREEZE = TFS - ZERODEGC)
      PARAMETER(ONEEM8  = 1.0E-8 )
C
C
C     -----------------------------------------------------------
C
      DO J = 1,L2
C
       IF  (.NOT. LAND(J)) THEN
C
C         1. CONVERT SLAB OCEAN TEMPERATURE AND CLIMATOLOGICAL SST
C            FROM KELVIN TO CELSIUS AND SET
C            VALUES OF SLABICE FOR ICE POINTS AND LAND POINTS.
C
C         2. CONVERT FROM EQUIVALENT ICE DEPTH
C            TO MEAN ACTUAL ICE DEPTH
C            NOTE THAT THE EXTRA PIECE OF ICE, HICEMIN, ADDED
C            IN AT THE END OF THE LAST SLAB TIMESTEP TO PREVENT
C            VERY SMALL ICE DEPTHS CAUSING THE ATMOSPHERE MODEL
C            TO FAIL IS NOW REMOVED
C
             IF (L_THERM) THEN
               TCLIMC(J) = TCLIM(J) - ZERODEGC
             ELSE
               IF ( HCLIM(J) .LT. ONEEM8 )  THEN
                 TCLIMC(J) = TCLIM(J) - ZERODEGC
               ELSE
                 TCLIMC(J) = TFREEZE
               ENDIF
             ENDIF
C
             IF (l_therm) THEN
               IF ( AICEATM(J) .LT. AICEMIN )   THEN
                  HICESLB(J)  = 0.0
               ELSE
                  HICESLB(J)  = AICEATM(J) * ( HICEATM(J) -
     +                          CONRATIO * ( HSNOWATM(J) / RHOSNOW ) )
     +                          - HICEMIN
               ENDIF
             ELSE
               IF ( AICEATM(J) .LT. AICEMIN )   THEN
                  SLABTEMP(J) = TSTARATM(J) - ZERODEGC
                  HICESLB(J)  = 0.0
               ELSE
                  SLABTEMP(J) = TFREEZE
                  HICESLB(J)  = AICEATM(J) * ( HICEATM(J) -
     +                          CONRATIO * ( HSNOWATM(J) / RHOSNOW ) )
     +                          - HICEMIN
               ENDIF
             ENDIF
C
C         3. CONVERT SNOWDEPTH FROM KG/M**2 TO M
C
             HSNOWSLB(J) = HSNOWATM(J) / RHOSNOW
C
C         4. INITIALISE THE ICE CONCENTRATIONS TO BE USED IN THE
C            SLAB MODEL. NB THE 'OLD VALES' INPUT FROM THE ATMOSPHERE
C            MUST BE SAVED AS THEY ARE USED BY THE ROUTINE
C            TRANSS2A, HENCE THE TWO ARRAYS!
C
             AICESLB(J) = AICEATM(J)
C
       ELSE
C
             SLABTEMP(J)  = RMDI
             TCLIMC(J)    = RMDI
C
             HSNOWSLB(J)  = RMDI
C
             HICESLB(J)   = RMDI
C
             AICESLB(J)   = RMDI
C
       ENDIF
      END DO
C
      RETURN
      END
