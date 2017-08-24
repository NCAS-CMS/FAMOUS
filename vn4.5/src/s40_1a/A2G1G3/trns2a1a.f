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
CLL    SUBROUTINE TRANSS2A
CLL    -------------------
CLL
CLL    THIS ROUTINE TRANSFERS DATA NEEDED FOR COUPLING
CLL    FROM THE SLAB OCEAN TO THE ATMOSPHERE,PERFORMING VARIOUS
CLL    MANIPULATIONS ON THE WAY.
CLL    THE SNOW DEPTH FIELD OVER LAND POINTS CONTAINS LAND SNOW
CLL    DEPTHS USED BY THE ATMOSPHERE MODEL AND SO IS NOT RESET
CLL    IT CAN BE COMPILED BY CFT77, BUT DOES
CLL    NOT CONFORM TO THE ANSI FORTRAN77 STANDARDS
CLL    BECAUSE OF IN_LINE COMMENTS AND THE USE OF ENDDO
CLL    CALLED BY: SLABCNTL
CLL    VERSION NUMBER 1.1
CLL    WRITTEN BY D L ROBERTS (14/1/91)
CLL    MODIFIED BY A.B.KEEN (02/02/93)
CLL    MODIFIED BY C.A.SENIOR (25/02/94)
CLL    REVIEWED BY W.INGRAM (01/03/93)
CLL    FOLLOWS DOCUMENTATION PAPER 3, VERSION 5 FOR STANDARDS.
CLL
CLLEND
C*L
C-----------------------------------------------------------------
      SUBROUTINE TRANSS2A(L1,L2,
     +                    LAND,
     +                    SLABTEMP,
     +                    TSTARATM,
     +                    AICESLB,
     +                    AICEATM,
     +                    HICESLB,
     +                    HICEATM,
     +                    HICEMIN,
     +                    HSNOWSLB,
     +                    HSNOWATM,
     +                    AICEMIN)
C
C     THE FLOW OF CONTROL IS STRAIGHTFORWARD.
C
      INTEGER
     + L1,              ! IN SIZE OF DATA VECTORS
     + L2               ! IN AMOUNT OF DATA TO BE PROCESSED
C
      LOGICAL LAND(L1) ! IN ATMOSPHERIC MODEL LAND-SEA
     +                  !       MASK (FALSE AT OCEAN POINTS).
C
      REAL
     + SLABTEMP(L1)   ! IN TEMPERATURE OF OCEAN SURFACE LAYER
     +,TSTARATM(L1)   ! INOUT SURFACE TEMPERATURE OF ATMOSPHERIC MODEL
     +,AICESLB(L1)    ! IN ICE CONCENTRATION FROM SLAB
     +,AICEATM(L1)    ! INOUT ICE CONCENTRATION IN ATMOS MODEL
     +,HICESLB(L1)    ! IN ICE DEPTH FROM SLAB
     +,HICEATM(L1)    ! OUT ICE DEPTH IN ATMOSPHERIC MODEL
     +,HSNOWSLB(L1)   ! IN SNOW DEPTH FROM SLAB
     +,HSNOWATM(L1)   ! INOUT SNOW DEPTH IN ATMOSPHERIC MODEL
C*
      REAL
     + AICEMIN        ! IN MIN ICE CONCENTRATION IF ICE PRESENT
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
C
C     LOCAL VARIABLES
C
      INTEGER
     + J                   ! LOOP COUNTER
      REAL
     + ONEEM8           ! SMALL +VE VALUE TO ELIMINATE ROUNDING
C
      PARAMETER(ONEEM8 = 1.0E-08 )
C
C
C
C     ----------------------------------------------------------
C
      DO J = 1,L2
C
          IF ( .NOT. LAND(J) ) THEN
C
CL          1. ICE DEPTH.
C
C     CONVERT FROM THE GRID BOX MEAN ACTUAL ICE DEPTH TO THE
C     EQUIVALENT ICE DEPTH AVERAGED OVER ICE AREA.
C     THIS PROCESS USES THE ICE CONCENTRATION AND SNOW DEPTH FIELDS.
C     NOTE THAT AN EXTRA PIECE OF ICE OF DEPTH HICEMIN IS ADDED
C     TO PREVENT VERY SMALL ICE DEPTHS OCCURING WHICH CAN CAUSE
C     FAILURE IN THE ATMOSPHERE MODEL
C
            IF ( AICESLB(J) .LT. ( AICEMIN - ONEEM8 ) ) THEN
              HICEATM(J) = 0.0
            ELSE
              HICEATM(J) = ( HICESLB(J) + HICEMIN )/ AICESLB(J)
     +                     + CONRATIO * HSNOWSLB(J)
            ENDIF
C
CL          2. SNOW DEPTH.
C
C     NOTE THAT THIS HAS TO BE CONVERTED FROM M TO KG/M**2.
C
C     THE SNOW DEPTH OVER LAND POINTS IS NOT SET TO 'RMDI'
C     BECAUSE THIS FIELD CONTAINS LAND SNOW DEPTHS USED BY THE
C     ATMOSPHERE MODEL. THIS IS NOT THE CASE FOR THE OTHER
C     VARIABLES PASSED TO THE ATMOSPHERE, WHICH ARE INITIALISED FOR
C     SAFETY!
C
            HSNOWATM(J) = HSNOWSLB(J) * RHOSNOW
C
CL          3. SEA SURFACE TEMPERATURE.
C
C     NOTE THAT THIS HAS TO BE CONVERTED FROM CELSIUS TO KELVIN.
C
C     AT SEA-ICE POINTS, THE GRID BOX MEAN SURFACE TEMPERATURE IS
C     ALTERED IN SUCH A WAY THAT THE SURFACE TEMPERATURE OF THE ICY
C     PORTION OF THE BOX IS THE SAME AS IT WAS AT THE END OF THE LAST
C     ATMOSPHERIC PHASE. (NOTE - THIS IS NON-CONSERVATIVE)
C     HOWEVER, IF ICE APPEARED DURING THE
C     MOST RECENT OCEAN PHASE, ITS TEMPERATURE IS INITIALISED AT THE
C     FREEZING POINT OF SEAWATER.
C     THIS CODE USES THE OLD VALUES OF ICE CONCENTRATION.
C
            IF ( AICESLB(J) .EQ. 0.0 ) THEN
              TSTARATM(J) = SLABTEMP(J) + ZERODEGC
            ELSEIF ( AICEATM(J) .GE.
     +                     ( AICEMIN - ONEEM8 ) ) THEN
              TSTARATM(J) = TFS + ( AICESLB(J) / AICEATM(J) )
     +                             *( TSTARATM(J) - TFS )
            ELSE
              TSTARATM(J) = TFS
            ENDIF
C
CL    SECTION 4: ICE CONCENTRATION.
C
C     FINALLY, UPDATE THE ICE CONCENTRATION THAT IS PASSED TO
C     THE ATMOSPHERE MODEL.
C
C
            AICEATM(J) = AICESLB(J)
C
          ELSE
C
            HICEATM(J) = RMDI
            AICEATM(J) = RMDI
          ENDIF
      END DO
C
      RETURN
      END
