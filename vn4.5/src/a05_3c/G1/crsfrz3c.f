C ******************************COPYRIGHT******************************
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
C
CLL  SUBROUTINE CRS_FRZL------------------------------------------
CLL
CLL  PURPOSE : CHANGE OF PHASE ROUTINE WHERE PRECIPITATION
CLL            CROSSES A MELTING OR FREEZING LEVEL
CLL
CLL            ADD LATENT HEATING
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
!LL   4.4   17/10/97  New version optimised for T3E.
!LL                   Single PE optimisations           D.Salmond
!LL   4.5   03/03/98  Insert missing brackets. Julie Gregory. 
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  LOGICAL COMPONENTS COVERED:
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE CRS_FRZL (NPNTS,RAIN,SNOW,THDD_KM1,EXKM1,
     &                     FLX_DD_KM1,BDDWT_KM1)
C
      IMPLICIT NONE
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
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

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
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
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS                ! IN VECTOR LENGTH
C
      INTEGER I                    ! LOOP COUNTER
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      LOGICAL BDDWT_KM1(NPNTS)     ! IN MASK FOR POINTS WHERE
                                   !    PRECIPITATION IS LIQUID
                                   !    IN LAYER K+1
C
      REAL EXKM1(NPNTS)            ! IN EXNER RATIO FOR LAYER K-1
C
      REAL FLX_DD_KM1(NPNTS)       ! IN MASS FLUX OF LAYER K-1 (PA/S)
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C----------------------------------------------------------------------
C
      REAL RAIN(NPNTS)             ! INOUT
                                   ! IN  AMOUNT OF FALLING RAIN
                                   !     DESCENDING FROM LAYER
                                   !     K-1 TO K-2 (KG/M**2/S)
                                   ! OUT UPDATED AMOUNT OF FALLING
                                   !     RAIN (KG/M**2/S)
C
      REAL SNOW(NPNTS)             ! INOUT
                                   ! IN  AMOUNT OF FALLING SNOW
                                   !     DESCENDING FROM LAYER
                                   !     K-1 TO K-2 (KG/M**2/S)
                                   ! OUT UPDATED AMOUNT OF FALLING
                                   !     SNOW (KG/M**2/S)
C
      REAL THDD_KM1(NPNTS)         ! INOUT
                                   ! IN  DOWNDRAUGHT POTENTIAL
                                   !     TEMPERATURE IN LAYER K-1 (K)
                                   ! OUT UPDATED DOWNDRAUGHT POTENTIAL
                                   !     TEMPERATURE IN LAYER K-1
                                   !     DUE TO CHANGE OF PHASE (K)
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C---------------------------------------------------------------------
C
      REAL FACTOR                  ! USED IN THE CALCULATION OF
                                   ! CHANGE OF PHASE OF FALLING
                                   ! PRECIPITATION
C
      REAL PRECIP_FRE              ! FREEZING PRECIPITATION
C
      REAL PRECIP_MELT             ! MELTING PRECIPITATION
C
CL
CL----------------------------------------------------------------------
CL  ADD LATENT HEATING WHERE PRECIP CROSSES A MELTING OR FREEZING LEVEL
CL
CL  UM DOCUMENTATION PAPER 27
CL  SECTION (11), EQUATION (42)
CL----------------------------------------------------------------------
CL
       DO I=1,NPNTS
C
        IF (.NOT.BDDWT_KM1(I).AND.RAIN(I).GT.0.0.AND.THDD_KM1(I)
     *       *EXKM1(I).LT.TM) THEN
C FREEZE
          FACTOR = (EXKM1(I)*CP*FLX_DD_KM1(I))/(LF*G)
          PRECIP_FRE = (TM/EXKM1(I)-THDD_KM1(I))* FACTOR
          PRECIP_FRE = MIN(RAIN(I),PRECIP_FRE)
          THDD_KM1(I) = THDD_KM1(I)+PRECIP_FRE/FACTOR
          RAIN(I) = RAIN(I)-PRECIP_FRE
          SNOW(I) = SNOW(I)+PRECIP_FRE
C
        ELSE IF (BDDWT_KM1(I).AND.SNOW(I).GT.0.0) THEN
C MELT
          FACTOR = (EXKM1(I)*CP*FLX_DD_KM1(I))/(LF*G)
          PRECIP_MELT = (THDD_KM1(I)-TM/EXKM1(I))*FACTOR
          PRECIP_MELT = MIN(SNOW(I),PRECIP_MELT)
          THDD_KM1(I) = THDD_KM1(I)-PRECIP_MELT/FACTOR
          RAIN(I) = RAIN(I)+PRECIP_MELT
          SNOW(I) = SNOW(I)-PRECIP_MELT
        END IF
      END DO
C
      RETURN
      END
C
