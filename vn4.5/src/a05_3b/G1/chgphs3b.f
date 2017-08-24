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
CLL  SUBROUTINE CHG_PHSE----------------------------------------------
CLL
CLL  PURPOSE : CHANGE OF PHASE ROUTINE FOR POINTS WHERE NO
CLL            DOWNDRAUGHT OCCURING
CLL
CLL            UPDATES POTENTIAL TEMPERATURE OF LAYER K
CLL            AS PRECIPITATION CHANGES PHASE IN SITU
CLL
CLL            ADD LATENT HEATING WHERE PRECIPITATION CROSSES A
CLL            MELTING OR FREEZING LEVEL
CLL            Limit the freezing and melting of convective
CLL            precipitation below cloud base so that the change
CLL            of phase cannot take the temperature to the
CLL            other side of the melting point of water (TM).
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  Code written by R.N.B.Smith (based on CHGPHS2A at vn4.2)  Feb.97
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  PROJECT TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE CHG_PHSE (NPNTS,K,RAIN,SNOW,DTHBYDT_KM1,
     &                     EXK,EXKM1,DELPKM1,THE_K,THE_KM1,
     &                     TIMESTEP,CCA)
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

      REAL CLDAREA                ! FRACTIONAL CLOUD AREA WHEN NO
                                  ! DD, USED IN EVAPORATION CALC
      PARAMETER (CLDAREA=1.0)
C
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS                ! IN VECTOR LENGTH
C
      INTEGER I                    ! LOOP COUNTER
C
      INTEGER K                    ! IN MODEL LAYER
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL EXK(NPNTS)              ! IN EXNER RATIO FOR LAYER K
C
      REAL EXKM1(NPNTS)            ! IN EXNER RATIO FOR LAYER K-1
C
      REAL DELPKM1(NPNTS)          ! IN PRESSURE DIFFERENCE ACROSS
                                   !    LAYER K-1 (PA)
C
      REAL THE_K(NPNTS)            ! IN POTENTIAL TEMPERATURE OF
                                   !    ENVIRONMENT IN LAYER K
C
      REAL THE_KM1(NPNTS)          ! IN POTENTIAL TEMPERATURE OF
                                   !    ENVIRONMENT IN LAYER K-1
C
      REAL TIMESTEP                ! IN Timestep in seconds
C
      REAL CCA(NPNTS)              ! IN Convective cloud area
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C----------------------------------------------------------------------
C
      REAL RAIN(NPNTS)             ! INOUT
                                   ! IN  AMOUNT OF FALLING RAIN
                                   !     (KG/M**2/S)
                                   ! OUT UPDATED AMOUNT OF FALLING
                                   !     RAIN (KG/M**2/S)
C
      REAL SNOW(NPNTS)             ! INOUT
                                   ! IN  AMOUNT OF FALLING SNOW
                                   !     (KG/M**2/S)
                                   ! OUT UPDATED AMOUNT OF FALLING
                                   !     SNOW (KG/M**2/S)
C
      REAL DTHBYDT_KM1(NPNTS)      ! INOUT
                                   ! IN  INCREMENT TO MODEL POTENTIAL
                                   !     TEMPERATURE IN LAYER K-1
                                   ! OUT UPDATED INCREMENT TO MODEL
                                   !     POTENTIAL TEMPERATURE IN LAYER
                                   !     K-1 DUE TO CHANGE OF PHASE
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C---------------------------------------------------------------------
C
      REAL FACTOR                  ! USED IN THE CALCULATION OF
                                   ! CHANGE OF PHASE OF FALLING
                                   ! PRECIPITATION
C
      REAL WPC                     ! AMOUNT OF PRECIPITATION WHICH CAN
C                                  ! CHANGE PHASE
C
      REAL CA                      ! LOCAL CLOUD AREA
C
      REAL THE_KM1_NEW             ! THE_KM1 UPDATED WITH ALL INCREMENTS
C                                  ! PRIOR TO THIS SUBROUTINE
C
      LOGICAL BPPNWT_K             ! MASK WHERE PRECIPITATION IS LIQUID
                                   ! IN LAYER K
C
      LOGICAL BPPNWT_KM1           ! MASK WHERE PRECIPITATION IS LIQUID
                                   ! IN LAYER K-1
C
CL
CL----------------------------------------------------------------------
CL  ADD LATENT HEATING WHERE PRECIP CROSSES A MELTING OR FREEZING LEVEL
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (11), EQUATION (42)
CL----------------------------------------------------------------------
CL
      DO I=1,NPNTS
        THE_KM1_NEW = THE_KM1(I) + TIMESTEP*DTHBYDT_KM1(I)
        BPPNWT_K = THE_K(I).GT.TM/EXK(I)
        BPPNWT_KM1 = THE_KM1_NEW .GT. TM/EXKM1(I)
        FACTOR = LF*G/(EXKM1(I)*CP*DELPKM1(I))
        CA = CLDAREA * CCA(I)
C
C FREEZE
C
        IF (.NOT.BPPNWT_KM1.AND.(BPPNWT_K.OR.RAIN(I).GT.0.0)) THEN
          WPC = MIN( RAIN(I),
     &               CA * (TM/EXKM1(I) - THE_KM1_NEW) /
     &                    (TIMESTEP * FACTOR) )
          DTHBYDT_KM1(I) = DTHBYDT_KM1(I) + WPC * FACTOR
          SNOW(I) = SNOW(I) + WPC
          RAIN(I) = RAIN(I) - WPC
        END IF
C
C MELT
C
        IF (BPPNWT_KM1.AND.(.NOT.BPPNWT_K.OR.SNOW(I).GT.0.0)) THEN
          WPC = MIN( SNOW(I),
     &               CA * (THE_KM1_NEW - TM/EXKM1(I)) /
     &                    (TIMESTEP * FACTOR) )
          DTHBYDT_KM1(I) = DTHBYDT_KM1(I) - WPC * FACTOR
          RAIN(I) = RAIN(I) + WPC
          SNOW(I) = SNOW(I) - WPC
        END IF
      END DO
C
      RETURN
      END
C
