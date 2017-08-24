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
C
!+ Subroutine to set the solar solar terms in a triple column.
!
! Method:
!       The Direct beam is calculated by propagating down through
!       the column. These direct fluxes are used to  the
!       source terms in each layer.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.2             24-05-96                Original Code
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TRIPLE_SOLAR_SOURCE(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &   , FLUX_INC_DIRECT
     &   , L_SCALE_SOLAR, ADJUST_SOLAR_KE
     &   , TRANS_0, SOURCE_COEFF
     &   , V11, V12, V13, V21, V22, V23, V31, V32, V33
     &   , FLUX_DIRECT
     &   , FLUX_DIRECT_GROUND
     &   , S_UP, S_DOWN
     &   , NPD_PROFILE, NPD_LAYER
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS
      INTEGER   !, INTENT(IN)
     &     NPD_PROFILE
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     COMDECKS INCLUDED
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET INTERNAL DIMENSIONS TIED TO ALGORITHMS,
!     MOSTLY FOR CLOUDS.
!
      INTEGER
     &     NPD_CLOUD_COMPONENT
!             NUMBER OF COMPONENTS OF CLOUDS
     &   , NPD_CLOUD_TYPE
!             NUMBER OF PERMITTED TYPES OF CLOUDS.
     &   , NPD_CLOUD_REPRESENTATION
!             NUMBER OF PERMITTED REPRESENTATIONS OF CLOUDS.
     &   , NPD_OVERLAP_COEFF
!             NUMBER OF OVERLAP COEFFICIENTS FOR CLOUDS
     &   , NPD_SOURCE_COEFF
!             NUMBER OF COEFFICIENTS FOR TWO-STREAM SOURCES
     &   , NPD_REGION
!             NUMBER OF REGIONS IN A LAYER
!
      PARAMETER(
     &     NPD_CLOUD_COMPONENT=4
     &   , NPD_CLOUD_TYPE=4
     &   , NPD_CLOUD_REPRESENTATION=4
     &   , NPD_OVERLAP_COEFF=18
     &   , NPD_SOURCE_COEFF=2
     &   , NPD_REGION=3
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET POINTERS TO SOURCE COEFFICIENTS.
!
      INTEGER
     &     IP_SCF_SOLAR_UP
!             POINTER TO SOURCE COEFICIENT FOR UPWARD SOLAR BEAM
     &   , IP_SCF_SOLAR_DOWN
!             POINTER TO SOURCE COEFICIENT FOR DOWNWARD SOLAR BEAM
     &   , IP_SCF_IR_1D
!             POINTER TO SOURCE COEFICIENT
!             FOR 1ST DIFFERENCE OF PLANCKIAN
     &   , IP_SCF_IR_2D
!             POINTER TO SOURCE COEFICIENT
!             FOR 2ND DIFFERENCE OF PLANCKIAN
!
      PARAMETER(
     &     IP_SCF_SOLAR_UP=1
     &   , IP_SCF_SOLAR_DOWN=2
     &   , IP_SCF_IR_1D=1
     &   , IP_SCF_IR_2D=2
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO DEFINE REFERENCE NUMBERS FOR REGIONS OF CLOUDS.
!
      INTEGER
     &     IP_REGION_CLEAR
!             REFERENCE NUMBER FOR CLEAR-SKY REGION
     &   , IP_REGION_STRAT
!             REFERENCE NUMBER FOR STRATIFORM CLOUDY REGION
     &   , IP_REGION_CONV
!             REFERENCE NUMBER FOR CONVECTIVE CLOUDY REGION
!
      PARAMETER(
     &     IP_REGION_CLEAR=1
     &   , IP_REGION_STRAT=2
     &   , IP_REGION_CONV=3
     &   )
!
!     ------------------------------------------------------------------
!
!
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , N_CLOUD_TOP
!             TOP CLOUDY LAYER
!
!     SPECIAL ARRAYS FOR EQUIVALENT EXTINCTION:
      LOGICAL   !, INTENT(IN)
     &     L_SCALE_SOLAR
!             SCALING APPLIED TO SOLAR FLUX
      REAL      !, INTENT(IN)
     &     ADJUST_SOLAR_KE(NPD_PROFILE, NPD_LAYER)
!             ADJUSTMENT TO SOLAR FLUXES WITH EQUIVALENT EXTINCTION
!
      REAL      !, INTENT(IN)
     &     FLUX_INC_DIRECT(NPD_PROFILE)
!             INCIDENT DIRECT SOLAR FLUX
!
!     OPTICAL PROPERTIES:
      REAL      !, INTENT(IN)
     &     TRANS_0(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             DIRECT TRANSMISSION
     &   , SOURCE_COEFF(NPD_PROFILE, NPD_LAYER
     &      , NPD_SOURCE_COEFF, NPD_REGION)
!             SOURCE COEFFICIENTS
!
!     ENERGY TRANSFER COEFFICIENTS:
      REAL      !, INTENT(IN)
     &     V11(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , V12(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , V13(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , V21(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , V22(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , V23(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , V31(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , V32(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , V33(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
!
!     CALCULATED DIRECT FLUX AND SOURCE TERMS:
      REAL      !, INTENT(OUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)
!             OVERALL DIRECT FLUX
     &   , FLUX_DIRECT_GROUND(NPD_PROFILE, NPD_REGION)
!             DIRECT FLUXES AT GROUND BENEATH EACH REGION
     &   , S_UP(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             UPWARD SOURCE FUNCTIONS
     &   , S_DOWN(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             DOWNWARD SOURCE FUNCTIONS
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     N_REGION
!             NUMBER OF REGIONS
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
     &   , K
!             LOOP VARIABLE
!
      REAL
     &     SOLAR_TOP(NPD_PROFILE, NPD_REGION)
!             SOLAR FLUXES AT TOP OF LAYER
     &   , SOLAR_BASE(NPD_PROFILE, NPD_REGION)
!             SOLAR FLUXES AT BASE OF LAYER
!
!
!     SET THE NUMBER OF REGIONS.
      N_REGION=3
!
!     THE CLEAR AND CLOUDY DIRECT FLUXES ARE CALCULATED SEPARATELY
!     AND ADDED TOGETHER TO FORM THE TOTAL DIRECT FLUX.
!
!     SET INCIDENT FLUXES.
      DO L=1, N_PROFILE
         FLUX_DIRECT(L, 0)=FLUX_INC_DIRECT(L)
      ENDDO
!
!     WITH EQUIVALENT EXTINCTION THE DIRECT SOLAR FLUX MUST BE
!     CORRECTED.
!
      IF (L_SCALE_SOLAR) THEN
!
         DO I=1, N_CLOUD_TOP-1
            DO L=1, N_PROFILE
               FLUX_DIRECT(L, I)
     &            =FLUX_DIRECT(L, I-1)*TRANS_0(L, I, IP_REGION_CLEAR)
     &            *ADJUST_SOLAR_KE(L, I)
               S_UP(L, I, IP_REGION_CLEAR)
     &            =SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP, IP_REGION_CLEAR)
     &            *FLUX_DIRECT(L, I-1)
               S_DOWN(L, I, IP_REGION_CLEAR)
     &            =(SOURCE_COEFF(L, I
     &            , IP_SCF_SOLAR_DOWN, IP_REGION_CLEAR)
     &            -TRANS_0(L, I, IP_REGION_CLEAR))*FLUX_DIRECT(L, I-1)
     &            +FLUX_DIRECT(L, I)
            ENDDO
         ENDDO
!
      ELSE
!
         DO I=1, N_CLOUD_TOP-1
            DO L=1, N_PROFILE
               FLUX_DIRECT(L, I)
     &            =FLUX_DIRECT(L, I-1)*TRANS_0(L, I, IP_REGION_CLEAR)
               S_UP(L, I, IP_REGION_CLEAR)
     &            =SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP, IP_REGION_CLEAR)
     &            *FLUX_DIRECT(L, I-1)
               S_DOWN(L, I, IP_REGION_CLEAR)
     &            =SOURCE_COEFF(L, I
     &            , IP_SCF_SOLAR_DOWN, IP_REGION_CLEAR)
     &            *FLUX_DIRECT(L, I-1)
            ENDDO
         ENDDO
!
      ENDIF
!
!
!
!     CLEAR AND CLOUDY REGION.
!     INITIALIZE PARTIAL FLUXES:
      DO L=1, N_PROFILE
         SOLAR_BASE(L, IP_REGION_CLEAR)=FLUX_DIRECT(L, N_CLOUD_TOP-1)
         SOLAR_BASE(L, IP_REGION_STRAT)=0.0E+00
         SOLAR_BASE(L, IP_REGION_CONV)=0.0E+00
      ENDDO
!
!
      DO I=N_CLOUD_TOP, N_LAYER
!
!        TRANSFER FLUXES ACROSS THE INTERFACE.
!
         DO L=1, N_PROFILE
            SOLAR_TOP(L, IP_REGION_CLEAR)
     &         =V11(L, I-1)*SOLAR_BASE(L, IP_REGION_CLEAR)
     &         +V12(L, I-1)*SOLAR_BASE(L, IP_REGION_STRAT)
     &         +V13(L, I-1)*SOLAR_BASE(L, IP_REGION_CONV)
            SOLAR_TOP(L, IP_REGION_STRAT)
     &         =V21(L, I-1)*SOLAR_BASE(L, IP_REGION_CLEAR)
     &         +V22(L, I-1)*SOLAR_BASE(L, IP_REGION_STRAT)
     &         +V23(L, I-1)*SOLAR_BASE(L, IP_REGION_CONV)
            SOLAR_TOP(L, IP_REGION_CONV)
     &         =V31(L, I-1)*SOLAR_BASE(L, IP_REGION_CLEAR)
     &         +V32(L, I-1)*SOLAR_BASE(L, IP_REGION_STRAT)
     &         +V33(L, I-1)*SOLAR_BASE(L, IP_REGION_CONV)
         ENDDO
!
!
!        PROPAGATE THE FLUXES THROUGH THE LAYER:
!
         IF (L_SCALE_SOLAR) THEN
!
            DO K=1, N_REGION
               DO L=1, N_PROFILE
                  SOLAR_BASE(L, K)
     &               =SOLAR_TOP(L, K)
     &               *TRANS_0(L, I, K)*ADJUST_SOLAR_KE(L, I)
                  S_UP(L, I, K)
     &               =SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP, K)
     &               *SOLAR_TOP(L, K)
                  S_DOWN(L, I, K)
     &               =(SOURCE_COEFF(L, I, IP_SCF_SOLAR_DOWN, K)
     &               -TRANS_0(L, I, K))*SOLAR_TOP(L, K)
     &               +SOLAR_BASE(L, K)
               ENDDO
            ENDDO
!
         ELSE
!
            DO K=1, N_REGION
               DO L=1, N_PROFILE
                  SOLAR_BASE(L, K)=SOLAR_TOP(L, K)
     &               *TRANS_0(L, I, K)
                  S_UP(L, I, K)
     &               =SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP, K)
     &               *SOLAR_TOP(L, K)
                  S_DOWN(L, I, K)
     &               =SOURCE_COEFF(L, I, IP_SCF_SOLAR_DOWN, K)
     &               *SOLAR_TOP(L, K)
               ENDDO
            ENDDO
!
         ENDIF
!
!
!        CALCULATE THE TOTAL DIRECT FLUX.
!
         DO L=1, N_PROFILE
            FLUX_DIRECT(L, I)=SOLAR_BASE(L, IP_REGION_CLEAR)
     &         +SOLAR_BASE(L, IP_REGION_STRAT)
     &         +SOLAR_BASE(L, IP_REGION_CONV)
         ENDDO
!
      ENDDO
!
!     PASS THE LAST VALUE AT THE BASE OF THE CLOUD OUT.
      DO K=1, N_REGION
         DO L=1, N_PROFILE
            FLUX_DIRECT_GROUND(L, K)=SOLAR_BASE(L, K)
         ENDDO
      ENDDO
!
!
!
      RETURN
      END
