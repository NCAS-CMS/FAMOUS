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
!+ Subroutine to set the solar source terms in a mixed column.
!
! Method:
!       The Direct beam is calculated by propagating down through
!       the column. These direct fluxes are used to define the
!       source terms in each layer.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             14-05-96                Cloudy direct flux
!                                               at ground returned
!                                               for new solver.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE MIXED_SOLAR_SOURCE(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &   , FLUX_INC_DIRECT
     &   , L_SCALE_SOLAR, ADJUST_SOLAR_KE
     &   , TRANS_0_FREE, SOURCE_COEFF_FREE
     &   , G_FF, G_FC, G_CF, G_CC
     &   , TRANS_0_CLOUD, SOURCE_COEFF_CLOUD
     &   , FLUX_DIRECT
     &   , FLUX_DIRECT_GROUND_CLOUD
     &   , S_UP_FREE, S_DOWN_FREE
     &   , S_UP_CLOUD, S_DOWN_CLOUD
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
!     CLEAR-SKY OPTICAL PROPERTIES:
      REAL      !, INTENT(IN)
     &     TRANS_0_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE DIRECT TRANSMISSION
     &   , SOURCE_COEFF_FREE(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             CLEAR-SKY SOURCE COEFFICIENTS
!
!     CLOUDY OPTICAL PROPERTIES:
      REAL      !, INTENT(IN)
     &     TRANS_0_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY TRANSMISSION
     &   , SOURCE_COEFF_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             CLOUDY REFLECTANCE
!
!     ENERGY TRANSFER COEFFICIENTS:
      REAL      !, INTENT(IN)
     &     G_FF(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , G_FC(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , G_CF(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , G_CC(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
!
!     CALCULATED DIRECT FLUX AND SOURCE TERMS:
      REAL      !, INTENT(OUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)
!             DIRECT FLUX
     &   , FLUX_DIRECT_GROUND_CLOUD(NPD_PROFILE)
!             DIRECT CLOUDY FLUX AT GROUND
     &   , S_UP_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE UPWARD SOURCE FUNCTION
     &   , S_DOWN_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE DOWNWARD SOURCE FUNCTION
     &   , S_UP_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY UPWARD SOURCE FUNCTION
     &   , S_DOWN_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY DOWNWARD SOURCE FUNCTION
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
      REAL
     &     SOLAR_TOP_FREE(NPD_PROFILE)
!             FREE SOLAR FLUX AT TOP OF LAYER
     &   , SOLAR_TOP_CLOUD(NPD_PROFILE)
!             CLOUDY SOLAR FLUX AT TOP OF LAYER
     &   , SOLAR_BASE_FREE(NPD_PROFILE)
!             FREE SOLAR FLUX AT BASE OF LAYER
     &   , SOLAR_BASE_CLOUD(NPD_PROFILE)
!             CLOUDY SOLAR FLUX AT BASE OF LAYER
!
!
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
     &            =FLUX_DIRECT(L, I-1)*TRANS_0_FREE(L, I)
     &            *ADJUST_SOLAR_KE(L, I)
               S_UP_FREE(L, I)=SOURCE_COEFF_FREE(L, I, IP_SCF_SOLAR_UP)
     &            *FLUX_DIRECT(L, I-1)
               S_DOWN_FREE(L, I)
     &            =(SOURCE_COEFF_FREE(L, I, IP_SCF_SOLAR_DOWN)
     &            -TRANS_0_FREE(L, I))*FLUX_DIRECT(L, I-1)
     &            +FLUX_DIRECT(L, I)
            ENDDO
         ENDDO
!
      ELSE
!
         DO I=1, N_CLOUD_TOP-1
            DO L=1, N_PROFILE
               FLUX_DIRECT(L, I)
     &            =FLUX_DIRECT(L, I-1)*TRANS_0_FREE(L, I)
               S_UP_FREE(L, I)=SOURCE_COEFF_FREE(L, I, IP_SCF_SOLAR_UP)
     &            *FLUX_DIRECT(L, I-1)
               S_DOWN_FREE(L, I)
     &            =SOURCE_COEFF_FREE(L, I, IP_SCF_SOLAR_DOWN)
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
         SOLAR_BASE_FREE(L)=FLUX_DIRECT(L, N_CLOUD_TOP-1)
         SOLAR_BASE_CLOUD(L)=0.0E+00
      ENDDO
!
!
      DO I=N_CLOUD_TOP, N_LAYER
!
!        TRANSFER FLUXES ACROSS THE INTERFACE. THE USE OF ONLY ONE
!        CLOUDY FLUX IMPLICITLY FORCES RANDOM OVERLAP OF DIFFERENT
!        SUBCLOUDS WITHIN THE CLOUDY PARTS OF THE LAYER.
!
         DO L=1, N_PROFILE
            SOLAR_TOP_CLOUD(L)=G_CC(L, I-1)*SOLAR_BASE_CLOUD(L)
     &         +G_FC(L, I-1)*SOLAR_BASE_FREE(L)
            SOLAR_TOP_FREE(L)=G_FF(L, I-1)*SOLAR_BASE_FREE(L)
     &         +G_CF(L, I-1)*SOLAR_BASE_CLOUD(L)
         ENDDO
!
!
!        PROPAGATE THE CLEAR AND CLOUDY FLUXES THROUGH THE LAYER:
!
         IF (L_SCALE_SOLAR) THEN
!
            DO L=1, N_PROFILE
               SOLAR_BASE_FREE(L)=SOLAR_TOP_FREE(L)
     &            *TRANS_0_FREE(L, I)*ADJUST_SOLAR_KE(L, I)
               SOLAR_BASE_CLOUD(L)=SOLAR_TOP_CLOUD(L)
     &            *TRANS_0_CLOUD(L, I)*ADJUST_SOLAR_KE(L, I)
               S_UP_FREE(L, I)=SOURCE_COEFF_FREE(L, I, IP_SCF_SOLAR_UP)
     &            *SOLAR_TOP_FREE(L)
               S_DOWN_FREE(L, I)
     &            =(SOURCE_COEFF_FREE(L, I, IP_SCF_SOLAR_DOWN)
     &            -TRANS_0_FREE(L, I))*SOLAR_TOP_FREE(L)
     &            +SOLAR_BASE_FREE(L)
               S_UP_CLOUD(L, I)
     &            =SOURCE_COEFF_CLOUD(L, I, IP_SCF_SOLAR_UP)
     &            *SOLAR_TOP_CLOUD(L)
               S_DOWN_CLOUD(L, I)
     &            =(SOURCE_COEFF_CLOUD(L, I, IP_SCF_SOLAR_DOWN)
     &            -TRANS_0_CLOUD(L, I))*SOLAR_TOP_CLOUD(L)
     &            +SOLAR_BASE_CLOUD(L)
            ENDDO
!
         ELSE
!
            DO L=1, N_PROFILE
               SOLAR_BASE_FREE(L)=SOLAR_TOP_FREE(L)
     &            *TRANS_0_FREE(L, I)
               SOLAR_BASE_CLOUD(L)=SOLAR_TOP_CLOUD(L)
     &            *TRANS_0_CLOUD(L, I)
               S_UP_FREE(L, I)=SOURCE_COEFF_FREE(L, I, IP_SCF_SOLAR_UP)
     &            *SOLAR_TOP_FREE(L)
               S_DOWN_FREE(L, I)
     &            =SOURCE_COEFF_FREE(L, I, IP_SCF_SOLAR_DOWN)
     &            *SOLAR_TOP_FREE(L)
               S_UP_CLOUD(L, I)
     &            =SOURCE_COEFF_CLOUD(L, I, IP_SCF_SOLAR_UP)
     &            *SOLAR_TOP_CLOUD(L)
               S_DOWN_CLOUD(L, I)
     &            =SOURCE_COEFF_CLOUD(L, I, IP_SCF_SOLAR_DOWN)
     &            *SOLAR_TOP_CLOUD(L)
            ENDDO
!
         ENDIF
!
!
!        CALCULATE THE TOTAL DIRECT FLUX.
!
         DO L=1, N_PROFILE
            FLUX_DIRECT(L, I)=SOLAR_BASE_FREE(L)+SOLAR_BASE_CLOUD(L)
         ENDDO
!
      ENDDO
!
!     PASS THE LAST VALUE AT THE BASE OF THE CLOUD OUT.
      DO L=1, N_PROFILE
         FLUX_DIRECT_GROUND_CLOUD(L)=SOLAR_BASE_CLOUD(L)
      ENDDO
!
!
!
      RETURN
      END
