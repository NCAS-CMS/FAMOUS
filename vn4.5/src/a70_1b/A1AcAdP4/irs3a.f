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
!+ Subroutine to calcaulate IR source function for differential flux.
!
! Method:
!       The linear contribution to the source function is proportional
!       to the absorption divided by the optical depth. A tolerance is
!       added to the optical depth to allow for the depth's being 0.
!       Corrections may also be made for cwa quadratic variation in the
!       temperature across the layer and for the effects of edges.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE IR_SOURCE(N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST
     &   , SOURCE_COEFF, DEL_PLANCK, L_IR_SOURCE_QUAD, DIFF_PLANCK_2
     &   , L_2_STREAM_CORRECT, PLANCK_SOURCE
     &   , GROUND_EMISSION, N_LAYER
     &   , TAU, TRANS
     &   , S_DOWN, S_UP
     &   , NPD_PROFILE, NPD_LAYER
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS.
      INTEGER   !, INTENT(IN)
     &     NPD_PROFILE
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     INCLUDE COMDECKS.
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
!     DUMMY VARIABLES.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , I_LAYER_FIRST
!             FIRST LAYER TO CONSIDER
     &   , I_LAYER_LAST
!             LAST LAYER TO CONSIDER
     &   , N_LAYER
!             NUMBER OF LAYERS
!
      LOGICAL   !, INTENT(IN)
     &     L_IR_SOURCE_QUAD
!             USE A QUADRATIC REPRESENTATION
     &   , L_2_STREAM_CORRECT
!             EDGE CORRECTION TO 2-STREAM
!
      REAL      !, INTENT(IN)
     &     SOURCE_COEFF(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             COEFFICIENTS FOR SOURCE TERMS
     &   , DEL_PLANCK(NPD_PROFILE, NPD_LAYER)
!             DIFFERENCE IN PLANCK FUNCTION ACROSS THE LAYER
     &   , DIFF_PLANCK_2(NPD_PROFILE, NPD_LAYER)
!             2x2ND DIFFERENCE OF PLANCKIAN
     &   , TAU(NPD_PROFILE, NPD_LAYER)
!             OPTCIAL DEPTH
     &   , TRANS(NPD_PROFILE, NPD_LAYER)
!             TRANSMISSION COEFFICIENT
     &   , PLANCK_SOURCE(NPD_PROFILE, 0: NPD_LAYER)
!             PLANCKIAN SOURCE FUNCTION
     &   , GROUND_EMISSION(NPD_PROFILE)
!             TOTAL FLUX EMITTED FROM GROUND
!
      REAL      !, INTENT(OUT)
     &     S_DOWN(NPD_PROFILE, NPD_LAYER)
!             UPWARD SOURCE FUNCTION
     &   , S_UP(NPD_PROFILE, NPD_LAYER)
!             UPWARD SOURCE FUNCTION
!
!
!     LOCAL VARIABLES.
!
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
      REAL
     &     TAUC(NPD_PROFILE, 0: NPD_LAYER)
!             CUMULATIVE OPTICAL DEPTH
     &   , PLANCK_AVE(NPD_PROFILE, 0: NPD_LAYER)
!             AVERAGE PLANCKIAN
     &   , DELTA_TAU_UP_TOP
!             OPTICAL DEPTH: SURF-TOP OF LAYER
     &   , DELTA_TAU_UP_BASE
!             OPTICAL DEPTH: SURF-BASE OF LAYER
!
!     FUNCTIONS CALLED:
      REAL
     &     E3_ACC01
!             THIRD EXPONENTIAL INTEGRAL TO 1%
      EXTERNAL
     &     E3_ACC01
!
!
!
!     MULTIPLY THE SOURCE COEFFICIENTS BY THE PLANCKIAN DIFFERENCES
!     TO THE ORDER REQUIRED.
!
      IF (L_IR_SOURCE_QUAD) THEN
!
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               S_UP(L, I)=SOURCE_COEFF(L, I, IP_SCF_IR_1D)
     &            *DEL_PLANCK(L, I)
     &            +SOURCE_COEFF(L, I, IP_SCF_IR_2D)
     &            *DIFF_PLANCK_2(L, I)
               S_DOWN(L, I)=-SOURCE_COEFF(L, I, IP_SCF_IR_1D)
     &            *DEL_PLANCK(L, I)
     &            +SOURCE_COEFF(L, I, IP_SCF_IR_2D)
     &            *DIFF_PLANCK_2(L, I)
            ENDDO
!
         ENDDO
!
      ELSE
!
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               S_UP(L, I)=SOURCE_COEFF(L, I, IP_SCF_IR_1D)
     &            *DEL_PLANCK(L, I)
               S_DOWN(L, I)=-S_UP(L, I)
            ENDDO
         ENDDO
!
      ENDIF
!
!
!     EDGE CORRECTIONS TO 2-STREAM EQUATIONS.
!
      IF (L_2_STREAM_CORRECT) THEN
!
         DO L=1, N_PROFILE
            TAUC(L, 0)=0.0E+00
         ENDDO
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               TAUC(L, I)=TAUC(L, I-1)+TAU(L, I)
               PLANCK_AVE(L, I)
     &            =0.5E+00*(PLANCK_SOURCE(L, I-1)+PLANCK_SOURCE(L, I))
            ENDDO
         ENDDO
!
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               DELTA_TAU_UP_TOP=TAUC(L, N_LAYER)-TAUC(L, I-1)
               DELTA_TAU_UP_BASE=TAUC(L, N_LAYER)-TAUC(L, I)
               S_UP(L, I)=S_UP(L, I)
     &            +2.0E+00*(GROUND_EMISSION(L)-PLANCK_AVE(L, I))
     &            *(E3_ACC01(DELTA_TAU_UP_TOP)
     &            -TRANS(L, I)*E3_ACC01(DELTA_TAU_UP_BASE))
               S_DOWN(L, I)=S_DOWN(L, I)
     &            +2.0E+00*PLANCK_AVE(L, I)
     &            *(TRANS(L, I)*E3_ACC01(TAUC(L, I-1))
     &            -E3_ACC01(TAUC(L, I)))
            ENDDO
         ENDDO

      ENDIF
!
!
!
      RETURN
      END
