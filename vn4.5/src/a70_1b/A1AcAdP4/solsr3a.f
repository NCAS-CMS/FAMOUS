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
!+ Subroutine to calculate the solar flux and source terms.
!
! Method:
!       Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             08-05-97                Formulation for
!                                               equivalent extinction
!                                               amended.
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLAR_SOURCE(N_PROFILE, N_LAYER
     &   , FLUX_INC_DIRECT
     &   , TRANS_0, SOURCE_COEFF
     &   , L_SCALE_SOLAR, ADJUST_SOLAR_KE
     &   , FLUX_DIRECT
     &   , S_DOWN, S_UP
     &   , NPD_PROFILE, NPD_LAYER
     &   )
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
!     COMDECKS INCLUDED.
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
     &   , N_LAYER
!             NUMBER OF LAYERS
!
      LOGICAL   !, INTENT(IN)
     &     L_SCALE_SOLAR
!             SCALING APPLIED TO SOLAR BEAM
!
      REAL      !, INTENT(IN)
     &     FLUX_INC_DIRECT(NPD_PROFILE)
!             INCIDENT SOLAR FLUX
     &   , TRANS_0(NPD_PROFILE, NPD_LAYER)
!             DIRECT TRANSMISSION COEFFICIENT
     &   , SOURCE_COEFF(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             REFLECTION COEFFICIENT
     &   , ADJUST_SOLAR_KE(NPD_PROFILE, NPD_LAYER)
!             ADJUSTMENT TO SOLAR FLUX
!
!
      REAL      !, INTENT(OUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)
!             DIRECT FLUX
     &   , S_DOWN(NPD_PROFILE, NPD_LAYER)
!             DOWNWARD SOURCE FUNCTION
     &   , S_UP(NPD_PROFILE, NPD_LAYER)
!             UPWARD SOURCE FUNCTION
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
!
!
      DO L=1, N_PROFILE
         FLUX_DIRECT(L, 0)=FLUX_INC_DIRECT(L)
      ENDDO
!
!     THE SOLAR FLUX MAY BE MULTIPLIED BY A SCALING FACTOR IF AN
!     EQUIVALENT EXTINCTION IS USED.
!
      IF (L_SCALE_SOLAR) THEN
!
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               FLUX_DIRECT(L, I)
     &            =FLUX_DIRECT(L, I-1)*TRANS_0(L, I)
     &            *ADJUST_SOLAR_KE(L, I)
               S_UP(L, I)=SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP)
     &            *FLUX_DIRECT(L, I-1)
               S_DOWN(L, I)=(SOURCE_COEFF(L, I, IP_SCF_SOLAR_DOWN)
     &            -TRANS_0(L, I))*FLUX_DIRECT(L, I-1)
     &            +FLUX_DIRECT(L, I)
            ENDDO
         ENDDO
!
      ELSE
!
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               FLUX_DIRECT(L, I)
     &            =FLUX_DIRECT(L, I-1)*TRANS_0(L, I)
               S_UP(L, I)=SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP)
     &            *FLUX_DIRECT(L, I-1)
               S_DOWN(L, I)=SOURCE_COEFF(L, I, IP_SCF_SOLAR_DOWN)
     &            *FLUX_DIRECT(L, I-1)
            ENDDO
         ENDDO
!
      ENDIF
!
!
      RETURN
      END
