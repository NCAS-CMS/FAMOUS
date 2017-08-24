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
!+ Subroutine to calculate transmission and reflection coefficients.
!
! Method:
!        Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             29-03-96                Half-precision
!                                               exponential introduced.
!                                               (J. M. Edwards)
!       4.2             Oct. 96     T3E migration: HF functions
!                                   replaced by T3E vec_lib functions
!                                               (S.J.Swarbrick)
!  4.5  12/05/98  Move constant exp(k*log( )) outside of loop.
!                                            RBarnes@ecmwf.int
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TRANS_SOURCE_COEFF(N_PROFILE
     &   , I_LAYER_FIRST, I_LAYER_LAST
     &   , ISOLIR, L_IR_SOURCE_QUAD
     &   , TAU, SUM, DIFF, LAMBDA, SEC_0
     &   , GAMMA_UP, GAMMA_DOWN
     &   , TRANS, REFLECT, TRANS_0, SOURCE_COEFF
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
!     MODULE TO SET FLAGS FOR DIFFERENT PORTIONS
!     OF THE SPECTRUM.
!
      INTEGER
     &     IP_SOLAR
!             SOLAR REGION
     &   , IP_INFRA_RED
!             INFRA-RED REGION
!
      PARAMETER(
     &     IP_SOLAR=1
     &   , IP_INFRA_RED=2
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE FOR SETTING MACHINE PRECISION.
!
      REAL
     &     TOL_MACHINE
!             MACHINE TOLERANCE
     &   , SQRT_TOL_MACHINE
!             SQRT OF MACHINE TOLERANCE
!
!
!     THE PRECISION SHOULD BE ABOUT 2/2^(SIZE OF SIGNIFICAND)
!
!     THE IEEE-FORMAT USES 53 BITS FOR THE SIGNIFICAND
!     IN DOUBLE PRECISION
!
!     THE CRAY FORMAT USES 47 BITS IN SINGLE PRECISION.
!
      PARAMETER(
     &     TOL_MACHINE=1.42E-14
     &   , SQRT_TOL_MACHINE=1.19E-7
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
!
!     ALGORITHMIC CONTROL
      LOGICAL   !, INTENT(IN)
     &     L_IR_SOURCE_QUAD
!             QUADRATIC SOURCE IN INFRA-RED
      INTEGER   !, INTENT(IN)
     &     ISOLIR
!             SPECTRAL REGION
!
!     OPTICAL PROPERTIES OF THE LAYER
      REAL      !, INTENT(IN)
     &     TAU(NPD_PROFILE, NPD_LAYER)
!             OPTICAL DEPTHS OF LAYERS
     &   , SUM(NPD_PROFILE, NPD_LAYER)
!             SUM OF ALPHA_1 AND ALPHA_2
     &   , DIFF(NPD_PROFILE, NPD_LAYER)
!             DIFFERENCE OF ALPHA_1 AND ALPHA_2
     &   , LAMBDA(NPD_PROFILE, NPD_LAYER)
!             LAMBDA
     &   , SEC_0(NPD_PROFILE)
!             SECANT OF SOLAR ZENITH ANGLE
     &   , GAMMA_UP(NPD_PROFILE, NPD_LAYER)
!             BASIC SOLAR COEFFICIENT FOR UPWARD RADIATION
     &   , GAMMA_DOWN(NPD_PROFILE, NPD_LAYER)
!             BASIC SOLAR COEFFICIENT FOR DOWNWARD RADIATION
!
!     TRANSMISSION AND REFLECTION COEFFICIENTS AND COEFFICIENTS FOR
!     SOURCE TERMS.
      REAL      !, INTENT(OUT)
     &     TRANS(NPD_PROFILE, NPD_LAYER)
!             DIFFUSE TRANSMISSION COEFFICIENT
     &   , REFLECT(NPD_PROFILE, NPD_LAYER)
!             DIFFUSE REFLECTION COEFFICIENT
     &   , TRANS_0(NPD_PROFILE, NPD_LAYER)
!             DIRECT TRANSMISSION COEFFICIENT
     &   , SOURCE_COEFF(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             SOURCE COEFFICIENTS
!
!
!     LOCAL VARIABLES
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
      REAL
     &     GAMMA
!             GAMMA
     &   , EXPONENTIAL
!             EXPONENTIAL OF SCALED OPTICAL DEPTH
      REAL  XLAMTAU(N_PROFILE,I_LAYER_LAST-I_LAYER_FIRST+1) !Workspace
      INTEGER n_input     ! No. of inputs for exp_v
!
!
!
!
!     DETERMINE THE DIFFUSE TRANSMISSION AND REFLECTION COEFFICIENTS.
!
      DO I=I_LAYER_FIRST, I_LAYER_LAST
         DO L=1, N_PROFILE
            XLAMTAU(L,I-I_LAYER_FIRST+1)=-LAMBDA(L,I)*TAU(L,I)
         ENDDO                                
      ENDDO     
      n_input=(I_LAYER_LAST-I_LAYER_FIRST+1)*N_PROFILE
      do I=1,I_LAYER_LAST-I_LAYER_FIRST+1
        do L=1,n_profile
          xlamtau(L,I)=exp(xlamtau(L,I))
        end do
      end do
!                                                                       
      DO I=I_LAYER_FIRST, I_LAYER_LAST                                  
         DO L=1, N_PROFILE                                              
            EXPONENTIAL=xlamtau(L,I-I_LAYER_FIRST+1)
            GAMMA=(SUM(L, I)-LAMBDA(L, I))
     &         /(SUM(L, I)+LAMBDA(L, I))
            TRANS(L, I)=(EXPONENTIAL*(1.0E+00-GAMMA**2)
     &         /(1.0E+00-(EXPONENTIAL*GAMMA)**2))
            REFLECT(L, I)=GAMMA*(1.0E+00-EXPONENTIAL**2)
     &         /(1.0E+00-(EXPONENTIAL*GAMMA)**2)
         ENDDO
      ENDDO
!
!
!
      IF (ISOLIR.EQ.IP_SOLAR) THEN
!
!        CALCULATE THE DIRECT TRANSMISSION AND THE SOURCE COEFFICIENTS
!        FOR THE SOLAR BEAM: IN THE SOLAR CASE THESE ARE
!        THE COEFFICIENTS WHICH WILL MULTIPLY THE DIRECT FLUX AT THE
!        TOP OF THE LAYER TO GIVE THE SOURCE TERMS FOR THE UPWARD
!        DIFFUSE FLUX AND THE TOTAL DOWNWARD FLUX.
!
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               TRANS_0(L, I)=EXP(-TAU(L, I)*SEC_0(L))
               SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP)
     &            =(GAMMA_UP(L, I)-REFLECT(L, I)
     &            *(1.0E+00+GAMMA_DOWN(L, I)))
     &            -GAMMA_UP(L, I)*TRANS(L, I)*TRANS_0(L, I)
               SOURCE_COEFF(L, I, IP_SCF_SOLAR_DOWN)=TRANS_0(L, I)
     &            *(1.0E+00+GAMMA_DOWN(L, I)
     &            -GAMMA_UP(L, I)*REFLECT(L, I))
     &            -(1.0E+00+GAMMA_DOWN(L, I))*TRANS(L, I)
            ENDDO
         ENDDO
!
!
      ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN
!
!        IN THE CASE OF INFRA-RED RADIATION, THE FIRST SOURCE
!        COEFFICIENT HOLDS THE MULTIPLIER FOR THE FIRST DIFFERENCE
!        OF THE PLANCK FUNCTION ACROSS THE LAYER, AND THE SECOND
!        THAT FOR THE SECOND DIFFERENCE.
!
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
!
!              A TOLERANCE IS ADDED TO THE NUMERATOR AND THE DENOMIATOR
!              TO AVOID ILL-CONDITIONING AT SMALL OPTICAL DEPTHS.
!
               SOURCE_COEFF(L, I, IP_SCF_IR_1D)=(1.0E+00-TRANS(L, I)
     &            +REFLECT(L, I)+SQRT_TOL_MACHINE)
     &            /(SQRT_TOL_MACHINE+TAU(L, I)*SUM(L, I))
!
            ENDDO
         ENDDO
!
!
         IF (L_IR_SOURCE_QUAD) THEN
!
!           QUADRATIC CORRECTION TO SOURCE FUNCTION.
!           THIS CORRECTION IS VERY ILL-CONDITIONED FOR
!           SMALL OPTICAL DEPTHS SO THE ASYMPTOTIC FORM IS THEN USED.
!
           EXPONENTIAL = EXP(3.3E-01*LOG(TOL_MACHINE))
            DO I=I_LAYER_FIRST, I_LAYER_LAST
               DO L=1, N_PROFILE
!
!                 USE A SEPARATE ASYMPTOTIC FORM WHEN THE OPTICAL
!                 DEPTH IS SMALL, MAKING THE TRANSITION WHEN THE
!                 OPTICAL DEPTH IS ROUGHLY EQUAL TO THE CUBE ROOT
!                 OF THE MACHINE'S PRECISION.
!
                  IF (TAU(L, I).GT.EXPONENTIAL) THEN
                     SOURCE_COEFF(L, I, IP_SCF_IR_2D)
     &                  =-2.0E+00*(1.0E+00-TRANS(L, I)-REFLECT(L, I)
     &                  +SQRT_TOL_MACHINE)
     &                  /(DIFF(L, I)*TAU(L, I)+SQRT_TOL_MACHINE)
                  ELSE
                     SOURCE_COEFF(L, I, IP_SCF_IR_2D)
     &                  =-2.0E+00+DIFF(L, I)*TAU(L, I)
                  ENDIF
!
                  SOURCE_COEFF(L, I, IP_SCF_IR_2D)
     &               =-(1.0E+00+REFLECT(L, I)+TRANS(L, I)
     &               +SOURCE_COEFF(L, I, IP_SCF_IR_2D))
     &               /(SUM(L, I)*TAU(L, I)+SQRT_TOL_MACHINE)
!
               ENDDO
            ENDDO
!
         ENDIF
!
      ENDIF
!
!
!
      RETURN
      END
