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
!+ Subroutine to calculate fluxes including only gaseous absorption.
!
! Method:
!       Transmission coefficients for each layer are calculated
!       from the gaseous absorption alone. Fluxes are propagated
!       upward or downward through the column using these
!       coefficients and source terms.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             29-03-96                Half-precision
!                                               exponentials introduced.
!                                               (J. M. Edwards)
!       4.2             Oct. 96     T3E migration: HF functions
!                                   replaced by T3E vec_lib function
!                                               (S.J.Swarbrick)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE MONOCHROMATIC_GAS_FLUX(N_PROFILE, N_LAYER
     &   , L_NET
     &   , TAU_GAS
     &   , ISOLIR, SEC_0, FLUX_INC_DIRECT, FLUX_INC_DOWN
     &   , DIFF_PLANCK, SOURCE_GROUND
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR
     &   , DIFFUSIVITY_FACTOR
     &   , FLUX_DIRECT, FLUX_DIFFUSE
     &   , NPD_PROFILE, NPD_LAYER
     &         )
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
!     INCLUDE COMDECKS
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
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , ISOLIR
!             SPECTRAL REGION
      LOGICAL   !, INTENT(IN)
     &     L_NET
!             CALCULATE NET FLUXES.
      REAL      !, INTENT(IN)
     &     TAU_GAS(NPD_PROFILE, NPD_LAYER)
!             GASEOUS OPTICAL DEPTHS
     &   , SEC_0(NPD_PROFILE)
!             SECANT OF ZENITH ANGLE
     &   , FLUX_INC_DIRECT(NPD_PROFILE)
!             INCIDENT DIRECT FLUX
     &   , FLUX_INC_DOWN(NPD_PROFILE)
!             INCIDENT DIFFUSE FLUX
     &   , SOURCE_GROUND(NPD_PROFILE)
!             SOURCE FUNCTION OF GROUND
     &   , DIFF_PLANCK(NPD_PROFILE, NPD_LAYER)
!             DIFFERENCE IN PLANCK FUNCTION
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)
!             DIFFUSE SURFACE ALBEDO
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)
!             DIRECT SURFACE ALBEDO
     &   , DIFFUSIVITY_FACTOR
!             DIFFUSIVITY FACTOR
      REAL      !, INTENT(OUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)
!             DIRECT FLUX
     &   , FLUX_DIFFUSE(NPD_PROFILE, 2*NPD_LAYER+2)
!             DIFFUSE FLUX
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
      REAL
     &     TRANS(N_PROFILE, N_LAYER)
!             TRANSMISSIVITIES
      REAL
     &     SOURCE_UP(NPD_PROFILE, NPD_LAYER)
!             UPWARD SOURCE FUNCTION
     &   , SOURCE_DOWN(NPD_PROFILE, NPD_LAYER)
!             DOWNWARD SOURCE FUNCTION
!
!
!
!                                                                       
      DO I=1, N_LAYER                                                   
         DO L=1, N_PROFILE                                              
            TRANS(L, I)=EXP(-DIFFUSIVITY_FACTOR*TAU_GAS(L, I))    
         ENDDO                                                          
      ENDDO                                                             
!
      IF (ISOLIR.EQ.IP_SOLAR) THEN
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               SOURCE_UP(L, I)=0.0E+00
               SOURCE_DOWN(L, I)=0.0E+00
            ENDDO
         ENDDO
      ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               SOURCE_UP(L, I)=(1.0E+00-TRANS(L, I)+SQRT_TOL_MACHINE)
     &            *DIFF_PLANCK(L, I)
     &            /(DIFFUSIVITY_FACTOR*TAU_GAS(L, I)+SQRT_TOL_MACHINE)
               SOURCE_DOWN(L, I)=-SOURCE_UP(L, I)
            ENDDO
         ENDDO
      ENDIF
!
!     THE DIRECT FLUX.
      IF (ISOLIR.EQ.IP_SOLAR) THEN
         DO L=1, N_PROFILE
            FLUX_DIRECT(L, 0)=FLUX_INC_DIRECT(L)
         ENDDO
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               FLUX_DIRECT(L, I)
     &            =FLUX_DIRECT(L, I-1)*EXP(-TAU_GAS(L, I)*SEC_0(L))
            ENDDO
         ENDDO
      ENDIF
!
!     DOWNWARD FLUXES.
      DO L=1, N_PROFILE
         FLUX_DIFFUSE(L, 2)=FLUX_INC_DOWN(L)
      ENDDO
      DO I=1, N_LAYER
         DO L=1, N_PROFILE
            FLUX_DIFFUSE(L, 2*I+2)=TRANS(L, I)*FLUX_DIFFUSE(L, 2*I)
     &         +SOURCE_DOWN(L, I)
         ENDDO
      ENDDO
!
!     UPWARD FLUXES.
      IF (ISOLIR.EQ.IP_SOLAR) THEN
         DO L=1, N_PROFILE
            FLUX_DIFFUSE(L, 2*N_LAYER+1)=SOURCE_GROUND(L)
     &         +ALBEDO_SURFACE_DIFF(L)*FLUX_DIFFUSE(L, 2*N_LAYER+2)
     &         +ALBEDO_SURFACE_DIR(L)*FLUX_DIRECT(L, N_LAYER)
         ENDDO
      ELSE
         DO L=1, N_PROFILE
            FLUX_DIFFUSE(L, 2*N_LAYER+1)=SOURCE_GROUND(L)
     &         +ALBEDO_SURFACE_DIFF(L)*FLUX_DIFFUSE(L, 2*N_LAYER+2)
         ENDDO
      ENDIF
      DO I=N_LAYER, 1, -1
         DO L=1, N_PROFILE
            FLUX_DIFFUSE(L, 2*I-1)=TRANS(L, I)*FLUX_DIFFUSE(L, 2*I+1)
     &         +SOURCE_UP(L, I)
         ENDDO
      ENDDO
!
!     REDUCE TO THE NET FLUX IF THIS IS REQUIRED.
      IF (L_NET) THEN
         DO I=0, N_LAYER
            DO L=1, N_PROFILE
               FLUX_DIFFUSE(L, I+1)=FLUX_DIFFUSE(L, 2*I+2)
     &            -FLUX_DIFFUSE(L, 2*I+1)
            ENDDO
         ENDDO
      ENDIF
!
!
      RETURN
      END
