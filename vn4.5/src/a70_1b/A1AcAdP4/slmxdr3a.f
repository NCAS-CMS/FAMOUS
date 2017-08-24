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
!+ Subroutine to solve for mixed fluxes scattering without a matrix.
!
! Method:
!       Gaussian elimination in an upward direction is employed to
!       determine effective albedos for lower levels of the atmosphere.
!       This allows a downward pass of back-substitution to be carried
!       out to determine the upward and downward fluxes.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.1             10-04-96                Original Code
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLVER_MIX_DIRECT(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &   , T, R, S_DOWN, S_UP
     &   , T_CLOUD, R_CLOUD, S_DOWN_CLOUD, S_UP_CLOUD
     &   , V11, V21, V12, V22
     &   , U11, U12, U21, U22
     &   , L_NET
     &   , FLUX_INC_DOWN
     &   , SOURCE_GROUND_FREE, SOURCE_GROUND_CLOUD, ALBEDO_SURFACE_DIFF
     &   , FLUX_TOTAL
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
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , N_CLOUD_TOP
!             TOPMOST CLOUDY LAYER
      LOGICAL   !, INTENT(IN)
     &     L_NET
!             FLAG FOR CALCULATION OF NET FLUXES
      REAL      !, INTENT(IN)
     &     T(NPD_PROFILE, NPD_LAYER)
!             CLEAR-SKY TRANSMISSION
     &   , R(NPD_PROFILE, NPD_LAYER)
!             CLEAR-SKY REFLECTION
     &   , S_DOWN(NPD_PROFILE, NPD_LAYER)
!             CLEAR-SKY DOWNWARD SOURCE FUNCTION
     &   , S_UP(NPD_PROFILE, NPD_LAYER)
!             CLEAR-SKY UPWARD SOURCE FUNCTION
     &   , T_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY TRANSMISSION
     &   , R_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY REFLECTION
     &   , S_DOWN_CLOUD(NPD_PROFILE, NPD_LAYER)
!             DOWNWARD CLOUDY SOURCE FUNCTION
     &   , S_UP_CLOUD(NPD_PROFILE, NPD_LAYER)
!             UPWARD CLOUDY SOURCE FUNCTION
      REAL      !, INTENT(IN)
     &     V11(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , V21(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , V12(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , V22(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , U11(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , U12(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , U21(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , U22(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
      REAL      !, INTENT(IN)
     &     FLUX_INC_DOWN(NPD_PROFILE)
!             INCIDENT FLUX
     &   , SOURCE_GROUND_FREE(NPD_PROFILE)
!             SOURCE FROM GROUND (CLEAR SKY)
     &   , SOURCE_GROUND_CLOUD(NPD_PROFILE)
!             SOURCE FROM GROUND (CLOUDY REGION)
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)
!             DIFFUSE ALBEDO
      REAL      !, INTENT(OUT)
     &     FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)
!             TOTAL FLUX
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
!     EFFECTIVE COUPLING ALBEDOS AND SOURCE FUNCTIONS:
      REAL
     &     ALPHA11(NPD_PROFILE, NPD_LAYER+1)
     &   , ALPHA22(NPD_PROFILE, NPD_LAYER+1)
     &   , ALPHA21(NPD_PROFILE, NPD_LAYER+1)
     &   , ALPHA12(NPD_PROFILE, NPD_LAYER+1)
     &   , G1(NPD_PROFILE, NPD_LAYER+1)
     &   , G2(NPD_PROFILE, NPD_LAYER+1)
!     TERMS FOR DOWNWARD PROPAGATION:
      REAL
     &     GAMMA11(NPD_PROFILE, NPD_LAYER)
     &   , GAMMA12(NPD_PROFILE, NPD_LAYER)
     &   , GAMMA21(NPD_PROFILE, NPD_LAYER)
     &   , GAMMA22(NPD_PROFILE, NPD_LAYER)
     &   , BETA11_INV(NPD_PROFILE, NPD_LAYER)
     &   , BETA21(NPD_PROFILE, NPD_LAYER)
     &   , BETA22_INV(NPD_PROFILE, NPD_LAYER)
     &   , H1(NPD_PROFILE, NPD_LAYER)
     &   , H2(NPD_PROFILE, NPD_LAYER)
!
!     AUXILAIRY NUMERICAL VARIABLES REQUIRED ONLY IN THE CURRENT LAYER:
      REAL
     &     THETA11
     &   , THETA12
     &   , THETA21
     &   , THETA22
     &   , LAMBDA
     &   , LAMBDA_BAR
!
!     TEMPORARY FLUXES
      REAL
     &     FLUX_DOWN_1(NPD_PROFILE, 0: NPD_LAYER)
!             DOWNWARD FLUXES OUTSIDE CLOUDS JUST BELOW I'TH LEVEL
     &   , FLUX_DOWN_2(NPD_PROFILE, 0: NPD_LAYER)
!             DOWNWARD FLUXES INSIDE CLOUDS JUST BELOW I'TH LEVEL
     &   , FLUX_UP_1(NPD_PROFILE, 0: NPD_LAYER)
!             UPWARD FLUXES OUTSIDE CLOUDS JUST ABOVE I'TH LEVEL
     &   , FLUX_UP_2(NPD_PROFILE, 0: NPD_LAYER)
!             UPWARD FLUXES INSIDE CLOUDS JUST ABOVE I'TH LEVEL
!
!
!
!     INITIALIZE AT THE BOTTOM OF THE COLUMN FOR UPWARD ELIMINATION.
      DO L=1, N_PROFILE
         ALPHA11(L, N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         ALPHA22(L, N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         ALPHA21(L, N_LAYER+1)=0.0E+00
         ALPHA12(L, N_LAYER+1)=0.0E+00
         G1(L, N_LAYER+1)=SOURCE_GROUND_FREE(L)
         G2(L, N_LAYER+1)=SOURCE_GROUND_CLOUD(L)
      ENDDO
!
!     UPWARD ELIMINATION THROUGH THE CLOUDY LAYERS.
      DO I=N_LAYER, N_CLOUD_TOP, -1
         DO L=1, N_PROFILE
!
            THETA11=ALPHA11(L, I+1)*V11(L, I)+ALPHA12(L, I+1)*V21(L, I)
            THETA12=ALPHA11(L, I+1)*V12(L, I)+ALPHA12(L, I+1)*V22(L, I)
            THETA21=ALPHA21(L, I+1)*V11(L, I)+ALPHA22(L, I+1)*V21(L, I)
            THETA22=ALPHA21(L, I+1)*V12(L, I)+ALPHA22(L, I+1)*V22(L, I)

            BETA21(L, I)=-THETA21*R(L, I)
            BETA22_INV(L, I)=1.0E+00/(1.0E+00-THETA22*R_CLOUD(L, I))
            GAMMA21(L, I)=THETA21*T(L, I)
            GAMMA22(L, I)=THETA22*T_CLOUD(L, I)
            H2(L, I)=G2(L, I+1)+THETA21*S_DOWN(L, I)
     &         +THETA22*S_DOWN_CLOUD(L, I)

            LAMBDA=THETA12*R_CLOUD(L, I)*BETA22_INV(L, I)
            BETA11_INV(L, I)=1.0E+00
     &         /(1.0E+00-THETA11*R(L, I)+LAMBDA*BETA21(L, I))
            GAMMA11(L, I)=THETA11*T(L, I)+LAMBDA*GAMMA21(L, I)
            GAMMA12(L, I)=THETA12*T_CLOUD(L, I)+LAMBDA*GAMMA22(L, I)
            H1(L, I)=G1(L, I+1)+THETA11*S_DOWN(L, I)
     &         +THETA12*S_DOWN_CLOUD(L, I)+LAMBDA*H2(L, I)

            LAMBDA=U22(L, I-1)*T_CLOUD(L, I)*BETA22_INV(L, I)
            LAMBDA_BAR=(U21(L, I-1)*T(L, I)+LAMBDA*BETA21(L, I))
     &         *BETA11_INV(L, I)
            ALPHA21(L, I)=U21(L, I-1)*R(L, I)+LAMBDA*GAMMA21(L, I)
     &         +LAMBDA_BAR*GAMMA11(L, I)
            ALPHA22(L, I)=U22(L, I-1)*R_CLOUD(L, I)
     &         +LAMBDA*GAMMA22(L, I)+LAMBDA_BAR*GAMMA12(L, I)
            G2(L, I)=U21(L, I-1)*S_UP(L, I)+U22(L, I-1)*S_UP_CLOUD(L, I)
     &         +LAMBDA*H2(L, I)+LAMBDA_BAR*H1(L, I)
!
            LAMBDA=U12(L, I-1)*T_CLOUD(L, I)*BETA22_INV(L, I)
            LAMBDA_BAR=(U11(L, I-1)*T(L, I)+LAMBDA*BETA21(L, I))
     &         *BETA11_INV(L, I)
            ALPHA11(L, I)=U11(L, I-1)*R(L, I)+LAMBDA*GAMMA21(L, I)
     &         +LAMBDA_BAR*GAMMA11(L, I)
            ALPHA12(L, I)=U12(L, I-1)*R_CLOUD(L, I)
     &         +LAMBDA*GAMMA22(L, I)+LAMBDA_BAR*GAMMA12(L, I)
            G1(L, I)=U11(L, I-1)*S_UP(L, I)+U12(L, I-1)*S_UP_CLOUD(L, I)
     &         +LAMBDA*H2(L, I)+LAMBDA_BAR*H1(L, I)
!
         ENDDO
      ENDDO
!
!     THE LAYER ABOVE THE CLOUD: ONLY ONE SET OF ALPHAS IS NOW NEEDED.
!
      I=N_CLOUD_TOP-1
      DO L=1, N_PROFILE
!
         IF (N_CLOUD_TOP.LT.N_LAYER) THEN
!           IF THERE IS NO CLOUD IN THE COLUMN THE V'S WILL NOT BE
!           ASSIGNED SO AN IF TEST IS REQUIRED.
            THETA11=ALPHA11(L, I+1)*V11(L, I)+ALPHA12(L, I+1)*V21(L, I)
         ELSE
            THETA11=ALPHA11(L, I+1)
         ENDIF
!
         BETA11_INV(L, I)=1.0E+00/(1.0E+00-THETA11*R(L, I))
         GAMMA11(L, I)=THETA11*T(L, I)
         H1(L, I)=G1(L, I+1)+THETA11*S_DOWN(L, I)
!
         LAMBDA=T(L, I)*BETA11_INV(L, I)
         ALPHA11(L, I)=R(L, I)+LAMBDA*GAMMA11(L, I)
         G1(L, I)=S_UP(L, I)+LAMBDA*H1(L, I)
!
      ENDDO
!
      DO I=N_CLOUD_TOP-2, 1, -1
         DO L=1, N_PROFILE
!
            BETA11_INV(L, I)=1.0E+00/(1.0E+00-ALPHA11(L, I+1)*R(L, I))
            GAMMA11(L, I)=ALPHA11(L, I+1)*T(L, I)
            H1(L, I)=G1(L, I+1)+ALPHA11(L, I+1)*S_DOWN(L, I)
!
            LAMBDA=T(L, I)*BETA11_INV(L, I)
            ALPHA11(L, I)=R(L, I)+LAMBDA*GAMMA11(L, I)
            G1(L, I)=S_UP(L, I)+LAMBDA*H1(L, I)
!
         ENDDO
      ENDDO
!
!
!     INITIALIZE FOR DOWNWARD BACK-SUBSTITUTION.
      DO L=1, N_PROFILE
         FLUX_TOTAL(L, 2)=FLUX_INC_DOWN(L)
         FLUX_TOTAL(L, 1)=ALPHA11(L, 1)*FLUX_TOTAL(L, 2)+G1(L, 1)
      ENDDO
!
!     SWEEP DOWNWARD THROUGH THE CLEAR-SKY REGION, FINDING THE DOWNWARD
!     FLUX AT THE TOP OF THE LAYER AND THE UPWARD FLUX AT THE BOTTOM.
      DO I=1, N_CLOUD_TOP-1
         DO L=1, N_PROFILE
            FLUX_TOTAL(L, 2*I+1)=(GAMMA11(L, I)*FLUX_TOTAL(L, 2*I)
     &         +H1(L, I))*BETA11_INV(L, I)
            FLUX_TOTAL(L, 2*I+2)=T(L, I)*FLUX_TOTAL(L, 2*I)
     &         +R(L, I)*FLUX_TOTAL(L, 2*I+1)+S_DOWN(L, I)
         ENDDO
      ENDDO
!
!     PASS INTO THE TOP CLOUDY LAYER. USE FLUX_DOWN_[1,2] TO HOLD,
!     PROVISIONALLY, THE DOWNWARD FLUXES JUST BELOW THE TOP OF THE
!     LAYER, THEN CALCULATE THE UPWARD FLUXES AT THE BOTTOM AND
!     FINALLY THE DOWNWARD FLUXES AT THE BOTTOM OF THE LAYER.
      I=N_CLOUD_TOP
      DO L=1, N_PROFILE
         FLUX_DOWN_1(L, I)=V11(L, I-1)*FLUX_TOTAL(L, 2*I)
         FLUX_DOWN_2(L, I)=V21(L, I-1)*FLUX_TOTAL(L, 2*I)
         FLUX_UP_1(L, I)=(GAMMA11(L, I)*FLUX_DOWN_1(L, I)
     &      +GAMMA12(L, I)*FLUX_DOWN_2(L, I)+H1(L, I))*BETA11_INV(L, I)
         FLUX_UP_2(L, I)=(GAMMA21(L, I)*FLUX_DOWN_1(L, I)
     &      +GAMMA22(L, I)*FLUX_DOWN_2(L, I)+H2(L, I)
     &      -BETA21(L, I)*FLUX_UP_1(L, I))*BETA22_INV(L, I)
         FLUX_DOWN_1(L, I)=T(L, I)*FLUX_DOWN_1(L, I)
     &      +R(L, I)*FLUX_UP_1(L, I)+S_DOWN(L, I)
         FLUX_DOWN_2(L, I)=T_CLOUD(L, I)*FLUX_DOWN_2(L, I)
     &      +R_CLOUD(L, I)*FLUX_UP_2(L, I)+S_DOWN_CLOUD(L, I)
      ENDDO
!
!     THE MAIN LOOP OF BACK-SUBSTITUTION. THE PROVISIONAL USE OF THE
!     DOWNWARD FLUXES IS AS ABOVE.
      DO I=N_CLOUD_TOP+1, N_LAYER
         DO L=1, N_PROFILE
            FLUX_DOWN_1(L, I)=V11(L, I-1)*FLUX_DOWN_1(L, I-1)
     &         +V12(L, I-1)*FLUX_DOWN_2(L, I-1)
            FLUX_DOWN_2(L, I)=V21(L, I-1)*FLUX_DOWN_1(L, I-1)
     &         +V22(L, I-1)*FLUX_DOWN_2(L, I-1)
            FLUX_UP_1(L, I)=(GAMMA11(L, I)*FLUX_DOWN_1(L, I)
     &         +GAMMA12(L, I)*FLUX_DOWN_2(L, I)+H1(L, I))
     &         *BETA11_INV(L, I)
            FLUX_UP_2(L, I)=(GAMMA21(L, I)*FLUX_DOWN_1(L, I)
     &         +GAMMA22(L, I)*FLUX_DOWN_2(L, I)
     &         -BETA21(L, I)*FLUX_UP_1(L, I)+H2(L, I))
     &         *BETA22_INV(L, I)
            FLUX_DOWN_1(L, I)=T(L, I)*FLUX_DOWN_1(L, I)
     &         +R(L, I)*FLUX_UP_1(L, I)+S_DOWN(L, I)
            FLUX_DOWN_2(L, I)=T_CLOUD(L, I)*FLUX_DOWN_2(L, I)
     &         +R_CLOUD(L, I)*FLUX_UP_2(L, I)+S_DOWN_CLOUD(L, I)
         ENDDO
      ENDDO
!
!
!     CALCULATE THE OVERALL FLUX.
      DO I=N_CLOUD_TOP, N_LAYER
         DO L=1, N_PROFILE
            FLUX_TOTAL(L, 2*I+1)=FLUX_UP_1(L, I)+FLUX_UP_2(L, I)
            FLUX_TOTAL(L, 2*I+2)=FLUX_DOWN_1(L, I)+FLUX_DOWN_2(L, I)
         ENDDO
      ENDDO
!
!     REDUCE TO NET FLUXES IF REQUIRED.
      IF (L_NET) THEN
         DO I=0, N_LAYER
            DO L=1, N_PROFILE
               FLUX_TOTAL(L, I+1)
     &            =FLUX_TOTAL(L, 2*I+2)-FLUX_TOTAL(L, 2*I+1)
            ENDDO
         ENDDO
      ENDIF
!
!
!
      RETURN
      END
