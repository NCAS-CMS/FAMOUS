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
!       4.2             10-04-96                Original Code
!                                               (J. M. Edwards)
!  4.5  12/06/98  Fujitsu specific mod to split large loop to
!                 get longer vector length. R.Barnes@ecmwf.int 
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLVER_TRIPLE(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &   , T, R, S_DOWN, S_UP
     &   , T_STRAT, R_STRAT, S_DOWN_STRAT, S_UP_STRAT
     &   , T_CONV, R_CONV, S_DOWN_CONV, S_UP_CONV
     &   , V11, V12, V13, V21, V22, V23, V31, V32, V33
     &   , U11, U12, U13, U21, U22, U23, U31, U32, U33
     &   , L_NET
     &   , FLUX_INC_DOWN
     &   , SOURCE_GROUND_FREE, SOURCE_GROUND_STRAT
     &   , SOURCE_GROUND_CONV, ALBEDO_SURFACE_DIFF
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
     &   , T_STRAT(NPD_PROFILE, NPD_LAYER)
!             STRATFIFORM TRANSMISSION
     &   , R_STRAT(NPD_PROFILE, NPD_LAYER)
!             STRATFIFORM REFLECTION
     &   , S_DOWN_STRAT(NPD_PROFILE, NPD_LAYER)
!             DOWNWARD STRATFIFORM SOURCE FUNCTION
     &   , S_UP_STRAT(NPD_PROFILE, NPD_LAYER)
!             UPWARD STRATFIFORM SOURCE FUNCTION
     &   , T_CONV(NPD_PROFILE, NPD_LAYER)
!             CONVECTIVE TRANSMISSION
     &   , R_CONV(NPD_PROFILE, NPD_LAYER)
!             CONVECTIVE REFLECTION
     &   , S_DOWN_CONV(NPD_PROFILE, NPD_LAYER)
!             DOWNWARD CONVECTIVE SOURCE FUNCTION
     &   , S_UP_CONV(NPD_PROFILE, NPD_LAYER)
!             UPWARD CONVECTIVE SOURCE FUNCTION
      REAL      !, INTENT(IN)
     &     V11(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
     &   , V12(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
     &   , V13(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
     &   , V21(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
     &   , V22(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
     &   , V23(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
     &   , V31(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
     &   , V32(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
     &   , V33(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
      REAL
     &     U11(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
     &   , U12(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
     &   , U13(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
     &   , U21(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
     &   , U22(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
     &   , U23(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
     &   , U31(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
     &   , U32(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
     &   , U33(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
      REAL      !, INTENT(IN)
     &     FLUX_INC_DOWN(NPD_PROFILE)
!             INCIDENT FLUX
     &   , SOURCE_GROUND_FREE(NPD_PROFILE)
!             SOURCE FROM GROUND (CLEAR SKY)
     &   , SOURCE_GROUND_STRAT(NPD_PROFILE)
!             SOURCE FROM GROUND (CLOUDY REGION)
     &   , SOURCE_GROUND_CONV(NPD_PROFILE)
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
     &   , ALPHA12(NPD_PROFILE, NPD_LAYER+1)
     &   , ALPHA13(NPD_PROFILE, NPD_LAYER+1)
     &   , ALPHA21(NPD_PROFILE, NPD_LAYER+1)
     &   , ALPHA22(NPD_PROFILE, NPD_LAYER+1)
     &   , ALPHA23(NPD_PROFILE, NPD_LAYER+1)
     &   , ALPHA31(NPD_PROFILE, NPD_LAYER+1)
     &   , ALPHA32(NPD_PROFILE, NPD_LAYER+1)
     &   , ALPHA33(NPD_PROFILE, NPD_LAYER+1)
     &   , G1(NPD_PROFILE, NPD_LAYER+1)
     &   , G2(NPD_PROFILE, NPD_LAYER+1)
     &   , G3(NPD_PROFILE, NPD_LAYER+1)
!     TERMS FOR DOWNWARD PROPAGATION:
      REAL
     &     GAMMA11(NPD_PROFILE, NPD_LAYER)
     &   , GAMMA12(NPD_PROFILE, NPD_LAYER)
     &   , GAMMA13(NPD_PROFILE, NPD_LAYER)
     &   , GAMMA21(NPD_PROFILE, NPD_LAYER)
     &   , GAMMA22(NPD_PROFILE, NPD_LAYER)
     &   , GAMMA23(NPD_PROFILE, NPD_LAYER)
     &   , GAMMA31(NPD_PROFILE, NPD_LAYER)
     &   , GAMMA32(NPD_PROFILE, NPD_LAYER)
     &   , GAMMA33(NPD_PROFILE, NPD_LAYER)
     &   , BETA11_INV(NPD_PROFILE, NPD_LAYER)
     &   , BETA21(NPD_PROFILE, NPD_LAYER)
     &   , BETA22_INV(NPD_PROFILE, NPD_LAYER)
     &   , BETA31(NPD_PROFILE, NPD_LAYER)
     &   , BETA32(NPD_PROFILE, NPD_LAYER)
     &   , BETA33_INV(NPD_PROFILE, NPD_LAYER)
     &   , H1(NPD_PROFILE, NPD_LAYER)
     &   , H2(NPD_PROFILE, NPD_LAYER)
     &   , H3(NPD_PROFILE, NPD_LAYER)
!
!     AUXILAIRY NUMERICAL VARIABLES REQUIRED ONLY IN THE CURRENT LAYER:
      REAL
     &     THETA11
     &   , THETA12
     &   , THETA13
     &   , THETA21
     &   , THETA22
     &   , THETA23
     &   , THETA31
     &   , THETA32
     &   , THETA33
     &   , LAMBDA3
     &   , LAMBDA2
     &   , LAMBDA1
     &   , LAMBDA
!
!     TEMPORARY FLUXES
      REAL
     &     FLUX_DOWN_1(NPD_PROFILE, 0: NPD_LAYER)
!             DOWNWARD FLUXES OUTSIDE CLOUDS JUST BELOW I'TH LEVEL
     &   , FLUX_DOWN_2(NPD_PROFILE, 0: NPD_LAYER)
!             DOWNWARD FLUXES INSIDE CLOUDS JUST BELOW I'TH LEVEL
     &   , FLUX_DOWN_3(NPD_PROFILE, 0: NPD_LAYER)
!             DOWNWARD FLUXES INSIDE CLOUDS JUST BELOW I'TH LEVEL
     &   , FLUX_UP_1(NPD_PROFILE, 0: NPD_LAYER)
!             UPWARD FLUXES OUTSIDE CLOUDS JUST ABOVE I'TH LEVEL
     &   , FLUX_UP_2(NPD_PROFILE, 0: NPD_LAYER)
!             UPWARD FLUXES INSIDE CLOUDS JUST ABOVE I'TH LEVEL
     &   , FLUX_UP_3(NPD_PROFILE, 0: NPD_LAYER)
!             UPWARD FLUXES INSIDE CLOUDS JUST ABOVE I'TH LEVEL
!
!
!
!     THIS ROUTINE IS SPECIFIC TO CASES OF THREE REGIONS AND IT IS
!     ASSUMED THAT 1 REPRESENTS CLEAR SKIES, 2 REPRESENTS STARTIFORM
!     CLOUDS AND 3 REPRESENTS CONVECTIVE CLOUD.
!
!     INITIALIZE AT THE BOTTOM OF THE COLUMN FOR UPWARD ELIMINATION.
      DO L=1, N_PROFILE
         ALPHA11(L, N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         ALPHA12(L, N_LAYER+1)=0.0E+00
         ALPHA13(L, N_LAYER+1)=0.0E+00
         ALPHA21(L, N_LAYER+1)=0.0E+00
         ALPHA22(L, N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         ALPHA23(L, N_LAYER+1)=0.0E+00
         ALPHA31(L, N_LAYER+1)=0.0E+00
         ALPHA32(L, N_LAYER+1)=0.0E+00
         ALPHA33(L, N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         G1(L, N_LAYER+1)=SOURCE_GROUND_FREE(L)
         G2(L, N_LAYER+1)=SOURCE_GROUND_STRAT(L)
         G3(L, N_LAYER+1)=SOURCE_GROUND_CONV(L)
      ENDDO
!
!     UPWARD ELIMINATION THROUGH THE CLOUDY LAYERS.
      DO I=N_LAYER, N_CLOUD_TOP, -1
         DO L=1, N_PROFILE
!
            THETA11=ALPHA11(L, I+1)*V11(L, I)+ALPHA12(L, I+1)*V21(L, I)
     &         +ALPHA13(L, I+1)*V31(L, I)
            THETA12=ALPHA11(L, I+1)*V12(L, I)+ALPHA12(L, I+1)*V22(L, I)
     &         +ALPHA13(L, I+1)*V32(L, I)
            THETA13=ALPHA11(L, I+1)*V13(L, I)+ALPHA12(L, I+1)*V23(L, I)
     &         +ALPHA13(L, I+1)*V33(L, I)
            THETA21=ALPHA21(L, I+1)*V11(L, I)+ALPHA22(L, I+1)*V21(L, I)
     &         +ALPHA23(L, I+1)*V31(L, I)
            THETA22=ALPHA21(L, I+1)*V12(L, I)+ALPHA22(L, I+1)*V22(L, I)
     &         +ALPHA23(L, I+1)*V32(L, I)
            THETA23=ALPHA21(L, I+1)*V13(L, I)+ALPHA22(L, I+1)*V23(L, I)
     &         +ALPHA23(L, I+1)*V33(L, I)
            THETA31=ALPHA31(L, I+1)*V11(L, I)+ALPHA32(L, I+1)*V21(L, I)
     &         +ALPHA33(L, I+1)*V31(L, I)
            THETA32=ALPHA31(L, I+1)*V12(L, I)+ALPHA32(L, I+1)*V22(L, I)
     &         +ALPHA33(L, I+1)*V32(L, I)
            THETA33=ALPHA31(L, I+1)*V13(L, I)+ALPHA32(L, I+1)*V23(L, I)
     &         +ALPHA33(L, I+1)*V33(L, I)

            BETA31(L, I)=-THETA31*R(L, I)
            BETA32(L, I)=-THETA32*R_STRAT(L, I)
            BETA33_INV(L, I)=1.0E+00/(1.0E+00-THETA33*R_CONV(L, I))
            GAMMA31(L, I)=THETA31*T(L, I)
            GAMMA32(L, I)=THETA32*T_STRAT(L, I)
            GAMMA33(L, I)=THETA33*T_CONV(L, I)
            H3(L, I)=G3(L, I+1)+THETA31*S_DOWN(L, I)
     &         +THETA32*S_DOWN_STRAT(L, I)
     &         +THETA33*S_DOWN_CONV(L, I)
!
            LAMBDA3=THETA23*R_CONV(L, I)*BETA33_INV(L, I)
            BETA22_INV(L, I)=1.0E+00
     &         /(1.0E+00-THETA22*R_STRAT(L, I)+LAMBDA3*BETA32(L, I))
            BETA21(L, I)=-THETA21*R(L, I)+LAMBDA3*BETA31(L, I)
            GAMMA21(L, I)=THETA21*T(L, I)+LAMBDA3*GAMMA31(L, I)
            GAMMA22(L, I)=THETA22*T_STRAT(L, I)+LAMBDA3*GAMMA32(L, I)
            GAMMA23(L, I)=THETA23*T_CONV(L, I)+LAMBDA3*GAMMA33(L, I)
            H2(L, I)=G2(L, I+1)+THETA21*S_DOWN(L, I)
     &         +THETA22*S_DOWN_STRAT(L, I)+THETA23*S_DOWN_CONV(L, I)
     &         +LAMBDA3*H3(L, I)
!
            LAMBDA3=THETA13*R_CONV(L, I)*BETA33_INV(L, I)
            LAMBDA2=(THETA12*R_STRAT(L, I)-LAMBDA3*BETA32(L, I))
     &         *BETA22_INV(L, I)
            BETA11_INV(L, I)=1.0E+00
     &         /(1.0E+00-THETA11*R(L, I)+LAMBDA3*BETA31(L, I)
     &         +LAMBDA2*BETA21(L, I))
            GAMMA11(L, I)=THETA11*T(L, I)+LAMBDA3*GAMMA31(L, I)
     &         +LAMBDA2*GAMMA21(L, I)
            GAMMA12(L, I)=THETA12*T_STRAT(L, I)+LAMBDA3*GAMMA32(L, I)
     &         +LAMBDA2*GAMMA22(L, I)
            GAMMA13(L, I)=THETA13*T_CONV(L, I)+LAMBDA3*GAMMA33(L, I)
     &         +LAMBDA2*GAMMA23(L, I)
            H1(L, I)=G1(L, I+1)+THETA11*S_DOWN(L, I)
     &         +THETA12*S_DOWN_STRAT(L, I)+THETA13*S_DOWN_CONV(L, I)
     &         +LAMBDA3*H3(L, I)+LAMBDA2*H2(L, I)
            LAMBDA3=U33(L, I-1)*T_CONV(L, I)*BETA33_INV(L, I)
            LAMBDA2=(U32(L, I-1)*T_STRAT(L, I)+LAMBDA3*BETA32(L, I))
     &         *BETA22_INV(L, I)
            LAMBDA1=(U31(L, I-1)*T(L, I)+LAMBDA3*BETA31(L, I)
     &         +LAMBDA2*BETA21(L, I))*BETA11_INV(L, I)
            ALPHA31(L, I)=U31(L, I-1)*R(L, I)+LAMBDA3*GAMMA31(L, I)
     &         +LAMBDA2*GAMMA21(L, I)+LAMBDA1*GAMMA11(L, I)
            ALPHA32(L, I)=U32(L, I-1)*R_STRAT(L, I)
     &         +LAMBDA3*GAMMA32(L, I)+LAMBDA2*GAMMA22(L, I)
     &         +LAMBDA1*GAMMA12(L, I)
            ALPHA33(L, I)=U33(L, I-1)*R_CONV(L, I)
     &         +LAMBDA3*GAMMA33(L, I)+LAMBDA2*GAMMA23(L, I)
     &         +LAMBDA1*GAMMA13(L, I)
            G3(L, I)=U31(L, I-1)*S_UP(L, I)+U32(L, I-1)*S_UP_STRAT(L, I)
     &         +U33(L, I-1)*S_UP_CONV(L, I)
     &         +LAMBDA3*H3(L, I)+LAMBDA2*H2(L, I)+LAMBDA1*H1(L, I)
!
            LAMBDA3=U23(L, I-1)*T_CONV(L, I)*BETA33_INV(L, I)
            LAMBDA2=(U22(L, I-1)*T_STRAT(L, I)+LAMBDA3*BETA32(L, I))
     &         *BETA22_INV(L, I)
            LAMBDA1=(U21(L, I-1)*T(L, I)+LAMBDA3*BETA31(L, I)
     &         +LAMBDA2*BETA21(L, I))*BETA11_INV(L, I)
            ALPHA21(L, I)=U21(L, I-1)*R(L, I)+LAMBDA3*GAMMA31(L, I)
     &         +LAMBDA2*GAMMA21(L, I)+LAMBDA1*GAMMA11(L, I)
            ALPHA22(L, I)=U22(L, I-1)*R_STRAT(L, I)
     &         +LAMBDA3*GAMMA32(L, I)+LAMBDA2*GAMMA22(L, I)
     &         +LAMBDA1*GAMMA12(L, I)
            ALPHA23(L, I)=U23(L, I-1)*R_CONV(L, I)
     &         +LAMBDA3*GAMMA33(L, I)+LAMBDA2*GAMMA23(L, I)
     &         +LAMBDA1*GAMMA13(L, I)
            G2(L, I)=U21(L, I-1)*S_UP(L, I)+U22(L, I-1)*S_UP_STRAT(L, I)
     &         +U23(L, I-1)*S_UP_CONV(L, I)
     &         +LAMBDA3*H3(L, I)+LAMBDA2*H2(L, I)+LAMBDA1*H1(L, I)
!
            LAMBDA3=U13(L, I-1)*T_CONV(L, I)*BETA33_INV(L, I)
            LAMBDA2=(U12(L, I-1)*T_STRAT(L, I)+LAMBDA3*BETA32(L, I))
     &         *BETA22_INV(L, I)
            LAMBDA1=(U11(L, I-1)*T(L, I)+LAMBDA3*BETA31(L, I)
     &         +LAMBDA2*BETA21(L, I))*BETA11_INV(L, I)
            ALPHA11(L, I)=U11(L, I-1)*R(L, I)+LAMBDA3*GAMMA31(L, I)
     &         +LAMBDA2*GAMMA21(L, I)+LAMBDA1*GAMMA11(L, I)
            ALPHA12(L, I)=U12(L, I-1)*R_STRAT(L, I)
     &         +LAMBDA3*GAMMA32(L, I)+LAMBDA2*GAMMA22(L, I)
     &         +LAMBDA1*GAMMA12(L, I)
            ALPHA13(L, I)=U13(L, I-1)*R_CONV(L, I)
     &         +LAMBDA3*GAMMA33(L, I)+LAMBDA2*GAMMA23(L, I)
     &         +LAMBDA1*GAMMA13(L, I)
            G1(L, I)=U11(L, I-1)*S_UP(L, I)+U12(L, I-1)*S_UP_STRAT(L, I)
     &         +U13(L, I-1)*S_UP_CONV(L, I)
     &         +LAMBDA3*H3(L, I)+LAMBDA2*H2(L, I)+LAMBDA1*H1(L, I)
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
     &         +ALPHA13(L, I+1)*V31(L, I)
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
            LAMBDA1=T(L, I)*BETA11_INV(L, I)
            ALPHA11(L, I)=R(L, I)+LAMBDA1*GAMMA11(L, I)
            G1(L, I)=S_UP(L, I)+LAMBDA1*H1(L, I)
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
!     PASS INTO THE TOP CLOUDY LAYER. USE FLUX_DOWN_[1,2,3] TO HOLD,
!     PROVISIONALLY, THE DOWNWARD FLUXES JUST BELOW THE TOP OF THE
!     LAYER, THEN CALCULATE THE UPWARD FLUXES AT THE BOTTOM AND
!     FINALLY THE DOWNWARD FLUXES AT THE BOTTOM OF THE LAYER.
      I=N_CLOUD_TOP
      DO L=1, N_PROFILE
         FLUX_DOWN_1(L, I)=V11(L, I-1)*FLUX_TOTAL(L, 2*I)
         FLUX_DOWN_2(L, I)=V21(L, I-1)*FLUX_TOTAL(L, 2*I)
         FLUX_DOWN_3(L, I)=V31(L, I-1)*FLUX_TOTAL(L, 2*I)
         FLUX_UP_1(L, I)=(GAMMA11(L, I)*FLUX_DOWN_1(L, I)
     &      +GAMMA12(L, I)*FLUX_DOWN_2(L, I)
     &      +GAMMA13(L, I)*FLUX_DOWN_3(L, I)
     &      +H1(L, I))*BETA11_INV(L, I)
         FLUX_UP_2(L, I)=(GAMMA21(L, I)*FLUX_DOWN_1(L, I)
     &      +GAMMA22(L, I)*FLUX_DOWN_2(L, I)
     &      +GAMMA23(L, I)*FLUX_DOWN_3(L, I)+H2(L, I)
     &      -BETA21(L, I)*FLUX_UP_1(L, I))*BETA22_INV(L, I)
         FLUX_UP_3(L, I)=(GAMMA31(L, I)*FLUX_DOWN_1(L, I)
     &      +GAMMA32(L, I)*FLUX_DOWN_2(L, I)
     &      +GAMMA33(L, I)*FLUX_DOWN_3(L, I)+H3(L, I)
     &      -BETA31(L, I)*FLUX_UP_1(L, I)-BETA32(L, I)*FLUX_UP_2(L, I))
     &      *BETA33_INV(L, I)
         FLUX_DOWN_1(L, I)=T(L, I)*FLUX_DOWN_1(L, I)
     &      +R(L, I)*FLUX_UP_1(L, I)+S_DOWN(L, I)
         FLUX_DOWN_2(L, I)=T_STRAT(L, I)*FLUX_DOWN_2(L, I)
     &      +R_STRAT(L, I)*FLUX_UP_2(L, I)+S_DOWN_STRAT(L, I)
         FLUX_DOWN_3(L, I)=T_CONV(L, I)*FLUX_DOWN_3(L, I)
     &      +R_CONV(L, I)*FLUX_UP_3(L, I)+S_DOWN_CONV(L, I)
      ENDDO
!
!     THE MAIN LOOP OF BACK-SUBSTITUTION. THE PROVISIONAL USE OF THE
!     DOWNWARD FLUXES IS AS ABOVE.
      DO I=N_CLOUD_TOP+1, N_LAYER
         DO L=1, N_PROFILE
            FLUX_DOWN_1(L, I)=V11(L, I-1)*FLUX_DOWN_1(L, I-1)
     &         +V12(L, I-1)*FLUX_DOWN_2(L, I-1)
     &         +V13(L, I-1)*FLUX_DOWN_3(L, I-1)
            FLUX_DOWN_2(L, I)=V21(L, I-1)*FLUX_DOWN_1(L, I-1)
     &         +V22(L, I-1)*FLUX_DOWN_2(L, I-1)
     &         +V23(L, I-1)*FLUX_DOWN_3(L, I-1)
            FLUX_DOWN_3(L, I)=V31(L, I-1)*FLUX_DOWN_1(L, I-1)
     &         +V32(L, I-1)*FLUX_DOWN_2(L, I-1)
     &         +V33(L, I-1)*FLUX_DOWN_3(L, I-1)
            FLUX_UP_1(L, I)=(GAMMA11(L, I)*FLUX_DOWN_1(L, I)
     &         +GAMMA12(L, I)*FLUX_DOWN_2(L, I)
     &         +GAMMA13(L, I)*FLUX_DOWN_3(L, I)
     &         +H1(L, I))*BETA11_INV(L, I)
            FLUX_UP_2(L, I)=(GAMMA21(L, I)*FLUX_DOWN_1(L, I)
     &         +GAMMA22(L, I)*FLUX_DOWN_2(L, I)
     &         +GAMMA23(L, I)*FLUX_DOWN_3(L, I)+H2(L, I)
     &         -BETA21(L, I)*FLUX_UP_1(L, I))*BETA22_INV(L, I)
            FLUX_UP_3(L, I)=(GAMMA31(L, I)*FLUX_DOWN_1(L, I)
     &         +GAMMA32(L, I)*FLUX_DOWN_2(L, I)
     &         +GAMMA33(L, I)*FLUX_DOWN_3(L, I)+H3(L, I)
     &         -BETA31(L, I)*FLUX_UP_1(L, I)
     &         -BETA32(L, I)*FLUX_UP_2(L, I))
     &         *BETA33_INV(L, I)
            FLUX_DOWN_1(L, I)=T(L, I)*FLUX_DOWN_1(L, I)
     &         +R(L, I)*FLUX_UP_1(L, I)+S_DOWN(L, I)
            FLUX_DOWN_2(L, I)=T_STRAT(L, I)*FLUX_DOWN_2(L, I)
     &         +R_STRAT(L, I)*FLUX_UP_2(L, I)+S_DOWN_STRAT(L, I)
            FLUX_DOWN_3(L, I)=T_CONV(L, I)*FLUX_DOWN_3(L, I)
     &         +R_CONV(L, I)*FLUX_UP_3(L, I)+S_DOWN_CONV(L, I)
         ENDDO
      ENDDO
!
!
!     CALCULATE THE OVERALL FLUX.
      DO I=N_CLOUD_TOP, N_LAYER
         DO L=1, N_PROFILE
            FLUX_TOTAL(L, 2*I+1)=FLUX_UP_1(L, I)+FLUX_UP_2(L, I)
     &         +FLUX_UP_3(L, I)
            FLUX_TOTAL(L, 2*I+2)=FLUX_DOWN_1(L, I)+FLUX_DOWN_2(L, I)
     &         +FLUX_DOWN_3(L, I)
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
