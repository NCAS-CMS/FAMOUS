C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1998, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!       4.5             11-06-98                Optimised version
!                                               (P. Burton)
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
     &     ALPHA(3,3,NPD_LAYER+1),G(3,NPD_LAYER+1)

!     TERMS FOR DOWNWARD PROPAGATION:
      REAL
     &     GAMMA(3,3,NPD_LAYER), H(3,NPD_LAYER)
     &   , BETA(3,3,NPD_LAYER)

!     AUXILAIRY NUMERICAL VARIABLES REQUIRED ONLY IN THE CURRENT LAYER:
      REAL
     &     THETA(3,3,0:NPD_LAYER)
     &   , LAMBDA(0:3,0:NPD_LAYER)
!
!     TEMPORARY FLUXES
      REAL
     &     FLUX_DOWN(3,0: NPD_LAYER)
     &   , FLUX_UP(3,0: NPD_LAYER)
     &   , FLUX_TEMP(6,0:NPD_LAYER)

!     New temporary arrays for the optimised version
      REAL
     &     UV(18,0:NPD_LAYER,NPD_PROFILE)
     &   , T_R(12,NPD_LAYER,NPD_PROFILE)


!     THIS ROUTINE IS SPECIFIC TO CASES OF THREE REGIONS AND IT IS
!     ASSUMED THAT 1 REPRESENTS CLEAR SKIES, 2 REPRESENTS STARTIFORM
!     CLOUDS AND 3 REPRESENTS CONVECTIVE CLOUD.

! Write U and V into a scalar array

      DO I = 0,N_LAYER
!*DIR$ CACHE_BYPASS UV,V11,V21,V31,V12
         DO L=1,N_PROFILE
            UV(1,I,L) = V11(L,I)
            UV(2,I,L) = V21(L,I)
            UV(3,I,L) = V31(L,I)
            UV(4,I,L) = V12(L,I)
         END DO
      END DO
      DO I=0,N_LAYER
!*DIR$ CACHE_BYPASS UV,V22,V32,V13,V23
         DO L=1,N_PROFILE
            UV(5,I,L) = V22(L,I)
            UV(6,I,L) = V32(L,I)
            UV(7,I,L) = V13(L,I)
            UV(8,I,L) = V23(L,I)
         END DO
      END DO
      DO I=0,N_LAYER
!*DIR$ CACHE_BYPASS UV,U11,U21,U31,V33
         DO L=1,N_PROFILE
            UV(9,I,L) = V33(L,I)
            UV(10,I,L) = U11(L,I)
            UV(11,I,L) = U21(L,I)
            UV(12,I,L) = U31(L,I)
         END DO
      END DO
      DO I=0,N_LAYER
!*DIR$ CACHE_BYPASS UV,U12,U22,U32,U13
         DO L=1,N_PROFILE
            UV(13,I,L) = U12(L,I)
            UV(14,I,L) = U22(L,I)
            UV(15,I,L) = U32(L,I)
            UV(16,I,L) = U13(L,I)
         END DO
      END DO
      DO I=0,N_LAYER
!*DIR$ CACHE_BYPASS UV,U23,U33
         DO L=1,N_PROFILE
            UV(17,I,L) = U23(L,I)
            UV(18,I,L) = U33(L,I)
         END DO
      END DO
      DO I=1,N_LAYER
!*DIR$ CACHE_BYPASS T_R,T,T_STRAT,T_CONV,R
         DO L=1,N_PROFILE
            T_R(1,I,L) = T(L,I)
            T_R(2,I,L) = T_STRAT(L,I)
            T_R(3,I,L) = T_CONV(L,I)
            T_R(4,I,L) = R(L,I)
         END DO
      END DO
      DO I=1,N_LAYER
!*DIR$ CACHE_BYPASS T_R,R_STRAT,R_CONV,S_DOWN,S_UP
         DO L=1,N_PROFILE
            T_R(5,I,L) = R_STRAT(L,I)
            T_R(6,I,L) = R_CONV(L,I)
            T_R(7,I,L) = S_DOWN(L,I)
            T_R(8,I,L) = S_UP(L,I)
         END DO
      END DO
      DO I=1,N_LAYER
!*DIR$ CACHE_BYPASS T_R,S_DOWN_STRAT,S_UP_STRAT,S_UP_CONV,S_DOWN_CONV
         DO L=1,N_PROFILE
            T_R(9,I,L) = S_DOWN_STRAT(L,I)
            T_R(10,I,L) = S_UP_STRAT(L,I)
            T_R(11,I,L) = S_DOWN_CONV(L,I)
            T_R(12,I,L) = S_UP_CONV(L,I)
         END DO
      END DO

!
!     INITIALIZE AT THE BOTTOM OF THE COLUMN FOR UPWARD ELIMINATION.
      DO L=1, N_PROFILE
         ALPHA(1,1,N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         ALPHA(1,2,N_LAYER+1)=0.0E+00
         ALPHA(1,3,N_LAYER+1)=0.0E+00
         ALPHA(2,1,N_LAYER+1)=0.0E+00
         ALPHA(2,2,N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         ALPHA(2,3,N_LAYER+1)=0.0E+00
         ALPHA(3,1,N_LAYER+1)=0.0E+00
         ALPHA(3,2,N_LAYER+1)=0.0E+00
         ALPHA(3,3,N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         G(1,N_LAYER+1)=SOURCE_GROUND_FREE(L)
         G(2,N_LAYER+1)=SOURCE_GROUND_STRAT(L)
         G(3,N_LAYER+1)=SOURCE_GROUND_CONV(L)
!
!     UPWARD ELIMINATION THROUGH THE CLOUDY LAYERS.


         DO I=N_LAYER, N_CLOUD_TOP, -1
            THETA(1,1,I)=ALPHA(1,1,I+1)*UV(1,I,L)
     &          +ALPHA(1,2,I+1)*UV(2,I,L)
     &          +ALPHA(1,3,I+1)*UV(3,I,L)
            THETA(1,2,I)=ALPHA(1,1,I+1)*UV(4,I,L)
     &           +ALPHA(1,2,I+1)*UV(5,I,L)
     &           +ALPHA(1,3,I+1)*UV(6,I,L)
            THETA(1,3,I)=ALPHA(1,1,I+1)*UV(7,I,L)
     &           +ALPHA(1,2,I+1)*UV(8,I,L)
     &           +ALPHA(1,3,I+1)*UV(9,I,L)
            THETA(2,1,I)=ALPHA(2,1,I+1)*UV(1,I,L)
     &            +ALPHA(2,2,I+1)*UV(2,I,L)
     &            +ALPHA(2,3,I+1)*UV(3,I,L)
            THETA(2,2,I)=ALPHA(2,1,I+1)*UV(4,I,L)
     &           +ALPHA(2,2,I+1)*UV(5,I,L)
     &           +ALPHA(2,3,I+1)*UV(6,I,L)
            THETA(2,3,I)=ALPHA(2,1,I+1)*UV(7,I,L)
     &           +ALPHA(2,2,I+1)*UV(8,I,L)
     &           +ALPHA(2,3,I+1)*UV(9,I,L)
            THETA(3,1,I)=ALPHA(3,1,I+1)*UV(1,I,L)
     &           +ALPHA(3,2,I+1)*UV(2,I,L)
     &           +ALPHA(3,3,I+1)*UV(3,I,L)
            THETA(3,2,I)=ALPHA(3,1,I+1)*UV(4,I,L)
     &           +ALPHA(3,2,I+1)*UV(5,I,L)
     &           +ALPHA(3,3,I+1)*UV(6,I,L)
            THETA(3,3,I)=ALPHA(3,1,I+1)*UV(7,I,L)
     &           +ALPHA(3,2,I+1)*UV(8,I,L)
     &           +ALPHA(3,3,I+1)*UV(9,I,L)

            BETA(3,1,I)=-THETA(3,1,I)*T_R(4,I,L)
            BETA(3,2,I)=-THETA(3,2,I)*T_R(5,I,L)
            BETA(3,3,I)=1.0E+00/(1.0E+00-THETA(3,3,I)*T_R(6,I,L))
            GAMMA(3,1,I)=THETA(3,1,I)*T_R(1,I,L)
            GAMMA(3,2,I)=THETA(3,2,I)*T_R(2,I,L)
            GAMMA(3,3,I)=THETA(3,3,I)*T_R(3,I,L)
            H(3,I)=G(3,I+1)+THETA(3,1,I)*T_R(7,I,L)
     &         +THETA(3,2,I)*T_R(9,I,L)
     &         +THETA(3,3,I)*T_R(11,I,L)
!
            LAMBDA(3,I)=THETA(2,3,I)*T_R(6,I,L)*BETA(3,3,I)
            BETA(2,2,I)=1.0E+00
     &         /(1.0E+00-THETA(2,2,I)*T_R(5,I,L)
     &           +LAMBDA(3,I)*BETA(3,2,I))
            BETA(2,1,I)=-THETA(2,1,I)*T_R(4,I,L)+LAMBDA(3,I)*BETA(3,1,I)
            GAMMA(2,1,I)=THETA(2,1,I)*T_R(1,I,L)
     &           +LAMBDA(3,I)*GAMMA(3,1,I)
            GAMMA(2,2,I)=THETA(2,2,I)*T_R(2,I,L)+
     &           LAMBDA(3,I)*GAMMA(3,2,I)
            GAMMA(2,3,I)=THETA(2,3,I)*T_R(3,I,L)+LAMBDA(3,I)*
     &           GAMMA(3,3,I)
            H(2,I)=G(2,I+1)+THETA(2,1,I)*T_R(7,I,L)
     &         +THETA(2,2,I)*T_R(9,I,L)
     &         +THETA(2,3,I)*T_R(11,I,L)
     &         +LAMBDA(3,I)*H(3,I)
!
            LAMBDA(3,I)=THETA(1,3,I)*T_R(6,I,L)*BETA(3,3,I)
            LAMBDA(2,I)=(THETA(1,2,I)*T_R(5,I,L)
     &         -LAMBDA(3,I)*BETA(3,2,I))*BETA(2,2,I)
            BETA(1,1,I)=1.0E+00
     &       /(1.0E+00-THETA(1,1,I)*T_R(4,I,L)+LAMBDA(3,I)*BETA(3,1,I)
     &         +LAMBDA(2,I)*BETA(2,1,I))
            GAMMA(1,1,I)=THETA(1,1,I)*T_R(1,I,L)
     &           +LAMBDA(3,I)*GAMMA(3,1,I)
     &          +LAMBDA(2,I)*GAMMA(2,1,I)
            GAMMA(1,2,I)=THETA(1,2,I)*T_R(2,I,L)
     &           +LAMBDA(3,I)*GAMMA(3,2,I)
     &           +LAMBDA(2,I)*GAMMA(2,2,I)
            GAMMA(1,3,I)=THETA(1,3,I)*T_R(3,I,L)
     &           +LAMBDA(3,I)*GAMMA(3,3,I)
     &           +LAMBDA(2,I)*GAMMA(2,3,I)
            H(1,I)=G(1,I+1)+THETA(1,1,I)*T_R(7,I,L)
     &         +THETA(1,2,I)*T_R(9,I,L)
     &         +THETA(1,3,I)*T_R(11,I,L)
     &         +LAMBDA(3,I)*H(3,I)+LAMBDA(2,I)*H(2,I)

!
            LAMBDA(3,I)=UV(18,I-1,L)*T_R(3,I,L)*BETA(3,3,I)
            LAMBDA(2,I)=(UV(15,I-1,L)*T_R(2,I,L)+LAMBDA(3,I)
     &         *BETA(3,2,I))*BETA(2,2,I)
            LAMBDA(1,I)=(UV(12,I-1,L)*T_R(1,I,L)
     &           +LAMBDA(3,I)*BETA(3,1,I)
     &         +LAMBDA(2,I)*BETA(2,1,I))*BETA(1,1,I)
            ALPHA(3,1,I)=UV(12,I-1,L)*T_R(4,I,L)
     &           +LAMBDA(3,I)*GAMMA(3,1,I)
     &           +LAMBDA(2,I)*GAMMA(2,1,I)+LAMBDA(1,I)*GAMMA(1,1,I)
            ALPHA(3,2,I)=UV(15,I-1,L)*T_R(5,I,L)
     &         +LAMBDA(3,I)*GAMMA(3,2,I)+LAMBDA(2,I)*GAMMA(2,2,I)
     &         +LAMBDA(1,I)*GAMMA(1,2,I)
            ALPHA(3,3,I)=UV(18,I-1,L)*T_R(6,I,L)
     &         +LAMBDA(3,I)*GAMMA(3,3,I)+LAMBDA(2,I)*GAMMA(2,3,I)
     &         +LAMBDA(1,I)*GAMMA(1,3,I)

!
              G(3,I)=UV(12,I-1,L)*T_R(8,I,L)
     &         + UV(15,I-1,L)*T_R(10,I,L)
     &         +UV(18,I-1,L)*T_R(12,I,L)
     &         +LAMBDA(3,I)*H(3,I)+LAMBDA(2,I)*H(2,I)+LAMBDA(1,I)*H(1,I)
!
            LAMBDA(3,I)=UV(17,I-1,l)*T_R(3,I,L)*BETA(3,3,I)
            LAMBDA(2,I)=(UV(14,I-1,l)*T_R(2,I,L)
     &           +LAMBDA(3,I)*BETA(3,2,I))*BETA(2,2,I)
            LAMBDA(1,I)=(UV(11,I-1,l)*T_R(1,I,L)
     &           +LAMBDA(3,I)*BETA(3,1,I)
     &           +LAMBDA(2,I)*BETA(2,1,I))*BETA(1,1,I)
            ALPHA(2,1,I)=UV(11,I-1,L)*T_R(4,I,L)
     &           +LAMBDA(3,I)*GAMMA(3,1,I)
     &           +LAMBDA(2,I)*GAMMA(2,1,I)+LAMBDA(1,I)*GAMMA(1,1,I)
            ALPHA(2,2,I)=UV(14,I-1,l)*T_R(5,I,L)
     &           +LAMBDA(3,I)*GAMMA(3,2,I)+LAMBDA(2,I)*GAMMA(2,2,I)
     &           +LAMBDA(1,I)*GAMMA(1,2,I)
            ALPHA(2,3,I)=UV(17,I-1,L)*T_R(6,I,L)
     &         +LAMBDA(3,I)*GAMMA(3,3,I)+LAMBDA(2,I)*GAMMA(2,3,I)
     &         +LAMBDA(1,I)*GAMMA(1,3,I)
            G(2,I)=UV(11,I-1,L)*T_R(8,I,L)
     &           +UV(14,I-1,L)*T_R(10,I,L)
     &           +UV(17,I-1,L)*T_R(12,I,L)
     &           +LAMBDA(3,I)*H(3,I)+LAMBDA(2,I)*H(2,I)
     &           +LAMBDA(1,I)*H(1,I)
!
            LAMBDA(3,I)=UV(16,I-1,L)*T_R(3,I,L)*BETA(3,3,I)
            LAMBDA(2,I)=(UV(13,I-1,L)*T_R(2,I,L)
     &         +LAMBDA(3,I)*BETA(3,2,I))
     &         *BETA(2,2,I)
            LAMBDA(1,I)=(UV(10,I-1,L)*T_R(1,I,L)
     &           +LAMBDA(3,I)*BETA(3,1,I)
     &           +LAMBDA(2,I)*BETA(2,1,I))*BETA(1,1,I)
            ALPHA(1,1,I)=UV(10,I-1,L)*T_R(4,I,L)
     &           +LAMBDA(3,I)*GAMMA(3,1,I)
     &           +LAMBDA(2,I)*GAMMA(2,1,I)+LAMBDA(1,I)*GAMMA(1,1,I)
            ALPHA(1,2,I)=UV(13,I-1,L)*T_R(5,I,L)
     &         +LAMBDA(3,I)*GAMMA(3,2,I)+LAMBDA(2,I)*GAMMA(2,2,I)
     &         +LAMBDA(1,I)*GAMMA(1,2,I)
            ALPHA(1,3,I)=UV(16,I-1,L)*T_R(6,I,L)
     &         +LAMBDA(3,I)*GAMMA(3,3,I)+LAMBDA(2,I)*GAMMA(2,3,I)
     &         +LAMBDA(1,I)*GAMMA(1,3,I)
            G(1,I)=UV(10,I-1,L)*T_R(8,I,L)
     &           +UV(13,I-1,L)*T_R(10,I,L)
     &           +UV(16,I-1,L)*T_R(12,I,L)
     &           +LAMBDA(3,I)*H(3,I)
     &           +LAMBDA(2,I)*H(2,I)+LAMBDA(1,I)*H(1,I)
!
         ENDDO
!
!     THE LAYER ABOVE THE CLOUD: ONLY ONE SET OF ALPHAS IS NOW NEEDED.
!
      I=N_CLOUD_TOP-1
!
         IF (N_CLOUD_TOP.LT.N_LAYER) THEN
!           IF THERE IS NO CLOUD IN THE COLUMN THE V'S WILL NOT BE
!           ASSIGNED SO AN IF TEST IS REQUIRED.
            THETA(1,1,I)=ALPHA(1,1,I+1)*UV(1,I,L)
     &          +ALPHA(1,2,I+1)*UV(2,I,L)
     &          +ALPHA(1,3,I+1)*UV(3,I,L)
         ELSE
            THETA(1,1,I)=ALPHA(1,1,I+1)
         ENDIF
!
         BETA(1,1,I)=1.0E+00/(1.0E+00-THETA(1,1,I)*T_R(4,I,L))
         GAMMA(1,1,I)=THETA(1,1,I)*T_R(1,I,L)
         H(1,I)=G(1,I+1)+THETA(1,1,I)*T_R(7,I,L)
!
         LAMBDA(0,I)=T_R(1,I,L)*BETA(1,1,I)
         ALPHA(1,1,I)=T_R(4,I,L)+LAMBDA(0,I)*GAMMA(1,1,I)
         G(1,I)=T_R(8,I,L)+LAMBDA(0,I)*H(1,I)
!
!
         DO I=N_CLOUD_TOP-2, 1, -1
!
            BETA(1,1,I)=1.0E+00/(1.0E+00-ALPHA(1,1,I+1)*T_R(4,I,L))
            GAMMA(1,1,I)=ALPHA(1,1,I+1)*T_R(1,I,L)
            H(1,I)=G(1,I+1)+ALPHA(1,1,I+1)*T_R(7,I,L)
!
            LAMBDA(1,I)=T_R(1,I,L)*BETA(1,1,I)
            ALPHA(1,1,I)=T_R(4,I,L)+LAMBDA(1,I)*GAMMA(1,1,I)
            G(1,I)=T_R(8,I,L)+LAMBDA(1,I)*H(1,I)
!
         ENDDO

!
!
!     INITIALIZE FOR DOWNWARD BACK-SUBSTITUTION.
         FLUX_TOTAL(L, 2)=FLUX_INC_DOWN(L)
         FLUX_TOTAL(L, 1)=ALPHA(1,1,1)*FLUX_TOTAL(L, 2)+G(1,1)
!
!     SWEEP DOWNWARD THROUGH THE CLEAR-SKY REGION, FINDING THE DOWNWARD
!     FLUX AT THE TOP OF THE LAYER AND THE UPWARD FLUX AT THE BOTTOM.
      DO I=1, N_CLOUD_TOP-1
            FLUX_TOTAL(L, 2*I+1)=(GAMMA(1,1,I)*FLUX_TOTAL(L, 2*I)
     &         +H(1,I))*BETA(1,1,I)
            FLUX_TOTAL(L, 2*I+2)=T_R(1,I,L)*FLUX_TOTAL(L, 2*I)
     &        +T_R(4,I,L)*FLUX_TOTAL(L, 2*I+1)+T_R(7,I,L)
      ENDDO
!
!     PASS INTO THE TOP CLOUDY LAYER. USE FLUX_DOWN_[1,2,3] TO HOLD,
!     PROVISIONALLY, THE DOWNWARD FLUXES JUST BELOW THE TOP OF THE
!     LAYER, THEN CALCULATE THE UPWARD FLUXES AT THE BOTTOM AND
!     FINALLY THE DOWNWARD FLUXES AT THE BOTTOM OF THE LAYER.
      I=N_CLOUD_TOP
         FLUX_TEMP(1,I)=UV(1,I-1,L)*FLUX_TOTAL(L, 2*I)
         FLUX_TEMP(2,I)=UV(2,I-1,L)*FLUX_TOTAL(L, 2*I)
         FLUX_TEMP(3,I)=UV(3,I-1,L)*FLUX_TOTAL(L, 2*I)
         FLUX_TEMP(4,I)=(GAMMA(1,1,I)*FLUX_TEMP(1,I)
     &      +GAMMA(1,2,I)*FLUX_TEMP(2,I)
     &      +GAMMA(1,3,I)*FLUX_TEMP(3,I)
     &      +H(1,I))*BETA(1,1,I)
         FLUX_TEMP(5,I)=(GAMMA(2,1,I)*FLUX_TEMP(1,I)
     &      +GAMMA(2,2,I)*FLUX_TEMP(2,I)
     &      +GAMMA(2,3,I)*FLUX_TEMP(3,I)+H(2,I)
     &      -BETA(2,1,I)*FLUX_TEMP(4,I))*BETA(2,2,I)
         FLUX_TEMP(6,I)=(GAMMA(3,1,I)*FLUX_TEMP(1,I)
     &      +GAMMA(3,2,I)*FLUX_TEMP(2,I)
     &      +GAMMA(3,3,I)*FLUX_TEMP(3,I)+H(3,I)
     &      -BETA(3,1,I)*FLUX_TEMP(4,I)-BETA(3,2,I)*FLUX_TEMP(5,I))
     &      *BETA(3,3,I)
         FLUX_TEMP(1,I)=T_R(1,I,L)*FLUX_TEMP(1,I)
     &      +T_R(4,I,L)*FLUX_TEMP(4,I)+T_R(7,I,L)
         FLUX_TEMP(2,I)=T_R(2,I,L)*FLUX_TEMP(2,I)
     &      +T_R(5,I,L)*FLUX_TEMP(5,I)+T_R(9,I,L)
         FLUX_TEMP(3,I)=T_R(3,I,L)*FLUX_TEMP(3,I)
     &      +T_R(6,I,L)*FLUX_TEMP(6,I)+T_R(11,I,L)
!
!     THE MAIN LOOP OF BACK-SUBSTITUTION. THE PROVISIONAL USE OF THE
!     DOWNWARD FLUXES IS AS ABOVE.
         DO I=N_CLOUD_TOP+1, N_LAYER
            FLUX_TEMP(1,I)=UV(1,I-1,L)*FLUX_TEMP(1,I-1)
     &           +UV(4,I-1,L)*FLUX_TEMP(2,I-1)
     &           +UV(7,I-1,L)*FLUX_TEMP(3,I-1)
            FLUX_TEMP(2,I)=UV(2,I-1,L)*FLUX_TEMP(1,I-1)
     &           +UV(5,I-1,L)*FLUX_TEMP(2,I-1)
     &           +UV(8,I-1,L)*FLUX_TEMP(3,I-1)
            FLUX_TEMP(3,I)=UV(3,I-1,l)*FLUX_TEMP(1,I-1)
     &           +UV(6,I-1,L)*FLUX_TEMP(2,I-1)
     &           +UV(9,I-1,L)*FLUX_TEMP(3,I-1)
            FLUX_TEMP(4,I)=(GAMMA(1,1,I)*FLUX_TEMP(1,I)
     &         +GAMMA(1,2,I)*FLUX_TEMP(2,I)
     &         +GAMMA(1,3,I)*FLUX_TEMP(3,I)
     &         +H(1,I))*BETA(1,1,I)
            FLUX_TEMP(5,I)=(GAMMA(2,1,I)*FLUX_TEMP(1,I)
     &         +GAMMA(2,2,I)*FLUX_TEMP(2,I)
     &         +GAMMA(2,3,I)*FLUX_TEMP(3,I)+H(2,I)
     &         -BETA(2,1,I)*FLUX_TEMP(4,I))*BETA(2,2,I)
            FLUX_TEMP(6,I)=(GAMMA(3,1,I)*FLUX_TEMP(1,I)
     &         +GAMMA(3,2,I)*FLUX_TEMP(2,I)
     &         +GAMMA(3,3,I)*FLUX_TEMP(3,I)+H(3,I)
     &         -BETA(3,1,I)*FLUX_TEMP(4,I)
     &         -BETA(3,2,I)*FLUX_TEMP(5,I))
     &         *BETA(3,3,I)
            FLUX_TEMP(1,I)=T_R(1,I,L)*FLUX_TEMP(1,I)
     &         +T_R(4,I,L)*FLUX_TEMP(4,I)+T_R(7,I,L)
            FLUX_TEMP(2,I)=T_R(2,I,L)*FLUX_TEMP(2,I)
     &         +T_R(5,I,L)*FLUX_TEMP(5,I)+T_R(9,I,L)
            FLUX_TEMP(3,I)=T_R(3,I,L)*FLUX_TEMP(3,I)
     &         +T_R(6,I,L)*FLUX_TEMP(6,I)+T_R(11,I,L)
         ENDDO


!
!
!     CALCULATE THE OVERALL FLUX.
         DO I=N_CLOUD_TOP, N_LAYER
            FLUX_TOTAL(L, 2*I+1)=FLUX_TEMP(4,I)+FLUX_TEMP(5,I)
     &         +FLUX_TEMP(6,I)
            FLUX_TOTAL(L, 2*I+2)=FLUX_TEMP(1,I)+FLUX_TEMP(2,I)
     &         +FLUX_TEMP(3,I)
         ENDDO

      END DO
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
      RETURN
      END
