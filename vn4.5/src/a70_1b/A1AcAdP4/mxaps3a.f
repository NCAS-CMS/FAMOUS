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
!+ Subroutine to solve for fluxes treating scattering approximately.
!
! Method:
!       The routine is applicable in the infra-red. Downward
!       differential fluxes are calculated first assuming that the
!       upward differential fluxes are 0. Upward fluxes are then
!       calculated using the previously calculated downward fluxes
!       in the reflected terms.
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
      SUBROUTINE MIX_APP_SCAT(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &   , T_FREE, R_FREE, S_DOWN_FREE, S_UP_FREE
     &   , T_CLOUD, R_CLOUD, S_DOWN_CLOUD, S_UP_CLOUD
     &   , G_FF, G_FC, G_CF, G_CC
     &   , B_FF, B_FC, B_CF, B_CC
     &   , L_NET
     &   , FLUX_INC_DOWN
     &   , SOURCE_GROUND, ALBEDO_SURFACE_DIFF
     &   , FLUX_DIFFUSE
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
     &     T_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE TRANSMISSION
     &   , R_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE REFLECTION
     &   , S_DOWN_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE DOWNWARD SOURCE FUNCTION
     &   , S_UP_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE UPWARD SOURCE FUNCTION
     &   , T_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY TRANSMISSION
     &   , R_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY REFLECTION
     &   , S_DOWN_CLOUD(NPD_PROFILE, NPD_LAYER)
!             DOWNWARD CLOUDY SOURCE FUNCTION
     &   , S_UP_CLOUD(NPD_PROFILE, NPD_LAYER)
!             UPWARD CLOUDY SOURCE FUNCTION
      REAL      !, INTENT(IN)
     &     G_FF(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , G_FC(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , G_CF(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , G_CC(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , B_FF(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , B_FC(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , B_CF(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
     &   , B_CC(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT
      REAL      !, INTENT(IN)
     &     FLUX_INC_DOWN(NPD_PROFILE)
!             INCIDENT DIFFUSE FLUX
     &   , SOURCE_GROUND(NPD_PROFILE)
!             SOURCE FROM GROUND
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)
!             DIFFUSE ALBEDO
      REAL      !, INTENT(OUT)
     &     FLUX_DIFFUSE(NPD_PROFILE, 2*NPD_LAYER+2)
!             DIFFUSE FLUX
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
      REAL
     &     FLUX_DOWN(NPD_PROFILE, 0: NPD_LAYER)
!             DOWNWARD FLUXES OUTSIDE CLOUDS JUST BELOW I'TH LEVEL
     &   , FLUX_DOWN_CLOUD(NPD_PROFILE, 0: NPD_LAYER)
!             DOWNWARD FLUXES INSIDE CLOUDS JUST BELOW I'TH LEVEL
     &   , FLUX_UP(NPD_PROFILE, 0: NPD_LAYER)
!             UPWARD FLUXES OUTSIDE CLOUDS JUST ABOVE I'TH LEVEL
     &   , FLUX_UP_CLOUD(NPD_PROFILE, 0: NPD_LAYER)
!             UPWARD FLUXES INSIDE CLOUDS JUST ABOVE I'TH LEVEL
     &   , FLUX_PROPAGATED
!             TEMPORARY PROPAGATED FLUX OUTSIDE CLOUD
     &   , FLUX_PROPAGATED_CLOUD
!             TEMPORARY PROPAGATED FLUX INSIDE CLOUD
     &   , FLUX_CLOUD_TOP(NPD_PROFILE)
!             TOTAL DOWNWARD FLUX AT TOP OF CLOUD
!
!
!
!     THE ARRAYS FLUX_DOWN AND FLUX_UP WILL EVENTUALLY CONTAIN THE TOTAL
!     FLUXES, BUT INITIALLY THEY ARE USED FOR THE CLEAR FLUXES.
!     NOTE THAT DOWNWARD FLUXES REFER TO VALUES JUST BELOW THE INTERFACE
!     AND UPWARD FLUXES TO VALUES JUST ABOVE IT.
!
!
!     DOWNWARD FLUX:
!
!     REGION ABOVE CLOUDS:
      DO L=1, N_PROFILE
         FLUX_DOWN(L, 0)=FLUX_INC_DOWN(L)
      ENDDO
      DO I=1, N_CLOUD_TOP-1
         DO L=1, N_PROFILE
            FLUX_DOWN(L, I)=T_FREE(L, I)*FLUX_DOWN(L, I-1)
     &         +S_DOWN_FREE(L, I)
         ENDDO
      ENDDO
      DO L=1, N_PROFILE
         FLUX_CLOUD_TOP(L)=FLUX_DOWN(L, N_CLOUD_TOP-1)
      ENDDO
!
!     REGION OF CLOUDS:
      DO L=1, N_PROFILE
         FLUX_DOWN(L, N_CLOUD_TOP-1)
     &      =G_FF(L, N_CLOUD_TOP-1)*FLUX_CLOUD_TOP(L)
         FLUX_DOWN_CLOUD(L, N_CLOUD_TOP-1)
     &      =G_FC(L, N_CLOUD_TOP-1)*FLUX_CLOUD_TOP(L)
      ENDDO
!
      DO I=N_CLOUD_TOP, N_LAYER-1
         DO L=1, N_PROFILE
!
!           PROPAGATE DOWNWARD FLUXES THROUGH THE LAYER.
            FLUX_PROPAGATED=T_FREE(L, I)*FLUX_DOWN(L, I-1)
     &         +S_DOWN_FREE(L, I)
            FLUX_PROPAGATED_CLOUD=T_CLOUD(L, I)*FLUX_DOWN_CLOUD(L, I-1)
     &         +S_DOWN_CLOUD(L, I)
!           TRANSFER DOWNWARD FLUXES ACROSS THE INTERFACE.
            FLUX_DOWN(L, I)
     &         =G_FF(L, I)*FLUX_PROPAGATED
     &         +G_CF(L, I)*FLUX_PROPAGATED_CLOUD
            FLUX_DOWN_CLOUD(L, I)
     &         =G_CC(L, I)*FLUX_PROPAGATED_CLOUD
     &         +G_FC(L, I)*FLUX_PROPAGATED
!
         ENDDO
      ENDDO
!
!     PROPAGATE ACROSS THE BOTTOM LAYER, BUT WITHOUT TRANSFERRING
!     ACROSS THE SURFACE AND FORM THE REFLECTED BEAMS.
      DO L=1, N_PROFILE
!        PROPAGATE DOWNWARD FLUXES THROUGH THE LAYER.
         FLUX_DOWN(L, N_LAYER)
     &      =T_FREE(L, N_LAYER)*FLUX_DOWN(L, N_LAYER-1)
     &      +S_DOWN_FREE(L, N_LAYER)
         FLUX_DOWN_CLOUD(L, N_LAYER)
     &      =T_CLOUD(L, N_LAYER)*FLUX_DOWN_CLOUD(L, N_LAYER-1)
     &      +S_DOWN_CLOUD(L, N_LAYER)
         FLUX_UP(L, N_LAYER)
     &      =ALBEDO_SURFACE_DIFF(L)*FLUX_DOWN(L, N_LAYER)
     &      +B_FF(L, N_LAYER)*SOURCE_GROUND(L)
         FLUX_UP_CLOUD(L, N_LAYER)
     &      =ALBEDO_SURFACE_DIFF(L)*FLUX_DOWN_CLOUD(L, N_LAYER)
     &      +B_CF(L, N_LAYER)*SOURCE_GROUND(L)
      ENDDO
!
!
!     CALCULATE THE UPWARD FLUXES USING THE PREVIOUS DOWNWARD FLUXES
!     TO APPROXIMATE THE SCATTERING TERM.
      DO I=N_LAYER, N_CLOUD_TOP, -1
         DO L=1, N_PROFILE
!
!           PROPAGATE UPWARD FLUXES THROUGH THE LAYER.
            FLUX_PROPAGATED=T_FREE(L, I)*FLUX_UP(L, I)+S_UP_FREE(L, I)
     &         +R_FREE(L, I)*FLUX_DOWN(L, I-1)
            FLUX_PROPAGATED_CLOUD=T_CLOUD(L, I)*FLUX_UP_CLOUD(L, I)
     &         +S_UP_CLOUD(L, I)+R_CLOUD(L, I)*FLUX_DOWN_CLOUD(L, I-1)
!           TRANSFER UPWARD FLUXES ACROSS THE INTERFACE.
            FLUX_UP(L, I-1)=B_FF(L, I-1)*FLUX_PROPAGATED
     &         +B_FC(L, I-1)*FLUX_PROPAGATED_CLOUD
            FLUX_UP_CLOUD(L, I-1)=B_CC(L, I-1)*FLUX_PROPAGATED_CLOUD
     &         +B_CF(L, I-1)*FLUX_PROPAGATED
!
         ENDDO
      ENDDO
!
!     CONTINUE THROUGH THE REGION ABOVE CLOUDS.
      DO I=N_CLOUD_TOP-1, 1, -1
         DO L=1, N_PROFILE
            FLUX_UP(L, I-1)=T_FREE(L, I)*FLUX_UP(L,I)+S_UP_FREE(L, I)
     &         +R_FREE(L, I)*FLUX_DOWN(L, I-1)
         ENDDO
      ENDDO
!
!
!
!     CALCULATE THE OVERALL FLUX.
      IF (L_NET) THEN
         DO I=0, N_CLOUD_TOP-2
            DO L=1, N_PROFILE
               FLUX_DIFFUSE(L, I+1)=FLUX_DOWN(L, I)-FLUX_UP(L, I)
            ENDDO
         ENDDO
         DO I=N_CLOUD_TOP-1, N_LAYER
            DO L=1, N_PROFILE
               FLUX_DIFFUSE(L, I+1)
     &            =FLUX_DOWN(L, I)+FLUX_DOWN_CLOUD(L, I)
     &            -FLUX_UP(L, I)-FLUX_UP_CLOUD(L, I)
            ENDDO
         ENDDO
      ELSE
         DO I=0, N_CLOUD_TOP-2
            DO L=1, N_PROFILE
               FLUX_DIFFUSE(L, 2*I+1)=FLUX_UP(L, I)
               FLUX_DIFFUSE(L, 2*I+2)=FLUX_DOWN(L, I)
            ENDDO
         ENDDO
         DO I=N_CLOUD_TOP-1, N_LAYER
            DO L=1, N_PROFILE
               FLUX_DIFFUSE(L, 2*I+1)=FLUX_UP(L, I)+FLUX_UP_CLOUD(L, I)
               FLUX_DIFFUSE(L, 2*I+2)=FLUX_DOWN(L, I)
     &            +FLUX_DOWN_CLOUD(L, I)
            ENDDO
         ENDDO
      ENDIF
!
!
      RETURN
      END
