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
!+ Subroutine to solve for triple overlaps with approximate scattering.
!
! Method:
!       The flux is propagated downwards, ignoring reflection terms.
!       Since the routine uses differential fluxes, this effectively
!       treats the upward flux as Planckian at this point. Upward
!       fluxes are calculated using the newly available approximate
!       downward fluxes in the reflected terms.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.2             10-04-96                Original Code
!                                               (J. M. Edwards)
!       4.5             27-05-98                Non-scientific, but
!                                               non-bit-comparable
!                                               change to the indexing
!                                               of the DO-loops at
!                                               SLVTPA3A.263 and 266
!                                               to use the explicit
!                                               surface index.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLVER_TRIPLE_APP_SCAT(N_PROFILE, N_LAYER, N_CLOUD_TOP
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
     &   , FLUX_PROPAG_1(NPD_PROFILE)
!             TEMPORARY FLUXES FOR PROPAGATION ACROSS LAYERS
     &   , FLUX_PROPAG_2(NPD_PROFILE)
!             TEMPORARY FLUXES FOR PROPAGATION ACROSS LAYERS
     &   , FLUX_PROPAG_3(NPD_PROFILE)
!             TEMPORARY FLUXES FOR PROPAGATION ACROSS LAYERS
!
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
         FLUX_TOTAL(L, 2)=FLUX_INC_DOWN(L)
      ENDDO
      DO I=1, N_CLOUD_TOP-1
         DO L=1, N_PROFILE
            FLUX_TOTAL(L, 2*I+2)=T(L, I)*FLUX_TOTAL(L, 2*I)
     &         +S_DOWN(L, I)
         ENDDO
      ENDDO
!
!     PASS INTO THE CLOUDY REGION. HERE, DOWNWARD FLUXES HOLD VALUES
!     JUST BELOW THE LEVEL AND UPWARD FLUXES THE VALUES JUST ABOVE IT.
!     THUS THE FLUXES IMPINGING ON THE LAYER ARE HELD.
      I=N_CLOUD_TOP-1
      DO L=1, N_PROFILE
         FLUX_DOWN_1(L, I)=V11(L, I)*FLUX_TOTAL(L, 2*I+2)
         FLUX_DOWN_2(L, I)=V21(L, I)*FLUX_TOTAL(L, 2*I+2)
         FLUX_DOWN_3(L, I)=V31(L, I)*FLUX_TOTAL(L, 2*I+2)
      ENDDO
!
      DO I=N_CLOUD_TOP, N_LAYER-1
         DO L=1, N_PROFILE
!
!           PROPAGTE THE FLUX ACROSS THE LAYER.
            FLUX_PROPAG_1(L)=T(L, I)*FLUX_DOWN_1(L, I-1)
     &         +S_DOWN(L, I)
            FLUX_PROPAG_2(L)=T_STRAT(L, I)*FLUX_DOWN_2(L, I-1)
     &         +S_DOWN_STRAT(L, I)
            FLUX_PROPAG_3(L)=T_CONV(L, I)*FLUX_DOWN_3(L, I-1)
     &         +S_DOWN_CONV(L, I)
!
!           TRANSFER ACROSS THE INTERFACE.
            FLUX_DOWN_1(L, I)=V11(L, I)*FLUX_PROPAG_1(L)
     &         +V12(L, I)*FLUX_PROPAG_2(L)
     &         +V13(L, I)*FLUX_PROPAG_3(L)
            FLUX_DOWN_2(L, I)=V21(L, I)*FLUX_PROPAG_1(L)
     &         +V22(L, I)*FLUX_PROPAG_2(L)
     &         +V23(L, I)*FLUX_PROPAG_3(L)
            FLUX_DOWN_3(L, I)=V31(L, I)*FLUX_PROPAG_1(L)
     &         +V32(L, I)*FLUX_PROPAG_2(L)
     &         +V33(L, I)*FLUX_PROPAG_3(L)
!
         ENDDO
      ENDDO
!
!     PROPAGATE ACROSS THE BOTTOM LAYER AND FORM THE REFLECTED BEAM.
!     WE DO NOT TRANSFER FLUXES ACROSS THE BOTTOM INTERFACE, SO AS
!     TO MAKE THE REFLECTION CONSISTENT BETWEEN REGIONS.
      DO L=1, N_PROFILE
!
!        PROPAGTE THE FLUX THROUGH THE LAYER.
         FLUX_DOWN_1(L, N_LAYER)
     &      =T(L, N_LAYER)*FLUX_DOWN_1(L, N_LAYER-1)
     &      +S_DOWN(L, N_LAYER)
         FLUX_DOWN_2(L, N_LAYER)
     &      =T_STRAT(L, N_LAYER)*FLUX_DOWN_2(L, N_LAYER-1)
     &      +S_DOWN_STRAT(L, N_LAYER)
         FLUX_DOWN_3(L, N_LAYER)
     &      =T_CONV(L, N_LAYER)*FLUX_DOWN_3(L, N_LAYER-1)
     &      +S_DOWN_CONV(L, N_LAYER)
!
!        REFLECT FROM THE SURFACE.
         FLUX_UP_1(L, N_LAYER)
     &      =ALBEDO_SURFACE_DIFF(L)*FLUX_DOWN_1(L, N_LAYER)
     &      +SOURCE_GROUND_FREE(L)
         FLUX_UP_2(L, N_LAYER)
     &      =ALBEDO_SURFACE_DIFF(L)*FLUX_DOWN_2(L, N_LAYER)
     &      +SOURCE_GROUND_STRAT(L)
         FLUX_UP_3(L, N_LAYER)
     &      =ALBEDO_SURFACE_DIFF(L)*FLUX_DOWN_3(L, N_LAYER)
     &      +SOURCE_GROUND_CONV(L)
!
!        PROPAGATE ACROSS THE BOTTOM LAYER.
         FLUX_PROPAG_1(L)
     &      =T(L, N_LAYER)*FLUX_UP_1(L, N_LAYER)+S_UP(L, N_LAYER)
     &      +R(L, N_LAYER)*FLUX_DOWN_1(L, N_LAYER-1)
         FLUX_PROPAG_2(L)
     &      =T_STRAT(L, N_LAYER)*FLUX_UP_2(L, N_LAYER)
     &      +S_UP_STRAT(L, N_LAYER)
     &      +R_STRAT(L, N_LAYER)*FLUX_DOWN_2(L, N_LAYER-1)
         FLUX_PROPAG_3(L)
     &      =T_CONV(L, N_LAYER)*FLUX_UP_3(L, N_LAYER)
     &      +S_UP_CONV(L, N_LAYER)
     &      +R_CONV(L, N_LAYER)*FLUX_DOWN_3(L, N_LAYER-1)
!
      ENDDO
!
!
!
!     WORK BACK UP THROUGH THE COLUMN ASSIGNING THE UPWARD FLUXES.
      DO I=N_LAYER-1, N_CLOUD_TOP, -1
         DO L=1, N_PROFILE
!
            FLUX_UP_1(L, I)=U11(L, I)*FLUX_PROPAG_1(L)
     &         +U12(L, I)*FLUX_PROPAG_2(L)
     &         +U13(L, I)*FLUX_PROPAG_3(L)
            FLUX_UP_2(L, I)=U21(L, I)*FLUX_PROPAG_1(L)
     &         +U22(L, I)*FLUX_PROPAG_2(L)
     &         +U23(L, I)*FLUX_PROPAG_3(L)
            FLUX_UP_3(L, I)=U31(L, I)*FLUX_PROPAG_1(L)
     &         +U32(L, I)*FLUX_PROPAG_2(L)
     &         +U33(L, I)*FLUX_PROPAG_3(L)
!
            FLUX_PROPAG_1(L)=T(L, I)*FLUX_UP_1(L, I)+S_UP(L, I)
     &         +R(L, I)*FLUX_DOWN_1(L, I-1)
            FLUX_PROPAG_2(L)=T_STRAT(L, I)*FLUX_UP_2(L, I)
     &         +S_UP_STRAT(L, I)+R_STRAT(L, I)*FLUX_DOWN_2(L, I-1)
            FLUX_PROPAG_3(L)=T_CONV(L, I)*FLUX_UP_3(L, I)
     &         +S_UP_CONV(L, I)+R_CONV(L, I)*FLUX_DOWN_3(L, I-1)
!
         ENDDO
      ENDDO
!
!     PROPAGATE INTO THE CLOUD-FREE REGION.
      I=N_CLOUD_TOP-1
      DO L=1, N_PROFILE
         FLUX_TOTAL(L, 2*I+1)=FLUX_PROPAG_1(L)+FLUX_PROPAG_2(L)
     &      +FLUX_PROPAG_3(L)
      ENDDO
!
!     CONTINUE THROUGH THE LAYERS ABOVE CLOUDS.
      DO I=N_CLOUD_TOP-1, 1, -1
         DO L=1, N_PROFILE
            FLUX_TOTAL(L, 2*I-1)=T(L, I)*FLUX_TOTAL(L, 2*I+1)
     &         +R(L, I)*FLUX_TOTAL(L, 2*I)+S_UP(L, I)
         ENDDO
      ENDDO
!
!     ASSIGN THE TOTAL FLUXES ON THE INTERMEDIATE CLOUDY LAYERS.
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
               FLUX_TOTAL(L, I+1)=FLUX_TOTAL(L, 2*I+2)
     &            -FLUX_TOTAL(L, 2*I-1)
            ENDDO
         ENDDO
      ENDIF

!
!
!
      RETURN
      END
