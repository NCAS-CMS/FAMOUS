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
!+ Subroutine to zero an array.
!
! Purpose:
!   The routine fills a 1-dimensional array with zeros.
!
! Method:
!   Straightforward.
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
      SUBROUTINE R2_ZERO_1D(N, X)
!
!
!
      IMPLICIT NONE
!
!
!     DUMMY ARGUMENTS
      INTEGER   !, INTENT(IN)
     &     N
!             LENGTH OF ARRAY
      REAL      !, INTENT(OUT)
     &     X(N)
!             ARRAY TO BE ZEROED
!
!     LOCAL VARIABLES
      INTEGER
     &     I
!             LOOP VARIABLE
!
!
!
      DO I=1, N
         X(I)=0.0E+00
      ENDDO
!
!
!
      RETURN
      END
!+ Subroutine to initialize diagnostics and coupling arrays.
!
! Purpose:
!   The coupling and diagnostic arrays are zeroed.
!
! Method:
!   Straightforward.
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
      SUBROUTINE R2_INIT_COUPLE_DIAG(N_PROFILE
     &   , SEA_FLUX
     &   , L_SURFACE_DOWN_FLUX, SURFACE_DOWN_FLUX
     &   , L_SURF_DOWN_CLR, SURF_DOWN_CLR
     &   , L_SURF_UP_CLR, SURF_UP_CLR
     &   , L_FLUX_BELOW_690NM_SURF, FLUX_BELOW_690NM_SURF
     &   , NPD_PROFILE
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     DUMMY ARGUMENTS
!
!     DIMENSIONS OF ARRAYS
      INTEGER   !, INTENT(IN)
     &     NPD_PROFILE
!             MAXIMUM NUMBER OF ATMOSPHERIC PROFILES
!
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF ATMOSPHERIC PROFILES
!
!     SWITCHES FOR DIAGNOSTICS:
      LOGICAL   !, INTENT(IN)
     &     L_FLUX_BELOW_690NM_SURF
!             FLUX BELOW 690NM AT SURFACE TO BE CALCULATED
     &   , L_SURFACE_DOWN_FLUX
!             DOWNWARD SURFACE FLUX REQUIRED
     &   , L_SURF_DOWN_CLR
!             CALCULATE DOWNWARD CLEAR FLUX
     &   , L_SURF_UP_CLR
!             CALCULATE UPWARD CLEAR FLUX
!
!     SURFACE FLUXES FOR COUPLING OR DIAGNOSTIC USE
      REAL      !, INTENT(OUT)
     &     SEA_FLUX(NPD_PROFILE)
!             NET DOWNWARD FLUX INTO SEA
     &   , SURFACE_DOWN_FLUX(NPD_PROFILE)
!             DOWNWARD FLUX AT SURFACE
     &   , SURF_DOWN_CLR(NPD_PROFILE)
!             CLEAR-SKY DOWNWARD FLUX AT SURFACE
     &   , SURF_UP_CLR(NPD_PROFILE)
!             CLEAR-SKY UPWARD FLUX AT SURFACE
     &   , FLUX_BELOW_690NM_SURF(NPD_PROFILE)
!             SURFACE FLUX BELOW 690NM
!
!
!
      CALL R2_ZERO_1D(N_PROFILE, SEA_FLUX)
!
      IF (L_SURFACE_DOWN_FLUX) THEN
         CALL R2_ZERO_1D(N_PROFILE, SURFACE_DOWN_FLUX)
      ENDIF
!
      IF (L_SURF_DOWN_CLR) THEN
         CALL R2_ZERO_1D(N_PROFILE, SURF_DOWN_CLR)
      ENDIF
!
      IF (L_SURF_UP_CLR) THEN
         CALL R2_ZERO_1D(N_PROFILE, SURF_UP_CLR)
      ENDIF
!
      IF (L_FLUX_BELOW_690NM_SURF) THEN
         CALL R2_ZERO_1D(N_PROFILE, FLUX_BELOW_690NM_SURF)
      ENDIF
!
!
!
      RETURN
      END
!+ Subroutine to calculate spectral diagnostics and coupling arrays.
!
! Purpose:
!   The coupling and diagnostic arrays are calculated.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                Formulation over
!                                               sea-ice revised.
!                                               Corrections to
!                                               some diagnostics.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_COUPLE_DIAG(N_PROFILE, L_NET, ISOLIR
     &   , ALBEDO_FIELD_DIFF, ALBEDO_FIELD_DIR
     &   , ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR
     &   , N_FRAC_ICE_POINT, I_FRAC_ICE_POINT, ICE_FRACTION
     &   , PLANCK_FREEZE_SEA
     &   , PLANCK_AIR_SURFACE, THERMAL_SOURCE_GROUND
     &   , FLUX_DOWN, FLUX_UP, FLUX_DIRECT
     &   , FLUX_DOWN_CLEAR, FLUX_UP_CLEAR, FLUX_DIRECT_CLEAR
     &   , WEIGHT_690NM
     &   , SEA_FLUX
     &   , L_SURFACE_DOWN_FLUX, SURFACE_DOWN_FLUX
     &   , L_SURF_DOWN_CLR, SURF_DOWN_CLR
     &   , L_SURF_UP_CLR, SURF_UP_CLR
     &   , L_FLUX_BELOW_690NM_SURF, FLUX_BELOW_690NM_SURF
     &   , NPD_PROFILE
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED
!     SPECTRAL REGIONS
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
!
!     DUMMY ARGUMENTS
!
!     DIMENSIONS OF ARRAYS
      INTEGER   !, INTENT(IN)
     &     NPD_PROFILE
!             MAXIMUM NUMBER OF ATMOSPHERIC PROFILES
!
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF ATMOSPHERIC PROFILES
     &   , ISOLIR
!             SPECTRAL REGION
!
!     LOGICAL SWITCHES FOR THE CODE
      LOGICAL   !, INTENT(IN)
     &     L_NET
!             FLAG FOR NET FLUXES
!
!     SWITCHES FOR DIAGNOSTICS:
      LOGICAL   !, INTENT(IN)
     &     L_FLUX_BELOW_690NM_SURF
!             FLUX BELOW 690NM AT SURFACE TO BE CALCULATED
     &   , L_SURFACE_DOWN_FLUX
!             DOWNWARD SURFACE FLUX REQUIRED
     &   , L_SURF_DOWN_CLR
!             CALCULATE DOWNWARD CLEAR FLUX
     &   , L_SURF_UP_CLR
!             CALCULATE UPWARD CLEAR FLUX
!
!     ALBEDOS
      REAL      !, INTENT(IN)
     &     ALBEDO_FIELD_DIFF(NPD_PROFILE)
!             DIFFUSE ALBEDO MEANED OVER GRID BOX
     &   , ALBEDO_FIELD_DIR(NPD_PROFILE)
!             DIRECT ALBEDO MEANED OVER GRID BOX
     &   , ALBEDO_SEA_DIFF(NPD_PROFILE)
!             DIFFUSE ALBEDO OF OPEN SEA
     &   , ALBEDO_SEA_DIR(NPD_PROFILE)
!             DIRECT ALBEDO MEANED OF OPEN SEA
!
      REAL      !, INTENT(IN)
     &     THERMAL_SOURCE_GROUND(NPD_PROFILE)
!             THERMAL SOURCE AT GROUND
     &   , PLANCK_AIR_SURFACE(NPD_PROFILE)
!             PLANCK FUNCTION AT NEAR-SURFACE AIR TEMPERATURE IN BAND
!
!     ARGUMENTS RELATING TO SEA ICE.
      INTEGER   !, INTENT(IN)
     &     N_FRAC_ICE_POINT
!             NUMBER OF POINTS WITH FRACTIONAL ICE COVER
     &   , I_FRAC_ICE_POINT(NPD_PROFILE)
!             INDICES OF POINTS WITH FRACTIONAL ICE COVER
      REAL  !, INTENT(IN)
     &     ICE_FRACTION(NPD_PROFILE)
!             ICE FRACTION
      REAL  !, INTENT(IN)
     &     PLANCK_FREEZE_SEA
!             PLANCK FUNCTION OVER FREEZING SEA
!
      REAL      !, INTENT(IN)
     &     WEIGHT_690NM
!             WEIGHTING APPLIED TO BAND FOR REGION BELOW 690 NM
!
!     CALCULATED FLUXES
      REAL      !, INTENT(IN)
     &     FLUX_DOWN(NPD_PROFILE)
!             TOTAL DOWNWARD OR NET FLUX AT SURFACE
     &   , FLUX_DIRECT(NPD_PROFILE)
!             DIRECT SOLAR FLUX AT SURFACE
     &   , FLUX_UP(NPD_PROFILE)
!             UPWARD FLUX AT SURFACE
     &   , FLUX_DOWN_CLEAR(NPD_PROFILE)
!             TOTAL CLEAR-SKY DOWNWARD OR NET FLUX AT SURFACE
     &   , FLUX_UP_CLEAR(NPD_PROFILE)
!             CLEAR-SKY UPWARD FLUX AT SURFACE
     &   , FLUX_DIRECT_CLEAR(NPD_PROFILE)
!             CLEAR-SKY DIRECT SOLAR FLUX AT SURFACE
!
!
!     SURFACE FLUXES FOR COUPLING OR DIAGNOSTIC USE
      REAL      !, INTENT(INOUT)
     &     SEA_FLUX(NPD_PROFILE)
!             NET DOWNWARD FLUX INTO SEA
     &   , SURFACE_DOWN_FLUX(NPD_PROFILE)
!             DOWNWARD FLUX AT SURFACE
     &   , SURF_DOWN_CLR(NPD_PROFILE)
!             CLEAR-SKY DOWNWARD FLUX AT SURFACE
     &   , SURF_UP_CLR(NPD_PROFILE)
!             CLEAR-SKY UPWARD FLUX AT SURFACE
     &   , FLUX_BELOW_690NM_SURF(NPD_PROFILE)
!             SURFACE FLUX BELOW 690NM
!
!
!     LOCAL VARIABLES
      INTEGER
     &     L
!             LOOP VARIABLE
!
!
!
!
!     DEPENDING ON THE SOLVER THE TOTAL FLUX AVAILABLE WILL BE EITHER
!     THE NET FLUX OR THE SEPARATE UPWARD AND DOWNWARD FLUXES, HENCE
!     EACH DIAGNOSTIC MUST BE ENFOLDED IN AN IF-TEST.
!
!     SINCE DIFFERENTIAL FLUXES ARE USED IN THE INFRA-RED APPROPRIATE
!     PLANCKIAN SOURCES MUST BE ADDED TO NON-NET FLUXES. A SLIGHTLY
!     INEFFICIENT FORM HAS BEEN USED IN THE NON-NET CASE, DERIVED IN
!     ANALOGY WITH THE CASE OF NET FLUXES SINCE THIS MATCHES THE USE
!     OF ARRAYS IN THE MAIN CODE.
!
      IF (L_NET) THEN
         IF (ISOLIR.EQ.IP_SOLAR) THEN
            DO L=1, N_PROFILE
               SEA_FLUX(L)=SEA_FLUX(L)+FLUX_DIRECT(L)
     &            *(ALBEDO_SEA_DIFF(L)-ALBEDO_SEA_DIR(L))
     &            +((1.0E+00-ALBEDO_SEA_DIFF(L))
     &            /(1.0E+00-ALBEDO_FIELD_DIFF(L)))
     &            *(FLUX_DOWN(L)-FLUX_DIRECT(L)
     &            *(ALBEDO_FIELD_DIFF(L)-ALBEDO_FIELD_DIR(L)))
            ENDDO
         ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN
            DO L=1, N_PROFILE
               SEA_FLUX(L)=SEA_FLUX(L)
     &            +(FLUX_DOWN(L)+THERMAL_SOURCE_GROUND(L))
     &            *(1.0E+00-ALBEDO_SEA_DIFF(L))
     &            /(1.0E+00-ALBEDO_FIELD_DIFF(L))
     &            -(1.0E+00-ALBEDO_SEA_DIFF(L))*PLANCK_FREEZE_SEA
            ENDDO
         ENDIF
      ELSE
         IF (ISOLIR.EQ.IP_SOLAR) THEN
            DO L=1, N_PROFILE
               SEA_FLUX(L)=SEA_FLUX(L)+FLUX_DOWN(L)-FLUX_UP(L)
            ENDDO
         ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN
            DO L=1, N_PROFILE
               SEA_FLUX(L)=SEA_FLUX(L)
     &            +(1.0E+00-ALBEDO_SEA_DIFF(L))
     &            *(FLUX_DOWN(L)+PLANCK_AIR_SURFACE(L)
     &            -PLANCK_FREEZE_SEA)
            ENDDO
         ENDIF
      ENDIF
!
      IF (L_SURFACE_DOWN_FLUX) THEN
         IF (L_NET) THEN
            IF (ISOLIR.EQ.IP_SOLAR) THEN
               DO L=1, N_PROFILE
                  SURFACE_DOWN_FLUX(L)=SURFACE_DOWN_FLUX(L)
     &               +(FLUX_DOWN(L)+FLUX_DIRECT(L)
     &               *(ALBEDO_FIELD_DIR(L)-ALBEDO_FIELD_DIFF(L)))
     &               /(1.0E+00-ALBEDO_FIELD_DIFF(L))
               ENDDO
            ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN
               DO L=1, N_PROFILE
                  SURFACE_DOWN_FLUX(L)=SURFACE_DOWN_FLUX(L)
     &               +(FLUX_DOWN(L)+THERMAL_SOURCE_GROUND(L))
     &               /(1.0E+00-ALBEDO_FIELD_DIFF(L))
               ENDDO
            ENDIF
         ELSE
            IF (ISOLIR.EQ.IP_SOLAR) THEN
               DO L=1, N_PROFILE
                  SURFACE_DOWN_FLUX(L)=SURFACE_DOWN_FLUX(L)
     &               +FLUX_DOWN(L)
               ENDDO
            ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN
               DO L=1, N_PROFILE
                  SURFACE_DOWN_FLUX(L)=SURFACE_DOWN_FLUX(L)
     &               +FLUX_DOWN(L)+PLANCK_AIR_SURFACE(L)
               ENDDO
            ENDIF
         ENDIF
      ENDIF
!
      IF (L_SURF_DOWN_CLR) THEN
         IF (L_NET) THEN
            IF (ISOLIR.EQ.IP_SOLAR) THEN
               DO L=1, N_PROFILE
                  SURF_DOWN_CLR(L)=SURF_DOWN_CLR(L)
     &               +(FLUX_DOWN_CLEAR(L)+FLUX_DIRECT_CLEAR(L)
     &               *(ALBEDO_FIELD_DIR(L)-ALBEDO_FIELD_DIFF(L)))
     &               /(1.0E+00-ALBEDO_FIELD_DIFF(L))
               ENDDO
            ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN
               DO L=1, N_PROFILE
                  SURF_DOWN_CLR(L)=SURF_DOWN_CLR(L)
     &               +(FLUX_DOWN_CLEAR(L)+THERMAL_SOURCE_GROUND(L))
     &               /(1.0E+00-ALBEDO_FIELD_DIFF(L))
               ENDDO
            ENDIF
         ELSE
            IF (ISOLIR.EQ.IP_SOLAR) THEN
               DO L=1, N_PROFILE
                  SURF_DOWN_CLR(L)=SURF_DOWN_CLR(L)
     &               +FLUX_DOWN_CLEAR(L)
               ENDDO
            ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN
               DO L=1, N_PROFILE
                  SURF_DOWN_CLR(L)=SURF_DOWN_CLR(L)
     &               +FLUX_DOWN_CLEAR(L)+PLANCK_AIR_SURFACE(L)
               ENDDO
            ENDIF
         ENDIF
      ENDIF
!
      IF (L_SURF_UP_CLR) THEN
         IF (L_NET) THEN
            IF (ISOLIR.EQ.IP_SOLAR) THEN
               DO L=1, N_PROFILE
                  SURF_UP_CLR(L)=SURF_UP_CLR(L)
     &               +((ALBEDO_FIELD_DIR(L)-ALBEDO_FIELD_DIFF(L))
     &               *FLUX_DIRECT_CLEAR(L)
     &               +ALBEDO_FIELD_DIFF(L)*FLUX_DOWN_CLEAR(L))
     &               /(1.0E+00-ALBEDO_FIELD_DIFF(L))
               ENDDO
            ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN
               DO L=1, N_PROFILE
                  SURF_UP_CLR(L)=SURF_UP_CLR(L)
     &               +(THERMAL_SOURCE_GROUND(L)+ALBEDO_FIELD_DIFF(L)
     &               *FLUX_DOWN_CLEAR(L))/(1.0E+00-ALBEDO_FIELD_DIFF(L))
               ENDDO
            ENDIF
         ELSE
            IF (ISOLIR.EQ.IP_SOLAR) THEN
               DO L=1, N_PROFILE
                  SURF_UP_CLR(L)=SURF_UP_CLR(L)
     &               +FLUX_UP_CLEAR(L)
               ENDDO
            ELSE IF (ISOLIR.EQ.IP_INFRA_RED) THEN
               DO L=1, N_PROFILE
                  SURF_UP_CLR(L)=SURF_UP_CLR(L)
     &               +FLUX_UP_CLEAR(L)+PLANCK_AIR_SURFACE(L)
               ENDDO
            ENDIF
         ENDIF
      ENDIF
!
!     THIS DIAGNOSTIC IS AVAILABLE ONLY IN THE SOLAR REGION.
      IF (L_FLUX_BELOW_690NM_SURF) THEN
         IF (ISOLIR.EQ.IP_SOLAR) THEN
            IF (L_NET) THEN
               DO L=1, N_PROFILE
                  FLUX_BELOW_690NM_SURF(L)=FLUX_BELOW_690NM_SURF(L)
     &               +WEIGHT_690NM*FLUX_DOWN(L)
               ENDDO
            ELSE
               DO L=1, N_PROFILE
                  FLUX_BELOW_690NM_SURF(L)=FLUX_BELOW_690NM_SURF(L)
     &               +WEIGHT_690NM*(FLUX_DOWN(L)-FLUX_UP(L))
               ENDDO
            ENDIF
         ENDIF
      ENDIF
!
!
!
      RETURN
      END
