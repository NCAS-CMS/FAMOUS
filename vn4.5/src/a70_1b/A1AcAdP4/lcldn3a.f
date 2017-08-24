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
!+ Function to determine whether densities are required for clouds.
!
! Method:
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                L_CLOUD_DENSITY set
!                                               as .FALSE. initially
!                                               to provide a default.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77 with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      FUNCTION L_CLOUD_DENSITY(N_CONDENSED, I_PHASE_CMP, L_CLOUD_CMP
     &   , I_CONDENSED_PARAM
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     INCLUDE COMDECKS
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
!     MODULE TO SET NUMBERS FOR ICE CLOUD SCHEMES.
!
      INTEGER
     &     NPD_ICE_CLOUD_FIT
!             NUMBER OF CLOUD FITTING SCHEMES
     &   , IP_SLINGO_SCHRECKER_ICE
!             PARAMETRIZATION OF SLINGO AND SCHRECKER.
     &   , IP_ICE_UNPARAMETRIZED
!             UNPARAMETRIZED ICE CRYSTAL DATA
     &   , IP_SUN_SHINE_VN2_VIS
!             SUN AND SHINE'S PARAMETRIZATION IN THE VISIBLE (VERSION 2)
     &   , IP_SUN_SHINE_VN2_IR
!             SUN AND SHINE'S PARAMETRIZATION IN THE IR (VERSION 2)
     &   , IP_ICE_ADT
!             SCHEME BASED ON ANOMALOUS DIFFRACTION THEORY
!             FOR ICE CRYSTALS

!
      PARAMETER(
     &     NPD_ICE_CLOUD_FIT=6
     &   , IP_SLINGO_SCHRECKER_ICE=1
     &   , IP_ICE_UNPARAMETRIZED=3
     &   , IP_SUN_SHINE_VN2_VIS=4
     &   , IP_SUN_SHINE_VN2_IR=5
     &   , IP_ICE_ADT=6
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET INDICES FOR PHASES.
!
      INTEGER
     &     IP_PHASE_WATER
!             LIQUID PHASE
     &   , IP_PHASE_ICE
!             ICE PHASE
!
      PARAMETER(
     &     IP_PHASE_WATER=1
     &   , IP_PHASE_ICE=2
     &   )
!
!     ------------------------------------------------------------------
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_CONDENSED
!             NUMBER OF TYPES OF CONDENSATE
     &   , I_PHASE_CMP(NPD_CLOUD_COMPONENT)
!             PHASES OF COMPONENTS
     &   , I_CONDENSED_PARAM(NPD_CLOUD_COMPONENT)
!             PARAMETRIZATIONS OF COMPONENTS
      LOGICAL   !, INTENT(IN)
     &     L_CLOUD_CMP(NPD_CLOUD_COMPONENT)
!             FLAGS FOR ENABLED COMPONENTS
      LOGICAL   !, INTENT(OUT)
     &     L_CLOUD_DENSITY
!             RETURNED FLAG FOR CALCULATING DENSITY
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     K
!             LOOP VARIABLE
!
!
      L_CLOUD_DENSITY=.FALSE.
!
!     DENSITIES MUST BE CALCULATED IF SUN & SHINE'S PARAMETRIZATIONS
!     ARE USED.
      DO K=1, N_CONDENSED
         L_CLOUD_DENSITY=L_CLOUD_DENSITY.OR.
     &      (L_CLOUD_CMP(K).AND.(I_PHASE_CMP(K).EQ.IP_PHASE_ICE).AND.
     &      ( (I_CONDENSED_PARAM(K).EQ.IP_SUN_SHINE_VN2_VIS).OR.
     &        (I_CONDENSED_PARAM(K).EQ.IP_SUN_SHINE_VN2_IR) ) )
      ENDDO
!
!
!
      RETURN
      END
