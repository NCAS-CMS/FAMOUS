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
!+ Subroutine to find single scattering propeties of all regions.
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
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SINGLE_SCATTERING_ALL(I_SCATTER_METHOD_BAND
!                       Atmospheric Propeties
     &   , N_PROFILE, N_LAYER, D_MASS
!                       Cloudy Propeties
     &   , L_CLOUD, N_CLOUD_TOP, N_CLOUD_TYPE
!                       Optical Propeties
     &   , K_GREY_TOT_FREE, K_EXT_SCAT_FREE
     &   , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD
     &   , K_GAS_ABS
!                       Single Scattering Propeties
     &   , TAU_FREE, OMEGA_FREE
     &   , TAU_CLOUD, OMEGA_CLOUD
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER
     &   )
!
!
      IMPLICIT NONE
!
!
!     SIZES OF ARRAYS.
      INTEGER   !, INTENT(IN)
     &     NPD_PROFILE
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     INCLUDE COMDECKS.
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
!
!
!     DUMMY VARIABLES.
      INTEGER   !, INTENT(IN)
     &     I_SCATTER_METHOD_BAND
!             TREATMENT OF SCATTERING IN THE BAND
!
!                       Atmospheric Properties
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
      REAL      !, INTENT(IN)
     &     D_MASS(NPD_PROFILE, NPD_LAYER)
!             MASS THICKNESS OF EACH LAYER
!
!                       Cldoudy Propeties
      LOGICAL   !, INTENT(IN)
     &     L_CLOUD
!             FLAG FOR CLOUDS
      INTEGER   !, INTENT(IN)
     &     N_CLOUD_TOP
!             TOPMOST CLOUDY LAYER
     &   , N_CLOUD_TYPE
!             NUMBER OF TYPES OF CLOUDS
!
!                       Optical Properties
      REAL      !, INTENT(IN)
     &     K_GREY_TOT_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE ABSORPTIVE EXTINCTION
     &   , K_EXT_SCAT_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE SCATTERING EXTINCTION
     &   , K_GREY_TOT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             STRATIFORM ABSORPTIVE EXTINCTION
     &   , K_EXT_SCAT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             STRATIFORM SCATTERING EXTINCTION
     &   , K_GAS_ABS(NPD_PROFILE, NPD_LAYER)
!             GASEOUS EXTINCTION
!
!                       Single Scattering Properties
      REAL      !, INTENT(OUT)
     &     TAU_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE OPTICAL DEPTH
     &   , OMEGA_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE SINGLE SCATTERING ALBEDO
     &   , TAU_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY OPTICAL DEPTH
     &   , OMEGA_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY SINGLE SCATTERING ALBEDO
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     K
!             LOOP VARIABLE
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     SINGLE_SCATTERING
!
!
!
!     CLEAR-SKY PROPERTIES:
!
      CALL SINGLE_SCATTERING(I_SCATTER_METHOD_BAND
     &   , N_PROFILE, N_LAYER, 1
     &   , D_MASS
     &   , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, K_GAS_ABS
     &   , TAU_FREE, OMEGA_FREE
     &   , NPD_PROFILE, NPD_LAYER
     &   )
!
      IF (L_CLOUD) THEN
         DO K=1, N_CLOUD_TYPE
            CALL SINGLE_SCATTERING(I_SCATTER_METHOD_BAND
     &         , N_PROFILE, N_LAYER, N_CLOUD_TOP
     &         , D_MASS
     &         , K_GREY_TOT_CLOUD(1, 1, K)
     &         , K_EXT_SCAT_CLOUD(1, 1, K)
     &         , K_GAS_ABS
     &         , TAU_CLOUD(1, 1, K), OMEGA_CLOUD(1, 1, K)
     &         , NPD_PROFILE, NPD_LAYER
     &         )
         ENDDO
      ENDIF
!
!
!
      RETURN
      END
