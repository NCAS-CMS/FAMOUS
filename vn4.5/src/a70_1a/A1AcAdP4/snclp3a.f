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
!+ Function to set number of cloudy parameters.
!
! Method:
!       Straightforward
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.5             18-05-98                Code for new parametr-
!                                               ization of droplets
!                                               included.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77 with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      FUNCTION SET_N_CLOUD_PARAMETER(I_SCHEME, I_COMPONENT
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
!     MODULE TO SET COMPONENTS OF CLOUDS.
!
      INTEGER
     &     IP_CLCMP_ST_WATER
!             STRATIFORM WATER DROPLETS
     &   , IP_CLCMP_ST_ICE
!             STRATIFORM ICE CRYSTALS
     &   , IP_CLCMP_CNV_WATER
!             CONVECTIVE WATER DROPLETS
     &   , IP_CLCMP_CNV_ICE
!             CONVECTIVE ICE CRYSTALS
!
      PARAMETER(
     &     IP_CLCMP_ST_WATER=1
     &   , IP_CLCMP_ST_ICE=2
     &   , IP_CLCMP_CNV_WATER=3
     &   , IP_CLCMP_CNV_ICE=4
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET NUMBERS FOR WATER CLOUD SCHEMES.
!
      INTEGER
     &     NPD_CLOUD_FIT
!             NUMBER OF CLOUD FITTING SCHEMES
     &   , IP_SLINGO_SCHRECKER
!             PARAMETRIZATION OF SLINGO-SCHRECKER
     &   , IP_ACKERMAN_STEPHENS
!             PARAMETRIZATION OF ACKERMAN & STEPHENS
     &   , IP_DROP_UNPARAMETRIZED
!             UNPARAMETRIZED DROPLET DATA
     &   , IP_DROP_PADE_2
!             PADE APPROXIMATION OF THE SECOND ORDER
!             (THIRD ORDER FOR THE EXTINCTION)
!
!
      PARAMETER(
     &     NPD_CLOUD_FIT=3
     &   , IP_SLINGO_SCHRECKER=1
     &   , IP_ACKERMAN_STEPHENS=2
     &   , IP_DROP_UNPARAMETRIZED=3
     &   , IP_DROP_PADE_2=5
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
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     I_SCHEME
!             PARAMETRIZATION SCHEME
     &   , I_COMPONENT
!             COMPONENT IN CLOUD
!
      INTEGER   !, INTENT(OUT)
     &     SET_N_CLOUD_PARAMETER
!             RETURNED NUMBER OF COEFFICIENTS IN PARAMETRIZATION
!
!
!
      IF ( (I_COMPONENT.EQ.IP_CLCMP_ST_WATER).OR.
     &     (I_COMPONENT.EQ.IP_CLCMP_CNV_WATER) ) THEN
!
         IF (I_SCHEME.EQ.IP_SLINGO_SCHRECKER) THEN
           SET_N_CLOUD_PARAMETER =6
         ELSE IF (I_SCHEME.EQ.IP_ACKERMAN_STEPHENS) THEN
            SET_N_CLOUD_PARAMETER=9
         ELSE IF (I_SCHEME.EQ.IP_DROP_PADE_2) THEN
            SET_N_CLOUD_PARAMETER=16
         ENDIF
!
      ELSE IF ( (I_COMPONENT.EQ.IP_CLCMP_ST_ICE).OR.
     &          (I_COMPONENT.EQ.IP_CLCMP_CNV_ICE) ) THEN
!
         IF (I_SCHEME.EQ.IP_SLINGO_SCHRECKER_ICE) THEN
            SET_N_CLOUD_PARAMETER=6
         ELSE IF (I_SCHEME.EQ.IP_ICE_ADT) THEN
            SET_N_CLOUD_PARAMETER=30
         ELSE IF (I_SCHEME.EQ.IP_SUN_SHINE_VN2_VIS) THEN
            SET_N_CLOUD_PARAMETER=6
         ELSE IF (I_SCHEME.EQ.IP_SUN_SHINE_VN2_IR) THEN
            SET_N_CLOUD_PARAMETER=0
         ENDIF
!
      ENDIF
!
!
!
      RETURN
      END
