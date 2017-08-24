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
!+ Subroutine to aggregate clouds into regions.
!
! Method:
!       The clouds in a layer are combined in groups to form regions
!       which will be considered as bulk entities in the solution of the
!       equation of transfer. The extents of these regions are also
!       determined.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       HADAM3          05-06-96                Original Code
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE AGGREGATE_CLOUD(IERR
     &   , N_PROFILE, N_LAYER, N_CLOUD_TOP
     &   , I_CLOUD, I_CLOUD_REPRESENTATION, N_CLOUD_TYPE
     &   , FRAC_CLOUD
     &   , I_REGION_CLOUD, FRAC_REGION
     &   , NPD_PROFILE, NPD_LAYER
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     DUMMY ARRAY SIZES
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
!     MODULE TO SET UNIT NUMBERS FOR STANDARD I/O.
!
      INTEGER
     &     IU_STDIN
!             UNIT NUMBER FOR STANDARD INPUT
     &   , IU_STDOUT
!             UNIT NUMBER FOR STANDARD OUTPUT
     *   , IU_ERR
!             UNIT NUMBER FOR ERROR MESSAGES
!
      PARAMETER(
     &     IU_STDIN=5
     &   , IU_STDOUT=6
     *   , IU_ERR=6
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET ERROR FLAGS IN THE RADIATION CODE.
!
      INTEGER
     &     I_NORMAL
!             ERROR FREE CONDITION
     &   , I_ERR_FATAL
!             FATAL ERROR: IMMEDIATE RETURN
     &   , I_ABORT_CALCULATION
!             CALCULATION ABORTED
     &   , I_MISSING_DATA
!             MISSING DATA ERROR: CONDITIONAL
     &   , I_ERR_IO
!             I/O ERROR
     &   , I_ERR_RANGE
!             INTERPOLATION RANGE ERROR
     &   , I_ERR_EXIST
!             EXISTENCE ERROR
!
      PARAMETER(
     &     I_NORMAL=0
     &   , I_ERR_FATAL=1
     &   , I_ABORT_CALCULATION=2
     &   , I_MISSING_DATA=3
     &   , I_ERR_IO=4
     &   , I_ERR_RANGE=5
     &   , I_ERR_EXIST=6
     &   )
!
!     ------------------------------------------------------------------
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
!     MODULE TO SET REPRESENTATIONS OF CLOUDS.
!
!     REPRESENTATIONS
      INTEGER
     &     IP_CLOUD_HOMOGEN
!             ALL COMPONENTS ARE MIXED HOMOGENEOUSLY
     &   , IP_CLOUD_ICE_WATER
!             ICE AND WATER CLOUDS ARE TREATED SEPARATELY
     &   , IP_CLOUD_CONV_STRAT
!             CLOUDS ARE DIVIDED INTO HOMOGENEOUSLY MIXED
!             STRATIFORM AND CONVECTIVE PARTS
     &   , IP_CLOUD_CSIW
!             CLOUDS DIVIDED INTO ICE AND WATER PHASES AND
!             INTO STRATIFORM AND CONVECTIVE COMPONENTS.
!
      PARAMETER(
     &     IP_CLOUD_HOMOGEN=1
     &   , IP_CLOUD_ICE_WATER=2
     &   , IP_CLOUD_CONV_STRAT=3
     &   , IP_CLOUD_CSIW=4
     &   )
!
!     TYPES OF CLOUDS:
      INTEGER
     &     NP_CLOUD_TYPE(NPD_CLOUD_REPRESENTATION)
!             NUMBER OF TYPE OF CLOUDS IN REPRESENTATION
!
      INTEGER
     &     IP_CLOUD_TYPE_MAP(NPD_CLOUD_COMPONENT
     &       , NPD_CLOUD_REPRESENTATION)
!            MAP OF COMPONENTS CONTRIBUTING TO TYPES OF CLOUDS
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET TYPES OF CLOUDS.
!
      INTEGER
     &     IP_CLOUD_TYPE_HOMOGEN
!             CLOUD COMPOSED OF MIXED WATER AND ICE
     &   , IP_CLOUD_TYPE_WATER
!             CLOUD COMPOSED ONLY OF WATER
     &   , IP_CLOUD_TYPE_ICE
!             CLOUD COMPOSED ONLY OF ICE
     &   , IP_CLOUD_TYPE_STRAT
!             MIXED-PHASE STRATIFORM CLOUD
     &   , IP_CLOUD_TYPE_CONV
!             MIXED-PHASE CONVECTIVE CLOUD
     &   , IP_CLOUD_TYPE_SW
!             STRATIFORM WATER CLOUD
     &   , IP_CLOUD_TYPE_SI
!             STRATIFORM ICE CLOUD
     &   , IP_CLOUD_TYPE_CW
!             CONVECTIVE WATER CLOUD
     &   , IP_CLOUD_TYPE_CI
!             CONVECTIVE ICE CLOUD
!
      PARAMETER(
     &     IP_CLOUD_TYPE_HOMOGEN=1
     &   , IP_CLOUD_TYPE_WATER=1
     &   , IP_CLOUD_TYPE_ICE=2
     &   , IP_CLOUD_TYPE_STRAT=1
     &   , IP_CLOUD_TYPE_CONV=2
     &   , IP_CLOUD_TYPE_SW=1
     &   , IP_CLOUD_TYPE_SI=2
     &   , IP_CLOUD_TYPE_CW=3
     &   , IP_CLOUD_TYPE_CI=4
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO DEFINE REFERENCE NUMBERS FOR CLOUD SCHEMES.
!
      INTEGER
     &     IP_CLOUD_MIX_MAX
!             MAXIMUM/RANDOM OVERLAP IN A MIXED COLUMN
     &   , IP_CLOUD_MIX_RANDOM
!             RANDOM OVERLAP IN A MIXED COLUMN
     &   , IP_CLOUD_COLUMN_MAX
!             MAXIMUM OVERLAP IN A COLUMN MODEL
     &   , IP_CLOUD_CLEAR
!             CLEAR COLUMN
     &   , IP_CLOUD_TRIPLE
!             MIXED COLUMN WITH SPLIT BETWEEN 
!             CONVECTIVE AND LAYER CLOUD.
!
      PARAMETER(
     &     IP_CLOUD_MIX_MAX=2
     &   , IP_CLOUD_MIX_RANDOM=4
     &   , IP_CLOUD_COLUMN_MAX=3
     &   , IP_CLOUD_CLEAR=5
     &   , IP_CLOUD_TRIPLE=6
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO DEFINE REFERENCE NUMBERS FOR REGIONS OF CLOUDS.
!
      INTEGER
     &     IP_REGION_CLEAR
!             REFERENCE NUMBER FOR CLEAR-SKY REGION
     &   , IP_REGION_STRAT
!             REFERENCE NUMBER FOR STRATIFORM CLOUDY REGION
     &   , IP_REGION_CONV
!             REFERENCE NUMBER FOR CONVECTIVE CLOUDY REGION
!
      PARAMETER(
     &     IP_REGION_CLEAR=1
     &   , IP_REGION_STRAT=2
     &   , IP_REGION_CONV=3
     &   )
!
!     ------------------------------------------------------------------
!
!     DUMMY VARIABLES.
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , N_CLOUD_TOP
!             TOPMOST CLOUDY LAYER
      INTEGER   !, INTENT(IN)
     &     I_CLOUD
!             CLOUD SCHEME USED
     &   , I_CLOUD_REPRESENTATION
!             REPRESENTATION OF CLOUDS USED
     &   , N_CLOUD_TYPE
!             NUMBER OF TYPES OF CLOUD
!
      REAL      !, INTENT(OUT)
     &     FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             FRACTIONS OF EACH TYPE OF CLOUD
!
      INTEGER   !, INTENT(OUT)
     &     I_REGION_CLOUD(NPD_CLOUD_TYPE)
!             REGIONS IN WHICH PARTICULAR TYPES OF CLOUD FALL
      REAL      !, INTENT(OUT)
     &     FRAC_REGION(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             FRACTIONS OF TOTAL CLOUD OCCUPIED BY EACH REGION
!
!
!     LOCAL VARIABLES
      INTEGER
     &     I
!            LOOP VARIABLE
     &   , L
!            LOOP VARIABLE
     &   , K
!            LOOP VARIABLE
!
!
!
      IF (I_CLOUD.EQ.IP_CLOUD_TRIPLE) THEN
!
         IF (I_CLOUD_REPRESENTATION.EQ.IP_CLOUD_CSIW) THEN
!
            DO K=1, N_CLOUD_TYPE
               IF (K.EQ.IP_CLOUD_TYPE_SW) THEN
                  I_REGION_CLOUD(K)=IP_REGION_STRAT
               ELSE IF (K.EQ.IP_CLOUD_TYPE_SI) THEN
                  I_REGION_CLOUD(K)=IP_REGION_STRAT
               ELSE IF (K.EQ.IP_CLOUD_TYPE_CW) THEN
                  I_REGION_CLOUD(K)=IP_REGION_CONV
               ELSE IF (K.EQ.IP_CLOUD_TYPE_CI) THEN
                  I_REGION_CLOUD(K)=IP_REGION_CONV
               ENDIF
            ENDDO
!
            DO I=N_CLOUD_TOP, N_LAYER
               DO L=1, N_PROFILE
                  FRAC_REGION(L, I, IP_REGION_STRAT)
     &               =FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SW)
     &               +FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SI)
                  FRAC_REGION(L, I, IP_REGION_CONV)
     &               =FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CW)
     &               +FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CI)
               ENDDO
            ENDDO
!
         ELSE IF (I_CLOUD_REPRESENTATION.EQ.IP_CLOUD_CONV_STRAT) THEN
!
            DO K=1, N_CLOUD_TYPE
               IF (K.EQ.IP_CLOUD_TYPE_STRAT) THEN
                  I_REGION_CLOUD(K)=IP_REGION_STRAT
               ELSE IF (K.EQ.IP_CLOUD_TYPE_CONV) THEN
                  I_REGION_CLOUD(K)=IP_REGION_CONV
               ENDIF
            ENDDO
!
            DO I=N_CLOUD_TOP, N_LAYER
               DO L=1, N_PROFILE
                  FRAC_REGION(L, I, IP_REGION_STRAT)
     &               =FRAC_CLOUD(L, I, IP_CLOUD_TYPE_STRAT)
                  FRAC_REGION(L, I, IP_REGION_CONV)
     &               =FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CONV)
               ENDDO
            ENDDO
!
!
         ELSE
            WRITE(IU_ERR, '(/A)')
     &         '*** ERROR: THIS REPRESENTATION OF CLOUDS IS NOT '
     &         //'COMPATIBLE WITH THE TRIPLE OVERLAP.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
!
      ENDIF
!
!
!
      RETURN
      END
