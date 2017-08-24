C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!+ Subroutine to calculate an observed effective radius.
!
! Purpose:
!   An effective radius as observed from above the cloud-top is
!   calculated.
!
! Method:
!   For each type of cloud containing water in any layer the effective
!   radius is weighted with the product of the area of the cloud and the
!   probability that light emitted from the cloud reaches the observing
!   instrument.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.4             24-01-97                Original Code
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_CLOUD_LEVEL_DIAG(IERR, N_PROFILE, NLEVS, NCLDS
     &   , I_GATHER
     &   , I_CLOUD, I_CLOUD_REPRESENTATION
     &   , W_CLOUD, FRAC_CLOUD
     &   , CONDENSED_MIX_RATIO, CONDENSED_RE
     &   , L_OBSERVED_RE, WEIGHTED_RE, SUM_WEIGHT_RE
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED.
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
!
!
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!     DIMENSIONS OF ARRAYS:
      INTEGER   !, INTENT(IN)
     &     NPD_FIELD
!             SIZE OF ARRAY OF ARRAYS PASSED FROM MAIN CODE
     &   , NPD_PROFILE
!             SIZE OF ARRAY OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     ACTUAL SIZES USED:
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , NLEVS
!             NUMBER OF ATMOSPHERIC LAYERS
     &   , NCLDS
!             NUMBER OF CLOUDY LEVELS
     &   , I_GATHER(NPD_FIELD)
!             LIST OF GATHERED POINTS
!
!     LOGICAL FLAGS FOR DIAGNOSTICS
      LOGICAL   !, INTENT(IN)
     &     L_OBSERVED_RE
!             FLAG TO ENABLE DIAGNOSIS OF EFFECTIVE RADIUS SEEN FROM
!             SPACE (N.B. THE ROUTINE IS AT PRESENT CALLED ONLY IF
!             THIS IS TRUE, BUT ITS PRESENCE HERE ALLOWS FOR POSSIBLE
!             FUTURE EXTENSION OF THE ROUTINE).
!
!     REPRESENTATION OF CLOUDS
      INTEGER   !, INTENT(IN)
     &     I_CLOUD_REPRESENTATION
!             REPRESENTATION OF CLOUDS
     &   , I_CLOUD
!             TREATMENT OF OVERLAPS
!
      REAL      !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             TOTAL AMOUNTS OF CLOUD
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             FRACTION OF TYPES OF CLOUD
     &   , CONDENSED_RE(NPD_PROFILE, 0: NPD_LAYER, NPD_CLOUD_COMPONENT)
!             EFFECTIVE RADII OF CLOUDY COMPONENTS
     &   , CONDENSED_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER
     &      , NPD_CLOUD_COMPONENT)
!             MASS MIXING RATIOS OF CONDENSED COMPONENTS
!
      REAL      !, INTENT(OUT)
     &     WEIGHTED_RE(NPD_FIELD)
!             WEIGHTED SUM OF EFFECTIVE RADIUS AND WEIGHTING FUNCTION
     &   , SUM_WEIGHT_RE(NPD_FIELD)
!             SUM OF WEIGHTS FOR EFFECTIVE RADIUS
!
!
!
!     LOCAL VARIABLES:
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
     &   , I_INV
!             INVERTED LOOP INDEX
      REAL
     &     TRANS_OVERLYING_SPACE(NPD_PROFILE)
!             PROBABILITY OF A PHOTON IN CLEAR AIR IN THE LEVEL ABOVE
!             THE CURRENT ONE REACHING SPACE
     &   , AREA_EXPOSED(NPD_PROFILE)
!             TOTAL AREA OF CLOUD IN THE CURRENT LAYER EXPOSED TO
!             CLEAR AIR IN THE LAYER ABOVE
     &   , AREA_EXPOSED_ST(NPD_PROFILE)
!             TOTAL AREA OF STRATIFORM CLOUD IN THE CURRENT LAYER
!             EXPOSED TO CLEAR AIR IN THE LAYER ABOVE
     &   , AREA_EXPOSED_CNV(NPD_PROFILE)
!             TOTAL AREA OF CONVECTIVE CLOUD IN THE CURRENT LAYER
!             EXPOSED TO CLEAR AIR IN THE LAYER ABOVE
     &   , AREA_CLEAR_ABOVE(NPD_PROFILE)
!             AREA OF THE CLEAR SKY REGION IN THE LAYER ABOVE
     &   , AREA_STRAT(NPD_PROFILE)
!             AREA OF STRATIFORM CLOUD IN THE CURRENT LAYER
     &   , AREA_STRAT_ABOVE(NPD_PROFILE)
!             AREA OF STRATIFORM CLOUD IN THE LAYER ABOVE
     &   , AREA_CONV(NPD_PROFILE)
!             AREA OF CONVECTIVE CLOUD IN THE CURRENT LAYER
     &   , AREA_CONV_ABOVE(NPD_PROFILE)
!             AREA OF CONVECTIVE CLOUD IN THE LAYER ABOVE
     &   , AREA_CLEAR_CLEAR(NPD_PROFILE)
!             AREA OF BOUNDARY WHERE CLEAR SKY OVERLIES CLEAR SKY
     &   , AREA_CLEAR(NPD_PROFILE)
!             AREA OF CLEAR SKY IN THE CURRENT LAYER
!             DOWN TO A LEVEL
     &   , AREA_UNCORRELATED(NPD_PROFILE)
!             UNCORRELATED REGION ON THE INTERFACE
     &   , WEIGHTED_RE_G(NPD_PROFILE)
!             WEIGHTED SUM OF EFFECTIVE RADIUS AND WEIGHTING FUNCTION
     &   , SUM_WEIGHT_RE_G(NPD_PROFILE)
!             SUM OF WEIGHTS FOR EFFECTIVE RADIUS
!
!     VARIABLES FOR GATHERING
      INTEGER
     &     N_LIST
!             NUMBER OF POINTS IN LIST
     &   , L_LIST(NPD_PROFILE)
!             INDICES OF POINTS IN LIST
!
!     INDICATOR FUNCTION
      REAL
     &     CHI_CNV(NPD_PROFILE)
!             CONVECTIVE INDICATOR FUNCTION
     &   , CHI_ST(NPD_PROFILE)
!             STRATIFORM INDICATOR FUNCTION
!
!
!
!
!     INITIALIZATION OF DIAGNOSTIC FIELDS.
!
      IF (L_OBSERVED_RE) THEN
         CALL R2_ZERO_1D(NPD_PROFILE, WEIGHTED_RE)
         CALL R2_ZERO_1D(NPD_PROFILE, SUM_WEIGHT_RE)
         DO L=1, N_PROFILE
            WEIGHTED_RE_G(L)=0.0E+00
            SUM_WEIGHT_RE_G(L)=0.0E+00
         ENDDO
      ENDIF
!
!     INITIALIZE THE TRANSMISION ABOVE CLOUDS.
      DO L=1, N_PROFILE
         TRANS_OVERLYING_SPACE(L)=1.0E+00
         AREA_CLEAR_ABOVE(L)=1.0E+00
      ENDDO
      IF (L_OBSERVED_RE.AND.(I_CLOUD.EQ.IP_CLOUD_TRIPLE)) THEN
         DO L=1, N_PROFILE
            AREA_STRAT_ABOVE(L)=0.0E+00
            AREA_CONV_ABOVE(L)=0.0E+00
         ENDDO
      ENDIF
!
!     STEP DOWN THROUGH THE ATMOSPHERE CALCULATING CONTRIBUTIONS TO
!     THE DIAGNOSTICS AND SUBSEQUENTLY ALLOWING FOR TRANSMISSION
!     THROUGH THE CURRENT LAYER.
!
      DO I=NCLDS, 1, -1
         I_INV=NLEVS+1-I
!
         DO L=1, N_PROFILE
            AREA_CLEAR(L)=1.0E+00-W_CLOUD(L, I_INV)
         ENDDO
!
!        CALCULATE THE LOCAL AREA OF CLOUD RADIATING INTO CLEAR AIR.
         IF (I_CLOUD.EQ.IP_CLOUD_MIX_RANDOM) THEN
            DO L=1, N_PROFILE
               AREA_EXPOSED(L)=W_CLOUD(L, I_INV)
     &            *AREA_CLEAR_ABOVE(L)
            ENDDO
         ELSE IF ( (I_CLOUD.EQ.IP_CLOUD_MIX_MAX).OR.
     &             (I_CLOUD.EQ.IP_CLOUD_TRIPLE) ) THEN
            DO L=1, N_PROFILE
               AREA_EXPOSED(L)=MAX(0.0E+00, (W_CLOUD(L, I_INV)
     &            +AREA_CLEAR_ABOVE(L)-1.0E+00))
            ENDDO
         ENDIF
!
!
!
!
         IF (L_OBSERVED_RE) THEN
!
!
!
            IF (I_CLOUD_REPRESENTATION.EQ.IP_CLOUD_CONV_STRAT) THEN
!
!
               IF ( (I_CLOUD.EQ.IP_CLOUD_MIX_MAX).OR.
     &              (I_CLOUD.EQ.IP_CLOUD_MIX_RANDOM) ) THEN
!
!                 IF THE OVERLAP OF CONVECTIVE CLOUD IS NOT ASSUMED
!                 TO BE COHERENT THE OVERALL EXPOSED AREA MAY BE
!                 PARTITIONED ACCORDING TO THE FRACTIONAL
!                 CONTRIBUTIONS OF CLOUD IN THE CURRENT LAYER.
!
                  DO L=1, N_PROFILE
                     AREA_EXPOSED_ST(L)=AREA_EXPOSED(L)
     &                  *FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_STRAT)
                     AREA_EXPOSED_CNV(L)=AREA_EXPOSED(L)
     &                  *FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_CONV)
                  ENDDO
!
               ELSE IF (I_CLOUD.EQ.IP_CLOUD_TRIPLE) THEN
!
!                 HERE, THE DIFFERENT TYPES OF CLOUDS OVERLAP
!                 COHERENTLY SO STRATIFORM CLOUD WILL BE EXPOSED
!                 ONLY IF THERE IS LESS STRATIFORM CLOUD IN THE
!                 LAYER ABOVE AND MORE CLEAR AIR IN THE LAYER ABOVE:
!                 UNDER THESE CONDITIONS THE NON-CORRELATED AREAS
!                 OVERLAP RANDOMLY.
!
                  DO L=1, N_PROFILE
                     AREA_STRAT(L)=W_CLOUD(L, I_INV)
     &                  *FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_STRAT)
                     AREA_EXPOSED_ST(L)=MAX(0.0E+00
     &                  , (AREA_STRAT(L)-AREA_STRAT_ABOVE(L)))
                     AREA_EXPOSED_ST(L)=MAX(0.0E+00, AREA_EXPOSED_ST(L)
     &                  *(AREA_CLEAR_ABOVE(L)-AREA_CLEAR(L)))
                     AREA_EXPOSED_CNV(L)
     &                  =AREA_EXPOSED(L)-AREA_EXPOSED_ST(L)
                  ENDDO
               ELSE
                  WRITE(IU_ERR, '(/A)')
     &               '*** ERROR: THE DIAGNOSTIC OF OBSERVED RE HAS NOT '
     &               //'BEEN IMPLEMENTED WITH THIS OVERLAP OPTION.'
                  IERR=I_ERR_FATAL
                  RETURN
               ENDIF
!
!              THE INDICATOR FUNCTIONS FOR LIQUID WATER IN
!              CONVECTIVE OR STRAIFORM CLOUDS ARE SET TO 1
!              IF THERE IS ANY LIQUID WATER AND TO 0 OTHERWISE.
               DO L=1, N_PROFILE
                  IF (CONDENSED_MIX_RATIO(1, I_INV, IP_CLCMP_CNV_WATER)
     &               .GT.0.0E+00) THEN
                     CHI_CNV(L)=1.0E+00
                  ELSE
                     CHI_CNV(L)=0.0E+00
                  ENDIF
                  IF (CONDENSED_MIX_RATIO(1, I_INV, IP_CLCMP_ST_WATER)
     &               .GT.0.0E+00) THEN
                     CHI_ST(L)=1.0E+00
                  ELSE
                     CHI_ST(L)=0.0E+00
                  ENDIF
                  CHI_ST(L)=0.0E+00
               ENDDO
!
!              INCLUDE CONTRIBUTIONS FROM CONVECTIVE AND STRATIFORM
!              WATER CLOUDS.
               DO L=1, N_PROFILE
                  WEIGHTED_RE_G(L)=WEIGHTED_RE_G(L)
     &               +TRANS_OVERLYING_SPACE(L)
     &               *(AREA_EXPOSED_CNV(L)*CHI_CNV(L)
     &               *CONDENSED_RE(L, I_INV, IP_CLCMP_CNV_WATER)
     &               +AREA_EXPOSED_ST(L)*CHI_ST(L)
     &               *CONDENSED_RE(L, I_INV, IP_CLCMP_ST_WATER))
                  SUM_WEIGHT_RE_G(L)=SUM_WEIGHT_RE_G(L)
     &               +TRANS_OVERLYING_SPACE(L)
     &               *(AREA_EXPOSED_CNV(L)*CHI_CNV(L)
     &               +AREA_EXPOSED_ST(L)*CHI_ST(L))
               ENDDO
!
            ELSE IF (I_CLOUD_REPRESENTATION.EQ.IP_CLOUD_CSIW) THEN
!
               IF ( (I_CLOUD.EQ.IP_CLOUD_MIX_MAX).OR.
     &              (I_CLOUD.EQ.IP_CLOUD_MIX_RANDOM) ) THEN
!
!                 IF THE OVERLAP OF CONVECTIVE CLOUD IS NOT ASSUMED
!                 TO BE COHERENT THE OVERALL EXPOSED AREA MAY BE
!                 PARTITIONED ACCORDING TO THE FRACTIONAL
!                 CONTRIBUTIONS OF CLOUD IN THE CURRENT LAYER.
!                 THE EXPOSED AREAS INCLUDE ONLY THE PARTS OF THE
!                 CLOUDS CONTAINING WATER DROPLETS.
!
                  DO L=1, N_PROFILE
                     AREA_EXPOSED_ST(L)=AREA_EXPOSED(L)
     &                  *FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_SW)
                     AREA_EXPOSED_CNV(L)=AREA_EXPOSED(L)
     &                  *FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_CW)
                  ENDDO
!
               ELSE IF (I_CLOUD.EQ.IP_CLOUD_TRIPLE) THEN
!
!                 HERE, THE DIFFERENT TYPES OF CLOUDS OVERLAP
!                 COHERENTLY SO STRATIFORM CLOUD WILL BE EXPOSED
!                 ONLY IF THERE IS LESS STRATIFORM CLOUD IN THE
!                 LAYER ABOVE AND MORE CLEAR AIR IN THE LAYER ABOVE:
!                 UNDER THESE CONDITIONS THE NON-CORRELATED AREAS
!                 OVERLAP RANDOMLY.
!                 THE ACTUAL EXPOSED AREAS OF CONVECTIVE OR
!                 STRATIFORM CLOUD MUST THEN BE WEIGHTED BY FACTORS
!                 REPRESENTING THE LIQUID PORTION OF EACH CLOUD, SINCE
!                 NOTHING IS RETRIEVED OVER ICE. (THE HORIZONTAL
!                 ARRANGEMENT OF ICE AND WATER WITHIN EITHER TYPE OF
!                 CLOUD IS RANDOM).
!
                  DO L=1, N_PROFILE
!
                     AREA_STRAT(L)=W_CLOUD(L, I_INV)
     &                  *(FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_SW)
     &                  +FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_SI))
                     AREA_CONV(L)=W_CLOUD(L, I_INV)
     &                  *(FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_CW)
     &                  +FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_CI))
                     AREA_UNCORRELATED(L)=1.0E+00
     &                  -MIN(AREA_CLEAR(L), AREA_CLEAR_ABOVE(L))
     &                  -MIN(AREA_STRAT(L), AREA_STRAT_ABOVE(L))
     &                  -MIN(AREA_CONV(L), AREA_CONV_ABOVE(L))
                     AREA_EXPOSED_ST(L)=MAX(0.0E+00
     &                  , (AREA_STRAT(L)-AREA_STRAT_ABOVE(L)))
                     IF (AREA_UNCORRELATED(L).GT.0.0E+00) THEN
                        AREA_EXPOSED_ST(L)
     &                     =MAX(0.0E+00, AREA_EXPOSED_ST(L)
     &                     *(AREA_CLEAR_ABOVE(L)-AREA_CLEAR(L)))
     &                     /AREA_UNCORRELATED(L)
                     ELSE
                        AREA_EXPOSED_ST(L)=0.0E+00
                     ENDIF
                     AREA_EXPOSED_CNV(L)
     &                  =AREA_EXPOSED(L)-AREA_EXPOSED_ST(L)
!
                     IF (FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_CW)
     &                  .GT.0.0E+00) THEN
                        AREA_EXPOSED_CNV(L)=AREA_EXPOSED_CNV(L)
     &                     /(1.0E+00
     &                     +FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_CI)
     &                     /FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_CW))
                     ELSE
                        AREA_EXPOSED_CNV(L)=0.0E+00
                     ENDIF
!
                     IF (FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_SW)
     &                  .GT.0.0E+00) THEN
                        AREA_EXPOSED_ST(L)=AREA_EXPOSED_ST(L)
     &                     /(1.0E+00
     &                     +FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_SI)
     &                     /FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_SW))
                     ELSE
                        AREA_EXPOSED_ST(L)=0.0E+00
                     ENDIF
!
                  ENDDO
               ELSE
                  WRITE(IU_ERR, '(/A)')
     &               '*** ERROR: THE DIAGNOSTIC OF OBSERVED RE HAS NOT '
     &               //'BEEN IMPLEMENTED WITH THIS OVERLAP OPTION.'
                  IERR=I_ERR_FATAL
                  RETURN
               ENDIF
!
!
               DO L=1, N_PROFILE
!
                  WEIGHTED_RE_G(L)=WEIGHTED_RE_G(L)
     &               +TRANS_OVERLYING_SPACE(L)
     &               *(AREA_EXPOSED_CNV(L)
     &               *CONDENSED_RE(L, I_INV, IP_CLCMP_CNV_WATER)
     &               +AREA_EXPOSED_ST(L)
     &               *CONDENSED_RE(L, I_INV, IP_CLCMP_ST_WATER))
                  SUM_WEIGHT_RE_G(L)=SUM_WEIGHT_RE_G(L)
     &               +TRANS_OVERLYING_SPACE(L)
     &               *(AREA_EXPOSED_CNV(L)+AREA_EXPOSED_ST(L))
               ENDDO
!
            ENDIF
!
!
         ENDIF
!
!
!
!        ADVANCE THE STORED QUANTITIES REFFERRING TO OVERLYING LAYERS.
!
!
!        THE TRANSMISSION TO SPACE CURRENTLY HOLDS THE PROBABILITY THAT
!        A PHOTON TRAVELLING UPWARDS IN THE CLEAR AIR IN THE LAYER ABOVE
!        WILL ESCAPE TO SPACE WITHOUT ENCOUNTERING A CLOUD. TO ADVANCE
!        THIS TO THE CURRENT LAYER IT MUST BE MULTIPLIED BY A FACTOR
!        REPRESENTING THE OVERLAP ASSUMPTION AT THE TOP OF THE PRESENT
!        LAYER.
!
         IF (I_CLOUD.EQ.IP_CLOUD_MIX_RANDOM) THEN
!
            DO L=1, N_PROFILE
               TRANS_OVERLYING_SPACE(L)=TRANS_OVERLYING_SPACE(L)
     &            *AREA_CLEAR_ABOVE(L)
            ENDDO
!
         ELSE IF ( (I_CLOUD.EQ.IP_CLOUD_MIX_MAX).OR.
     &             (I_CLOUD.EQ.IP_CLOUD_TRIPLE) ) THEN
!
            DO L=1, N_PROFILE
               AREA_CLEAR_CLEAR(L)=MIN(AREA_CLEAR(L)
     &            , AREA_CLEAR_ABOVE(L))
               IF (AREA_CLEAR(L).GT.0.0E+00) THEN
                  TRANS_OVERLYING_SPACE(L)=TRANS_OVERLYING_SPACE(L)
     &               *AREA_CLEAR_CLEAR(L)/AREA_CLEAR(L)
               ELSE
                  TRANS_OVERLYING_SPACE(L)=0.0E+00
               ENDIF
            ENDDO
!
         ENDIF
!
!        ADVANCE THE AREAS OF CLOUD.
         DO L=1, N_PROFILE
            AREA_CLEAR_ABOVE(L)=AREA_CLEAR(L)
         ENDDO
         IF (I_CLOUD_REPRESENTATION.EQ.IP_CLOUD_CONV_STRAT) THEN
            DO L=1, N_PROFILE
               AREA_STRAT_ABOVE(L)=W_CLOUD(L, I_INV)
     &            *FRAC_CLOUD(L, I_INV, IP_CLOUD_TYPE_STRAT)
            ENDDO
         ELSE IF (I_CLOUD_REPRESENTATION.EQ.IP_CLOUD_CSIW) THEN
            DO L=1, N_PROFILE
               AREA_STRAT_ABOVE(L)=AREA_STRAT(L)
               AREA_CONV_ABOVE(L)=AREA_CONV(L)
            ENDDO
         ENDIF
!
      ENDDO
!
!
!
      IF (L_OBSERVED_RE) THEN
!        SCATTER THE DIAGNOSTICS BACK TO THE OUTPUT ARRAYS AND CONVERT
!        TO MICRONS (TO AVOID FIELDS BEING CORRUPTED BY PACKING).
         DO L=1, N_PROFILE
            WEIGHTED_RE(I_GATHER(L))=1.0E+06*WEIGHTED_RE_G(L)
            SUM_WEIGHT_RE(I_GATHER(L))=SUM_WEIGHT_RE_G(L)
         ENDDO
      ENDIF
!
!
!
      RETURN
      END
