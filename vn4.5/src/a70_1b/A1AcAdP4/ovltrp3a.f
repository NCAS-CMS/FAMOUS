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
!+ Subroutine to find energy transfer coefficients for triple overlap.
!
! Method:
!       Energy transfer coefficients for upward and downward radiation
!       at the edges of the layers are calculated assuming maximal
!       overlap of regions of the same nature and random overlap of
!       regions of a different nature.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.2             24-05-96                Original Code
!       4.3             20-02-97                Vector searching
!                                               routine WHENFGT
!                                               replaced by IF-tests.
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE OVERLAP_TRIPLE(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &   , W_CLOUD, W_FREE, FRAC_REGION
     &   , CLOUD_OVERLAP
     &   , NPD_PROFILE, NPD_LAYER
     &   )
!
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
!     MODULE FOR SETTING MACHINE PRECISION.
!
      REAL
     &     TOL_MACHINE
!             MACHINE TOLERANCE
     &   , SQRT_TOL_MACHINE
!             SQRT OF MACHINE TOLERANCE
!
!
!     THE PRECISION SHOULD BE ABOUT 2/2^(SIZE OF SIGNIFICAND)
!
!     THE IEEE-FORMAT USES 53 BITS FOR THE SIGNIFICAND
!     IN DOUBLE PRECISION
!
!     THE CRAY FORMAT USES 47 BITS IN SINGLE PRECISION.
!
      PARAMETER(
     &     TOL_MACHINE=1.42E-14
     &   , SQRT_TOL_MACHINE=1.19E-7
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE FOR SETTING MACHINE-DEPENDENT TOLERANCES.
!     (THE COMDECK PRMCH3A MUST ALWAYS BE INCLUDED BEFORE THIS COMDECK.)
!
      REAL
     &     TOL_DIV
!             TOLERANCE FOR DIVISION
     &   , TOL_TEST
!             TOLERANCE FOR TESTING EQUALITY
!
      PARAMETER(
     &     TOL_DIV=3.2E+01*TOL_MACHINE
     &   , TOL_TEST=1.6E+01*TOL_MACHINE
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
!     MODULE TO SET POINTERS IN CLOUD_OVERLAP.
!
!     NOTE THAT SEVERAL POINTERS ARE IDENTICAL SINCE ONLY CERTAIN
!     GROUPS OF COEFFICIENTS ARE RELEVANT TO A PARTICULAR SCHEME.
!
      INTEGER
     &     IP_CLOVLP_GFF
!             POINTER TO GAMMA-FREE-FREE
     &   , IP_CLOVLP_GFC
!             POINTER TO GAMMA-FREE-CLOUD
     &   , IP_CLOVLP_GCF
!             POINTER TO GAMMA-CLOUD-FREE
     &   , IP_CLOVLP_GCC
!             POINTER TO GAMMA-CLOUD-CLOUD
     &   , IP_CLOVLP_BFF
!             POINTER TO BETA-FREE-FREE
     &   , IP_CLOVLP_BFC
!             POINTER TO BETA-FREE-CLOUD
     &   , IP_CLOVLP_BCF
!             POINTER TO BETA-CLOUD-FREE
     &   , IP_CLOVLP_BCC
!             POINTER TO BETA-CLOUD-CLOUD
     &   , IP_CLOVLP_GFM
!             POINTER TO GAMMA_F-
     &   , IP_CLOVLP_GFP
!             POINTER TO GAMMA_F+
     &   , IP_CLOVLP_BFM
!             POINTER TO BETA_F-
     &   , IP_CLOVLP_BFP
!             POINTER TO BETA_F+
     &   , IP_CLOVLP_GM
!             POINTER TO GAMMA_-
     &   , IP_CLOVLP_GP
!             POINTER TO GAMMA_+
     &   , IP_CLOVLP_BM
!             POINTER TO BETA_-
     &   , IP_CLOVLP_BP
!             POINTER TO BETA_+
!
!     POINTERS FOR TRIPLE OVERLAPS:
      INTEGER
     &     IP_CLOVLP_V11
     &   , IP_CLOVLP_V12
     &   , IP_CLOVLP_V13
     &   , IP_CLOVLP_V21
     &   , IP_CLOVLP_V22
     &   , IP_CLOVLP_V23
     &   , IP_CLOVLP_V31
     &   , IP_CLOVLP_V32
     &   , IP_CLOVLP_V33
     &   , IP_CLOVLP_U11
     &   , IP_CLOVLP_U12
     &   , IP_CLOVLP_U13
     &   , IP_CLOVLP_U21
     &   , IP_CLOVLP_U22
     &   , IP_CLOVLP_U23
     &   , IP_CLOVLP_U31
     &   , IP_CLOVLP_U32
     &   , IP_CLOVLP_U33
!
      PARAMETER(
     &     IP_CLOVLP_GFF=1
     &   , IP_CLOVLP_GFC=2
     &   , IP_CLOVLP_GCF=3
     &   , IP_CLOVLP_GCC=4
     &   , IP_CLOVLP_BFF=5
     &   , IP_CLOVLP_BFC=6
     &   , IP_CLOVLP_BCF=7
     &   , IP_CLOVLP_BCC=8
     &   , IP_CLOVLP_GFM=5
     &   , IP_CLOVLP_GFP=6
     &   , IP_CLOVLP_BFM=7
     &   , IP_CLOVLP_BFP=8
     &   , IP_CLOVLP_GM=5
     &   , IP_CLOVLP_GP=6
     &   , IP_CLOVLP_BM=7
     &   , IP_CLOVLP_BP=8
     &   )
!
      PARAMETER(
     &     IP_CLOVLP_V11=1
     &   , IP_CLOVLP_V12=2
     &   , IP_CLOVLP_V13=3
     &   , IP_CLOVLP_V21=4
     &   , IP_CLOVLP_V22=5
     &   , IP_CLOVLP_V23=6
     &   , IP_CLOVLP_V31=7
     &   , IP_CLOVLP_V32=8
     &   , IP_CLOVLP_V33=9
     &   , IP_CLOVLP_U11=10
     &   , IP_CLOVLP_U12=11
     &   , IP_CLOVLP_U13=12
     &   , IP_CLOVLP_U21=13
     &   , IP_CLOVLP_U22=14
     &   , IP_CLOVLP_U23=15
     &   , IP_CLOVLP_U31=16
     &   , IP_CLOVLP_U32=17
     &   , IP_CLOVLP_U33=18
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
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , N_CLOUD_TOP
!             TOPMOST CLOUDY LAYER
      REAL      !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUD AMOUNTS
     &   , FRAC_REGION(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             FRACTIONS OF TOTAL CLOUD AMOUNT OCCUPIED BY
!             DIFFERENT REGIONS
!
      REAL      !, INTENT(OUT)
     &     W_FREE(NPD_PROFILE, NPD_LAYER)
!             CLOUD-FREE AMOUNTS
     &   , CLOUD_OVERLAP(NPD_PROFILE, 0: NPD_LAYER, NPD_OVERLAP_COEFF)
!             COEFFICIENTS FOR TRANSFER OF ENERGY AT INTERFACE
!
!
!     LOCAL ARGUMENTS.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
     &   , K
!             LOOP VARIABLE
!
!
!     FIXED LOCAL VALUES:
      INTEGER
     &     N_REGION
!             NUMBER OF REGIONS
      REAL
     &     AREA_LOWER(NPD_PROFILE, NPD_REGION)
!             AREAS OF REGIONS IN LOWER LAYER
     &   , AREA_UPPER(NPD_PROFILE, NPD_REGION)
!             AREAS OF REGIONS IN LOWER LAYER
     &   , AREA_OVERLAP(NPD_PROFILE, NPD_REGION, NPD_REGION)
!             AREAS OF REGIONS IN LOWER LAYER
     &   , AREA_RANDOM(NPD_PROFILE)
!             AREAS OF RANDOM OVERLAP
!
!
!
!     SET THE FREE FRACTIONS IN EACH LAYER.
      DO I=1, N_CLOUD_TOP-1
         DO L=1, N_PROFILE
            W_FREE(L, I)=1.0E+00
         ENDDO
      ENDDO
      DO I=N_CLOUD_TOP, N_LAYER
         DO L=1, N_PROFILE
            W_FREE(L, I)=1.0E+00-W_CLOUD(L, I)
         ENDDO
      ENDDO
!
!
!
!     SET NUMBER OF REGIONS FOR POSSIBLE FUTURE EXTENSION.
      N_REGION=3
!
!     WE CONSIDER EACH BOUNDARY IN TURN, COMPARING THE FRACTIONS
!     OF EACH REGION IN THE LAYERS ABOVE AND BELOW THE BOUNDARY.
!
!     INITIALIZE FOR THE LAYER ABOVE THE CLOUDS.
      DO L=1, N_PROFILE
         AREA_UPPER(L, IP_REGION_CLEAR)=1.0E+00
         AREA_UPPER(L, IP_REGION_STRAT)=0.0E+00
         AREA_UPPER(L, IP_REGION_CONV)=0.0E+00
      ENDDO

      DO I=N_CLOUD_TOP-1, N_LAYER
!
!        SET AREAS OF THE REGIONS IN THE LOWER LAYER.
         IF (I.LT.N_LAYER) THEN
            DO L=1, N_PROFILE
               AREA_LOWER(L, IP_REGION_CLEAR)=W_FREE(L, I+1)
               AREA_LOWER(L, IP_REGION_STRAT)=W_CLOUD(L, I+1)
     &            *FRAC_REGION(L, I+1, IP_REGION_STRAT)
               AREA_LOWER(L, IP_REGION_CONV)=W_CLOUD(L, I+1)
     &            *FRAC_REGION(L, I+1, IP_REGION_CONV)
            ENDDO
         ELSE
            DO L=1, N_PROFILE
               AREA_LOWER(L, IP_REGION_CLEAR)=1.0E+00
               AREA_LOWER(L, IP_REGION_STRAT)=0.0E+00
               AREA_LOWER(L, IP_REGION_CONV)=0.0E+00
            ENDDO
         ENDIF
!
!        SET THE AREAS OF OVERLAP BETWEEN LIKE REGIONS.
         DO K=1, N_REGION
            DO L=1, N_PROFILE
               AREA_OVERLAP(L, K, K)=MIN(AREA_LOWER(L, K)
     &            , AREA_UPPER(L, K))
            ENDDO
         ENDDO
!
!        FIND THE AREAS OF OVERLAP BETWEEN UNLIKE REGIONS. THE OVERLAP
!        BETWEEN UNLIKE REGIONS IS ASSUMED TO BE RANDOM. THESE AREAS
!        ARE SET EQUAL TO 0 FOR THE CASE WHERE THERE IS NO SUCH AREA
!        AND ARE RESET WHEN SUCH AN AREA IS PRESENT.
         DO K=1, N_REGION
            DO J=1, K-1
               DO L=1, N_PROFILE
                  AREA_OVERLAP(L, K, J)=0.0E+00
                  AREA_OVERLAP(L, J, K)=0.0E+00
               ENDDO
            ENDDO
         ENDDO
         DO L=1, N_PROFILE
            AREA_RANDOM(L)=1.0E+00
         ENDDO
         DO K=1, N_REGION
            DO L=1, N_PROFILE
               AREA_RANDOM(L)=AREA_RANDOM(L)-AREA_OVERLAP(L, K, K)
            ENDDO
         ENDDO
         DO K=1, N_REGION
            DO J=1, K-1
               DO L=1, N_PROFILE
                  IF (AREA_RANDOM(L).GT.TOL_DIV) THEN
                     AREA_OVERLAP(L, K, J)
     &                  =(AREA_UPPER(L, K)-AREA_OVERLAP(L, K, K))
     &                  *(AREA_LOWER(L, J)-AREA_OVERLAP(L, J, J))
     &                  /AREA_RANDOM(L)
                     AREA_OVERLAP(L, J, K)
     &                  =(AREA_UPPER(L, J)-AREA_OVERLAP(L, J, J))
     &                  *(AREA_LOWER(L, K)-AREA_OVERLAP(L, K, K))
     &                  /AREA_RANDOM(L)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
!
!        NOW PROCEED TO FIND THE ENERGY TRANSFER COEFFICIENTS
!        BETWEEN THE VARIOUS REGIONS.
!
!        COEFFICIENTS FOR THE DOWNWARD TRANSFER OF ENERGY:
!
!        TO AVOID DIVISION BY 0 WE INITIALIZE TO DEFAULT VALUES
!        AND RESET.
         DO L=1, N_PROFILE
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V11)=1.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V21)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V31)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V12)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V22)=1.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V32)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V13)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V23)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_V33)=1.0E+00
         ENDDO
!
!        TRANSFER FROM CLEAR-SKY REGION:
         DO L=1, N_PROFILE
            IF (AREA_UPPER(L, IP_REGION_CLEAR).GT.TOL_DIV) THEN
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V11)
     &            =AREA_OVERLAP(L, IP_REGION_CLEAR, IP_REGION_CLEAR)
     &            /AREA_UPPER(L, IP_REGION_CLEAR)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V21)
     &            =AREA_OVERLAP(L, IP_REGION_CLEAR, IP_REGION_STRAT)
     &            /AREA_UPPER(L, IP_REGION_CLEAR)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V31)
     &            =AREA_OVERLAP(L, IP_REGION_CLEAR, IP_REGION_CONV)
     &            /AREA_UPPER(L, IP_REGION_CLEAR)
            ENDIF
         ENDDO
!
!        TRANSFER FROM STRATIFORM REGION:
         DO L=1, N_PROFILE
            IF (AREA_UPPER(L, IP_REGION_STRAT).GT.TOL_DIV) THEN
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V12)
     &            =AREA_OVERLAP(L, IP_REGION_STRAT, IP_REGION_CLEAR)
     &            /AREA_UPPER(L, IP_REGION_STRAT)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V22)
     &            =AREA_OVERLAP(L, IP_REGION_STRAT, IP_REGION_STRAT)
     &            /AREA_UPPER(L, IP_REGION_STRAT)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V32)
     &            =AREA_OVERLAP(L, IP_REGION_STRAT, IP_REGION_CONV)
     &            /AREA_UPPER(L, IP_REGION_STRAT)
            ENDIF
         ENDDO
!
!        TRANSFER FROM CONVECTIVE REGION:
         DO L=1, N_PROFILE
            IF (AREA_UPPER(L, IP_REGION_CONV).GT.TOL_DIV) THEN
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V13)
     &            =AREA_OVERLAP(L, IP_REGION_CONV, IP_REGION_CLEAR)
     &            /AREA_UPPER(L, IP_REGION_CONV)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V23)
     &            =AREA_OVERLAP(L, IP_REGION_CONV, IP_REGION_STRAT)
     &            /AREA_UPPER(L, IP_REGION_CONV)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_V33)
     &            =AREA_OVERLAP(L, IP_REGION_CONV, IP_REGION_CONV)
     &            /AREA_UPPER(L, IP_REGION_CONV)
            ENDIF
         ENDDO
!
!
!        TRANSFER COEFFICIENTS FOR UPWARD FLOW OF ENERGY:
!
!        TO AVOID DIVISION BY 0 WE INITIALIZE TO DEFAULT VALUES
!        AND RESET.
         DO L=1, N_PROFILE
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U11)=1.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U21)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U31)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U12)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U22)=1.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U32)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U13)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U23)=0.0E+00
            CLOUD_OVERLAP(L, I, IP_CLOVLP_U33)=1.0E+00
         ENDDO
!
!        TRANSFER FROM CLEAR-SKY REGION:
         DO L=1, N_PROFILE
            IF (AREA_LOWER(L, IP_REGION_CLEAR).GT.TOL_DIV) THEN
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U11)
     &            =AREA_OVERLAP(L, IP_REGION_CLEAR, IP_REGION_CLEAR)
     &            /AREA_LOWER(L, IP_REGION_CLEAR)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U21)
     &            =AREA_OVERLAP(L, IP_REGION_STRAT, IP_REGION_CLEAR)
     &            /AREA_LOWER(L, IP_REGION_CLEAR)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U31)
     &            =AREA_OVERLAP(L, IP_REGION_CONV, IP_REGION_CLEAR)
     &            /AREA_LOWER(L, IP_REGION_CLEAR)
            ENDIF
         ENDDO
!
!        TRANSFER FROM STRATIFORM REGION:
         DO L=1, N_PROFILE
            IF (AREA_LOWER(L, IP_REGION_STRAT).GT.TOL_DIV) THEN
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U12)
     &            =AREA_OVERLAP(L, IP_REGION_CLEAR, IP_REGION_STRAT)
     &            /AREA_LOWER(L, IP_REGION_STRAT)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U22)
     &            =AREA_OVERLAP(L, IP_REGION_STRAT, IP_REGION_STRAT)
     &            /AREA_LOWER(L, IP_REGION_STRAT)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U32)
     &            =AREA_OVERLAP(L, IP_REGION_CONV, IP_REGION_STRAT)
     &            /AREA_LOWER(L, IP_REGION_STRAT)
            ENDIF
         ENDDO
!
!        TRANSFER FROM CONVECTIVE REGION:
         DO L=1, N_PROFILE
            IF (AREA_LOWER(L, IP_REGION_CONV).GT.TOL_DIV) THEN
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U13)
     &            =AREA_OVERLAP(L, IP_REGION_CLEAR, IP_REGION_CONV)
     &            /AREA_LOWER(L, IP_REGION_CONV)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U23)
     &            =AREA_OVERLAP(L, IP_REGION_STRAT, IP_REGION_CONV)
     &            /AREA_LOWER(L, IP_REGION_CONV)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_U33)
     &            =AREA_OVERLAP(L, IP_REGION_CONV, IP_REGION_CONV)
     &            /AREA_LOWER(L, IP_REGION_CONV)
            ENDIF
         ENDDO
!
!        REASSIGN THE FRACTIONS IN THE UPPER LAYER TO STEP DOWN
!        THROUGH THE ATMOSPHERE.
         IF (I.LT.N_LAYER) THEN
            DO K=1, N_REGION
               DO L=1, N_PROFILE
                  AREA_UPPER(L, K)=AREA_LOWER(L, K)
               ENDDO
            ENDDO
         ENDIF
!
      ENDDO
!
!
!
      RETURN
      END
