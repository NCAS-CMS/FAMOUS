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
!+ Subroutine to set pointers to types of clouds
!
! Method:
!       The types of condensate included are examined. Their phases
!       are set and depending on the representation of clouds adopted
!       it is determined to which type of cloud they contribute.
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
      SUBROUTINE SET_CLOUD_POINTER(IERR
     &   , N_CONDENSED, TYPE_CONDENSED, I_CLOUD_REPRESENTATION
     &   , L_DROP, L_ICE
     &   , I_PHASE_CMP, I_CLOUD_TYPE, L_CLOUD_CMP
     &   )
!
!
!
      IMPLICIT NONE
!
!
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
!
!     DUMMY VARIABLES.
      INTEGER   !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER   !, INTENT(IN)
     &     N_CONDENSED
!             NUMBER OF CONDENSED COMPONENTS
     &   , TYPE_CONDENSED(NPD_CLOUD_COMPONENT)
!             TYPES OF COMPONENTS
     &   , I_CLOUD_REPRESENTATION
!             REPRESENTATION OF CLOUDS USED
      LOGICAL   !, INTENT(IN)
     &     L_DROP
!             FLAG FOR INCLUSION OF DROPLETS
     &   , L_ICE
!             FLAG FOR INCLUSION OF ICE CRYSTALS
!
      INTEGER   !, INTENT(OUT)
     &     I_PHASE_CMP(NPD_CLOUD_COMPONENT)
!             PHASES OF COMPONENTS
     &   , I_CLOUD_TYPE(NPD_CLOUD_COMPONENT)
!             TYPES OF CLOUD TO WHICH EACH COMPONENT CONTRIBUTES
      LOGICAL   !, INTENT(OUT)
     &     L_CLOUD_CMP(NPD_CLOUD_COMPONENT)
!             LOGICAL SWITCHES TO INCLUDE COMPONENTS
!
!
!     LOCAL VARIABLES
      INTEGER
     &     K
!            LOOP VARIABLE
!
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET PROPERTIES OF REPRESENTATIONS OF CLOUDS.
!
      DATA
     &    NP_CLOUD_TYPE(IP_CLOUD_HOMOGEN)/1/
     &  , NP_CLOUD_TYPE(IP_CLOUD_ICE_WATER)/2/
     &  , NP_CLOUD_TYPE(IP_CLOUD_CONV_STRAT)/2/
     &  , NP_CLOUD_TYPE(IP_CLOUD_CSIW)/4/
!
!
!     THE ARRAY IP_CLOUD_TYPE_MAP INDICATES TO WHICH TYPE OF CLOUD
!     EACH COMPONENT BELONGS IN A PARTICULAR REPRESENTATION. AN
!     ENTRY OF 0 INDICATES THAT THAT COMPONENT SHOULD NOT BE
!     PRESENT IN THE REPRESENTATION.
!
      DATA
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_HOMOGEN)
     &       /IP_CLOUD_TYPE_HOMOGEN/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_HOMOGEN)
     &       /IP_CLOUD_TYPE_HOMOGEN/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_HOMOGEN)
     &       /0/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_HOMOGEN)
     &       /0/
      DATA
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_ICE_WATER)
     &       /IP_CLOUD_TYPE_WATER/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_ICE_WATER)
     &       /IP_CLOUD_TYPE_ICE/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_ICE_WATER)
     &       /0/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_ICE_WATER)
     &       /0/
      DATA
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_CONV_STRAT)
     &       /IP_CLOUD_TYPE_STRAT/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_CONV_STRAT)
     &       /IP_CLOUD_TYPE_CONV/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_CONV_STRAT)
     &       /IP_CLOUD_TYPE_STRAT/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_CONV_STRAT)
     &       /IP_CLOUD_TYPE_CONV/
      DATA
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_CSIW)
     &       /IP_CLOUD_TYPE_SW/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_CSIW)
     &       /IP_CLOUD_TYPE_SI/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_CSIW)
     &       /IP_CLOUD_TYPE_CW/
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_CSIW)
     &       /IP_CLOUD_TYPE_CI/
!
!     ------------------------------------------------------------------
!
!
!
      DO K=1, N_CONDENSED
!
         I_CLOUD_TYPE(K)=IP_CLOUD_TYPE_MAP(TYPE_CONDENSED(K)
     &      , I_CLOUD_REPRESENTATION)
!
!        CHECK FOR 0 FLAGGING ILLEGAL TYPES.
         IF (I_CLOUD_TYPE(K).EQ.0) THEN
            WRITE(IU_ERR, '(/A)')
     &         '*** ERROR: A COMPONENT IS NOT COMPATIBLE WITH THE'
     &         //'REPRESENTATION OF CLOUDS SELECTED.'
            IERR=I_ERR_FATAL
            RETURN
         ENDIF
!
         IF (TYPE_CONDENSED(K).EQ.IP_CLCMP_ST_WATER) THEN
!
            I_PHASE_CMP(K)=IP_PHASE_WATER
            L_CLOUD_CMP(K)=L_DROP
!
         ELSE IF (TYPE_CONDENSED(K).EQ.IP_CLCMP_ST_ICE) THEN
!
            I_PHASE_CMP(K)=IP_PHASE_ICE
            L_CLOUD_CMP(K)=L_ICE
!
         ELSE IF (TYPE_CONDENSED(K).EQ.IP_CLCMP_CNV_WATER) THEN
!
            I_PHASE_CMP(K)=IP_PHASE_WATER
            L_CLOUD_CMP(K)=L_DROP
!
         ELSE IF (TYPE_CONDENSED(K).EQ.IP_CLCMP_CNV_ICE) THEN
!
            I_PHASE_CMP(K)=IP_PHASE_ICE
            L_CLOUD_CMP(K)=L_ICE
!
         ENDIF
!
      ENDDO
!
!
!
      RETURN
      END
