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
!+ Subroutine to calculate two-stream coefficients in the regions.
!
! Method:
!       The coefficients for each region are determined and
!       averaged.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.2             15-05-96                Original Code
!                                               (J. M. Edwards)
!       4.3             20-02-97                Calls to vector
!                                               searching routine
!                                               WHENFGT removed.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Unnecessary final
!                                               dimensions removed
!                                               from arrays.
!                                               (J. M. Edwards)
!LL  4.5  27/04/98  Add Fujitsu vectorization directive. 
!LL                                           RBarnes@ecmwf.int
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
! Fujitsu directive to encourage vectorization for whole routine
!OCL NOVREC
      SUBROUTINE TWO_COEFF_REGION(IERR
     &   , N_PROFILE, N_LAYER, N_CLOUD_TOP
     &   , I_2STREAM, L_IR_SOURCE_QUAD, N_SOURCE_COEFF
     &   , N_CLOUD_TYPE, FRAC_CLOUD
     &   , I_REGION_CLOUD, FRAC_REGION
     &   , ASYMMETRY_FREE, OMEGA_FREE, TAU_FREE
     &   , ASYMMETRY_CLOUD, OMEGA_CLOUD, TAU_CLOUD
     &   , ISOLIR, SEC_0
     &   , TRANS, REFLECT, TRANS_0
     &   , SOURCE_COEFF
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
!
!
!     DUMMY ARGUMENTS.
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
     &   , ISOLIR
!             SPECTRAL REGION
     &   , N_CLOUD_TYPE
!             NUMBER OF TYPES OF CLOUDS
     &   , I_2STREAM
!             TWO STREAM SCHEME
     &   , N_SOURCE_COEFF
!             NUMBER OF SOURCE COEFFICIENTS
!
      INTEGER   !, INTENT(IN)
     &     I_REGION_CLOUD(NPD_CLOUD_TYPE)
!             REGIONS IN WHICH TYPES OF CLOUDS FALL
!
      LOGICAL   !, INTENT(IN)
     &     L_IR_SOURCE_QUAD
!             USE A QUADRATIC SOURCE IN THE INFRA-RED
!
!     OPTICAL PROPERTIES OF LAYER:
      REAL      !, INTENT(IN)
     &     FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             FRACTIONS OF DIFFERENT TYPES OF CLOUDS
     &   , FRAC_REGION(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             FRACTIONS OF TOTAL CLOUD OCCUPIED BY EACH REGION
     &   , ASYMMETRY_FREE(NPD_PROFILE, NPD_LAYER)
!             CLEAR-SKY ASYMMETRY FACTOR
     &   , OMEGA_FREE(NPD_PROFILE, NPD_LAYER)
!             CLEAR-SKY ALBEDO OF SINGLE SCATTERING
     &   , TAU_FREE(NPD_PROFILE, NPD_LAYER)
!             CLEAR-SKY OPTICAL DEPTH
     &   , ASYMMETRY_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             ASYMMETRY FACTOR
     &   , OMEGA_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             ALBEDO OF SINGLE SCATTERING
     &   , TAU_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             OPTICAL DEPTH
!
!     SOLAR BEAM
      REAL      !, INTENT(IN)
     &     SEC_0(NPD_PROFILE)
!             SECANT OF ZENITH ANGLE
!
!
!     COEFFICIENTS IN THE TWO-STREAM EQUATIONS:
      REAL      !, INTENT(OUT)
     &     TRANS(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             DIFFUSE TRANSMISSION COEFFICIENT
     &   , REFLECT(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             DIFFUSE REFLECTION COEFFICIENT
     &   , TRANS_0(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             DIRECT TRANSMISSION COEFFICIENT
     &   , SOURCE_COEFF(NPD_PROFILE, NPD_LAYER
     &      , NPD_SOURCE_COEFF, NPD_REGION)
!             SOURCE COEFFICIENTS IN TWO-STREAM EQUATIONS
!
!     LOCAL VARIABLES.
      INTEGER
     &     N_REGION
!             NUMBER OF REGIONS
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
     &   , K
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
     &   , I_REGION
!             LOOP VARIABLE OVER REGIONS
!
!     COEFFICIENTS IN THE TWO-STREAM EQUATIONS:
      REAL      !, INTENT(OUT)
     &     TRANS_TEMP(NPD_PROFILE, NPD_LAYER)
!             TEMPORARY DIFFUSE TRANSMISSION COEFFICIENT
     &   , REFLECT_TEMP(NPD_PROFILE, NPD_LAYER)
!             TEMPORARY DIFFUSE REFLECTION COEFFICIENT
     &   , TRANS_0_TEMP(NPD_PROFILE, NPD_LAYER)
!             TEMPORARY DIRECT TRANSMISSION COEFFICIENT
     &   , SOURCE_COEFF_TEMP(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)
!             TEMPORARY SOURCE COEFFICIENTS IN TWO-STREAM EQUATIONS
!
!     VARIABLES FOR GATHERING:
      INTEGER
     &     N_LIST
!             NUMBER OF POINTS IN LIST
     &   , L_LIST(NPD_PROFILE)
!             LIST OF COLLECTED POINTS
     &   , LL
      REAL
     &     TAU_GATHERED(NPD_PROFILE, NPD_LAYER)
!             GATHERED OPTICAL DEPTH
     &   , OMEGA_GATHERED(NPD_PROFILE, NPD_LAYER)
!             GATHERED ALEBDO OF SINGLE SCATTERING
     &   , ASYMMETRY_GATHERED(NPD_PROFILE, NPD_LAYER)
!             GATHERED ASYMMETRY
     &   , SEC_0_GATHERED(NPD_PROFILE)
!             GATHERED ASYMMETRY
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     TWO_COEFF
!
!     CRAY DIRECTIVES FOR THE WHOLE ROUTINE:
!     POINTS ARE NOT REPEATED IN THE INDEXING ARRAY, SO IT IS SAFE
!     TO VECTORIZE OVER INDIRECTLY ADDRESSED ARRAYS.
Cfpp$ NODEPCHK R
!
!
!
!     FOR THE TRIPLE OVERLAP THE NUMBER OF REGIONS IS 3.
      N_REGION=3
!
!     DETERMINE THE OPTICAL PROPERTIES OF THE CLEAR-SKY REGIONS OF
!     THE LAYERS.
!
      CALL TWO_COEFF(IERR
     &   , N_PROFILE, 1, N_LAYER
     &   , I_2STREAM, L_IR_SOURCE_QUAD
     &   , ASYMMETRY_FREE, OMEGA_FREE, TAU_FREE
     &   , ISOLIR, SEC_0
     &   , TRANS(1, 1, IP_REGION_CLEAR)
     &   , REFLECT(1, 1, IP_REGION_CLEAR)
     &   , TRANS_0(1, 1, IP_REGION_CLEAR)
     &   , SOURCE_COEFF(1, 1, 1, IP_REGION_CLEAR)
     &   , NPD_PROFILE, NPD_LAYER
     &   )
      IF (IERR.NE.I_NORMAL) RETURN
!
!
!     NOW DEAL WITH CLOUDS.
!
!     INITIALIZE THE FULL ARRAYS FOR CLOUDY REGIONS.
!
      DO I_REGION=1, N_REGION
         IF (I_REGION.NE.IP_REGION_CLEAR) THEN
            DO I=N_CLOUD_TOP, N_LAYER
               DO L=1, N_PROFILE
                  TRANS(L, I, I_REGION)=0.0E+00
                  REFLECT(L, I, I_REGION)=0.0E+00
               ENDDO
            ENDDO
            DO J=1, N_SOURCE_COEFF
               DO I=N_CLOUD_TOP, N_LAYER
                  DO L=1, N_PROFILE
                     SOURCE_COEFF(L, I, J, I_REGION)=0.0E+00
                  ENDDO
               ENDDO
            ENDDO
!
            IF (ISOLIR.EQ.IP_SOLAR) THEN
               DO I=N_CLOUD_TOP, N_LAYER
                  DO L=1, N_PROFILE
                     TRANS_0(L, I, I_REGION)=0.0E+00
                  ENDDO
               ENDDO
            ENDIF
!
         ENDIF
!
      ENDDO
!
!
!
!     CONSIDER EACH TYPE OF CLOUD IN TURN, CHECKING WHICH REGION IT
!     CONTRUBUTES TO AND FORM WEIGHTED SUMS OF CLOUD PROPERTIES.
!
      DO K=1, N_CLOUD_TYPE
!
!
!        SET THE REGION IN WHICH CLOUDS OF THIS TYPE ARE INCLUDED.
         I_REGION=I_REGION_CLOUD(K)
!
         DO I=N_CLOUD_TOP, N_LAYER
!
!           FORM A LIST OF POINTS WHERE CLOUD OF THIS TYPE EXISTS
!           ON THIS ROW FOR GATHERING.
            N_LIST=0
            DO L=1, N_PROFILE
               IF (FRAC_CLOUD(L, I, K).GT.0.0E+00) THEN
                  N_LIST=N_LIST+1
                  L_LIST(N_LIST)=L
               ENDIF
            ENDDO
!
!
            IF (N_LIST.GT.0) THEN
!
!              GATHER THE OPTICAL PROPERTIES. THOUGH WE CONSIDER ONLY
!              ONE LAYER AT A TIME THE LOWER ROUTINES WILL OPERATE ON
!              ARRAYS WITH VERTICAL STRUCTURE, SO THE GATHERED ARRAYS
!              ARE TWO-DIMENSIONAL.
!
               DO L=1, N_LIST
                  TAU_GATHERED(L, I)
     &              =TAU_CLOUD(L_LIST(L), I, K)
                  OMEGA_GATHERED(L, I)
     &              =OMEGA_CLOUD(L_LIST(L), I, K)
                  ASYMMETRY_GATHERED(L, I)
     &              =ASYMMETRY_CLOUD(L_LIST(L), I, K)
               ENDDO
               IF (ISOLIR.EQ.IP_SOLAR) THEN
                  DO L=1, N_LIST
                     SEC_0_GATHERED(L)=SEC_0(L_LIST(L))
                  ENDDO
               ENDIF
!
!
               CALL TWO_COEFF(IERR
     &            , N_LIST, I, I
     &            , I_2STREAM, L_IR_SOURCE_QUAD
     &            , ASYMMETRY_GATHERED, OMEGA_GATHERED
     &            , TAU_GATHERED
     &            , ISOLIR, SEC_0_GATHERED
     &            , TRANS_TEMP, REFLECT_TEMP, TRANS_0_TEMP
     &            , SOURCE_COEFF_TEMP
     &            , NPD_PROFILE, NPD_LAYER
     &            )
               IF (IERR.NE.I_NORMAL) RETURN
!
!
               DO L=1, N_LIST
                  LL=L_LIST(L)
                  TRANS(LL, I, I_REGION)=TRANS(LL, I, I_REGION)
     &               +FRAC_CLOUD(LL, I, K)*TRANS_TEMP(L, I)
                  REFLECT(LL, I, I_REGION)=REFLECT(LL, I, I_REGION)
     &               +FRAC_CLOUD(LL, I, K)*REFLECT_TEMP(L, I)
               ENDDO
               DO J=1, N_SOURCE_COEFF
                  DO L=1, N_LIST
                     LL=L_LIST(L)
                     SOURCE_COEFF(LL, I, J, I_REGION)
     &                  =SOURCE_COEFF(LL, I, J, I_REGION)
     &                  +FRAC_CLOUD(LL, I, K)
     &                  *SOURCE_COEFF_TEMP(L, I, J)
                  ENDDO
               ENDDO
               IF (ISOLIR.EQ.IP_SOLAR) THEN
                  DO L=1, N_LIST
                     LL=L_LIST(L)
                     TRANS_0(LL, I, I_REGION)=TRANS_0(LL, I, I_REGION)
     &                  +FRAC_CLOUD(LL, I, K)*TRANS_0_TEMP(L, I)
                  ENDDO
               ENDIF
!
            ENDIF
!
         ENDDO
      ENDDO
!
!
!     FINALLY, SCALE THE WEIGHTED SUMS BY THE CLOUD FRACTIONS.
      DO I_REGION=1, N_REGION
         IF (I_REGION.NE.IP_REGION_CLEAR) THEN
            DO I=N_CLOUD_TOP, N_LAYER
!
!              GATHER POINTS WITHIN THIS REGION.
               N_LIST=0
               DO L=1,N_PROFILE
                  IF (FRAC_REGION(L, I, I_REGION).GT.0.0E+00) THEN
                     N_LIST=N_LIST+1
                     L_LIST(N_LIST)=L
                  ENDIF
               ENDDO
               DO L=1, N_LIST
                  LL=L_LIST(L)
                  TRANS(LL, I, I_REGION)=TRANS(LL, I, I_REGION)
     &               /FRAC_REGION(LL, I, I_REGION)
                  REFLECT(LL, I, I_REGION)=REFLECT(LL, I, I_REGION)
     &               /FRAC_REGION(LL, I, I_REGION)
               ENDDO
               DO J=1, N_SOURCE_COEFF
                  DO L=1, N_LIST
                     LL=L_LIST(L)
                     SOURCE_COEFF(LL, I, J, I_REGION)
     &                  =SOURCE_COEFF(LL, I, J, I_REGION)
     &                  /FRAC_REGION(LL, I, I_REGION)
                  ENDDO
               ENDDO
               IF (ISOLIR.EQ.IP_SOLAR) THEN
                  DO L=1, N_LIST
                     LL=L_LIST(L)
                     TRANS_0(LL, I, I_REGION)=TRANS_0(LL, I, I_REGION)
     &                  /FRAC_REGION(LL, I, I_REGION)
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO
!
!
!
      RETURN
      END
