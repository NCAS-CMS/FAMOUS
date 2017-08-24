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
!+ Subroutine to find maximally overlapped energy transfer coefficients.
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
!       4.2             Nov. 96   T3E migration: CALL WHENFGT replaced
!                                  by portable fortran code.
!                                                S.J.Swarbrick
!       4.5             18-05-98                Reference to obsolete
!                                               solvers removed.
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
      SUBROUTINE OVERLAP_MIX_MAXIMUM(N_PROFILE, N_LAYER, N_CLOUD_TOP
     &   , ISOLIR, I_SOLVER
     &   , W_CLOUD, W_FREE
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
!     MODULE TO DEFINE REFERENCE NUMBERS FOR SOLVERS.
!
      INTEGER
     &     IP_SOLVER_PENTADIAGONAL
!             PENTADIAGONAL SCHEME
     &   , IP_SOLVER_MIX_11
!             MIXED COLUMN SCHEME USING FULL ENDEKADIAGONAL MATRIX
     &   , IP_SOLVER_MIX_APP_SCAT
!             MIXED COLUMN SCHEME WITH APPROXIMATE SCATTERING
     &   , IP_SOLVER_MIX_DIRECT
!             DIRECT MIXED COLUMN SCHEME FOR FULL FLUXES 
     &   , IP_SOLVER_HOMOGEN_DIRECT
!             DIRECT SOLVER FOR A HOMOGENEOUS COLUMN
     &   , IP_SOLVER_TRIPLE
!             DIRECT SOLVER FOR TRIPLE COLUMN
     &   , IP_SOLVER_TRIPLE_APP_SCAT
!             DIRECT SOLVER FOR TRIPLE COLUMN APPROXIMATING SCATTERING
!
      PARAMETER(
     &     IP_SOLVER_PENTADIAGONAL=1
     &   , IP_SOLVER_MIX_11=6
     &   , IP_SOLVER_MIX_APP_SCAT=9
     &   , IP_SOLVER_MIX_DIRECT=11
     &   , IP_SOLVER_HOMOGEN_DIRECT=13
     &   , IP_SOLVER_TRIPLE=14
     &   , IP_SOLVER_TRIPLE_APP_SCAT=15
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
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , N_CLOUD_TOP
!             TOPMOST CLOUDY LAYER
     &   , ISOLIR
!             SPECTRAL REGION
     &   , I_SOLVER
!             SOLVER TO BE USED
      REAL      !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUD AMOUNTS
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
     &   , L
!             LOOP VARIABLE
     &   , K
!             LOOP VARIABLE
     &   , KK
!             LOOP VARIABLE
     &   , N_INDEX
!             NUMBER OF INDICES IN WHEN TEST
     &   , INDEX(NPD_PROFILE)
!             INDICES OF WHEN TEST
      REAL
     &     G_FF(NPD_PROFILE, 0: NPD_LAYER)
!             BASIC TRANSFER COEFFICIENT
     &   , G_FC(NPD_PROFILE, 0: NPD_LAYER)
!             BASIC TRANSFER COEFFICIENT
     &   , G_CF(NPD_PROFILE, 0: NPD_LAYER)
!             BASIC TRANSFER COEFFICIENT
     &   , G_CC(NPD_PROFILE, 0: NPD_LAYER)
!             BASIC TRANSFER COEFFICIENT
     &   , B_FF(NPD_PROFILE, 0: NPD_LAYER)
!             BASIC TRANSFER COEFFICIENT
     &   , B_FC(NPD_PROFILE, 0: NPD_LAYER)
!             BASIC TRANSFER COEFFICIENT
     &   , B_CF(NPD_PROFILE, 0: NPD_LAYER)
!             BASIC TRANSFER COEFFICIENT
     &   , B_CC(NPD_PROFILE, 0: NPD_LAYER)
!             BASIC TRANSFER COEFFICIENT
!
!     SUBROUTINES CALLED:
!
!     CRAY DIRECTIVES FOR THE WHOLE ROUTINE:
!     POINTS ARE NOT REPEATED IN THE INDEXING ARRAY, SO IT IS SAFE
!     TO VECTORIZE OVER INDIRECTLY ADDRESSED ARRAYS.
Cfpp$ NODEPCHK R
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
!     EVALUATE THE EXTENT OF OVERLAP BETWEEN LAYERS OF CLOUD
!     AT THE INTERFACE BETWEEN THE ITH AND (I+1)ST LAYER ON THE
!     ASSUMPTION OF MAXIMUM OVERLAP BETWEEN ADJACENT LAYERS.
!     THE TOP AND BOTTOM BOUNDARIES ARE EXCEPTIONAL.
!
!     HERE IS IS EASIEST TO FIND THE BASIC TRANSFER COEFFICIENTS
!     AND DEDUCE ALL OTHERS FROM THEM.
!
!     COEFFICIENTS ARE NOT REQUIRED IF THERE ARE NO CLOUDY LAYERS,
!     IN WHICH CASE N_CLOUD_TOP WILL BE SET TO N_LAYER+1.
      IF (N_CLOUD_TOP.LE.N_LAYER) THEN
         DO L=1, N_PROFILE
!
!           TOPMOST CLOUDY BOUNDARY:
            G_FF(L, N_CLOUD_TOP-1)=W_FREE(L, N_CLOUD_TOP)
            G_FC(L, N_CLOUD_TOP-1)=W_CLOUD(L, N_CLOUD_TOP)
            G_CF(L, N_CLOUD_TOP-1)=0.0E+00
            G_CC(L, N_CLOUD_TOP-1)=1.0E+00
            B_FF(L, N_CLOUD_TOP-1)=1.0E+00
            B_FC(L, N_CLOUD_TOP-1)=1.0E+00
            B_CF(L, N_CLOUD_TOP-1)=0.0E+00
            B_CC(L, N_CLOUD_TOP-1)=0.0E+00
!
!           BOTTOM BOUNDARY:
            G_FF(L, N_LAYER)=1.0E+00
            G_FC(L, N_LAYER)=0.0E+00
            G_CF(L, N_LAYER)=1.0E+00
            G_CC(L, N_LAYER)=0.0E+00
            B_FF(L, N_LAYER)=W_FREE(L, N_LAYER)
            B_FC(L, N_LAYER)=0.0E+00
            B_CF(L, N_LAYER)=W_CLOUD(L, N_LAYER)
            B_CC(L, N_LAYER)=1.0E+00
!
         ENDDO
      ENDIF
!
      DO I=N_CLOUD_TOP, N_LAYER-1
!
!        SET DEFAULT VALUES FOR COMPLETELY CLEAR OR CLOUDY LAYERS
!        AND OVERWRITE IN NORMAL CASES.
!
!        GAMMAS:
!
         DO L=1, N_PROFILE
            G_CC(L, I)=1.0E+00
            G_FF(L, I)=1.0E+00
         ENDDO
!
         N_INDEX=0
         DO L   =1,N_PROFILE
           IF (W_CLOUD(L,I).GT.TOL_DIV) THEN
             N_INDEX =N_INDEX+1
             INDEX(N_INDEX)=L
           END IF
         END DO
!
         DO K=1, N_INDEX
            KK=INDEX(K)
            G_CC(KK, I)=MIN(1.0E+00, W_CLOUD(KK, I+1)/W_CLOUD(KK, I))
         ENDDO
!
         N_INDEX=0
         DO L   =1,N_PROFILE
           IF (W_FREE(L,I).GT.TOL_DIV) THEN
             N_INDEX =N_INDEX+1
             INDEX(N_INDEX)=L
           END IF
         END DO
!
         DO K=1, N_INDEX
            KK=INDEX(K)
            G_FF(KK, I)=MIN(1.0E+00, W_FREE(KK, I+1)/W_FREE(KK, I))
         ENDDO
!        INFER THE OTHER VALUES FROM GENERIC RELATIONS.
         DO L=1, N_PROFILE
            G_FC(L, I)=1.0E+00-G_FF(L, I)
            G_CF(L, I)=1.0E+00-G_CC(L, I)
         ENDDO
!
!
!        BETAS:
!
         DO L=1, N_PROFILE
            B_CC(L, I)=1.0E+00
            B_FF(L, I)=1.0E+00
         ENDDO
!
         N_INDEX=0
         DO L   =1,N_PROFILE
           IF (W_CLOUD(L,I+1).GT.TOL_DIV) THEN
             N_INDEX =N_INDEX+1
             INDEX(N_INDEX)=L
           END IF
         END DO
!
         DO K=1, N_INDEX
            KK=INDEX(K)
            B_CC(KK, I)=MIN(1.0E+00, W_CLOUD(KK, I)/W_CLOUD(KK, I+1))
         ENDDO
!
         N_INDEX=0
         DO L   =1,N_PROFILE
           IF (W_FREE(L,I+1).GT.TOL_DIV) THEN
             N_INDEX =N_INDEX+1
             INDEX(N_INDEX)=L
           END IF
         END DO
!
         DO K=1, N_INDEX
            KK=INDEX(K)
            B_FF(KK, I)=MIN(1.0E+00, W_FREE(KK, I)/W_FREE(KK, I+1))
         ENDDO
!        INFER THE OTHER VALUSE FROM GENERIC RELATIONS.
         DO L=1, N_PROFILE
            B_FC(L, I)=1.0E+00-B_CC(L, I)
            B_CF(L, I)=1.0E+00-B_FF(L, I)
         ENDDO
!
      ENDDO
!
!     NOW CALCULATE THE OVERLAPPED ARRAYS USING THE BASIC COEFFICIENTS.
!
!     IN THE SOLAR REGION COEFFICIENTS FOR DOWNWARD COUPLING OF THE
!     FLUXES ARE REQUIRED. THESE COEFFICIENTS ARE ALSO NEEDED FOR
!     INFRA-RED CALCULATIONS WITH APPROXIMATE SCATTERING, OR
!     WHENEVER THE DIRECT SOLVER IS USED.
!
      IF ( (I_SOLVER.EQ.IP_SOLVER_MIX_DIRECT).OR.
     &     (ISOLIR.EQ.IP_SOLAR).OR.
     &     ( (ISOLIR.EQ.IP_INFRA_RED).AND.
     &     (I_SOLVER.EQ.IP_SOLVER_MIX_APP_SCAT) ) ) THEN
!
         DO I=N_CLOUD_TOP-1, N_LAYER
            DO L=1, N_PROFILE
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GFF)=G_FF(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GFC)=G_FC(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GCF)=G_CF(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GCC)=G_CC(L, I)
            ENDDO
         ENDDO
!
      ENDIF
!
!     WITH APPROXIMATE SCATTERING IN THE LONGWAVE, OR WHENEVER
!     A DIRECT SOLVER IS USED, THE CORRESPONDING
!     UPWARD COEFFICIENTS ARE NEEDED.
!
      IF ( (I_SOLVER.EQ.IP_SOLVER_MIX_DIRECT).OR.
     &     ( (ISOLIR.EQ.IP_INFRA_RED).AND.
     &     (I_SOLVER.EQ.IP_SOLVER_MIX_APP_SCAT) ) ) THEN
!
         DO I=N_CLOUD_TOP-1, N_LAYER
            DO L=1, N_PROFILE
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BFF)=B_FF(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BFC)=B_FC(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BCF)=B_CF(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BCC)=B_CC(L, I)
            ENDDO
         ENDDO
      ENDIF
!
!
      IF (I_SOLVER.EQ.IP_SOLVER_MIX_11) THEN
!
         DO I=N_CLOUD_TOP-1, N_LAYER
            DO L=1, N_PROFILE
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GM)
     &            =G_FF(L, I)+G_CC(L, I)-G_CF(L, I)-G_FC(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GP)
     &            =G_FF(L, I)-G_CC(L, I)+G_CF(L, I)-G_FC(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BM)
     &            =B_FF(L, I)+B_CC(L, I)-B_CF(L, I)-B_FC(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BP)
     &            =B_FF(L, I)-B_CC(L, I)-B_CF(L, I)+B_FC(L, I)
            ENDDO
         ENDDO
!
!
      ENDIF
!
!
!
      RETURN
      END
