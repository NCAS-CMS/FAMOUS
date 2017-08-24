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
!+ Subroutine to split the atmosphere into maximally overlapped columns.
!
! Method:
!       The layers are first ranked in order of increasing cloudiness.
!       This operation cannot be vectorized and is done for one profile
!       at a time. The areal extent of each column and the logical
!       cloud mask are then set.
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
      SUBROUTINE SPLIT_MAXIMUM(N_PROFILE, N_LAYER
     &   , W_CLOUD
     &   , N_COLUMN, AREA_COLUMN, L_COLUMN
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN
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
     &   , NPD_COLUMN
!             NUMBER OF COLUMNS PER POINT
!
!     INCLUDE COMDECKS
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
!
!     DUMMY ARGUMENTS
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
      INTEGER   !, INTENT(INOUT)
     &     N_COLUMN(NPD_PROFILE)
!             NUMBER OF COLUMNS
      LOGICAL   !, INTENT(IN)
     &     L_COLUMN(NPD_PROFILE, NPD_LAYER, NPD_COLUMN)
!             ARRAY OF TYPES
      REAL      !, INTENT(IN)
     &     AREA_COLUMN(NPD_PROFILE, NPD_COLUMN)
!             AREA OF EACH COLUMN
     &   , W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUD AMOUNTS
!
!     LOCAL ARGUMENTS
      INTEGER
     &     IRANK(NPD_LAYER)
!             ARRAY TO RANK COLUMNS BY W
     &   , I
!             LOOP VARIABLE
     &   , K
!             LOOP VARIBLE
     &   , L
!             LOOP VARIBLE
      REAL
     &     W_CLOUD_SINGLE(NPD_LAYER)
!             CLOUD AMOUNTS FOR SINGLE PROFILE
     &   , W
!             SINGLE CLOUD AMOUNT
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     RANK
!
!
!
      DO L=1, N_PROFILE
!        GATHER THE CLOUD AMOUNTS FOR ONE PROFILE
         DO I=1, N_LAYER
            W_CLOUD_SINGLE(I)=W_CLOUD(L, I)
         ENDDO
!
!        FIRST FORM THE VECTOR IRANK, RANKING THE LAYERS IN ORDER OF
!        INCREASING CLOUD CONTENT.
         CALL RANK(N_LAYER
     &      , W_CLOUD_SINGLE, IRANK
     &      )
!
!        PASS THROUGH ALL THE COLUMNS SETTING L_COLUMN EQUAL TO .FALSE.
!        IF THE COLUMN IS ACTUALLY CLEAR ON THAT LEVEL. THE ASSUMPTION
!        OF MAXIMUM OVERLAP IS USED HERE.
         N_COLUMN(L)=1
         W=0.0E+00
         I=1
30       IF (I.LE.N_LAYER) THEN
            IF ( W_CLOUD_SINGLE(IRANK(I)).LT.(W+TOL_TEST) ) THEN
               I=I+1
               GOTO 30
            ELSE
               DO K=1, I-1
                  L_COLUMN(L, IRANK(K), N_COLUMN(L))=.FALSE.
               ENDDO
               DO K=I, N_LAYER
                  L_COLUMN(L, IRANK(K), N_COLUMN(L))=.TRUE.
               ENDDO
               AREA_COLUMN(L, N_COLUMN(L))
     &            =W_CLOUD_SINGLE(IRANK(I))-W
               W=W_CLOUD_SINGLE(IRANK(I))
            ENDIF
            IF (W.LT.W_CLOUD_SINGLE(IRANK(N_LAYER))-TOL_TEST) THEN
               N_COLUMN(L)=N_COLUMN(L)+1
               GOTO 30
            ENDIF
         ENDIF
!
!        THERE IS A TOTALLY CLEAR COLUMN UNLESS AT LEAST ONE LAYER IS
!        TOTALLY CLOUDY.
         IF ((1.0E+00-W).GT.TOL_TEST) THEN
!           INCREMENT THE NUMBER OF COLUMNS IF THE FIRST IS NOT BLANK.
            IF (W.GE.TOL_TEST) THEN
               N_COLUMN(L)=N_COLUMN(L)+1
            ENDIF
            DO K=1, N_LAYER
               L_COLUMN(L, IRANK(K), N_COLUMN(L))=.FALSE.
            ENDDO
            AREA_COLUMN(L, N_COLUMN(L))=1.0E+00-W
         ENDIF
!
      ENDDO
!
!
      RETURN
      END
