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
!+ Subroutine to rank amounts of cloudiness in acolumn.
!
! Method:
!       This routine uses a heap sorting algorithm taken from
!       "Numerical Recipes" by D. E. Knuth to form an array of
!       pointers to layers containing increasing amounts of cloud.
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
      SUBROUTINE RANK(N_LAYER
     &   , W_CLOUD, IRANK
     &   , NPD_LAYER
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS.
      INTEGER   !, INTENT(IN)
     &     NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     DUMMY ARGUMENTS
      INTEGER   !, INTENT(IN)
     &     N_LAYER
!             NUMBER OF LAYERS
      REAL      !, INTENT(IN)
     &     W_CLOUD(NPD_LAYER)
!             CLOUD AMOUNT
      INTEGER   !, INTENT(OUT)
     &     IRANK(NPD_LAYER)
!             RANK ARRAY FOR CLOUDS
!
!
!     LOCAL ARGUMENTS
      INTEGER
     &     I
!             LOOP VARIBALE
     &   , J
!             LOOP VARIBALE
     &   , K
!             LOOP VARIBALE
     &   , IR
!             LOOP VARIBALE
     &   , I_RANK_TEMPORARY
!             TEMPORARY RANK VALUE
      REAL
     &     W
!             SINGLE CLOUD AMOUNT
!
!
!     FORM AN ARRAY RANKING THE AMOUNTS OF CLOUDINESS IN EACH LAYER.
      DO I=1, N_LAYER
         IRANK(I)=I
      ENDDO
!
      K=N_LAYER/2+1
      IR=N_LAYER
20    CONTINUE
         IF (K.GT.1) THEN
            K=K-1
            I_RANK_TEMPORARY=IRANK(K)
            W=W_CLOUD(I_RANK_TEMPORARY)
         ELSE
            I_RANK_TEMPORARY=IRANK(IR)
            W=W_CLOUD(I_RANK_TEMPORARY)
            IRANK(IR)=IRANK(1)
            IR=IR-1
            IF (IR.EQ.1) THEN
               IRANK(1)=I_RANK_TEMPORARY
               RETURN
            ENDIF
         ENDIF
         I=K
         J=K+K
30       IF (J.LE.IR) THEN
            IF (J.LT.IR) THEN
               IF (W_CLOUD(IRANK(J)).LT.W_CLOUD(IRANK(J+1))) THEN
                  J=J+1
               ENDIF
            ENDIF
            IF (W.LT.W_CLOUD(IRANK(J))) THEN
               IRANK(I)=IRANK(J)
               I=J
               J=J+J
            ELSE
               J=IR+1
            ENDIF
         GOTO 30
         ENDIF
         IRANK(I)=I_RANK_TEMPORARY
      GOTO 20
!
!
!
      END
