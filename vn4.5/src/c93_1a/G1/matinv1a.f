C ******************************COPYRIGHT******************************
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
CLL  SUBROUTINE MATINV--------------------------------------------
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL vn4.2  07/11/96:Allow this routine to be used in C93_2A (Farnon)
CLL
CLL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL  VERSION 2, DATED 18/01/90
CLL
CLL  SYSTEM TASK: ????
CLL
CLL  PURPOSE: INVERTS AN N BY N MATRIX BY ROW REDUCTION
CLL
CLL  DOCUMENTATION:  ???
CLL
CLLEND-------------------------------------------------------------
C
C*L ARGUMENTS:-----------------------------------------------------
      SUBROUTINE MATINV
     & (C,B,N)
C
      IMPLICIT NONE
C
      INTEGER
     *  N         ! IN  SIZE OF MATRIX TO BE INVERTED
C
      REAL
     *  C(N,N)    ! IN  MATRIX TO BE INVERTED
     *, B(N,N)    ! OUT INVERTED MATRIX
C*-----------------------------------------------------------------
C
C*L WORKSPACE USAGE------------------------------------------------
      REAL
     *  A(N,N)
C*-----------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED------------------------------------
C     NONE
C*-----------------------------------------------------------------
C
C------------------------------------------------------------------
C DEFINE LOCAL VARIABLES
C------------------------------------------------------------------
      INTEGER
     *  M1,M2,J,M  ! LOOP COUNTERS
     *, L
C
      REAL
     *  AMM   !  = A(M,M)
     *, AMJ   !  = A(M,J)
C
C
C------------------------------------------------------------------
CL  1. SET UP EXTENDED MATRIX (A|B)=(C|I) TO ROW REDUCE
C------------------------------------------------------------------
      DO 400 M1=1,N
        DO 401 M2=1,N
          IF (M1.EQ.M2) THEN
            B(M1,M2)=1
          ELSE
            B(M1,M2)=0
          ENDIF
          A(M1,M2)=C(M1,M2)
  401   CONTINUE
  400 CONTINUE
C
C------------------------------------------------------------------
CL 2. FOR EACH COLUMN M FIND THE FIRST NON-ZERO ENTRY.
CL    ADD ROW WITH NON-ZERO ENTRY TO ROW M SO THAT A(M,M) IS NOT 0.
CL    IF COLUMN M CONSISTS OF ZEROS THEN MATRIX CANNOT BE INVERTED,
CL    SO STOP.
C------------------------------------------------------------------
      DO 2 M=1,N
        L=M
  200   CONTINUE
        IF (A(M,M).EQ.0.) THEN
          IF (L.EQ.N) STOP 5
          L=L+1
          DO 402 M1=1,N
            A(M1,M)=A(M1,M)+A(M1,L)
            B(M1,M)=B(M1,M)+B(M1,L)
  402     CONTINUE
        ENDIF
        IF (A(M,M).EQ.0.) GOTO 200
C
C------------------------------------------------------------------
CL 3. DIVIDE ROW BY A(M,M) IN ORDER TO OBTAIN A(M,M)=1
C------------------------------------------------------------------
        AMM=A(M,M)
        DO 403 M1=1,N
          A(M1,M)=A(M1,M)/AMM
          B(M1,M)=B(M1,M)/AMM
  403   CONTINUE
C
C------------------------------------------------------------------
CL 4. IN COLUMN M MAKE EVERY VALUE 0 (EXCEPT A(M,M)) BY SUBTRACTING
CL    ROW M*THE APPROPRIATE VALUE.
C------------------------------------------------------------------
        DO 3 J=1,N
          IF (J.NE.M) THEN
            AMJ=A(M,J)
            DO 404 M1=1,N
              A(M1,J)=A(M1,J)-AMJ*A(M1,M)
              B(M1,J)=B(M1,J)-AMJ*B(M1,M)
  404       CONTINUE
          ENDIF
 3      CONTINUE
 2    CONTINUE
C
CL  EXTENDED MATRIX (A|B)=(I|C**(-1)) WHERE C**(-1)=INVERSE OF C
      RETURN
      END
