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
C*LL
CLL   SUBROUTINE ROWSWAP
CLL   ------------------
CLL
CLL   THIS ROUTINE INVERTS THE ORDER OF THE ROWS IN A TWO-DIMENSIONAL
CLL   ARRAY, FACILITATING COUPLING AN OCEAN MODEL IN WHICH THE FIRST
CLL   ROW IS AT THE SOUTH EDGE WITH AN ATMOSPHERE MODEL IN WHICH THE
CLL   ROWS START AT THE NORTH EDGE. IT IS CALLED BY TRANSO2A AND BY
CLL   TRANSA2O.
CLL   A WORK ARRAY IS USED IN ORDER THAT THE ROUTINE CAN BE CALLED
CLL   WITH THE SAME FIELD AS INPUT AND OUTPUT ARGUMENTS, IF REQUIRED.
CLL
CLL   ROUTINE WRITTEN BY D.L.ROBERTS
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL   THIS IS SYSTEM COMPONENT S193 (PART OF TASK D2).
CLL
CLL   THIS ROUTINE CAN BE COMPILED BY CFT77 BUT DOES NOT CONFORM TO
CLL   FORTRAN77 STANDARDS, BECAUSE OF THE INLINE COMMENTS. IT FOLLOWS
CLL   VERSION 1 OF DOCUMENTATION PAPER NO. 3.
CLL
CLLEND
C
      SUBROUTINE ROWSWAP(UP,DOWN,ICOLS,JROWS)
C     ---------------------------------------
C
      IMPLICIT NONE
C*L
      INTEGER ICOLS          ! IN    NUMBER OF COLUMNS IN FIELD
      INTEGER JROWS          ! IN    NUMBER OF ROWS IN FIELD
      REAL UP(ICOLS,JROWS)   ! IN    INPUT ARRAY
      REAL DOWN(ICOLS,JROWS) ! OUT   OUTPUT ARRAY
C*
      INTEGER
     + I,J                   !       LOOP COUNTERS
C
      REAL WORK(ICOLS,JROWS) !       WORK ARRAY
C
      DO 200 J = 1,JROWS
        DO 100 I = 1,ICOLS
          WORK(I,J) = UP(I,J)
100     CONTINUE
200   CONTINUE
C
      DO 400 J = 1,JROWS
        DO 300 I = 1,ICOLS
          DOWN(I,J) = WORK(I,JROWS+1-J)
300     CONTINUE
400   CONTINUE
C
      RETURN
      END
