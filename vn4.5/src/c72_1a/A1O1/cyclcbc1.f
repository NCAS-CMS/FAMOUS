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
CLL   SUBROUTINE CYCLICBC ---------------------------------------------
CLL   -------------------
CLL
CLL   THIS ROUTINE COPIES THE FIRST TWO COLUMNS OF A TWO-DIMENSIONAL
CLL   ARRAY TO THE LAST TWO COLUMNS, OVERWRITING ANY DATA THAT HAPPEN
CLL   TO BE IN THOSE COLUMNS. THE MOTIVATION FOR THIS IS THAT THE OCEAN
CLL   MODEL HAS TWO SUCH DUPLICATE COLUMNS WHEN IT IS WORKING WITH A
CLL   DOMAIN WITH CYCLICALLY CONTINUOUS EAST-WEST BOUNDARIES (SUCH AS
CLL   A GLOBAL MODEL OR A FRAM-TYPE CONFIGURATION).
CLL   THIS ROUTINE IS CALLED FROM TRANSA2O.
CLL
CLL   ROUTINE WRITTEN BY D.L.ROBERTS
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL Programming standard :
CLL   This routine can be compiled by cft77 but does not conform to
CLL   Fortran77 standards, because of the inline comments. it follows
CLL   version 1 of Documentation paper no. 3.
CLL
CLL Logical components covered : S194
CLL
CLL Project task : D2
CLL
CLL External documentation: Unified Model documentation paper No:
CLL                         Version:
CLL
CLLEND -----------------------------------------------------------------
C
      SUBROUTINE CYCLICBC(ARRAY,ICOLS,JROWS)
C     --------------------------------------
C
      IMPLICIT NONE
C*L
      INTEGER ICOLS           ! IN  TOTAL NUMBER OF COLUMNS IN FIELD,
     +                        !     INCLUDING THE DUPLICATED COLUMNS.
      INTEGER JROWS           ! IN  NUMBER OF ROWS IN FIELD.
      REAL ARRAY(ICOLS,JROWS) ! IN OUT ARRAY TO BE OPERATED ON.
C*
      INTEGER
     + ICOLSM1,               !   THE PENULTIMATE COLUMN.
     + J                      !   LOOP COUNTER.
C
      ICOLSM1 = ICOLS - 1
C
      DO 100 J = 1,JROWS
        ARRAY(ICOLSM1,J) = ARRAY(1,J)
        ARRAY(ICOLS,J) = ARRAY(2,J)
100   CONTINUE
C
      RETURN
      END
