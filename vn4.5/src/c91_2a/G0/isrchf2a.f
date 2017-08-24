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
CLL Function ISRCHFLT
CLL
CLL Purpose:  Portable version of Cray library function to find the
CLL           first real array element in relation to a real target.
CLL
CLL Tested under compiler:   fort77
CLL Tested under OS version: HP-UX A.08.07
CLL
CLL  Model            Modification history :
CLL version  Date
CLL  3.2   16/07/93   New deck. Tracey Smith.
CLL  3.3   22/09/93   Improved comments  Tracey Smith
CLL  4.0   06/12/95   Fixed incorrect loop bounds  P.Burton
CLL  4.4   24/04/97   Returns the index of the first value less than
CLL                   TARGET within the field of selected points.
CLL                                                 Ian Edmond
CLL
CLL Programming Standard: UM Doc Paper 3, version 5 (08/12/92)
CLL
      INTEGER FUNCTION ISRCHFLT(N,ARRAY,INC,TARGET)
      IMPLICIT NONE
      INTEGER
     &  N               ! IN number of elements to be searched
     & ,INC             ! IN increment between elements in the array
     & ,I               ! loop counter
      REAL
     &  ARRAY(1+(N-1)*INC) ! IN array to be searched
     & ,TARGET          ! IN real value to be searched for
      ISRCHFLT=N+1
      IF(N.LE.0) THEN
        ISRCHFLT=0
      ELSE
        DO 100 I=1,N
          IF(ARRAY(1+INC*(I-1)).LT.TARGET) THEN
            ISRCHFLT=I
            RETURN
          END IF
 100    CONTINUE
      END IF
      END
