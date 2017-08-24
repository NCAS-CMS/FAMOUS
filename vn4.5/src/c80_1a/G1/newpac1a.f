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
CLL  SUBROUTINE NEWPACK:--------------------------------------
CLL
CLL  Purpose: Packing codes stored in LOOKUP(21,K) & LOOKUP(39,K)
CLL           are changed from pre vn2.8 values to
CLL           specification required at release 2.8
CLL
CLL  Written by A. Dickinson 28/08/92
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.2    13/07/93 Tidyied up integer declararions
CLL
CLL Programming standard :
CLL
CLL Logical components covered :
CLL
CLL Project task :
CLL
CLL  Documentation: UM Documentation Paper F3
CLL
CLLEND -----------------------------------------------------------------
C
      SUBROUTINE NEWPACK
     1(LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP)

      IMPLICIT NONE

      INTEGER
     1 LEN1_LOOKUP
     1,LEN2_LOOKUP
     1,LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)


      INTEGER
     1 N1
     1,N2
     1,N3
     1,K

      DO K=1,LEN2_LOOKUP
        N1=0
        N2=0
        N3=0
        IF(LOOKUP(21,K).EQ.-2)N1=2
C Ocean field packed using index array
        IF(LOOKUP(21,K).GT.9.AND.LOOKUP(21,K).LT.100)THEN
          N2=1
          N3=LOOKUP(21,K)-10
        ENDIF
C Ocean field compressed using bit mask
        IF(LOOKUP(21,K).GT.99)THEN
          N2=2
          N3=LOOKUP(21,K)-100
        ENDIF
C Real field stored at land pts
        IF(LOOKUP(39,K).EQ.4)THEN
          LOOKUP(39,K)=1
          N2=2
          N3=1
        ENDIF
C Integer field stored at land pts
        IF(LOOKUP(39,K).EQ.5)THEN
          LOOKUP(39,K)=2
          N2=2
          N3=1
        ENDIF

        LOOKUP(21,K)=100*N3+10*N2+N1
        ENDDO

      RETURN
      END
