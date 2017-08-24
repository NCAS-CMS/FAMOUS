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
CLL SUBROUTINE EXPT_ENC--------------------------------------------
CLL
CLL     Given a valid five character experiment RUN_ID, an INTEGER*4
CLL   code number is generated. The valid experiment name characters
CLL   are A-Z (uppercase), 0-9. Each letter of the experiment name is
CLL   stored as a 6 bit number. The last letter is converted to the
CLL   6 least significant bits of the integer code, the next letter
CLL   the next 6 lsb's, etc. Hence 30 bits are used in all.
CLL     The lookup table is capable of holding 64 elements. This
CLL   number cannot be exceeded if the code number is to remain
CLL   INTEGER*4. Similarly, the experiment RUN_ID length (5 chars)
CLL   cannot be exceeded.
CLL     Subroutine called from PP_HEAD.
CLL
CLL   A. Brady <- programmer of some or all of previous code or changes
CLL
CLL    Model            Modification history from model version 3.0:
CLL   version  Date
CLL
CLL   Programming standard:
CLL
CLL   Logical components covered:
CLL
CLL   Project TASK:
CLL
CLL   External documentation:
CLL
CLLEND-------------------------------------------------------------

C*L  INTERFACE and ARGUMENTS:--------------------------------------

      SUBROUTINE EXPT_ENC(EXPTSTRG
     &  ,EXPTCODE
     &  ,ICODE
     &  ,CMESSAGE)
C*-----------------------------------------------------------------

      IMPLICIT NONE

      INTEGER       ICODE       !OUT  Return code: successful=0
      CHARACTER*80  CMESSAGE    !OUT  Error message if ICODE > 0

      CHARACTER*5   EXPTSTRG    !IN   Experiment name string. Length
                                !     must equal parameter STRSIZE
      INTEGER       EXPTCODE    !OUT  Experiment code integer

C     Define local variables

      LOGICAL       TEST
      CHARACTER*1   LETTER
      INTEGER       I,J,NEWNUM,LETNUM,NSTRINGS,NBITS,STRSIZE

      PARAMETER(NSTRINGS=36,
     &  NBITS=6,
     &  STRSIZE=5)

      CHARACTER*1   USTRINGS    ! Upper case strings
      CHARACTER*1   LSTRINGS    ! Lower case strings

      DIMENSION     USTRINGS(0:NSTRINGS-1)
      DIMENSION     LSTRINGS(0:NSTRINGS-1)

      DATA USTRINGS/'A','B','C','D','E','F','G','H','I','J',
     &  'K','L','M','N','O','P','Q','R','S','T',
     &  'U','V','W','X','Y','Z','0','1','2','3',
     &  '4','5','6','7','8','9'/

      DATA LSTRINGS/'a','b','c','d','e','f','g','h','i','j',
     &  'k','l','m','n','o','p','q','r','s','t',
     &  'u','v','w','x','y','z','0','1','2','3',
     &  '4','5','6','7','8','9'/

C     Begin main

      EXPTCODE=0
      LETNUM=STRSIZE

C     Loop over letters in EXPTSTRG
      DO 20 I=0,STRSIZE-1
        TEST=.FALSE.
        READ(EXPTSTRG(LETNUM:LETNUM),"(A1)")LETTER

C       Loop over letters in lookup table USTRINGS/LSTRINGS
        DO 10 J=0,NSTRINGS-1
          IF ((LETTER.EQ.USTRINGS(J)).OR.(LETTER.EQ.LSTRINGS(J))) THEN
            TEST=.TRUE.
            NEWNUM=J*(2**(I*NBITS))
            EXPTCODE=EXPTCODE+NEWNUM
C           Exit loop as we have found the code for this letter
            GOTO 15
          ENDIF
 10     CONTINUE
 15     CONTINUE

C       Check experiment name is valid
        IF (.NOT.TEST) THEN
          ICODE=99
          CMESSAGE='EXPT_ENC: Invalid letter in expt name (RUN_ID)'
          RETURN
        ENDIF
        LETNUM=LETNUM-1
 20   CONTINUE

      RETURN
      END
