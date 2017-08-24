CLL  Routine: EREPORT --------------------------------------------------
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
CLL
CLL  Purpose: Reports error exit code and message at end of model run.
CLL
CLL  Tested under compiler:   cft77
CLL  Tested under OS version: UNICOS 5.1
CLL
CLL  Author:   T.C.Johns
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL  3.4  10/10/94  Minor simplification of output format statements.
CLL                 R.Rawlins.
CLL  4.4  15/10/97  Added code to print the error message to
CLL                 stderr, and call abort in case all the
CLL                 PE's have not detected the error condition.
CLL                   Author: Bob Carruthers, Cray Research
CLL  4.5  08/07/98  Print only the leading non-blank 
CLL                 characters in 'cmessage'
CLL                   Author: Bob Carruthers, Cray Research
CLL  4.5  26/08/98  Changed resetting of ICODE  (A Van der Wal)
CLL
CLL  Programming standard: UM Doc Paper 1, version 1 (15/1/90)
CLL
CLL  Logical components covered: C0
CLL
CLL  Project task: C0
CLL
CLL  External documentation: On-line UM document C0 - The top-level
CLL                          control system
CLL
CLL  -------------------------------------------------------------------
C*L  Interface and arguments: ------------------------------------------
C
      SUBROUTINE EREPORT (ICODE,CMESSAGE)
      IMPLICIT NONE
      INTEGER ICODE            ! In - Error code from model
      CHARACTER*256 CMESSAGE   ! In - Error message from model

      integer get_char_len   ! Returns the length of the string,
                             ! excluding trailing blanks
C
C*----------------------------------------------------------------------
C  Local variables
C
CL----------------------------------------------------------------------
CL 1. Write informative message summarising completion state of model
CL
      WRITE(6,1000)
      IF (ICODE.LT.0) WRITE(6,1010) ICODE,
     2 CMESSAGE(1:get_char_len(cmessage))
      IF (ICODE.GT.0) WRITE(6,1020) ICODE,
     2 CMESSAGE(1:get_char_len(cmessage))
      WRITE(6,1000)
CL 1.1  Reset error code
      IF (ICODE.LT.0) ICODE=0
C
 1000 FORMAT(" ****************************************",
     &       "*****************************************")
 1010 FORMAT(" Model completed with warning code - ",I4,
     &       " Routine and message:-",(/A80))
 1020 FORMAT(" Model aborted with error code - ",I4,
     &       " Routine and message:-",(/A80))
      RETURN
CL----------------------------------------------------------------------
      END
