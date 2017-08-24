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
CLL  SUBROUTINE POSERROR---------------------------------------
CLL
CLL  Purpose:
CLL           Prints out a message when position of a data block as
CLL           pointed to by fixed length header differs from actual
CLL           position in model dump.
CLL
CLL  Written by A. Dickinson 29/12/89
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL   4.4    15/10/97 Added code to print the error message to
CLL                   stderr, and call abort in case all the
CLL                   PE's have not detected the error condition.
CLL                     Author: Bob Carruthers, Cray Research
CLL
CLL  Programming standard:
CLL           Unified Model Documentation Paper No 3
CLL           Version No 1 15/1/90
CLL
CLL  System component: E4
CLL
CLL  System task: F3
CLL
CLL  Documentation:
CLL           None
CLL------------------------------------------------------------
C*L Arguments:-------------------------------------------------
      SUBROUTINE POSERROR(STRING,START_BLOCK,HEAD_POS,HEAD_ADDRESS)

      IMPLICIT NONE

      INTEGER
     * START_BLOCK  !IN Actual position of data block
     *,HEAD_POS     !IN Position in FIXHD of pointer
     *,HEAD_ADDRESS !IN Position in file pointed to by FIXHD(HEAD_POS)

      CHARACTER*(80) STRING  !IN Description of block

C -------------------------------------------------------------
C Workspace usage:---------------------------------------------
C None
C -------------------------------------------------------------
C*L External subroutines called:-------------------------------
C None
C*-------------------------------------------------------------

CL Internal structure: none

      WRITE(6,'('' ******FATAL ERROR WHEN READING MODEL DUMP******'')')
      WRITE(6,'('' Conflict between start position of '',A)')STRING
      WRITE(6,'('' block and pointer in fixed length header: FIXHD('',
     *I3,'') ='',I9)')HEAD_POS,HEAD_ADDRESS
      WRITE(6,'('' Current position in file ='',I9,'' words in'')')
     *START_BLOCK
      WRITE(6,'('' ***********************************************'')')

      RETURN
      END

