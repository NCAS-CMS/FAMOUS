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
CLL  Subroutine: DEL_HIST -------------------------------------------
CLL
CLL  Purpose: delete a history file -called if problems writing
CLL           out partial sums
CLL
CLL  Tested under compiler:   cft77
CLL  Tested under OS version: UNICOS 6.1.5A
CLL
CLL  Author:   R A Stratton
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL   3.1  05/02/93    Portable Fortran unit no assigns
CLL                    Author: A. Dickinson    Reviewer: R. Stratton
CLL
CLL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
CLL
CLL  Logical components covered: H
CLL
CLL  Project task: H
CLL
CLL  External documentation:
CLL
CLLEND------------------------------------------------------------------
C
C*L  Interface and arguments: ------------------------------------------
      SUBROUTINE DEL_HIST(HUNIT)
C
      IMPLICIT NONE
C
      INTEGER
     &       HUNIT    ! IN  - history file unit number
C
C NOTE no error code and error message are passed back from this routine
C as this routine is only called in the event of a non-zero error
C code in U_MODEL.
C*----------------------------------------------------------------------
C  Common blocks
C
C  Subroutines called
C
      EXTERNAL GET_FILE
C
C  Local variables
C
      INTEGER   ICODE        ! error from open and close
      CHARACTER*80 CMESSAGE  ! error message this routine only
      CHARACTER*80 FILENAME
C
CL----------------------------------------------------------------------
CL 1. Open  history file
CL
      CALL GET_FILE(HUNIT,FILENAME,80,ICODE)
        OPEN(HUNIT,FILE=FILENAME,FORM='UNFORMATTED',IOSTAT=ICODE)
      IF (ICODE.NE.0) THEN
        CMESSAGE='DELHIST: failure to open history file prior to its del
     &etion'
        WRITE(6,*)CMESSAGE
      ENDIF
CL
CL 2. Close history file deleting it in the process.
CL
      CLOSE(HUNIT,IOSTAT=ICODE,STATUS='DELETE')
      IF (ICODE.NE.0) THEN
        CMESSAGE='DELHIST: failure to delete history file'
        WRITE(6,*)CMESSAGE
      ENDIF

      RETURN
CL----------------------------------------------------------------------
      END
C
