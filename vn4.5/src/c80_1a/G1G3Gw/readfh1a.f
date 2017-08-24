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
CLL  SUBROUTINE READ_FLH --------------------------------------
CLL
CLL  Written by D. Robinson 17/06/92
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL
CLL  Programming standard: Unified Model Documentation Paper No 3
CLL                        Version No 4  5/2/92
CLL
CLL  System component: R30
CLL
CLL  System task: F3
CLL
CLL  Purpose:
CLL           Reads in the fixed length header from file attached to
CLL           unit NFTIN.
CLL
CLL  Documentation:
CLL           Unified Model Documentation Paper No F3
CLL           Version No 5 9/2/90
CLL
CLL------------------------------------------------------------
C*L Arguments:-------------------------------------------------
      SUBROUTINE READ_FLH (NFTIN,FIXHD,LEN_FIXHD,
     *           ICODE,CMESSAGE)

      IMPLICIT NONE

      INTEGER
     * NFTIN            ! IN  Unit no of dump
     *,LEN_FIXHD        ! IN  Length of fixed length header
     *,FIXHD(LEN_FIXHD) ! OUT Fixed length header

      INTEGER  ICODE  ! OUT Return code; successful=0, error > 0
      CHARACTER*(80)
     * CMESSAGE       !OUT Error message if ICODE > 0

C Local arrays:------------------------------------------------
C None
C -------------------------------------------------------------
C External subroutines called:---------------------------------
      EXTERNAL IOERROR,BUFFIN
C Local variables:---------------------------------------------
      INTEGER LEN_IO
      REAL A
C -------------------------------------------------------------

      ICODE=0
      CMESSAGE=' '

CL 1. Buffer in fixed length header record

      CALL BUFFIN (NFTIN,FIXHD(1),LEN_FIXHD,LEN_IO,A)

CL 2. Check for I/O errors
      IF(A.NE.-1.0.OR.LEN_IO.NE.LEN_FIXHD)THEN
        CALL IOERROR('buffer in of fixed length header',A,LEN_IO
     *               ,LEN_FIXHD)
        CMESSAGE='READ_FLH: I/O error'
        ICODE=1
      ENDIF

      RETURN
      END
