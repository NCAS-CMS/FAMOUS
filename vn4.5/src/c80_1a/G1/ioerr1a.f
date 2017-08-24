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
CLL  SUBROUTINE IOERROR----------------------------------------
CLL
CLL  Purpose: Prints out a message after using buffer in/out when
CLL           either a return code < 0.0 is encountered
CLL           by UNIT function or value returned by LENGTH
CLL           differs from length of I/O request.
CLL
CLL  Written by A. Dickinson
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL   4.1    12/06/96 Break up write statement. D. Robinson.
CLL   4.4    15/10/97 Added code to print the error message to
CLL                   stderr, and call abort in case all the
CLL                   PE's have not detected the error condition.
CLL                     Author: Bob Carruthers, Cray Research
CLL   4.5    08/07/98 Print only the leading non-blank 
CLL                   characters in 'string'
CLL                     Author: Bob Carruthers, Cray Research
CLL
CLL  Programming standard: Unified Model Documentation Paper No 3
CLL                        Version No 1 15/1/90
CLL
CLL  Logical component: E4
CLL
CLL  System task: F3
CLL
CLL  Documentation: CFT77 reference manual SR-0018 C  Page 9-3
CLL------------------------------------------------------------
C*L Arguments:-------------------------------------------------
      SUBROUTINE IOERROR(STRING,ERROR,LEN_IO,LEN_IO_REQ)

      IMPLICIT NONE

      INTEGER
     * LEN_IO  ! Number of 64-bit words transferred as registered
     *         ! by LENGTH function
     *,LEN_IO_REQ  ! Number of 64-bit words requested for
     *         ! transfer via BUFFER IN/OUT

      CHARACTER*(80) STRING ! User provided character string

      REAL
     * ERROR   ! Error code returned by UNIT function

      integer get_char_len   ! Returns the length of the string,
                             ! excluding trailing blanks

C -------------------------------------------------------------
C Workspace usage:---------------------------------------------
C None
C -------------------------------------------------------------
C*L External subroutines called:-------------------------------
C None
C*-------------------------------------------------------------

CL Internal structure: none

      WRITE(6,'('' **FATAL ERROR WHEN READING/WRITING MODEL DUMP**'')')
      WRITE(6,'('' '',A)') STRING(1:get_char_len(string))
      WRITE(6,'('' Error code = '',F6.2)') ERROR
      WRITE(6,'('' Length requested            = '',I9)') LEN_IO_REQ
      WRITE(6,'('' Length actually transferred = '',I9)') LEN_IO
      WRITE(6,'(''  Fatal error codes are as follows:'')')
      WRITE(6,'('' -1.0 Mismatch between actual and requested data'',
     *          '' length'')')
      WRITE(6,'(''  0.0 End-of-file was read'')')
      WRITE(6,'(''  1.0 Error occurred during read'')')
      WRITE(6,'(''  2.0 Other disk malfunction'')')
      WRITE(6,'('' 3.0 File does not exist'')')
      WRITE(6,'('' ***********************************************'')')

      RETURN
      END
