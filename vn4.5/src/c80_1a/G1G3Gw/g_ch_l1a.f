C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.
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
CLL  Routine: G_CH_L1A ------------------------------------------------
CLL
CLL  Purpose: To find the number of characters in a Fortran string
CLL           which are terminated by blanks.
CLL
CLL  Author:  Bob Carruthers, Cray Research.   Date: 18 September 1997
CLL           from the original code in UMSHELL1 by Paul Burton.
CLL
CLL Version Date      Modification history
CLL  4.5    09/10/98  Added def UTILHIST to top def line to allow
CLL                   history executables to be built. K Rogers
CLL
CLL  -------------------------------------------------------------------
C*L  Interface and arguments: ------------------------------------------
      INTEGER FUNCTION GET_CHAR_LEN(string)
! finds the length of the contents of character variable string

      IMPLICIT NONE

      CHARACTER *(*) string  ! IN : string to find length of
      INTEGER real_len
      LOGICAL found_end

      found_end=.FALSE.
      real_len=LEN(string)

      DO WHILE ( .NOT. found_end )

        IF (string(real_len:real_len) .NE. " ") THEN
          found_end=.TRUE.
        ELSE
          real_len=real_len-1
        ENDIF

      ENDDO

      GET_CHAR_LEN=real_len

      RETURN
      END
