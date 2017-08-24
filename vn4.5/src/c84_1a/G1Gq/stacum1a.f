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
CLL  Routine: STACCUM  -------------------------------------------------
CLL
CLL  Purpose: Accumulates fields within STASH (temporal service routine)
CLL
CLL  Tested under compiler:   cft77
CLL  Tested under OS version: UNICOS 5.1
CLL
CLL  Author:   S.Tett/T.Johns
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL
CLL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
CLL
CLL  Logical components covered: D721
CLL
CLL  Project task: D7
CLL
CLL  External documentation:
CLL    Unified Model Doc Paper C4 - Storage handling and diagnostic
CLL                                 system (STASH)
CLL
C*L  Interface and arguments: ------------------------------------------
      SUBROUTINE STACCUM(fieldin,result,size,masking,amdi)
C
      IMPLICIT NONE
C
      INTEGER size              ! IN size of fieldin and result.
      REAL fieldin(size)        ! IN  input field to be processed
      REAL result(size)         ! OUT where accum is done.
      LOGICAL masking           ! IN true if masked (ie. MDI possible)
      REAL amdi                 ! IN missing data indicator
C* ---------------------------------------------------------------------
C
C  Local variables
C
      INTEGER i                 ! loop count
CL----------------------------------------------------------------------
CL 1.1 loop over array size, if either result or fieldin is amdi, set
CL    result to amdi else accumulate it in result.
CL
      IF (masking) THEN
        DO i=1,size
          IF ((result(i).ne.amdi).and.(fieldin(i).ne.amdi)) THEN
            result(i)=result(i)+fieldin(i)
          ELSE
            result(i)=amdi
          ENDIF
        ENDDO
      ELSE
CL
CL 1.2 loop over array size, accumulate result without checking for
CL     missing data
CL
        DO i=1,size
          result(i)=result(i)+fieldin(i)
        ENDDO
      ENDIF
C
      RETURN
      END
