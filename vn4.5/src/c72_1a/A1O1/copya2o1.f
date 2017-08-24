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
      SUBROUTINE COPYA2O(IMAX,JMAX,ODIN,INVERT,IMT,WANT,OMASK,ODOUT)
CLL
CLL   SUBROUTINE COPYA2O ----------------------------------------------
CLL
CLL   Auxiliary to TRANSA2O, used for transfer of fields from
CLL   atmosphere to ocean when the model grids are congruent, and as a
CLL   preliminary step when they are not congruent to remove missing
CLL   data and invert the rows.
CLL
CLL   Transfer field ODIN on the input grid to ODOUT on the output grid
CLL   by straight copying, for the case where the fields arecongruent.
CLL   The matrices may have different first dimension. If (INVERT), the
CLL   rows are in opposite orders. Only those points .EQV.WANT (where
CLL   WANT=.F. for sea) in the output mask are copied. Thus, this
CLL   routine is appropriate when the source field defines values every-
CLL   where, but the output field does not require a value everywhere.
CLL
CLL   WRITTEN BY J M GREGORY (1.7.91) (Extracted from TRANSA2O)
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL Programming standard :
CLL   Follows documentation paper 3, version 1 for standards.
CLL
CLL Logical components covered :
CLL
CLL Project task :
CLL
CLL Documentation:
CLL
CLLEND -------------------------------------------------------------
C*L
      INTEGER
     & IMAX                  !IN First dimension of ODIN
     &,JMAX                  !IN Second dimension of ODIN,ODOUT
     &,IMT                   !IN First dimension of ODOUT
C
      REAL
     & ODIN(IMAX,JMAX)       !IN Input field
     &,ODOUT(IMT,JMAX)       !INOUT Output field
C
      LOGICAL
     & INVERT                !IN Row inversion is required
     &,WANT                  !IN Logical value of selected points
     &,OMASK(IMT,JMAX)       !IN Output mask
C*
      INTEGER
     & I,J,JI,JO             ! Loop indices
C
      IF (INVERT) THEN
        DO 50 JI = 1,JMAX
          JO=JMAX-JI+1
          DO 45 I = 1,IMAX
            IF (OMASK(I,JO).EQV.WANT) ODOUT(I,JO)=ODIN(I,JI)
45        CONTINUE
50      CONTINUE
      ELSE
        DO 60 J = 1,JMAX
          DO 55 I = 1,IMAX
            IF (OMASK(I,J).EQV.WANT) ODOUT(I,J) = ODIN(I,J)
55        CONTINUE
60      CONTINUE
      ENDIF
C
      RETURN
      END
