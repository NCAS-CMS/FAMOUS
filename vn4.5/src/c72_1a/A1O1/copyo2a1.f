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
      SUBROUTINE COPYO2A(IMT,JMAX,ODIN,WANT,OMASK
     &,ZERO,INVERT,IMAX,ODOUT)
CLL   SUBROUTINE COPYO2A --------------------------------------------
CLL
CLL        Purpose:
CLL   Auxiliary to TRANSO2A, used for transfer of fields from
CLL   ocean to atmosphere when the model grids are congruent, and as a
CLL   preliminary step when they are not congruent to remove missing
CLL   data and extra columns and invert the rows.
CLL
CLL   Transfer field ODIN on the input grid to ODOUT on the output grid
CLL   by straight copying, for the case where the grids are congruent.
CLL   The matrices may have different first dimension. If (INVERT), the
CLL   rows are in oppsite orders. Only those points .EQV.WANT (where
CLL   WANT=.FALSE. for sea) in the source mask are copied. If (ZERO),
CLL   the rest are set to zero in the output field; otherwise they are
CLL   not changed. The mask should be the same way up as ODOUT, opposite
CLL   to ODIN if (INVERT). This routine is appropriate where the target
CLL   field could define a value everywhere, but the source field does
CLL   not supply one everywhere.
CLL
CLL   WRITTEN BY J M GREGORY (31.5.91) (Extracted from TRANSO2A)
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL PROGRAMMING STANDARD :
CLL   FOLLOWS DOCUMENTATION PAPER 3, VERSION 1 FOR STANDARDS.
CLL
CLL LOGICAL COMPONENTS COVERED :
CLL
CLL PROJECT TASK :
CLL
CLL EXTERNAL DOCUMENTATION: UNIFIED MODEL DOCUMENTATION PAPER NO
CLL
CLLEND -----------------------------------------------------------------
CLL
C*L
      INTEGER
     & IMT                    !IN First dimension of ODIN
     &,JMAX                   !IN Second dimension of ODIN,ODOUT
     &,IMAX                   !IN First dimension of ODOUT
C
      REAL
     & ODIN(IMT,JMAX)         !IN Input field
     &,RMDI                   !IN Missing data value in input field
     &,ODOUT(IMAX,JMAX)       !INOUT Output field
C
      LOGICAL
     & WANT                   !IN Mark of required input points
     &,OMASK(IMT,JMAX)        !IN Input mask
     &,ZERO                   !IN Missing data to be replaced by zero
     &,INVERT                 !IN Row inversion is required
C*
      INTEGER
     & I,J,JI,JO              ! Loop indices
C
      IF (INVERT) THEN
        DO 50 JI = 1,JMAX
          JO=JMAX-JI+1
          DO 45 I = 1,IMAX
            IF (OMASK(I,JO).EQV.WANT) THEN
              ODOUT(I,JO) = ODIN(I,JI)
            ELSEIF (ZERO) THEN
              ODOUT(I,JO) = 0.0
            ENDIF
45        CONTINUE
50      CONTINUE
      ELSE
        DO 60 J = 1,JMAX
          DO 55 I = 1,IMAX
            IF (OMASK(I,J).EQV.WANT) THEN
              ODOUT(I,J) = ODIN(I,J)
            ELSEIF (ZERO) THEN
              ODOUT(I,J) = 0.0
            ENDIF
55        CONTINUE
60      CONTINUE
      ENDIF
C
      RETURN
      END
