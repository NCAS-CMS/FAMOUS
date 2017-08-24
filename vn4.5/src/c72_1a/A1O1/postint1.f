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
      SUBROUTINE POST_H_INT(NCOASTAL,INDEXS,INDEXI,NS,FIELDS,SMINT
     &,NI,FIELDI,TPOINT,RMDI,NT,FIELDT)
CLL
CLL   Subroutine POST_H_INT
CLL
CLL       Purpose
CLL   POST_H_INT performs coastal adjustment and transfers the field
CLL   onto the target grid.
CLL
CLL   In the first step, values at points on the source grid (the
CLL   input to H_INT) are used to replace output values at coastal
CLL   points identified by COAST_AJ, provided the source grid point
CLL   selected by COAST_AJ is a sea point on that grid. Otherwise
CLL   missing data is substituted. In the second step, the final
CLL   set of values are written to the appropriate locations in the
CLL   target grid. (The interpolation routine H_INT does not use a
CLL   target grid explicitly, but a set of target points instead.)
CLL   Missing data values are not transferred.
CLL
CLL   WRITTEN BY J M GREGORY (3.6.91)
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL   FOLLOWS DOCUMENTATION PAPER 3, VERSION 1 FOR STANDARDS.
CLL
C*L
      IMPLICIT NONE
C
      INTEGER
     + NCOASTAL             !IN NUMBER OF COASTAL POINTS
     +,INDEXS(*)            !IN coastal points on source grid
     +,INDEXI(*)            !IN coastal points in target sea set
     &,NS                   !IN number of points in source grid
     &,SMINT(NS)            !IN land-sea mask on source grid
     &,NI                   !IN number of target sea points
     &,TPOINT(NI)           !IN list of sea points in target grid
     &,NT                   !IN number of points in target grid
C
      REAL
     & FIELDS(NS)           !IN field on source grid
     &,FIELDI(NI)           !INOUT values at target sea points
     &,RMDI                 !IN missing data indicator
     &,FIELDT(NT)           !INOUT field on target grid
C*
      INTEGER
     + K                    ! Loop index
C
CDIR$ IVDEP
      DO K=1,NCOASTAL
        IF (SMINT(INDEXS(K)).EQ.0) THEN
          FIELDI(INDEXI(K))=FIELDS(INDEXS(K))

        ENDIF
      ENDDO
      DO K=1,NI
        IF (FIELDI(K).NE.RMDI) FIELDT(TPOINT(K))=FIELDI(K)
      ENDDO
C
      RETURN
      END
